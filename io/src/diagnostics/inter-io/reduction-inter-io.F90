!> Reduction inter IO action which will perform reductions between IO servers. This is not as trivial as calling the MPI function
!! as it is nondeterministic when messages will arrive and hence when one reduction on a process and a reduction on another
!! should be called.
module reduction_inter_io_mod
  use datadefn_mod, only : DEFAULT_PRECISION, DOUBLE_PRECISION, STRING_LENGTH
  use configuration_parser_mod, only : io_configuration_type, io_configuration_inter_communication_description
  use collections_mod, only : map_type, hashmap_type, list_type, iterator_type, mapentry_type, c_get_generic, c_get_string, &
       c_add_string, c_put_generic, c_remove, c_is_empty, c_contains, c_free, c_generic_at, c_get_iterator, c_has_next, &
       c_next_mapentry, c_next_string
  use conversions_mod, only : conv_to_string
  use forthread_mod, only : forthread_mutex_init, forthread_mutex_lock, forthread_mutex_trylock, &
       forthread_mutex_unlock, forthread_mutex_destroy, forthread_rwlock_rdlock, forthread_rwlock_wrlock, &
       forthread_rwlock_unlock, forthread_rwlock_init, forthread_rwlock_destroy, forthread_rwlock_trywrlock
  use threadpool_mod, only : check_thread_status, threadpool_start_thread
  use logging_mod, only : LOG_ERROR, log_log
  use inter_io_specifics_mod, only : handle_completion, register_inter_io_communication, find_inter_io_from_name, &
       package_inter_io_communication_message, unpackage_inter_io_communication_message
  use mpi, only : MPI_DOUBLE_PRECISION, MPI_INT, MPI_ANY_SOURCE, MPI_REQUEST_NULL, MPI_STATUS_IGNORE, MPI_CHARACTER, MPI_BYTE
  use mpi_communication_mod, only : lock_mpi, unlock_mpi, wait_for_mpi_request
  implicit none

#ifndef TEST_MODE
  private
#endif

  !< The types of reduction operator that are supported
  integer, parameter :: MEAN=1, MIN=2, MAX=3, SUM=4
  integer, parameter :: MY_INTER_IO_TAG=1, PERFORM_CLEAN_EVERY=200
  character(len=*), parameter :: MY_INTER_IO_NAME="reductioninterio"

  type reduction_progress_type
     real(DEFAULT_PRECISION), dimension(:), allocatable :: values
     character, dimension(:), allocatable :: send_buffer
     character(len=STRING_LENGTH) :: field_name
     integer :: contributed_moncs, contributed_io_servers, timestep, reduction_operator, async_handle, mutex, root
     procedure(handle_completion), pointer, nopass :: completion_procedure
  end type reduction_progress_type  

  integer, volatile :: reduction_progress_rwlock, inter_io_description_mutex, clean_progress_mutex, &
       reduction_count_mutex, previous_clean_reduction_count, reduction_count
  type(hashmap_type), volatile :: reduction_progresses
  logical, volatile :: initialised=.false.

  public init_reduction_inter_io, check_reduction_inter_io_for_completion, finalise_reduction_inter_io, &
       perform_inter_io_reduction, get_reduction_operator
contains

  !> Initialises the reduction action
  !! @param io_configuration The IO server configuration
  subroutine init_reduction_inter_io(io_configuration)
    type(io_configuration_type), intent(inout) :: io_configuration

    if (.not. initialised) then
      initialised=.true.
      previous_clean_reduction_count=0
      reduction_count=0
      call check_thread_status(forthread_rwlock_init(reduction_progress_rwlock, -1))
      call check_thread_status(forthread_mutex_init(inter_io_description_mutex, -1))
      call check_thread_status(forthread_mutex_init(clean_progress_mutex, -1))
      call check_thread_status(forthread_mutex_init(reduction_count_mutex, -1))
      call register_inter_io_communication(io_configuration, MY_INTER_IO_TAG, handle_recv_data_from_io_server, MY_INTER_IO_NAME)
    end if
  end subroutine init_reduction_inter_io  

  !> Handles the receiving of data from some other IO server. This is issued call back style within a thread
  !! to handle that data
  !! @param io_configuration Configuration state of the IO server
  !! @param inter_io_index Index of the inter IO communication description
  subroutine handle_recv_data_from_io_server(io_configuration, data_buffer, inter_io_index)
    type(io_configuration_type), intent(inout) :: io_configuration
    character, dimension(:), intent(inout) :: data_buffer
    integer, intent(in) :: inter_io_index    

    call handle_process_recv_from_other_IO_server(io_configuration, io_configuration%inter_io_communications(inter_io_index), &
         io_configuration%my_io_rank, data_buffer, io_configuration%number_of_io_servers)
  end subroutine handle_recv_data_from_io_server

  !> Checks this action for completion, when all are completed then the IO server can shutdown as this is called
  !! once all MONC processes have deregistered
  !! @param io_configuration Configuration state of the IO server
  !! @returns Whether the action has completed
  logical function check_reduction_inter_io_for_completion(io_configuration)
    type(io_configuration_type), intent(inout) :: io_configuration

    check_reduction_inter_io_for_completion=check_and_clean_progress(io_configuration%my_io_rank)
  end function check_reduction_inter_io_for_completion  

  !> Finalises the reduction action, waiting for all outstanding requests and then freeing data
  !! @param io_configuration Configuration state of the IO server
  subroutine finalise_reduction_inter_io(io_configuration)
    type(io_configuration_type), intent(inout) :: io_configuration

    type(reduction_progress_type) :: progress
    type(iterator_type) :: iterator

    if (initialised) then
      call check_thread_status(forthread_rwlock_rdlock(reduction_progress_rwlock))
      if (.not. c_is_empty(reduction_progresses)) then
        iterator=c_get_iterator(reduction_progresses)
        do while (c_has_next(iterator))
          progress=retrieve_reduction_progress(c_next_mapentry(iterator))
          if (progress%async_handle /= MPI_REQUEST_NULL) then
            call wait_for_mpi_request(progress%async_handle)
          end if
        end do
      end if
      call check_thread_status(forthread_rwlock_unlock(reduction_progress_rwlock))
      call check_thread_status(forthread_rwlock_destroy(reduction_progress_rwlock))
      call check_thread_status(forthread_mutex_destroy(inter_io_description_mutex))
      call check_thread_status(forthread_mutex_destroy(clean_progress_mutex))
      call check_thread_status(forthread_mutex_destroy(reduction_count_mutex))
      initialised=.false.
    end if
  end subroutine finalise_reduction_inter_io

  !> Actually handles the processing for this data wrt the vertical reduction
  !! @param io_configuration Configuration of the IO server
  !! @param field_values The values to communicate
  !! @param field_size Number of elements to communicate
  !! @param reduction_field_name Field name that the reduction will be performed over
  !! @param reduction_op The reduction operator to use
  !! @param root The root IO server process
  !! @param timestep The timestep this is issued at
  !! @param completion_procedure Callback completion procedure
  subroutine perform_inter_io_reduction(io_configuration, field_values, field_size, reduction_field_name, reduction_op, &
       root, timestep, completion_procedure)
    type(io_configuration_type), intent(inout) :: io_configuration
    real(kind=DOUBLE_PRECISION), dimension(:) :: field_values
    integer, intent(in) :: field_size, reduction_op, root, timestep    
    character(len=*), intent(in) :: reduction_field_name
    procedure(handle_completion) :: completion_procedure

    type(reduction_progress_type), pointer :: reduction_progress
    logical :: collective_values_new
        
    call clean_progress(io_configuration%my_io_rank)
    reduction_progress=>find_or_add_reduction_progress(timestep, reduction_op, root, reduction_field_name, completion_procedure)

    call check_thread_status(forthread_mutex_lock(reduction_progress%mutex))
    reduction_progress%contributed_moncs=reduction_progress%contributed_moncs+1

    collective_values_new=.not. allocated(reduction_progress%values)
    if (collective_values_new) allocate(reduction_progress%values(field_size))

    call integrate_io_server_collective_values(reduction_op, reduction_progress, field_values, field_size, collective_values_new)
    if (reduction_progress%contributed_moncs == io_configuration%number_of_moncs) then
      reduction_progress%contributed_io_servers=reduction_progress%contributed_io_servers+1      
      call handle_local_moncs_completed_collective(io_configuration, reduction_progress)
    else
      call check_thread_status(forthread_mutex_unlock(reduction_progress%mutex))
    end if    
  end subroutine perform_inter_io_reduction

  subroutine clean_progress(myrank)
    integer, intent(in) :: myrank

    logical :: cc_dummy

    call check_thread_status(forthread_mutex_lock(reduction_count_mutex))  
    reduction_count=reduction_count+1
    if (previous_clean_reduction_count + PERFORM_CLEAN_EVERY .lt. reduction_count) then      
      previous_clean_reduction_count=reduction_count
      call check_thread_status(forthread_mutex_unlock(reduction_count_mutex))        
      cc_dummy=check_and_clean_progress(myrank)
    else
      call check_thread_status(forthread_mutex_unlock(reduction_count_mutex))      
    end if
  end subroutine clean_progress

  !> Checks all the reduction progresses and will remove any that have completed. This is designed to be
  !! called from an IO server other than 0 (the master IO server) and it checks if the outstanding async
  !! send handle has completed. Checking on the master IO server or checking any progress that is not currently
  !! sending is fine and will not impact the correctness (but obviously the progress wont be freed)
  logical function check_and_clean_progress(myrank)
    integer, intent(in) :: myrank

    integer :: i, entries, completed, ierr, num_to_remove, have_lock
    type(list_type) :: entries_to_remove
    type(iterator_type) :: iterator
    type(mapentry_type) :: mapentry
    type(reduction_progress_type), pointer :: specific_reduction_progress
    character(len=STRING_LENGTH) :: entry_key
    class(*), pointer :: generic
    logical :: destroy_lock

    check_and_clean_progress=.true.
    have_lock=forthread_mutex_trylock(clean_progress_mutex)
    if (have_lock == 0) then
      call check_thread_status(forthread_rwlock_rdlock(reduction_progress_rwlock))
      iterator=c_get_iterator(reduction_progresses)
      do while (c_has_next(iterator))
        mapentry=c_next_mapentry(iterator)
        destroy_lock=.false.
        specific_reduction_progress=>retrieve_reduction_progress(mapentry)
        if (myrank /= specific_reduction_progress%root) then
          call check_thread_status(forthread_mutex_lock(specific_reduction_progress%mutex))
          if (specific_reduction_progress%async_handle /= MPI_REQUEST_NULL) then
            call wait_for_mpi_request(specific_reduction_progress%async_handle)
            !if (completed == 1) then
            if (specific_reduction_progress%async_handle == MPI_REQUEST_NULL) then
              if (allocated(specific_reduction_progress%send_buffer)) deallocate(specific_reduction_progress%send_buffer)
              destroy_lock=.true.
              call c_add_string(entries_to_remove, mapentry%key)
            else
              check_and_clean_progress=.false.
            end if
          end if
          call check_thread_status(forthread_mutex_unlock(specific_reduction_progress%mutex))
          if (destroy_lock) call check_thread_status(forthread_mutex_destroy(specific_reduction_progress%mutex))
        end if
      end do
      call check_thread_status(forthread_rwlock_unlock(reduction_progress_rwlock))

      if (.not. c_is_empty(entries_to_remove)) then
        call check_thread_status(forthread_rwlock_wrlock(reduction_progress_rwlock))
        iterator=c_get_iterator(entries_to_remove)
        do while (c_has_next(iterator))
          entry_key=c_next_string(iterator)
          generic=>c_get_generic(reduction_progresses, entry_key)
          call c_remove(reduction_progresses, entry_key)
          deallocate(generic)
        end do
        call check_thread_status(forthread_rwlock_unlock(reduction_progress_rwlock))
      end if
      call c_free(entries_to_remove)
      call check_thread_status(forthread_mutex_unlock(clean_progress_mutex))
    end if
  end function check_and_clean_progress

  !> Integrates the collective values from another IO server into the currently stored values
  !! @param reduction_op The reduction operator to perform
  !! @param reduction_progress The progress data type which is updated
  !! @param single_server_values The values from the IO server which need to be integrated
  !! @param dim_one_size Size in first dimension
  !! @param target_index The index where we are putting the values into the current value array of reduction progress
  !! @param collective_values_empty Whether the collective values is empty
  subroutine integrate_io_server_collective_values(reduction_op, reduction_progress, single_server_values, &
       number_elements, collective_values_empty)
    integer, intent(in) :: reduction_op, number_elements
    logical, intent(in) :: collective_values_empty
    type(reduction_progress_type), intent(inout) :: reduction_progress
    real(kind=DOUBLE_PRECISION), dimension(:), intent(in) :: single_server_values

    integer :: k

    if (collective_values_empty) then
      reduction_progress%values=single_server_values
    else
      if (reduction_op == MEAN .or. reduction_op == SUM) then
        reduction_progress%values=reduction_progress%values+single_server_values       
      else if (reduction_op == MIN .or. reduction_op == MAX) then
        do k=1, number_elements
          if (reduction_op == MIN) then
            if (single_server_values(k) .lt. reduction_progress%values(k)) &
                 reduction_progress%values(k)=single_server_values(k)
          else if (reduction_op == MAX) then
            if (single_server_values(k) .gt. reduction_progress%values(k)) &
                 reduction_progress%values(k)=single_server_values(k)
          end if
        end do
      end if
    end if
  end subroutine integrate_io_server_collective_values

  !> Handles the case where the local MONC processes have completed their collective operation for a specific reduction
  !! and, for this IO server, it either needs to send its value to the master IO server or, if it is the master, check
  !! for completion
  !! @param io_configuration Configuration state of the IO server
  !! @param reduction_progress The specific reduction progress data item that represents this reduction
  !! @param z_size Size in Z
  !! @param reduction_progress_location Location in the reduction progresses list that this single progress item resides at
  subroutine handle_local_moncs_completed_collective(io_configuration, reduction_progress)
    type(io_configuration_type), intent(inout) :: io_configuration
    type(reduction_progress_type), intent(inout) :: reduction_progress

    integer :: ierr, inter_io_comm_index
    
    if (io_configuration%my_io_rank == reduction_progress%root .and. &
         reduction_progress%contributed_io_servers == io_configuration%number_of_io_servers) then
      call handle_collective_completed(io_configuration, reduction_progress)
    else
      if (io_configuration%my_io_rank /= reduction_progress%root) then
        inter_io_comm_index=find_inter_io_from_name(io_configuration, MY_INTER_IO_NAME)

        reduction_progress%send_buffer=package_inter_io_communication_message(reduction_progress%field_name, &
             reduction_progress%timestep, reduction_progress%values, reduction_progress%reduction_operator)
        call lock_mpi()
        call mpi_isend(reduction_progress%send_buffer, size(reduction_progress%send_buffer), &
             MPI_BYTE, reduction_progress%root, &
             io_configuration%inter_io_communications(inter_io_comm_index)%message_tag, &
             io_configuration%io_communicator, reduction_progress%async_handle, ierr)
        call unlock_mpi()
        ! Deallocate the current value as this is finished with and has been packed into the send buffer
        if (allocated(reduction_progress%values)) deallocate(reduction_progress%values)
      end if
      call check_thread_status(forthread_mutex_unlock(reduction_progress%mutex))
    end if
  end subroutine handle_local_moncs_completed_collective

  !> Handles the data received from another IO server, locates the correct reduction progress, appends the information
  !! and then checks for & deals with the situation where that reduction is completed
  !! @param number_io_servers The total number of IO servers
  !! @param z_size Number of levels in the vertical
  subroutine handle_process_recv_from_other_IO_server(io_configuration, inter_io_comm, myrank, data_buffer, number_io_servers)
    type(io_configuration_type), intent(inout) :: io_configuration
    type(io_configuration_inter_communication_description), intent(inout) :: inter_io_comm
    character, dimension(:), intent(inout) :: data_buffer
    integer, intent(in) :: number_io_servers, myrank

    type(reduction_progress_type), pointer :: reduction_progress
    character(len=STRING_LENGTH) :: field_name
    integer :: timestep, reduction_op
    real(kind=DEFAULT_PRECISION), dimension(:), allocatable :: field_values
    logical :: collective_values_new

    call unpackage_inter_io_communication_message(data_buffer, field_name, timestep, field_values, reduction_op)
    reduction_progress=>find_or_add_reduction_progress(timestep, reduction_op, myrank, field_name)

    call check_thread_status(forthread_mutex_lock(reduction_progress%mutex))
    collective_values_new=.not. allocated(reduction_progress%values)
    if (collective_values_new) allocate(reduction_progress%values(size(field_values)))

    reduction_progress%contributed_io_servers=reduction_progress%contributed_io_servers+1    
    call integrate_io_server_collective_values(reduction_op, reduction_progress, &
         field_values, size(field_values), collective_values_new)
    if (reduction_progress%contributed_io_servers == number_io_servers) then
      call handle_collective_completed(io_configuration, reduction_progress)
      deallocate(field_values)
      return
    end if
    deallocate(field_values)
    call check_thread_status(forthread_mutex_unlock(reduction_progress%mutex))
  end subroutine handle_process_recv_from_other_IO_server

  !> Handles the situation where collective communication for a specific reduction has completed across all IO servers
  !! @param reduction_progress The reduction progress data type
  !! @param number_io_servers The total number of IO servers
  subroutine handle_collective_completed(io_configuration, reduction_progress)
    type(io_configuration_type), intent(inout) :: io_configuration
    type(reduction_progress_type), intent(inout) :: reduction_progress
    
    if (reduction_progress%reduction_operator == MEAN) then
      reduction_progress%values=reduction_progress%values/io_configuration%number_of_global_moncs
    end if    
    call reduction_progress%completion_procedure(io_configuration, reduction_progress%values, &
         reduction_progress%field_name, reduction_progress%timestep)
    call check_thread_status(forthread_mutex_unlock(reduction_progress%mutex))
    call check_thread_status(forthread_mutex_destroy(reduction_progress%mutex))
    if (allocated(reduction_progress%values)) deallocate(reduction_progress%values)
    call remove_reduction_progress(reduction_progress)
  end subroutine handle_collective_completed

  !> Finds or adds a specific reduction progress based upon the timestep and reduction operator. If none can be found
  !! then a new progress is added in. With new progresses this procedure will initialise them
  !! @param timestep The timestep to match
  !! @param reduction_operator The reduction operator to match
  !! @param field_name The name of the field that the reduction type represents
  !! @param num_vectors The number of reduction vectors (items to reduce) to be stored
  !! @returns A reduction progress data object
  function find_or_add_reduction_progress(timestep, reduction_operator, root, field_name, completion_procedure)
    integer, intent(in) :: timestep, reduction_operator, root
    type(reduction_progress_type), pointer :: find_or_add_reduction_progress
    character(len=*), intent(in) :: field_name
    procedure(handle_completion), optional :: completion_procedure

    class(*), pointer :: generic
    type(reduction_progress_type), pointer :: new_progress

    find_or_add_reduction_progress=>find_reduction_progress(timestep, reduction_operator, field_name)
    if (.not. associated(find_or_add_reduction_progress)) then
      call check_thread_status(forthread_rwlock_wrlock(reduction_progress_rwlock))
      find_or_add_reduction_progress=>find_reduction_progress(timestep, reduction_operator, field_name, issue_read_lock=.false.)
      if (.not. associated(find_or_add_reduction_progress)) then
        allocate(new_progress)
        call check_thread_status(forthread_mutex_init(new_progress%mutex, -1))
        new_progress%timestep=timestep
        new_progress%reduction_operator=reduction_operator
        new_progress%contributed_moncs=0
        new_progress%contributed_io_servers=0
        new_progress%root=root
        new_progress%async_handle=MPI_REQUEST_NULL
        new_progress%field_name=field_name
        if (present(completion_procedure)) then
          new_progress%completion_procedure=>completion_procedure
        else
          new_progress%completion_procedure=>null()
        end if        
        generic=>new_progress
        call c_put_generic(reduction_progresses, generate_reduction_key(field_name, timestep, reduction_operator), &
             generic, .false.)
        find_or_add_reduction_progress=>new_progress
      end if
      call check_thread_status(forthread_rwlock_unlock(reduction_progress_rwlock))
    end if      
    if (.not. associated(find_or_add_reduction_progress%completion_procedure) .and. present(completion_procedure)) then
      find_or_add_reduction_progress%completion_procedure=>completion_procedure
    end if
  end function find_or_add_reduction_progress

  !> Locates a specific reduction progress based upon the timestep, operator and field name
  !! @param timestep The timestep to search for
  !! @param reduction_operator The reduction operator to search for
  !! @param field_name The field name which must match
  !! @param reduction_progress_location Optional location which is set to be the index of the matching progress item
  !! @returns Pointer to the reduction progress or null if none is found
  function find_reduction_progress(timestep, reduction_operator, field_name, issue_read_lock)
    integer, intent(in) :: timestep, reduction_operator    
    logical, intent(in), optional :: issue_read_lock
    type(reduction_progress_type), pointer :: find_reduction_progress
    character(len=*), intent(in) :: field_name

    class(*), pointer :: generic
    logical :: do_read_lock

    if (present(issue_read_lock)) then
      do_read_lock=issue_read_lock      
    else
      do_read_lock=.true.
    end if

    if (do_read_lock) call check_thread_status(forthread_rwlock_rdlock(reduction_progress_rwlock))
    generic=>c_get_generic(reduction_progresses, generate_reduction_key(field_name, timestep, reduction_operator))     
    if (do_read_lock) call check_thread_status(forthread_rwlock_unlock(reduction_progress_rwlock))
    if (associated(generic)) then
      select type(generic)
      type is (reduction_progress_type)      
        find_reduction_progress=>generic
      end select
    else
      find_reduction_progress=>null()
    end if
  end function find_reduction_progress

  !> Removes a specific reduction progress
  !! @param reduction_progress The reduction progress to remove from the list
  subroutine remove_reduction_progress(reduction_progress)
    type(reduction_progress_type), intent(in) :: reduction_progress
        
    class(*), pointer :: generic
    character(len=STRING_LENGTH) :: specific_key

    specific_key=generate_reduction_key(reduction_progress%field_name, reduction_progress%timestep,&
         reduction_progress%reduction_operator)
    call check_thread_status(forthread_rwlock_wrlock(reduction_progress_rwlock))
    generic=>c_get_generic(reduction_progresses, specific_key)
    call c_remove(reduction_progresses, specific_key)
    call check_thread_status(forthread_rwlock_unlock(reduction_progress_rwlock))
    if (associated(generic)) deallocate(generic)
  end subroutine remove_reduction_progress

  !> Generates the lookup key that is used for the map storage of reduction progresses
  !! @param field_name The field name
  !! @param timestep The timestep
  !! @param reduction_operator The reduction operator
  character(len=STRING_LENGTH) function generate_reduction_key(field_name, timestep, reduction_operator)
    character(len=*), intent(in) :: field_name
    integer, intent(in) :: timestep, reduction_operator

    generate_reduction_key=trim(field_name)//"#"//trim(conv_to_string(timestep))//"#"// trim(conv_to_string(reduction_operator))
  end function generate_reduction_key

  !> Helper function to retrieve the reduction progress from a mapentry
  !! @param mapentry The map entry to retrieve from
  !! @returns The progress data object in the map entry or null if none is found
  function retrieve_reduction_progress(mapentry)
    type(mapentry_type), intent(in) :: mapentry
    type(reduction_progress_type), pointer :: retrieve_reduction_progress

    class(*), pointer :: generic

    generic=>c_get_generic(mapentry)

    if (associated(generic)) then
      select type(generic)
      type is (reduction_progress_type)      
        retrieve_reduction_progress=>generic
      end select
    else
      retrieve_reduction_progress=>null()
    end if 
  end function retrieve_reduction_progress

  !> Given the map of action attributes this procedure will identify the reduction operator that has been
  !! selected by the configuration
  !! @param action_attributes Action attributes from the IO server configuration
  !! @returns The reduction operator
  integer function get_reduction_operator(op_string)
    character(len=*), intent(in) :: op_string
    
    if (op_string .eq. "mean") then
      get_reduction_operator=MEAN
    else if (op_string .eq. "min") then
      get_reduction_operator=MIN
    else if (op_string .eq. "max") then
      get_reduction_operator=MAX
    else if (op_string .eq. "sum") then
      get_reduction_operator=SUM
    else
      call log_log(LOG_ERROR, "Reduction operator '"//trim(op_string)//"' not recognised")
    end if
  end function get_reduction_operator  
end module reduction_inter_io_mod
