!> Broadcast inter IO communication which sends a value from one IO server to all others. This tracks field name and timestep
!! and only issues one call (and one results call to completion) for that combination
module broadcast_inter_io_mod
  use datadefn_mod, only : DEFAULT_PRECISION, DOUBLE_PRECISION, STRING_LENGTH
  use configuration_parser_mod, only : io_configuration_type, io_configuration_inter_communication_description
  use collections_mod, only : hashmap_type, list_type, iterator_type, mapentry_type, c_add_string, c_remove, c_free, &
       c_get_generic, c_get_string, c_put_generic, c_generic_at, c_get_iterator, c_has_next, c_next_mapentry, c_next_string, &
       c_is_empty
  use conversions_mod, only : conv_to_string
  use forthread_mod, only : forthread_mutex_init, forthread_mutex_lock, forthread_mutex_trylock, forthread_mutex_unlock, &
       forthread_mutex_destroy, forthread_rwlock_rdlock, forthread_rwlock_wrlock, forthread_rwlock_unlock, &
       forthread_rwlock_init, forthread_rwlock_destroy
  use threadpool_mod, only : check_thread_status, threadpool_start_thread
  use threadpool_mod, only : check_thread_status
  use inter_io_specifics_mod, only : handle_completion, register_inter_io_communication, find_inter_io_from_name, &
       package_inter_io_communication_message, unpackage_inter_io_communication_message
  use mpi, only : MPI_DOUBLE_PRECISION, MPI_INT, MPI_ANY_SOURCE, MPI_REQUEST_NULL, MPI_STATUSES_IGNORE, MPI_CHARACTER, MPI_BYTE
  implicit none

#ifndef TEST_MODE
  private
#endif

  integer, parameter :: MY_INTER_IO_TAG=2, PERFORM_CLEAN_EVERY=200
  character(len=*), parameter :: MY_INTER_IO_NAME="bcastinterio"

  !< Type keeping track of broadcast statuses
  type inter_io_broadcast
     logical :: handled
     integer :: mutex
     integer, dimension(:), allocatable :: send_requests
     character, dimension(:), allocatable :: send_buffer
     real(kind=DEFAULT_PRECISION), dimension(:), allocatable :: cached_values
     procedure(handle_completion), pointer, nopass :: completion_procedure
  end type inter_io_broadcast

  type(hashmap_type), volatile :: broadcast_statuses
  integer, volatile :: broadcast_statuses_rwlock, inter_io_description_mutex, clean_progress_mutex, &
       bcast_count_mutex, bcast_clean_reduction_count, bcast_count
  logical, volatile :: initialised=.false.

  public init_broadcast_inter_io, perform_inter_io_broadcast, finalise_broadcast_inter_io, check_broadcast_inter_io_for_completion
contains

  !> Initialises the broadcast inter IO functionality
  !! @param io_configuration The IO server configuration
  subroutine init_broadcast_inter_io(io_configuration)
    type(io_configuration_type), intent(inout) :: io_configuration

    if (.not. initialised) then
      initialised=.true.
      bcast_count=0
      bcast_clean_reduction_count=0
      call check_thread_status(forthread_rwlock_init(broadcast_statuses_rwlock, -1))
      call check_thread_status(forthread_mutex_init(inter_io_description_mutex, -1))
      call check_thread_status(forthread_mutex_init(clean_progress_mutex, -1))
      call check_thread_status(forthread_mutex_init(bcast_count_mutex, -1))
      call register_inter_io_communication(io_configuration, MY_INTER_IO_TAG, handle_recv_data_from_io_server, MY_INTER_IO_NAME)
    end if
  end subroutine init_broadcast_inter_io

  !> Finalises the broadcast inter IO functionality
  subroutine finalise_broadcast_inter_io()
    integer :: ierr
    type(inter_io_broadcast), pointer :: broadcast_item
    type(iterator_type) :: iterator

    if (initialised) then
      call check_thread_status(forthread_rwlock_rdlock(broadcast_statuses_rwlock))
      if (.not. c_is_empty(broadcast_statuses)) then
        iterator=c_get_iterator(broadcast_statuses)
        do while (c_has_next(iterator))
          broadcast_item=>retrieve_broadcast_item(c_next_mapentry(iterator))
          call check_thread_status(forthread_mutex_lock(broadcast_item%mutex))
          if (allocated(broadcast_item%send_requests)) then
            call mpi_waitall(size(broadcast_item%send_requests), broadcast_item%send_requests, MPI_STATUSES_IGNORE, ierr)
            deallocate(broadcast_item%send_requests)
            if (allocated(broadcast_item%send_buffer)) deallocate(broadcast_item%send_buffer)
          end if
          call check_thread_status(forthread_mutex_unlock(broadcast_item%mutex))
          call check_thread_status(forthread_mutex_destroy(broadcast_item%mutex))
        end do
      end if
      call check_thread_status(forthread_rwlock_unlock(broadcast_statuses_rwlock))
      call check_thread_status(forthread_rwlock_destroy(broadcast_statuses_rwlock))
      call check_thread_status(forthread_mutex_destroy(inter_io_description_mutex))
      call check_thread_status(forthread_mutex_destroy(clean_progress_mutex))
      call check_thread_status(forthread_mutex_destroy(bcast_count_mutex))
      initialised=.false.
    end if
  end subroutine finalise_broadcast_inter_io

  !> Checks the statuses for broadcast completion and returns whether they are all finished or not
  !! @param io_configuration The IO server configuration
  !! @returns Whether the broadcast inter IO has finished (i.e. no messages in transit)
  logical function check_broadcast_inter_io_for_completion(io_configuration)
    type(io_configuration_type), intent(inout) :: io_configuration

    type(inter_io_broadcast), pointer :: broadcast_item
    type(iterator_type) :: iterator

    check_broadcast_inter_io_for_completion=.true.
    call check_thread_status(forthread_rwlock_rdlock(broadcast_statuses_rwlock))
    if (.not. c_is_empty(broadcast_statuses)) then
      iterator=c_get_iterator(broadcast_statuses)
      do while (c_has_next(iterator))
        broadcast_item=>retrieve_broadcast_item(c_next_mapentry(iterator))
        if (.not. broadcast_item%handled) then
          check_broadcast_inter_io_for_completion=.false.
          exit
        end if
      end do
    end if
    call check_thread_status(forthread_rwlock_unlock(broadcast_statuses_rwlock))
  end function check_broadcast_inter_io_for_completion

  !> Handles receiving data from another IO server and processing this. If the request has already been registered (with a 
  !! callback) then this simply calls out. Otherwise it has to cache the data and awaits a thread calling the broadcast
  !! to call out to the callback
  !! @param io_configuration The IO server configuration
  !! @param data_buffer Data received from other IO server
  !! @param inter_io_index Index of the inter IO communication description
  subroutine handle_recv_data_from_io_server(io_configuration, data_buffer, inter_io_index)
    type(io_configuration_type), intent(inout) :: io_configuration
    character, dimension(:), intent(inout) :: data_buffer
    integer, intent(in) :: inter_io_index    

    type(inter_io_broadcast), pointer :: broadcast_item
    real(kind=DEFAULT_PRECISION), dimension(:), allocatable :: data_values
    character(len=STRING_LENGTH) :: field_name
    integer :: timestep

    call unpackage_inter_io_communication_message(data_buffer, field_name, timestep, data_values)

    broadcast_item=>find_or_add_broadcast_item(field_name, timestep)
    
    if (associated(broadcast_item%completion_procedure)) then
      call check_thread_status(forthread_mutex_lock(broadcast_item%mutex))
      broadcast_item%handled=.true.
      call check_thread_status(forthread_mutex_unlock(broadcast_item%mutex))
      call broadcast_item%completion_procedure(io_configuration, data_values, field_name, timestep)      
    else
      call check_thread_status(forthread_mutex_lock(broadcast_item%mutex))
      allocate(broadcast_item%cached_values(size(data_values)), source=data_values)
      broadcast_item%cached_values=data_values
      call check_thread_status(forthread_mutex_unlock(broadcast_item%mutex))
    end if
    if (allocated(data_values)) deallocate(data_values)
  end subroutine handle_recv_data_from_io_server
  
  !> Performs an inter IO broadcast of data from the root to all other IO servers. Note that this is on the IO server (and not
  !! MONC level) so might require some translation between the user's logical view and this view. Broadcasts are only issued once
  !! for a specific field_name and timestep pair.
  !! @param io_configuration Configuration of the IO server
  !! @param field_values The values to communicate
  !! @param field_size Number of elements to communicate
  !! @param field_name Field name that the reduction will be performed over
  !! @param root The root IO server process
  !! @param timestep The timestep this is issued at
  !! @param completion_procedure Callback completion procedure
  subroutine perform_inter_io_broadcast(io_configuration, field_values, field_size, field_name, root, &
       timestep, completion_procedure)
    type(io_configuration_type), intent(inout) :: io_configuration
    real(kind=DOUBLE_PRECISION), dimension(:) :: field_values
    integer, intent(in) :: field_size, root, timestep    
    character(len=*), intent(in) :: field_name
    procedure(handle_completion) :: completion_procedure
    
    type(inter_io_broadcast), pointer :: broadcast_item
    integer :: inter_io_comm_index, i, ierr

    call clean_broadcast_progress_if_needed()
    inter_io_comm_index=find_inter_io_from_name(io_configuration, MY_INTER_IO_NAME)
    broadcast_item=>find_or_add_broadcast_item(field_name, timestep, completion_procedure)    
    
    call check_thread_status(forthread_mutex_lock(broadcast_item%mutex))
    if (io_configuration%my_io_rank == root .and. .not. broadcast_item%handled) then
      broadcast_item%handled=.true.

      allocate(broadcast_item%send_requests(io_configuration%number_of_io_servers))
      broadcast_item%send_buffer=package_inter_io_communication_message(field_name, timestep, field_values)

      do i=0, io_configuration%number_of_io_servers-1
        if (i .ne. io_configuration%my_io_rank) then
          call mpi_isend(broadcast_item%send_buffer, size(broadcast_item%send_buffer), MPI_BYTE, i, &
               io_configuration%inter_io_communications(inter_io_comm_index)%message_tag, &
               io_configuration%io_communicator, broadcast_item%send_requests(i+1), ierr)
        else
          broadcast_item%send_requests(i+1)=MPI_REQUEST_NULL
        end if
      end do
      ! Still call the completion procedure on the root      
      call completion_procedure(io_configuration, field_values, field_name, timestep)
    else
      if (allocated(broadcast_item%cached_values) .and. .not. broadcast_item%handled) then
        broadcast_item%handled=.true.       
        call completion_procedure(io_configuration, broadcast_item%cached_values, field_name, timestep)
        if (allocated(broadcast_item%cached_values)) deallocate(broadcast_item%cached_values)
      end if
    end if
    call check_thread_status(forthread_mutex_unlock(broadcast_item%mutex))
  end subroutine perform_inter_io_broadcast

  !> Calls out to do a broadcast progress clean if needed (i.e. every n steps.)
  subroutine clean_broadcast_progress_if_needed()
    call check_thread_status(forthread_mutex_lock(bcast_count_mutex))
    bcast_count=bcast_count+1
    if (bcast_clean_reduction_count + PERFORM_CLEAN_EVERY .lt. bcast_count) then
      bcast_clean_reduction_count=bcast_count
      call check_thread_status(forthread_mutex_unlock(bcast_count_mutex))   
      call clean_broadcast_progress()
    else
      call check_thread_status(forthread_mutex_unlock(bcast_count_mutex))   
    end if
  end subroutine clean_broadcast_progress_if_needed  

  !> Performs a clean of the broadcast progresses that no longer need to be stored
  subroutine clean_broadcast_progress()
    type(inter_io_broadcast), pointer :: specific_broadcast_item_at_index
    integer :: completion_flag, ierr, num_to_remove, have_lock
    character(len=STRING_LENGTH) :: entry_key
    type(list_type) :: entries_to_remove
    logical :: destroy_lock
    type(iterator_type) :: iterator
    type(mapentry_type) :: mapentry
    class(*), pointer :: generic

    have_lock=forthread_mutex_trylock(clean_progress_mutex)
    if (have_lock==0) then
      call check_thread_status(forthread_rwlock_rdlock(broadcast_statuses_rwlock))
      iterator=c_get_iterator(broadcast_statuses)
      do while (c_has_next(iterator))
        destroy_lock=.false.
        mapentry=c_next_mapentry(iterator)
        specific_broadcast_item_at_index=>retrieve_broadcast_item(mapentry)
        call check_thread_status(forthread_mutex_lock(specific_broadcast_item_at_index%mutex))
        if (allocated(specific_broadcast_item_at_index%send_requests)) then
          call mpi_testall(size(specific_broadcast_item_at_index%send_requests), specific_broadcast_item_at_index%send_requests, &
               completion_flag, MPI_STATUSES_IGNORE, ierr)
          if (completion_flag == 1) then
            deallocate(specific_broadcast_item_at_index%send_requests)
            if (allocated(specific_broadcast_item_at_index%send_buffer)) deallocate(specific_broadcast_item_at_index%send_buffer)
            call c_add_string(entries_to_remove, mapentry%key)
            destroy_lock=.true.
          end if
        else
          if (specific_broadcast_item_at_index%handled) then
            if (allocated(specific_broadcast_item_at_index%cached_values)) then
              deallocate(specific_broadcast_item_at_index%cached_values)
            end if
            call c_add_string(entries_to_remove, mapentry%key)
            destroy_lock=.true.
          end if
        end if
        call check_thread_status(forthread_mutex_unlock(specific_broadcast_item_at_index%mutex))
        if (destroy_lock) call check_thread_status(forthread_mutex_destroy(specific_broadcast_item_at_index%mutex))
      end do      
      call check_thread_status(forthread_rwlock_unlock(broadcast_statuses_rwlock))

      if (.not. c_is_empty(entries_to_remove)) then
        call check_thread_status(forthread_rwlock_wrlock(broadcast_statuses_rwlock))
        iterator=c_get_iterator(entries_to_remove)
        do while (c_has_next(iterator))
          entry_key=c_next_string(iterator)
          generic=>c_get_generic(broadcast_statuses, entry_key)
          call c_remove(broadcast_statuses, entry_key)
          deallocate(generic)
        end do
        call check_thread_status(forthread_rwlock_unlock(broadcast_statuses_rwlock))
      end if
      call c_free(entries_to_remove)
      call check_thread_status(forthread_mutex_unlock(clean_progress_mutex))
    end if    
  end subroutine clean_broadcast_progress  

  !> Locates and returns or adds and returns a specific broadcast item representing a timestep and field
  !! @param field_name The field name this represents
  !! @param timestep The timestep this represents
  !! @param completion_procedure The (optional) completion procedure which is called once values are received
  !! @returns The existing or new broadcast item
  function find_or_add_broadcast_item(field_name, timestep, completion_procedure)
    character(len=*), intent(in) :: field_name
    integer, intent(in) :: timestep
    procedure(handle_completion), optional :: completion_procedure
    type(inter_io_broadcast), pointer :: find_or_add_broadcast_item

    class(*), pointer :: generic

    find_or_add_broadcast_item=>find_broadcast_item(field_name, timestep, .true.)
    if (.not. associated(find_or_add_broadcast_item)) then
      call check_thread_status(forthread_rwlock_wrlock(broadcast_statuses_rwlock))
      find_or_add_broadcast_item=>find_broadcast_item(field_name, timestep, .false.)
      if (.not. associated(find_or_add_broadcast_item)) then
        allocate(find_or_add_broadcast_item)
        if (present(completion_procedure)) then
          find_or_add_broadcast_item%completion_procedure=>completion_procedure
        else
          find_or_add_broadcast_item%completion_procedure=>null()
        end if
        find_or_add_broadcast_item%handled=.false.
        call check_thread_status(forthread_mutex_init(find_or_add_broadcast_item%mutex, -1))
        generic=>find_or_add_broadcast_item
        call c_put_generic(broadcast_statuses, trim(field_name)//"#"//conv_to_string(timestep), generic, .false.)
      end if
      call check_thread_status(forthread_rwlock_unlock(broadcast_statuses_rwlock))
    end if
  end function find_or_add_broadcast_item
  
  !> Finds a specific broadcast item or null if none is found
  !! @param field_name Corresponding field name to find
  !! @param timestep Corresponding timestep to find
  !! @param do_read_lock Whether to issue a read lock or not
  !! @returns The corresponding broadcast status item or null if none is found
  function find_broadcast_item(field_name, timestep, do_read_lock)
    character(len=*), intent(in) :: field_name
    integer, intent(in) :: timestep
    logical, intent(in) :: do_read_lock
    type(inter_io_broadcast), pointer :: find_broadcast_item

    class(*), pointer :: generic

    if (do_read_lock) call check_thread_status(forthread_rwlock_rdlock(broadcast_statuses_rwlock))
    generic=>c_get_generic(broadcast_statuses, trim(field_name)//"#"//conv_to_string(timestep))
    if (do_read_lock) call check_thread_status(forthread_rwlock_unlock(broadcast_statuses_rwlock))

    if (associated(generic)) then
      select type(generic)
        type is (inter_io_broadcast)
          find_broadcast_item=>generic
      end select
    else
      find_broadcast_item=>null()
    end if    
  end function find_broadcast_item

  !> Locates a broadcast item within a mapentry or null if none exists
  !! @param mapentry The map entry to use for this retrieval
  !! @returns The broadcast status or null if none exists
  function retrieve_broadcast_item(mapentry)   
    type(mapentry_type), intent(in) :: mapentry
    type(inter_io_broadcast), pointer :: retrieve_broadcast_item

    class(*), pointer :: generic

    generic=>c_get_generic(mapentry)

    if (associated(generic)) then
      select type(generic)
        type is (inter_io_broadcast)
          retrieve_broadcast_item=>generic
      end select
    else
      retrieve_broadcast_item=>null()
    end if    
  end function retrieve_broadcast_item
end module broadcast_inter_io_mod
