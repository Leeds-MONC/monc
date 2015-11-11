!> All reduction, which does a reduce and then broadcasts the data to all IO servers
module allreduction_inter_io_mod
  use datadefn_mod, only : DEFAULT_PRECISION, DOUBLE_PRECISION, STRING_LENGTH
  use configuration_parser_mod, only : io_configuration_type
  use collections_mod, only : hashmap_type, c_get_generic, c_put_generic, c_remove, c_is_empty
  use conversions_mod, only : conv_to_string
  use forthread_mod, only : forthread_rwlock_rdlock, forthread_rwlock_wrlock, forthread_rwlock_unlock, &
       forthread_rwlock_init, forthread_rwlock_destroy
  use threadpool_mod, only : check_thread_status
  use inter_io_specifics_mod, only : handle_completion
  use reduction_inter_io_mod, only : init_reduction_inter_io, check_reduction_inter_io_for_completion, &
       finalise_reduction_inter_io, perform_inter_io_reduction, get_reduction_operator
  use broadcast_inter_io_mod, only : init_broadcast_inter_io, perform_inter_io_broadcast, finalise_broadcast_inter_io
  implicit none

#ifndef TEST_MODE
  private
#endif

  type allreduce_type
     integer :: root
     procedure(handle_completion), pointer, nopass :: completion_procedure
  end type allreduce_type  

  logical, volatile :: initialised=.false.

  integer, volatile :: allreduce_rwlock
  type(hashmap_type), volatile :: allreduce_types

  public init_allreduction_inter_io, finalise_allreduction_inter_io, perform_inter_io_allreduction, &
       check_allreduction_inter_io_for_completion
contains

  !> Initialises the all reduction inter IO functionality
  !! @param io_configuration The IO server configuration
  subroutine init_allreduction_inter_io(io_configuration)
    type(io_configuration_type), intent(inout) :: io_configuration

    if (.not. initialised) then
      initialised=.true.
      call check_thread_status(forthread_rwlock_init(allreduce_rwlock, -1))
    end if
    call init_reduction_inter_io(io_configuration)
    call init_broadcast_inter_io(io_configuration)
  end subroutine init_allreduction_inter_io

  !> Finalises the all reduction inter IO functionality
  !! @param io_configuration The IO server configuration
  subroutine finalise_allreduction_inter_io(io_configuration)
    type(io_configuration_type), intent(inout) :: io_configuration

    if (initialised) then
      initialised=.false.
      call check_thread_status(forthread_rwlock_destroy(allreduce_rwlock))
    end if
    call finalise_reduction_inter_io(io_configuration)
    call finalise_broadcast_inter_io()
  end subroutine finalise_allreduction_inter_io

  !> Determines whether this all reduction inter IO functionality has completed or not
  !! @param io_configuration The IO server configuration
  !! @returns Whether the inter IO communication has completed or not
  logical function check_allreduction_inter_io_for_completion(io_configuration)
    type(io_configuration_type), intent(inout) :: io_configuration

    call check_thread_status(forthread_rwlock_rdlock(allreduce_rwlock))
    check_allreduction_inter_io_for_completion=c_is_empty(allreduce_types)
    call check_thread_status(forthread_rwlock_unlock(allreduce_rwlock))
  end function check_allreduction_inter_io_for_completion
  
  !> Performs the all reduction inter IO reduction
  !! @param io_configuration Configuration of the IO server
  !! @param field_values The values to communicate
  !! @param field_size Number of elements to communicate
  !! @param reduction_field_name Field name that the reduction will be performed over
  !! @param reduction_op The reduction operator to use
  !! @param root The root IO server process
  !! @param timestep The timestep this is issued at
  !! @param completion_procedure Callback completion procedure
  subroutine perform_inter_io_allreduction(io_configuration, field_values, field_size, field_name, reduction_op, root, &
       timestep, completion_procedure)
    type(io_configuration_type), intent(inout) :: io_configuration
    real(kind=DOUBLE_PRECISION), dimension(:) :: field_values
    integer, intent(in) :: field_size, reduction_op, root, timestep    
    character(len=*), intent(in) :: field_name
    procedure(handle_completion) :: completion_procedure

    if (io_configuration%my_io_rank .eq. root) then
      call add_allreduce_information_if_needed(field_name, timestep, root, completion_procedure)
    end if
    call perform_inter_io_reduction(io_configuration, field_values, field_size, field_name, reduction_op, &
       root, timestep, internal_reduction_completion_procedure)
    if (io_configuration%my_io_rank .ne. root) then
      ! None root processes issue a broadcast
      call perform_inter_io_broadcast(io_configuration, field_values, size(field_values), field_name, root, &
           timestep, completion_procedure)
    end if    
  end subroutine perform_inter_io_allreduction  

  !> Internal completion, called after the reduce has completed (on root) and calls out to broadcast
  !! @param io_configuration Configuration of the IO server
  !! @param values Array of values resulting from the communication
  !! @param field_name Corresponding field name
  !! @param timestep Corresponding timestep
  subroutine internal_reduction_completion_procedure(io_configuration, values, field_name, timestep)
    type(io_configuration_type), intent(inout) :: io_configuration
    real(DEFAULT_PRECISION), dimension(:) :: values
    character(len=*) :: field_name
    integer :: timestep

    type(allreduce_type), pointer :: allreduce_information

    allreduce_information=>find_allreduce_information(field_name, timestep, .true.)    
    call perform_inter_io_broadcast(io_configuration, values, size(values), field_name, &
         allreduce_information%root, timestep, allreduce_information%completion_procedure)
    call remove_allreduce_information(field_name, timestep, .true.)
  end subroutine internal_reduction_completion_procedure

  !> Adds an all reduce information to the status if it does not exist
  !! @param field_name The corresponding field name
  !! @param timestep The corresponding timestep
  !! @param root The root process
  !! @param completion_procedure Procedure to call on communication completion
  subroutine add_allreduce_information_if_needed(field_name, timestep, root, completion_procedure)
    character(len=*), intent(in) :: field_name
    integer, intent(in) :: timestep, root
    procedure(handle_completion) :: completion_procedure

    type(allreduce_type), pointer :: allreduce_information
    class(*), pointer :: generic

    allreduce_information=>find_allreduce_information(field_name, timestep, .true.)
    if (.not. associated(allreduce_information)) then
      call check_thread_status(forthread_rwlock_wrlock(allreduce_rwlock))
      allreduce_information=>find_allreduce_information(field_name, timestep, .false.)
      if (.not. associated(allreduce_information)) then
        allocate(allreduce_information)
        allreduce_information%completion_procedure=>completion_procedure
        allreduce_information%root=root
        generic=>allreduce_information
        call c_put_generic(allreduce_types, trim(field_name)//"#"//conv_to_string(timestep), generic, .true.)
      end if
      call check_thread_status(forthread_rwlock_unlock(allreduce_rwlock))
    end if
  end subroutine add_allreduce_information_if_needed
  
  !> Finds an all reduce status information based on the field name and timestep, or returns null if none is found
  !! @param field_name The corresponding field name
  !! @param timestep The corresponding timestep
  !! @param dolock Whether to issue a read lock or not
  !! @returns The corresponding all reduce status or null if none is found
  function find_allreduce_information(field_name, timestep, dolock)
    character(len=*), intent(in) :: field_name
    integer, intent(in) :: timestep
    logical, intent(in) :: dolock
    type(allreduce_type), pointer :: find_allreduce_information

    class(*), pointer :: generic

    if (dolock) call check_thread_status(forthread_rwlock_rdlock(allreduce_rwlock))
    generic=>c_get_generic(allreduce_types, trim(field_name)//"#"//conv_to_string(timestep))
    if (dolock) call check_thread_status(forthread_rwlock_unlock(allreduce_rwlock))
    
    if (associated(generic)) then
      select type(generic)
        type is (allreduce_type)
          find_allreduce_information=>generic
      end select      
    else
      find_allreduce_information=>null()
    end if
  end function find_allreduce_information

  !> Removes an all reduce status information based on the field name and timestep
  !! @param field_name The corresponding field name
  !! @param timestep The corresponding timestep
  !! @param dolock Whether to issue a write lock or not
  subroutine remove_allreduce_information(field_name, timestep, dolock)
    character(len=*), intent(in) :: field_name
    integer, intent(in) :: timestep
    logical, intent(in) :: dolock

    if (dolock) call check_thread_status(forthread_rwlock_wrlock(allreduce_rwlock))
    call c_remove(allreduce_types, trim(field_name)//"#"//conv_to_string(timestep))
    if (dolock) call check_thread_status(forthread_rwlock_unlock(allreduce_rwlock))
  end subroutine remove_allreduce_information  
end module allreduction_inter_io_mod
