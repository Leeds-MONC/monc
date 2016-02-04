!> Performs time averaged, time manipulation and only returns a value if the output frequency determines one should be
module timeaveraged_time_manipulation_mod
  use datadefn_mod, only : DEFAULT_PRECISION, STRING_LENGTH
  use configuration_parser_mod, only : io_configuration_type
  use collections_mod, only : hashmap_type, c_put_generic, c_get_generic
  use forthread_mod, only : forthread_rwlock_rdlock, forthread_rwlock_wrlock, forthread_rwlock_unlock, &
       forthread_rwlock_init, forthread_rwlock_destroy, forthread_mutex_init, forthread_mutex_lock, &
       forthread_mutex_unlock, forthread_mutex_destroy
  use threadpool_mod, only : check_thread_status
  use configuration_parser_mod, only : data_values_type
  implicit none

#ifndef TEST_MODE
  private
#endif

  !> The completed time averaged values
  type time_averaged_completed_type
     character(len=STRING_LENGTH) :: field_name
     real(kind=DEFAULT_PRECISION) :: start_time, previous_time, previous_output_time
     integer :: mutex
     logical :: empty_values
     real(kind=DEFAULT_PRECISION), dimension(:), allocatable :: time_averaged_values
  end type time_averaged_completed_type

  type(hashmap_type), volatile :: timeaveraged_values
  integer, volatile :: timeaveraged_value_rw_lock

  public init_time_averaged_manipulation, finalise_time_averaged_manipulation, perform_timeaveraged_time_manipulation
contains

   !> Initialises the reduction action
  subroutine init_time_averaged_manipulation()
    call check_thread_status(forthread_rwlock_init(timeaveraged_value_rw_lock, -1))
  end subroutine init_time_averaged_manipulation  

  !> Finalises the reduction action, waiting for all outstanding requests and then freeing data
  !! @param io_configuration Configuration state of the IO server
  subroutine finalise_time_averaged_manipulation()
    call check_thread_status(forthread_rwlock_destroy(timeaveraged_value_rw_lock))
  end subroutine finalise_time_averaged_manipulation

  !> Performs the time averaged manipulation and only returns values if these are to be stored (i.e. past an output frequency)
  !! @param instant_values The instantaneous values to work with
  !! @param output_frequency The output frequency configuration option
  !! @param field_name The field name
  !! @param timestep The timestep
  !! @param time The model time
  !! @returns An allocated array of reals if data is to be stored, otherwise this is unallocated
  type(data_values_type) function perform_timeaveraged_time_manipulation(instant_values, output_frequency, &
       field_name, timestep, time)
    real(kind=DEFAULT_PRECISION), dimension(:), intent(in) :: instant_values
    real, intent(in) :: output_frequency, time
    character(len=*), intent(in) :: field_name
    integer, intent(in) :: timestep

    type(time_averaged_completed_type), pointer :: timeaveraged_value

    timeaveraged_value=>find_or_add_timeaveraged_value(timestep, field_name)

    call check_thread_status(forthread_mutex_lock(timeaveraged_value%mutex))
    call time_average(timeaveraged_value, instant_values, time)

    if (time-timeaveraged_value%previous_output_time .ge. output_frequency) then
      timeaveraged_value%previous_output_time=time
      allocate(perform_timeaveraged_time_manipulation%values(size(timeaveraged_value%time_averaged_values)))
      perform_timeaveraged_time_manipulation%values=timeaveraged_value%time_averaged_values
      timeaveraged_value%time_averaged_values=0.0_DEFAULT_PRECISION
      timeaveraged_value%start_time=time
      timeaveraged_value%previous_time=time
      timeaveraged_value%empty_values=.true.
    end if
    call check_thread_status(forthread_mutex_unlock(timeaveraged_value%mutex))
  end function perform_timeaveraged_time_manipulation

  !> Does the time averaging itself
  !! @param timeaveraged_value The time averaged value to update
  !! @param instant_values The instant values to integrate in
  !! @param time The model time
  subroutine time_average(timeaveraged_value, instant_values, time)
    type(time_averaged_completed_type), intent(inout) :: timeaveraged_value
    real(kind=DEFAULT_PRECISION), dimension(:), intent(in) :: instant_values
    real, intent(in) :: time

    integer :: i
    real(kind=DEFAULT_PRECISION) :: timeav, timedg, combined_add    

    timeav=time-timeaveraged_value%start_time
    timedg=time-timeaveraged_value%previous_time
    combined_add=timeav+timedg

    if (.not. allocated(timeaveraged_value%time_averaged_values)) then
      allocate(timeaveraged_value%time_averaged_values(size(instant_values)))
      timeaveraged_value%time_averaged_values=0.0_DEFAULT_PRECISION
    end if
    
    if (timeaveraged_value%empty_values) then
      timeaveraged_value%empty_values=.false.
      timeaveraged_value%time_averaged_values=instant_values
    else
      do i=1, size(instant_values)
        timeaveraged_value%time_averaged_values(i)=(timeav*timeaveraged_value%time_averaged_values(i)+&
             timedg*instant_values(i)) / combined_add
      end do
    end if
    
    timeaveraged_value%previous_time=time
  end subroutine time_average  

  !> Retrieves or creates (and retrieves) a time averaged value based upon the information provided
  !! @param timestep The corresponding timestep
  !! @param field_name The corresponding field name
  !! @returns A matching or new time averaged value
  function find_or_add_timeaveraged_value(timestep, field_name)
    integer, intent(in) :: timestep
    character(len=*), intent(in) :: field_name
    type(time_averaged_completed_type), pointer :: find_or_add_timeaveraged_value
    
    class(*), pointer :: generic
    type(time_averaged_completed_type), pointer :: new_entry

    find_or_add_timeaveraged_value=>find_timeaveraged_value(field_name)
    if (.not. associated(find_or_add_timeaveraged_value)) then
      call check_thread_status(forthread_rwlock_wrlock(timeaveraged_value_rw_lock))
      find_or_add_timeaveraged_value=>find_timeaveraged_value(field_name, .false.)
      if (.not. associated(find_or_add_timeaveraged_value)) then
        allocate(new_entry)
        new_entry%field_name=field_name
        new_entry%start_time=0.0_DEFAULT_PRECISION
        new_entry%previous_time=0.0_DEFAULT_PRECISION
        new_entry%empty_values=.true.
        new_entry%previous_output_time=0.0_DEFAULT_PRECISION
        call check_thread_status(forthread_mutex_init(new_entry%mutex, -1))
        generic=>new_entry
        call c_put_generic(timeaveraged_values, field_name, generic, .false.)
        find_or_add_timeaveraged_value=>new_entry
      end if
      call check_thread_status(forthread_rwlock_unlock(timeaveraged_value_rw_lock))
    end if
  end function find_or_add_timeaveraged_value  

  !> Finds a time averaged value based upon its field name
  !! @param field_name The corresponding field name
  !! @param issue_read_lock Optional flag whether we should issue a read lock during this operation, if omitted then issues lock
  !! @returns The matching time averaged value or null if none is found
  function find_timeaveraged_value(field_name, issue_read_lock)
    character(len=*), intent(in) :: field_name
    type(time_averaged_completed_type), pointer :: find_timeaveraged_value
    logical, intent(in), optional :: issue_read_lock

    class(*), pointer :: generic
    logical :: do_read_lock

    if (present(issue_read_lock)) then
      do_read_lock=issue_read_lock      
    else
      do_read_lock=.true.
    end if        

    if (do_read_lock) call check_thread_status(forthread_rwlock_rdlock(timeaveraged_value_rw_lock))
    generic=>c_get_generic(timeaveraged_values, field_name)    
    if (do_read_lock) call check_thread_status(forthread_rwlock_unlock(timeaveraged_value_rw_lock))
    if (associated(generic)) then
      select type(generic)
      type is (time_averaged_completed_type)      
        find_timeaveraged_value=>generic
      end select
    else
      find_timeaveraged_value=>null()
    end if
  end function find_timeaveraged_value
end module timeaveraged_time_manipulation_mod
