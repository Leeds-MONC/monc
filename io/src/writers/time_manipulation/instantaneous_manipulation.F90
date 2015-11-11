!> Performs instantaneous time manipulation and only returns a value if the output frequency determines one should be
module instantaneous_time_manipulation_mod
  use datadefn_mod, only : DEFAULT_PRECISION, STRING_LENGTH
  use collections_mod, only : hashmap_type, c_put_real, c_get_real, c_contains
  use conversions_mod, only : conv_single_real_to_double
  use forthread_mod, only : forthread_mutex_init, forthread_mutex_lock, forthread_mutex_unlock, forthread_mutex_destroy
  use threadpool_mod, only : check_thread_status 
  use configuration_parser_mod, only : data_values_type
  implicit none

#ifndef TEST_MODE
  private
#endif

  integer, volatile ::  existing_instantaneous_writes_mutex
  type(hashmap_type), volatile :: existing_instantaneous_writes

  public init_instantaneous_manipulation, finalise_instantaneous_manipulation, perform_instantaneous_time_manipulation
contains

  !> Initialises the instantaneous time manipulation
  subroutine init_instantaneous_manipulation()
    call check_thread_status(forthread_mutex_init(existing_instantaneous_writes_mutex, -1))
  end subroutine init_instantaneous_manipulation

  !> Finalises the instantaneous time manipulation
  subroutine finalise_instantaneous_manipulation()
    call check_thread_status(forthread_mutex_destroy(existing_instantaneous_writes_mutex))
  end subroutine finalise_instantaneous_manipulation

  !> Performs the instantaneous time manipulation and returns data only if this is to be written to the 
  !! storage. Internally a state is maintained which tracks when the write was last done to allow for flexibility in the
  !! time criteria.
  !! @param instant_values The instantaneous values to work with
  !! @param output_frequency The output frequency configuration option
  !! @param field_name The field name
  !! @param timestep The timestep
  !! @param time The model time
  !! @returns An allocated array of reals if data is to be stored, otherwise this is unallocated
  type(data_values_type) function perform_instantaneous_time_manipulation(instant_values, output_frequency, &
       field_name, timestep, time)
    real(kind=DEFAULT_PRECISION), dimension(:), intent(in) :: instant_values
    real, intent(in) :: output_frequency, time
    character(len=*), intent(in) :: field_name
    integer, intent(in) :: timestep

    if (deduce_whether_to_issue_values(field_name, output_frequency, time)) then
      allocate(perform_instantaneous_time_manipulation%values(size(instant_values)))
      perform_instantaneous_time_manipulation%values=instant_values
    end if    
  end function perform_instantaneous_time_manipulation

  !> Determines whether to issue values for write or not. This depends on the time and output frequency
  !! @param field_name The field name that we are manipulating
  !! @param output_frequency Configured output time frequency
  !! @param time The current model time
  !! @returns Whether or not one should issue the instantaneous values to write
  logical function deduce_whether_to_issue_values(field_name, output_frequency, time)
    character(len=*), intent(in) :: field_name
    real, intent(in) :: output_frequency, time

    real :: previous_time_write, time_difference

    call check_thread_status(forthread_mutex_lock(existing_instantaneous_writes_mutex))
    if (c_contains(existing_instantaneous_writes, field_name)) then
      previous_time_write=real(c_get_real(existing_instantaneous_writes, field_name))
      time_difference=time-previous_time_write
    else
      ! Rethink this as only works if started at time=0
      time_difference=time
    end if
    if (time_difference .ge. output_frequency) then
      deduce_whether_to_issue_values=.true.
      call c_put_real(existing_instantaneous_writes, field_name, conv_single_real_to_double(time))
    else
      deduce_whether_to_issue_values=.false.
    end if
    call check_thread_status(forthread_mutex_unlock(existing_instantaneous_writes_mutex))
  end function deduce_whether_to_issue_values  
end module instantaneous_time_manipulation_mod
