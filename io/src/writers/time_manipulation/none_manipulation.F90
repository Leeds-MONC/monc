!> Performs no time manipulation and returns the value, basically a no-op
module none_time_manipulation_mod
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

  public perform_none_time_manipulation, is_none_time_manipulation_ready_to_write
contains

  logical function is_none_time_manipulation_ready_to_write(latest_time, output_frequency, write_time, &
       latest_timestep, write_timestep)
    real, intent(in) :: latest_time, output_frequency, write_time
    integer, intent(in) :: latest_timestep, write_timestep

    is_none_time_manipulation_ready_to_write=latest_timestep .ge. write_timestep
  end function is_none_time_manipulation_ready_to_write  

  !> Performs no time manipulation and returns data
  !! @param instant_values The instantaneous values to work with
  !! @param output_frequency The output frequency configuration option
  !! @param field_name The field name
  !! @param timestep The timestep
  !! @param time The model time
  !! @param time_basis True for diagnostics interval in time coordinates, False for timestep coordinates
  !! @returns An allocated array of reals if data is to be stored, otherwise this is unallocated
  type(data_values_type) function perform_none_time_manipulation(instant_values, output_frequency, &
       field_name, timestep, time, time_basis)
    real(kind=DEFAULT_PRECISION), dimension(:), intent(in) :: instant_values
    real, intent(in) :: output_frequency
    real(kind=DEFAULT_PRECISION), intent(in) :: time
    character(len=*), intent(in) :: field_name
    integer, intent(in) :: timestep
    logical, intent(in) :: time_basis

    integer :: i

    allocate(perform_none_time_manipulation%values(size(instant_values)))
    do i=1,size(instant_values)
      perform_none_time_manipulation%values(i)=instant_values(i)
    end do
  end function perform_none_time_manipulation 
end module none_time_manipulation_mod
