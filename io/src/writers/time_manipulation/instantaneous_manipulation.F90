!> Performs instantaneous time manipulation and only returns a value if the output frequency determines one should be
module instantaneous_time_manipulation_mod
  use datadefn_mod, only : DEFAULT_PRECISION, STRING_LENGTH
  use collections_mod, only : hashmap_type, mapentry_type, iterator_type, c_get_iterator, c_next_mapentry, c_put_real, &
       c_get_real, c_contains, c_size, c_has_next
  use conversions_mod, only : conv_single_real_to_double
  use forthread_mod, only : forthread_mutex_init, forthread_mutex_lock, forthread_mutex_unlock, forthread_mutex_destroy
  use threadpool_mod, only : check_thread_status 
  use configuration_parser_mod, only : data_values_type
  use data_utils_mod, only : unpack_scalar_integer_from_bytedata, unpack_scalar_string_from_bytedata, &
       unpack_scalar_dp_real_from_bytedata
  use io_server_client_mod, only : pack_scalar_field
  implicit none

#ifndef TEST_MODE
  private
#endif

  integer, volatile ::  existing_instantaneous_writes_mutex
  type(hashmap_type), volatile :: existing_instantaneous_writes

  public init_instantaneous_manipulation, finalise_instantaneous_manipulation, perform_instantaneous_time_manipulation, &
       is_instantaneous_time_manipulation_ready_to_write, serialise_instantaneous_state, unserialise_instantaneous_state, &
       prepare_to_serialise_instantaneous_state
contains

  !> Initialises the instantaneous time manipulation
  subroutine init_instantaneous_manipulation()
    call check_thread_status(forthread_mutex_init(existing_instantaneous_writes_mutex, -1))
  end subroutine init_instantaneous_manipulation

  !> Finalises the instantaneous time manipulation
  subroutine finalise_instantaneous_manipulation()
    call check_thread_status(forthread_mutex_destroy(existing_instantaneous_writes_mutex))
  end subroutine finalise_instantaneous_manipulation

  logical function is_instantaneous_time_manipulation_ready_to_write(latest_time, output_frequency, write_time, &
       latest_timestep, write_timestep)
    real, intent(in) :: latest_time, output_frequency, write_time
    integer, intent(in) :: latest_timestep, write_timestep

    is_instantaneous_time_manipulation_ready_to_write=latest_time + output_frequency .gt. write_time
  end function is_instantaneous_time_manipulation_ready_to_write

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
    real, intent(in) :: output_frequency
    real(kind=DEFAULT_PRECISION), intent(in) :: time
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
    real, intent(in) :: output_frequency
    real(kind=DEFAULT_PRECISION), intent(in) :: time

    real(kind=DEFAULT_PRECISION) :: previous_time_write, time_difference

    call check_thread_status(forthread_mutex_lock(existing_instantaneous_writes_mutex))
    if (c_contains(existing_instantaneous_writes, field_name)) then
      previous_time_write=c_get_real(existing_instantaneous_writes, field_name)
      time_difference=(aint(time*10000000.0)-aint(previous_time_write*10000000.0))/10000000.0
    else
      ! Rethink this as only works if started at time=0
      time_difference=time
    end if
    if (time_difference .ge. output_frequency) then
      deduce_whether_to_issue_values=.true.
      call c_put_real(existing_instantaneous_writes, field_name, time)
    else
      deduce_whether_to_issue_values=.false.
    end if
    call check_thread_status(forthread_mutex_unlock(existing_instantaneous_writes_mutex))
  end function deduce_whether_to_issue_values

  !> Prepares to serialise the instantaneous state, both determines the byte storage size and issues any locks
  !! @returns The number of bytes needed to store the serialised state
  integer(kind=8) function prepare_to_serialise_instantaneous_state()
    real(kind=DEFAULT_PRECISION) :: a
    
    call check_thread_status(forthread_mutex_lock(existing_instantaneous_writes_mutex))
    prepare_to_serialise_instantaneous_state=kind(prepare_to_serialise_instantaneous_state) + &
         ((STRING_LENGTH + kind(a)) * c_size(existing_instantaneous_writes))
  end function prepare_to_serialise_instantaneous_state  

  !> Will serialise the state of this manipulator so that it can be later restarted. Any locks issued during preparation
  !! are released here
  !! @param byte_data The state of this manipulator is written into here
  subroutine serialise_instantaneous_state(byte_data)
    character, dimension(:), allocatable, intent(inout) :: byte_data

    integer :: current_data_point
    type(mapentry_type) :: map_entry
    type(iterator_type) :: iterator
    
    current_data_point=1
    current_data_point=pack_scalar_field(byte_data, current_data_point, c_size(existing_instantaneous_writes))
    iterator=c_get_iterator(existing_instantaneous_writes)
    do while (c_has_next(iterator))
      map_entry=c_next_mapentry(iterator)
      current_data_point=pack_scalar_field(byte_data, current_data_point, string_value=map_entry%key)
      current_data_point=pack_scalar_field(byte_data, current_data_point, double_real_value=c_get_real(map_entry))
    end do
    call check_thread_status(forthread_mutex_unlock(existing_instantaneous_writes_mutex))
  end subroutine serialise_instantaneous_state

  !> Unpacks some serialised byte data to initialise this manipulator to some previous state
  !! @param byte_data The byte data to unpack and reinitialise from
  subroutine unserialise_instantaneous_state(byte_data)
    character, dimension(:), allocatable, intent(in) :: byte_data

    integer :: current_data_point, number_entries, i

    current_data_point=1
    number_entries=unpack_scalar_integer_from_bytedata(byte_data, current_data_point)
    if (number_entries .gt. 0) then
      do i=1, number_entries
        call c_put_real(existing_instantaneous_writes, unpack_scalar_string_from_bytedata(byte_data, current_data_point), &
             unpack_scalar_dp_real_from_bytedata(byte_data, current_data_point))        
      end do
    end if    
  end subroutine unserialise_instantaneous_state
end module instantaneous_time_manipulation_mod
