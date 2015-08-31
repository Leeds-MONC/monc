!> Writer types which are shared across writing functionality
module writer_types_mod
  use collections_mod, only : queue_type, map_type
  use datadefn_mod, only : DEFAULT_PRECISION, STRING_LENGTH
  use configuration_parser_mod, only : data_values_type
  implicit none

#ifndef TEST_MODE
  private
#endif

  abstract interface
     !> Time manipulation interface which is implemented by the instantaneous and time averaged manipulations
     type(data_values_type) function perform_time_manipulation(instant_values, output_frequency, field_name, timestep, time)
       import DEFAULT_PRECISION, data_values_type
       real(kind=DEFAULT_PRECISION), dimension(:), intent(in) :: instant_values
       real, intent(in) :: output_frequency, time
       character(len=*), intent(in) :: field_name
       integer, intent(in) :: timestep
     end function perform_time_manipulation
  end interface

  !< Pending writes which will be dealt with sequentially
  type pending_write_type
     integer :: timestep
     real :: write_time
  end type pending_write_type  

  !< Field values stored for multiple MONCS (collective fields)
  type write_field_collective_values_type
     type(map_type) :: monc_values
  end type write_field_collective_values_type  

  !< The field type, many of these make up a specific writer
  type writer_field_type
     character(len=STRING_LENGTH) :: field_name, dim_size_defns(4), units
     procedure(perform_time_manipulation), pointer, nopass :: time_manipulation
     integer :: time_manipulation_type, values_mutex, dimensions, field_type, data_type, timestep_frequency, &
          max_timeseries_points_per_dump, actual_dim_size(4)
     real :: output_frequency, previous_write_time, previous_tracked_write_point
     logical :: collective_write, pending_to_write, enabled
     type(map_type) :: values_to_write
     logical :: duplicate_field_name
  end type writer_field_type
  
  !< A writer which will write to a file and contains many fields.
  type writer_type
     character(len=STRING_LENGTH) :: filename
     type(writer_field_type), dimension(:), allocatable :: contents
     integer :: trigger_and_write_mutex, write_timestep, num_fields_to_write, num_fields_to_write_mutex, pending_writes_mutex
     real :: write_time_frequency, previous_write_time, latest_pending_write_time, write_time
     logical :: currently_writing
     type(queue_type) :: pending_writes
  end type writer_type

  public writer_type, writer_field_type, write_field_collective_values_type, pending_write_type, &
       perform_time_manipulation
end module writer_types_mod
