!> Writer types which are shared across writing functionality. Also includes serialisation functionality for these types
module writer_types_mod
  use collections_mod, only : mapentry_type, list_type, iterator_type, queue_type, map_type, hashmap_type, c_get_generic, &
       c_size, c_next_mapentry, c_has_next, c_get_iterator, c_put_generic, c_get_string, c_put_string, c_next_generic, &
       c_push_generic
  use datadefn_mod, only : DEFAULT_PRECISION, STRING_LENGTH
  use configuration_parser_mod, only : data_values_type
  use io_server_client_mod, only : pack_scalar_field, pack_array_field
  use data_utils_mod, only : unpack_scalar_integer_from_bytedata, unpack_scalar_real_from_bytedata
  use logging_mod, only : LOG_ERROR, log_log
  use forthread_mod, only : forthread_mutex_lock, forthread_mutex_unlock
  use threadpool_mod, only : check_thread_status
  implicit none

#ifndef TEST_MODE
  private
#endif

  abstract interface
     !> Time manipulation interface which is implemented by the instantaneous and time averaged manipulations
     type(data_values_type) function perform_time_manipulation(instant_values, output_frequency, field_name, timestep, &
                                                               time, time_basis)
       import DEFAULT_PRECISION, data_values_type
       real(kind=DEFAULT_PRECISION), dimension(:), intent(in) :: instant_values
       real, intent(in) :: output_frequency
       real(kind=DEFAULT_PRECISION), intent(in) :: time
       character(len=*), intent(in) :: field_name
       integer, intent(in) :: timestep
       logical, intent(in) :: time_basis
     end function perform_time_manipulation

     logical function is_field_ready_to_write(latest_time, output_frequency, write_time, latest_timestep, write_timestep)
       real, intent(in) :: latest_time, output_frequency, write_time
       integer, intent(in) :: latest_timestep, write_timestep
     end function is_field_ready_to_write     
  end interface

  !< Pending writes which will be dealt with sequentially
  type pending_write_type
     integer :: timestep
     real :: write_time
     logical :: terminated_write
  end type pending_write_type  

  !< Field values stored for multiple MONCS (collective fields)
  type write_field_collective_values_type
     type(map_type) :: monc_values
  end type write_field_collective_values_type

  !< Describes a collective contiguous block of data
  type write_field_collective_descriptor_type
     integer, dimension(:), allocatable :: absolute_start, count
     integer :: split_dim
     type(list_type) :: specific_monc_info
  end type write_field_collective_descriptor_type

  !< Specific information for a corresponding MONC process
  type write_field_collective_monc_info_type
     integer :: relative_dimension_start, counts(3), monc_location, monc_source
  end type write_field_collective_monc_info_type  

  !< The field type, many of these make up a specific writer
  type writer_field_type
     character(len=STRING_LENGTH) :: field_name, field_namespace, dim_size_defns(4), units
     procedure(perform_time_manipulation), pointer, nopass :: time_manipulation
     procedure(is_field_ready_to_write), pointer, nopass :: ready_to_write
     integer :: time_manipulation_type, values_mutex, dimensions, field_type, data_type, timestep_frequency, &
          actual_dim_size(4), latest_timestep_values, max_num_collective_writes, max_num_collective_writes_request_handle
     real :: output_frequency, previous_write_time, previous_tracked_write_point
     logical :: collective_write, collective_initialised, collective_contiguous_optimisation, &
          pending_to_write, enabled, expected_here, issue_write
     type(map_type) :: values_to_write
     type(list_type) :: collective_descriptors
     logical :: duplicate_field_name, prognostic_field, diagnostic_field
  end type writer_field_type
  
  !< A writer which will write to a file and contains many fields.
  type writer_type
     character(len=STRING_LENGTH) :: filename, title
     type(writer_field_type), dimension(:), allocatable :: contents
     integer :: trigger_and_write_mutex, write_timestep, previous_write_timestep, num_fields_to_write,              &
          num_fields_to_write_mutex, pending_writes_mutex, write_timestep_frequency, latest_pending_write_timestep, &
          write_precision
     real :: write_time_frequency, previous_write_time, latest_pending_write_time, write_time, defined_write_time
     logical :: write_on_model_time, contains_io_status_dump, write_on_terminate, include_in_io_state_write
     type(queue_type) :: pending_writes
  end type writer_type

  !< Represents the dimension information associated with a Q field that is written collectively
  type collective_q_field_representation_type
     integer, dimension(:), allocatable :: dimensions
  end type collective_q_field_representation_type

  !< Represents the dimension information associated with a tracer field that is written collectively
  type collective_tracer_representation_type
     integer, dimension(:), allocatable :: dimensions
  end type collective_tracer_representation_type

  type netcdf_diagnostics_timeseries_type
     integer :: netcdf_dim_id, netcdf_var_id, num_entries
     real :: last_write_point
     logical :: variable_written
  end type netcdf_diagnostics_timeseries_type

  !< Keeps track of a specific diagnostics NetCDF file
  type netcdf_diagnostics_type
     integer :: ncid, mutex, key_value_dim_id, string_dim_id
     type(map_type) :: dimension_to_id
     type(hashmap_type) :: variable_to_id, timeseries_dimension
     type(writer_type), pointer :: corresponding_writer_entry
     logical :: termination_write
  end type netcdf_diagnostics_type

  public writer_type, writer_field_type, write_field_collective_values_type, pending_write_type, &
       perform_time_manipulation, collective_q_field_representation_type,  collective_tracer_representation_type, &
       netcdf_diagnostics_timeseries_type, &
       netcdf_diagnostics_type, serialise_writer_type, unserialise_writer_type, serialise_data_values_type, &
       unserialise_data_values_type, write_field_collective_descriptor_type, write_field_collective_monc_info_type, &
       prepare_to_serialise_data_values_type, prepare_to_serialise_writer_type
contains

  !> Prepares to serialise the writer type by issuing locks and determining the size of serialised bytes needed
  !! @param writer_to_serialise Writer type to serialise
  !! @returns The number of bytes needed to store the serialised state
  integer(kind=8) function prepare_to_serialise_writer_type(writer_to_serialise)
    type(writer_type), intent(inout) :: writer_to_serialise

    integer :: i

    prepare_to_serialise_writer_type=(kind(writer_to_serialise%write_timestep) * 6) + &
         (kind(writer_to_serialise%previous_write_time) * 4) + &
         c_size(writer_to_serialise%pending_writes) * (kind(writer_to_serialise%write_timestep) + &
         kind(writer_to_serialise%previous_write_time))
    call check_thread_status(forthread_mutex_lock(writer_to_serialise%num_fields_to_write_mutex))

    if (size(writer_to_serialise%contents) .gt. 0) then
      do i=1, size(writer_to_serialise%contents)
        prepare_to_serialise_writer_type=prepare_to_serialise_writer_type+&
             prepare_to_serialise_writer_field_type(writer_to_serialise%contents(i))+kind(prepare_to_serialise_writer_type)
      end do
    end if
  end function prepare_to_serialise_writer_type

  !> Serialises a specific writer type into byte data (for storage or transmission.) Releases any locks issued during preparation
  !! @param writer_to_serialise The writer type to serialise
  !! @param byte_data The byte data which will be packed with the serialised data
  subroutine serialise_writer_type(writer_to_serialise, byte_data)
    type(writer_type), intent(inout) :: writer_to_serialise
    character, dimension(:), allocatable, intent(inout) :: byte_data

    integer :: prev_pt, i, current_data_point

    type(iterator_type) :: iterator
    class(*), pointer :: generic

    current_data_point=1

    current_data_point=pack_scalar_field(byte_data, current_data_point, writer_to_serialise%write_timestep)
    current_data_point=pack_scalar_field(byte_data, current_data_point, writer_to_serialise%num_fields_to_write)
    call check_thread_status(forthread_mutex_unlock(writer_to_serialise%num_fields_to_write_mutex))
    current_data_point=pack_scalar_field(byte_data, current_data_point, writer_to_serialise%previous_write_timestep)
    current_data_point=pack_scalar_field(byte_data, current_data_point, writer_to_serialise%latest_pending_write_timestep)
    current_data_point=pack_scalar_field(byte_data, current_data_point, single_real_value=writer_to_serialise%previous_write_time)
    current_data_point=pack_scalar_field(byte_data, current_data_point, &
         single_real_value=writer_to_serialise%latest_pending_write_time)
    current_data_point=pack_scalar_field(byte_data, current_data_point, single_real_value=writer_to_serialise%write_time)
    current_data_point=pack_scalar_field(byte_data, current_data_point, single_real_value=writer_to_serialise%defined_write_time)
    current_data_point=pack_scalar_field(byte_data, current_data_point, c_size(writer_to_serialise%pending_writes))
    current_data_point=pack_scalar_field(byte_data, current_data_point, size(writer_to_serialise%contents))

    if (c_size(writer_to_serialise%pending_writes) .gt. 0) then
      iterator=c_get_iterator(writer_to_serialise%pending_writes)
      do while (c_has_next(iterator))
        generic=>c_next_generic(iterator)
        select type(generic)
          type is (pending_write_type)
            current_data_point=pack_scalar_field(byte_data, current_data_point, generic%timestep)
            current_data_point=pack_scalar_field(byte_data, current_data_point, single_real_value=generic%write_time)
        end select        
      end do
    end if
    
    if (size(writer_to_serialise%contents) .gt. 0) then
      do i=1, size(writer_to_serialise%contents)
        prev_pt=current_data_point
        current_data_point=current_data_point+kind(current_data_point)
        call serialise_writer_field_type(writer_to_serialise%contents(i), byte_data, current_data_point)
        prev_pt=pack_scalar_field(byte_data, prev_pt, (current_data_point-prev_pt)-kind(current_data_point))
      end do
    end if   
  end subroutine serialise_writer_type
  
  !> Unserialises some byte data into the writer in order to recreate the state of the writer
  !! @param writer_to_unserialise The writer to unserialise and fill in
  !! @param byte_data The raw byte data to read from
  subroutine unserialise_writer_type(writer_to_unserialise, byte_data)
    type(writer_type), intent(inout) :: writer_to_unserialise
    character, dimension(:), intent(in) :: byte_data

    integer :: current_data_point, byte_size, expected_pending_writes, expected_contents, i
    type(pending_write_type), pointer :: pwt
    class(*), pointer :: generic

    current_data_point=1
    writer_to_unserialise%write_timestep=unpack_scalar_integer_from_bytedata(byte_data, current_data_point)
    writer_to_unserialise%num_fields_to_write=unpack_scalar_integer_from_bytedata(byte_data, current_data_point)
    writer_to_unserialise%previous_write_timestep=unpack_scalar_integer_from_bytedata(byte_data, current_data_point)
    writer_to_unserialise%latest_pending_write_timestep=unpack_scalar_integer_from_bytedata(byte_data, current_data_point)
    writer_to_unserialise%previous_write_time=unpack_scalar_real_from_bytedata(byte_data, current_data_point)
    writer_to_unserialise%latest_pending_write_time=unpack_scalar_real_from_bytedata(byte_data, current_data_point)
    writer_to_unserialise%write_time=unpack_scalar_real_from_bytedata(byte_data, current_data_point)
    writer_to_unserialise%defined_write_time=unpack_scalar_real_from_bytedata(byte_data, current_data_point)
    expected_pending_writes=unpack_scalar_integer_from_bytedata(byte_data, current_data_point)
    expected_contents=unpack_scalar_integer_from_bytedata(byte_data, current_data_point)

    if (expected_contents .ne. size(writer_to_unserialise%contents)) then
      call log_log(LOG_ERROR, "Expected number of writer entry fields in the checkpoint does not match the configured number")
    end if

    if (expected_pending_writes .gt. 0) then
      do i=1, expected_pending_writes
        allocate(pwt)        
        pwt%timestep=unpack_scalar_integer_from_bytedata(byte_data, current_data_point)
        pwt%write_time=unpack_scalar_real_from_bytedata(byte_data, current_data_point)
        generic=>pwt
        call c_push_generic(writer_to_unserialise%pending_writes, generic, .false.)
      end do      
    end if    

    if (expected_contents .gt. 0) then
      do i=1, expected_contents
        byte_size=unpack_scalar_integer_from_bytedata(byte_data, current_data_point)
        call unserialise_writer_field_type(writer_to_unserialise%contents(i), &
             byte_data(current_data_point:current_data_point+byte_size-1))
        current_data_point=current_data_point+byte_size
      end do
    end if
  end subroutine unserialise_writer_type

  !> Prepares to serialise a specific writer field, both determines the data size and issues any locks
  !! @param writer_field_to_serialise The writer field type to prepare for serialisation
  !! @returns The number of bytes needed to store the serialised state
  integer(kind=8) function prepare_to_serialise_writer_field_type(writer_field_to_serialise)
    type(writer_field_type), intent(inout) :: writer_field_to_serialise

    type(iterator_type) :: iterator
    class(*), pointer :: generic
    type(mapentry_type) :: map_entry

    call check_thread_status(forthread_mutex_lock(writer_field_to_serialise%values_mutex))
    prepare_to_serialise_writer_field_type=(kind(writer_field_to_serialise%latest_timestep_values) * 2) + &
         (kind(writer_field_to_serialise%previous_write_time) * 2)

    iterator=c_get_iterator(writer_field_to_serialise%values_to_write)
    do while (c_has_next(iterator))
      map_entry=c_next_mapentry(iterator)
      generic=>c_get_generic(map_entry)
      if (associated(generic)) then
        select type(generic)
        type is (data_values_type)
          prepare_to_serialise_writer_field_type=prepare_to_serialise_writer_field_type+&
               prepare_to_serialise_data_values_type(generic)+&
               (kind(prepare_to_serialise_writer_field_type)*3)+len(trim(map_entry%key))          
        type is (write_field_collective_values_type)
          prepare_to_serialise_writer_field_type=prepare_to_serialise_writer_field_type+&
               prepare_to_serialise_collective_values_type(generic)+&
               (kind(prepare_to_serialise_writer_field_type)*3)+len(trim(map_entry%key))          
        class default
          call log_log(LOG_ERROR, "Unknown data type in writer field type")
        end select        
      end if
    end do
  end function prepare_to_serialise_writer_field_type

  !> Serialises a specific writer field type for storage or transmission. This releases any locks issued during preparation
  !! @param writer_field_to_serialise The writer field type to serialise
  !! @param byte_data The byte data to pack with the serialised data
  !! @param current_data_point The current write point in the byte data, is updated during call so represents next point on return
  subroutine serialise_writer_field_type(writer_field_to_serialise, byte_data, current_data_point)
    type(writer_field_type), intent(inout) :: writer_field_to_serialise
    character, dimension(:), allocatable, intent(inout) :: byte_data
    integer, intent(inout) :: current_data_point

    integer :: prev_pt, byte_size, entry_type
    class(*), pointer :: generic
    type(mapentry_type) :: map_entry
    type(iterator_type) :: iterator
        
    current_data_point=pack_scalar_field(byte_data, current_data_point, writer_field_to_serialise%latest_timestep_values)
    current_data_point=pack_scalar_field(byte_data, current_data_point, &
         single_real_value=writer_field_to_serialise%previous_write_time)
    current_data_point=pack_scalar_field(byte_data, current_data_point, &
         single_real_value=writer_field_to_serialise%previous_tracked_write_point)
    
    current_data_point=pack_scalar_field(byte_data, current_data_point, c_size(writer_field_to_serialise%values_to_write))

    iterator=c_get_iterator(writer_field_to_serialise%values_to_write)
    do while (c_has_next(iterator))
      map_entry=c_next_mapentry(iterator)
      generic=>c_get_generic(map_entry)
      if (associated(generic)) then
        current_data_point=pack_scalar_field(byte_data, current_data_point, len(trim(map_entry%key)))
        byte_data(current_data_point:current_data_point+len(trim(map_entry%key))-1) = transfer(trim(map_entry%key), &
             byte_data(current_data_point:current_data_point+len(trim(map_entry%key))-1))          
        current_data_point=current_data_point+len(trim(map_entry%key))
        prev_pt=current_data_point
        current_data_point=current_data_point+(kind(current_data_point)*2)
        select type(generic)
        type is (data_values_type)
          call serialise_data_values_type(generic, byte_data, current_data_point)
          entry_type=1
        type is (write_field_collective_values_type)
          call serialise_collective_values_type(generic, byte_data, current_data_point)
          entry_type=2
        class default
          call log_log(LOG_ERROR, "Unknown data type in writer field type")
        end select        
        prev_pt=pack_scalar_field(byte_data, prev_pt, (current_data_point-(kind(current_data_point)*2)) - prev_pt)
        prev_pt=pack_scalar_field(byte_data, prev_pt, entry_type)        
      end if
    end do
    call check_thread_status(forthread_mutex_unlock(writer_field_to_serialise%values_mutex))
  end subroutine serialise_writer_field_type

  !> Unserialises byte data into a writer field type
  !! @param writer_field_to_unserialise The writer field to fill in
  !! @param byte_data The raw byte data to read from
  subroutine unserialise_writer_field_type(writer_field_to_unserialise, byte_data)
    type(writer_field_type), intent(inout) :: writer_field_to_unserialise
    character, dimension(:), intent(in) :: byte_data

    integer :: current_data_point, number_values_stored, i, byte_size, key_size, entry_type
    character(len=STRING_LENGTH) :: value_key
    class(*), pointer :: generic

    current_data_point=1
    writer_field_to_unserialise%latest_timestep_values=unpack_scalar_integer_from_bytedata(byte_data, current_data_point)    
    writer_field_to_unserialise%previous_write_time=unpack_scalar_real_from_bytedata(byte_data, current_data_point)
    writer_field_to_unserialise%previous_tracked_write_point=unpack_scalar_real_from_bytedata(byte_data, current_data_point)
    number_values_stored=unpack_scalar_integer_from_bytedata(byte_data, current_data_point)

    if (number_values_stored .gt. 0) then
      do i=1, number_values_stored
        key_size=unpack_scalar_integer_from_bytedata(byte_data, current_data_point)
        value_key=transfer(byte_data(current_data_point:current_data_point+key_size-1), value_key)
        value_key(key_size+1:)=" "
        current_data_point=current_data_point+key_size
        byte_size=unpack_scalar_integer_from_bytedata(byte_data, current_data_point)
        entry_type=unpack_scalar_integer_from_bytedata(byte_data, current_data_point)      
        if (entry_type == 1) then
          generic=>unserialise_data_values_type(byte_data(current_data_point:current_data_point+byte_size-1))
        else if (entry_type == 2) then
          generic=>unserialise_collective_values_type(byte_data(current_data_point:current_data_point+byte_size-1))
        else
          call log_log(LOG_ERROR, "Unknown entry type in writer field type serialisation bytes")
        end if
        call c_put_generic(writer_field_to_unserialise%values_to_write, value_key, generic, .false.)
        current_data_point=current_data_point+byte_size
      end do
    end if
  end subroutine unserialise_writer_field_type

  !> Prepares to serialise a specific collective value, both determines the required byte storate size and issues any locks
  !! @param collective_values_to_serialise The collective values to prepare for serialisation
  !! @returns The number of bytes needed to store the serialised state
  integer(kind=8) function prepare_to_serialise_collective_values_type(collective_values_to_serialise)
    type(write_field_collective_values_type), intent(inout) :: collective_values_to_serialise

    class(*), pointer :: generic
    type(mapentry_type) :: map_entry
    type(iterator_type) :: iterator

    prepare_to_serialise_collective_values_type=kind(prepare_to_serialise_collective_values_type)
    if (c_size(collective_values_to_serialise%monc_values) .gt. 0) then
      iterator=c_get_iterator(collective_values_to_serialise%monc_values)
      do while (c_has_next(iterator))
        map_entry=c_next_mapentry(iterator)
        generic=>c_get_generic(map_entry)
        if (associated(generic)) then
          select type(generic)
          type is (data_values_type)
            prepare_to_serialise_collective_values_type=prepare_to_serialise_collective_values_type+&
                 prepare_to_serialise_data_values_type(generic)+&
                 (kind(prepare_to_serialise_collective_values_type)*2)+len(trim(map_entry%key))         
          class default
            call log_log(LOG_ERROR, "Unknown data type in collective values type")
          end select
        end if
      end do
    end if
  end function prepare_to_serialise_collective_values_type

  !> Serialises collective values. This releases any locks issued during preparation
  !! @param collective_values_to_serialise The collective values to serialise
  !! @param byte_data The byte data which will be packed with the serialised byte code
  !! @param current_data_point The current write point in the byte data, is updated during call so represents next point on return
  subroutine serialise_collective_values_type(collective_values_to_serialise, byte_data, current_data_point)
    type(write_field_collective_values_type), intent(inout) :: collective_values_to_serialise
    character, dimension(:), allocatable, intent(inout) :: byte_data
    integer, intent(inout) :: current_data_point

    integer :: prev_pt
    class(*), pointer :: generic
    type(mapentry_type) :: map_entry
    type(iterator_type) :: iterator

    current_data_point=pack_scalar_field(byte_data, current_data_point, c_size(collective_values_to_serialise%monc_values))

    if (c_size(collective_values_to_serialise%monc_values) .gt. 0) then
      iterator=c_get_iterator(collective_values_to_serialise%monc_values)
      do while (c_has_next(iterator))
        map_entry=c_next_mapentry(iterator)
        generic=>c_get_generic(map_entry)
        if (associated(generic)) then
          select type(generic)
          type is (data_values_type)
            current_data_point=pack_scalar_field(byte_data, current_data_point, len(trim(map_entry%key)))
            byte_data(current_data_point:current_data_point+len(trim(map_entry%key))-1) = transfer(trim(map_entry%key), &
                 byte_data(current_data_point:current_data_point+len(trim(map_entry%key))-1))          
            current_data_point=current_data_point+len(trim(map_entry%key))

            prev_pt=current_data_point
            current_data_point=current_data_point+kind(current_data_point)
            call serialise_data_values_type(generic, byte_data, current_data_point)            
            
            prev_pt=pack_scalar_field(byte_data, prev_pt, (current_data_point-kind(current_data_point))-prev_pt)
            class default
            call log_log(LOG_ERROR, "Unknown data type in collective values type")
          end select
        end if
      end do
    end if
  end subroutine serialise_collective_values_type

  !> Unserialsies collective values contained in some data
  !! @param byte_data The byte data to unserialise
  !! @returns The collective values type which represents this byte data
  function unserialise_collective_values_type(byte_data)
    character, dimension(:), intent(in) :: byte_data
    type(write_field_collective_values_type), pointer :: unserialise_collective_values_type

    integer :: current_data_point, number_entries, i, key_size, byte_size
    character(len=STRING_LENGTH) :: value_key
    class(*), pointer :: generic

    allocate(unserialise_collective_values_type)
    
    current_data_point=1
    number_entries=unpack_scalar_integer_from_bytedata(byte_data, current_data_point)

    if (number_entries .gt. 0) then
      do i=1, number_entries
        key_size=unpack_scalar_integer_from_bytedata(byte_data, current_data_point)
        value_key=transfer(byte_data(current_data_point:current_data_point+key_size-1), value_key)
        value_key(key_size+1:)=" "
        current_data_point=current_data_point+key_size
        byte_size=unpack_scalar_integer_from_bytedata(byte_data, current_data_point)
        generic=>unserialise_data_values_type(byte_data(current_data_point:current_data_point+byte_size-1))
        call c_put_generic(unserialise_collective_values_type%monc_values, value_key, generic, .false.)
        current_data_point=current_data_point+byte_size
      end do
    end if
  end function unserialise_collective_values_type

  !> Prepares to serialise a specific data values type, both determines the byte size required and also issues any locks
  !! @param data_values_to_serialise The data values to prepare for serialisation
  !! @returns The number of bytes needed to store the serialised state
  integer(kind=8) function prepare_to_serialise_data_values_type(data_values_to_serialise)
    type(data_values_type), intent(inout) :: data_values_to_serialise

    integer :: values_size

    prepare_to_serialise_data_values_type=kind(prepare_to_serialise_data_values_type) * 8
    
    if (allocated(data_values_to_serialise%values)) then
      prepare_to_serialise_data_values_type=prepare_to_serialise_data_values_type+&
           (kind(data_values_to_serialise%values)*size(data_values_to_serialise%values))
    else if (allocated(data_values_to_serialise%string_values)) then
      prepare_to_serialise_data_values_type=prepare_to_serialise_data_values_type+&
           (size(data_values_to_serialise%string_values) * STRING_LENGTH)
    else
      prepare_to_serialise_data_values_type=prepare_to_serialise_data_values_type+&
           prepare_to_serialise_string_map(data_values_to_serialise%map_values)
    end if        
  end function prepare_to_serialise_data_values_type
  
  !> Serialises some data values to store or transmit. This releases any locks issued during preparation
  !! @param data_values_to_serialise The data values to serialise
  !! @param byte_data The byte data which will be packaged with the serialised byte code
  !! @param current_data_point The current write point in the byte data, is updated during call so represents next point on return
  subroutine serialise_data_values_type(data_values_to_serialise, byte_data, current_data_point)
    type(data_values_type), intent(inout) :: data_values_to_serialise
    character, dimension(:), allocatable, intent(inout) :: byte_data
    integer, intent(inout) :: current_data_point

    integer :: values_size, prev_pt
    character, dimension(:), allocatable :: dvt_byte_data, temp

    if (allocated(data_values_to_serialise%values)) then
      values_size=kind(data_values_to_serialise%values)*size(data_values_to_serialise%values)
    else if (allocated(data_values_to_serialise%string_values)) then
      values_size=size(data_values_to_serialise%string_values) * STRING_LENGTH
    else
      values_size=0
    end if
    
    current_data_point=pack_scalar_field(byte_data, current_data_point, data_values_to_serialise%data_type)
    current_data_point=pack_scalar_field(byte_data, current_data_point, data_values_to_serialise%dimensions)
    current_data_point=pack_array_field(byte_data, current_data_point, data_values_to_serialise%dim_sizes)
    if (allocated(data_values_to_serialise%values)) then
      current_data_point=pack_scalar_field(byte_data, current_data_point, 1)
      current_data_point=pack_scalar_field(byte_data, current_data_point, size(data_values_to_serialise%values))
      current_data_point=pack_array_field(byte_data, current_data_point, real_array_1d=data_values_to_serialise%values)
    else if (allocated(data_values_to_serialise%string_values)) then
      current_data_point=pack_scalar_field(byte_data, current_data_point, 2)
      current_data_point=pack_scalar_field(byte_data, current_data_point, &
           size(data_values_to_serialise%string_values) * STRING_LENGTH)
      byte_data(current_data_point:current_data_point+(size(data_values_to_serialise%string_values) * STRING_LENGTH)-1) = &
           transfer(data_values_to_serialise%string_values, &
           byte_data(current_data_point:current_data_point+(size(data_values_to_serialise%string_values) * STRING_LENGTH)-1))
    else
      current_data_point=pack_scalar_field(byte_data, current_data_point, 3)      
      call serialise_string_map(data_values_to_serialise%map_values, byte_data, current_data_point)      
    end if
  end subroutine serialise_data_values_type

  !> Unserialises some byte data into data values
  !! @param byte_data The byte data which is read from
  !! @returns The unserialised and allocated data values
  function unserialise_data_values_type(byte_data)
    character, dimension(:), intent(in) :: byte_data
    type(data_values_type), pointer :: unserialise_data_values_type

    integer :: current_data_point, i, values_size, byte_size, values_type

    allocate(unserialise_data_values_type)
    current_data_point=1
    unserialise_data_values_type%data_type=unpack_scalar_integer_from_bytedata(byte_data, current_data_point)
    unserialise_data_values_type%dimensions=unpack_scalar_integer_from_bytedata(byte_data, current_data_point)
    do i=1, 4
      unserialise_data_values_type%dim_sizes(i)=unpack_scalar_integer_from_bytedata(byte_data, current_data_point)
    end do
    values_type=unpack_scalar_integer_from_bytedata(byte_data, current_data_point)
    if (values_type == 1) then
      values_size=unpack_scalar_integer_from_bytedata(byte_data, current_data_point)
      allocate(unserialise_data_values_type%values(values_size))
      byte_size=values_size*kind(unserialise_data_values_type%values)
      unserialise_data_values_type%values=transfer(byte_data(current_data_point:current_data_point+byte_size-1), &
           unserialise_data_values_type%values)
    else if (values_type == 2) then
      values_size=unpack_scalar_integer_from_bytedata(byte_data, current_data_point)
      allocate(unserialise_data_values_type%string_values(values_size))
      byte_size=values_size*STRING_LENGTH
      unserialise_data_values_type%string_values=transfer(byte_data(current_data_point:current_data_point+byte_size-1), &
           unserialise_data_values_type%string_values)              
    else if (values_type == 3) then
      byte_size=unpack_scalar_integer_from_bytedata(byte_data, current_data_point)
      unserialise_data_values_type%map_values=unserialise_string_map(byte_data(current_data_point:current_data_point+byte_size-1))
    else
      call log_log(LOG_ERROR, "Unknown values type in data values serialisation bytes")
    end if      
  end function unserialise_data_values_type

  !> Prepares a map for serialisation, both determines the size of storage required and also issues any locks
  !! @param map_to_serialise The map which will be prepared for serialisation
  !! @returns The number of bytes needed to store the serialised state
  integer(kind=8) function prepare_to_serialise_string_map(map_to_serialise)
    type(map_type), intent(inout) :: map_to_serialise

    type(mapentry_type) :: map_entry
    type(iterator_type) :: iterator
    character(len=STRING_LENGTH) :: str_value

    prepare_to_serialise_string_map=kind(prepare_to_serialise_string_map)
    iterator=c_get_iterator(map_to_serialise)
    do while (c_has_next(iterator))
      map_entry=c_next_mapentry(iterator)
      str_value=c_get_string(map_entry)
      prepare_to_serialise_string_map=prepare_to_serialise_string_map+len(trim(map_entry%key)) + &
           len(trim(str_value)) + (kind(prepare_to_serialise_string_map) * 2)
    end do
  end function prepare_to_serialise_string_map  

  !> Serialises a string map, where the values are assumed to be strings. This releases any locks issued during preparation
  !! @param map_to_serialise The map which will be serialised
  !! @param byte_data The byte data representation of this map, this is allocated here
  !! @param current_data_point The current write point in the byte data, is updated during call so represents next point on return
  subroutine serialise_string_map(map_to_serialise, byte_data, current_data_point)
    type(map_type), intent(inout) :: map_to_serialise
    character, dimension(:), allocatable, intent(inout) :: byte_data
    integer, intent(inout) :: current_data_point

    type(mapentry_type) :: map_entry
    type(iterator_type) :: iterator
    character(len=STRING_LENGTH) :: str_value

    current_data_point=pack_scalar_field(byte_data, current_data_point, c_size(map_to_serialise))

    iterator=c_get_iterator(map_to_serialise)
    do while (c_has_next(iterator))
      map_entry=c_next_mapentry(iterator)
      str_value=c_get_string(map_entry)
            
      current_data_point=pack_scalar_field(byte_data, current_data_point, len(trim(map_entry%key)))
      byte_data(current_data_point:current_data_point+len(trim(map_entry%key))-1) = transfer(trim(map_entry%key), &
           byte_data(current_data_point:current_data_point+len(trim(map_entry%key))-1))          
      current_data_point=current_data_point+len(trim(map_entry%key))
      current_data_point=pack_scalar_field(byte_data, current_data_point, len(trim(str_value)))
      if (len(trim(str_value)) .gt. 0) then
        byte_data(current_data_point:current_data_point+len(trim(str_value))-1) = transfer(trim(str_value), &
             byte_data(current_data_point:current_data_point+len(trim(str_value))-1))
      end if
      current_data_point=current_data_point+len(trim(str_value))
    end do
  end subroutine serialise_string_map

  !> Inflates some byte data into a string map
  !! @param byte_data The byte data to unserialise
  !! @returns The inflated map representation
  function unserialise_string_map(byte_data)
    character, dimension(:), intent(in) :: byte_data
    type(map_type), pointer :: unserialise_string_map

    integer :: current_data_point, number_entries, i, key_size, value_size
    character(len=STRING_LENGTH) :: value_key, value_value

    allocate(unserialise_string_map)
    current_data_point=1
    number_entries=unpack_scalar_integer_from_bytedata(byte_data, current_data_point)

    if (number_entries .gt. 0) then
      do i=1, number_entries
        key_size=unpack_scalar_integer_from_bytedata(byte_data, current_data_point)
        value_key=transfer(byte_data(current_data_point:current_data_point+key_size-1), value_key)
        value_key(key_size+1:)=" "
        current_data_point=current_data_point+key_size
        value_size=unpack_scalar_integer_from_bytedata(byte_data, current_data_point)
        if (value_size .gt. 0) then
          value_value=transfer(byte_data(current_data_point:current_data_point+value_size-1), value_value)
        else
          value_value=""
        end if
        current_data_point=current_data_point+value_size
        call c_put_string(unserialise_string_map, value_key, value_value)
      end do
    end if
  end function unserialise_string_map  
end module writer_types_mod
