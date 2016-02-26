!> This federates over the writing of diagnostic and prognostic data to the file system. It also manages the time manipulation
!! of fields and groups.
module writer_federator_mod
  use datadefn_mod, only : DEFAULT_PRECISION, STRING_LENGTH
  use configuration_parser_mod, only : TIME_AVERAGED_TYPE, INSTANTANEOUS_TYPE, GROUP_TYPE, FIELD_TYPE, io_configuration_type, &
       io_configuration_field_type, io_configuration_diagnostic_field_type, io_configuration_data_definition_type, &
       data_values_type, get_data_value_by_field_name, get_diagnostic_field_configuration, get_prognostic_field_configuration, &
       get_monc_location
  use instantaneous_time_manipulation_mod, only : init_instantaneous_manipulation, finalise_instantaneous_manipulation, &
       perform_instantaneous_time_manipulation
  use timeaveraged_time_manipulation_mod, only : init_time_averaged_manipulation, finalise_time_averaged_manipulation, &
       perform_timeaveraged_time_manipulation
  use collections_mod, only : queue_type, list_type, map_type, hashmap_type, hashset_type, iterator_type, mapentry_type, &
       c_contains, c_size, c_get_string, c_get_generic, c_get_integer, c_add_string, c_free, c_put_real, c_put_generic, &
       c_key_at, c_is_empty, c_remove, c_push_generic, c_pop_generic, c_real_at, c_get_real, c_get_iterator, &
       c_has_next, c_next_mapentry, c_next_string, c_get_real, c_put_integer
  use conversions_mod, only : conv_to_string, conv_single_real_to_double, conv_to_integer, conv_to_real
  use io_server_client_mod, only : ARRAY_FIELD_TYPE
  use forthread_mod, only : forthread_mutex_init, forthread_mutex_lock, forthread_mutex_unlock, forthread_mutex_destroy, &
       forthread_rwlock_rdlock, forthread_rwlock_wrlock, forthread_rwlock_unlock, forthread_rwlock_init, forthread_rwlock_destroy
  use threadpool_mod, only : check_thread_status
  use logging_mod, only : LOG_DEBUG, LOG_ERROR, LOG_WARN, log_log, log_master_log, log_get_logging_level, log_is_master
  use writer_types_mod, only : writer_type, writer_field_type, write_field_collective_values_type, pending_write_type, &
       collective_q_field_representation_type
  use netcdf_filetype_writer_mod, only : initialise_netcdf_filetype, finalise_netcdf_filetype, define_netcdf_file, &
       write_variable, close_netcdf_file  
  use global_callback_inter_io_mod, only : perform_global_callback
  use data_utils_mod, only : get_scalar_integer_from_monc, get_scalar_real_from_monc, is_field_present 
  implicit none

#ifndef TEST_MODE
  private
#endif

  type(writer_type), volatile, dimension(:), allocatable :: writer_entries
  type(hashset_type), volatile :: used_field_names, q_field_names
  type(hashmap_type), volatile :: time_points, q_field_splits, collective_q_field_dims

  integer, volatile :: time_points_rwlock

  public initialise_writer_federator, finalise_writer_federator, provide_ordered_field_to_writer_federator, &
       check_writer_for_trigger, issue_actual_write, is_field_used_by_writer_federator, inform_writer_federator_fields_present, &
       inform_writer_federator_time_point, provide_q_field_names_to_writer_federator, is_field_split_on_q
contains

  !> Initialises the write federator and configures it based on the user configuration. Also initialises the time manipulations
  !! @param io_configuration The IO server configuration
  subroutine initialise_writer_federator(io_configuration, diagnostic_generation_frequency)
    type(io_configuration_type), intent(inout) :: io_configuration
    type(hashmap_type), intent(inout) :: diagnostic_generation_frequency

    integer :: i, j, number_contents, current_field_index
    type(hashset_type) :: writer_field_names, duplicate_field_names
    
    call check_thread_status(forthread_rwlock_init(time_points_rwlock, -1))

    call init_time_averaged_manipulation()
    call init_instantaneous_manipulation()
    call initialise_netcdf_filetype()
    
    allocate(writer_entries(io_configuration%number_of_writers))    
    do i=1, io_configuration%number_of_writers
      current_field_index=0
      number_contents=io_configuration%file_writers(i)%number_of_contents
      allocate(writer_entries(i)%contents(get_total_number_writer_fields(io_configuration, i)))
      writer_entries(i)%filename=io_configuration%file_writers(i)%file_name
      call check_thread_status(forthread_mutex_init(writer_entries(i)%trigger_and_write_mutex, -1))
      call check_thread_status(forthread_mutex_init(writer_entries(i)%num_fields_to_write_mutex, -1))
      call check_thread_status(forthread_mutex_init(writer_entries(i)%pending_writes_mutex, -1))
      writer_entries(i)%write_time_frequency=io_configuration%file_writers(i)%write_time_frequency
      writer_entries(i)%previous_write_time=0
      writer_entries(i)%defined_write_time=io_configuration%file_writers(i)%write_time_frequency
      writer_entries(i)%latest_pending_write_time=0
      writer_entries(i)%currently_writing=.false.
      do j=1, number_contents
        if (io_configuration%file_writers(i)%contents(j)%facet_type == GROUP_TYPE) then
          current_field_index=add_group_of_fields_to_writer_entry(io_configuration, i, j, current_field_index, &
               writer_field_names, duplicate_field_names, diagnostic_generation_frequency)
        else if (io_configuration%file_writers(i)%contents(j)%facet_type == FIELD_TYPE) then
          current_field_index=current_field_index+add_field_to_writer_entry(io_configuration, &
           i, j, current_field_index, io_configuration%file_writers(i)%contents(j)%facet_name, writer_field_names, &
           duplicate_field_names, diagnostic_generation_frequency)
        end if
      end do
      if (.not. c_is_empty(duplicate_field_names)) call handle_duplicate_field_names(writer_entries(i), duplicate_field_names)
      call c_free(writer_field_names)
      call c_free(duplicate_field_names)
    end do
  end subroutine initialise_writer_federator

  !> Finalises the write federator and the manipulations
  subroutine finalise_writer_federator()
    call check_thread_status(forthread_rwlock_destroy(time_points_rwlock))
    call finalise_time_averaged_manipulation()
    call finalise_instantaneous_manipulation()
    call finalise_netcdf_filetype()
  end subroutine finalise_writer_federator

  subroutine inform_writer_federator_time_point(io_configuration, source, data_id, data_dump)
    type(io_configuration_type), intent(inout) :: io_configuration
    integer, intent(in) :: source, data_id
    character, dimension(:), intent(in) :: data_dump

    real(kind=DEFAULT_PRECISION) :: time
    integer :: timestep
    character(len=STRING_LENGTH) :: timestep_key

    if (is_field_present(io_configuration, source, data_id, "time") .and. &
         is_field_present(io_configuration, source, data_id, "timestep")) then      
      time=get_scalar_real_from_monc(io_configuration, source, data_id, data_dump, "time")
      timestep=get_scalar_integer_from_monc(io_configuration, source, data_id, data_dump, "timestep")

      timestep_key=conv_to_string(timestep)

      call check_thread_status(forthread_rwlock_rdlock(time_points_rwlock))
      if (.not. c_contains(time_points, timestep_key)) then
        call check_thread_status(forthread_rwlock_unlock(time_points_rwlock))
        call check_thread_status(forthread_rwlock_wrlock(time_points_rwlock))
        if (.not. c_contains(time_points, timestep_key)) then
          call c_put_real(time_points, timestep_key, time)
        end if        
      end if
      call check_thread_status(forthread_rwlock_unlock(time_points_rwlock))      
    end if
  end subroutine inform_writer_federator_time_point

  !> Informs the writer federator that specific fields are present and should be reflected in the diagnostics output
  !! @param field_names The set of field names that are present
  subroutine inform_writer_federator_fields_present(io_configuration, field_names)
    type(io_configuration_type), intent(inout) :: io_configuration
    type(hashset_type), intent(inout) :: field_names

    type(iterator_type) :: iterator
    character(len=STRING_LENGTH) :: specific_name
    integer :: i, number_q_fields

    iterator=c_get_iterator(used_field_names)
    do while (c_has_next(iterator))
      specific_name=c_next_string(iterator)
      if (c_contains(field_names, specific_name)) then
        call enable_specific_field_by_name(specific_name)      
      end if
    end do
    iterator=c_get_iterator(q_field_names)
    do while (c_has_next(iterator))
      specific_name=c_next_string(iterator)
      if (c_contains(field_names, specific_name)) then
        number_q_fields=c_get_integer(io_configuration%dimension_sizing, "qfields")
        do i=1, number_q_fields
          if (c_size(io_configuration%q_field_names) .ge. i) then
            call enable_specific_field_by_name(trim(specific_name)//"_"//trim(c_get_string(io_configuration%q_field_names, i)))
          else
            call enable_specific_field_by_name(trim(specific_name)//"_udef"//trim(conv_to_string(i)))
          end if
        end do
      end if
    end do    
  end subroutine inform_writer_federator_fields_present  

  !> Determines whether a field is used by the writer federator or not
  !! @param field_name The field name to check whether it is being used or not
  !! @returns Whether this field is used or not
  logical function is_field_used_by_writer_federator(field_name)
    character(len=*), intent(in) :: field_name

    is_field_used_by_writer_federator=c_contains(used_field_names, field_name)
  end function is_field_used_by_writer_federator

  !> Determines whether a field is split on Q or not
  !! @param field_name The field name to check whether it is being used or not
  !! @returns Whether this field is used or not (and then further split to be constituient parts of Q)
  logical function is_field_split_on_q(field_name)
    character(len=*), intent(in) :: field_name

    is_field_split_on_q=c_contains(q_field_names, field_name)
  end function is_field_split_on_q  

  !> Enables a specific field by its name, this will locate all the fields with this name and enable them
  !! @param field_name The name of the field to enable
  subroutine enable_specific_field_by_name(field_name)
    character(len=*), intent(in) :: field_name

    logical :: continue_search
    integer :: writer_index, contents_index

    continue_search=.true.
    writer_index=1
    contents_index=0
    do while (continue_search)
      contents_index=contents_index+1
      continue_search=get_next_applicable_writer_entry(field_name, writer_index, contents_index)
      if (continue_search) then
        writer_entries(writer_index)%contents(contents_index)%enabled=.true.
      end if      
    end do    
  end subroutine enable_specific_field_by_name

  !> Provides the Q field names to the write federator, this is required as on initialisation we don't know what these are and
  !! only when MONC register do they inform the IO server of the specifics
  !! @param q_field_names An ordered list of Q field names
  subroutine provide_q_field_names_to_writer_federator(q_provided_field_names)
    type(list_type), intent(inout) :: q_provided_field_names

    type(iterator_type) :: iterator, q_field_iterator
    logical :: continue_search
    integer :: writer_index, contents_index, i
    character(len=STRING_LENGTH) :: search_field, field_name, specific_name

    iterator=c_get_iterator(q_field_names)
    do while (c_has_next(iterator))
      specific_name=c_next_string(iterator)
      q_field_iterator=c_get_iterator(q_provided_field_names)
      i=1
      do while (c_has_next(q_field_iterator))
        search_field=trim(specific_name)//"_udef"//trim(conv_to_string(i))
        field_name=trim(specific_name)//"_"//trim(c_next_string(q_field_iterator))        
        continue_search=.true.
        writer_index=1
        contents_index=0
        do while (continue_search)
          contents_index=contents_index+1
          continue_search=get_next_applicable_writer_entry(search_field, writer_index, contents_index)
          if (continue_search) then
            writer_entries(writer_index)%contents(contents_index)%field_name=field_name
          end if
        end do
        i=i+1
        call c_add_string(used_field_names, field_name)
        call c_remove(used_field_names, search_field)
      end do     
    end do    
  end subroutine provide_q_field_names_to_writer_federator

  !> Provides fields (either diagnostics or prognostics) to the write federator which will action these as appropriate. This will
  !! split Q fields up if appropriate
  !! @param io_configuration The IO server configuration
  !! @param field_name The field name to write (if appropriate)
  !! @param field_values The field values to write (if appropriate)
  !! @param timestep Corresponding MONC timestep
  !! @param time Corresponding MONC model time
  !! @param source Optional MONC source for the communicated fields
  subroutine provide_ordered_field_to_writer_federator(io_configuration, field_name, field_values, &
       timestep, time, source)
    type(io_configuration_type), intent(inout) :: io_configuration
    character(len=*), intent(in) :: field_name
    integer, intent(in) :: timestep, source
    real(kind=DEFAULT_PRECISION), dimension(:), intent(in) :: field_values
    real(kind=DEFAULT_PRECISION), intent(in) :: time

    type(iterator_type) :: iterator
    integer :: individual_size, index
    
    if (c_contains(used_field_names, field_name)) then
      call provide_ordered_single_field_to_writer_federator(io_configuration, field_name, field_values, &
           timestep, time, source)
    else if (c_contains(q_field_names, field_name)) then
      if (c_contains(q_field_splits, field_name)) then
        individual_size=c_get_integer(q_field_splits, field_name)
      else if (source .gt. -1) then
        individual_size=get_size_of_collective_q(io_configuration, field_name, source)
      else
        call log_log(LOG_WARN, "Can not find Q split field in Q field names or collective field names with source, ignoring")
        return
      end if
      iterator=c_get_iterator(io_configuration%q_field_names)
      index=1
      do while (c_has_next(iterator))
        call provide_ordered_single_field_to_writer_federator(io_configuration, &
             trim(field_name)//"_"//trim(c_next_string(iterator)), field_values(index:index+individual_size-1), &
             timestep, time, source)
        index=index+individual_size
      end do
    end if
  end subroutine provide_ordered_field_to_writer_federator

  !> Retrieves the data size for each Q entry of a collective Q field for the specific source MONC that has sent data
  !! @param io_configuration The IO server configuration
  !! @param field_name The field name to write (if appropriate)
  !! @param source  MONC source for the communicated fields
  !! @returns The size (elements) per Q split field
  integer function get_size_of_collective_q(io_configuration, field_name, source)
    type(io_configuration_type), intent(inout) :: io_configuration
    character(len=*), intent(in) :: field_name
    integer, intent(in) :: source

    class(*), pointer :: generic
    integer :: i, monc_index

    get_size_of_collective_q=1
    monc_index=get_monc_location(io_configuration, source)
    generic=>c_get_generic(collective_q_field_dims, field_name)
    select type(generic)
      type is(collective_q_field_representation_type)
        do i=1, size(generic%dimensions)
          get_size_of_collective_q=&
               get_size_of_collective_q*io_configuration%registered_moncs(monc_index)%local_dim_sizes(generic%dimensions(i))
        end do        
    end select    
  end function get_size_of_collective_q  

  !> Provides a single ordered field, i.e. Q fields have been split by this point
  !! @param io_configuration The IO server configuration
  !! @param field_name The field name to write (if appropriate)
  !! @param field_values The field values to write (if appropriate)
  !! @param timestep Corresponding MONC timestep
  !! @param time Corresponding MONC model time
  !! @param source Optional MONC source for the communicated fields
  subroutine provide_ordered_single_field_to_writer_federator(io_configuration, field_name, field_values, &
       timestep, time, source)
    type(io_configuration_type), intent(inout) :: io_configuration
    character(len=*), intent(in) :: field_name
    integer, intent(in) :: timestep, source
    real(kind=DEFAULT_PRECISION), dimension(:), intent(in) :: field_values
    real(kind=DEFAULT_PRECISION), intent(in) :: time

    integer :: writer_index, contents_index
    logical :: continue_search
    type(data_values_type), pointer :: result_values
    type(hashmap_type) :: typed_result_values
    class(*), pointer :: generic

    continue_search=.true.
    writer_index=1
    contents_index=0
    if (c_contains(used_field_names, field_name)) then
      do while (continue_search)
        contents_index=contents_index+1
        continue_search=get_next_applicable_writer_entry(field_name, writer_index, contents_index)
        if (continue_search) then
          if (.not. writer_entries(writer_index)%contents(contents_index)%enabled) then
            call log_log(LOG_WARN, "Received data for previously un-enabled field '"//&
                 writer_entries(writer_index)%contents(contents_index)%field_name//"'")
          end if          
          writer_entries(writer_index)%contents(contents_index)%enabled=.true.
          if (.not. c_contains(typed_result_values, conv_to_string(&
               writer_entries(writer_index)%contents(contents_index)%time_manipulation_type))) then
            allocate(result_values)
            if (writer_entries(writer_index)%contents(contents_index)%collective_write .and. source .gt. -1) then
              result_values=writer_entries(writer_index)%contents(contents_index)%time_manipulation(field_values, &
                   writer_entries(writer_index)%contents(contents_index)%output_frequency, &
                   trim(field_name)//"#"//conv_to_string(source), timestep, time)
            else
              result_values=writer_entries(writer_index)%contents(contents_index)%time_manipulation(field_values, &
                   writer_entries(writer_index)%contents(contents_index)%output_frequency, &
                   field_name, timestep, time)
            end if
            generic=>result_values
            call c_put_generic(typed_result_values, conv_to_string(&
                 writer_entries(writer_index)%contents(contents_index)%time_manipulation_type), generic, .false.)
          else
            result_values=>get_data_value_by_field_name(typed_result_values, conv_to_string(&
                 writer_entries(writer_index)%contents(contents_index)%time_manipulation_type))            
          end if
          if (allocated(result_values%values)) then
            if (log_get_logging_level() .ge. LOG_DEBUG) then
              call log_log(LOG_DEBUG, "[WRITE FED VALUE STORE] Storing value for field "//trim(field_name)//" ts="//&
                   trim(conv_to_string(timestep))// " t="//trim(conv_to_string(time)))
            end if
            call check_thread_status(forthread_mutex_lock(writer_entries(writer_index)%contents(contents_index)%values_mutex))
            if (writer_entries(writer_index)%contents(contents_index)%collective_write .and.  source .gt. -1) then
              call write_collective_write_value(result_values, writer_index, contents_index, source, conv_to_string(time))
            else
              call c_put_generic(writer_entries(writer_index)%contents(contents_index)%values_to_write, conv_to_string(time), &
                   generic, .false.)
            end if
            call check_thread_status(forthread_mutex_unlock(writer_entries(writer_index)%contents(contents_index)%values_mutex))
            if (writer_entries(writer_index)%contents(contents_index)%pending_to_write) then
              call determine_if_outstanding_field_can_be_written(io_configuration, writer_entries(writer_index), &
                   writer_entries(writer_index)%contents(contents_index))
            end if                        
          end if
        end if
      end do
    end if
    call c_free(typed_result_values)
  end subroutine provide_ordered_single_field_to_writer_federator

  !> Writes the collective values, this is held differently to independent values which are written directly - instead
  !! here we need to store the values for each MONC hence a specific type is used instead
  !! @param result_values The data values to store
  !! @param writer_index The file writer index
  !! @param contents_index The contents index
  !! @param source The MONC process id
  !! @param lookup_key The values lookup key
  subroutine write_collective_write_value(result_values, writer_index, contents_index, source, lookup_key)
    integer, intent(in) :: writer_index, contents_index, source
    type(data_values_type), pointer :: result_values
    character(len=*), intent(in) :: lookup_key

    class(*), pointer :: generic
    type(write_field_collective_values_type), pointer :: stored_monc_values

    if (c_contains(writer_entries(writer_index)%contents(contents_index)%values_to_write, lookup_key)) then
      generic=>c_get_generic(writer_entries(writer_index)%contents(contents_index)%values_to_write, lookup_key)
      select type(generic)
        type is(write_field_collective_values_type)
          stored_monc_values=>generic
      end select      
    else
      allocate(stored_monc_values)
      generic=>stored_monc_values
      call c_put_generic(writer_entries(writer_index)%contents(contents_index)%values_to_write, lookup_key, generic, .false.)
    end if
    generic=>result_values
    call c_put_generic(stored_monc_values%monc_values, conv_to_string(source), generic, .false.)
  end subroutine write_collective_write_value  

  !> For a specific field wil determine and handle any outstanding fields writes until an outstanding write
  !! can not be performed or the outstanding list is empty
  !! @param specific_field The specific field that we are concerned with
  subroutine determine_if_outstanding_field_can_be_written(io_configuration, writer_entry, specific_field)
    type(io_configuration_type), intent(inout) :: io_configuration
    type(writer_type), intent(inout) :: writer_entry
    type(writer_field_type), intent(inout) :: specific_field

    logical :: field_write_success, do_close_num_fields

    if (specific_field%pending_to_write) then    
      call determine_if_field_can_be_written(io_configuration, writer_entry, specific_field, writer_entry%write_timestep, &
           writer_entry%write_time, writer_entry%previous_write_time, field_write_success)
      if (field_write_success) then
        if (log_get_logging_level() .ge. LOG_DEBUG) then
          call log_log(LOG_DEBUG, "Flushed outstanding field ts="//conv_to_string(writer_entry%write_timestep)//&
               " write time="//conv_to_string(writer_entry%write_time))
        end if
        call check_thread_status(forthread_mutex_lock(writer_entry%num_fields_to_write_mutex))
        writer_entry%num_fields_to_write=writer_entry%num_fields_to_write-1
        do_close_num_fields=writer_entry%num_fields_to_write == 0
        call check_thread_status(forthread_mutex_unlock(writer_entry%num_fields_to_write_mutex))
        if (do_close_num_fields) then
          call close_diagnostics_file(io_configuration, writer_entry, writer_entry%write_timestep, writer_entry%write_time)
        end if    
      end if
    end if
  end subroutine determine_if_outstanding_field_can_be_written  

  !> Determines if a file can be written to its overarching write representation. If so then a write is issued, otherwise
  !! an outstanding write point is registered which will be checked frequency to do a write later on
  !! @param specific_field The specific field we are checking and going to write if possible
  !! @param timestep The current timestep that we are at for this write
  !! @param write_time The current time that we are at for this write
  !! @param field_written An optional output logical representing whether a write was performed or not
  subroutine determine_if_field_can_be_written(io_configuration, writer_entry, specific_field, &
       timestep, write_time, previous_write_time, field_written)
    type(io_configuration_type), intent(inout) :: io_configuration
    type(writer_type), intent(inout) :: writer_entry
    type(writer_field_type), intent(inout) :: specific_field
    integer, intent(in) :: timestep
    real, intent(in) :: write_time, previous_write_time
    logical, intent(out), optional :: field_written

    real :: value_to_test, largest_value_found
    integer :: num_matching
    type(iterator_type) :: iterator
    type(mapentry_type) :: map_entry
        
    num_matching=0
    largest_value_found=0.0
    call check_thread_status(forthread_mutex_lock(specific_field%values_mutex))
    if (.not. c_is_empty(specific_field%values_to_write)) then
      iterator=c_get_iterator(specific_field%values_to_write)      
      do while (c_has_next(iterator))
        map_entry=c_next_mapentry(iterator)
        value_to_test=conv_to_real(map_entry%key)
        if (value_to_test .le. write_time .and. value_to_test .gt. previous_write_time) then        
          num_matching=num_matching+1
          if (largest_value_found .lt. value_to_test) largest_value_found=value_to_test
        end if
      end do
    end if
    if (largest_value_found + specific_field%output_frequency .gt. write_time)  then
      if (num_matching .gt. 0) call write_variable(io_configuration, specific_field, writer_entry%filename, timestep, write_time)
      specific_field%previous_write_time=write_time
      specific_field%pending_to_write=.false.
      if (present(field_written)) field_written=.true.
    else
      if (log_get_logging_level() .ge. LOG_DEBUG) then
          call log_log(LOG_DEBUG, "Setting outstanding field ts="//conv_to_string(writer_entry%write_timestep)//&
               " write time="//conv_to_string(writer_entry%write_time)//" prev="//conv_to_string(previous_write_time)//&
               " largest entry="//conv_to_string(largest_value_found)//" num matching="//conv_to_string(num_matching))
        end if
      specific_field%pending_to_write=.true.
      if (present(field_written)) field_written=.false.
    end if
    call check_thread_status(forthread_mutex_unlock(specific_field%values_mutex))
  end subroutine determine_if_field_can_be_written

  !> Checks all writer entries for any trigger fires and issues the underlying file storage
  !! @param io_configuration Configuration of the IO server
  !! @param source The source PID of the MONC process
  !! @param data_id The ID of the data definition that is represented by the dump
  !! @param data_dump The data we have received from MONC
  subroutine check_writer_for_trigger(io_configuration, source, data_id, data_dump)
    type(io_configuration_type), intent(inout) :: io_configuration
    integer, intent(in) :: source, data_id
    character, dimension(:), intent(in) :: data_dump

    integer :: i, timestep
    real(kind=DEFAULT_PRECISION) :: time

    if (is_field_present(io_configuration, source, data_id, "timestep") .and. &
         is_field_present(io_configuration, source, data_id, "time")) then
      timestep=get_scalar_integer_from_monc(io_configuration, source, data_id, data_dump, "timestep")
      time=get_scalar_real_from_monc(io_configuration, source, data_id, data_dump, "time")

      do i=1, size(writer_entries)
        call check_writer_trigger(io_configuration, i, timestep, real(time, kind=4))
      end do
    end if
  end subroutine check_writer_for_trigger

  !> Checks a writer trigger and issues a file creation along with field write if the conditions (time or timestep) are met.
  !! This will either create and write to the file or store a pending state if one is already open (required due to 
  !! NetCDF/HDF5 limitations with thread safety and parallel access.)
  !! @param io_configuration The IO server configuration
  !! @param writer_entry_index Index of the writer we are concerned with
  !! @param timestep The corresponding timestep
  !! @param time The corresponding model time
  subroutine check_writer_trigger(io_configuration, writer_entry_index, timestep, time)
    type(io_configuration_type), intent(inout) :: io_configuration
    integer, intent(in) :: writer_entry_index, timestep
    real, intent(in) :: time

    real :: time_difference
    
    call check_thread_status(forthread_mutex_lock(writer_entries(writer_entry_index)%trigger_and_write_mutex))
    time_difference=time-writer_entries(writer_entry_index)%latest_pending_write_time
    if (time_difference .ge. writer_entries(writer_entry_index)%write_time_frequency) then
      writer_entries(writer_entry_index)%latest_pending_write_time=time      
      if (writer_entries(writer_entry_index)%currently_writing) then
        call check_thread_status(forthread_mutex_unlock(writer_entries(writer_entry_index)%trigger_and_write_mutex))
        call register_pending_file_write(writer_entry_index, timestep, time)
      else
        writer_entries(writer_entry_index)%currently_writing=.true.
        call check_thread_status(forthread_mutex_unlock(writer_entries(writer_entry_index)%trigger_and_write_mutex))
        call issue_actual_write(io_configuration, writer_entries(writer_entry_index), timestep, time)
      end if
    else
      call check_thread_status(forthread_mutex_unlock(writer_entries(writer_entry_index)%trigger_and_write_mutex))
    end if    
  end subroutine check_writer_trigger

  !> Issues the actual file creation, write of available fields and closure if all completed. 
  !! @param io_configuration The IO server configuration
  !! @param writer_entry_index Index of the writer we are concerned with
  !! @param timestep The corresponding timestep
  !! @param time The corresponding model time
  subroutine issue_actual_write(io_configuration, writer_entry, timestep, time)
    type(io_configuration_type), intent(inout) :: io_configuration
    type(writer_type), intent(inout) :: writer_entry
    integer, intent(in) :: timestep
    real, intent(in) :: time

    integer :: j, total_outstanding, num_written, total_flds
    logical :: empty_contents_here, field_written
    type(map_type) :: applicable_time_points
    
    writer_entry%write_time=time
    writer_entry%write_timestep=timestep
    applicable_time_points=extract_applicable_time_points(writer_entry%previous_write_time, time)
    call define_netcdf_file(io_configuration, writer_entry, timestep, time, applicable_time_points)
    call c_free(applicable_time_points)
    empty_contents_here=.true.
    total_outstanding=0
    total_flds=0
    num_written=0
    call check_thread_status(forthread_mutex_lock(writer_entry%num_fields_to_write_mutex))
    do j=1, size(writer_entry%contents)      
      if (.not. c_is_empty(writer_entry%contents(j)%values_to_write)) then
        empty_contents_here=.false.
        total_flds=total_flds+1
        call determine_if_field_can_be_written(io_configuration, writer_entry, writer_entry%contents(j), timestep, time, &
             writer_entry%contents(j)%previous_write_time, field_written)
        if (.not. field_written) then
          total_outstanding=total_outstanding+1
        else
          num_written=num_written+1
        end if
      end if
    end do
    writer_entry%num_fields_to_write=total_outstanding
    call check_thread_status(forthread_mutex_unlock(writer_entry%num_fields_to_write_mutex))
    if (log_get_logging_level() .ge. LOG_DEBUG) then
      call log_log(LOG_DEBUG, "Started write for NetCDF file, timestep= "//trim(conv_to_string(timestep))&
           //" total="//trim(conv_to_string(total_flds))//" written="//trim(conv_to_string(num_written))//&
           " outstanding="//trim(conv_to_string(total_outstanding)))
    end if
    if (empty_contents_here .or. total_outstanding == 0) then
      call close_diagnostics_file(io_configuration, writer_entry, timestep, time)
    end if    
  end subroutine issue_actual_write

  !> Extracts the applicable time points from the overall map that lie within a specific range
  !! @param start_time The start time where values must be greater than
  !! @param end_time The end time where values must be equal or less than
  !! @returns The sorted (based on timestep) list of MONC communication time points that correspond to this
  type(map_type) function extract_applicable_time_points(start_time, end_time)
    real, intent(in) :: start_time, end_time

    real :: time_entry
    type(iterator_type) :: iterator
    type(mapentry_type) :: map_entry

    call check_thread_status(forthread_rwlock_rdlock(time_points_rwlock))
    iterator=c_get_iterator(time_points)
    do while (c_has_next(iterator))
      map_entry=c_next_mapentry(iterator)
      time_entry=real(c_get_real(map_entry))
      if (time_entry .gt. start_time .and. time_entry .le. end_time) then
        call c_put_real(extract_applicable_time_points, map_entry%key, conv_single_real_to_double(time_entry))
      end if
    end do
    call check_thread_status(forthread_rwlock_unlock(time_points_rwlock))
    extract_applicable_time_points=sort_applicable_time_points(extract_applicable_time_points)
  end function extract_applicable_time_points

  !> Sorts the time points based upon their timestep, smallest to largest. Note that this is a bubble sort and as such
  !! inefficient, so would be good to change to something else but works OK for now
  !! @param unsorted_timepoints The unsorted timepoints, which is changed by destroying the map
  !! @returns The sorted version of the input map based upon timestep
  type(map_type) function sort_applicable_time_points(unsorted_timepoints)
    type(map_type), intent(inout) :: unsorted_timepoints

    integer :: i, entries, specific_ts, smallest_ts
    character(len=STRING_LENGTH) :: smallest_key
    real(kind=DEFAULT_PRECISION) :: rvalue
    type(iterator_type) :: iterator
    type(mapentry_type) :: map_entry

    entries=c_size(unsorted_timepoints)
    do i=1, entries
      smallest_key=""
      iterator=c_get_iterator(unsorted_timepoints)
      do while (c_has_next(iterator))
        map_entry=c_next_mapentry(iterator)
        specific_ts=conv_to_integer(map_entry%key)
        if (len_trim(smallest_key) == 0 .or. smallest_ts .gt. specific_ts) then
          smallest_ts=specific_ts
          smallest_key=map_entry%key
          rvalue=c_get_real(map_entry)
        end if        
      end do
      call c_put_real(sort_applicable_time_points, smallest_key, rvalue)
      call c_remove(unsorted_timepoints, smallest_key)
    end do
    call c_free(unsorted_timepoints)
  end function sort_applicable_time_points  

  !> Closes the diagnostics file, this is done via a global callback to issue the closes synchronously (collective
  !! operation)
  !! @param io_configuration The IO server configuration
  !! @param writer_entry The writer entry
  !! @param timestep The file write timestep
  !! @param time The file write time
  subroutine close_diagnostics_file(io_configuration, writer_entry, timestep, time)
    type(io_configuration_type), intent(inout) :: io_configuration
    type(writer_type), intent(inout) :: writer_entry
    integer, intent(in) :: timestep
    real, intent(in) :: time

    if (log_get_logging_level() .ge. LOG_DEBUG) then
      call log_log(LOG_DEBUG, "Issue close for NetCDF file at timestep "//trim(conv_to_string(timestep)))
    end if
    call perform_global_callback(io_configuration, writer_entry%filename, timestep, handle_close_diagnostics_globalcallback)
  end subroutine close_diagnostics_file

  !> Call back for the inter IO reduction which actually does the NetCDF file closing which is a 
  !! collective (synchronous) operation. Calls out to the NetCDF code to do the call and then checks the list of
  !! pending file writes to process any others that are waiting in the queue
  !! @param io_configuration The IO server configuration
  !! @param values The inter IO resulting values, we don't care about these
  !! @param field_name The field name that is being communicated
  !! @param timestep The write timestep
  subroutine handle_close_diagnostics_globalcallback(io_configuration, values, field_name, timestep)
    type(io_configuration_type), intent(inout) :: io_configuration
    real(DEFAULT_PRECISION), dimension(:) :: values
    character(len=STRING_LENGTH) :: field_name
    integer :: timestep

    class(*), pointer :: generic
    type(writer_type), pointer :: writer_entry

    writer_entry=>close_netcdf_file(io_configuration, field_name, timestep)

    writer_entry%previous_write_time=writer_entry%write_time   
    writer_entry%defined_write_time=writer_entry%defined_write_time+writer_entry%write_time_frequency
    call check_thread_status(forthread_mutex_lock(writer_entry%pending_writes_mutex))
    if (.not. c_is_empty(writer_entry%pending_writes)) then      
      generic=>c_pop_generic(writer_entry%pending_writes)    
      call check_thread_status(forthread_mutex_unlock(writer_entry%pending_writes_mutex))
      select type(generic)
      type is (pending_write_type)
        if (log_get_logging_level() .ge. LOG_DEBUG) then
          call log_log(LOG_DEBUG, "Chain to next pending entry ts= "//trim(conv_to_string(generic%timestep)))
        end if
        call issue_actual_write(io_configuration, writer_entry, generic%timestep, generic%write_time)
        deallocate(generic)
      end select
    else
      call check_thread_status(forthread_mutex_unlock(writer_entry%pending_writes_mutex))
      writer_entry%currently_writing=.false.
      if (log_get_logging_level() .ge. LOG_DEBUG) then
        call log_log(LOG_DEBUG, "No more pending entries to chain to at ts= "//trim(conv_to_string(timestep)))
      end if
    end if
  end subroutine handle_close_diagnostics_globalcallback
  
  !> Registers a pending file write which will be actioned later on
  !! @param writer_entry_index Index of the writer entry
  !! @param timestep The timestep that the pending write represents
  !! @param time The time of the pending write
  subroutine register_pending_file_write(writer_entry_index, timestep, time)
    integer, intent(in) :: writer_entry_index, timestep
    real, intent(in) :: time

    type(pending_write_type), pointer :: pending_write
    class(*), pointer :: generic

    allocate(pending_write)
    pending_write%write_time=time
    pending_write%timestep=timestep

    generic=>pending_write
    call check_thread_status(forthread_mutex_lock(writer_entries(writer_entry_index)%pending_writes_mutex))
    call c_push_generic(writer_entries(writer_entry_index)%pending_writes, generic, .false.)
    call check_thread_status(forthread_mutex_unlock(writer_entries(writer_entry_index)%pending_writes_mutex))
  end subroutine register_pending_file_write  

  !> Retrieves the index of the next writer which uses a specific field. If none is found then returns false, otherwise true
  !! @param field_name The field name to search for
  !! @param writer_start_point The index that we start searching from
  !! @param result_writer The index of the (next) writer that requires this field is written in here
  !! @param result_contents The index of the contents in the next writer that corresponds to this field
  !! @returns Whether or not a next entry has been found
  logical function get_next_applicable_writer_entry(field_name, writer_index_point, contents_index_point)
    character(len=*), intent(in) :: field_name
    integer, intent(inout) :: writer_index_point, contents_index_point

    integer :: i, j

    if (writer_index_point .le. size(writer_entries)) then    
      do i=writer_index_point, size(writer_entries)
        if (contents_index_point .le. size(writer_entries(i)%contents)) then
          do j=contents_index_point, size(writer_entries(i)%contents)           
            if (writer_entries(i)%contents(j)%field_name==field_name) then
              writer_index_point=i
              contents_index_point=j
              get_next_applicable_writer_entry=.true.
              return
            end if
          end do
        end if
        contents_index_point=1
      end do
    end if
    get_next_applicable_writer_entry=.false.
  end function get_next_applicable_writer_entry  

  !> Determines the total number of fields that make up a writer entry, this is all the fields of the groups that make up
  !! this writer and individual fields specified too
  !! @param io_configuration The IO server configuration
  !! @param writer_entry_index The index of the writer that we are inquiring about
  !! @returns The total number of fields that make up this writer
  integer function get_total_number_writer_fields(io_configuration, writer_entry_index)
    type(io_configuration_type), intent(inout) :: io_configuration
    integer, intent(in) :: writer_entry_index

    integer :: i, number_contents, group_index, number_q_fields

    get_total_number_writer_fields=0
    number_q_fields=c_get_integer(io_configuration%dimension_sizing, "qfields")

    number_contents=io_configuration%file_writers(writer_entry_index)%number_of_contents
    do i=1, number_contents
      if (io_configuration%file_writers(writer_entry_index)%contents(i)%facet_type == GROUP_TYPE) then          
        group_index=get_index_of_group(io_configuration, &
             io_configuration%file_writers(writer_entry_index)%contents(i)%facet_name)
        if (group_index == 0) call log_log(LOG_ERROR, "Can not find group '"//trim(&
             io_configuration%file_writers(writer_entry_index)%contents(i)%facet_name)//"'")
        get_total_number_writer_fields=get_total_number_writer_fields+&
             get_group_number_of_fields(io_configuration, io_configuration%groups(group_index)%members, number_q_fields)
      else if (io_configuration%file_writers(writer_entry_index)%contents(i)%facet_type == FIELD_TYPE) then
        get_total_number_writer_fields=get_total_number_writer_fields+get_field_number_of_fields(io_configuration, &
             io_configuration%file_writers(writer_entry_index)%contents(i)%facet_name, number_q_fields)
      end if
    end do
  end function get_total_number_writer_fields

  !> Retrieves the number of fields within a group of fields
  !! @param io_configuration The IO server configuration
  !! @param group_members The members of the group
  !! @param num_q_fields The number of Q fields
  !! @returns The number of fields that make up this group
  integer function get_group_number_of_fields(io_configuration, group_members, num_q_fields)
    type(io_configuration_type), intent(inout) :: io_configuration
    type(list_type) :: group_members
    integer, intent(in) :: num_q_fields

    type(iterator_type) :: iterator
    character(len=STRING_LENGTH) :: field_name
    
    get_group_number_of_fields=0
    iterator=c_get_iterator(group_members)
    do while (c_has_next(iterator))
      field_name=c_next_string(iterator)
      get_group_number_of_fields=get_group_number_of_fields+get_field_number_of_fields(io_configuration, field_name, num_q_fields)
    end do
  end function get_group_number_of_fields

  !> Retrieves the number of fields that make up this field, if it is a Q field then it will be split into many subfields
  !! hence it is not a simple 1-1 mapping
  !! @param io_configuration The IO server configuration
  !! @param field_name The name of the field
  !! @param num_q_fields The number of Q fields
  !! @returns The number of fields that make up this field
  integer function get_field_number_of_fields(io_configuration, field_name, num_q_fields)
    type(io_configuration_type), intent(inout) :: io_configuration
    character(len=STRING_LENGTH), intent(in) :: field_name
    integer, intent(in) :: num_q_fields

    type(io_configuration_field_type) :: prognostic_field_configuration
    type(io_configuration_data_definition_type) :: prognostic_containing_data_defn
    type(io_configuration_diagnostic_field_type) :: diagnostic_field_configuration
    
    if (get_diagnostic_field_configuration(io_configuration, field_name, diagnostic_field_configuration)) then
      if (diagnostic_field_configuration%field_type == ARRAY_FIELD_TYPE) then
        if (diagnostic_field_configuration%dim_size_defns(diagnostic_field_configuration%dimensions) .eq. "qfields") then
          get_field_number_of_fields=num_q_fields
          return
        end if
      end if
      get_field_number_of_fields=1
    else if (get_prognostic_field_configuration(io_configuration, field_name, &
         prognostic_field_configuration, prognostic_containing_data_defn)) then
      if (prognostic_field_configuration%field_type == ARRAY_FIELD_TYPE) then
        if (prognostic_field_configuration%dim_size_defns(prognostic_field_configuration%dimensions) .eq. "qfields") then
          get_field_number_of_fields=num_q_fields
          return
        end if
      end if
        get_field_number_of_fields=1
    end if
  end function get_field_number_of_fields

  !> Adds a group of fields to a writer entry, groups are expanded out into individual fields, each inherit the properties
  !! of the group
  !! @param io_configuration The IO server configuration
  !! @param writer_entry_index Index of the writer entry that we are dealing with
  !! @param facet_index Index of the facet (group) in the IO server configuration
  !! @param current_field_index The current field index in this internal module representation of the structure
  !! @returns The next field index to write to
  integer function add_group_of_fields_to_writer_entry(io_configuration, writer_entry_index, facet_index, current_field_index, &
       writer_field_names, duplicate_field_names, diagnostic_generation_frequency)
    type(io_configuration_type), intent(inout) :: io_configuration
    integer, intent(in) :: writer_entry_index, facet_index, current_field_index
    type(hashset_type), intent(inout) :: writer_field_names, duplicate_field_names
    type(hashmap_type), intent(inout) :: diagnostic_generation_frequency

    integer :: group_index
    character(len=STRING_LENGTH) :: field_name
    type(iterator_type) :: iterator

    add_group_of_fields_to_writer_entry=current_field_index
    group_index=get_index_of_group(io_configuration, &
         io_configuration%file_writers(writer_entry_index)%contents(facet_index)%facet_name)
    if (group_index == 0) then
      call log_log(LOG_ERROR, "Can not find group '"//&
           trim(io_configuration%file_writers(writer_entry_index)%contents(facet_index)%facet_name)//"' in the configuration")
    end if
    iterator=c_get_iterator(io_configuration%groups(group_index)%members)
    do while (c_has_next(iterator))
      field_name=c_next_string(iterator)
      add_group_of_fields_to_writer_entry=add_group_of_fields_to_writer_entry+add_field_to_writer_entry(io_configuration, &
           writer_entry_index, facet_index, add_group_of_fields_to_writer_entry, field_name, writer_field_names, &
           duplicate_field_names, diagnostic_generation_frequency)
    end do    
  end function add_group_of_fields_to_writer_entry

  !> Adds a field to the writer entry, this will split the Q fields. However at initialisation we don't know what the Q
  !! fields are called, hence place a marker which will be replaced later on
  !! @param io_configuration The IO server configuration
  !! @param writer_entry_index Index of the writer entry that we are dealing with
  !! @param io_config_facet_index Index of the facet (group) in the IO server configuration
  !! @param my_facet_index The current field index in this internal module representation of the structure
  !! @param field_name The name of the field that we are constructing
  !! @param writer_field_names The field names in the writer (for duplication checking)
  !! @param duplicate_field_names Duplicate field names in the wrier, for duplication checking
  !! @param diagnostic_generation_frequency Generation frequency of the diagnostics
  !! @returns Location for next field to be written to
  integer function add_field_to_writer_entry(io_configuration, writer_entry_index, io_config_facet_index, &
       my_facet_index, field_name, writer_field_names, duplicate_field_names, diagnostic_generation_frequency)
    type(io_configuration_type), intent(inout) :: io_configuration
    integer, intent(in) :: writer_entry_index, io_config_facet_index, my_facet_index
    character(len=*), intent(in) :: field_name
    type(hashset_type), intent(inout) :: writer_field_names, duplicate_field_names
    type(hashmap_type), intent(inout) :: diagnostic_generation_frequency

    integer :: i, number_q_fields, tot_size
    type(io_configuration_field_type) :: prognostic_field_configuration
    type(io_configuration_data_definition_type) :: prognostic_containing_data_defn
    type(io_configuration_diagnostic_field_type) :: diagnostic_field_configuration
    type(collective_q_field_representation_type), pointer :: collective_q_field
    class(*), pointer :: generic
    
    if (get_diagnostic_field_configuration(io_configuration, field_name, diagnostic_field_configuration)) then
      if (diagnostic_field_configuration%field_type == ARRAY_FIELD_TYPE) then
        if (diagnostic_field_configuration%dim_size_defns(diagnostic_field_configuration%dimensions) .eq. "qfields") then  
          number_q_fields=c_get_integer(io_configuration%dimension_sizing, "qfields")
          do i=1, number_q_fields
            call add_specific_field_to_writer_entry(io_configuration, writer_entry_index, io_config_facet_index, &
                 my_facet_index+i, trim(field_name)//"_udef"//trim(conv_to_string(i)), writer_field_names, &
                 duplicate_field_names, c_get_integer(diagnostic_generation_frequency, field_name), diagnostic_field_configuration)
          end do
          tot_size=1
          do i=1, writer_entries(writer_entry_index)%contents(my_facet_index+number_q_fields)%dimensions
            tot_size=tot_size*writer_entries(writer_entry_index)%contents(my_facet_index+number_q_fields)%actual_dim_size(i)
          end do
          call c_put_integer(q_field_splits, field_name, tot_size)
          add_field_to_writer_entry=number_q_fields
          return
        end if
      end if
      call add_specific_field_to_writer_entry(io_configuration, writer_entry_index, io_config_facet_index, &
           my_facet_index+1, field_name, writer_field_names, duplicate_field_names, &
           c_get_integer(diagnostic_generation_frequency, field_name), diagnostic_field_configuration)
      add_field_to_writer_entry=1      
    else if (get_prognostic_field_configuration(io_configuration, field_name, &
         prognostic_field_configuration, prognostic_containing_data_defn)) then
      if (prognostic_field_configuration%field_type == ARRAY_FIELD_TYPE) then
        if (prognostic_field_configuration%dim_size_defns(prognostic_field_configuration%dimensions) .eq. "qfields") then  
          number_q_fields=c_get_integer(io_configuration%dimension_sizing, "qfields")
          do i=1, number_q_fields
            call add_specific_field_to_writer_entry(io_configuration, writer_entry_index, io_config_facet_index, &
                 my_facet_index+i, trim(field_name)//"_udef"//trim(conv_to_string(i)), writer_field_names, &
                 duplicate_field_names, prognostic_containing_data_defn%frequency, &
                 prognostic_field_configuration=prognostic_field_configuration)
          end do
          if (prognostic_field_configuration%collective) then
            allocate(collective_q_field)
            allocate(collective_q_field%dimensions(&
                 writer_entries(writer_entry_index)%contents(my_facet_index+number_q_fields)%dimensions))
            do i=1, writer_entries(writer_entry_index)%contents(my_facet_index+number_q_fields)%dimensions
              if (trim(writer_entries(writer_entry_index)%contents(my_facet_index+number_q_fields)%dim_size_defns(i)) == "z") then
                collective_q_field%dimensions(i)=1
              else if (trim(writer_entries(writer_entry_index)%contents(my_facet_index+number_q_fields)%dim_size_defns(i)) &
                   == "y") then
                collective_q_field%dimensions(i)=2
              else if (trim(writer_entries(writer_entry_index)%contents(my_facet_index+number_q_fields)%dim_size_defns(i)) &
                   == "x") then
                collective_q_field%dimensions(i)=3
              end if
            end do
            generic=>collective_q_field
            call c_put_generic(collective_q_field_dims, field_name, generic, .false.)
          else
            tot_size=1
            do i=1, writer_entries(writer_entry_index)%contents(my_facet_index+number_q_fields)%dimensions
              tot_size=tot_size*writer_entries(writer_entry_index)%contents(my_facet_index+number_q_fields)%actual_dim_size(i)
            end do
            call c_put_integer(q_field_splits, field_name, tot_size)
          end if
          call c_add_string(q_field_names, field_name)
          add_field_to_writer_entry=number_q_fields
          return
        end if
      end if
      call add_specific_field_to_writer_entry(io_configuration, writer_entry_index, io_config_facet_index, &
           my_facet_index+1, field_name, writer_field_names, duplicate_field_names, prognostic_containing_data_defn%frequency, &
           prognostic_field_configuration=prognostic_field_configuration)
      add_field_to_writer_entry=1
    end if
  end function add_field_to_writer_entry  

  !> Adds a specific field and its information to a writer entry
  !! @param io_configuration The IO server configuration
  !! @param writer_entry_index Index of the writer entry that we are dealing with
  !! @param io_config_facet_index Index of the facet (group) in the IO server configuration
  !! @param my_facet_index The current field index in this internal module representation of the structure
  !! @param field_name The name of the field that we are constructing
  !! @param writer_field_names The field names in the writer (for duplication checking)
  !! @param duplicate_field_names Duplicate field names in the wrier, for duplication checking
  !! @param timestep_frequency Timestepping frequency
  !! @param diagnostic_field_configuration The diagnostic field configuration (optional)
  !! @param prognostic_field_configuration The prognostic field configuration (optional)
  subroutine add_specific_field_to_writer_entry(io_configuration, writer_entry_index, io_config_facet_index, &
       my_facet_index, field_name, writer_field_names, duplicate_field_names, timestep_frequency, &
       diagnostic_field_configuration, prognostic_field_configuration)
    type(io_configuration_type), intent(inout) :: io_configuration
    integer, intent(in) :: writer_entry_index, io_config_facet_index, my_facet_index, timestep_frequency
    character(len=*), intent(in) :: field_name
    type(hashset_type), intent(inout) :: writer_field_names, duplicate_field_names
    type(io_configuration_diagnostic_field_type), intent(inout), optional :: diagnostic_field_configuration
    type(io_configuration_field_type), intent(inout), optional :: prognostic_field_configuration
    
    integer :: i

    writer_entries(writer_entry_index)%contents(my_facet_index)%field_name=field_name

    call c_add_string(used_field_names, field_name)
    
    if (.not. c_contains(writer_field_names, field_name)) then
      call c_add_string(writer_field_names, writer_entries(writer_entry_index)%contents(my_facet_index)%field_name)
    else
      call c_add_string(duplicate_field_names, writer_entries(writer_entry_index)%contents(my_facet_index)%field_name)
    end if

    if (io_configuration%file_writers(writer_entry_index)%contents(io_config_facet_index)%time_manipulation_type == &
         INSTANTANEOUS_TYPE) then
      writer_entries(writer_entry_index)%contents(my_facet_index)%time_manipulation=>perform_instantaneous_time_manipulation
    else if (io_configuration%file_writers(writer_entry_index)%contents(io_config_facet_index)%time_manipulation_type == &
         TIME_AVERAGED_TYPE) then        
      writer_entries(writer_entry_index)%contents(my_facet_index)%time_manipulation=>perform_timeaveraged_time_manipulation
    end if
    writer_entries(writer_entry_index)%contents(my_facet_index)%time_manipulation_type=&
         io_configuration%file_writers(writer_entry_index)%contents(io_config_facet_index)%time_manipulation_type
    writer_entries(writer_entry_index)%contents(my_facet_index)%output_frequency=&
         io_configuration%file_writers(writer_entry_index)%contents(io_config_facet_index)%output_time_frequency
    writer_entries(writer_entry_index)%contents(my_facet_index)%previous_write_time=0.0
    writer_entries(writer_entry_index)%contents(my_facet_index)%previous_tracked_write_point=0.0   
    writer_entries(writer_entry_index)%contents(my_facet_index)%duplicate_field_name=.false.
    writer_entries(writer_entry_index)%contents(my_facet_index)%pending_to_write=.false.
    writer_entries(writer_entry_index)%contents(my_facet_index)%enabled=.false.

    if (present(diagnostic_field_configuration)) then
      writer_entries(writer_entry_index)%contents(my_facet_index)%timestep_frequency=timestep_frequency
      writer_entries(writer_entry_index)%contents(my_facet_index)%dimensions=diagnostic_field_configuration%dimensions
      writer_entries(writer_entry_index)%contents(my_facet_index)%data_type=diagnostic_field_configuration%data_type
      writer_entries(writer_entry_index)%contents(my_facet_index)%field_type=diagnostic_field_configuration%field_type
      writer_entries(writer_entry_index)%contents(my_facet_index)%dim_size_defns=diagnostic_field_configuration%dim_size_defns
      writer_entries(writer_entry_index)%contents(my_facet_index)%units=diagnostic_field_configuration%units
      writer_entries(writer_entry_index)%contents(my_facet_index)%collective_write=diagnostic_field_configuration%collective
    else if (present(prognostic_field_configuration)) then
      writer_entries(writer_entry_index)%contents(my_facet_index)%timestep_frequency=timestep_frequency
      writer_entries(writer_entry_index)%contents(my_facet_index)%data_type=prognostic_field_configuration%data_type
      writer_entries(writer_entry_index)%contents(my_facet_index)%field_type=prognostic_field_configuration%field_type
      writer_entries(writer_entry_index)%contents(my_facet_index)%units=prognostic_field_configuration%units
      writer_entries(writer_entry_index)%contents(my_facet_index)%dimensions=prognostic_field_configuration%dimensions
      writer_entries(writer_entry_index)%contents(my_facet_index)%collective_write=prognostic_field_configuration%collective
      if (prognostic_field_configuration%field_type == ARRAY_FIELD_TYPE) then
        if (prognostic_field_configuration%dimensions .gt. 0) then
          writer_entries(writer_entry_index)%contents(my_facet_index)%dimensions=prognostic_field_configuration%dimensions
          writer_entries(writer_entry_index)%contents(my_facet_index)%dim_size_defns=&
               prognostic_field_configuration%dim_size_defns
        else
          call log_log(LOG_ERROR, "The writing prognostic field '"//trim(field_name)//"' configuration must have dimensions")
        end if
      end if
    else
      call log_log(LOG_ERROR, "A diagnostic or prognostic configuration for the field '"//trim(field_name)//"' was not found")
    end if
    if (writer_entries(writer_entry_index)%contents(my_facet_index)%dimensions .gt. 0) then
      if (writer_entries(writer_entry_index)%contents(my_facet_index)%dim_size_defns(&
           writer_entries(writer_entry_index)%contents(my_facet_index)%dimensions) .eq. "qfields") then
        writer_entries(writer_entry_index)%contents(my_facet_index)%dimensions=&
             writer_entries(writer_entry_index)%contents(my_facet_index)%dimensions-1
      end if      
      do i=1, writer_entries(writer_entry_index)%contents(my_facet_index)%dimensions
        writer_entries(writer_entry_index)%contents(my_facet_index)%actual_dim_size(i)=c_get_integer(&
             io_configuration%dimension_sizing, writer_entries(writer_entry_index)%contents(my_facet_index)%dim_size_defns(i))
      end do      
    end if    
    call check_thread_status(forthread_mutex_init(writer_entries(writer_entry_index)%contents(my_facet_index)%values_mutex, -1))
  end subroutine add_specific_field_to_writer_entry

  !> Marks duplicate field names in a writer entry as duplicates so that the NetCDF layer can then deal with this by issuing
  !! unique names
  !! @param writer_entry The writer state
  !! @param duplicate_field_names A hashset of duplicate field names which are marked
  subroutine handle_duplicate_field_names(writer_entry, duplicate_field_names)
    type(writer_type), intent(inout) :: writer_entry
    type(hashset_type), intent(inout) :: duplicate_field_names

    integer :: i

    do i=1, size(writer_entry%contents)
      if (c_contains(duplicate_field_names, writer_entry%contents(i)%field_name)) then
        writer_entry%contents(i)%duplicate_field_name=.true.
      end if      
    end do    
  end subroutine handle_duplicate_field_names  

  !> Searches the IO server configuration for a group with a specific name and returns the index to that group or 0 if no
  !! corresponding group is found
  !! @param io_configuration The IO server configuration
  !! @param group_name The group name to look up
  !! @returns The index to this group or 0 if no group is found
  integer function get_index_of_group(io_configuration, group_name)
    type(io_configuration_type), intent(inout) :: io_configuration
    character(len=*), intent(in) :: group_name

    integer :: i, entries

    entries=io_configuration%number_of_groups
    do i=1, entries
      if (io_configuration%groups(i)%name == group_name) then
        get_index_of_group=i
        return
      end if
    end do
    get_index_of_group=0
  end function get_index_of_group  
end module writer_federator_mod
