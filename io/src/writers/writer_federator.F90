!> This federates over the writing of diagnostic and prognostic data to the file system. It also manages the time manipulation
!! of fields and groups.
module writer_federator_mod
  use datadefn_mod, only : DEFAULT_PRECISION, STRING_LENGTH
  use configuration_parser_mod, only : TIME_AVERAGED_TYPE, INSTANTANEOUS_TYPE, NONE_TYPE, GROUP_TYPE, FIELD_TYPE, IO_STATE_TYPE, &
       io_configuration_type, io_configuration_field_type, io_configuration_diagnostic_field_type, &
       io_configuration_data_definition_type, data_values_type, get_data_value_by_field_name, get_diagnostic_field_configuration,&
       get_prognostic_field_configuration, get_monc_location
  use none_time_manipulation_mod, only : perform_none_time_manipulation, is_none_time_manipulation_ready_to_write
  use instantaneous_time_manipulation_mod, only : init_instantaneous_manipulation, finalise_instantaneous_manipulation, &
       perform_instantaneous_time_manipulation, is_instantaneous_time_manipulation_ready_to_write
  use timeaveraged_time_manipulation_mod, only : init_time_averaged_manipulation, finalise_time_averaged_manipulation, &
       perform_timeaveraged_time_manipulation, is_time_averaged_time_manipulation_ready_to_write
  use collections_mod, only : queue_type, list_type, map_type, hashmap_type, hashset_type, iterator_type, mapentry_type, &
       c_contains, c_size, c_get_string, c_get_generic, c_get_integer, c_add_string, c_free, c_put_real, c_put_generic, &
       c_key_at, c_is_empty, c_remove, c_push_generic, c_pop_generic, c_real_at, c_get_real, c_get_iterator, &
       c_has_next, c_next_mapentry, c_next_string, c_get_real, c_put_integer, c_add_generic, c_add_string
  use conversions_mod, only : conv_to_string, conv_single_real_to_double, conv_to_integer, conv_to_real
  use io_server_client_mod, only : ARRAY_FIELD_TYPE
  use forthread_mod, only : forthread_mutex_init, forthread_mutex_lock, forthread_mutex_unlock, forthread_mutex_destroy, &
       forthread_rwlock_rdlock, forthread_rwlock_wrlock, forthread_rwlock_unlock, forthread_rwlock_init, forthread_rwlock_destroy
  use threadpool_mod, only : check_thread_status
  use logging_mod, only : LOG_DEBUG, LOG_ERROR, LOG_WARN, log_log, log_master_log, log_get_logging_level, log_is_master
  use writer_types_mod, only : writer_type, writer_field_type, write_field_collective_values_type, pending_write_type, &
       collective_q_field_representation_type, write_field_collective_descriptor_type, write_field_collective_monc_info_type
  use netcdf_filetype_writer_mod, only : initialise_netcdf_filetype, finalise_netcdf_filetype, define_netcdf_file, &
       write_variable, close_netcdf_file, store_io_server_state, get_writer_entry_from_netcdf
  use global_callback_inter_io_mod, only : perform_global_callback
  use data_utils_mod, only : get_scalar_integer_from_monc, get_scalar_real_from_monc, get_scalar_logical_from_monc, &
       is_field_present
  use io_server_client_mod, only : ARRAY_FIELD_TYPE, MAP_FIELD_TYPE, DOUBLE_DATA_TYPE, STRING_DATA_TYPE
  use io_server_state_writer_mod, only : is_io_server_state_writer_ready
  use io_server_state_reader_mod, only : reactivate_writer_federator_state
  use grids_mod, only : Z_INDEX, Y_INDEX, X_INDEX
  use mpi, only : MPI_INT, MPI_MAX
  use mpi_communication_mod, only : lock_mpi, unlock_mpi
  implicit none

#ifndef TEST_MODE
  private
#endif

  type(writer_type), volatile, dimension(:), allocatable :: writer_entries
  type(hashset_type), volatile :: used_field_names, q_field_names
  type(hashmap_type), volatile :: time_points, q_field_splits, collective_q_field_dims

  integer, volatile :: time_points_rwlock, collective_contiguous_initialisation_mutex, currently_writing_mutex
  logical, volatile :: currently_writing

  public initialise_writer_federator, finalise_writer_federator, provide_ordered_field_to_writer_federator, &
       check_writer_for_trigger, issue_actual_write, is_field_used_by_writer_federator, inform_writer_federator_fields_present, &
       inform_writer_federator_time_point, provide_q_field_names_to_writer_federator, is_field_split_on_q
contains

  !> Initialises the write federator and configures it based on the user configuration. Also initialises the time manipulations
  !! @param io_configuration The IO server configuration
  subroutine initialise_writer_federator(io_configuration, diagnostic_generation_frequency, continuation_run)
    type(io_configuration_type), intent(inout) :: io_configuration
    type(hashmap_type), intent(inout) :: diagnostic_generation_frequency
    logical, intent(in) :: continuation_run

    integer :: i, j, number_contents, current_field_index
    type(hashset_type) :: writer_field_names, duplicate_field_names
    
    call check_thread_status(forthread_rwlock_init(time_points_rwlock, -1))
    call check_thread_status(forthread_mutex_init(collective_contiguous_initialisation_mutex, -1))
    call check_thread_status(forthread_mutex_init(currently_writing_mutex, -1))    

    currently_writing=.false.

    call init_time_averaged_manipulation()
    call init_instantaneous_manipulation()
    call initialise_netcdf_filetype()
    
    allocate(writer_entries(io_configuration%number_of_writers))    
    do i=1, io_configuration%number_of_writers
      current_field_index=0
      number_contents=io_configuration%file_writers(i)%number_of_contents
      allocate(writer_entries(i)%contents(get_total_number_writer_fields(io_configuration, i)))
      writer_entries(i)%filename=io_configuration%file_writers(i)%file_name
      writer_entries(i)%title=io_configuration%file_writers(i)%title
      writer_entries(i)%write_on_terminate=io_configuration%file_writers(i)%write_on_terminate
      writer_entries(i)%include_in_io_state_write=io_configuration%file_writers(i)%include_in_io_state_write
      call check_thread_status(forthread_mutex_init(writer_entries(i)%trigger_and_write_mutex, -1))
      call check_thread_status(forthread_mutex_init(writer_entries(i)%num_fields_to_write_mutex, -1))
      call check_thread_status(forthread_mutex_init(writer_entries(i)%pending_writes_mutex, -1))
      writer_entries(i)%write_on_model_time=io_configuration%file_writers(i)%write_on_model_time
      if (writer_entries(i)%write_on_model_time) then
        writer_entries(i)%write_timestep_frequency=0
        writer_entries(i)%write_time_frequency=io_configuration%file_writers(i)%write_time_frequency
      else
        writer_entries(i)%write_time_frequency=0
        writer_entries(i)%write_timestep_frequency=io_configuration%file_writers(i)%write_timestep_frequency
      end if
      writer_entries(i)%previous_write_time=0
      writer_entries(i)%defined_write_time=io_configuration%file_writers(i)%write_time_frequency
      writer_entries(i)%latest_pending_write_time=0
      writer_entries(i)%latest_pending_write_timestep=0      
      writer_entries(i)%contains_io_status_dump=.false.
      do j=1, number_contents
        if (io_configuration%file_writers(i)%contents(j)%facet_type == GROUP_TYPE) then
          current_field_index=add_group_of_fields_to_writer_entry(io_configuration, i, j, current_field_index, &
               writer_field_names, duplicate_field_names, diagnostic_generation_frequency)
        else if (io_configuration%file_writers(i)%contents(j)%facet_type == FIELD_TYPE) then
          current_field_index=current_field_index+add_field_to_writer_entry(io_configuration, &
           i, j, current_field_index, io_configuration%file_writers(i)%contents(j)%facet_name, "", writer_field_names, &
           duplicate_field_names, diagnostic_generation_frequency)
        else if (io_configuration%file_writers(i)%contents(j)%facet_type == IO_STATE_TYPE) then
          writer_entries(i)%contains_io_status_dump=.true.
        end if
      end do
      if (.not. c_is_empty(duplicate_field_names)) call handle_duplicate_field_names(writer_entries(i), duplicate_field_names)
      call c_free(writer_field_names)
      call c_free(duplicate_field_names)
    end do
    if (continuation_run) then
      call reactivate_writer_federator_state(io_configuration, writer_entries, time_points)
    end if    
  end subroutine initialise_writer_federator

  !> Finalises the write federator and the manipulations
  subroutine finalise_writer_federator()
    call check_thread_status(forthread_rwlock_destroy(time_points_rwlock))
    call check_thread_status(forthread_mutex_destroy(collective_contiguous_initialisation_mutex))
    call check_thread_status(forthread_mutex_destroy(currently_writing_mutex))    
    call finalise_time_averaged_manipulation()
    call finalise_instantaneous_manipulation()
    call finalise_netcdf_filetype()
  end subroutine finalise_writer_federator

  subroutine inform_writer_federator_time_point(io_configuration, source, data_id, data_dump)
    type(io_configuration_type), intent(inout) :: io_configuration
    integer, intent(in) :: source, data_id
    character, dimension(:), allocatable, intent(in) :: data_dump

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
  subroutine inform_writer_federator_fields_present(io_configuration, field_names, diag_field_names_and_roots)
    type(io_configuration_type), intent(inout) :: io_configuration
    type(hashset_type), intent(inout), optional :: field_names
    type(hashmap_type), intent(inout), optional :: diag_field_names_and_roots

    type(iterator_type) :: iterator
    character(len=STRING_LENGTH) :: specific_name
    integer :: i, number_q_fields, expected_io
    logical :: field_found, expected_here, diagnostics_mode

    iterator=c_get_iterator(used_field_names)
    expected_io=-1
    do while (c_has_next(iterator))
      specific_name=c_next_string(iterator)
      if (present(field_names)) then
        field_found=c_contains(field_names, specific_name)
        diagnostics_mode=.false.
      else if (present(diag_field_names_and_roots)) then
        field_found=c_contains(diag_field_names_and_roots, specific_name)
        if (field_found) expected_io=c_get_integer(diag_field_names_and_roots, specific_name)
        diagnostics_mode=.true.
      else
        field_found=.false.
      end if
      if (field_found) then
        expected_here=expected_io == -1 .or. expected_io == io_configuration%my_io_rank
        call enable_specific_field_by_name(specific_name, diagnostics_mode, expected_here)      
      end if
    end do
    iterator=c_get_iterator(q_field_names)
    do while (c_has_next(iterator))
      specific_name=c_next_string(iterator)
      if (present(field_names)) then
        field_found=c_contains(field_names, specific_name)
        diagnostics_mode=.false.
      else if (present(diag_field_names_and_roots)) then
        field_found=c_contains(diag_field_names_and_roots, specific_name)
        if (field_found) expected_io=c_get_integer(diag_field_names_and_roots, specific_name)
        diagnostics_mode=.true.
      else
        field_found=.false.
      end if
      if (field_found) then
        expected_here=expected_io == -1 .or. expected_io == io_configuration%my_io_rank
        number_q_fields=c_get_integer(io_configuration%dimension_sizing, "qfields")
        do i=1, number_q_fields
          if (c_size(io_configuration%q_field_names) .ge. i) then
            call enable_specific_field_by_name(trim(specific_name)//"_"//trim(c_get_string(io_configuration%q_field_names, i)), &
                 diagnostics_mode, expected_here)
          else
            call enable_specific_field_by_name(trim(specific_name)//"_udef"//trim(conv_to_string(i)), &
                 diagnostics_mode, expected_here)
          end if
        end do
      end if
    end do    
  end subroutine inform_writer_federator_fields_present  

  !> Determines whether a field is used by the writer federator or not
  !! @param field_name The field name to check whether it is being used or not
  !! @returns Whether this field is used or not
  logical function is_field_used_by_writer_federator(field_name, field_namespace)
    character(len=*), intent(in) :: field_name, field_namespace

    integer :: writer_index, contents_index

    writer_index=1
    contents_index=1
    if (c_contains(used_field_names, field_name)) then      
      is_field_used_by_writer_federator=get_next_applicable_writer_entry(field_name, field_namespace, writer_index, contents_index)
    else
      is_field_used_by_writer_federator=.false.
    end if
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
  subroutine enable_specific_field_by_name(field_name, diagnostics_mode, expected_here)
    character(len=*), intent(in) :: field_name
    logical, intent(in) :: diagnostics_mode
    logical, intent(in), optional :: expected_here

    logical :: continue_search
    integer :: writer_index, contents_index

    continue_search=.true.
    writer_index=1
    contents_index=0
    do while (continue_search)
      contents_index=contents_index+1
      continue_search=get_next_applicable_writer_entry(field_name, writer_index_point=writer_index, &
           contents_index_point=contents_index)
      if (continue_search) then
        if ((writer_entries(writer_index)%contents(contents_index)%diagnostic_field .and. diagnostics_mode) .or. &
             (writer_entries(writer_index)%contents(contents_index)%prognostic_field .and. .not. diagnostics_mode)) then
          writer_entries(writer_index)%contents(contents_index)%enabled=.true.
          if (present(expected_here)) then
            writer_entries(writer_index)%contents(contents_index)%expected_here=expected_here
          end if
        end if
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
          continue_search=get_next_applicable_writer_entry(search_field, writer_index_point=writer_index, &
               contents_index_point=contents_index)
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

  subroutine provide_ordered_field_to_writer_federator(io_configuration, field_name, field_namespace, field_values, &
       timestep, time, source)
    type(io_configuration_type), intent(inout) :: io_configuration
    character(len=*), intent(in) :: field_name, field_namespace
    integer, intent(in) :: timestep, source
    type(data_values_type), target :: field_values
    real(kind=DEFAULT_PRECISION), intent(in) :: time

    integer :: writer_index, contents_index
    logical :: continue_search
    type(data_values_type), pointer :: result_values
    type(hashmap_type) :: typed_result_values
    class(*), pointer :: generic

    if (field_values%data_type == DOUBLE_DATA_TYPE) then
      call provide_ordered_field_to_writer_federator_real_values(io_configuration, field_name, field_namespace, &
           field_values%values, timestep, time, source)
    else if (field_values%data_type == STRING_DATA_TYPE) then
      continue_search=.true.
      writer_index=1
      contents_index=0
      allocate(result_values, source=field_values)
      generic=>result_values
      if (c_contains(used_field_names, field_name)) then
        do while (continue_search)
          contents_index=contents_index+1
          continue_search=get_next_applicable_writer_entry(field_name, field_namespace, writer_index, contents_index)
          if (continue_search) then
            if (.not. writer_entries(writer_index)%contents(contents_index)%enabled) then
              call log_log(LOG_WARN, "Received data for previously un-enabled field '"//&
                   writer_entries(writer_index)%contents(contents_index)%field_name//"'")
            end if
            writer_entries(writer_index)%contents(contents_index)%enabled=.true.          
            writer_entries(writer_index)%contents(contents_index)%latest_timestep_values=timestep
            if (log_get_logging_level() .ge. LOG_DEBUG) then
              call log_log(LOG_DEBUG, "[WRITE FED VALUE STORE] Storing value for field "//trim(field_name)//" ts="//&
                   trim(conv_to_string(timestep))// " t="//trim(conv_to_string(time)))
            end if
            call check_thread_status(forthread_mutex_lock(writer_entries(writer_index)%contents(contents_index)%values_mutex))           
            call c_put_generic(writer_entries(writer_index)%contents(contents_index)%values_to_write, conv_to_string(time), &
                 generic, .false.)
            call check_thread_status(forthread_mutex_unlock(writer_entries(writer_index)%contents(contents_index)%values_mutex))
            if (writer_entries(writer_index)%contents(contents_index)%pending_to_write) then
              call determine_if_outstanding_field_can_be_written(io_configuration, writer_entries(writer_index), &
                   writer_entries(writer_index)%contents(contents_index))
            end if
          end if
        end do
      end if
    end if
  end subroutine provide_ordered_field_to_writer_federator  

  !> Provides fields (either diagnostics or prognostics) to the write federator which will action these as appropriate. This will
  !! split Q fields up if appropriate
  !! @param io_configuration The IO server configuration
  !! @param field_name The field name to write (if appropriate)
  !! @param field_values The field values to write (if appropriate)
  !! @param timestep Corresponding MONC timestep
  !! @param time Corresponding MONC model time
  !! @param source Optional MONC source for the communicated fields
  subroutine provide_ordered_field_to_writer_federator_real_values(io_configuration, field_name, field_namespace, field_values, &
       timestep, time, source)
    type(io_configuration_type), intent(inout) :: io_configuration
    character(len=*), intent(in) :: field_name, field_namespace
    integer, intent(in) :: timestep, source
    real(kind=DEFAULT_PRECISION), dimension(:), intent(in) :: field_values
    real(kind=DEFAULT_PRECISION), intent(in) :: time

    type(iterator_type) :: iterator
    integer :: individual_size, index
    
    if (c_contains(used_field_names, field_name)) then
      call provide_ordered_single_field_to_writer_federator(io_configuration, field_name, field_namespace, field_values, &
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
             trim(field_name)//"_"//trim(c_next_string(iterator)), field_namespace, field_values(index:index+individual_size-1), &
             timestep, time, source)
        index=index+individual_size
      end do
    end if
  end subroutine provide_ordered_field_to_writer_federator_real_values

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
  subroutine provide_ordered_single_field_to_writer_federator(io_configuration, field_name, field_namespace, field_values, &
       timestep, time, source)
    type(io_configuration_type), intent(inout) :: io_configuration
    character(len=*), intent(in) :: field_name, field_namespace
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
        continue_search=get_next_applicable_writer_entry(field_name, field_namespace, writer_index, contents_index)
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
            writer_entries(writer_index)%contents(contents_index)%latest_timestep_values=timestep
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
           writer_entry%previous_write_timestep, writer_entry%write_time, writer_entry%previous_write_time, field_write_success)
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
       timestep, previous_write_timestep, write_time, previous_write_time, field_written)
    type(io_configuration_type), intent(inout) :: io_configuration
    type(writer_type), intent(inout) :: writer_entry
    type(writer_field_type), intent(inout) :: specific_field
    integer, intent(in) :: timestep, previous_write_timestep
    real, intent(in) :: write_time, previous_write_time
    logical, intent(out), optional :: field_written

    real :: value_to_test, largest_value_found
    integer :: num_matching
    logical :: entry_beyond_this_write
    type(iterator_type) :: iterator
    type(mapentry_type) :: map_entry
    type(write_field_collective_values_type), pointer :: multi_monc_entries
    class(*), pointer :: generic
        
    num_matching=0
    largest_value_found=0.0
    entry_beyond_this_write=.false.
    call check_thread_status(forthread_mutex_lock(specific_field%values_mutex))
    if (.not. c_is_empty(specific_field%values_to_write)) then
      iterator=c_get_iterator(specific_field%values_to_write)
      do while (c_has_next(iterator))
        map_entry=c_next_mapentry(iterator)
        value_to_test=conv_to_real(map_entry%key)
        if (specific_field%collective_write) then
          generic=>c_get_generic(map_entry)
          select type(generic)
          type is(write_field_collective_values_type)
            multi_monc_entries=>generic
          end select
          if (c_size(multi_monc_entries%monc_values) .ne. io_configuration%number_of_moncs) cycle
        end if
        if (value_to_test .gt. write_time) entry_beyond_this_write=.true.
        if (value_to_test .le. write_time .and. value_to_test .gt. previous_write_time) then        
          num_matching=num_matching+1
          if (largest_value_found .lt. value_to_test) largest_value_found=value_to_test
        end if
      end do
    end if

    if (num_matching .gt. 0 .and. (specific_field%ready_to_write(largest_value_found, specific_field%output_frequency, write_time, &
         specific_field%latest_timestep_values, timestep) .or. entry_beyond_this_write)) then
      if (.not. specific_field%collective_write .or. .not. specific_field%collective_contiguous_optimisation) then
        if (specific_field%issue_write) then
          call write_variable(io_configuration, specific_field, writer_entry%filename, timestep, write_time)
        end if
        specific_field%previous_write_time=writer_entry%write_time
      end if      
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
    character, dimension(:), allocatable, intent(in) :: data_dump

    integer :: i, timestep
    real(kind=DEFAULT_PRECISION) :: time
    logical :: terminated

    if (is_field_present(io_configuration, source, data_id, "timestep") .and. &
         is_field_present(io_configuration, source, data_id, "time")) then
      timestep=get_scalar_integer_from_monc(io_configuration, source, data_id, data_dump, "timestep")
      time=get_scalar_real_from_monc(io_configuration, source, data_id, data_dump, "time")

      if (is_field_present(io_configuration, source, data_id, "terminated")) then
        terminated=get_scalar_logical_from_monc(io_configuration, source, data_id, data_dump, "terminated")
      else
        terminated=.false.
      end if
      do i=1, size(writer_entries)
        call check_writer_trigger(io_configuration, i, timestep, real(time, kind=4), terminated)
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
  subroutine check_writer_trigger(io_configuration, writer_entry_index, timestep, time, terminated)
    type(io_configuration_type), intent(inout) :: io_configuration
    integer, intent(in) :: writer_entry_index, timestep
    real, intent(in) :: time
    logical, intent(in) :: terminated

    real :: time_difference
    integer :: i
    logical :: issue_write, issue_terminated_write
    
    call check_thread_status(forthread_mutex_lock(writer_entries(writer_entry_index)%trigger_and_write_mutex))
    issue_terminated_write=writer_entries(writer_entry_index)%write_on_terminate .and. terminated
    if (writer_entries(writer_entry_index)%write_on_model_time) then
      time_difference=time-writer_entries(writer_entry_index)%latest_pending_write_time
      issue_write=time_difference .ge. writer_entries(writer_entry_index)%write_time_frequency      
    else
      if (writer_entries(writer_entry_index)%write_timestep_frequency .gt. 0) then
        issue_write=writer_entries(writer_entry_index)%latest_pending_write_timestep .ne. timestep .and. &
             mod(timestep, writer_entries(writer_entry_index)%write_timestep_frequency) == 0
      else
        issue_write=.false.
      end if
      issue_terminated_write=issue_terminated_write .and. &
           writer_entries(writer_entry_index)%latest_pending_write_timestep .ne. timestep
    end if

    if (issue_write .or. issue_terminated_write) then
      writer_entries(writer_entry_index)%latest_pending_write_time=time
      writer_entries(writer_entry_index)%latest_pending_write_timestep=timestep

      call check_thread_status(forthread_mutex_lock(currently_writing_mutex))    

      if (currently_writing) then
        call check_thread_status(forthread_mutex_unlock(currently_writing_mutex))
        call check_thread_status(forthread_mutex_unlock(writer_entries(writer_entry_index)%trigger_and_write_mutex))
        call register_pending_file_write(writer_entry_index, timestep, time, &
             writer_entries(writer_entry_index)%write_on_terminate .and. terminated)
      else
        currently_writing=.true.
        call check_thread_status(forthread_mutex_unlock(currently_writing_mutex))
        call check_thread_status(forthread_mutex_unlock(writer_entries(writer_entry_index)%trigger_and_write_mutex))
        call issue_actual_write(io_configuration, writer_entries(writer_entry_index), timestep, time, &
             writer_entries(writer_entry_index)%write_on_terminate .and. terminated)
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
  subroutine issue_actual_write(io_configuration, writer_entry, timestep, time, terminated_write)
    type(io_configuration_type), intent(inout) :: io_configuration
    type(writer_type), intent(inout) :: writer_entry
    integer, intent(in) :: timestep
    real, intent(in) :: time
    logical, intent(in) :: terminated_write

    integer :: i, j, total_outstanding, num_written, total_flds
    logical :: field_written
    type(map_type) :: applicable_time_points

    call check_thread_status(forthread_mutex_lock(collective_contiguous_initialisation_mutex))
    do i=1, size(writer_entry%contents)
      if (writer_entry%contents(i)%enabled .and. writer_entry%contents(i)%collective_write) then
        if (.not. writer_entry%contents(i)%collective_initialised) then
          call determine_collective_type_and_optimise_if_possible(io_configuration, writer_entry%contents(i))
        end if
      end if
    end do
    call check_thread_status(forthread_mutex_unlock(collective_contiguous_initialisation_mutex))
    
    writer_entry%write_time=time
    writer_entry%write_timestep=timestep
    applicable_time_points=extract_applicable_time_points(writer_entry%previous_write_time, time)    
    call define_netcdf_file(io_configuration, writer_entry, timestep, time, applicable_time_points, terminated_write)
    call c_free(applicable_time_points)
    total_outstanding=0
    total_flds=0
    num_written=0
    call check_thread_status(forthread_mutex_lock(writer_entry%num_fields_to_write_mutex))
    do j=1, size(writer_entry%contents)  
      if (writer_entry%contents(j)%enabled .and. writer_entry%contents(j)%expected_here) then
        total_flds=total_flds+1
        call determine_if_field_can_be_written(io_configuration, writer_entry, writer_entry%contents(j), timestep, &
             writer_entry%previous_write_timestep, time, writer_entry%contents(j)%previous_write_time, field_written)
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
    if (total_outstanding == 0) then
      call close_diagnostics_file(io_configuration, writer_entry, timestep, time)
    end if    
  end subroutine issue_actual_write

  !> Cleans out old timepoints which are no longer going to be of any relavence to the file writing. This ensures that we
  !! don't have lots of stale points that need to be processed and searched beyond but are themselves pointless
  subroutine clean_time_points()
    real :: time_entry
    type(iterator_type) :: iterator
    type(mapentry_type) :: map_entry
    type(list_type) :: removed_entries
    integer :: i
    logical :: remove_timepoint

    call check_thread_status(forthread_rwlock_rdlock(time_points_rwlock))
    iterator=c_get_iterator(time_points)
    do while (c_has_next(iterator))
      map_entry=c_next_mapentry(iterator)
      time_entry=real(c_get_real(map_entry))
      remove_timepoint=.true.
      do i=1, size(writer_entries)
        if (writer_entries(i)%previous_write_time .le. time_entry) then
          remove_timepoint=.false.
          if (writer_entries(i)%write_on_terminate) then
            !print *, "Ignore CTP ", time_entry, writer_entries(i)%previous_write_time
          end if
          exit
        end if
      end do
      if (remove_timepoint) call c_add_string(removed_entries, map_entry%key)
    end do
    iterator=c_get_iterator(removed_entries)
    do while (c_has_next(iterator))
      call c_remove(time_points, c_next_string(iterator))
    end do
    call check_thread_status(forthread_rwlock_unlock(time_points_rwlock))
    call c_free(removed_entries)
  end subroutine clean_time_points


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
    
    type(writer_type), pointer :: writer_entry
    integer :: i
    logical :: terminated, done_chain_run

    writer_entry=>get_writer_entry_from_netcdf(field_name, timestep, terminated)    

    do i=1, size(writer_entry%contents)
      if (writer_entry%contents(i)%enabled .and. writer_entry%contents(i)%collective_write .and. &
           writer_entry%contents(i)%collective_contiguous_optimisation) then        
        call check_thread_status(forthread_mutex_lock(writer_entry%contents(i)%values_mutex))
        call write_variable(io_configuration, writer_entry%contents(i), writer_entry%filename, timestep, writer_entry%write_time)
        writer_entry%contents(i)%previous_write_time=writer_entry%write_time
        call check_thread_status(forthread_mutex_unlock(writer_entry%contents(i)%values_mutex))
      end if      
    end do

    writer_entry%previous_write_time=writer_entry%write_time
    writer_entry%previous_write_timestep=writer_entry%write_timestep
    writer_entry%defined_write_time=writer_entry%defined_write_time+writer_entry%write_time_frequency

    call clean_time_points()

    if (writer_entry%contains_io_status_dump) then
      if (.not. terminated) then
        do while (.not. is_io_server_state_writer_ready(timestep))
        end do
      end if
      call check_thread_status(forthread_rwlock_rdlock(time_points_rwlock))
      call store_io_server_state(io_configuration, writer_entries, time_points, writer_entry, timestep)
      call check_thread_status(forthread_rwlock_unlock(time_points_rwlock))
    end if

    call close_netcdf_file(io_configuration, field_name, timestep)

    done_chain_run=.false.
    do i=1, size(writer_entries)
      if (writer_entries(i)%filename .ne. writer_entry%filename) then
        done_chain_run=check_for_and_issue_chain_write(io_configuration, writer_entries(i))
        if (done_chain_run) exit
      end if
    end do
    if (.not. done_chain_run) done_chain_run=check_for_and_issue_chain_write(io_configuration, writer_entry)

    if (.not. done_chain_run) then
      call check_thread_status(forthread_mutex_lock(currently_writing_mutex))
      currently_writing=.false.
      call check_thread_status(forthread_mutex_unlock(currently_writing_mutex))
      if (log_get_logging_level() .ge. LOG_DEBUG) then
        call log_log(LOG_DEBUG, "No more pending entries to chain to at ts= "//trim(conv_to_string(timestep)))
      end if
    end if
  end subroutine handle_close_diagnostics_globalcallback

  !> Will check whether there are any pending writes and if so will issue a chain write for this
  !! @param io_configuration The IO server configuration
  !! @param writer_entry The specific writer entry
  !! @returns Whether a chain write was issued or not
  logical function check_for_and_issue_chain_write(io_configuration, writer_entry)
    type(io_configuration_type), intent(inout) :: io_configuration
    type(writer_type), intent(inout) :: writer_entry

    class(*), pointer :: generic

    call check_thread_status(forthread_mutex_lock(writer_entry%pending_writes_mutex))
    if (.not. c_is_empty(writer_entry%pending_writes)) then
      check_for_and_issue_chain_write=.true.
      generic=>c_pop_generic(writer_entry%pending_writes)    
      call check_thread_status(forthread_mutex_unlock(writer_entry%pending_writes_mutex))
      select type(generic)
      type is (pending_write_type)
        if (log_get_logging_level() .ge. LOG_DEBUG) then
          call log_log(LOG_DEBUG, "Chain to next pending entry ts= "//trim(conv_to_string(generic%timestep)))
        end if
        call issue_actual_write(io_configuration, writer_entry, generic%timestep, &
             generic%write_time, generic%terminated_write)
        deallocate(generic)
      end select
    else
      check_for_and_issue_chain_write=.false.
      call check_thread_status(forthread_mutex_unlock(writer_entry%pending_writes_mutex))
    end if    
  end function check_for_and_issue_chain_write
  
  !> Registers a pending file write which will be actioned later on
  !! @param writer_entry_index Index of the writer entry
  !! @param timestep The timestep that the pending write represents
  !! @param time The time of the pending write
  subroutine register_pending_file_write(writer_entry_index, timestep, time, terminated_write)
    integer, intent(in) :: writer_entry_index, timestep
    real, intent(in) :: time
    logical, intent(in) :: terminated_write

    type(pending_write_type), pointer :: pending_write
    class(*), pointer :: generic

    allocate(pending_write)
    pending_write%write_time=time
    pending_write%timestep=timestep
    pending_write%terminated_write=terminated_write

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
  logical function get_next_applicable_writer_entry(field_name, field_namespace, writer_index_point, contents_index_point)
    character(len=*), intent(in) :: field_name
    character(len=*), intent(in), optional :: field_namespace
    integer, intent(inout) :: writer_index_point, contents_index_point

    integer :: i, j

    if (writer_index_point .le. size(writer_entries)) then    
      do i=writer_index_point, size(writer_entries)
        if (contents_index_point .le. size(writer_entries(i)%contents)) then
          do j=contents_index_point, size(writer_entries(i)%contents)           
            if (writer_entries(i)%contents(j)%field_name==field_name) then
              if (present(field_namespace)) then                
                if (writer_entries(i)%contents(j)%field_namespace .ne. field_namespace) cycle
              end if
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
             get_group_number_of_fields(io_configuration, io_configuration%groups(group_index)%members, number_q_fields, &
             io_configuration%groups(group_index)%namespace)
      else if (io_configuration%file_writers(writer_entry_index)%contents(i)%facet_type == FIELD_TYPE) then
        ! NSE
        get_total_number_writer_fields=get_total_number_writer_fields+get_field_number_of_fields(io_configuration, &
             io_configuration%file_writers(writer_entry_index)%contents(i)%facet_name, "", number_q_fields)
      end if
    end do
  end function get_total_number_writer_fields

  !> Retrieves the number of fields within a group of fields
  !! @param io_configuration The IO server configuration
  !! @param group_members The members of the group
  !! @param num_q_fields The number of Q fields
  !! @returns The number of fields that make up this group
  integer function get_group_number_of_fields(io_configuration, group_members, num_q_fields, namespace)
    type(io_configuration_type), intent(inout) :: io_configuration
    type(list_type) :: group_members
    integer, intent(in) :: num_q_fields
    character(len=STRING_LENGTH), intent(in) :: namespace

    type(iterator_type) :: iterator
    character(len=STRING_LENGTH) :: field_name
    
    get_group_number_of_fields=0
    iterator=c_get_iterator(group_members)
    do while (c_has_next(iterator))
      field_name=c_next_string(iterator)
      get_group_number_of_fields=get_group_number_of_fields+get_field_number_of_fields(io_configuration, field_name, namespace, &
           num_q_fields)
    end do
  end function get_group_number_of_fields

  !> Retrieves the number of fields that make up this field, if it is a Q field then it will be split into many subfields
  !! hence it is not a simple 1-1 mapping
  !! @param io_configuration The IO server configuration
  !! @param field_name The name of the field
  !! @param num_q_fields The number of Q fields
  !! @returns The number of fields that make up this field
  integer function get_field_number_of_fields(io_configuration, field_name, field_namespace, num_q_fields)
    type(io_configuration_type), intent(inout) :: io_configuration
    character(len=STRING_LENGTH), intent(in) :: field_name, field_namespace
    integer, intent(in) :: num_q_fields

    type(io_configuration_field_type) :: prognostic_field_configuration
    type(io_configuration_data_definition_type) :: prognostic_containing_data_defn
    type(io_configuration_diagnostic_field_type) :: diagnostic_field_configuration
    
    if (get_diagnostic_field_configuration(io_configuration, field_name, field_namespace, diagnostic_field_configuration)) then
      if (diagnostic_field_configuration%field_type == ARRAY_FIELD_TYPE) then
        if (diagnostic_field_configuration%dim_size_defns(diagnostic_field_configuration%dimensions) .eq. "qfields") then
          get_field_number_of_fields=num_q_fields
          return
        end if
      end if
      get_field_number_of_fields=1
    else if (get_prognostic_field_configuration(io_configuration, field_name, field_namespace, &
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
           writer_entry_index, facet_index, add_group_of_fields_to_writer_entry, field_name, &
           io_configuration%groups(group_index)%namespace, writer_field_names, duplicate_field_names, &
           diagnostic_generation_frequency)
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
       my_facet_index, field_name, field_namespace, writer_field_names, duplicate_field_names, diagnostic_generation_frequency)
    type(io_configuration_type), intent(inout) :: io_configuration
    integer, intent(in) :: writer_entry_index, io_config_facet_index, my_facet_index
    character(len=*), intent(in) :: field_name, field_namespace
    type(hashset_type), intent(inout) :: writer_field_names, duplicate_field_names
    type(hashmap_type), intent(inout) :: diagnostic_generation_frequency

    integer :: i, number_q_fields, tot_size
    type(io_configuration_field_type) :: prognostic_field_configuration
    type(io_configuration_data_definition_type) :: prognostic_containing_data_defn
    type(io_configuration_diagnostic_field_type) :: diagnostic_field_configuration
    type(collective_q_field_representation_type), pointer :: collective_q_field
    class(*), pointer :: generic
    
    if (get_diagnostic_field_configuration(io_configuration, field_name, field_namespace, diagnostic_field_configuration)) then
      if (diagnostic_field_configuration%field_type == ARRAY_FIELD_TYPE) then
        if (diagnostic_field_configuration%dim_size_defns(diagnostic_field_configuration%dimensions) .eq. "qfields") then  
          number_q_fields=c_get_integer(io_configuration%dimension_sizing, "qfields")
          do i=1, number_q_fields
            call add_specific_field_to_writer_entry(io_configuration, writer_entry_index, io_config_facet_index, &
                 my_facet_index+i, trim(field_name)//"_udef"//trim(conv_to_string(i)), field_namespace, writer_field_names, &
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
           my_facet_index+1, field_name, field_namespace, writer_field_names, duplicate_field_names, &
           c_get_integer(diagnostic_generation_frequency, field_name), diagnostic_field_configuration)
      add_field_to_writer_entry=1      
    else if (get_prognostic_field_configuration(io_configuration, field_name, field_namespace, &
         prognostic_field_configuration, prognostic_containing_data_defn)) then
      if (prognostic_field_configuration%field_type == ARRAY_FIELD_TYPE) then
        if (prognostic_field_configuration%dim_size_defns(prognostic_field_configuration%dimensions) .eq. "qfields") then  
          number_q_fields=c_get_integer(io_configuration%dimension_sizing, "qfields")
          do i=1, number_q_fields
            call add_specific_field_to_writer_entry(io_configuration, writer_entry_index, io_config_facet_index, &
                 my_facet_index+i, trim(field_name)//"_udef"//trim(conv_to_string(i)), field_namespace, writer_field_names, &
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
           my_facet_index+1, field_name, field_namespace, writer_field_names, duplicate_field_names, &
           prognostic_containing_data_defn%frequency, prognostic_field_configuration=prognostic_field_configuration)
      add_field_to_writer_entry=1
    else
      call log_log(LOG_ERROR, "Field '"//trim(field_name)//&
           "' configured for file write but can not find this as a prognostic or diagnostic definition")
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
       my_facet_index, field_name, field_namespace, writer_field_names, duplicate_field_names, timestep_frequency, &
       diagnostic_field_configuration, prognostic_field_configuration)
    type(io_configuration_type), intent(inout) :: io_configuration
    integer, intent(in) :: writer_entry_index, io_config_facet_index, my_facet_index, timestep_frequency
    character(len=*), intent(in) :: field_name, field_namespace
    type(hashset_type), intent(inout) :: writer_field_names, duplicate_field_names
    type(io_configuration_diagnostic_field_type), intent(inout), optional :: diagnostic_field_configuration
    type(io_configuration_field_type), intent(inout), optional :: prognostic_field_configuration
    
    integer :: i

    writer_entries(writer_entry_index)%contents(my_facet_index)%field_name=field_name
    writer_entries(writer_entry_index)%contents(my_facet_index)%field_namespace=field_namespace

    call c_add_string(used_field_names, field_name)
    
    if (.not. c_contains(writer_field_names, field_name)) then
      call c_add_string(writer_field_names, writer_entries(writer_entry_index)%contents(my_facet_index)%field_name)
    else
      call c_add_string(duplicate_field_names, writer_entries(writer_entry_index)%contents(my_facet_index)%field_name)
    end if

    if (io_configuration%file_writers(writer_entry_index)%contents(io_config_facet_index)%time_manipulation_type == &
         INSTANTANEOUS_TYPE) then
      writer_entries(writer_entry_index)%contents(my_facet_index)%time_manipulation=>perform_instantaneous_time_manipulation
      writer_entries(writer_entry_index)%contents(my_facet_index)%ready_to_write=>is_instantaneous_time_manipulation_ready_to_write
    else if (io_configuration%file_writers(writer_entry_index)%contents(io_config_facet_index)%time_manipulation_type == &
         TIME_AVERAGED_TYPE) then        
      writer_entries(writer_entry_index)%contents(my_facet_index)%time_manipulation=>perform_timeaveraged_time_manipulation
      writer_entries(writer_entry_index)%contents(my_facet_index)%ready_to_write=>is_time_averaged_time_manipulation_ready_to_write
    else if (io_configuration%file_writers(writer_entry_index)%contents(io_config_facet_index)%time_manipulation_type == &
         NONE_TYPE) then        
      writer_entries(writer_entry_index)%contents(my_facet_index)%time_manipulation=>perform_none_time_manipulation
      writer_entries(writer_entry_index)%contents(my_facet_index)%ready_to_write=>is_none_time_manipulation_ready_to_write
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
    writer_entries(writer_entry_index)%contents(my_facet_index)%expected_here=.true.
    writer_entries(writer_entry_index)%contents(my_facet_index)%prognostic_field=.false.
    writer_entries(writer_entry_index)%contents(my_facet_index)%diagnostic_field=.false.

    if (present(diagnostic_field_configuration)) then
      writer_entries(writer_entry_index)%contents(my_facet_index)%timestep_frequency=timestep_frequency
      writer_entries(writer_entry_index)%contents(my_facet_index)%dimensions=diagnostic_field_configuration%dimensions
      writer_entries(writer_entry_index)%contents(my_facet_index)%data_type=diagnostic_field_configuration%data_type
      writer_entries(writer_entry_index)%contents(my_facet_index)%field_type=diagnostic_field_configuration%field_type
      writer_entries(writer_entry_index)%contents(my_facet_index)%dim_size_defns=diagnostic_field_configuration%dim_size_defns
      writer_entries(writer_entry_index)%contents(my_facet_index)%units=diagnostic_field_configuration%units
      writer_entries(writer_entry_index)%contents(my_facet_index)%collective_write=diagnostic_field_configuration%collective
      writer_entries(writer_entry_index)%contents(my_facet_index)%collective_initialised=.false.
      writer_entries(writer_entry_index)%contents(my_facet_index)%issue_write=.true.
      writer_entries(writer_entry_index)%contents(my_facet_index)%diagnostic_field=.true.
    else if (present(prognostic_field_configuration)) then
      writer_entries(writer_entry_index)%contents(my_facet_index)%timestep_frequency=timestep_frequency
      writer_entries(writer_entry_index)%contents(my_facet_index)%data_type=prognostic_field_configuration%data_type
      writer_entries(writer_entry_index)%contents(my_facet_index)%field_type=prognostic_field_configuration%field_type
      writer_entries(writer_entry_index)%contents(my_facet_index)%units=prognostic_field_configuration%units
      writer_entries(writer_entry_index)%contents(my_facet_index)%dimensions=prognostic_field_configuration%dimensions
      writer_entries(writer_entry_index)%contents(my_facet_index)%collective_write=prognostic_field_configuration%collective
      writer_entries(writer_entry_index)%contents(my_facet_index)%collective_initialised=.false.
      writer_entries(writer_entry_index)%contents(my_facet_index)%prognostic_field=.true.
      if (.not. prognostic_field_configuration%collective) then
        writer_entries(writer_entry_index)%contents(my_facet_index)%issue_write=io_configuration%my_io_rank==0
      else
        writer_entries(writer_entry_index)%contents(my_facet_index)%issue_write=.true.
      end if      
      if (prognostic_field_configuration%field_type == ARRAY_FIELD_TYPE .or. &
           prognostic_field_configuration%field_type == MAP_FIELD_TYPE) then
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

  !> Determines whether it can optimise a specific collective field. If the field fits into certain limited parameters
  !! then it will optimise it. These parameters are very common, hence most fields can be optimised. Basically, it is looking
  !! to contiguous blocks of data from different MONCs so that the number of writes to the NetCDF file is limited
  !! @param io_configuration The IO server configuration
  !! @param field_to_write_information Description of the the specific field that we are looking at here
  subroutine determine_collective_type_and_optimise_if_possible(io_configuration, field_to_write_information)
    type(io_configuration_type), intent(inout) :: io_configuration
    type(writer_field_type), intent(inout) :: field_to_write_information

    if (field_to_write_information%dimensions .eq. 3 .and. &
         get_dimension_identifier(field_to_write_information%dim_size_defns(1)) == Z_INDEX .and. &
         get_dimension_identifier(field_to_write_information%dim_size_defns(2)) == Y_INDEX .and. &
         get_dimension_identifier(field_to_write_information%dim_size_defns(3)) == X_INDEX) then
      field_to_write_information%collective_contiguous_optimisation=.true.
      call initialise_contiguous_data_regions(io_configuration, field_to_write_information)
    else
      field_to_write_information%collective_contiguous_optimisation=.false.
    end if
    field_to_write_information%collective_initialised=.true.
  end subroutine determine_collective_type_and_optimise_if_possible  

  !> Will initialise the collective data regions that form contiguous blocks within the data. This is quite an expensive
  !! operation so only done once for each field, but has the potential for very significant performance advantages for the
  !! fields that match it
  !! @param io_configuration The IO server configuration
  !! @param field_to_write_information The specific field that is being looked at and identified for contiguous data
  subroutine initialise_contiguous_data_regions(io_configuration, field_to_write_information)
    type(io_configuration_type), intent(inout) :: io_configuration
    type(writer_field_type), intent(inout) :: field_to_write_information

    integer :: start(field_to_write_information%dimensions, io_configuration%number_of_moncs), &
         count(field_to_write_information%dimensions, io_configuration%number_of_moncs), &
         common_starters(io_configuration%number_of_moncs), num_common, num_current_contents, active_dim, other_dim, &
         j, k, i, dim_identifier, number_distinct_writes, start_blocks(io_configuration%number_of_moncs), ierr, &
         count_blocks(io_configuration%number_of_moncs), current_contents(io_configuration%number_of_moncs), &
         monc_write_start_offset_per_dim(field_to_write_information%dimensions,io_configuration%number_of_moncs)
    logical :: processed(io_configuration%number_of_moncs)

    type(write_field_collective_descriptor_type), pointer :: collective_descriptor
    type(write_field_collective_monc_info_type), pointer :: specific_monc_collective
    class(*), pointer :: generic  

    processed=.false.
    number_distinct_writes=0
    do j=1, io_configuration%number_of_moncs
      do k=1, field_to_write_information%dimensions
        dim_identifier=get_dimension_identifier(field_to_write_information%dim_size_defns(k))
        start(k, j)=io_configuration%registered_moncs(j)%local_dim_starts(dim_identifier)
        count(k, j)=io_configuration%registered_moncs(j)%local_dim_sizes(dim_identifier)
      end do
    end do

    do j=1, io_configuration%number_of_moncs
      if (.not. processed(j)) then
        call get_common_starts(Y_INDEX, start(Y_INDEX, j), start, common_starters, num_common)
        if (num_common == 0) then
          call get_common_starts(X_INDEX, start(X_INDEX, j), start, common_starters, num_common)
          if (num_common .gt. 0) then
            active_dim=X_INDEX
            other_dim=Y_INDEX
          end if
        else
          active_dim=Y_INDEX
          other_dim=X_INDEX
        end if
          number_distinct_writes=number_distinct_writes+1
          allocate(collective_descriptor)
          allocate(collective_descriptor%absolute_start(field_to_write_information%dimensions), &
               collective_descriptor%count(field_to_write_information%dimensions))
          start_blocks(number_distinct_writes)=start(other_dim, j)
          count_blocks(number_distinct_writes)=count(other_dim, j)
          num_current_contents=1
          current_contents(num_current_contents)=j
          processed(j)=.true.
          if (num_common .gt. 0) then
            do k=1, num_common
              do i=1, num_common
                if (.not. processed(common_starters(i)) .and. count(active_dim, j) == count(active_dim, i)) then
                  if (start(other_dim, common_starters(i)) .lt. start_blocks(number_distinct_writes) .and. &
                       start(other_dim, common_starters(i)) + count(other_dim, common_starters(i)) &
                       == start_blocks(number_distinct_writes)) then
                    start_blocks(number_distinct_writes)=start(other_dim, common_starters(i))
                    count_blocks(number_distinct_writes)=count_blocks(number_distinct_writes)+count(other_dim, common_starters(i))
                    processed(common_starters(i))=.true.                    
                    num_current_contents=num_current_contents+1
                    current_contents(num_current_contents)=common_starters(i)
                  else if (start(other_dim, common_starters(i)) .gt. start_blocks(number_distinct_writes) .and. &
                       start_blocks(number_distinct_writes) + count_blocks(number_distinct_writes) &
                       == start(other_dim, common_starters(i))) then
                    count_blocks(number_distinct_writes)=count_blocks(number_distinct_writes)+count(other_dim, common_starters(i))
                    processed(common_starters(i))=.true.                    
                    num_current_contents=num_current_contents+1
                    current_contents(num_current_contents)=common_starters(i)
                  end if
                end if
              end do
            end do
          end if
          collective_descriptor%absolute_start=start(:,j)
          collective_descriptor%count=count(:,j)
          collective_descriptor%absolute_start(other_dim)=start_blocks(number_distinct_writes)
          collective_descriptor%count(other_dim)=count_blocks(number_distinct_writes)
          collective_descriptor%split_dim=other_dim
          if (num_current_contents .gt. 0) then
            do k=1, num_current_contents
              allocate(specific_monc_collective)                          
              specific_monc_collective%relative_dimension_start=(start(other_dim,current_contents(k))-&
                   start_blocks(number_distinct_writes)) + 1              
              specific_monc_collective%counts=count(:, current_contents(k))
              specific_monc_collective%monc_location=current_contents(k)
              specific_monc_collective%monc_source=io_configuration%registered_moncs(current_contents(k))%source_id
              generic=>specific_monc_collective
              call c_add_generic(collective_descriptor%specific_monc_info, generic, .false.)
            end do
          end if
          generic=>collective_descriptor
          call c_add_generic(field_to_write_information%collective_descriptors, generic, .false.)
      end if
    end do
    call lock_mpi()
    call mpi_iallreduce(number_distinct_writes, field_to_write_information%max_num_collective_writes, 1, MPI_INT, MPI_MAX, &
         io_configuration%io_communicator, field_to_write_information%max_num_collective_writes_request_handle, ierr)
    call unlock_mpi()
  end subroutine initialise_contiguous_data_regions

  !> Retrieves the number of common starting points that match a specific input value
  !! @param dim The dimension we are searching in
  !! @param val The value that we are looking to match against
  !! @param vals All the other dimensions information to search against
  !! @param common_starters The common starting point locations
  !! @param num_common The number of common starting points identified
  subroutine get_common_starts(dim, val, vals, common_starters, num_common)
    integer, intent(in) :: dim, val
    integer, dimension(:,:), intent(in) :: vals
    integer, dimension(:), intent(out) :: common_starters
    integer, intent(out) :: num_common

    integer :: i

    num_common=0
    do i=1, size(vals, 2)
      if (vals(dim, i) == val) then
        num_common=num_common+1
        common_starters(num_common)=i
      end if      
    end do    
  end subroutine get_common_starts

  !> Translates a dimension name to its numeric corresponding identifier
  !! @param dim_name The name of the dimension to look up
  !! @param is_auto_dimension Optional parameter determining whether dimension is auto or not
  !! @returns Corresponding identifier
  integer function get_dimension_identifier(dim_name, is_auto_dimension)
    character(len=*), intent(in) :: dim_name
    logical, intent(out), optional :: is_auto_dimension

    integer :: dash_idx
    logical :: is_modified_size

    dash_idx=index(dim_name, "_")
    dash_idx=dash_idx-1
    is_modified_size=dash_idx .ne. -1
    if (.not. is_modified_size) dash_idx=len_trim(dim_name)

    if (dim_name(:dash_idx) .eq. "z" .or. dim_name(:dash_idx) .eq. "zn") then
      get_dimension_identifier=Z_INDEX
    else if (dim_name(:dash_idx) .eq. "y") then
      get_dimension_identifier=Y_INDEX
    else if (dim_name(:dash_idx) .eq. "x") then
      get_dimension_identifier=X_INDEX
    else
      get_dimension_identifier=-1
    end if

    if (present(is_auto_dimension)) is_auto_dimension=is_modified_size
  end function get_dimension_identifier
end module writer_federator_mod
