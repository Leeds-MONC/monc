!> The writer field manager will manage aspects of the fields being provided to the writer federator. There are two major
!! aspects to this, firstly extraction of fields from the data sent directly from MONC and then passed on if these
!! fields are needed. Secondly, this will queue up the data to ensure that it is sent to the writer federator in strict
!! order and that only one piece of a fields data is in the writer federator at any one time
module writer_field_manager_mod
  use datadefn_mod, only : DEFAULT_PRECISION, STRING_LENGTH
  use configuration_parser_mod, only : io_configuration_type, data_values_type, io_configuration_field_type       
  use collections_mod, only : hashmap_type, map_type, iterator_type, mapentry_type, c_size, c_get_iterator, c_has_next, &
       c_next_mapentry, c_contains, c_get_string, c_get_generic, c_put_generic, c_is_empty, c_remove
  use conversions_mod, only : conv_to_string, conv_to_real
  use forthread_mod, only : forthread_mutex_init, forthread_mutex_lock, forthread_mutex_unlock, forthread_mutex_destroy, &
       forthread_rwlock_rdlock, forthread_rwlock_wrlock, forthread_rwlock_unlock, forthread_rwlock_init, forthread_rwlock_destroy
  use threadpool_mod, only : check_thread_status
  use data_utils_mod, only : get_scalar_integer_from_monc, get_scalar_real_from_monc, get_scalar_logical_from_monc, &
       is_field_present, get_array_double_from_monc, get_array_integer_from_monc, get_string_from_monc, get_map_from_monc, &
       unpack_scalar_integer_from_bytedata, unpack_scalar_dp_real_from_bytedata
  use io_server_client_mod, only : DOUBLE_DATA_TYPE, INTEGER_DATA_TYPE, STRING_DATA_TYPE, SCALAR_FIELD_TYPE, MAP_FIELD_TYPE
  use logging_mod, only : LOG_WARN, LOG_ERROR, log_log
  use writer_federator_mod, only : is_field_used_by_writer_federator, provide_ordered_field_to_writer_federator, &
       is_field_split_on_q, is_field_split_on_tracer 
  use writer_types_mod, only : prepare_to_serialise_data_values_type, serialise_data_values_type, unserialise_data_values_type
  use optionsdatabase_mod, only : options_get_logical
  use io_server_client_mod, only : pack_scalar_field
  use io_server_state_writer_mod, only : set_serialise_write_field_manager_state
  use io_server_state_reader_mod, only : reactivate_writer_field_manager_state
  implicit none

#ifndef TEST_MODE
  private
#endif

  type field_ordering_value_type
     character(len=STRING_LENGTH) :: field_name, field_namespace
     integer :: timestep, frequency, source
     type(data_values_type), allocatable :: field_values
     real(kind=DEFAULT_PRECISION) :: time
  end type field_ordering_value_type

  type field_ordering_type
     type(hashmap_type) :: timestep_to_value
     integer :: access_mutex, last_timestep_access, frequency, last_time_access
  end type field_ordering_type

  interface provide_field_to_writer_federator
     module procedure provide_field_to_writer_federator_rvalues_src, provide_field_to_writer_federator_rvalues_nosrc, &
          provide_field_to_writer_federator_src, provide_field_to_writer_federator_nosrc
  end interface provide_field_to_writer_federator

  integer, volatile :: field_lock
  type(hashmap_type), volatile :: field_orderings
  logical :: time_basis
  real(kind=DEFAULT_PRECISION) :: model_initial_time

  public initialise_writer_field_manager, finalise_writer_field_manager, provide_monc_data_to_writer_federator, &
       provide_field_to_writer_federator, is_write_field_manager_up_to_date
contains

  !> Initialises the writer field manager
  !! @param io_configuration Configuration of the IO server
  !! @param continuation_run Whether or not this is a continuation run
  subroutine initialise_writer_field_manager(io_configuration, continuation_run, reconfig_initial_time)
    type(io_configuration_type), intent(inout) :: io_configuration
    logical, intent(in) :: continuation_run
    real(kind=DEFAULT_PRECISION), intent(in) :: reconfig_initial_time

    model_initial_time = reconfig_initial_time
    time_basis = options_get_logical(io_configuration%options_database,"time_basis")

    call check_thread_status(forthread_rwlock_init(field_lock, -1))
    call set_serialise_write_field_manager_state(serialise_field_manager_state, prepare_to_serialise_field_manager_state, &
         is_write_field_manager_up_to_date)
    if (continuation_run) then
      call reactivate_writer_field_manager_state(io_configuration, unserialise_field_manager_state)
    end if
  end subroutine initialise_writer_field_manager

  !> Finalises the writer field manager
  subroutine finalise_writer_field_manager()
    call check_thread_status(forthread_rwlock_destroy(field_lock))
  end subroutine finalise_writer_field_manager

  !> Determines whether the state of the write field manager is up to date with respect to the timestep that has been
  !! provided, checks each entry to see if it is missing data or not
  !! @param timestep The timestep to check against
  !! @returns Whether the write field manager is up to date or not
  logical function is_write_field_manager_up_to_date(timestep)
    integer, intent(in) :: timestep

    type(iterator_type) :: iterator
    type(mapentry_type) :: mapentry
    class(*), pointer :: generic

    is_write_field_manager_up_to_date=.true.
    call check_thread_status(forthread_rwlock_rdlock(field_lock))
    iterator=c_get_iterator(field_orderings)
    do while (c_has_next(iterator))
      mapentry=c_next_mapentry(iterator)
      generic=>c_get_generic(mapentry)
      select type(generic)
        type is (field_ordering_type)
          if (generic%last_timestep_access .lt. timestep) then
            is_write_field_manager_up_to_date=.false.
            exit
          end if
      end select      
    end do    
    call check_thread_status(forthread_rwlock_unlock(field_lock))
  end function is_write_field_manager_up_to_date  

  !> Data communicated from MONC is provided to this write federator and then included if the configuration has selected
  !! to use that MONC field (as opposed to a field produced by the diagnostics federator)
  !! @param io_configuration Configuration of the IO server
  !! @param source The source PID of the MONC process
  !! @param data_id The ID of the data definition that is represented by the dump
  !! @param data_dump The data we have received from MONC
  subroutine provide_monc_data_to_writer_federator(io_configuration, source, data_id, data_dump)
    type(io_configuration_type), intent(inout) :: io_configuration
    integer, intent(in) :: source, data_id
    character, dimension(:), allocatable, intent(in) :: data_dump

    integer :: i, num_fields, timestep
    character(len=STRING_LENGTH) :: field_name, field_namespace
    real(kind=DEFAULT_PRECISION) :: time
    type(data_values_type) :: monc_value
    logical :: terminated_case

    if (is_field_present(io_configuration, source, data_id, "timestep") .and. &
         is_field_present(io_configuration, source, data_id, "time")) then
      timestep=get_scalar_integer_from_monc(io_configuration, source, data_id, data_dump, "timestep")
      time=get_scalar_real_from_monc(io_configuration, source, data_id, data_dump, "time")

      if (is_field_present(io_configuration, source, data_id, "terminated")) then
        terminated_case=get_scalar_logical_from_monc(io_configuration, source, data_id, data_dump, "terminated") .and. &
             io_configuration%data_definitions(data_id)%send_on_terminate
      else
        terminated_case=.false.
      end if

      num_fields=io_configuration%data_definitions(data_id)%number_of_data_fields

      do i=1, num_fields
        field_name=io_configuration%data_definitions(data_id)%fields(i)%name
        field_namespace=io_configuration%data_definitions(data_id)%fields(i)%namespace
        if (is_field_present(io_configuration, source, data_id, field_name) .and. &
             (is_field_used_by_writer_federator(field_name, field_namespace) .or. &
             is_field_split_on_q(field_name) .or. is_field_split_on_tracer(field_name))) then
          monc_value=get_value_from_monc_data(io_configuration, source, data_id, data_dump, field_name, field_namespace)
          call provide_field_to_writer_federator_src(io_configuration, field_name, field_namespace, monc_value, timestep, time, &
               io_configuration%data_definitions(data_id)%frequency, source, terminated_case)
          !<<<< CONSIDER WHY DEALLOC IS COMMENTED OUT
          !deallocate(monc_value)
        end if
      end do
    else
      call log_log(LOG_WARN, "Can not run pass MONC fields to writer federator without a time and timestep")
    end if
  end subroutine provide_monc_data_to_writer_federator

  !> Retrieves a value from the communicated MONC data. If this was an integer then converts to a real
  !! @param io_configuration Configuration of the IO server
  !! @param source The source PID of the MONC process
  !! @param data_id The ID of the data definition that is represented by the dump
  !! @param data_dump The data we have received from MONC
  !! @param field_name The field to retrieve
  !! @returns Data value wrapper containing the values and meta-data
  function get_value_from_monc_data(io_configuration, source, data_id, data_dump, field_name, field_namespace)
    type(io_configuration_type), intent(inout) :: io_configuration
    integer, intent(in) :: source, data_id
    character, dimension(:), allocatable, intent(in) :: data_dump
    character(len=*), intent(in) :: field_name, field_namespace
    type(data_values_type) :: get_value_from_monc_data
    type(map_type) :: retrieved_map
    type(iterator_type) :: iterator
    type(mapentry_type) :: map_entry

    integer :: field_data_type, field_field_type, i
    integer, dimension(:), allocatable :: int_values
    real(kind=DEFAULT_PRECISION), dimension(:), allocatable :: real_values

    call get_type_of_field(io_configuration%data_definitions(data_id)%fields, field_name, field_namespace, &
         field_field_type, field_data_type)
    if (field_data_type == 0) then
      call log_log(LOG_ERROR, "No data type for field '"//trim(field_name)//"'")
    end if

    if (field_data_type == DOUBLE_DATA_TYPE) then
      real_values=get_array_double_from_monc(io_configuration, source, data_id, data_dump, field_name)
      allocate(get_value_from_monc_data%values(size(real_values)))
      do i=1,size(real_values)
        get_value_from_monc_data%values(i)=real_values(i)
      end do
      get_value_from_monc_data%data_type=DOUBLE_DATA_TYPE
      deallocate(real_values)
    else if (field_data_type == INTEGER_DATA_TYPE) then
      get_value_from_monc_data%data_type=DOUBLE_DATA_TYPE
      int_values=get_array_integer_from_monc(io_configuration, source, data_id, data_dump, field_name)
      allocate(get_value_from_monc_data%values(size(int_values)))
      do i=1, size(int_values)
        get_value_from_monc_data%values(i)=conv_to_real(int_values(i))
      end do
      deallocate(int_values)
    else if (field_data_type == STRING_DATA_TYPE) then
      get_value_from_monc_data%data_type=STRING_DATA_TYPE
      if (field_field_type == SCALAR_FIELD_TYPE) then
        allocate(get_value_from_monc_data%string_values(1))
        get_value_from_monc_data%string_values(1)=get_string_from_monc(io_configuration, source, data_id, data_dump, field_name)
      else if (field_field_type == MAP_FIELD_TYPE) then        
        get_value_from_monc_data%map_values=get_map_from_monc(io_configuration, source, data_id, data_dump, field_name)        
      end if
    end if
  end function get_value_from_monc_data

  !> Retrieves the data type of a field or 0 if the field was not found
  !! @param fields Array of fields to search
  !! @param field_name The name of the field to locate
  !! @returns The data type of this field or 0 if the field was not found
  subroutine get_type_of_field(fields, field_name, field_namespace, field_type, data_type)
    type(io_configuration_field_type), dimension(:), intent(in) :: fields
    character(len=*), intent(in) :: field_name, field_namespace
    integer, intent(out) :: field_type, data_type

    integer :: i

    do i=1, size(fields)
      if (fields(i)%name .eq. field_name .and. fields(i)%namespace .eq. field_namespace) then
        data_type=fields(i)%data_type
        field_type=fields(i)%field_type
        return
      end if      
    end do
    data_type=0
    field_type=0
  end subroutine get_type_of_field

  subroutine provide_field_to_writer_federator_rvalues_src(io_configuration, field_name, field_namespace, &
       field_values, timestep, time, frequency, source)
    type(io_configuration_type), intent(inout) :: io_configuration
    character(len=*), intent(in) :: field_name, field_namespace
    integer, intent(in) :: timestep, frequency, source
    real(kind=DEFAULT_PRECISION), dimension(:), intent(in) :: field_values
    real(kind=DEFAULT_PRECISION), intent(in) :: time

    type(data_values_type) :: packaged_data
    
    allocate(packaged_data%values(size(field_values)))
    packaged_data%values=field_values
    packaged_data%data_type=DOUBLE_DATA_TYPE
    call provide_field_to_writer_federator_src(io_configuration, field_name, field_namespace, packaged_data, &
         timestep, time, frequency, source)
  end subroutine provide_field_to_writer_federator_rvalues_src

  subroutine provide_field_to_writer_federator_rvalues_nosrc(io_configuration, field_name, field_namespace, &
       field_values, timestep, time, frequency)
    type(io_configuration_type), intent(inout) :: io_configuration
    character(len=*), intent(in) :: field_name, field_namespace
    integer, intent(in) :: timestep, frequency
    real(kind=DEFAULT_PRECISION), dimension(:), intent(in) :: field_values
    real(kind=DEFAULT_PRECISION), intent(in) :: time

    call provide_field_to_writer_federator_rvalues_src(io_configuration, field_name, field_namespace, &
         field_values, timestep, time, frequency, -1)
  end subroutine provide_field_to_writer_federator_rvalues_nosrc

  !> Provides a field to the write federator with no source (a none collective diagnostic)
  !! @param io_configuration The IO server configuration
  !! @param field_name The field name
  !! @param field_values The fields values
  !! @param timestep The corresponding model timestep
  !! @param time Corresponding MONC model time
  !! @param frequency Configured sampling frequency
  subroutine provide_field_to_writer_federator_nosrc(io_configuration, field_name, field_namespace, &
       field_values, timestep, time, frequency)
    type(io_configuration_type), intent(inout) :: io_configuration
    character(len=*), intent(in) :: field_name, field_namespace
    integer, intent(in) :: timestep, frequency
    type(data_values_type), intent(inout) :: field_values
    real(kind=DEFAULT_PRECISION), intent(in) :: time

    call provide_field_to_writer_federator_src(io_configuration, field_name, field_namespace, &
         field_values, timestep, time, frequency, -1)
  end subroutine provide_field_to_writer_federator_nosrc
      
  !> Provides a field to the write federator (a collective diagnostic or prognostic)
  !! @param io_configuration The IO server configuration
  !! @param field_name The field name
  !! @param field_values The fields values
  !! @param timestep The corresponding model timestep
  !! @param time Corresponding MONC model time
  !! @param frequency Configured sampling frequency
  !! @param source The MONC source ID
  subroutine provide_field_to_writer_federator_src(io_configuration, field_name, field_namespace, &
       field_values, timestep, time, frequency, source, terminated_case)
    type(io_configuration_type), intent(inout) :: io_configuration
    character(len=*), intent(in) :: field_name, field_namespace
    integer, intent(in) :: timestep, frequency, source
    type(data_values_type), intent(inout) :: field_values
    real(kind=DEFAULT_PRECISION), intent(in) :: time
    logical, intent(in), optional :: terminated_case

    type(field_ordering_type), pointer :: field_ordering
    class(*), pointer :: generic
    logical :: this_is_termination, ready_to_provide

    if (present(terminated_case)) then
      this_is_termination=terminated_case
    else
      this_is_termination=.false.
    end if

    field_ordering=>get_or_add_field_ordering(field_name, field_namespace, frequency, source)

    if (time_basis) then
      ready_to_provide = (nint(time) == field_ordering%last_time_access + frequency)
    else
      ready_to_provide = (timestep == field_ordering%last_timestep_access + frequency)
    end if

    call check_thread_status(forthread_mutex_lock(field_ordering%access_mutex))
    if (ready_to_provide .or. this_is_termination) then
      if (.not. this_is_termination) then
        field_ordering%last_timestep_access=timestep
        field_ordering%last_time_access=nint(time)
      end if
      call provide_ordered_field_to_writer_federator(io_configuration, field_name, field_namespace, &
           field_values, timestep, time, source)
      if (allocated(field_values%values)) deallocate(field_values%values)
    else
      generic=>generate_value_container(field_name, field_namespace, field_values, timestep, time, frequency, source)
      if (time_basis) then
        call c_put_generic(field_ordering%timestep_to_value, conv_to_string(nint(time)), generic, .false.)
      else
        call c_put_generic(field_ordering%timestep_to_value, conv_to_string(timestep), generic, .false.)
      end if
    end if
    call process_queued_items(io_configuration, field_ordering)
    call check_thread_status(forthread_mutex_unlock(field_ordering%access_mutex))
  end subroutine provide_field_to_writer_federator_src

  !> Processes queued up items for a specific field's ordering. This will send any available fields to the write federator
  !! with a guaranteed ordering and clean up the memory associated with them
  !! @param io_configuration The IO server configuration
  !! @param field_ordering The field ordering status associated with this field
  subroutine process_queued_items(io_configuration, field_ordering)
    type(io_configuration_type), intent(inout) :: io_configuration
    type(field_ordering_type), intent(inout) :: field_ordering

    integer :: next_timestep
    type(field_ordering_value_type), pointer :: field_ordering_value_at_timestep

    do while (.not. c_is_empty(field_ordering%timestep_to_value))
      if (time_basis) then
        next_timestep=field_ordering%last_time_access + field_ordering%frequency
      else
        next_timestep=field_ordering%last_timestep_access + field_ordering%frequency
      end if
      if (c_contains(field_ordering%timestep_to_value, conv_to_string(next_timestep))) then
        field_ordering_value_at_timestep=>get_field_ordering_value_at_timestep(field_ordering%timestep_to_value, next_timestep)
        call c_remove(field_ordering%timestep_to_value, conv_to_string(next_timestep))
        field_ordering%last_timestep_access=field_ordering_value_at_timestep%timestep
        field_ordering%last_time_access=nint(field_ordering_value_at_timestep%time)
        call provide_ordered_field_to_writer_federator(io_configuration, field_ordering_value_at_timestep%field_name, &
             field_ordering_value_at_timestep%field_namespace, field_ordering_value_at_timestep%field_values, &
             field_ordering_value_at_timestep%timestep, field_ordering_value_at_timestep%time, &
             field_ordering_value_at_timestep%source)
        if (allocated(field_ordering_value_at_timestep%field_values)) then
          if (allocated(field_ordering_value_at_timestep%field_values%values)) &
               deallocate(field_ordering_value_at_timestep%field_values%values)
          deallocate(field_ordering_value_at_timestep%field_values)
        end if
        deallocate(field_ordering_value_at_timestep)
      else
        exit
      end if
    end do
  end subroutine process_queued_items

  !> Retrieves a specific field ordering value at the corresponding timestep or null if none is found
  !! @param collection The map of timesteps to field ordering values that we are looking up
  !! @param timestep The timestep to look up
  !! @returns The corresponding field ordering or null if none is found
  function get_field_ordering_value_at_timestep(collection, timestep)
    type(hashmap_type), intent(inout) :: collection
    integer, intent(in) :: timestep
    type(field_ordering_value_type), pointer :: get_field_ordering_value_at_timestep

    class(*), pointer :: generic

    generic=>c_get_generic(collection, conv_to_string(timestep))
    if (associated(generic)) then
      select type(generic)
        type is (field_ordering_value_type)
          get_field_ordering_value_at_timestep=>generic
      end select      
    else
      get_field_ordering_value_at_timestep=>null()
    end if    
  end function get_field_ordering_value_at_timestep  

  !> Generates the field value container which is then filled in with appropriate values and added into the overall field ordering
  !! data structure
  !! @param field_name The field name
  !! @param field_values The fields values
  !! @param timestep The corresponding MONC model timestep
  !! @param time The corresponding MONC model time
  !! @param frequency Configured sampling frequency
  !! @param source MONC source pid
  !! @returns A value container which can be added into the field ordering data structure
  function generate_value_container(field_name, field_namespace, field_values, timestep, time, frequency, source)
    character(len=*), intent(in) :: field_name, field_namespace
    integer, intent(in) :: timestep, frequency, source
    type(data_values_type), intent(in) :: field_values
    real(kind=DEFAULT_PRECISION), intent(in) :: time
    type(field_ordering_value_type), pointer :: generate_value_container

    allocate(generate_value_container)
    generate_value_container%field_name=field_name
    generate_value_container%field_namespace=field_namespace
    generate_value_container%timestep=timestep
    generate_value_container%frequency=frequency
    generate_value_container%time=time
    generate_value_container%source=source
    allocate(generate_value_container%field_values, source=field_values)
  end function generate_value_container  

  !> Retrieves or adds ordering for a specific field (and MONC source)
  !! @param field_name The field name
  !! @param frequency The sampling frequency of this field
  !! @param source MONC source PID
  !! @returns Either the existing or a newly created field ordering
  function get_or_add_field_ordering(field_name, field_namespace, frequency, source)
    character(len=*), intent(in) :: field_name, field_namespace
    integer, intent(in) :: frequency, source
    type(field_ordering_type), pointer :: get_or_add_field_ordering

    class(*), pointer :: generic
    character(len=STRING_LENGTH) :: entry_key

    if (source .gt. -1) then
      entry_key=trim(field_name)//"#"//trim(field_namespace)//"#"//trim(conv_to_string(source))
    else
      entry_key=field_name
    end if

    get_or_add_field_ordering=>get_field_ordering(entry_key, .true.)
    if (.not. associated(get_or_add_field_ordering)) then
      call check_thread_status(forthread_rwlock_wrlock(field_lock))
      get_or_add_field_ordering=>get_field_ordering(entry_key, .false.)
      if (.not. associated(get_or_add_field_ordering)) then
        allocate(get_or_add_field_ordering)
        get_or_add_field_ordering%last_timestep_access=0
        get_or_add_field_ordering%last_time_access = merge(                                   &
                    nint(model_initial_time - mod(real(model_initial_time),real(frequency))), &
                    0,                                                                        &
                    frequency .gt. 0 )
        get_or_add_field_ordering%frequency=frequency
        call check_thread_status(forthread_mutex_init(get_or_add_field_ordering%access_mutex, -1))
        generic=>get_or_add_field_ordering
        call c_put_generic(field_orderings, entry_key, generic, .false.)
      end if
      call check_thread_status(forthread_rwlock_unlock(field_lock))
    end if
  end function get_or_add_field_ordering
  
  !> Retrieves a field ordering based upon the name or null if none can be found
  !! @param field_name The name of the field to look up
  !! @param do_lock Whether to issue a read lock or not
  !! @returns The corresponding field ordering data structure or null if none is found
  function get_field_ordering(field_name, do_lock)
    character(len=*), intent(in) :: field_name
    logical, intent(in) :: do_lock
    type(field_ordering_type), pointer :: get_field_ordering

    class(*), pointer :: generic

    if (do_lock) call check_thread_status(forthread_rwlock_rdlock(field_lock))
    generic=>c_get_generic(field_orderings, field_name)
    if (do_lock) call check_thread_status(forthread_rwlock_unlock(field_lock))
    if (associated(generic)) then
      select type(generic)
        type is (field_ordering_type)
          get_field_ordering=>generic
      end select      
    else
      get_field_ordering=>null()
    end if
  end function get_field_ordering

  !> Prepares to serialise the field manager state, both determines storage needed and also issues any locks
  !! @returns The number of bytes needed to store the serialised state
  integer(kind=8) function prepare_to_serialise_field_manager_state()
    type(mapentry_type) :: map_entry
    type(iterator_type) :: iterator
    class(*), pointer :: generic

    call check_thread_status(forthread_rwlock_rdlock(field_lock))
    prepare_to_serialise_field_manager_state=kind(prepare_to_serialise_field_manager_state)

    iterator=c_get_iterator(field_orderings)
    do while (c_has_next(iterator))
      map_entry=c_next_mapentry(iterator)
      generic=>c_get_generic(map_entry)
      if (associated(generic)) then
        select type(generic)
        type is (field_ordering_type)
          prepare_to_serialise_field_manager_state=prepare_to_serialise_field_manager_state+&
               prepare_to_serialise_specific_field_ordering(generic)+&
               (kind(prepare_to_serialise_field_manager_state)*2)+len(trim(map_entry%key))
        end select
      end if
    end do
  end function prepare_to_serialise_field_manager_state

  !> Serialises the current field manager, releases any locks issued during preparation
  !> @param byte_data Packaged field manager serialised state.
  subroutine serialise_field_manager_state(byte_data)
    character, dimension(:), allocatable, intent(inout) :: byte_data

    integer :: current_data_point, prev_pt
    type(mapentry_type) :: map_entry
    type(iterator_type) :: iterator
    class(*), pointer :: generic
 
    current_data_point=1
    current_data_point=pack_scalar_field(byte_data, current_data_point, c_size(field_orderings))

    iterator=c_get_iterator(field_orderings)
    do while (c_has_next(iterator))
      map_entry=c_next_mapentry(iterator)
      generic=>c_get_generic(map_entry)
      if (associated(generic)) then
        select type(generic)
        type is (field_ordering_type)
          current_data_point=pack_scalar_field(byte_data, current_data_point, len(trim(map_entry%key)))
          byte_data(current_data_point:current_data_point+len(trim(map_entry%key))-1) = transfer(trim(map_entry%key), &
               byte_data(current_data_point:current_data_point+len(trim(map_entry%key))-1))
          current_data_point=current_data_point+len(trim(map_entry%key))

          prev_pt=current_data_point
          current_data_point=current_data_point+kind(current_data_point)
          call serialise_specific_field_ordering(generic, byte_data, current_data_point)          
          prev_pt=pack_scalar_field(byte_data, prev_pt, (current_data_point-kind(current_data_point))-prev_pt)
        end select
      end if
    end do    
    
    call check_thread_status(forthread_rwlock_unlock(field_lock))
  end subroutine serialise_field_manager_state

  !> Unserialses from some byte data into the state
  !! @param byte_data Serialised byte data to inflate and initialise from
  subroutine unserialise_field_manager_state(byte_data)
    character, dimension(:), intent(in) :: byte_data

    integer :: current_data_point, number_of_entries, i, byte_size, key_size
    character(len=STRING_LENGTH) :: value_key
    class(*), pointer :: generic

    current_data_point=1
    number_of_entries=unpack_scalar_integer_from_bytedata(byte_data, current_data_point)

    if (number_of_entries .gt. 0) then
      do i=1, number_of_entries
        key_size=unpack_scalar_integer_from_bytedata(byte_data, current_data_point)
        value_key=transfer(byte_data(current_data_point:current_data_point+key_size-1), value_key)
        value_key(key_size+1:)=" "
        current_data_point=current_data_point+key_size
        byte_size=unpack_scalar_integer_from_bytedata(byte_data, current_data_point)        
        generic=>unserialise_specific_field_ordering(byte_data(current_data_point:current_data_point+byte_size-1))
        call c_put_generic(field_orderings, value_key, generic, .false.)
        current_data_point=current_data_point+byte_size
      end do
    end if
  end subroutine unserialise_field_manager_state

  !> Prepares to serialise a specific field ordering, both determines storage size and issues locks
  !! @param specific_field_ordering The field ordering to prepare for serialisation
  !! @returns The number of bytes needed to store the serialised state
  integer(kind=8) function prepare_to_serialise_specific_field_ordering(specific_field_ordering)
    type(field_ordering_type), intent(inout) :: specific_field_ordering

    type(mapentry_type) :: map_entry
    type(iterator_type) :: iterator
    class(*), pointer :: generic

    call check_thread_status(forthread_mutex_lock(specific_field_ordering%access_mutex))

    prepare_to_serialise_specific_field_ordering=kind(prepare_to_serialise_specific_field_ordering) * 3

    iterator=c_get_iterator(specific_field_ordering%timestep_to_value)
    do while (c_has_next(iterator))
      map_entry=c_next_mapentry(iterator)
      generic=>c_get_generic(map_entry)
      if (associated(generic)) then
        select type(generic)
        type is (field_ordering_value_type)
          prepare_to_serialise_specific_field_ordering=prepare_to_serialise_specific_field_ordering+&
               prepare_to_serialise_field_ordering_value(generic)+(kind(prepare_to_serialise_specific_field_ordering)*2)+&
               len(trim(map_entry%key))
        end select
      end if
    end do
  end function prepare_to_serialise_specific_field_ordering  

  !> Serialises a specific fields ordering and releases any locks issued during preparation
  !! @param specific_field_ordering The field ordering to serialise
  !! @param byte_data Byte data which contains the packaged state of the field ordering
  !! @param current_data_point The current write point in the byte data, is updated during call so represents next point on return
  subroutine serialise_specific_field_ordering(specific_field_ordering, byte_data, current_data_point)
    type(field_ordering_type), intent(inout) :: specific_field_ordering
    character, dimension(:), allocatable, intent(inout) :: byte_data
    integer, intent(inout) :: current_data_point

    integer :: prev_pt
    type(mapentry_type) :: map_entry
    type(iterator_type) :: iterator
    class(*), pointer :: generic
    
    current_data_point=pack_scalar_field(byte_data, current_data_point, specific_field_ordering%last_timestep_access)
    current_data_point=pack_scalar_field(byte_data, current_data_point, specific_field_ordering%last_time_access)
    current_data_point=pack_scalar_field(byte_data, current_data_point, specific_field_ordering%frequency)
    current_data_point=pack_scalar_field(byte_data, current_data_point, c_size(specific_field_ordering%timestep_to_value))

    iterator=c_get_iterator(specific_field_ordering%timestep_to_value)
    do while (c_has_next(iterator))
      map_entry=c_next_mapentry(iterator)
      generic=>c_get_generic(map_entry)
      if (associated(generic)) then
        select type(generic)
        type is (field_ordering_value_type)
          current_data_point=pack_scalar_field(byte_data, current_data_point, len(trim(map_entry%key)))
          byte_data(current_data_point:current_data_point+len(trim(map_entry%key))-1) = transfer(trim(map_entry%key), &
               byte_data(current_data_point:current_data_point+len(trim(map_entry%key))-1))
          current_data_point=current_data_point+len(trim(map_entry%key))

          prev_pt=current_data_point
          current_data_point=current_data_point+kind(current_data_point)
          call serialise_field_ordering_value(generic, byte_data, current_data_point)          
          prev_pt=pack_scalar_field(byte_data, prev_pt, (current_data_point-kind(current_data_point))-prev_pt)
        end select
      end if
    end do
    
    call check_thread_status(forthread_mutex_unlock(specific_field_ordering%access_mutex))
  end subroutine serialise_specific_field_ordering
  
  !> Unserialises some field ordering
  !! @param byte_data The serialised representation which we work from
  !! @returns The initialised field ordering represented by the provided byte data
  function unserialise_specific_field_ordering(byte_data)
    character, dimension(:), intent(in) :: byte_data
    type(field_ordering_type), pointer :: unserialise_specific_field_ordering

    integer :: current_data_point, number_of_values, byte_size, i, key_size
    character(len=STRING_LENGTH) :: value_key
    class(*), pointer :: generic

    allocate(unserialise_specific_field_ordering)
    
    current_data_point=1
    unserialise_specific_field_ordering%last_timestep_access=unpack_scalar_integer_from_bytedata(byte_data, current_data_point)
    unserialise_specific_field_ordering%last_time_access=unpack_scalar_integer_from_bytedata(byte_data, current_data_point)
    unserialise_specific_field_ordering%frequency=unpack_scalar_integer_from_bytedata(byte_data, current_data_point)
    number_of_values=unpack_scalar_integer_from_bytedata(byte_data, current_data_point)

    if (number_of_values .gt. 0) then
      do i=1, number_of_values
        key_size=unpack_scalar_integer_from_bytedata(byte_data, current_data_point)
        value_key=transfer(byte_data(current_data_point:current_data_point+key_size-1), value_key)
        value_key(key_size+1:)=" "
        current_data_point=current_data_point+key_size
        byte_size=unpack_scalar_integer_from_bytedata(byte_data, current_data_point)
        generic=>unserialise_field_ordering_value(byte_data(current_data_point:current_data_point+byte_size-1))
        call c_put_generic(unserialise_specific_field_ordering%timestep_to_value, value_key, generic, .false.)
        current_data_point=current_data_point+byte_size
      end do
    end if
    call check_thread_status(forthread_mutex_init(unserialise_specific_field_ordering%access_mutex, -1))
  end function unserialise_specific_field_ordering

  !> Prepares to serialise a specific field ordering value, determines both the storage size and issues locks
  !! @param specific_field_value The specific value to prepare for serialisation
  !! @returns The number of bytes needed to store the serialised state
  integer(kind=8) function prepare_to_serialise_field_ordering_value(specific_field_value)
    type(field_ordering_value_type), intent(inout) :: specific_field_value

    prepare_to_serialise_field_ordering_value=prepare_to_serialise_data_values_type(specific_field_value%field_values)
    prepare_to_serialise_field_ordering_value=prepare_to_serialise_field_ordering_value+&
         (kind(specific_field_value%timestep) * 6) + len(trim(specific_field_value%field_name)) + &
         len(trim(specific_field_value%field_namespace)) + kind(specific_field_value%time)
  end function prepare_to_serialise_field_ordering_value

  !> Serialises a field ordering value and releases any issued locks
  !! @param specific_field_value The specific value to serialise
  !! @param byte_data Target byte data whihc is allocated here
  !! @param current_data_point The current write point in the byte data, is updated during call so represents next point on return
  subroutine serialise_field_ordering_value(specific_field_value, byte_data, current_data_point)
    type(field_ordering_value_type), intent(inout) :: specific_field_value
    character, dimension(:), allocatable, intent(inout) :: byte_data
    integer, intent(inout) :: current_data_point

    integer :: prev_pt    
    
    current_data_point=pack_scalar_field(byte_data, current_data_point, specific_field_value%timestep)
    current_data_point=pack_scalar_field(byte_data, current_data_point, specific_field_value%frequency)
    current_data_point=pack_scalar_field(byte_data, current_data_point, specific_field_value%source)
    current_data_point=pack_scalar_field(byte_data, current_data_point, len(trim(specific_field_value%field_name)))
    byte_data(current_data_point:current_data_point+len(trim(specific_field_value%field_name))-1) = &
         transfer(trim(specific_field_value%field_name), byte_data(current_data_point:current_data_point+&
         len(trim(specific_field_value%field_name))-1))
    current_data_point=current_data_point+len(trim(specific_field_value%field_name))
    current_data_point=pack_scalar_field(byte_data, current_data_point, len(trim(specific_field_value%field_namespace)))
    byte_data(current_data_point:current_data_point+len(trim(specific_field_value%field_namespace))-1) = &
         transfer(trim(specific_field_value%field_namespace), byte_data(current_data_point:current_data_point+&
         len(trim(specific_field_value%field_namespace))-1))
    current_data_point=current_data_point+len(trim(specific_field_value%field_namespace))
    current_data_point=pack_scalar_field(byte_data, current_data_point, double_real_value=specific_field_value%time)

    prev_pt=current_data_point
    current_data_point=current_data_point+kind(current_data_point)
    call serialise_data_values_type(specific_field_value%field_values, byte_data, current_data_point)    
    prev_pt=pack_scalar_field(byte_data, prev_pt, (current_data_point-kind(current_data_point))-prev_pt)
  end subroutine serialise_field_ordering_value

  !> Unseralises some field ordering from its byte representation
  !! @param byte_data The serialised representation
  !! @returns The corresponding initialised and completed ordering value
  function unserialise_field_ordering_value(byte_data)
    character, dimension(:), intent(in) :: byte_data
    type(field_ordering_value_type), pointer :: unserialise_field_ordering_value

    integer :: current_data_point, byte_size, str_size

    allocate(unserialise_field_ordering_value)
    current_data_point=1
    unserialise_field_ordering_value%timestep=unpack_scalar_integer_from_bytedata(byte_data, current_data_point)
    unserialise_field_ordering_value%frequency=unpack_scalar_integer_from_bytedata(byte_data, current_data_point)
    unserialise_field_ordering_value%source=unpack_scalar_integer_from_bytedata(byte_data, current_data_point)
    str_size=unpack_scalar_integer_from_bytedata(byte_data, current_data_point)
    unserialise_field_ordering_value%field_name=transfer(byte_data(current_data_point:current_data_point+str_size-1), &
         unserialise_field_ordering_value%field_name)
    unserialise_field_ordering_value%field_name(str_size+1:)=" "
    current_data_point=current_data_point+str_size
    str_size=unpack_scalar_integer_from_bytedata(byte_data, current_data_point)
    unserialise_field_ordering_value%field_namespace=transfer(byte_data(current_data_point:current_data_point+str_size-1), &
         unserialise_field_ordering_value%field_namespace)
    unserialise_field_ordering_value%field_namespace(str_size+1:)=" "
    current_data_point=current_data_point+str_size
    unserialise_field_ordering_value%time=unpack_scalar_dp_real_from_bytedata(byte_data, current_data_point)
    byte_size=unpack_scalar_integer_from_bytedata(byte_data, current_data_point)        
    unserialise_field_ordering_value%field_values=unserialise_data_values_type(byte_data(current_data_point:&
         current_data_point+byte_size-1))
  end function unserialise_field_ordering_value
end module writer_field_manager_mod
