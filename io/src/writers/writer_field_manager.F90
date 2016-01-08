!> The writer field manager will manage aspects of the fields being provided to the writer federator. There are two major
!! aspects to this, firstly extraction of fields from the data sent directly from MONC and then passed on if these
!! fields are needed. Secondly, this will queue up the data to ensure that it is sent to the writer federator in strict
!! order and that only one piece of a fields data is in the writer federator at any one time
module writer_field_manager_mod
  use datadefn_mod, only : DEFAULT_PRECISION, STRING_LENGTH
  use configuration_parser_mod, only : io_configuration_type, io_configuration_field_type
  use collections_mod, only : hashmap_type, c_contains, c_get_generic, c_put_generic, c_is_empty, c_remove
  use conversions_mod, only : conv_to_string, conv_to_real
  use forthread_mod, only : forthread_mutex_init, forthread_mutex_lock, forthread_mutex_unlock, forthread_mutex_destroy, &
       forthread_rwlock_rdlock, forthread_rwlock_wrlock, forthread_rwlock_unlock, forthread_rwlock_init, forthread_rwlock_destroy
  use threadpool_mod, only : check_thread_status
  use data_utils_mod, only : get_scalar_integer_from_monc, get_scalar_real_from_monc, is_field_present, &
       get_array_double_from_monc, get_array_integer_from_monc
  use io_server_client_mod, only : DOUBLE_DATA_TYPE, INTEGER_DATA_TYPE
  use logging_mod, only : LOG_WARN, LOG_ERROR, log_log
  use writer_federator_mod, only : is_field_used_by_writer_federator, provide_ordered_field_to_writer_federator, &
       is_field_split_on_q
  implicit none

#ifndef TEST_MODE
  private
#endif

  type field_value_type
     character(len=STRING_LENGTH) :: field_name
     integer :: timestep, frequency, source
     real(kind=DEFAULT_PRECISION), dimension(:), allocatable :: field_values
     real(kind=DEFAULT_PRECISION) :: time
  end type field_value_type  

  type field_ordering_type
     type(hashmap_type) :: timestep_to_value
     integer :: access_mutex, last_timestep_access, frequency
  end type field_ordering_type  

  interface provide_field_to_writer_federator
     module procedure provide_field_to_writer_federator_src, provide_field_to_writer_federator_nosrc
  end interface provide_field_to_writer_federator

  integer :: field_lock
  type(hashmap_type) :: field_orderings

  public initialise_writer_field_manager, finalise_writer_field_manager, provide_monc_data_to_writer_federator, &
       provide_field_to_writer_federator
contains

  !> Initialises the writer field manager
  !! @param io_configuration Configuration of the IO server
  subroutine initialise_writer_field_manager(io_configuration)
    type(io_configuration_type), intent(inout) :: io_configuration

    call check_thread_status(forthread_rwlock_init(field_lock, -1))
  end subroutine initialise_writer_field_manager

  !> Finalises the writer field manager
  subroutine finalise_writer_field_manager()
    call check_thread_status(forthread_rwlock_destroy(field_lock))
  end subroutine finalise_writer_field_manager

  !> Data communicated from MONC is provided to this write federator and then included if the configuration has selected
  !! to use that MONC field (as opposed to a field produced by the diagnostics federator)
  !! @param io_configuration Configuration of the IO server
  !! @param source The source PID of the MONC process
  !! @param data_id The ID of the data definition that is represented by the dump
  !! @param data_dump The data we have received from MONC
  subroutine provide_monc_data_to_writer_federator(io_configuration, source, data_id, data_dump)
    type(io_configuration_type), intent(inout) :: io_configuration
    integer, intent(in) :: source, data_id
    character, dimension(:), allocatable :: data_dump

    integer :: i, num_fields, timestep
    character(len=STRING_LENGTH) :: field_name
    real(kind=DEFAULT_PRECISION) :: time
    real(kind=DEFAULT_PRECISION), dimension(:), allocatable :: monc_value

    if (is_field_present(io_configuration, source, data_id, "timestep") .and. &
         is_field_present(io_configuration, source, data_id, "time")) then
      timestep=get_scalar_integer_from_monc(io_configuration, source, data_id, data_dump, "timestep")
      time=get_scalar_real_from_monc(io_configuration, source, data_id, data_dump, "time")

      num_fields=io_configuration%data_definitions(data_id)%number_of_data_fields

      do i=1, num_fields
        field_name=io_configuration%data_definitions(data_id)%fields(i)%name
        if (is_field_present(io_configuration, source, data_id, field_name) .and. &
             (is_field_used_by_writer_federator(field_name) .or. is_field_split_on_q(field_name))) then
          monc_value=get_value_from_monc_data(io_configuration, source, data_id, data_dump, field_name)
          call provide_field_to_writer_federator_src(io_configuration, field_name, monc_value, timestep, time, &
               io_configuration%data_definitions(data_id)%frequency, source)
          deallocate(monc_value)
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
  function get_value_from_monc_data(io_configuration, source, data_id, data_dump, field_name)
    type(io_configuration_type), intent(inout) :: io_configuration
    integer, intent(in) :: source, data_id
    character, dimension(:), allocatable :: data_dump
    character(len=*), intent(in) :: field_name
    real(kind=DEFAULT_PRECISION), dimension(:), allocatable :: get_value_from_monc_data

    integer :: field_data_type, i
    integer, dimension(:), allocatable :: int_values

    field_data_type=get_datatype_of_field(io_configuration%data_definitions(data_id)%fields, field_name)
    if (field_data_type == 0) then
      call log_log(LOG_ERROR, "No data type for field '"//trim(field_name)//"'")
    end if

    if (field_data_type == DOUBLE_DATA_TYPE) then
      get_value_from_monc_data=get_array_double_from_monc(io_configuration, source, data_id, data_dump, field_name)
    else if (field_data_type == INTEGER_DATA_TYPE) then
      int_values=get_array_integer_from_monc(io_configuration, source, data_id, data_dump, field_name)
      allocate(get_value_from_monc_data(size(int_values)))
      do i=1, size(int_values)
        get_value_from_monc_data(i)=conv_to_real(int_values(i))
      end do
      deallocate(int_values)
    end if
  end function get_value_from_monc_data

  !> Retrieves the data type of a field or 0 if the field was not found
  !! @param fields Array of fields to search
  !! @param field_name The name of the field to locate
  !! @returns The data type of this field or 0 if the field was not found
  integer function get_datatype_of_field(fields, field_name)
    type(io_configuration_field_type), dimension(:), intent(in) :: fields
    character(len=*), intent(in) :: field_name

    integer :: i

    do i=1, size(fields)
      if (fields(i)%name .eq. field_name) then
        get_datatype_of_field=fields(i)%data_type
        return
      end if      
    end do
    get_datatype_of_field=0
  end function get_datatype_of_field

  !> Provides a field to the write federator with no source (a none collective diagnostic)
  !! @param io_configuration The IO server configuration
  !! @param field_name The field name
  !! @param field_values The fields values
  !! @param timestep The corresponding model timestep
  !! @param time Corresponding MONC model time
  !! @param frequency Configured sampling frequency
  subroutine provide_field_to_writer_federator_nosrc(io_configuration, field_name, field_values, timestep, time, frequency)
    type(io_configuration_type), intent(inout) :: io_configuration
    character(len=*), intent(in) :: field_name
    integer, intent(in) :: timestep, frequency
    real(kind=DEFAULT_PRECISION), dimension(:), intent(in) :: field_values
    real(kind=DEFAULT_PRECISION), intent(in) :: time

    call provide_field_to_writer_federator_src(io_configuration, field_name, field_values, timestep, time, frequency, -1)
  end subroutine provide_field_to_writer_federator_nosrc
      
  !> Provides a field to the write federator (a collective diagnostic or prognostic)
  !! @param io_configuration The IO server configuration
  !! @param field_name The field name
  !! @param field_values The fields values
  !! @param timestep The corresponding model timestep
  !! @param time Corresponding MONC model time
  !! @param frequency Configured sampling frequency
  !! @param source The MONC source ID
  subroutine provide_field_to_writer_federator_src(io_configuration, field_name, field_values, timestep, time, frequency, source)
    type(io_configuration_type), intent(inout) :: io_configuration
    character(len=*), intent(in) :: field_name
    integer, intent(in) :: timestep, frequency, source
    real(kind=DEFAULT_PRECISION), dimension(:), intent(in) :: field_values
    real(kind=DEFAULT_PRECISION), intent(in) :: time

    type(field_ordering_type), pointer :: field_ordering
    class(*), pointer :: generic

    field_ordering=>get_or_add_field_ordering(field_name, frequency, source)
    call check_thread_status(forthread_mutex_lock(field_ordering%access_mutex))
    if (timestep == field_ordering%last_timestep_access + frequency) then
      call provide_ordered_field_to_writer_federator(io_configuration, field_name, field_values, timestep, time, source)
      field_ordering%last_timestep_access=timestep
    else
      generic=>generate_value_container(field_name, field_values, timestep, time, frequency, source)
      call c_put_generic(field_ordering%timestep_to_value, conv_to_string(timestep), generic, .false.)
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
    type(field_value_type), pointer :: field_ordering_value_at_timestep

    do while (.not. c_is_empty(field_ordering%timestep_to_value))
      next_timestep=field_ordering%last_timestep_access + field_ordering%frequency
      if (c_contains(field_ordering%timestep_to_value, conv_to_string(next_timestep))) then
        field_ordering_value_at_timestep=>get_field_ordering_value_at_timestep(field_ordering%timestep_to_value, next_timestep)
        call c_remove(field_ordering%timestep_to_value, conv_to_string(next_timestep))
        field_ordering%last_timestep_access=next_timestep
        call provide_ordered_field_to_writer_federator(io_configuration, field_ordering_value_at_timestep%field_name, &
             field_ordering_value_at_timestep%field_values, field_ordering_value_at_timestep%timestep, &
             field_ordering_value_at_timestep%time, field_ordering_value_at_timestep%source)
        if (allocated(field_ordering_value_at_timestep%field_values)) deallocate(field_ordering_value_at_timestep%field_values)
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
    type(field_value_type), pointer :: get_field_ordering_value_at_timestep

    class(*), pointer :: generic

    generic=>c_get_generic(collection, conv_to_string(timestep))
    if (associated(generic)) then
      select type(generic)
        type is (field_value_type)
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
  function generate_value_container(field_name, field_values, timestep, time, frequency, source)
    character(len=*), intent(in) :: field_name
    integer, intent(in) :: timestep, frequency, source
    real(kind=DEFAULT_PRECISION), dimension(:), intent(in) :: field_values
    real(kind=DEFAULT_PRECISION), intent(in) :: time
    type(field_value_type), pointer :: generate_value_container

    allocate(generate_value_container)
    generate_value_container%field_name=field_name
    generate_value_container%timestep=timestep
    generate_value_container%frequency=frequency
    generate_value_container%time=time
    generate_value_container%source=source
    allocate(generate_value_container%field_values(size(field_values)), source=field_values)
  end function generate_value_container  

  !> Retrieves or adds ordering for a specific field (and MONC source)
  !! @param field_name The field name
  !! @param frequency The sampling frequency of this field
  !! @param source MONC source PID
  !! @returns Either the existing or a newly created field ordering
  function get_or_add_field_ordering(field_name, frequency, source)
    character(len=*), intent(in) :: field_name
    integer, intent(in) :: frequency, source
    type(field_ordering_type), pointer :: get_or_add_field_ordering

    class(*), pointer :: generic
    character(len=STRING_LENGTH) :: entry_key

    if (source .gt. -1) then
      entry_key=trim(field_name)//"#"//trim(conv_to_string(source))
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
end module writer_field_manager_mod
