!> Performs time averaged, time manipulation and only returns a value if the output frequency determines one should be
module timeaveraged_time_manipulation_mod
  use datadefn_mod, only : DEFAULT_PRECISION, STRING_LENGTH
  use configuration_parser_mod, only : io_configuration_type
  use collections_mod, only : mapentry_type, iterator_type, hashmap_type, c_get_generic, c_size, &
       c_next_mapentry, c_has_next, c_get_iterator, c_put_generic
  use forthread_mod, only : forthread_rwlock_rdlock, forthread_rwlock_wrlock, forthread_rwlock_unlock, &
       forthread_rwlock_init, forthread_rwlock_destroy, forthread_mutex_init, forthread_mutex_lock, &
       forthread_mutex_unlock, forthread_mutex_destroy
  use threadpool_mod, only : check_thread_status
  use configuration_parser_mod, only : data_values_type
  use io_server_client_mod, only : pack_scalar_field, pack_array_field
  use data_utils_mod, only : unpack_scalar_integer_from_bytedata, unpack_scalar_dp_real_from_bytedata, &
       unpack_scalar_logical_from_bytedata

  implicit none

  real(kind=DEFAULT_PRECISION) :: model_initial_time

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

  public init_time_averaged_manipulation, finalise_time_averaged_manipulation, perform_timeaveraged_time_manipulation, &
       is_time_averaged_time_manipulation_ready_to_write, serialise_time_averaged_state, unserialise_time_averaged_state, &
       prepare_to_serialise_time_averaged_state
contains  

   !> Initialises the reduction action
  subroutine init_time_averaged_manipulation(reconfig_initial_time)
    real(kind=DEFAULT_PRECISION), intent(in) :: reconfig_initial_time
    
    model_initial_time = reconfig_initial_time
    call check_thread_status(forthread_rwlock_init(timeaveraged_value_rw_lock, -1))
  end subroutine init_time_averaged_manipulation  

  !> Finalises the reduction action, waiting for all outstanding requests and then freeing data
  !! @param io_configuration Configuration state of the IO server
  subroutine finalise_time_averaged_manipulation()
    call check_thread_status(forthread_rwlock_destroy(timeaveraged_value_rw_lock))
  end subroutine finalise_time_averaged_manipulation

  logical function is_time_averaged_time_manipulation_ready_to_write(latest_time, output_frequency, write_time, &
       latest_timestep, write_timestep)
    real, intent(in) :: latest_time, output_frequency, write_time
    integer, intent(in) :: latest_timestep, write_timestep

    is_time_averaged_time_manipulation_ready_to_write=latest_time + output_frequency .ge. write_time
  end function is_time_averaged_time_manipulation_ready_to_write

  !> Performs the time averaged manipulation and only returns values if these are to be stored (i.e. past an output frequency)
  !! @param instant_values The instantaneous values to work with
  !! @param output_frequency The output frequency configuration option
  !! @param field_name The field name
  !! @param timestep The timestep
  !! @param time The model time
  !! @param time_basis True for diagnostics interval in time coordinates, False for timestep coordinates
  !! @returns An allocated array of reals if data is to be stored, otherwise this is unallocated
  type(data_values_type) function perform_timeaveraged_time_manipulation(instant_values, output_frequency, &
       field_name, timestep, time, time_basis)
    real(kind=DEFAULT_PRECISION), dimension(:), intent(in) :: instant_values
    real, intent(in) :: output_frequency
    real(kind=DEFAULT_PRECISION), intent(in) :: time
    character(len=*), intent(in) :: field_name
    integer, intent(in) :: timestep
    logical, intent(in) :: time_basis

    type(time_averaged_completed_type), pointer :: timeaveraged_value
    logical :: select_value

    timeaveraged_value=>find_or_add_timeaveraged_value(timestep, field_name)

    call check_thread_status(forthread_mutex_lock(timeaveraged_value%mutex))
    call time_average(timeaveraged_value, instant_values, time)

    ! time_basis requires regular-interval entries.  Timestep requires time .ge. time+previous_output_time
    if (time_basis) then
      select_value = mod(nint(time), nint(output_frequency)) == 0
    else
      ! these are 'reals' to be consistent with other such calculations elsewhere
      select_value = real(time)-real(timeaveraged_value%previous_output_time) .ge. real(output_frequency)
    end if

    if (select_value) then
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
    real(kind=DEFAULT_PRECISION), intent(in) :: time

    integer :: i
    real(kind=DEFAULT_PRECISION) :: timeav, timedg, combined_add    

    timedg = time - timeaveraged_value%previous_time
    timeav = time - timeaveraged_value%start_time - timedg

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

  !> Prepares to serialise the time averaged state values. Both determines the storage size required and also issue locks
  !! @returns The number of bytes needed to store the serialised state
  integer(kind=8) function prepare_to_serialise_time_averaged_state()
    type(mapentry_type) :: map_entry
    type(iterator_type) :: iterator
    class(*), pointer :: generic
    
    call check_thread_status(forthread_rwlock_rdlock(timeaveraged_value_rw_lock))
    
    prepare_to_serialise_time_averaged_state=kind(prepare_to_serialise_time_averaged_state)
    iterator=c_get_iterator(timeaveraged_values)
    do while (c_has_next(iterator))
       map_entry=c_next_mapentry(iterator)
       generic=>c_get_generic(map_entry)
      if (associated(generic)) then
        select type(generic)
        type is (time_averaged_completed_type)
          prepare_to_serialise_time_averaged_state=prepare_to_serialise_time_averaged_state+&
               prepare_to_serialise_time_averaged_completed_value(generic)+&
               (kind(prepare_to_serialise_time_averaged_state)*2)+len(trim(map_entry%key))
        end select
      end if      
    end do
  end function prepare_to_serialise_time_averaged_state

  !> Serialises the state of this manipulator so that it can be restarted later on. Releases any locks issue during preparation.
  !! @param byte_data The byte data that represents the serialised state
  subroutine serialise_time_averaged_state(byte_data)
    character, dimension(:), allocatable, intent(inout) :: byte_data

    integer :: current_data_point, prev_pt
    type(mapentry_type) :: map_entry
    type(iterator_type) :: iterator
    class(*), pointer :: generic
    
    current_data_point=1
    current_data_point=pack_scalar_field(byte_data, current_data_point, c_size(timeaveraged_values))

    iterator=c_get_iterator(timeaveraged_values)
    do while (c_has_next(iterator))
       map_entry=c_next_mapentry(iterator)
       generic=>c_get_generic(map_entry)
      if (associated(generic)) then
        select type(generic)
        type is (time_averaged_completed_type)
          current_data_point=pack_scalar_field(byte_data, current_data_point, len(trim(map_entry%key)))
          byte_data(current_data_point:current_data_point+len(trim(map_entry%key))-1) = transfer(trim(map_entry%key), &
               byte_data(current_data_point:current_data_point+len(trim(map_entry%key))-1))
          current_data_point=current_data_point+len(trim(map_entry%key))

          prev_pt=current_data_point
          current_data_point=current_data_point+kind(current_data_point)
          call serialise_time_averaged_completed_value(generic, byte_data, current_data_point)          
          prev_pt=pack_scalar_field(byte_data, prev_pt, (current_data_point-kind(current_data_point)) - prev_pt)
        end select
      end if      
    end do
    call check_thread_status(forthread_rwlock_unlock(timeaveraged_value_rw_lock))
  end subroutine serialise_time_averaged_state

  !> Unserialises some byte data to initialise the state from some previous version
  !! @param byte_data The byte data to read from and initialise from
  subroutine unserialise_time_averaged_state(byte_data)
    character, dimension(:), intent(in) :: byte_data

    integer :: current_data_point, number_entries, i, key_size, byte_size
    character(len=STRING_LENGTH) :: value_key
    class(*), pointer :: generic

    current_data_point=1
    number_entries=unpack_scalar_integer_from_bytedata(byte_data, current_data_point)
    if (number_entries .gt. 0) then
      do i=1, number_entries
        key_size=unpack_scalar_integer_from_bytedata(byte_data, current_data_point)
        value_key=transfer(byte_data(current_data_point:current_data_point+key_size-1), value_key)
        value_key(key_size+1:)=" "
        current_data_point=current_data_point+key_size
        byte_size=unpack_scalar_integer_from_bytedata(byte_data, current_data_point)
        generic=>unserialise_time_averaged_completed_value(byte_data(current_data_point:current_data_point+byte_size-1))
        call c_put_generic(timeaveraged_values, value_key, generic, .false.)
        current_data_point=current_data_point+byte_size
      end do
    end if
  end subroutine unserialise_time_averaged_state

  !> Prepares to serialise a time averaged completed value, both determines the storage size and also issue any locks
  !! @param time_av_value The time averaged completed value to prepare for serialisation
  !! @returns The number of bytes needed to store the serialised state
  integer(kind=8) function prepare_to_serialise_time_averaged_completed_value(time_av_value)
    type(time_averaged_completed_type), intent(inout) :: time_av_value

    call check_thread_status(forthread_mutex_lock(time_av_value%mutex))

    prepare_to_serialise_time_averaged_completed_value=(kind(time_av_value%start_time) * 3) + kind(time_av_value%empty_values) + &
         (size(time_av_value%time_averaged_values) * kind(time_av_value%time_averaged_values)) + &
         (kind(prepare_to_serialise_time_averaged_completed_value) * 2) + len(time_av_value%field_name)
  end function prepare_to_serialise_time_averaged_completed_value

  !> Serialises a specific time averaged completed value, releases any locks issued during preparation
  !! @param time_av_value The time averaged completed value to serialise
  !! @param byte_data Resulting byte data that holds a representation of this
  !! @param current_data_point The current write point in the byte data, is updated during call so represents next point on return
  subroutine serialise_time_averaged_completed_value(time_av_value, byte_data, current_data_point)
    type(time_averaged_completed_type), intent(inout) :: time_av_value
    character, dimension(:), allocatable, intent(inout) :: byte_data
    integer, intent(inout) :: current_data_point

    integer :: i

    current_data_point=pack_scalar_field(byte_data, current_data_point, double_real_value=time_av_value%start_time)
    current_data_point=pack_scalar_field(byte_data, current_data_point, double_real_value=time_av_value%previous_time)
    current_data_point=pack_scalar_field(byte_data, current_data_point, double_real_value=time_av_value%previous_output_time)
    current_data_point=pack_scalar_field(byte_data, current_data_point, logical_value=time_av_value%empty_values)
    current_data_point=pack_scalar_field(byte_data, current_data_point, len(trim(time_av_value%field_name)))
    byte_data(current_data_point:current_data_point+len(trim(time_av_value%field_name))-1) = transfer(&
         trim(time_av_value%field_name), byte_data(current_data_point:current_data_point+len(trim(time_av_value%field_name))-1))
    current_data_point=current_data_point+len(trim(time_av_value%field_name))
    current_data_point=pack_scalar_field(byte_data, current_data_point, size(time_av_value%time_averaged_values))
    current_data_point=pack_array_field(byte_data, current_data_point, real_array_1d=time_av_value%time_averaged_values)
    call check_thread_status(forthread_mutex_unlock(time_av_value%mutex))
  end subroutine serialise_time_averaged_completed_value

  !> Will create a specific time averaged completed value based upon the provided serialised data
  !! @param byte_data The serialised byte data to read and initialise from
  function unserialise_time_averaged_completed_value(byte_data)
    character, dimension(:), intent(in) :: byte_data
    type(time_averaged_completed_type), pointer :: unserialise_time_averaged_completed_value

    integer :: current_data_point, i, values_size, byte_size, str_size

    allocate(unserialise_time_averaged_completed_value)
    current_data_point=1
    unserialise_time_averaged_completed_value%start_time=unpack_scalar_dp_real_from_bytedata(byte_data, current_data_point)
    unserialise_time_averaged_completed_value%previous_time=unpack_scalar_dp_real_from_bytedata(byte_data, current_data_point)
    unserialise_time_averaged_completed_value%previous_output_time=&
         unpack_scalar_dp_real_from_bytedata(byte_data, current_data_point)
    unserialise_time_averaged_completed_value%empty_values=unpack_scalar_logical_from_bytedata(byte_data, current_data_point)
    str_size=unpack_scalar_integer_from_bytedata(byte_data, current_data_point)
    unserialise_time_averaged_completed_value%field_name=&
         transfer(byte_data(current_data_point:current_data_point+str_size-1), &
         unserialise_time_averaged_completed_value%field_name)
    unserialise_time_averaged_completed_value%field_name(str_size+1:)=" "
    current_data_point=current_data_point+str_size
    values_size=unpack_scalar_integer_from_bytedata(byte_data, current_data_point)
    allocate(unserialise_time_averaged_completed_value%time_averaged_values(values_size))
    byte_size=values_size*kind(unserialise_time_averaged_completed_value%time_averaged_values)
    unserialise_time_averaged_completed_value%time_averaged_values=transfer(byte_data(current_data_point:&
         current_data_point+byte_size-1), unserialise_time_averaged_completed_value%time_averaged_values)
    call check_thread_status(forthread_mutex_init(unserialise_time_averaged_completed_value%mutex, -1))
  end function unserialise_time_averaged_completed_value

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
        new_entry%start_time=model_initial_time
        new_entry%previous_time=model_initial_time
        new_entry%empty_values=.true.
        new_entry%previous_output_time=model_initial_time
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
