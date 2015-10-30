!> Manages the options database. Contains administration functions and deduce runtime options from the command line.
!!
!! The key-value character length limit for each option is 64 characters, the
!! value length limit is 44 characters and key length 20 characters.
!! Note that the options database should be entirely agnostic of where or now the database is stored (in our
!! case in the state.)
module optionsdatabase_mod
  use datadefn_mod, only : DEFAULT_PRECISION, STRING_LENGTH, LONG_STRING_LENGTH
  use collections_mod, only : list_type, map_type, c_size, c_get, c_contains, c_value_at, c_key_at, c_put, c_add, c_remove
  use conversions_mod, only : conv_to_generic, conv_to_logical, conv_to_integer, conv_to_real, conv_is_logical, conv_is_integer, &
       conv_is_real, conv_to_string
  use logging_mod, only: LOG_ERROR, log_log
  implicit none

#ifndef TEST_MODE
  private
#endif
  integer, parameter :: LOGICAL_TYPE=0,&  !< Type of logical value data
                        INTEGER_TYPE=1,&  !< Type of integer value data
                        REAL_TYPE=2,&     !< Type of real value data
                        STRING_TYPE=3     !< Type of string value data

  ! Extra key size for a specific array element, the lookup key is allocated with the length of the key plus this
  integer, parameter :: ARRAY_APPEND_SIZE = 10

  !> Generic add interface for adding different types of data to the databases
  interface options_add
    module procedure options_add_integer, options_add_real, options_add_logical, options_add_string
  end interface options_add

  public load_command_line_into_options_database, options_has_key, options_get_logical, options_get_integer, &
    options_get_string, options_get_real, options_add, options_size, options_key_at, options_value_at, &
    options_get_array_size, options_get_integer_array, options_get_real_array, &
    options_get_string_array, options_get_logical_array, options_remove_key

  contains

  !> Returns the number of entries in the options database
  !! @param options_database The options database
  !! @returns The number of entries
  integer function options_size(options_database)
    type(map_type), intent(inout) :: options_database

    options_size = c_size(options_database)
  end function options_size

  !> Returns the ith key in the options database
  !! @param options_database The options database
  !! @param i The index to retrieve the key at
  !! @returns The key at index i
  character(len=STRING_LENGTH) function options_key_at(options_database, i)
    type(map_type), intent(inout) :: options_database
    integer, intent(in) :: i

    options_key_at = c_key_at(options_database, i)
  end function options_key_at

  !> Returns the value at index in the database
  !! @param options_database The options database
  !! @param i The index to retrieve the value at
  !! @returns The value at index i
  function options_value_at(options_database, i)
    type(map_type), intent(inout) :: options_database
    integer, intent(in) :: i
    class(*), pointer :: options_value_at

    options_value_at => c_value_at(options_database, i)
  end function options_value_at

  !> Determines whether a specific key is in the database
  !! @param options_database The options database
  !! @param key The key to find
  !! @returns Whether they key is in the database or not
  logical function options_has_key(options_database, key)
    type(map_type), intent(inout) :: options_database
    character(len=*), intent(in) :: key

    options_has_key = c_contains(options_database, trim(key))
    if (.not. options_has_key) then
      options_has_key = c_contains(options_database, trim(key)//"a_size")
    end if
  end function options_has_key

  !> Retrieves a real value from the database that matches the provided key
  !! @param options_database The options database
  !! @param key The key to search for and return the matching value
  !! @param index Optional array index to look up an array value
  !! @returns The matching real value
  real(kind=DEFAULT_PRECISION) function options_get_real(options_database, key, index)
    type(map_type), intent(inout) :: options_database
    character(len=*), intent(in) :: key
    integer, intent(in), optional :: index

    class(*), pointer :: raw_data
    character(len=len(key)+ARRAY_APPEND_SIZE) :: lookup_key

    if (present(index)) then
      lookup_key=get_options_array_key(key, index)
    else
      lookup_key=key
    end if
    call check_options_key_exists(options_database, lookup_key)
    raw_data=>c_get(options_database, lookup_key)
    if (associated(raw_data)) options_get_real = conv_to_real(raw_data, .false.)
  end function options_get_real

  !> Retrieves an entire (or subset) real array
  !! @param options_database The options database
  !! @param key The key to search for
  !! @param array_data The array data to write into
  !! @param from Optional starting index
  !! @param to Optional end index
  subroutine options_get_real_array(options_database, key, array_data, from, to)
    type(map_type), intent(inout) :: options_database
    character(len=*), intent(in) :: key
    real(kind=DEFAULT_PRECISION), dimension(:), intent(inout) :: array_data
    integer, intent(in), optional :: from, to

    integer :: num_elements, i, start, end

    num_elements=options_get_array_size(options_database, key)
    if (num_elements .gt. 0 .and. options_has_key(options_database, trim(key)//"a_size")) then
      if (present(from)) then
        start=from
      else
        start=1
      end if
      if (present(to)) then
        end=to
      else
        end=num_elements
      end if
      do i=start, num_elements
        array_data((i-start)+1)=options_get_real(options_database, key, i)
      end do
    else
      if (options_has_key(options_database, key)) then
        if (present(from)) then
          if (from .gt. 1) return          
        end if
        if (present(to)) then
          if (to .lt. 1) return
        end if
        array_data(1)=options_get_real(options_database, key)
      end if
    end if
  end subroutine options_get_real_array  

  !> Retrieves a logical value from the database that matches the provided key
  !! @param options_database The options database
  !! @param key The key to search for and return the matching value
  !! @param index Optional array index to look up an array value
  !! @returns The matching logical value
  logical function options_get_logical(options_database, key, index)
    type(map_type), intent(inout) :: options_database
    character(len=*), intent(in) :: key
    integer, intent(in), optional :: index

    class(*), pointer :: raw_data
    character(len=len(key)+ARRAY_APPEND_SIZE) :: lookup_key

    if (present(index)) then
      lookup_key=get_options_array_key(key, index)
    else
      lookup_key=key
    end if
    call check_options_key_exists(options_database, lookup_key)
    raw_data=>c_get(options_database, lookup_key)
    if (associated(raw_data)) options_get_logical = conv_to_logical(raw_data, .false.)
  end function options_get_logical

  !> Retrieves an entire (or subset) logical array
  !! @param options_database The options database
  !! @param key The key to search for
  !! @param array_data The array data to write into
  !! @param from Optional starting index
  !! @param to Optional end index
  subroutine options_get_logical_array(options_database, key, array_data, from, to)
    type(map_type), intent(inout) :: options_database
    character(len=*), intent(in) :: key
    logical, dimension(:), intent(inout) :: array_data
    integer, intent(in), optional :: from, to

    integer :: num_elements, i, start, end

    num_elements=options_get_array_size(options_database, key)
    if (num_elements .gt. 0 .and. options_has_key(options_database, trim(key)//"a_size")) then
      if (present(from)) then
        start=from
      else
        start=1
      end if
      if (present(to)) then
        end=to
      else
        end=num_elements
      end if
      do i=start, num_elements
        array_data((i-start)+1)=options_get_logical(options_database, key, i)
      end do
    else
      if (options_has_key(options_database, key)) then
        if (present(from)) then
          if (from .gt. 1) return          
        end if
        if (present(to)) then
          if (to .lt. 1) return
        end if
        array_data(1)=options_get_logical(options_database, key)
      end if
    end if
  end subroutine options_get_logical_array  

  !> Retrieves an integer value from the database that matches the provided key
  !! @param options_database The options database
  !! @param key The key to search for and return the matching value
  !! @param index Optional array index to look up an array value
  !! @returns The matching integer value
  integer function options_get_integer(options_database, key, index)
    type(map_type), intent(inout) :: options_database
    character(len=*), intent(in) :: key
    integer, intent(in), optional :: index

    class(*), pointer :: raw_data
    character(len=len(key)+ARRAY_APPEND_SIZE) :: lookup_key

    if (present(index)) then
      lookup_key=get_options_array_key(key, index)
    else
      lookup_key=key
    end if
    call check_options_key_exists(options_database, lookup_key)
    raw_data=>c_get(options_database, lookup_key)
    if (associated(raw_data)) options_get_integer = conv_to_integer(raw_data, .false.)
  end function options_get_integer

  !> Retrieves an entire (or subset) integer array
  !! @param options_database The options database
  !! @param key The key to search for
  !! @param array_data The array data to write into
  !! @param from Optional starting index
  !! @param to Optional end index
  subroutine options_get_integer_array(options_database, key, array_data, from, to)
    type(map_type), intent(inout) :: options_database
    character(len=*), intent(in) :: key
    integer, dimension(:), intent(inout) :: array_data
    integer, intent(in), optional :: from, to

    integer :: num_elements, i, start, end

    num_elements=options_get_array_size(options_database, key)
    if (num_elements .gt. 0 .and. options_has_key(options_database, trim(key)//"a_size")) then
      if (present(from)) then
        start=from
      else
        start=1
      end if
      if (present(to)) then
        end=to
      else
        end=num_elements
      end if
      do i=start, num_elements
        array_data((i-start)+1)=options_get_integer(options_database, key, i)
      end do
    else
      if (options_has_key(options_database, key)) then
        if (present(from)) then
          if (from .gt. 1) return          
        end if
        if (present(to)) then
          if (to .lt. 1) return
        end if
        array_data(1)=options_get_integer(options_database, key)
      end if
    end if
  end subroutine options_get_integer_array  

  !> Retrieves a string value from the database that matches the provided key
  !! @param options_database The options database
  !! @param key The key to search for and return the matching value
  !! @param index Optional array index to look up an array value
  !! @returns The matching string value
  character(len=STRING_LENGTH) function options_get_string(options_database, key, index)
    type(map_type), intent(inout) :: options_database
    character(len=*), intent(in) :: key
    integer, intent(in), optional :: index

    class(*), pointer :: raw_data
    character(len=len(key)+ARRAY_APPEND_SIZE) :: lookup_key

    if (present(index)) then
      lookup_key=get_options_array_key(key, index)
    else
      lookup_key=key
    end if
    call check_options_key_exists(options_database, lookup_key)
    raw_data=>c_get(options_database, lookup_key)
    if (associated(raw_data)) options_get_string = conv_to_string(raw_data, .false., STRING_LENGTH)
  end function options_get_string

  !> Retrieves an entire (or subset) string array
  !! @param options_database The options database
  !! @param key The key to search for
  !! @param array_data The array data to write into
  !! @param from Optional starting index
  !! @param to Optional end index
  subroutine options_get_string_array(options_database, key, array_data, from, to)
    type(map_type), intent(inout) :: options_database
    character(len=*), intent(in) :: key
    character(len=STRING_LENGTH), dimension(:), intent(inout) :: array_data
    integer, intent(in), optional :: from, to

    integer :: num_elements, i, start, end

    num_elements=options_get_array_size(options_database, key)
    if (num_elements .gt. 0 .and. options_has_key(options_database, trim(key)//"a_size")) then
      if (present(from)) then
        start=from
      else
        start=1
      end if
      if (present(to)) then
        end=to
      else
        end=num_elements
      end if
      do i=start, num_elements
        array_data((i-start)+1)=options_get_string(options_database, key, i)
      end do
    else
      if (options_has_key(options_database, key)) then
        if (present(from)) then
          if (from .gt. 1) return          
        end if
        if (present(to)) then
          if (to .lt. 1) return
        end if
        array_data(1)=options_get_string(options_database, key)
      end if
    end if
  end subroutine options_get_string_array  

  !> Gets the size of the array held in the options database corresponding to a specific key
  !! @param options_database The options database
  !! @param key The key to look up
  !! @returns The number of elements in the array, 0 means zero elements (key is not in the database or corresponds to scalar.)
  integer function options_get_array_size(options_database, key)
    type(map_type), intent(inout) :: options_database
    character(len=*), intent(in) :: key

    class(*), pointer :: raw_data

    if (options_has_key(options_database, trim(key)//"a_size")) then
      options_get_array_size=options_get_integer(options_database, trim(key)//"a_size")
    else if (options_has_key(options_database, trim(key))) then
      options_get_array_size=1
    else
      options_get_array_size=0
    end if    
  end function options_get_array_size

  !> Loads in the command line arguments and stores them in the options database
  !! @returns map_type of option-value pairs
  subroutine load_command_line_into_options_database(options_database)
    type(map_type), intent(inout) :: options_database

    integer :: i, arguments, equals_posn, type_of_config
    character(len=LONG_STRING_LENGTH) :: specific_arg

    arguments = command_argument_count()

    do i=1,arguments
      call get_command_argument(i, value = specific_arg)
      specific_arg = trim(specific_arg)
      equals_posn=index(specific_arg, "=")  ! Get the position of the equals character
      if (index(specific_arg, "--") .eq. 1) then
        if (equals_posn .gt. 0) then
          type_of_config = get_argument_value_type(specific_arg(equals_posn+1 : len(specific_arg)))
        else
          ! If no equals then it is a logical switch (to true)
          type_of_config = LOGICAL_TYPE
        end if        
        call add_specific_option_key_value_pair(type_of_config, options_database, specific_arg)
      end if
    end do
  end subroutine load_command_line_into_options_database

  !> Removes a specific key from the options database, if it is an array then the entire array is removed
  !! @param options_database The options database
  !! @param key The key to locate and remove
  subroutine options_remove_key(options_database, key)
    type(map_type), intent(inout) :: options_database
    character(len=*), intent(in) :: key

    integer :: array_size, i
    
    if (options_has_key(options_database, key)) then          
      array_size=options_get_array_size(options_database, key)
      if (array_size .gt. 0) then
        do i=1,array_size
          if (options_has_key(options_database, get_options_array_key(key, i))) then
            call c_remove(options_database, get_options_array_key(key, i))
          end if
        end do
        call c_remove(options_database, trim(key)//"a_size")
      end if
    else 
      call c_remove(options_database, key)
    end if
  end subroutine options_remove_key  

  !--------------------------------------------------------------------------
  ! Private procedures acting as helpers to the functionality of the registry
  !--------------------------------------------------------------------------

  !> Determines whether a specific options key exists in the database or not, if it doesn't then this results in a log
  !! error being issued
  !! @param options_database The options database
  !! @param key The key to search for and return the matching value
  subroutine check_options_key_exists(options_database, key)
    type(map_type), intent(inout) :: options_database
    character(len=*), intent(in) :: key

    if (.not. options_has_key(options_database, key)) then
      call log_log(LOG_ERROR, "No configuration option with key "//trim(key)//" present")
    end if    
  end subroutine check_options_key_exists  

  !> Adds a real value to the options database with a specific key
  !! @param options_database The options database
  !! @param key The key to use
  !! @param value The real value to add to the database
  !! @param do_not_replace Optional flag whether to ignore existing values or replace them
  !! @param array_index Optional array index which specifies which location in the array to write to
  subroutine options_add_real(options_database, key, value, do_not_replace, array_index)
    type(map_type), intent(inout) :: options_database
    character(len=*), intent(in) :: key
    real, intent(in) :: value
    integer, intent(in), optional :: array_index
    logical, intent(in), optional :: do_not_replace

    class(*), pointer :: raw_data
    integer :: temp_size
    character(len=len(key)+ARRAY_APPEND_SIZE) :: lookup_key

    if (present(do_not_replace)) then
      if (do_not_replace) then
        if (present(array_index)) then
          if (options_has_key(options_database, trim(key)//"a_"//conv_to_string(array_index))) return
        else
          if (options_has_key(options_database, key)) return
        end if
      end if      
    end if

    raw_data => conv_to_generic(value, .true.)

    if (present(array_index)) then
      if (options_has_key(options_database, trim(key)//"a_size")) then
        temp_size=options_get_integer(options_database, trim(key)//"a_size")
        if (temp_size .lt. array_index) temp_size=temp_size+1
      else
        temp_size=1
      end if
      call options_add(options_database, trim(key)//"a_size", temp_size)      
      lookup_key=get_options_array_key(key, array_index)
    else
      lookup_key=key
    end if
    call c_put(options_database, lookup_key, raw_data)
  end subroutine options_add_real

  !> Adds a logical value to the options database with a specific key
  !! @param options_database The options database
  !! @param key The key to use
  !! @param value The logical value to add to the database
  !! @param do_not_replace Optional flag whether to ignore existing values or replace them
  !! @param array_index Optional array index which specifies which location in the array to write to
  subroutine options_add_logical(options_database, key, value, do_not_replace, array_index)
    type(map_type), intent(inout) :: options_database
    character(len=*), intent(in) :: key
    logical, intent(in) :: value
    integer, intent(in), optional :: array_index
    logical, intent(in), optional :: do_not_replace

    class(*), pointer :: raw_data
    integer :: temp_size
    character(len=len(key)+ARRAY_APPEND_SIZE) :: lookup_key

    if (present(do_not_replace)) then
      if (do_not_replace) then
        if (present(array_index)) then
          if (options_has_key(options_database, trim(key)//"a_"//conv_to_string(array_index))) return
        else
          if (options_has_key(options_database, key)) return
        end if
      end if      
    end if

    raw_data => conv_to_generic(value, .true.)

    if (present(array_index)) then
      if (options_has_key(options_database, trim(key)//"a_size")) then
        temp_size=options_get_integer(options_database, trim(key)//"a_size")
        if (temp_size .lt. array_index) temp_size=temp_size+1
      else
        temp_size=1
      end if
      call options_add(options_database, trim(key)//"a_size", temp_size)
      lookup_key=get_options_array_key(key, array_index)
    else
      lookup_key=key
    end if
    call c_put(options_database, lookup_key, raw_data)    
  end subroutine options_add_logical

  !> Adds a string value to the options database with a specific key
  !! @param options_database The options database
  !! @param key The key to use
  !! @param value The string value to add to the database
  !! @param do_not_replace Optional flag whether to ignore existing values or replace them
  !! @param array_index Optional array index which specifies which location in the array to write to
  subroutine options_add_string(options_database, key, value, do_not_replace, array_index)
    type(map_type), intent(inout) :: options_database
    character(len=*), intent(in) :: key, value
    integer, intent(in), optional :: array_index
    logical, intent(in), optional :: do_not_replace

    character(len=STRING_LENGTH) :: value_to_store
    class(*), pointer :: raw_data
    integer :: temp_size
    character(len=len(key)+ARRAY_APPEND_SIZE) :: lookup_key

    if (present(do_not_replace)) then
      if (do_not_replace) then
        if (present(array_index)) then
          if (options_has_key(options_database, trim(key)//"a_"//conv_to_string(array_index))) return
        else
          if (options_has_key(options_database, key)) return
        end if
      end if      
    end if

    value_to_store=value
    raw_data => conv_to_generic(value_to_store, .true.)

    if (present(array_index)) then
      if (options_has_key(options_database, trim(key)//"a_size")) then
        temp_size=options_get_integer(options_database, trim(key)//"a_size")
        if (temp_size .lt. array_index) temp_size=temp_size+1
      else
        temp_size=1
      end if
      call options_add(options_database, trim(key)//"a_size", temp_size)
      lookup_key=get_options_array_key(key, array_index)
    else
      lookup_key=key
    end if
    call c_put(options_database, lookup_key, raw_data)
  end subroutine options_add_string

  !> Adds an integer value to the options database with a specific key
  !! @param options_database The options database
  !! @param key The key to use
  !! @param value The integer value to add to the database
  !! @param do_not_replace Optional flag whether to ignore existing values or replace them
  !! @param array_index Optional array index which specifies which location in the array to write to
  recursive subroutine options_add_integer(options_database, key, value, do_not_replace, array_index)
    type(map_type), intent(inout) :: options_database
    character(len=*), intent(in) :: key
    integer, intent(in) :: value
    integer, intent(in), optional :: array_index
    logical, intent(in), optional :: do_not_replace

    class(*), pointer :: raw_data
    integer :: temp_size
    character(len=len(key)+ARRAY_APPEND_SIZE) :: lookup_key

    if (present(do_not_replace)) then
      if (do_not_replace) then
        if (present(array_index)) then
          if (options_has_key(options_database, trim(key)//"a_"//conv_to_string(array_index))) return
        else
          if (options_has_key(options_database, key)) return
        end if
      end if      
    end if

    raw_data => conv_to_generic(value, .true.)

    if (present(array_index)) then
      if (options_has_key(options_database, trim(key)//"a_size")) then
        temp_size=options_get_integer(options_database, trim(key)//"a_size")
        if (temp_size .lt. array_index) temp_size=temp_size+1
      else
        temp_size=1
      end if
      call options_add(options_database, trim(key)//"a_size", temp_size)
      lookup_key=get_options_array_key(key, array_index)
    else
      lookup_key=key
    end if
    call c_put(options_database, lookup_key, raw_data)    
  end subroutine options_add_integer

  !> Gets a key corresponding to the correct options key and index combination
  !! @param key The key
  !! @param index The index
  !! @returns The key to store the indexed item in
  character(len=STRING_LENGTH) function get_options_array_key(key, index)
    character(len=*), intent(in) :: key
    integer, intent(in) :: index

    get_options_array_key=trim(key)//"a_"//trim(conv_to_string(index))
  end function get_options_array_key  

  !> Given a specific value this will determine the type of data
  !!
  !! It basically wraps around the conversion test utilities, but will select a real if
  !! the string contains a dot and an integer if not. If no conversions_mod can be achieved then
  !! the value is assumed to be of type string
  !! @param specific_value String value to deduce the type of
  !! @returns Integer representing the type as per the module parameter type integers
  integer function get_argument_value_type(specific_value)
    character(len=*), intent(in) :: specific_value
    integer :: dot_posn, exponent_posn

    if (conv_is_logical(specific_value)) then
      get_argument_value_type = LOGICAL_TYPE
      return
      else
    end if

    dot_posn = index(specific_value,".")
    exponent_posn = index(specific_value,"e")
    if (dot_posn .eq. 0 .and. exponent_posn .eq. 0) then
      ! No dot, therefore an integer or string
      if (conv_is_integer(specific_value)) then
        get_argument_value_type = INTEGER_TYPE
        return
      end if
    else
      ! A dot is present, therefore a real or string
      if (conv_is_real(specific_value)) then
        get_argument_value_type = REAL_TYPE
        return
      end if
    end if
    get_argument_value_type = STRING_TYPE ! Default string type
  end function get_argument_value_type

  !> This will add a specific option key value pair to the options map_type.
  !!
  !! As per the map_type semantics, if the user supplies an option then this will override
  !! any existing default value.
  !! @param type_of_config The type of value to add to the configuration map_type
  !! @param parse_options The key-value pairs of configuration options
  !! @param specific_arg The argument string (key=value or just key)
  subroutine add_specific_option_key_value_pair(type_of_config, parse_options, specific_arg)
    integer, intent(in) :: type_of_config
    character(len=*), intent(in) :: specific_arg
    type(map_type), intent(inout) :: parse_options

    integer :: equals_posn
    equals_posn=index(specific_arg, "=")  ! Get the position of the equals character

    if (type_of_config == LOGICAL_TYPE) then
      if (equals_posn .gt. 0) then
        call set_options_logical_value(parse_options, specific_arg(3:equals_posn-1), &
          conv_to_logical(specific_arg(equals_posn+1:len(specific_arg))))
      else
        ! Handle the case where there is just a key (no value)
        call set_options_logical_value(parse_options, specific_arg(3:len(specific_arg)), .true.)
      end if
    else if (type_of_config == INTEGER_TYPE) then
      call set_options_integer_value(parse_options, specific_arg(3:equals_posn-1), &
        conv_to_integer(specific_arg(equals_posn+1:len(specific_arg))))
    else if (type_of_config == REAL_TYPE) then
      call set_options_real_value(parse_options, specific_arg(3:equals_posn-1), &
        conv_to_real(specific_arg(equals_posn+1:len(specific_arg))))
    else if (type_of_config == STRING_TYPE) then
      call set_options_string_value(parse_options, specific_arg(3:equals_posn-1), specific_arg(equals_posn+1:len(specific_arg)))
    end if
  end subroutine add_specific_option_key_value_pair

  !> A helper procedure to set a specific logical value
  subroutine set_options_logical_value(optionmap_type, key, value)
    type(map_type), intent(inout) :: optionmap_type
    character(len=*), intent(in) :: key
    logical, intent(in) :: value

    class(*), pointer :: generic_data

    generic_data=>conv_to_generic(value, .true.)
    call c_put(optionmap_type, key, generic_data)
  end subroutine set_options_logical_value

  !> A helper procedure to set a specific real value
  subroutine set_options_real_value(optionmap_type, key, value)
    type(map_type), intent(inout) :: optionmap_type
    character(len=*), intent(in) :: key
    real, intent(in) :: value

    class(*), pointer :: generic_data

    generic_data=>conv_to_generic(value, .true.)
    call c_put(optionmap_type, key, generic_data)
  end subroutine set_options_real_value

  !> A helper procedure to set a specific integer value
  subroutine set_options_integer_value(optionmap_type, key, value)
    type(map_type), intent(inout) :: optionmap_type
    character(len=*), intent(in) :: key
    integer, intent(in) :: value

    class(*), pointer :: generic_data

    generic_data=>conv_to_generic(value, .true.)
    call c_put(optionmap_type, key, generic_data)
  end subroutine set_options_integer_value

  !> A helper procedure to set a specific string value
  subroutine set_options_string_value(optionmap_type, key, value)
    type(map_type), intent(inout) :: optionmap_type
    character(len=*), intent(in) :: key
    character(len=*), intent(in) :: value

    character(len=STRING_LENGTH) :: write_value
    class(*), pointer :: generic_data

    ! We do an assignment here to force the value to be 44 characters, otherwise might underflow
    write_value = value
    generic_data=>conv_to_generic(write_value, .true.)
    call c_put(optionmap_type, key, generic_data)
  end subroutine set_options_string_value
end module optionsdatabase_mod
