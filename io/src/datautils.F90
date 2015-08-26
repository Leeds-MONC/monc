!> Contains functionality for managing and extracting data from the raw data dumps that the IO
!! server receives from a MONC process during model dumping
module data_utils_mod
  use datadefn_mod, only : DOUBLE_PRECISION, SINGLE_PRECISION, STRING_LENGTH
  use conversions_mod, only : conv_to_integer, conv_to_string, conv_to_generic, conv_is_integer
  use collections_mod, only : map_type, c_get, c_put, c_contains
  use io_server_client_mod, only : INTEGER_DATA_TYPE, BOOLEAN_DATA_TYPE, STRING_DATA_TYPE, FLOAT_DATA_TYPE, &
       DOUBLE_DATA_TYPE
  use configuration_parser_mod, only : io_configuration_type
  use logging_mod, only : LOG_ERROR, log_log
  implicit none

#ifndef TEST_MODE
  private
#endif

  public is_field_present, get_map, get_map_from_monc, get_scalar_integer, get_scalar_integer_from_monc, get_scalar_real, &
       get_scalar_real_from_monc, get_array_double, get_array_double_from_monc, get_3darray_double, &
       get_3darray_double_from_monc, get_4darray_double, get_4darray_double_from_monc, get_2darray_double, &
       get_2darray_double_from_monc, get_array_integer, get_array_integer_from_monc, &
       get_string, get_string_from_monc, get_field_size, get_action_attribute_string, get_action_attribute_logical, &
       get_action_attribute_integer
contains

  !> Retrieves the name of a field from the attributes specified in the configuration
  !! @param action_attributes Action attributes from the IO server configuration
  !! @param Name of the field to perform the collective operation over
  character(len=STRING_LENGTH) function get_action_attribute_string(action_attributes, field_name)
    type(map_type), intent(inout) :: action_attributes
    character(len=*), intent(in) :: field_name

    class(*), pointer :: generic

    if (.not. c_contains(action_attributes, field_name)) call log_log(LOG_ERROR, &
         "You must provide the field name in the collective operation configuration")

    generic=>c_get(action_attributes, field_name)
    get_action_attribute_string=conv_to_string(generic, .false., STRING_LENGTH)
  end function get_action_attribute_string

  !> Retrieves the name of a field from the attributes specified in the configuration
  !! @param action_attributes Action attributes from the IO server configuration
  !! @param Name of the field to perform the collective operation over
  integer function get_action_attribute_integer(action_attributes, field_name)
    type(map_type), intent(inout) :: action_attributes
    character(len=*), intent(in) :: field_name

    character(len=STRING_LENGTH) :: str_val

    str_val=get_action_attribute_string(action_attributes, field_name)
    if (.not. conv_is_integer(str_val)) call log_log(LOG_ERROR, "Can not convert string '"//trim(str_val)//"' to an integer")
    get_action_attribute_integer=conv_to_integer(str_val)
  end function get_action_attribute_integer

  !> Retrieves a logical value from the attribute which corresponds to a specific key
  !! @param action_attributes The attributes provided to the action
  !! @param name The name of the attribute to look up
  logical function get_action_attribute_logical(action_attributes, field_name)
    type(map_type), intent(inout) :: action_attributes
    character(len=*), intent(in) :: field_name

    class(*), pointer :: generic

    if (c_contains(action_attributes, field_name)) then
      generic=>c_get(action_attributes, field_name)
      get_action_attribute_logical=conv_to_string(generic, .false., STRING_LENGTH) .eq. "true"
    else
      get_action_attribute_logical=.false.
    end if
  end function get_action_attribute_logical

  !! Allows one to check if an optional field is present in the data being provided by a MONC
  !! process or not
  !! @param io_configuration Configuration of the IO server
  !! @param source PID of the MONC process
  !! @param data_id The ID of the data definition that is represented by the dump
  !! @param key Key of the field to retrieve
  !! @returns Whether the field is present or not
  logical function is_field_present(io_configuration, source, data_id, key)
    type(io_configuration_type), intent(inout) :: io_configuration
    integer, intent(in) :: source, data_id
    character(len=*), intent(in) :: key

    integer :: start_index, end_index, monc_location
    class(*), pointer :: generic

    generic=>c_get(io_configuration%monc_to_index, conv_to_string(source))
    monc_location=conv_to_integer(generic, .false.)

    generic=>c_get(io_configuration%registered_moncs(monc_location)%field_start_locations(data_id), key)
    if (.not. associated(generic)) then
      is_field_present=.false.
      return
    end if
    start_index=conv_to_integer(generic, .false.)
    generic=>c_get(io_configuration%registered_moncs(monc_location)%field_end_locations(data_id), key)
    if (.not. associated(generic)) then
      is_field_present=.false.
      return
    end if
    end_index=conv_to_integer(generic, .false.)

    is_field_present = end_index .gt. start_index
  end function is_field_present

  !> Retrieves the size of a field from the data definition
  !! @param field_starts Field starting locations
  !! @param field_ends Field ending locations
  !! @param key The key (name) of the field to check for
  !! @param data_type The data type of the field
  !! @returns The number of elements that this field contains (of type data_type)
  integer function get_field_size(field_starts, field_ends, key, data_type)
    type(map_type), intent(inout) :: field_starts, field_ends
    character(len=*), intent(in) :: key
    integer, intent(in) :: data_type

    integer :: start_index, end_index, element_size
    real(kind=DOUBLE_PRECISION) :: dreal
    real(kind=SINGLE_PRECISION) :: sreal
    class(*), pointer :: generic

    generic=>c_get(field_starts, key)
    start_index=conv_to_integer(generic, .false.)
    generic=>c_get(field_ends, key)
    end_index=conv_to_integer(generic, .false.)

    if (data_type == INTEGER_DATA_TYPE) then
      element_size=kind(start_index)
    else if (data_type == DOUBLE_DATA_TYPE) then
      element_size=kind(dreal)
    else if (data_type == FLOAT_DATA_TYPE) then
      element_size=kind(sreal)
    else if (data_type == STRING_DATA_TYPE) then
      element_size=STRING_LENGTH
    end if
    get_field_size=((end_index-start_index)+1)/element_size
  end function get_field_size

  !> Retrieves a map data structure with key->value pairs, each of which are strings
  !! @param field_starts Field starting locations
  !! @param field_ends Field ending locations
  !! @param data_dump Raw data that MONC has sent to us
  !! @param key Key of the field to retrieve
  !! @returns Corresponding map and contents
  type(map_type) function get_map(field_starts, field_ends, data_dump, key)
    type(map_type), intent(inout) :: field_starts, field_ends
    character, dimension(:) :: data_dump
    character(len=*), intent(in) :: key
    
    integer :: start_index, end_index, elements, i
    class(*), pointer :: generic
    character(len=STRING_LENGTH) :: retrieved1, retrieved2

    if (.not. c_contains(field_starts, key) .or. .not. c_contains(field_ends, key)) &
         call log_log(LOG_ERROR, "Field name `"//key//"` not found in the data definition")

    generic=>c_get(field_starts, key)
    start_index=conv_to_integer(generic, .false.)
    generic=>c_get(field_ends, key)
    end_index=conv_to_integer(generic, .false.)

    elements = (end_index+1 - start_index) / (STRING_LENGTH*2)

    do i=1, elements
      retrieved1=transfer(data_dump(start_index:start_index+STRING_LENGTH-1), retrieved1)
      start_index=start_index+STRING_LENGTH
      retrieved2=transfer(data_dump(start_index:start_index+STRING_LENGTH-1), retrieved2)
      start_index=start_index+STRING_LENGTH
      call c_put(get_map, retrieved1, conv_to_generic(retrieved2, .true.))
    end do    
  end function get_map  

  !> Retrieves a map data structure with key->value pairs, each of which are strings
  !! @param io_configuration Configuration of the IO server
  !! @param source PID of the MONC process
  !! @param data_dump Raw data that MONC has sent to us
  !! @param data_id The ID of the data definition that is represented by the dump
  !! @param key Key of the field to retrieve
  !! @returns Corresponding map and contents
  type(map_type) function get_map_from_monc(io_configuration, source, data_id, data_dump, key)
    type(io_configuration_type), intent(inout) :: io_configuration
    integer, intent(in) :: source, data_id
    character, dimension(:) :: data_dump
    character(len=*), intent(in) :: key

    integer :: monc_location
    class(*), pointer :: generic

    generic=>c_get(io_configuration%monc_to_index, conv_to_string(source))
    monc_location=conv_to_integer(generic, .false.)

    get_map_from_monc=get_map(io_configuration%registered_moncs(monc_location)%field_start_locations(data_id), &
         io_configuration%registered_moncs(monc_location)%field_end_locations(data_id), data_dump, key)
  end function get_map_from_monc

  !> Retrieves a string from the data dump
  !! @param field_starts Field starting locations
  !! @param field_ends Field ending locations
  !! @param data_dump Raw data that MONC has sent to us
  !! @param key Key of the field to retrieve
  !! @returns Corresponding string
  function get_string(field_starts, field_ends, data_dump, key)
    type(map_type), intent(inout) :: field_starts, field_ends
    character, dimension(:) :: data_dump
    character(len=*), intent(in) :: key
    character(len=STRING_LENGTH) :: get_string

    integer :: start_index, end_index
    class(*), pointer :: generic

    if (.not. c_contains(field_starts, key) .or. .not. c_contains(field_ends, key)) &
         call log_log(LOG_ERROR, "Field name `"//key//"` not found in the data definition")

    generic=>c_get(field_starts, key)
    start_index=conv_to_integer(generic, .false.)
    generic=>c_get(field_ends, key)
    end_index=conv_to_integer(generic, .false.)

    get_string=transfer(data_dump(start_index:end_index), get_string)
  end function get_string

  !> Retrieves a string from the data dump
  !! @param io_configuration Configuration of the IO server
  !! @param source PID of the MONC process
  !! @param data_dump Raw data that MONC has sent to us
  !! @param data_id The ID of the data definition that is represented by the dump
  !! @param key Key of the field to retrieve
  !! @returns Corresponding string
  function get_string_from_monc(io_configuration, source, data_id, data_dump, key)
    type(io_configuration_type), intent(inout) :: io_configuration
    integer, intent(in) :: source, data_id
    character, dimension(:) :: data_dump
    character(len=*), intent(in) :: key
    character(len=STRING_LENGTH) :: get_string_from_monc

    integer :: monc_location
    class(*), pointer :: generic

    generic=>c_get(io_configuration%monc_to_index, conv_to_string(source))
    monc_location=conv_to_integer(generic, .false.)

    get_string_from_monc=get_string(io_configuration%registered_moncs(monc_location)%field_start_locations(data_id), &
         io_configuration%registered_moncs(monc_location)%field_end_locations(data_id), data_dump, key)
  end function get_string_from_monc  

  !> Retrieves a single integer element (scalar) from the data dump
  !! @param field_starts Field starting locations
  !! @param field_ends Field ending locations
  !! @param data_dump Raw data that MONC has sent to us
  !! @param key Key of the field to retrieve
  !! @returns Corresponding integer
  integer function get_scalar_integer(field_starts, field_ends, data_dump, key)
    type(map_type), intent(inout) :: field_starts, field_ends
    character, dimension(:) :: data_dump
    character(len=*), intent(in) :: key

    integer :: start_index, end_index
    class(*), pointer :: generic

    if (.not. c_contains(field_starts, key) .or. .not. c_contains(field_ends, key)) &
         call log_log(LOG_ERROR, "Field name `"//key//"` not found in the data definition")

    generic=>c_get(field_starts, key)
    start_index=conv_to_integer(generic, .false.)
    generic=>c_get(field_ends, key)
    end_index=conv_to_integer(generic, .false.)

    get_scalar_integer=transfer(data_dump(start_index:end_index), get_scalar_integer)
  end function get_scalar_integer

  !> Retrieves a single integer element (scalar) from the data dump
  !! @param io_configuration Configuration of the IO server
  !! @param source PID of the MONC process
  !! @param data_dump Raw data that MONC has sent to us
  !! @param data_id The ID of the data definition that is represented by the dump
  !! @param key Key of the field to retrieve
  !! @returns Corresponding integer
  integer function get_scalar_integer_from_monc(io_configuration, source, data_id, data_dump, key)
    type(io_configuration_type), intent(inout) :: io_configuration
    integer, intent(in) :: source, data_id
    character, dimension(:) :: data_dump
    character(len=*), intent(in) :: key

    integer :: monc_location
    class(*), pointer :: generic

    generic=>c_get(io_configuration%monc_to_index, conv_to_string(source))
    monc_location=conv_to_integer(generic, .false.)

    get_scalar_integer_from_monc=get_scalar_integer(&
         io_configuration%registered_moncs(monc_location)%field_start_locations(data_id), &
         io_configuration%registered_moncs(monc_location)%field_end_locations(data_id), data_dump, key)
  end function get_scalar_integer_from_monc

  !> Retreives a scalar real with a corresponding key from the raw data dump
  !! @param field_starts Field starting locations
  !! @param field_ends Field ending locations
  !! @param data_dump Raw data that MONC has sent to us
  !! @param key Key of the field to retrieve
  !! @returns Real value of the corresponding field
  real(kind=DOUBLE_PRECISION) function get_scalar_real(field_starts, field_ends, data_dump, key)
    type(map_type), intent(inout) :: field_starts, field_ends
    character, dimension(:) :: data_dump
    character(len=*), intent(in) :: key

    integer :: start_index, end_index
    class(*), pointer :: generic

    if (.not. c_contains(field_starts, key) .or. .not. c_contains(field_ends, key)) &
         call log_log(LOG_ERROR, "Field name `"//key//"` not found in the data definition")

    generic=>c_get(field_starts, key)
    start_index=conv_to_integer(generic, .false.)
    generic=>c_get(field_ends, key)
    end_index=conv_to_integer(generic, .false.)

    get_scalar_real=transfer(data_dump(start_index:end_index), get_scalar_real)
  end function get_scalar_real

  !> Retreives a scalar real with a corresponding key from the raw data dump
  !! @param io_configuration Configuration of the IO server
  !! @param source PID of the MONC process
  !! @param data_id The ID of the data definition that is represented by the dump
  !! @param data_dump Raw data that MONC has sent to us
  !! @param key Key of the field to retrieve
  !! @returns Real value of the corresponding field
  real(kind=DOUBLE_PRECISION) function get_scalar_real_from_monc(io_configuration, source, data_id, data_dump, key)
    type(io_configuration_type), intent(inout) :: io_configuration
    integer, intent(in) :: source, data_id
    character, dimension(:) :: data_dump
    character(len=*), intent(in) :: key

    integer :: monc_location
    class(*), pointer :: generic

    generic=>c_get(io_configuration%monc_to_index, conv_to_string(source))
    monc_location=conv_to_integer(generic, .false.)

    get_scalar_real_from_monc=get_scalar_real(&
         io_configuration%registered_moncs(monc_location)%field_start_locations(data_id), &
         io_configuration%registered_moncs(monc_location)%field_end_locations(data_id), data_dump, key)
  end function get_scalar_real_from_monc

  !> Retreives an array of doubles with a corresponding key from the raw data dump. The size depends on the configuration
  !! that the MONC process supplied to the IO server when it registered which informs of the size of a field
  !! @param field_starts Field starting locations
  !! @param field_ends Field ending locations
  !! @param data_dump Raw data that MONC has sent to us
  !! @param key Key of the field to retrieve
  !! @returns Double array values of the corresponding field
  function get_array_double(field_starts, field_ends, data_dump, key)
    type(map_type), intent(inout) :: field_starts, field_ends
    character, dimension(:) :: data_dump
    character(len=*), intent(in) :: key
    real(kind=DOUBLE_PRECISION), dimension(:), allocatable :: get_array_double

    integer :: start_index, end_index, elements
    class(*), pointer :: generic    

    if (.not. c_contains(field_starts, key) .or. .not. c_contains(field_ends, key)) &
         call log_log(LOG_ERROR, "Field name `"//key//"` not found in the data definition")

    generic=>c_get(field_starts, key)
    start_index=conv_to_integer(generic, .false.)
    generic=>c_get(field_ends, key)
    end_index=conv_to_integer(generic, .false.)

    elements = (end_index - start_index) / kind(get_array_double)

    allocate(get_array_double(elements))
    get_array_double=transfer(data_dump(start_index:end_index), get_array_double)
  end function get_array_double

  !> Retreives an array of doubles with a corresponding key from the raw data dump. The size depends on the configuration
  !! that the MONC process supplied to the IO server when it registered which informs of the size of a field
  !! @param io_configuration Configuration of the IO server
  !! @param source PID of the MONC process
  !! @param data_id The ID of the data definition that is represented by the dump
  !! @param data_dump Raw data that MONC has sent to us
  !! @param key Key of the field to retrieve
  !! @returns Double array values of the corresponding field
  function get_array_double_from_monc(io_configuration, source, data_id, data_dump, key)
    type(io_configuration_type), intent(inout) :: io_configuration
    integer, intent(in) :: source, data_id
    character, dimension(:) :: data_dump
    character(len=*), intent(in) :: key
    real(kind=DOUBLE_PRECISION), dimension(:), allocatable :: get_array_double_from_monc

    integer :: monc_location
    class(*), pointer :: generic

    generic=>c_get(io_configuration%monc_to_index, conv_to_string(source))
    monc_location=conv_to_integer(generic, .false.)

    get_array_double_from_monc=get_array_double(&
         io_configuration%registered_moncs(monc_location)%field_start_locations(data_id), &
         io_configuration%registered_moncs(monc_location)%field_end_locations(data_id), data_dump, key)
  end function get_array_double_from_monc

  !> Retreives an array of integers with a corresponding key from the raw data dump. The size depends on the configuration
  !! that the MONC process supplied to the IO server when it registered which informs of the size of a field
  !! @param field_starts Field starting locations
  !! @param field_ends Field ending locations
  !! @param data_dump Raw data that MONC has sent to us
  !! @param key Key of the field to retrieve
  !! @returns Integer array values of the corresponding field
  function get_array_integer(field_starts, field_ends, data_dump, key)
    type(map_type), intent(inout) :: field_starts, field_ends
    character, dimension(:) :: data_dump
    character(len=*), intent(in) :: key
    integer, dimension(:), allocatable :: get_array_integer

    integer :: start_index, end_index, elements
    class(*), pointer :: generic

    if (.not. c_contains(field_starts, key) .or. .not. c_contains(field_ends, key)) &
         call log_log(LOG_ERROR, "Field name `"//key//"` not found in the data definition")

    generic=>c_get(field_starts, key)
    start_index=conv_to_integer(generic, .false.)
    generic=>c_get(field_ends, key)
    end_index=conv_to_integer(generic, .false.)

    elements = (end_index - start_index) / kind(get_array_integer)

    allocate(get_array_integer(elements))
    get_array_integer=transfer(data_dump(start_index:end_index), get_array_integer)
  end function get_array_integer  

  !> Retreives an array of integers with a corresponding key from the raw data dump. The size depends on the configuration
  !! that the MONC process supplied to the IO server when it registered which informs of the size of a field
  !! @param io_configuration Configuration of the IO server
  !! @param source PID of the MONC process
  !! @param data_id The ID of the data definition that is represented by the dump
  !! @param data_dump Raw data that MONC has sent to us
  !! @param key Key of the field to retrieve
  !! @returns Integer array values of the corresponding field
  function get_array_integer_from_monc(io_configuration, source, data_id, data_dump, key)
    type(io_configuration_type), intent(inout) :: io_configuration
    integer, intent(in) :: source, data_id
    character, dimension(:) :: data_dump
    character(len=*), intent(in) :: key
    integer, dimension(:), allocatable :: get_array_integer_from_monc

    integer ::  monc_location
    class(*), pointer :: generic

    generic=>c_get(io_configuration%monc_to_index, conv_to_string(source))
    monc_location=conv_to_integer(generic, .false.)

    get_array_integer_from_monc=get_array_integer(&
         io_configuration%registered_moncs(monc_location)%field_start_locations(data_id), &
         io_configuration%registered_moncs(monc_location)%field_end_locations(data_id), data_dump, key)   
  end function get_array_integer_from_monc

  !> Retreives a 2D array of doubles with a corresponding key from the raw data dump. The size depends on the configuration
  !! that the MONC process supplied to the IO server when it registered which informs of the size of a field
  !!
  !! Note that this replies in the F2008 `contiguous` keyword, but it is supported by the target compilers so we use it
  !! @param io_configuration Configuration of the IO server
  !! @param source PID of the MONC process
  !! @param data_id The ID of the data definition that is represented by the dump
  !! @param data_dump Raw data that MONC has sent to us
  !! @param key Key of the field to retrieve
  !! @param target_data The data that we are writing into
  !! @param size1 Size in first dimension
  !! @param size2 Size in second dimension
  subroutine get_2darray_double_from_monc(io_configuration, source, data_id, data_dump, key, target_data, size1, size2)
    type(io_configuration_type), intent(inout) :: io_configuration
    integer, intent(in) :: source, data_id, size1, size2
    character, dimension(:) :: data_dump
    character(len=*), intent(in) :: key
    real(kind=DOUBLE_PRECISION), dimension(:,:), pointer, contiguous, intent(inout) :: target_data

    integer :: monc_location
    class(*), pointer :: generic

    generic=>c_get(io_configuration%monc_to_index, conv_to_string(source))
    monc_location=conv_to_integer(generic, .false.)

    call get_2darray_double(io_configuration%registered_moncs(monc_location)%field_start_locations(data_id), &
         io_configuration%registered_moncs(monc_location)%field_end_locations(data_id), data_dump, key, target_data, &
         size1, size2)    
  end subroutine get_2darray_double_from_monc

  !> Retreives a 2D array of doubles with a corresponding key from the raw data dump. The size depends on the configuration
  !! that the MONC process supplied to the IO server when it registered which informs of the size of a field
  !!
  !! Note that this replies in the F2008 `contiguous` keyword, but it is supported by the target compilers so we use it
  !! @param field_starts Field starting locations
  !! @param field_ends Field ending locations
  !! @param data_dump Raw data that MONC has sent to us
  !! @param key Key of the field to retrieve
  !! @param target_data The data that we are writing into
  !! @param size1 Size in first dimension
  !! @param size2 Size in second dimension
  subroutine get_2darray_double(field_starts, field_ends, data_dump, key, target_data, size1, size2)
    type(map_type), intent(inout) :: field_starts, field_ends
    integer, intent(in) ::  size1, size2
    character, dimension(:) :: data_dump
    character(len=*), intent(in) :: key
    real(kind=DOUBLE_PRECISION), dimension(:,:), pointer, contiguous, intent(inout) :: target_data

    integer :: start_index, end_index, element_size
    real(kind=DOUBLE_PRECISION), dimension(:), pointer :: temp_data
    class(*), pointer :: generic

    ! Pointer bounds remapping as transfer needs 1D array but for performance don't want to allocate another array and copy using reshape
    temp_data(1:size1*size2)=>target_data

    if (.not. c_contains(field_starts, key) .or. .not. c_contains(field_ends, key)) &
         call log_log(LOG_ERROR, "Field name `"//key//"` not found in the data definition")

    generic=>c_get(field_starts, key)
    start_index=conv_to_integer(generic, .false.)
    generic=>c_get(field_ends, key)
    end_index=conv_to_integer(generic, .false.)

    element_size=(end_index-start_index) / kind(target_data)

    temp_data=transfer(data_dump(start_index:end_index), temp_data)
  end subroutine get_2darray_double

  !> Retreives a 3D array of doubles with a corresponding key from the raw data dump. The size depends on the configuration
  !! that the MONC process supplied to the IO server when it registered which informs of the size of a field
  !!
  !! Note that this replies in the F2008 `contiguous` keyword, but it is supported by the target compilers so we use it
  !! @param field_starts Field starting locations
  !! @param field_ends Field ending locations
  !! @param data_dump Raw data that MONC has sent to us
  !! @param key Key of the field to retrieve
  !! @param target_data The data that we are writing into
  !! @param size1 Size in first dimension
  !! @param size2 Size in second dimension
  !! @param size3 Size in third dimension
  subroutine get_3darray_double(field_starts, field_ends, data_dump, key, target_data, size1, size2, size3)
    type(map_type), intent(inout) :: field_starts, field_ends
    integer, intent(in) ::  size1, size2, size3
    character, dimension(:) :: data_dump
    character(len=*), intent(in) :: key
    real(kind=DOUBLE_PRECISION), dimension(:,:,:), pointer, contiguous, intent(inout) :: target_data
    
    integer :: start_index, end_index, element_size
    real(kind=DOUBLE_PRECISION), dimension(:), pointer :: temp_data
    class(*), pointer :: generic

    ! Pointer bounds remapping as transfer needs 1D array but for performance don't want to allocate another array and copy using reshape
    temp_data(1:size1*size2*size3)=>target_data

    if (.not. c_contains(field_starts, key) .or. .not. c_contains(field_ends, key)) &
         call log_log(LOG_ERROR, "Field name `"//key//"` not found in the data definition")

    generic=>c_get(field_starts, key)
    start_index=conv_to_integer(generic, .false.)
    generic=>c_get(field_ends, key)
    end_index=conv_to_integer(generic, .false.)

    element_size=(end_index-start_index) / kind(target_data)
    
    temp_data=transfer(data_dump(start_index:end_index), temp_data)
  end subroutine get_3darray_double  

  !> Retreives a 3D array of doubles with a corresponding key from the raw data dump. The size depends on the configuration
  !! that the MONC process supplied to the IO server when it registered which informs of the size of a field
  !!
  !! Note that this replies in the F2008 `contiguous` keyword, but it is supported by the target compilers so we use it
  !! @param io_configuration Configuration of the IO server
  !! @param source PID of the MONC process
  !! @param data_id The ID of the data definition that is represented by the dump
  !! @param data_dump Raw data that MONC has sent to us
  !! @param key Key of the field to retrieve
  !! @param target_data The data that we are writing into
  !! @param size1 Size in first dimension
  !! @param size2 Size in second dimension
  !! @param size3 Size in third dimension
  subroutine get_3darray_double_from_monc(io_configuration, source, data_id, data_dump, key, target_data, size1, size2, size3)
    type(io_configuration_type), intent(inout) :: io_configuration
    integer, intent(in) :: source, data_id, size1, size2, size3
    character, dimension(:) :: data_dump
    character(len=*), intent(in) :: key
    real(kind=DOUBLE_PRECISION), dimension(:,:,:), pointer, contiguous, intent(inout) :: target_data

    integer :: monc_location
    class(*), pointer :: generic

    generic=>c_get(io_configuration%monc_to_index, conv_to_string(source))
    monc_location=conv_to_integer(generic, .false.)

    call get_3darray_double(io_configuration%registered_moncs(monc_location)%field_start_locations(data_id), &
         io_configuration%registered_moncs(monc_location)%field_end_locations(data_id), data_dump, key, target_data, &
         size1, size2, size3)    
  end subroutine get_3darray_double_from_monc

  !> Retreives a 4D array of doubles with a corresponding key from the raw data dump. The size depends on the configuration
  !! that the MONC process supplied to the IO server when it registered which informs of the size of a field.
  !!
  !! Note that this replies in the F2008 `contiguous` keyword, but it is supported by the target compilers so we use it
  !! @param field_starts Field starting locations
  !! @param field_ends Field ending locations
  !! @param data_dump Raw data that MONC has sent to us
  !! @param key Key of the field to retrieve
  !! @param target_data The data that we are writing into
  !! @param size1 Size in first dimension
  !! @param size2 Size in second dimension
  !! @param size3 Size in third dimension
  !! @param size4 Size in fourth dimension
  subroutine get_4darray_double(field_starts, field_ends, data_dump, key, target_data, size1, size2, size3, size4)
    type(map_type), intent(inout) :: field_starts, field_ends
    integer, intent(in) :: size1, size2, size3, size4
    character, dimension(:) :: data_dump
    character(len=*), intent(in) :: key
    real(kind=DOUBLE_PRECISION), dimension(:,:,:,:), pointer, contiguous, intent(inout) :: target_data

    integer :: start_index, end_index, element_size
    real(kind=DOUBLE_PRECISION), dimension(:), pointer :: temp_data
    class(*), pointer :: generic
    
    ! Pointer bounds remapping as transfer needs 1D array but for performance don't want to allocate another array and copy using reshape
    temp_data(1:size1*size2*size3*size4)=>target_data     

    if (.not. c_contains(field_starts, key) .or. .not. c_contains(field_ends, key)) &
         call log_log(LOG_ERROR, "Field name `"//key//"` not found in the data definition")

    generic=>c_get(field_starts, key)
    start_index=conv_to_integer(generic, .false.)
    generic=>c_get(field_ends, key)
    end_index=conv_to_integer(generic, .false.)

    element_size=(end_index-start_index) / kind(target_data)
    
    temp_data=transfer(data_dump(start_index:end_index), temp_data)
  end subroutine get_4darray_double  
  
  !> Retreives a 4D array of doubles with a corresponding key from the raw data dump. The size depends on the configuration
  !! that the MONC process supplied to the IO server when it registered which informs of the size of a field.
  !!
  !! Note that this replies in the F2008 `contiguous` keyword, but it is supported by the target compilers so we use it
  !! @param io_configuration Configuration of the IO server
  !! @param source PID of the MONC process
  !! @param data_id The ID of the data definition that is represented by the dump
  !! @param data_dump Raw data that MONC has sent to us
  !! @param key Key of the field to retrieve
  !! @param target_data The data that we are writing into
  !! @param size1 Size in first dimension
  !! @param size2 Size in second dimension
  !! @param size3 Size in third dimension
  !! @param size4 Size in fourth dimension
  subroutine get_4darray_double_from_monc(io_configuration, source, data_id, data_dump, key, target_data, size1, &
       size2, size3, size4)
    type(io_configuration_type), intent(inout) :: io_configuration
    integer, intent(in) :: source, data_id, size1, size2, size3, size4
    character, dimension(:) :: data_dump
    character(len=*), intent(in) :: key
    real(kind=DOUBLE_PRECISION), dimension(:,:,:,:), pointer, contiguous, intent(inout) :: target_data

    integer :: monc_location
    class(*), pointer :: generic
    
    generic=>c_get(io_configuration%monc_to_index, conv_to_string(source))
    monc_location=conv_to_integer(generic, .false.)

    call get_4darray_double(io_configuration%registered_moncs(monc_location)%field_start_locations(data_id), &
         io_configuration%registered_moncs(monc_location)%field_end_locations(data_id), data_dump, key, target_data, &
         size1, size2, size3, size4)    
  end subroutine get_4darray_double_from_monc  
end module data_utils_mod
