!> Coarsens a field by selecting data with a specific period in any number of dimensions
module fieldcoarsener_operator_mod
  use datadefn_mod, only : DEFAULT_PRECISION, STRING_LENGTH
  use configuration_parser_mod, only : io_configuration_type, io_configuration_field_type, io_configuration_registered_monc_type, &
       data_values_type, get_data_value_by_field_name, get_prognostic_field_configuration
  use data_utils_mod, only : get_action_attribute_string, get_action_attribute_integer
  use grids_mod, only : Z_INDEX, Y_INDEX, X_INDEX
  use collections_mod, only : hashmap_type, list_type, map_type, c_add_string, c_get_integer
  use conversions_mod, only : conv_to_integer
  use logging_mod, only : LOG_ERROR, log_log
  implicit none

#ifndef TEST_MODE
  private
#endif

  public perform_fieldcoarsener_operator, fieldcoarsener_operator_get_required_fields, fieldcoarsener_operator_get_auto_size
contains

  !> Retrieves the size of an auto dimension based upon the work that will be completed here
  !! @param io_configuration The IO server configuration
  !! @param auto_dimension String of the auto dimension to look up
  !! @param action_attributes The XML configuration attributes
  !! @returns The size of an auto dimension
  integer function fieldcoarsener_operator_get_auto_size(io_configuration, auto_dimension, action_attributes)
    type(io_configuration_type), intent(inout) :: io_configuration
    character(len=*), intent(in) :: auto_dimension
    type(map_type), intent(inout) :: action_attributes

    integer :: auto_dim_id, index_of_match, entire_dim_size
    integer, dimension(:), allocatable :: dimensions_to_slice, indexes_to_slice

    auto_dim_id=convert_dimension_str_to_id(auto_dimension)
    call get_dimensions_and_indexes_to_slice(get_action_attribute_string(action_attributes, "dimension"), &
           get_action_attribute_string(action_attributes, "period"), dimensions_to_slice, indexes_to_slice)
    index_of_match=locate_dimension(auto_dim_id, dimensions_to_slice)
    if (index_of_match .ne. 0) then
      entire_dim_size=get_entire_dimension_size(io_configuration, auto_dimension)
      fieldcoarsener_operator_get_auto_size=ceiling(real(entire_dim_size)/indexes_to_slice(index_of_match))
    else
      fieldcoarsener_operator_get_auto_size=-1
    end if
    deallocate(dimensions_to_slice, indexes_to_slice)
  end function fieldcoarsener_operator_get_auto_size  

  !> Performs the field coarsener operator on a specific field
  !! @param io_configuration The IO server configuration
  !! @param field_values The field values
  !! @param action_attributes Attributes associated with the running of this operator
  !! @param source_monc_location The location in configuration data of MONC settings
  !! @param source_monc Process ID of the MONC that sent this data
  !! @param operator_result_values The resulting value or left unallocated if none are appropriate
  subroutine perform_fieldcoarsener_operator(io_configuration, field_values, action_attributes, source_monc_location, &
       source_monc, operator_result_values)
    type(io_configuration_type), intent(inout) :: io_configuration
    type(hashmap_type), intent(inout) :: field_values
    type(map_type), intent(inout) :: action_attributes
    integer, intent(in) :: source_monc_location, source_monc
    real(kind=DEFAULT_PRECISION), dimension(:), allocatable, intent(inout) :: operator_result_values

    character(len=STRING_LENGTH) :: field_to_slice
    integer, dimension(:), allocatable :: dimensions_to_slice, indexes_to_slice, dim_weights, dim_periods, dim_starts, indexes
    integer :: number_dims, sliced_size, i, j, source_dim
    type(io_configuration_field_type) :: corresponding_field_definition
    type(data_values_type), pointer :: field_local_values

    field_to_slice=get_action_attribute_string(action_attributes, "field")
    ! NSE
    if (get_prognostic_field_configuration(io_configuration, field_to_slice, "", corresponding_field_definition)) then
      call get_dimensions_and_indexes_to_slice(get_action_attribute_string(action_attributes, "dimension"), &
           get_action_attribute_string(action_attributes, "period"), dimensions_to_slice, indexes_to_slice)

      call determine_dimension_bounds(corresponding_field_definition, io_configuration%registered_moncs(source_monc_location), &
             dimensions_to_slice, indexes_to_slice, dim_weights, dim_periods, dim_starts, number_dims, sliced_size)

      field_local_values=>get_data_value_by_field_name(field_values, field_to_slice)

      allocate(operator_result_values(sliced_size), indexes(number_dims))
      indexes=dim_starts
      do i=1, sliced_size
        source_dim=1
        do j=1, number_dims
          source_dim=source_dim+(indexes(j)-1)*dim_weights(j)
        end do
        operator_result_values(i)=field_local_values%values(source_dim)
        indexes(1)=indexes(1)+dim_periods(1)
        do j=1, number_dims
          if (indexes(j) .gt. io_configuration%registered_moncs(source_monc_location)%local_dim_sizes(j)) then
            indexes(j)=1
            if (j .lt. number_dims) indexes(j+1)=indexes(j+1)+dim_periods(j+1)
          end if
        end do
      end do
      deallocate(dimensions_to_slice, indexes_to_slice, dim_weights, dim_periods, indexes, dim_starts)
    end if
  end subroutine perform_fieldcoarsener_operator

  !> Locates a dimension in a list of dimensions or 0 if none can be found
  !! @param dimension_id The id to look up
  !! @param list_of_dims The list of dimensions to check
  !! @param The index that this dimension resides at or 0 if it is not found
  integer function locate_dimension(dimension_id, list_of_dims)
    integer, intent(in) :: dimension_id, list_of_dims(:)

    integer :: i

    do i=1, size(list_of_dims)
      if (list_of_dims(i) == dimension_id) then
        locate_dimension=i
        return
      end if
    end do
    locate_dimension=0
  end function locate_dimension  

  !> Determines the dimension bounds which are used in the slicing
  !! @param corresponding_field_definition The definition of the field that is being sliced
  !! @param registered_monc_info This specific MONC process
  !! @param dimensions_to_slice Numeric code for the dimensions that should be sliced
  !! @param indexes_to_slice The indexes that we are slicing upon
  !! @param dim_weights Output weight for each dimension when dealing with a 1D data
  !! @param dim_periods The periods that we jump over in the data to grab the next source point
  !! @param dim_starts The local start point in data
  !! @param number_dims Output number of dimensions that the field comprises of
  !! @param sliced_size Output overall 1D data size of the sliced array
  subroutine determine_dimension_bounds(corresponding_field_definition, registered_monc_info, dimensions_to_slice, &
       indexes_to_slice, dim_weights, dim_periods, dim_starts, number_dims, sliced_size)
    type(io_configuration_field_type), intent(in) :: corresponding_field_definition
    type(io_configuration_registered_monc_type), intent(in) :: registered_monc_info
    integer, dimension(:), intent(in) :: dimensions_to_slice, indexes_to_slice
    integer, dimension(:), allocatable, intent(out) :: dim_weights, dim_periods, dim_starts
    integer, intent(out) :: number_dims, sliced_size

    integer :: i, j, dimension_id, amount_to_add, located_dim
    integer, dimension(:), allocatable :: dim_sizes
    logical, dimension(:), allocatable :: found_slice_field

    number_dims=corresponding_field_definition%dimensions
    allocate(dim_sizes(number_dims), dim_weights(number_dims), dim_periods(number_dims), dim_starts(number_dims),&
         found_slice_field(size(dimensions_to_slice)))
    found_slice_field=.false.

    do i=1, number_dims
      dimension_id=convert_dimension_str_to_id(corresponding_field_definition%dim_size_defns(i))
      located_dim=locate_dimension(dimension_id, dimensions_to_slice)      
      if (located_dim .gt. 0) then        
        dim_periods(i)=indexes_to_slice(located_dim)
        if (registered_monc_info%local_dim_starts(dimension_id) == 1) then
          dim_starts(i)=1
        else
          dim_starts(i)=dim_periods(i) - mod(registered_monc_info%local_dim_starts(dimension_id)-2, dim_periods(i))
        end if
        dim_sizes(i)=ceiling(real(registered_monc_info%local_dim_sizes(dimension_id) - (dim_starts(i)-1))/&
             real(indexes_to_slice(located_dim)))
        found_slice_field(located_dim)=.true.
      else
        dim_sizes(i)=registered_monc_info%local_dim_sizes(dimension_id)
        dim_periods(i)=1
        dim_starts(i)=1
      end if
    end do
    do i=1, size(found_slice_field)
      if (.not. found_slice_field(i)) call log_log(LOG_ERROR, "Can not find a dimension to slice in provided field")
    end do
    sliced_size=0
    do i=1, number_dims
      amount_to_add=1
      dim_weights(i)=1
      do j=1, i
        if (j .lt. i) dim_weights(i)=dim_weights(i)*registered_monc_info%local_dim_sizes(j)
        amount_to_add=amount_to_add*dim_sizes(j)
      end do
      sliced_size=sliced_size+amount_to_add
    end do
    deallocate(dim_sizes, found_slice_field)
  end subroutine determine_dimension_bounds

  !> Retrieves the dimensions and indexes to slice from the strings provided in configuration
  !! @param str_dim_to_slice The string dimensions that we are going to slice
  !! @param str_index_to_slice The string indexes that we are going to slice
  !! @param dimensions_to_slice Numeric corresponding codes for the dimensions to slice
  !! @param indexes_to_slice Numeric corresponding codes for the indexes to slice
  subroutine get_dimensions_and_indexes_to_slice(str_dim_to_slice, str_index_to_slice, dimensions_to_slice, indexes_to_slice)
    character(len=*), intent(in) :: str_dim_to_slice, str_index_to_slice
    integer, dimension(:), allocatable, intent(out) :: dimensions_to_slice, indexes_to_slice

    integer :: num_dims, num_indexes, i, dim_idx, index_idx, idx
    num_dims=get_occurances_of_character(str_dim_to_slice, ",")+1
    num_indexes=get_occurances_of_character(str_index_to_slice, ",")+1
    if (num_dims .ne. num_indexes) then
      call log_log(LOG_ERROR, "Coarsening number of dimensions and indexes must match")
    end if
    allocate(dimensions_to_slice(num_dims), indexes_to_slice(num_indexes))
    dim_idx=1
    index_idx=1
    do i=1, num_dims
      idx=index(str_dim_to_slice(dim_idx:), ",")
      idx=idx-1
      if (idx == -1) then
        idx=len_trim(str_dim_to_slice)
      else
        idx=idx+(dim_idx-1)
      end if
      dimensions_to_slice(i)=convert_dimension_str_to_id(str_dim_to_slice(dim_idx:idx))
      dim_idx=idx+2
      
      idx=index(str_index_to_slice(index_idx:), ",")
      idx=idx-1
      if (idx == -1) then
        idx=len_trim(str_index_to_slice)
      else
        idx=idx+(index_idx-1)
      end if
      indexes_to_slice(i)=conv_to_integer(str_index_to_slice(index_idx:idx))
      index_idx=idx+2
    end do
  end subroutine get_dimensions_and_indexes_to_slice 

  !> Retrieves the number of occurances of a character in a source string
  !! @param source_str The source string to search
  !! @param search_char The character that we are searching for
  !! @returns The number of occurances
  integer function get_occurances_of_character(source_str, search_char)
    character(len=*), intent(in) :: source_str, search_char

    integer :: i, n, idx

    get_occurances_of_character=0
    n=len_trim(source_str)
    i=1
    do while (i .le. n)
      idx=index(source_str(i:), search_char)
      if (idx == 0) exit
      i=i+idx+1
      get_occurances_of_character=get_occurances_of_character+1
    end do    
  end function get_occurances_of_character 

  !> Retrieves a list of the required fields for running this operator
  !! @param action_attributes The attributes which configure the operator
  !! @returns A list of required fields before the operator can run
  type(list_type) function fieldcoarsener_operator_get_required_fields(action_attributes)
    type(map_type), intent(inout) :: action_attributes

    character(len=STRING_LENGTH) :: field_to_slice

    field_to_slice=get_action_attribute_string(action_attributes, "field")
    call c_add_string(fieldcoarsener_operator_get_required_fields, field_to_slice)
  end function fieldcoarsener_operator_get_required_fields

  !> Converts a dimension string to the corresponding numeric ID
  !! @param dim_str The dimension string
  !! @returns The corresponding numeric dimension ID
  integer function convert_dimension_str_to_id(dim_str)
    character(len=*), intent(in) :: dim_str

    if (dim_str .eq. "x") then
      convert_dimension_str_to_id=X_INDEX
    else if (dim_str .eq. "y") then
      convert_dimension_str_to_id=Y_INDEX
    else if (dim_str .eq. "z") then
      convert_dimension_str_to_id=Z_INDEX
    else
      call log_log(LOG_ERROR, "Can not translate dimension "//trim(dim_str))
    end if
  end function convert_dimension_str_to_id  

  !> Looks up the global size of a dimension based upon its name
  !! @param io_configuration The IO server configuration
  !! @param dimension_name The name of the dimension that we are looking up
  !! @returns The global size of this dimension
  integer function get_entire_dimension_size(io_configuration, dimension_name)
    type(io_configuration_type), intent(inout) :: io_configuration
    character(len=STRING_LENGTH), intent(in) :: dimension_name

    get_entire_dimension_size=c_get_integer(io_configuration%dimension_sizing, dimension_name)
  end function get_entire_dimension_size
end module fieldcoarsener_operator_mod
