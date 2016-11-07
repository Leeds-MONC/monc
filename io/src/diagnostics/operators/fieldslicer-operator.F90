!> Slices a field based upon the selected dimension and index
module fieldslicer_operator_mod
  use datadefn_mod, only : DEFAULT_PRECISION, STRING_LENGTH
  use configuration_parser_mod, only : io_configuration_type, io_configuration_field_type, io_configuration_registered_monc_type, &
       data_values_type, get_data_value_by_field_name, get_prognostic_field_configuration
  use data_utils_mod, only : get_action_attribute_string, get_action_attribute_integer
  use grids_mod, only : Z_INDEX, Y_INDEX, X_INDEX
  use collections_mod, only : hashmap_type, list_type, map_type, c_add_string
  use logging_mod, only : LOG_ERROR, log_log
  implicit none

#ifndef TEST_MODE
  private
#endif

  public perform_fieldslicer_operator, fieldslicer_operator_get_required_fields
contains
  !> Performs the actual field slicing
  !! @param io_configuration Configuration of the IO server  
  !! @param field_values The field values
  !! @param action_attributes Attributes associated with the running of this operator
  !! @param source_monc_location The location in configuration data of MONC settings
  !! @param source_monc Process ID of the MONC that sent this data
  !! @param operator_result_values The resulting value or left unallocated if none are appropriate
  subroutine perform_fieldslicer_operator(io_configuration, field_values, action_attributes, source_monc_location, &
       source_monc, operator_result_values)
    type(io_configuration_type), intent(inout) :: io_configuration
    type(hashmap_type), intent(inout) :: field_values
    type(map_type), intent(inout) :: action_attributes
    integer, intent(in) :: source_monc_location, source_monc
    real(kind=DEFAULT_PRECISION), dimension(:), allocatable, intent(inout) :: operator_result_values

    character(len=STRING_LENGTH) :: field_to_slice
    integer :: i, j, dimension_to_slice, index_to_slice, number_dims, sliced_size, source_dim
    integer, dimension(:), allocatable :: dim_starts, dim_ends, dim_weights, indexes
    type(data_values_type), pointer :: field_local_values
    type(io_configuration_field_type) :: corresponding_field_definition

    field_to_slice=get_action_attribute_string(action_attributes, "field")
    ! NSE
    if (get_prognostic_field_configuration(io_configuration, field_to_slice, "", corresponding_field_definition)) then      
      dimension_to_slice=get_dimension_to_slice(action_attributes)
      index_to_slice=get_action_attribute_integer(action_attributes, "index")

      if (io_configuration%registered_moncs(source_monc_location)%local_dim_starts(dimension_to_slice) .le. index_to_slice .and.&
           io_configuration%registered_moncs(source_monc_location)%local_dim_ends(dimension_to_slice) .ge. index_to_slice) then
        field_local_values=>get_data_value_by_field_name(field_values, field_to_slice)      

        call determine_dimension_bounds(corresponding_field_definition, io_configuration%registered_moncs(source_monc_location), &
             dimension_to_slice, index_to_slice, dim_starts, dim_ends, dim_weights, number_dims, sliced_size)

        allocate(operator_result_values(sliced_size), indexes(number_dims))
        indexes=dim_starts
        do i=1, sliced_size
          source_dim=1
          do j=1, number_dims
            source_dim=source_dim+(indexes(j)-1)*dim_weights(j)
          end do
          operator_result_values(i)=field_local_values%values(source_dim)
          indexes(1)=indexes(1)+1
          do j=1, number_dims
            if (indexes(j) .gt. dim_ends(j)) then
              indexes(j)=dim_starts(j)
              if (j .lt. number_dims) indexes(j+1)=indexes(j+1)+1
            end if
          end do
        end do
        deallocate(dim_starts, dim_ends, dim_weights, indexes)
      end if
    end if    
  end subroutine perform_fieldslicer_operator 

  !> Determines the bounds for each dimension which is specified via the configuration file
  !! @param corresponding_field_definition The definition of the field that is being sliced
  !! @param registered_monc_info This specific MONC process
  !! @param dimension_to_slice Numeric code for the dimension that should be sliced
  !! @param index_to_slice The index that we are slicing upon
  !! @param dim_starts Output start index for each dimension
  !! @param dim_ends Output end index for each dimension
  !! @param dim_weights Output weight for each dimension when dealing with a 1D data
  !! @param number_dims Output number of dimensions that the field comprises of
  !! @param sliced_size Output overall 1D data size of the sliced array
  subroutine determine_dimension_bounds(corresponding_field_definition, registered_monc_info, dimension_to_slice, &
       index_to_slice, dim_starts, dim_ends, dim_weights, number_dims, sliced_size)
    type(io_configuration_field_type), intent(in) :: corresponding_field_definition
    type(io_configuration_registered_monc_type), intent(in) :: registered_monc_info
    integer, intent(in) :: dimension_to_slice, index_to_slice
    integer, dimension(:), allocatable, intent(out) :: dim_starts, dim_ends, dim_weights
    integer, intent(out) :: number_dims, sliced_size

    integer :: i, j, dimension_id, amount_to_add
    integer, dimension(:), allocatable :: dim_sizes
    logical :: found_slice_field

    number_dims=corresponding_field_definition%dimensions
    allocate(dim_sizes(number_dims), dim_starts(number_dims), dim_ends(number_dims), dim_weights(number_dims))
    found_slice_field=.false.

    do i=1, number_dims
      dimension_id=convert_dimension_str_to_id(corresponding_field_definition%dim_size_defns(i))
      if (dimension_id==dimension_to_slice) then
        dim_sizes(i)=1
        dim_starts(i)=index_to_slice
        dim_ends(i)=index_to_slice
        found_slice_field=.true.
      else
        dim_sizes(i)=registered_monc_info%local_dim_sizes(dimension_id)
        dim_starts(i)=1
        dim_ends(i)=registered_monc_info%local_dim_sizes(dimension_id)
      end if
    end do
    if (.not. found_slice_field) call log_log(LOG_ERROR, "Can not find dimension to slice in provided field")
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
    deallocate(dim_sizes)
  end subroutine determine_dimension_bounds

  !> Retrieves a list of the required fields for running this operator
  !! @param action_attributes The attributes which configure the operator
  !! @returns A list of required fields before the operator can run
  type(list_type) function fieldslicer_operator_get_required_fields(action_attributes)
    type(map_type), intent(inout) :: action_attributes

    character(len=STRING_LENGTH) :: field_to_slice

    field_to_slice=get_action_attribute_string(action_attributes, "field")
    call c_add_string(fieldslicer_operator_get_required_fields, field_to_slice)
  end function fieldslicer_operator_get_required_fields

  !> Retrieves the integer index of the dimension to slice
  !! @param action_attributes The attributes which configure the operator
  integer function get_dimension_to_slice(action_attributes)
    type(map_type), intent(inout) :: action_attributes

    get_dimension_to_slice=convert_dimension_str_to_id(get_action_attribute_string(action_attributes, "dimension"))
  end function get_dimension_to_slice

  !> Converts a dimension string to the corresponding numeric ID
  !! @param dim_str The dimension string
  !! @returns The corresponding numeric dimension ID
  integer function convert_dimension_str_to_id(dim_str)
    character(len=STRING_LENGTH), intent(in) :: dim_str

    if (dim_str .eq. "x") then
      convert_dimension_str_to_id=X_INDEX
    else if (dim_str .eq. "y") then
      convert_dimension_str_to_id=Y_INDEX
    else if (dim_str .eq. "z") then
      convert_dimension_str_to_id=Z_INDEX
    end if
  end function convert_dimension_str_to_id  
end module fieldslicer_operator_mod
