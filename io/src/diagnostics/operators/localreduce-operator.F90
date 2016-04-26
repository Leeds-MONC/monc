!> Performs a local reduction, reducing a local array into a single scalar value
module localreduce_operator_mod
  use datadefn_mod, only : DEFAULT_PRECISION, STRING_LENGTH
  use configuration_parser_mod, only : io_configuration_type, data_values_type, get_data_value_by_field_name
  use data_utils_mod, only : get_action_attribute_string
  use collections_mod, only : hashmap_type, list_type, map_type, c_add_string  
  implicit none

#ifndef TEST_MODE
  private
#endif

  public perform_localreduce_operator, localreduce_operator_get_required_fields
contains

  !> Executes this local reduction operator
  !! @param io_configuration Configuration of the IO server  
  !! @param field_values The field values
  !! @param action_attributes Attributes associated with the running of this operator
  !! @param source_monc_location Location of the source MONC
  !! @param source_monc The source MONC
  !! @param operator_result_values The operator resulting (scalar) value
  subroutine perform_localreduce_operator(io_configuration, field_values, action_attributes, source_monc_location, &
       source_monc, operator_result_values)
    type(io_configuration_type), intent(inout) :: io_configuration
    type(hashmap_type), intent(inout) :: field_values
    type(map_type), intent(inout) :: action_attributes
    integer, intent(in) :: source_monc_location, source_monc
    real(kind=DEFAULT_PRECISION), dimension(:), allocatable, intent(inout) :: operator_result_values

    character(len=STRING_LENGTH) :: field_to_reduce, reduction_operator
    type(data_values_type), pointer :: field_local_values

    allocate(operator_result_values(1))

    field_to_reduce=get_action_attribute_string(action_attributes, "field")
    reduction_operator=get_action_attribute_string(action_attributes, "operator")
    field_local_values=>get_data_value_by_field_name(field_values, field_to_reduce)
    operator_result_values=do_local_reduction(field_local_values%values, reduction_operator)
  end subroutine perform_localreduce_operator

  !> Does the actual local reduction, translating the array into a vector based upon the operator
  !! @param data The array data to reduce
  !! @param reduction_operator The operator to apply for this reduction
  !! @returns The resulting scalar value
  real(kind=DEFAULT_PRECISION) function do_local_reduction(data, reduction_operator)
    real(kind=DEFAULT_PRECISION), dimension(:), intent(in) :: data
    character(len=STRING_LENGTH), intent(in) :: reduction_operator

    if (reduction_operator .eq. "max") then
      do_local_reduction=maxval(data)
    else if (reduction_operator .eq. "min") then
      do_local_reduction=minval(data)
    else if (reduction_operator .eq. "sum") then
      do_local_reduction=sum(data)
    else if (reduction_operator .eq. "mean") then
      do_local_reduction=sum(data)/size(data)
    end if    
  end function do_local_reduction  

  !> Retrieves the list of fields needed by this operator for a specific configuration
  !! @param action_attributes The attributes which configure the operator
  !! @returns A list of required fields before the operator can run
  type(list_type) function localreduce_operator_get_required_fields(action_attributes)
    type(map_type), intent(inout) :: action_attributes

    character(len=STRING_LENGTH) :: field_to_reduce

    field_to_reduce=get_action_attribute_string(action_attributes, "field")
    call c_add_string(localreduce_operator_get_required_fields, field_to_reduce)
  end function localreduce_operator_get_required_fields
end module localreduce_operator_mod
