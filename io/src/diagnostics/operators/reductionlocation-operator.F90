module reductionlocation_operator_mod
  use datadefn_mod, only : DEFAULT_PRECISION, STRING_LENGTH
  use configuration_parser_mod, only : io_configuration_type, data_values_type, get_data_value_by_field_name
  use data_utils_mod, only : get_action_attribute_string
  use collections_mod, only : hashmap_type, list_type, map_type, hashmap_type, c_get, c_add
  use conversions_mod, only : conv_to_generic, generic_to_double_real
  implicit none

#ifndef TEST_MODE
  private
#endif

  public perform_reductionlocation_operator, reductionlocation_operator_get_required_fields
contains

  function perform_reductionlocation_operator(io_configuration, field_values, action_attributes)
    type(io_configuration_type), intent(inout) :: io_configuration
    type(hashmap_type), intent(inout) :: field_values
    type(map_type), intent(inout) :: action_attributes
    real(kind=DEFAULT_PRECISION), dimension(:), allocatable :: perform_reductionlocation_operator

    character(len=STRING_LENGTH) :: location_local, val, val_local
    type(data_values_type), pointer :: val_local_values, val_values, location_local_values
    integer :: i

    call extract_tripplet_variables(get_action_attribute_string(action_attributes, "input"), location_local, val, val_local)

    val_local_values=>get_data_value_by_field_name(field_values, val_local)
    val_values=>get_data_value_by_field_name(field_values, val)
    location_local_values=>get_data_value_by_field_name(field_values, location_local)

    allocate(perform_reductionlocation_operator(size(location_local_values%values)))
    do i=1, size(val_local_values%values)
      if (val_local_values%values(i) .eq. val_values%values(i)) then
        perform_reductionlocation_operator(i)=location_local_values%values(i)
      else
        perform_reductionlocation_operator(i)=-1.0_DEFAULT_PRECISION
      end if
    end do
  end function perform_reductionlocation_operator

  type(list_type) function reductionlocation_operator_get_required_fields(action_attributes)
    type(map_type), intent(inout) :: action_attributes

    character(len=STRING_LENGTH) :: location_local, val, val_local

    call extract_tripplet_variables(get_action_attribute_string(action_attributes, "input"), location_local, val, val_local)

    call c_add(reductionlocation_operator_get_required_fields, conv_to_generic(location_local, .true.))
    call c_add(reductionlocation_operator_get_required_fields, conv_to_generic(val, .true.))
    call c_add(reductionlocation_operator_get_required_fields, conv_to_generic(val_local, .true.))
  end function reductionlocation_operator_get_required_fields

  subroutine extract_tripplet_variables(field_str, location_local, val, val_local)
    character(len=*), intent(in) :: field_str
    character(len=*), intent(out) :: location_local, val, val_local

    character :: c
    character(len=STRING_LENGTH) :: artefacts(3)
    integer :: i, field_length, start_point, index_loc

    field_length=len(trim(field_str))

    start_point=1
    index_loc=1
    do i=1, field_length
      c=field_str(i:i)
      if (c .eq. ",") then
        artefacts(index_loc)=trim(adjustl(field_str(start_point: i-1)))
        start_point=i+1
        index_loc=index_loc+1
      end if
    end do
    if (start_point .lt. i) artefacts(index_loc)=trim(adjustl(field_str(start_point: i-1)))
    location_local=artefacts(1)
    val=artefacts(2)
    val_local=artefacts(3)
  end subroutine extract_tripplet_variables  
end module reductionlocation_operator_mod
