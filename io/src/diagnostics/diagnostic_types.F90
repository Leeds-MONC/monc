module diagnostic_types_mod
  use datadefn_mod, only : DEFAULT_PRECISION, STRING_LENGTH
  use collections_mod, only : map_type, list_type, c_get_generic
  use operator_mod, only : perform_activity
  use configuration_parser_mod, only : io_configuration_misc_item_type

  !< A diagnostic activity which is executed at some point with an input and returns an output
  type diagnostics_activity_type
     integer :: activity_type, communication_operator, root
     real(kind=DEFAULT_PRECISION) :: result
     type(list_type) :: required_fields
     type(map_type) :: activity_attributes
     character(len=STRING_LENGTH) :: result_name, activity_name, uuid
     procedure(perform_activity), pointer, nopass :: operator_procedure
  end type diagnostics_activity_type

 public diagnostics_activitity_type

contains

  !> Retrieves a misc action from the parsed user XML configuration at a specific index
  !! @param action_members The members to extract from
  !! @param index The index to look up
  !! @returns The misc item at this index or null if none is found
  function get_misc_action_at_index(action_members, index)
    type(list_type), intent(inout) :: action_members
    integer, intent(in) :: index
    type(io_configuration_misc_item_type), pointer :: get_misc_action_at_index

    class(*), pointer :: generic

    generic=>c_get_generic(action_members, index)
    if (associated(generic)) then
      select type(generic)
      type is(io_configuration_misc_item_type)
        get_misc_action_at_index=>generic
      end select
    else
      get_misc_action_at_index=>null()
    end if
  end function get_misc_action_at_index

end module
