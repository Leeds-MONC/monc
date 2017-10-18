module diagnostic_types_mod
  use datadefn_mod, only : DEFAULT_PRECISION, STRING_LENGTH
  use collections_mod, only : map_type, list_type, c_get_generic, iterator_type, c_has_next, &
     c_next_generic, c_get_iterator
  use operator_mod, only : perform_activity
  use configuration_parser_mod, only : io_configuration_misc_item_type

  !< A diagnostic which is a name and then the list of activities require to be executed
  type diagnostics_type
     character(len=STRING_LENGTH) :: diagnostic_name, diagnostic_namespace, uuid
     type(list_type) :: activities
     integer :: generation_timestep_frequency
     logical :: collective
  end type diagnostics_type  

  !< A diagnostic activity which is executed at some point with an input and returns an output
  type diagnostics_activity_type
     integer :: activity_type, communication_operator, root
     real(kind=DEFAULT_PRECISION) :: result
     type(list_type) :: required_fields
     type(map_type) :: activity_attributes
     character(len=STRING_LENGTH) :: result_name, activity_name, uuid
     procedure(perform_activity), pointer, nopass :: operator_procedure
  end type diagnostics_activity_type

  !< The type of activity
  integer, parameter :: OPERATOR_TYPE=1, REDUCTION_TYPE=2, BROADCAST_TYPE=3, ALLREDUCTION_TYPE=4, PERFORM_CLEAN_EVERY=100

 public diagnostics_activitity_type, diagnostics_type, get_diagnostic_activity_by_result_name, &
    retrieve_next_activity, OPERATOR_TYPE, REDUCTION_TYPE, BROADCAST_TYPE, ALLREDUCTION_TYPE, PERFORM_CLEAN_EVERY

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

  !> Retrives a diagnostic activity based upon its result name or null if none is found
  !! @param result_name The name of the result we are looking up
  !! @param diagnostic_entry_index The diagnostic index that we are concerned with
  !! @returns The corresponding activity or null if none is found
  function get_diagnostic_activity_by_result_name(diagnostic_activities, result_name)
    type(list_type) :: diagnostic_activities
    character(len=STRING_LENGTH), intent(in) :: result_name
    type(diagnostics_activity_type), pointer :: get_diagnostic_activity_by_result_name

    type(iterator_type) :: iterator

    iterator=c_get_iterator(diagnostic_activities)
    do while (c_has_next(iterator))
      get_diagnostic_activity_by_result_name=>retrieve_next_activity(iterator)
      if (get_diagnostic_activity_by_result_name%result_name == result_name) then
        return
      end if      
    end do
    get_diagnostic_activity_by_result_name=>null()
  end function get_diagnostic_activity_by_result_name

  !> Retrieves the next activity in a collection being iterated over by an iterator
  !! @param iterator The iterator we are using to iterate over the collection
  !! @returns The next activity or null if none is found
  function retrieve_next_activity(iterator)
    type(iterator_type), intent(inout) :: iterator
    type(diagnostics_activity_type), pointer :: retrieve_next_activity

    class(*), pointer :: generic

    generic=>c_next_generic(iterator)

    if (associated(generic)) then
      select type(generic)
        type is(diagnostics_activity_type)
          retrieve_next_activity=>generic
      end select      
    else
      retrieve_next_activity=>null()
    end if
  end function retrieve_next_activity  

end module
