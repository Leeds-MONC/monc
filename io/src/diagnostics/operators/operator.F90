!> Operator federator which manages the different operators which are available. Operators take in any number of scalar reals
!! and output a single scalar real.
module operator_mod
  use datadefn_mod, only : DEFAULT_PRECISION
  use configuration_parser_mod, only : io_configuration_type
  use collections_mod, only : hashmap_type, map_type, list_type
  use arithmetic_operator_mod, only : initialise_arithmetic_operator, finalise_arithmetic_operator, &
       arithmetic_operator_get_required_fields, perform_arithmetic_operator
  use reductionlocation_operator_mod, only : perform_reductionlocation_operator, reductionlocation_operator_get_required_fields
  use localreduce_operator_mod, only : perform_localreduce_operator, localreduce_operator_get_required_fields
  use fieldslicer_operator_mod, only : perform_fieldslicer_operator, fieldslicer_operator_get_required_fields
  use fieldcoarsener_operator_mod, only : perform_fieldcoarsener_operator, fieldcoarsener_operator_get_required_fields, &
       fieldcoarsener_operator_get_auto_size
  use logging_mod, only : LOG_ERROR, log_log
  implicit none

#ifndef TEST_MODE
  private
#endif

  abstract interface
     type(list_type) function get_operator_required_fields_interface(action_attributes)
       import :: list_type, map_type
       type(map_type), intent(inout) :: action_attributes
     end function get_operator_required_fields_interface
     subroutine perform_activity(io_configuration, field_values, action_attributes, source_monc_location, &
          source_monc, operator_result_values)       
       import DEFAULT_PRECISION, hashmap_type, map_type, io_configuration_type
       type(io_configuration_type), intent(inout) :: io_configuration
       type(hashmap_type), intent(inout) :: field_values
       type(map_type), intent(inout) :: action_attributes
       integer, intent(in) :: source_monc_location, source_monc
       real(kind=DEFAULT_PRECISION), dimension(:), allocatable, intent(inout) :: operator_result_values
     end subroutine perform_activity
  end interface

  public perform_activity, initialise_operators, finalise_operators, get_operator_required_fields, &
       get_operator_perform_procedure, get_operator_auto_size
contains

  !> Initialises any operators that require initialisation
  subroutine initialise_operators()
    call initialise_arithmetic_operator()
  end subroutine initialise_operators

  !> Finalises any operators that require finalisation
  subroutine finalise_operators()
    call finalise_arithmetic_operator()
  end subroutine finalise_operators 

  !> Retrieves the operator execution procedure of an operator with a specific name
  !! @param operator_name The name of the operator to retrieve the execution procedure of
  !! @returns The execution procedure which can be called to run the operator
  function get_operator_perform_procedure(operator_name)
    character(len=*), intent(in) :: operator_name
    procedure(perform_activity), pointer :: get_operator_perform_procedure

    if (trim(operator_name) .eq. "arithmetic") then
      get_operator_perform_procedure=>perform_arithmetic_operator
    else if (trim(operator_name) .eq. "localreduce") then
      get_operator_perform_procedure=>perform_localreduce_operator
    else if (trim(operator_name) .eq. "reductionlocation") then
      get_operator_perform_procedure=>perform_reductionlocation_operator
    else if (trim(operator_name) .eq. "field_slicer") then
      get_operator_perform_procedure=>perform_fieldslicer_operator
    else if (trim(operator_name) .eq. "field_coarsener") then
      get_operator_perform_procedure=>perform_fieldcoarsener_operator
    else
      call log_log(LOG_ERROR, "Operator '"//trim(operator_name)//"' not found so ignoring")
    end if    
  end function get_operator_perform_procedure  
  
  !> Retrieves the list of fields required by an operator before it can run
  !! @param operator_name The name of the operator to retrieve the execution procedure of
  !! @param action_attributes The attributes of the action as defined by the configuration
  type(list_type) function get_operator_required_fields(operator_name, action_attributes)
    character(len=*), intent(in) :: operator_name
    type(map_type), intent(inout) :: action_attributes

    if (trim(operator_name) .eq. "arithmetic") then
      get_operator_required_fields=arithmetic_operator_get_required_fields(action_attributes)
    else if (trim(operator_name) .eq. "localreduce") then
      get_operator_required_fields=localreduce_operator_get_required_fields(action_attributes)
    else if (trim(operator_name) .eq. "reductionlocation") then
      get_operator_required_fields=reductionlocation_operator_get_required_fields(action_attributes)
    else if (trim(operator_name) .eq. "field_slicer") then
      get_operator_required_fields=fieldslicer_operator_get_required_fields(action_attributes)
    else if (trim(operator_name) .eq. "field_coarsener") then
      get_operator_required_fields=fieldcoarsener_operator_get_required_fields(action_attributes)
    else
      call log_log(LOG_ERROR, "Operator '"//trim(operator_name)//"' not found so ignoring")
    end if    
  end function get_operator_required_fields

  integer function get_operator_auto_size(io_configuration, operator_name, auto_dimension, action_attributes)
    type(io_configuration_type), intent(inout) :: io_configuration
    character(len=*), intent(in) :: operator_name, auto_dimension
    type(map_type), intent(inout) :: action_attributes

    if (trim(operator_name) .eq. "field_coarsener") then
      get_operator_auto_size=fieldcoarsener_operator_get_auto_size(io_configuration, auto_dimension, action_attributes)
    else
      get_operator_auto_size=-1
    end if    
  end function get_operator_auto_size  
end module operator_mod
