!> Clears the source terms at the start of a timestep to then be populated by items in the dynamics group
module clearsourceterms_mod
  use datadefn_mod, only : DEFAULT_PRECISION
  use state_mod, only : model_state_type
  use monc_component_mod, only : COMPONENT_ARRAY_FIELD_TYPE, COMPONENT_DOUBLE_DATA_TYPE, component_descriptor_type, &
       component_field_value_type, component_field_information_type
implicit none

#ifndef TEST_MODE
  private
#endif

  public clearsourceterms_get_descriptor
contains
  !> Provides a description of this component for the core to register
  !! @returns The descriptor containing registration information for this component
  type(component_descriptor_type) function clearsourceterms_get_descriptor()
    clearsourceterms_get_descriptor%name="clearsourceterms"
    clearsourceterms_get_descriptor%version=0.1
    clearsourceterms_get_descriptor%timestep=>timestep_callback
  end function clearsourceterms_get_descriptor

  !> Timestep callback which simply clears the source terms
  !! @param current_state The current model state
  subroutine timestep_callback(current_state)
    type(model_state_type), target, intent(inout) :: current_state

    integer :: i

#ifdef U_ACTIVE
    current_state%su%data=0.0_DEFAULT_PRECISION
#endif
#ifdef V_ACTIVE
    current_state%sv%data=0.0_DEFAULT_PRECISION
#endif
#ifdef W_ACTIVE
    current_state%sw%data=0.0_DEFAULT_PRECISION
#endif

    if (current_state%sth%active) then
      current_state%sth%data=0.0_DEFAULT_PRECISION
    end if

    do i=1, current_state%number_q_fields
      current_state%sq(i)%data=0.0_DEFAULT_PRECISION
    end do
    
    if (current_state%n_tracers .gt. 0) then
      do i=1, current_state%n_tracers
        current_state%stracer(i)%data=0.0_DEFAULT_PRECISION
      end do
      if (current_state%traj_tracer_index >0 .and. current_state%reinit_tracer) then
        current_state%reinit_tracer=.false. 
      end if
    end if
  end subroutine timestep_callback  
end module clearsourceterms_mod
