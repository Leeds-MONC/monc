!> Determines the current stepping direction, which can be either forward or centred. This is mainly for field stepping,
!! which is u, v, w fields but also scalars as well which is th and q.
module steppingdirection_mod
  use monc_component_mod, only : component_descriptor_type
  use state_mod, only : model_state_type, CENTRED_STEPPING, FORWARD_STEPPING
  implicit none

#ifndef TEST_MODE
  private
#endif

  public steppingdirection_get_descriptor

contains

  !> Returns the descriptor of this component
  !! @returns The stepping direction component descriptor
  type(component_descriptor_type) function steppingdirection_get_descriptor()
    steppingdirection_get_descriptor%name="stepping_direction"
    steppingdirection_get_descriptor%version=0.1
    steppingdirection_get_descriptor%initialisation=>initialisation_callback
    steppingdirection_get_descriptor%timestep=>timestep_callback
  end function steppingdirection_get_descriptor

  !> Sets the scalar stepping on initialisation. This does not change throughout the model run so we can safely set it here
  !! @param current_state The current model state
  subroutine initialisation_callback(current_state)
    type(model_state_type), target, intent(inout) :: current_state

    ! TODO - allow configuration to use forward stepping for scalars (th and q) and momentum (u,v,w)
    current_state%scalar_stepping = CENTRED_STEPPING
    current_state%momentum_stepping = CENTRED_STEPPING
  end subroutine initialisation_callback  

  !> Determines whether we are forward or centre stepping
  !!
  !! This is important as forward stepping will effectively ignore the Z field
  !! and performs no smoothing on the field. Centre stepping uses the Z field and
  !! the Robert filter for smoothing.
  !! @param current_state The current model state
  subroutine timestep_callback(current_state)
    type(model_state_type), target, intent(inout) :: current_state
    
    current_state%field_stepping = CENTRED_STEPPING
    if (current_state%timestep .eq. current_state%start_timestep) current_state%field_stepping = FORWARD_STEPPING
  end subroutine timestep_callback  
end module steppingdirection_mod
