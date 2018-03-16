!> Dummy stub when not compiling with CASIM microphysics
module casim_profile_dgs_mod
  use monc_component_mod, only : component_descriptor_type
  implicit none

#ifndef TEST_MODE
  private
#endif

  public casim_profile_dgs_get_descriptor
contains

  !> Provides the descriptor back to the caller and is used in component registration
  !! @returns The termination check component descriptor
  type(component_descriptor_type) function casim_profile_dgs_get_descriptor()
    casim_profile_dgs_get_descriptor%name="casim_profile_dgs"
    casim_profile_dgs_get_descriptor%version=0.0
  end function casim_profile_dgs_get_descriptor

end module casim_profile_dgs_mod
