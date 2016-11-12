!> Dummy stub when not compiling with PETSc iterative solver
module petsc_solver_mod
  use monc_component_mod, only : component_descriptor_type
  implicit none

#ifndef TEST_MODE
  private
#endif

  public petsc_solver_get_descriptor
contains

  !> Provides the descriptor back to the caller and is used in component registration
  !! @returns The termination check component descriptor
  type(component_descriptor_type) function petsc_solver_get_descriptor()
    petsc_solver_get_descriptor%name="petsc_solver"
    petsc_solver_get_descriptor%version=0.0
  end function petsc_solver_get_descriptor
end module petsc_solver_mod
