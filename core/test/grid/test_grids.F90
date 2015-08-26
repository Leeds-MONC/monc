! Tests the logging_mod utility functions
module test_grid_mod

  use fruit, only : assert_equals, assert_not_equals, assert_true
  use datadefn_mod, only : DEFAULT_PRECISION
  use grids_mod, only :  X_INDEX, Y_INDEX, Z_INDEX, PRIMAL_GRID, DUAL_GRID, &
       global_grid_type, local_grid_type
  implicit none
  
contains
  
  ! Test index parameters
  subroutine test_XYZ_INDEX
    call assert_equals(1, Z_INDEX, "Z_INDEX is 1")
    call assert_equals(2, Y_INDEX, "Y_INDEX is 2")
    call assert_equals(3, X_INDEX, "X_INDEX is 3")
  end subroutine test_XYZ_INDEX

  ! Test type parameters
  subroutine test_type_parameters
    call assert_equals(1, PRIMAL_GRID, "PRIMAL_GRID is 1")
    call assert_equals(2, DUAL_GRID, "DUAL_GRID is 2")
  end subroutine test_type_parameters

  ! Test dimensions status of global_grid
  subroutine test_global_active_dimensions
    type(global_grid_type) :: global_grid
    logical, dimension(3) :: expected = (/ .true., .true., .true. /) 
    integer :: i
    do i = 1,3
       call assert_not_equals(expected(i), global_grid%active(i), &
            "Test the status of the dimensions after creating a global_grid_type object")
    end do
end subroutine test_global_active_dimensions

  ! Test dimensions status of local_grid
  subroutine test_local_active_dimensions
    type(local_grid_type) :: local_grid
    logical, dimension(3) :: expected = (/ .true., .true., .true. /) 
    integer :: i
    do i = 1,3
       call assert_not_equals(expected(i), local_grid%active(i), &
            "Test the status of the dimensions after creating a local_grid_type object")
    end do
end subroutine test_local_active_dimensions

end module test_grid_mod

  

  ! Driver for maths_mod utility tests
  program test_grid_driver
    use fruit, only : init_fruit, run_test_case, fruit_summary
    use test_grid_mod, only : test_XYZ_INDEX, test_type_parameters, test_global_active_dimensions, &
         test_local_active_dimensions
    
    implicit none
    
    call init_fruit
    call run_test_case(test_XYZ_INDEX, "Test expected XYZ index values")
    call run_test_case(test_type_parameters, "Test expected grid type parameters values")
    call run_test_case(test_global_active_dimensions, "Test dimension status when difining global_type object")
    call run_test_case(test_local_active_dimensions, "Test dimension status when difining local_type object")
    call fruit_summary
  end program test_grid_driver
