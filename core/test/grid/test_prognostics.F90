! Tests the logging_mod utility functions
module test_prognostic_mod
  use fruit, only : assert_equals, assert_not_equals, assert_true
  use mpi, only : MPI_REQUEST_NULL
  use datadefn_mod, only : DEFAULT_PRECISION
  use prognostics_mod, only : prognostic_field_type, prognostic_field_ptr_type, &
       get_field_interpolation_index
  use grids_mod, only : PRIMAL_GRID
  implicit none
  
contains
  
  ! Test async default value
  subroutine test_async_flux_value
    type(prognostic_field_type)  :: field
    call assert_equals(MPI_REQUEST_NULL, field%async_flux_handle , "Test async_flux_handle")    
  end subroutine test_async_flux_value


  ! Test get_field_interpolation_index 
  subroutine test_get_field_interpolation_index
    type(prognostic_field_type) :: field
    logical, dimension(3)  :: get_field_index
    integer :: i
    get_field_index = get_field_interpolation_index(field)
    do i=1,3
       call assert_equals(.false., get_field_index(i) , "Test interpolation index is false")    
    enddo
    
    field%grid(1) = PRIMAL_GRID
    get_field_index = get_field_interpolation_index(field)
    do i=1,3
       if (i == 1) then
          call assert_equals(.true., get_field_index(i) , " Should be true")    
       else
          call assert_equals(.false., get_field_index(i) , "Should be false")    
       endif
    enddo
    
    field%grid(1) = 10
    field%grid(2) = 20
    field%grid(3) = 30
    get_field_index = get_field_interpolation_index(field)
    do i=1,3
       call assert_equals(.false., get_field_index(i) , "Should be false")    
    enddo
    
  end subroutine test_get_field_interpolation_index


  
 end module test_prognostic_mod

  

  ! Driver for maths_mod utility tests
  program test_prognostic_driver
    use fruit, only : init_fruit, run_test_case, fruit_summary
    use maths_mod, only : random
    use test_prognostic_mod, only : test_async_flux_value, test_get_field_interpolation_index
    
    implicit none
    
    call init_fruit
    call run_test_case(test_async_flux_value, "Test async default value")
    call run_test_case(test_get_field_interpolation_index, "Test get_interpolation_index whith empty grid")
    call fruit_summary
  end program test_prognostic_driver
