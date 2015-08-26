! Tests the logging_mod utility functions
module test_interpolation_mod
  use fruit, only : assert_equals, assert_not_equals, &
       assert_true
  use datadefn_mod, only : DEFAULT_PRECISION
  use interpolation_mod, only : piecewise_linear_1d, interpolate_point_linear_1d
  implicit none
  
contains
  
  ! Test random produces a non zero value
  subroutine test_piecewise_zero_nnodes
    
    real(kind=DEFAULT_PRECISION) :: zvals(0), vals(2)
    real(kind=DEFAULT_PRECISION) :: zgrid(2)
    real(kind=DEFAULT_PRECISION) :: field(2)
    integer :: i
    do i=1,2
       vals(i)=10.0_DEFAULT_PRECISION
       zgrid(i)=10.0_DEFAULT_PRECISION
       field(i)=10.0_DEFAULT_PRECISION
    enddo
    call  piecewise_linear_1d(zvals, vals, zgrid, field)
    
    do i=1,2
       call assert_equals(10.0_DEFAULT_PRECISION, field(i), "Test no change in field")
    enddo
  end subroutine  test_piecewise_zero_nnodes
  
 ! Test value obtain when the product result is zero
  subroutine test_piecewise_zgrid_eq_zvals
    
    real(kind=DEFAULT_PRECISION) :: zvals(2), vals(2)
    real(kind=DEFAULT_PRECISION) :: zgrid(2)
    real(kind=DEFAULT_PRECISION) :: field(2)
    integer :: i
    zgrid(1) = 10.0_DEFAULT_PRECISION
    zgrid(2) = 9.0_DEFAULT_PRECISION
    do i=1,2
       vals(i)=50.0_DEFAULT_PRECISION
       zvals(i)=11.0_DEFAULT_PRECISION
    enddo
    zvals(1) = 10.0_DEFAULT_PRECISION
    call  piecewise_linear_1d(zvals, vals, zgrid, field)
    
    call assert_equals(50.0_DEFAULT_PRECISION, field(1), "Test  field(k) = vals(nn-1) when (zgrid(k)=zvals(nn-1)")
    
  end subroutine  test_piecewise_zgrid_eq_zvals
  
 ! Test random produces a non zero value
  subroutine test_interpolate_linear
    
    real(kind=DEFAULT_PRECISION) :: zvals(2), vals(2)
    integer :: i
    real(kind=DEFAULT_PRECISION) ::  z,f1,f2
    character(12) :: extrapolate
    
    extrapolate = 'linear'
    z = 0.0_DEFAULT_PRECISION
    do i=1,2
       vals(i)=50.0_DEFAULT_PRECISION
       zvals(i)=11.0_DEFAULT_PRECISION
    enddo
    zvals(1) = 10.0_DEFAULT_PRECISION
    call interpolate_point_linear_1d(zvals, vals, z, f1)
    call interpolate_point_linear_1d(zvals, vals, z, f2, extrapolate)
    
    call assert_equals(f1, f2, "Test linear is used by default")  
  end subroutine  test_interpolate_linear
  
! Test random produces a non zero value
  subroutine test_no_interpolation
    
    real(kind=DEFAULT_PRECISION) :: zvals(2), vals(2)
    integer :: i
    real(kind=DEFAULT_PRECISION) ::  z,f
    character(12) :: extrapolate
    
    extrapolate = 'linear'
    z = 0.0_DEFAULT_PRECISION
    f = 7.0_DEFAULT_PRECISION
    do i=1,2
       vals(i)=50.0_DEFAULT_PRECISION
       zvals(i)=11.0_DEFAULT_PRECISION
    enddo
    zvals(1) = 0.0_DEFAULT_PRECISION
    call interpolate_point_linear_1d(zvals, vals, z, f)

    
    call assert_equals(7.0_DEFAULT_PRECISION, f, "Test linear is used by default")  
  end subroutine  test_no_interpolation

  
end module test_interpolation_mod

  

  ! Driver for maths_mod utility tests
  program test_interpolation_driver
    use fruit, only : init_fruit, run_test_case, fruit_summary
    use maths_mod, only : random
    use test_interpolation_mod, only : test_piecewise_zero_nnodes,test_piecewise_zgrid_eq_zvals,&
         test_interpolate_linear, test_no_interpolation
    
    implicit none
    
    call init_fruit
    call run_test_case(test_piecewise_zero_nnodes, "Test piecewise when zvals is empty")
    call run_test_case(test_piecewise_zgrid_eq_zvals, "Test piecewise when zgrid=zvals(nn-1)")
    call run_test_case(test_interpolate_linear, "Test results with linear extrapolation")
    call run_test_case(test_no_interpolation, "Test there is no interpolation if z=zvals(1)=zvals(nn)")
    call fruit_summary
  end program test_interpolation_driver
