! Tests the logging_mod utility functions
module test_saturation_mod
  use fruit, only : assert_equals, assert_not_equals, assert_true
  use datadefn_mod, only : DEFAULT_PRECISION
  use saturation_mod, only :  qsaturation, dqwsatdt
  implicit none
  
contains
  
  ! Test qsaturation when P=T=0.0
  subroutine test_qsaturation_zero
    real(kind=DEFAULT_PRECISION) :: res
    real(kind=DEFAULT_PRECISION),parameter :: expected = -0.62203306596824359_DEFAULT_PRECISION
    res = qsaturation(0.0_DEFAULT_PRECISION, 0.0_DEFAULT_PRECISION)
    call assert_equals(expected, res, "Result when P=0 and T=0")
   end subroutine test_qsaturation_zero

   ! Test dqwsatdt when SATURATION=T=0.0
  subroutine test_dqwsatdt_zero
    real(kind=DEFAULT_PRECISION) :: res
    real(kind=DEFAULT_PRECISION),parameter :: expected = 0.0_DEFAULT_PRECISION
    res = dqwsatdt(0.0_DEFAULT_PRECISION, 0.0_DEFAULT_PRECISION)
    call assert_equals(expected, res, "Result when SAT=0 and T=0 is 0")
   end subroutine test_dqwsatdt_zero

   ! Test dqwsatdt result is 0 when T=qsa3
   subroutine test_dqwsatdt_temp_qsa3
    real(kind=DEFAULT_PRECISION) :: res
    real(kind=DEFAULT_PRECISION),parameter :: expected = 0.0_DEFAULT_PRECISION
    res = dqwsatdt(0.0_DEFAULT_PRECISION,35.86_DEFAULT_PRECISION)
    call assert_equals(expected, res, "Result when T=qsa3 should be 0")
    write(*,*)res
   end subroutine test_dqwsatdt_temp_qsa3
   

 end module test_saturation_mod

  

  ! Driver for maths_mod utility tests
  program test_saturation_driver
    use fruit, only : init_fruit, run_test_case, fruit_summary
    use test_saturation_mod, only : test_qsaturation_zero,test_dqwsatdt_zero,test_dqwsatdt_temp_qsa3
    
    implicit none
    
    call init_fruit
    call run_test_case(test_qsaturation_zero, "Test the result of qsaturation with T=P=0")
    call run_test_case(test_dqwsatdt_zero, "Test the result of  dqwsatdt with SAT=P=0")
    call run_test_case(test_dqwsatdt_temp_qsa3, "Test the result of  dqwsatdt when T=qsa3")
    call fruit_summary
  end program test_saturation_driver
