! Tests the logging_mod utility functions
module test_maths_mod
  use fruit, only : assert_equals, assert_not_equals, &
       assert_true
  use datadefn_mod, only : DEFAULT_PRECISION
  use maths_mod, only : random
  implicit none
  
contains
  
  ! Test random produces a non zero value
  subroutine test_random
    
    integer :: idum 
    real :: res

    res = 0.0 
    idum = 0
    res = random(idum)
    
    call assert_not_equals(res, 0.0, "Result is not 0")
    
  end subroutine test_random
  
  ! Test random produces different values 
  subroutine test_randomness
    
    integer :: idum 
    real :: res1, res2
    
    res1 = 0.0 
    res2 = 0.0
    idum = -4
    res1 = random(idum)
    res2 = random(idum)
    call assert_not_equals(res1, res2, &
         "Result are not equal")
      
  end subroutine test_randomness
  
end module test_maths_mod

  

  ! Driver for maths_mod utility tests
  program test_maths_mod_driver
    use fruit, only : init_fruit, run_test_case, fruit_summary
    use maths_mod, only : random
    use test_maths_mod, only : test_random, test_randomness
    
    implicit none
    
    call init_fruit
    call run_test_case(test_random, "Test maths_mod random function produces output")
    call run_test_case(test_randomness, "Test maths_mod randomness")
    call fruit_summary
  end program test_maths_mod_driver
