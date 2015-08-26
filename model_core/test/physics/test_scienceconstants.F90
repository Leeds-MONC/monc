! Tests the logging_mod utility functions
module test_science_constants_mod
  use fruit, only : assert_equals, assert_not_equals, assert_true
  use datadefn_mod, only : DEFAULT_PRECISION
  use monc_mod, only: load_model_configuration
  use science_constants_mod, only :  smallp, von_karman_constant, z0, z0th, alphah, &
       betam, betah, gammam, gammah, pi, surface_vapour_mixing_ratio, cp, rlvap, &
       rlvap_over_cp, r, r_over_cp, G, convective_limit, ratio_mol_wts, &
       rlargep, initialise_science_constants, seconds_in_a_day
  
  use state_mod, only : model_state_type
  implicit none
  
contains
  
  ! Test constants
  subroutine test_science_constants
  
    type(model_state_type) :: state
    call load_model_configuration(state%options_database)
    call initialise_science_constants(state)
  

    call assert_equals(rlvap_over_cp,rlvap/cp, "Test rlvap_over is what is expected")
    call assert_equals(r_over_cp,r/cp, "Test rlvap_over is what is expected")
   end subroutine test_science_constants


 end module test_science_constants_mod

  

  ! Driver for maths_mod utility tests
  program test_science_constants_driver
    use fruit, only : init_fruit, run_test_case, fruit_summary
    use test_science_constants_mod, only : test_science_constants
    
    implicit none
    
    call init_fruit
    call run_test_case(test_science_constants, "Test the initialization of science constants")
    call fruit_summary
  end program test_science_constants_driver
