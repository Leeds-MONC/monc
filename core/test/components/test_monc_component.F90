! Tests the component registry functionality
module test_monc_component_mod
  use fruit, only : assert_equals, add_fail, assert_true, &
       assert_false, set_unit_name,  assert_not_equals
  use conversions_mod, only : conv_to_generic, conv_to_string, conv_to_integer
  use state_mod, only : model_state_type
  use monc_component_mod, only :  INIT_PRIORITY_INDEX,TIMESTEP_PRIORITY_INDEX,&
       FINALISATION_PRIORITY_INDEX, component_descriptor_type
  use collections_mod, only : map_type, c_size, c_value_at, c_key_at
  implicit none



contains
  ! Test 
  subroutine test_parameters
    call assert_not_equals(2,INIT_PRIORITY_INDEX,"Test INIT_PRIORITY_INDEX is not real")
    call assert_equals(1,INIT_PRIORITY_INDEX,"Test INIT_PRIORITY_INDEX is value expected")
    call assert_equals(2,TIMESTEP_PRIORITY_INDEX,"Test TIMESTEP_PRIORITY_INDEX is value expected")
    call assert_equals(5,FINALISATION_PRIORITY_INDEX,"Test FINALISATION_PRIORITY_INDEX is value expected")
  end subroutine test_parameters
  
end module test_monc_component_mod


! The driver for testing the registry
program test_monc_component_driver
  use fruit, only : init_fruit, run_test_case, fruit_summary
  use test_monc_component_mod, only :  test_parameters
  
  implicit none

  call init_fruit
  call run_test_case(test_parameters, "Testing parameters")
  
  call fruit_summary
end program test_monc_component_driver

