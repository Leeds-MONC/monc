! Tests the logging_mod utility functions
module test_monc_mod
  use fruit, only : assert_equals, assert_not_equals, assert_true
  use state_mod, only : model_state_type
  use datadefn_mod, only : DEFAULT_PRECISION
  use monc_mod, only : load_model_configuration, fill_registry_with_components,&
       display_registed_components
  use optionsdatabase_mod, only : options_get_string, options_get_real
  use fftsolver_mod, only : fftsolver_get_descriptor
  use monc_component_mod, only : component_descriptor_type
  use collections_mod, only : list_type, c_add, map_type, c_size, c_key_at,&
       c_value_at, c_get
  use registry_mod, only : init_registry,  get_all_registered_components
  implicit none
  
contains
  
  subroutine add_component(component_descriptions, single_description)
    type(list_type), intent(inout) :: component_descriptions
    type(component_descriptor_type), intent(in) :: single_description

    class(*), pointer :: raw_data
    allocate(raw_data, source=single_description)
    call c_add(component_descriptions, raw_data)
  end subroutine add_component
  

  ! Test random produces a non zero value
  subroutine test_load_model
    type(model_state_type) :: current_state
    real :: z0
    
    call load_model_configuration(current_state%options_database)
    z0 = options_get_real(current_state%options_database, "z0")
    call assert_equals(2.0e-3,z0,"Test z0 has been read properly")
  end subroutine test_load_model
  
 ! Test random produces a non zero value
  subroutine test_fill_registry_components
    type(model_state_type) :: current_state
    type(list_type) :: component_descriptions
    type(map_type) :: registered_components
   
    registered_components = get_all_registered_components()
    call assert_equals(0,c_size(registered_components),"Test there are not components in the registry")
    ! add the component to the database
    call add_component(component_descriptions, fftsolver_get_descriptor())
    call fill_registry_with_components(current_state%options_database, component_descriptions)
    
    registered_components = get_all_registered_components()
    call assert_equals(1,c_size(registered_components),"Test there is only 1 component in the registry")
  end subroutine test_fill_registry_components
  
end module test_monc_mod

  

  ! Driver for maths_mod utility tests
program test_monc_driver
  use fruit, only : init_fruit, run_test_case, fruit_summary
  use test_monc_mod, only : test_load_model, test_fill_registry_components
  
  implicit none
    
  call init_fruit
  call run_test_case(test_load_model, "Test loading the configuration model to options_database")
  call run_test_case(test_fill_registry_components, "Test filling registry with components")
  call fruit_summary
end program test_monc_driver

