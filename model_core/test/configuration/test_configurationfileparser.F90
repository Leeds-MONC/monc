! Unit tests for the options database functionality
module test_configuration_file_parser_mod
  use collections_mod, only : map_type, c_get
  use conversions_mod, only : conv_to_logical, conv_to_integer, conv_to_real, &
       conv_is_logical, conv_is_integer, conv_is_real, conv_to_string
  use optionsdatabase_mod, only : options_add, options_get_string, options_has_key, &
       options_get_array_size, options_remove_key
  use configuration_file_parser_mod, only :  parse_configuration_file
  use fruit, only : assert_equals, assert_true, assert_false
  
  implicit none
  
contains
  
  ! Test some modules are read correctly
  subroutine test_parse_configuration_file
    type(map_type) :: options_database
      class(*), pointer :: raw_data

      call parse_configuration_file(options_database, "user_config")
      
      raw_data=>c_get(options_database, "viscosity_enabled")
      call assert_false(conv_to_logical(raw_data, .false.), "Test viscosity is not enable")
      
      raw_data=>c_get(options_database, "rhobous")
      call assert_equals(1.0, conv_to_real(raw_data, .false.), "Test setup configuration values")
      
      raw_data=>c_get(options_database, "passive_q")
      call assert_false(conv_to_logical(raw_data, .false.), "Test logical value is false when missing")
      
      raw_data=>c_get(options_database, "number_q_fields")
      call assert_equals(5, conv_to_integer(raw_data, .false.), "Test integer values read")

      raw_data=>c_get(options_database, "checkpoint_enable_write")
      call assert_true(conv_to_logical(raw_data, .false.), "Test result when entry is not defined in the user file but in  the global file")

    end subroutine test_parse_configuration_file

end module test_configuration_file_parser_mod

program test_configuration_file_parser_driver
  use fruit, only : init_fruit, run_test_case, fruit_summary
  use test_configuration_file_parser_mod, only : test_parse_configuration_file
  implicit none

  call init_fruit
  call run_test_case(test_parse_configuration_file, "Test values from user file")
  call fruit_summary
end program test_configuration_file_parser_driver
