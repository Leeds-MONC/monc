! Unit tests for the options database functionality
module test_optionsdatabase_mod
  use optionsdatabase_mod, only : LOGICAL_TYPE, INTEGER_TYPE, REAL_TYPE, STRING_TYPE, get_argument_value_type, &
       add_specific_option_key_value_pair, add_specific_option_key_value_pair
  use collections_mod, only : map_type, c_get
  use conversions_mod, only : conv_to_logical, conv_to_integer, conv_to_real, conv_to_string
  use fruit, only : assert_equals, assert_true, assert_false
  implicit none

  contains

  ! Tests the procedure that applies types to specific values
  subroutine test_get_argument_value_type()

    call assert_equals(LOGICAL_TYPE, get_argument_value_type("true"), "True is logical")
    call assert_equals(LOGICAL_TYPE, get_argument_value_type("false"), "False is logical")
    call assert_equals(INTEGER_TYPE, get_argument_value_type("12345"), "Simple number is integer")
    call assert_equals(INTEGER_TYPE, get_argument_value_type("+98"), "Positive number is integer")
    call assert_equals(INTEGER_TYPE, get_argument_value_type("-55232"), "Negative number is integer")
    call assert_equals(INTEGER_TYPE, get_argument_value_type("0"), "Zero is integer")
    call assert_equals(REAL_TYPE, get_argument_value_type("1.2"), "Simple floating point is real")
    call assert_equals(REAL_TYPE, get_argument_value_type("1.9867e6"), "Exponent floating point is real")
    call assert_equals(REAL_TYPE, get_argument_value_type("72.54e-6"), "Small exponent floating point is real")
    call assert_equals(REAL_TYPE, get_argument_value_type("-9983.2324"), "Negative floating point is real")
    call assert_equals(STRING_TYPE, get_argument_value_type("hello"), "Other characters are string")
    call assert_equals(STRING_TYPE, get_argument_value_type("1234abc"), "Number and then non numeric characters are string")
  end subroutine test_get_argument_value_type

  ! Test the addition of logical k-v pairs
  subroutine test_add_logical_key_value_pairs()
    type(map_type) :: options_database
    class(*), pointer :: raw_data

    call add_specific_option_key_value_pair(LOGICAL_TYPE, options_database, "--test1")
    call add_specific_option_key_value_pair(LOGICAL_TYPE, options_database, "--test2=false")
    call add_specific_option_key_value_pair(LOGICAL_TYPE, options_database, "--test3=true")

    raw_data=>c_get(options_database, "test1")
    call assert_true(conv_to_logical(raw_data, .false.), "Logical no value is true")
    raw_data=>c_get(options_database, "test2")
    call assert_false(conv_to_logical(raw_data, .false.), "Logical false values")
    raw_data=>c_get(options_database, "test3")
    call assert_true(conv_to_logical(raw_data, .false.), "Logical true values")
  end subroutine test_add_logical_key_value_pairs

  ! Test the addition of integer k-v pairs
  subroutine test_add_integer_key_value_pairs()
    type(map_type) :: options_database
    class(*), pointer :: raw_data

    call add_specific_option_key_value_pair(INTEGER_TYPE, options_database, "--test1=9542")
    call add_specific_option_key_value_pair(INTEGER_TYPE, options_database, "--test2=-1234")
    call add_specific_option_key_value_pair(INTEGER_TYPE, options_database, "--test3=+82")
    call add_specific_option_key_value_pair(INTEGER_TYPE, options_database, "--test4=0")

    raw_data=>c_get(options_database, "test1")
    call assert_equals(9542, conv_to_integer(raw_data, .false.), "Simple integer value added")
    raw_data=>c_get(options_database, "test2")
    call assert_equals(-1234, conv_to_integer(raw_data, .false.), "Negative integer value added")
    raw_data=>c_get(options_database, "test3")
    call assert_equals(82, conv_to_integer(raw_data, .false.), "Positive integer value added")
    raw_data=>c_get(options_database, "test4")
    call assert_equals(0, conv_to_integer(raw_data, .false.), "Zero integer value added")

  end subroutine test_add_integer_key_value_pairs

  ! Test the addition of real k-v pairs
  subroutine test_add_real_key_value_pairs()
    type(map_type) :: options_database
    class(*), pointer :: raw_data

    call add_specific_option_key_value_pair(REAL_TYPE, options_database, "--test1=9542.342")
    call add_specific_option_key_value_pair(REAL_TYPE, options_database, "--test2=-1.33e9")
    call add_specific_option_key_value_pair(REAL_TYPE, options_database, "--test3=987.232e-6")
    call add_specific_option_key_value_pair(REAL_TYPE, options_database, "--test4=-98.72")

    raw_data=>c_get(options_database, "test1")
    call assert_equals(9542.342, conv_to_real(raw_data, .false.), "Simple float value added")
    raw_data=>c_get(options_database, "test2")
    call assert_equals(-1.33e9, conv_to_real(raw_data, .false.), "Negative with exponent value added")
    raw_data=>c_get(options_database, "test3")
    call assert_equals(987.232e-6, conv_to_real(raw_data, .false.), &
      "Positive with negative exponent value added")
      raw_data=>c_get(options_database, "test4")
    call assert_equals(-98.72, conv_to_real(raw_data, .false.), "Negative float value added")

  end subroutine test_add_real_key_value_pairs

  ! Tests the addition of string k-v pairs
  subroutine test_add_string_key_value_pairs()
    type(map_type) :: options_database
    class(*), pointer :: raw_data

    call add_specific_option_key_value_pair(STRING_TYPE, options_database, "--test1=abc")
    call add_specific_option_key_value_pair(STRING_TYPE, options_database, "--test2=HeLlO")

    raw_data=>c_get(options_database, "test1")
    call assert_equals("abc", conv_to_string(raw_data, .false., 44), "Simple string")
    raw_data=>c_get(options_database, "test2")
    call assert_equals("HeLlO", conv_to_string(raw_data, .false., 44), "Case sensitive string")

  end subroutine test_add_string_key_value_pairs
end module test_optionsdatabase_mod

program test_optionsdatabase_driver
  use fruit, only : init_fruit, run_test_case, fruit_summary
  use test_optionsdatabase_mod, only : test_get_argument_value_type, test_add_logical_key_value_pairs, &
       test_add_integer_key_value_pairs, test_add_real_key_value_pairs, test_add_string_key_value_pairs

  implicit none

  call init_fruit
  call run_test_case(test_get_argument_value_type, "Test argument value type")
  call run_test_case(test_add_logical_key_value_pairs, "Test logical key value addition")
  call run_test_case(test_add_integer_key_value_pairs, "Test integer key value addition")
  call run_test_case(test_add_real_key_value_pairs, "Test real key value addition")
  call run_test_case(test_add_string_key_value_pairs, "Test string key value addition")
  call fruit_summary
end program test_optionsdatabase_driver
