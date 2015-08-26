! Tests the collections_mod utility functions
module test_naming_conventions_mod
  use fruit, only : assert_equals, add_fail, assert_false, assert_true
  use naming_conventions_mod, only :  k_per_day, k_per_second, m_per_second_per_day, m_per_second_per_second, &
     kg_per_kg_per_second, kg_per_kg_per_day, g_per_kg_per_day, g_per_kg_per_second
  implicit none

contains
  ! Tests the value of the different variables
  subroutine test_conventions_value


    call assert_equals(k_per_day,'K/day', "k_per_day")
    call assert_equals(k_per_second,'K/s', "k_per_second")
    call assert_equals(m_per_second_per_day,'m/s/day', "m_per_second_per_day")
    call assert_equals(m_per_second_per_second,'m/s/s', "m_per_second_per_second")
    call assert_equals(kg_per_kg_per_second,'kg/kg/s', "kg_per_kg_per_second")
    call assert_equals(kg_per_kg_per_day,'kg/kg/day', "kg_per_kg_per_day")
    call assert_equals(g_per_kg_per_day,'g/kg/day', "g_per_kg_per_day")
    call assert_equals(g_per_kg_per_second,'g/kg/s', "g_per_kg_per_second")

  end subroutine test_conventions_value
end module test_naming_conventions_mod

program test_naming_conventions
  use fruit, only : init_fruit, run_test_case, fruit_summary
  use test_naming_conventions_mod, only : test_conventions_value

  implicit none

  call init_fruit
  call run_test_case(test_conventions_value, "Values")
 
  call fruit_summary
end program test_naming_conventions
