! Unit tests for the options database functionality
module test_communication_types_mod
  use communication_types_mod, only : field_data_wrapper_type, neighbour_description_type,&
       halo_communication_type
  use fruit, only : assert_equals, assert_true, assert_false
  implicit none

  contains

   ! Test intializated values
  subroutine test_halo_default_size
    type(neighbour_description_type) :: neighbour
    
    call assert_equals(0, neighbour%halo_pages, "Test halo_pages")
    call assert_equals(0, neighbour%halo_corners, "Test halo_pages")
  end subroutine test_halo_default_size

  subroutine test_halo_communication_default_values
    type(halo_communication_type) :: halocomm
    
    call assert_false(halocomm%initialised, "Test initialised")
    call assert_false(halocomm%involve_corners, "Test involve_corners")
    call assert_false(halocomm%swap_in_progress, "Test swap_in_progress")
  end subroutine test_halo_communication_default_values
  
end module test_communication_types_mod

program test_communication_types_driver
  use fruit, only : init_fruit, run_test_case, fruit_summary
  use test_communication_types_mod, only : test_halo_default_size, &
       test_halo_communication_default_values

  implicit none

  call init_fruit
  call run_test_case(test_halo_default_size, "Test initial halo values")
  call run_test_case(test_halo_communication_default_values, "Test ")
  call fruit_summary
end program test_communication_types_driver
