! Tests the checkpointing functionality. Uses a dummy NetCDF implementation for unit testing
module test_checkpointer_mod
  use dummy_netcdf_mod, only : integer_array_wrapper_type, nf90_get_var_integer, dummy_netcdf_reset, variable_ids, &
       variable_data, global_attributes
  use checkpointer_write_checkpoint_mod, only : checkpoint_title, define_velocity_variable, write_out_misc_variables, &
       write_out_global_attributes
  use checkpointer_read_checkpoint_mod, only : does_field_exist
  use checkpointer_common_mod, only : CREATED_ATTRIBUTE_KEY, TITLE_ATTRIBUTE_KEY
  use collections_mod, only : c_get, c_contains, c_put
  use conversions_mod, only : conv_to_generic, conv_to_integer, conv_to_string
  use fruit, only : assert_equals, assert_true, assert_false 
  use state_mod, only : model_state_type
  implicit none

contains

  ! Tests the writing of global attributes
  subroutine test_write_out_global_attributes()
    class(*), pointer :: raw_data

    call dummy_netcdf_reset()
    call write_out_global_attributes(1)

    call assert_true(c_contains(global_attributes, TITLE_ATTRIBUTE_KEY), "Checkpoint title exists")
    raw_data=>c_get(global_attributes, TITLE_ATTRIBUTE_KEY)
    call assert_equals(CHECKPOINT_TITLE, trim(conv_to_string(raw_data, .true., 20)), "Checkpoint title correct")
    call assert_true(c_contains(global_attributes, CREATED_ATTRIBUTE_KEY), "Created attribute exists")
  end subroutine test_write_out_global_attributes

  ! Tests writing out misc variables
  subroutine test_write_out_misc_variables()
    type(model_state_type) :: state_mod
    integer, dimension(1) :: test_int
    integer :: ret

    call dummy_netcdf_reset()
    state_mod%timestep=100
    call write_out_misc_variables(state_mod, 1, 1, 2, 3, 4, 5)

    ret = nf90_get_var_integer(1, 1, test_int)
    call assert_equals(test_int(1), state_mod%timestep, "Timestep correct")
  end subroutine test_write_out_misc_variables

  ! Tests defining of velocity variables of different dimensions
  subroutine test_define_velocity_variable()
    integer :: id
    class(*), pointer :: raw_data, raw_id

    call dummy_netcdf_reset()
    call define_velocity_variable(1, 10, field_name="A",field_id=id)
    call assert_equals(1, id)
    raw_data => c_get(variable_data, "A")
    raw_id => c_get(variable_ids, "A")
    select type(raw_data)
      type is (integer_array_wrapper_type)
        call assert_equals(id, conv_to_integer(raw_id, .false.))
        call assert_equals(1, raw_data%size, "Size expected")
        call assert_equals(10, raw_data%data(1), "Dimension correct")
    end select
    call define_velocity_variable(1, 11, 20, field_name="B",field_id=id)
    call assert_equals(2, id)
    raw_data => c_get(variable_data, "B")
    raw_id => c_get(variable_ids, "B")
    select type(raw_data)
      type is (integer_array_wrapper_type)
        call assert_equals(id, conv_to_integer(raw_id, .false.))
        call assert_equals(2, raw_data%size, "Size expected")
        call assert_equals(11, raw_data%data(1), "Dimension correct")
        call assert_equals(20, raw_data%data(2), "Dimension correct")
    end select
    call define_velocity_variable(1, 11, 21, 30, field_name="C",field_id=id)
    call assert_equals(3, id)
    raw_data => c_get(variable_data, "C")
    raw_id => c_get(variable_ids, "C")
    select type(raw_data)
      type is (integer_array_wrapper_type)
        call assert_equals(id, conv_to_integer(raw_id, .false.))
        call assert_equals(3, raw_data%size, "Size expected")
        call assert_equals(11, raw_data%data(1), "Dimension correct")
        call assert_equals(21, raw_data%data(2), "Dimension correct")
        call assert_equals(30, raw_data%data(3), "Dimension correct")
    end select
  end subroutine test_define_velocity_variable

  ! Tests the field exists subroutine
  subroutine test_does_field_exist()
    logical :: exists
    integer :: dummy=1

    class(*), pointer :: raw_data

    raw_data=>conv_to_generic(dummy, .true.)
    call c_put(variable_ids, "BCD", raw_data)

    exists = does_field_exist(1, "ABC")
    call assert_false(exists, "Field does not exist")
    exists = does_field_exist(1, "BCD")
    call assert_true(exists, "Field does exist")
  end subroutine test_does_field_exist

end module test_checkpointer_mod

! Driver for checkpoint tests
program test_checkpointer_driver
  use fruit, only : init_fruit, run_test_case, fruit_summary
  use test_checkpointer_mod, only : test_write_out_global_attributes, test_write_out_misc_variables, &
       test_define_velocity_variable, test_does_field_exist

  implicit none

  call init_fruit
  call run_test_case(test_write_out_global_attributes, "Test writing of global attributes")
  call run_test_case(test_write_out_misc_variables, "Test writing of misc variables")
  call run_test_case(test_define_velocity_variable, "Test writing of velocity variables")
  call run_test_case(test_does_field_exist, "Test whether field exists or not")
  call fruit_summary
end program test_checkpointer_driver
