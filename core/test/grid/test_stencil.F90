! Tests the logging_mod utility functions
module test_stencil_mod
  use fruit, only : assert_equals, assert_not_equals, assert_true
  use datadefn_mod, only : DEFAULT_PRECISION
  use stencil_mod, only : create_stencil, interpolate_to_dual, free_stencil,grid_stencil_type,&
       calculate_interpolated_cell_value
  use grids_mod, only : local_grid_type, X_INDEX, Y_INDEX, Z_INDEX
  use prognostics_mod, only : prognostic_field_type, prognostic_field_ptr_type, &
       get_field_interpolation_index
  use state_mod, only : model_state_type, parallel_state_type, FORWARD_STEPPING
  implicit none
  
contains
  
  ! Test random produces a non zero value
  subroutine test_create_stencil
    type(model_state_type), target :: current_state
    type(grid_stencil_type) :: star_stencil
    type(prognostic_field_ptr_type), dimension(3)  :: fields
    integer, dimension(3, 2) :: sizes
    type(prognostic_field_type), dimension(:), allocatable :: interpolated_fields
    integer :: i,num_fields,max_y_point
    
    call assert_equals(0,star_stencil%nfields,"Test a stencil has been created")
    current_state%global_grid%size(Z_INDEX) = 2
    num_fields = 15
    current_state%local_grid%size(Y_INDEX) = 3
    current_state%local_grid%halo_size(Y_INDEX) = 4

    fields(num_fields)%ptr => current_state%u
    sizes(num_fields,:) = (/ 2, 2 /)
    
    allocate(interpolated_fields(0:num_fields))    
    allocate(interpolated_fields(num_fields)%data(current_state%global_grid%size(Z_INDEX),&
         -1:3, -1:3))
    max_y_point = (current_state%local_grid%size(Y_INDEX)+ &
         current_state%local_grid%halo_size(Y_INDEX) *2)-1
    
    interpolated_fields(num_fields)%active=.true.
    star_stencil = create_stencil(current_state%local_grid, fields, num_fields, 3, &
         sizes,.true., .false.)
    call assert_equals(max_y_point,star_stencil%max_y_point,"Test max Y point within the stencil")
    call assert_equals(-1,star_stencil%max_x_point,"Test max X point within the stencil")
    call assert_equals(num_fields+1,size(interpolated_fields),"Test size of interpolated_fields")
    call assert_equals(num_fields,star_stencil%nfields,"Test a stencil has been created")
  end subroutine test_create_stencil

  ! Test calculate interpolated cell value
  subroutine test_calculate_interpolated_cell_value
    type(model_state_type), target :: current_state
    type(grid_stencil_type) :: star_stencil
    type(prognostic_field_ptr_type), dimension(3)  :: fields
    integer, dimension(3, 2) :: sizes
    type(prognostic_field_type), dimension(:), allocatable :: interpolated_fields
    integer :: i,num_fields,max_y_point
    logical, dimension(3) :: interpolate_in_dimension
    real(kind=DEFAULT_PRECISION) :: r_value
    real(kind=DEFAULT_PRECISION) :: e_value

    
    current_state%global_grid%size(Z_INDEX) = 2
    num_fields = 15
    current_state%local_grid%size(Z_INDEX)= 0
    current_state%local_grid%size(Y_INDEX) = 3
    current_state%local_grid%halo_size(Y_INDEX) = 4

    fields(num_fields)%ptr => current_state%u
    sizes(num_fields,:) = (/ 2, 2 /)
    
    allocate(interpolated_fields(0:num_fields))    
    allocate(interpolated_fields(num_fields)%data(current_state%global_grid%size(Z_INDEX),&
         -1:3, -1:3))
    max_y_point = (current_state%local_grid%size(Y_INDEX)+ &
         current_state%local_grid%halo_size(Y_INDEX) *2)-1
    
    interpolated_fields(num_fields)%active=.true.
    star_stencil = create_stencil(current_state%local_grid, fields, num_fields, 3, &
         sizes,.true., .false.)
    e_value = 10_DEFAULT_PRECISION
    star_stencil%interpolated_fields(1,1)%data(1,1,1) = e_value
   
    call calculate_interpolated_cell_value(current_state%local_grid, 1, 1, (/.true.,.true.,.true./),&
         star_stencil%fields(1), star_stencil%interpolated_fields(1,1)%data(:,1,1))  
    
    r_value =  star_stencil%interpolated_fields(1,1)%data(1,1,1)
    
    call assert_equals(e_value, r_value,"Test")
  end subroutine test_calculate_interpolated_cell_value

  end module test_stencil_mod

  

  ! Driver for maths_mod utility tests
  program test_stencil_driver
    use fruit, only : init_fruit, run_test_case, fruit_summary
    use maths_mod, only : random
    use test_stencil_mod, only : test_create_stencil,test_calculate_interpolated_cell_value
    
    implicit none
    
    call init_fruit
    call run_test_case(test_create_stencil, "Test maths_mod random function produces output")
    call run_test_case(test_calculate_interpolated_cell_value, &
         "Test the calculation of interpolated cell values")

    call fruit_summary
  end program test_stencil_driver
