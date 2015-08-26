! Unit tests for the options database functionality
module test_halo_communication_mod
  use halo_communication_mod, only : copy_buffer_to_field, copy_field_to_buffer, copy_corner_to_buffer, &
       copy_buffer_to_corner, perform_local_data_copy_for_field, init_halo_communication, &
       finalise_halo_communication, initiate_nonblocking_halo_swap, complete_nonblocking_halo_swap,&
       blocking_halo_swap, get_single_field_per_halo_cell, get_number_communication_requests, &
       determine_halo_corner_size, determine_halo_corner_element_sizes, get_pid_neighbour_location,&
       has_pid_already_been_seen, retrieve_same_neighbour_information,perform_local_data_copy_for_dimension
  use collections_mod, only : map_type, c_get
  use conversions_mod, only : conv_to_logical, conv_to_integer, conv_to_real, conv_to_string
  use communication_types_mod, only:  halo_communication_type, neighbour_description_type, &
       field_data_wrapper_type
  use grids_mod, only : local_grid_type, X_INDEX, Y_INDEX, Z_INDEX
  use fruit, only : assert_equals, assert_true, assert_false, assert_not_equals
  use datadefn_mod, only : DEFAULT_PRECISION
  implicit none

  contains


  ! Test the number of requests
  subroutine test_get_number_communication_requests
    type(neighbour_description_type), dimension(:), allocatable :: halo_swap_neigh
    integer :: i, requests
    allocate(halo_swap_neigh(6))    
    
    do i=1,5
       halo_swap_neigh(i)%recv_size = i
       halo_swap_neigh(i)%recv_corner_size = i
    enddo
    
    requests = get_number_communication_requests(halo_swap_neigh, 5)
    call assert_equals(10, requests, "Test number of requests")
  end subroutine test_get_number_communication_requests

  ! Test the size of halo corner
  subroutine test_determine_halo_corner
    type(local_grid_type) :: local_grid
    integer :: halo_corner_size

    local_grid%halo_size(X_INDEX)=2
    local_grid%halo_size(Y_INDEX)=2
    local_grid%size(Z_INDEX) = 10
    
    halo_corner_size = determine_halo_corner_size(local_grid)

    call assert_equals(10*2*2, halo_corner_size, "Test halo_size")
  end subroutine test_determine_halo_corner

  ! Test the size of halo corner elements
  subroutine test_determine_halo_corner_elements
    type(local_grid_type) :: local_grid
    integer :: elemNumber
    allocate( local_grid%corner_neighbours(4,2) )
    local_grid%corner_neighbours(1,1) = 5
    local_grid%corner_neighbours(2,2) = 5 
    local_grid%size(Z_INDEX) = 10
    elemNumber = determine_halo_corner_element_sizes(local_grid, 5)    
    
    !4 elements per corner*10 in Z direction
    call assert_equals(10*4, elemNumber, "Test number elements halo_size")
    ! if not PID elemNumber should be 0
    elemNumber = determine_halo_corner_element_sizes(local_grid, 15)    
    call assert_equals(0, elemNumber, "Test number elements halo_size")
  end subroutine test_determine_halo_corner_elements

  ! Test get neighbour pid
  subroutine test_get_pid_neighbour_location
    type(local_grid_type) :: local_grid
    type(neighbour_description_type), dimension(:), allocatable :: halo_swap_neigh
    integer :: i, pids

    allocate(halo_swap_neigh(6))    
    do i=1,5
       halo_swap_neigh(i)%pid = i
    enddo
    halo_swap_neigh(2)%pid=1
    pids = get_pid_neighbour_location(halo_swap_neigh, 1, 5)
    
    call assert_equals(1, pids, "Test unique pid location")
    pids = get_pid_neighbour_location(halo_swap_neigh, 10, 5)
    call assert_equals(-1, pids, "Test pid location")
  end subroutine test_get_pid_neighbour_location
  
  ! Test has_pid_already_been_seen which actually tells if a list contains a pid
  subroutine test_pid_been_seen
    integer :: neigh_pids(8)
    logical :: seen
    neigh_pids(1) = -1
    neigh_pids(2) = 2
    neigh_pids(3) = 0

    seen  = has_pid_already_been_seen(neigh_pids, 1)
    call  assert_false( seen, "Test not seen if pid=-1" )
    neigh_pids(1) = 0
    seen  = has_pid_already_been_seen(neigh_pids, 2)
    call  assert_true( seen, "Test seen when there is not -1 if pid=2" )
    seen  = has_pid_already_been_seen(neigh_pids, 0)
    call  assert_true( seen, "Test seen when there is not - 1 if pid=0" )
    
  end subroutine test_pid_been_seen
  
  ! Test if we have same neighbourd L,R,U,D
  subroutine test_retrieve_same_neighbour_information
    type(local_grid_type) :: local_grid
    logical, dimension(3) :: retrieve_same_neigh_info
    integer :: i
    allocate( local_grid%neighbours(3,4) )
    local_grid%halo_size(Y_INDEX) = 2
    local_grid%halo_size(X_INDEX) = 2

    local_grid%neighbours(Y_INDEX,1) = 2
    local_grid%neighbours(Y_INDEX,2) = 2
    local_grid%neighbours(Y_INDEX,3) = 8
    local_grid%neighbours(Y_INDEX,4) = 8

    local_grid%neighbours(X_INDEX,1) = 4
    local_grid%neighbours(X_INDEX,2) = 4
    local_grid%neighbours(X_INDEX,3) = 6
    local_grid%neighbours(X_INDEX,4) = 6

    retrieve_same_neigh_info = retrieve_same_neighbour_information(local_grid)

    call assert_true(retrieve_same_neigh_info(1),"Test Z")
    call assert_false(retrieve_same_neigh_info(2),"Test Y")
    call assert_false(retrieve_same_neigh_info(3),"Test Z")


    local_grid%neighbours(Y_INDEX,3) = 2
    local_grid%neighbours(Y_INDEX,4) = 2

    retrieve_same_neigh_info = retrieve_same_neighbour_information(local_grid)

    call assert_true(retrieve_same_neigh_info(1),"Test Z")
    call assert_true(retrieve_same_neigh_info(2),"Test Y")
    call assert_false(retrieve_same_neigh_info(3),"Test Z")
    
  end subroutine test_retrieve_same_neighbour_information

  ! Test a local data copy for a specific dimension
  subroutine test_perform_local_data_copy_for_dimension
    type(local_grid_type) :: local_grid
    real(kind=DEFAULT_PRECISION), dimension(10,20,30) :: field_data
    integer :: i,j,k
    
     allocate( local_grid%neighbours(3,4) )
     ! start
     local_grid%local_domain_start_index(Z_INDEX) = 1
     local_grid%local_domain_start_index(Y_INDEX) = 11
     local_grid%local_domain_start_index(X_INDEX) = 21
     ! end
     local_grid%local_domain_end_index(Z_INDEX) = 10
     local_grid%local_domain_end_index(Y_INDEX) = 20
     local_grid%local_domain_end_index(X_INDEX) = 30
     
     do j=1,3
        do i =1,4
           local_grid%neighbours(j,i) = i*j
        enddo
     enddo
     do i=1,10
        do j=1,20
           do k=1,30
              field_data(i,j,k) = i*j*k
           enddo
        enddo
     enddo
  
     call assert_not_equals(field_data(1,2,21), field_data(1,10,21)&
          , "Test fields are not equal before calling")
     call perform_local_data_copy_for_dimension(X_INDEX, 1, 1, local_grid, &
          field_data)
  
   end subroutine test_perform_local_data_copy_for_dimension
   
 end module test_halo_communication_mod

program test_halo_communication_driver
  use fruit, only : init_fruit, run_test_case, fruit_summary
  use test_halo_communication_mod, only : test_get_number_communication_requests, &
       test_determine_halo_corner, test_determine_halo_corner_elements,&
       test_get_pid_neighbour_location, test_pid_been_seen,&
       test_retrieve_same_neighbour_information,test_perform_local_data_copy_for_dimension
  implicit none

  call init_fruit
  call run_test_case(test_get_number_communication_requests, "Test requests")
  call run_test_case(test_determine_halo_corner, "Test size of halo corners")
  call run_test_case(test_determine_halo_corner_elements, "Test number elements in halo corners")
  call run_test_case(test_get_pid_neighbour_location, "Test pid location")
  call run_test_case(test_pid_been_seen, "Test pid been seen")
  call run_test_case(test_retrieve_same_neighbour_information, "Test retrieve neighbour info")
  call run_test_case(test_perform_local_data_copy_for_dimension, "Test local data copies")
  call fruit_summary
end program test_halo_communication_driver
