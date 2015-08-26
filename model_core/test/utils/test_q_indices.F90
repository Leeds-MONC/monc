! Tests the collections_mod utility functions
module test_q_indices_mod
  use fruit, only : assert_equals, add_fail, assert_false, assert_true
  use q_indices_mod, only :  q_metadata_type, q_indices_add, get_indices_descriptor, get_max_number_q_indices, &
       get_number_active_q_indices, set_q_index
  implicit none

contains

  ! Tests q_indices_add - ensures that a new entry is added only if that entry does not 
  ! exist already
  subroutine test_q_indices_add
    character(7), parameter :: name='First'
    integer :: numActive,index
    index =  q_indices_add(name)
    numActive =  get_number_active_q_indices()

    call assert_equals(1, numActive, "First element")
    ! add the same element for second time
    index =  q_indices_add(name)
    numActive =  get_number_active_q_indices()
    call assert_equals(1, numActive, "Repeated element")    
  end subroutine test_q_indices_add

  ! Tests set_q_index - ensures that a new index is set a given index and name
  subroutine test_q_indices_set
    type(q_metadata_type) :: descriptor
    character(7), parameter :: name='Set'
    integer :: numActive,new_numActive,index
    index =  q_indices_add(name)
    numActive =  get_number_active_q_indices()

    ! add the same element for second time
    call set_q_index(index,"New")
    new_numActive =  get_number_active_q_indices()
    call assert_equals(new_numActive, numActive, "Overwritten element")   
 
    ! we expect to have overwritten the previous descriptor
    descriptor = get_indices_descriptor(index)
    call assert_equals("New",descriptor%name,"Check name of set element")
  end subroutine test_q_indices_set

  ! Tests limits - if an element is introduced after the n_maxqs limit, 
  ! it is not taken into account as an active index
  subroutine test_q_indices_limits
    integer :: numActive,new_numActive,index
    !number of elements active before
    numActive =  get_number_active_q_indices()
    !set new value out of range
    call set_q_index(101,"Limit")
    !number of elements active after
    new_numActive =  get_number_active_q_indices()
    !check the number of active elements has not changed
    call assert_equals(new_numActive, numActive, "Out of range element")

  end subroutine test_q_indices_limits




end module test_q_indices_mod




! Driver for the q_indices_mod unit tests
program test_q_indices_driver
  use fruit, only : init_fruit, run_test_case, fruit_summary
  use test_q_indices_mod, only : test_q_indices_add, test_q_indices_set, test_q_indices_limits

  implicit none

  call init_fruit
  call run_test_case(test_q_indices_add, "Elements addition")
  call run_test_case(test_q_indices_set, "Elements set")
  call run_test_case(test_q_indices_limits, "Elements limits")
  call fruit_summary
end program test_q_indices_driver
