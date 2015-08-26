! Tests the collections_mod utility functions
module test_collections_mod
  use fruit, only : assert_equals, add_fail, assert_false, assert_true
  use collections_mod, only : map_type, stack_type, queue_type, list_type, c_put, c_free, c_get, c_key_at, c_contains, c_size, &
       c_is_empty, c_pop, c_value_at, c_add, c_remove, c_insert, c_push, c_remove
  implicit none

contains
  ! Tests the map_type put - ensures that key-value pairs are being added correctly and are updated
  ! when additional put requests with existing keys are issued
  subroutine test_map_type_put
    type(map_type) :: my_map
    integer :: i, j, v
    class(*), pointer :: iptr

    call assert_equals(.true., c_is_empty(my_map), "map_type is empty")

    do i=1,10
      allocate(iptr, source=i)
      call c_put(my_map, str(i), iptr)
      call assert_equals(i, c_size(my_map), "map_type size incremented with put")
    end do

    do j=1,10
      v = j * 10
      allocate(iptr, source=v)
      call c_put(my_map, str(j), iptr)
      call assert_equals(10, c_size(my_map), "map_type size remains unchanged after duplicate key")
    end do

    call c_free(my_map)
  end subroutine test_map_type_put

  ! Tests that if we put many identical keys into a map_type then it will just update the value rather than
  ! add new key-value pairs
  subroutine test_map_type_put_unique
    type(map_type) :: my_map
    integer :: i
    class(*), pointer :: iptr, data

    call assert_equals(.true., c_is_empty(my_map), "map_type is empty")

    do i=1,10
      allocate(iptr, source=i)
      call c_put(my_map, "A", iptr)
      data=>c_get(my_map, "A")
      select type(data)
        type is (integer)
          call assert_equals(i, data, "Value of key-value pair is correct")
        class default
          call add_fail("Unknown type")
      end select
      call assert_equals(1, c_size(my_map), "Size of may is one due to unique key")
    end do

    call c_free(my_map)
  end subroutine test_map_type_put_unique

    subroutine test_map_type_pointers
    type(map_type) :: my_map
    integer, target :: i
    class(*), pointer :: iptr, data
    i=20
    iptr=>i
    call c_put(my_map, "A", iptr)
    i=50
    data=>c_get(my_map, "A")
    select type(data)
      type is (integer)
      call assert_equals(50, data, "Value of entry has been changed through modifying original variable")
      class default
      call add_fail("Unknown type")
    end select
    call c_free(my_map)
  end subroutine test_map_type_pointers

  ! Tests the map_type get functionality, will put in a number of key-value pairs, check that the map_type contains
  ! the key added, gets the key and checks the value along with ensuring the key is held where we expect it.
  ! Then it will modify the values associated with keys and check that these are represented correctly
  subroutine test_map_type_get
    type(map_type) :: my_map
    integer :: i, j, v, x
    class(*), pointer :: iptr, data

    call assert_equals(.true., c_is_empty(my_map), "map_type is empty")

    do i=1,10
      v = i * 10
      call assert_false(c_contains(my_map, str(i)), "map_type does not contain the key before put")
      allocate(iptr, source=v)
      call c_put(my_map, str(i), iptr)
      call assert_true(c_contains(my_map, str(i)), "map_type contains the key after put")
      data => c_get(my_map, str(i))
      select type (data)
        type is (integer)
          call assert_equals(i*10, data, "Value of entry is consistent")
      end select
      call assert_equals(str(i), c_key_at(my_map, i), "Key at location i is consistent")
    end do

    do j=1,10
      x = j * 100
      allocate(iptr, source=x)
      call c_put(my_map, str(j), iptr)
      data => c_get(my_map, str(j))
      select type (data)
        type is (integer)
          call assert_equals(j*100, data, "Value modified due to duplicate key")
        class default
          call add_fail("Unknown type")
      end select
      data => c_value_at(my_map, j)
      select type (data)
        type is (integer)
          call assert_equals(j*100, data, "Value at returned correct value at i")
        class default
          call add_fail("Unknown type")
      end select
    end do

    call c_free(my_map)
  end subroutine test_map_type_get

  ! Will add in a load of key-value pairs into the map_type and then test removing them and ensure that the map_type
  ! does not contain the key-value pair once the removal has completed.
  subroutine test_map_type_remove
    type(map_type) :: my_map
    integer :: i, j, v
    class(*), pointer :: iptr, data

    call assert_equals(.true., c_is_empty(my_map), "map_type is empty")

    do i=1,10
      v=i*10
      allocate(iptr, source=v)
      call c_put(my_map, str(i), iptr)
      data=>c_get(my_map, str(i))
      select type(data)
        type is (integer)
          call assert_equals(i*10, data, "Value of key is consistent")
        class default
          call add_fail("Unknown type")
      end select
    end do

    do j=1,10
      call assert_true(c_contains(my_map, str(11-j)), "map_type contains the key pre-removal")
      call c_remove(my_map, str(11-j))
      call assert_false(c_contains(my_map, str(11-j)), "map_type does not contain the key post-removal")
      call assert_equals(10-j, c_size(my_map), "Size of map_type is consistent post-removal")
    end do
    call assert_equals(.true., c_is_empty(my_map), "map_type is empty at the end")

    call c_free(my_map)
  end subroutine test_map_type_remove

  ! Helper function to convert an integer into a string (for map_type keying in a loop)
  character(len=15) function str(k)
    integer, intent(in) :: k
    write (str, *) k
    str = adjustl(str)
  end function str

  subroutine test_stack_type_pointers
    type(stack_type) :: my_stack
    integer, target :: i
    class(*), pointer :: iptr, data
    i=20
    iptr=>i
    call c_push(my_stack, iptr)
    i=50
    data=>c_pop(my_stack)
    select type(data)
      type is (integer)
      call assert_equals(50, data, "Stack data modified through changing original variable")
      class default
      call add_fail("Unknown type")
    end select

    call c_free(my_stack)
  end subroutine test_stack_type_pointers

  ! Tests a stack_type push and pop - ensures that it is working in LIFO order
  subroutine test_stack_type_push_pop
    type(stack_type) :: my_stack
    integer :: i,j
    class(*), pointer :: iptr, data

    call assert_equals(.true., c_is_empty(my_stack), "Stack is empty")

    do i=1,10
      allocate(iptr, source=i)
      call c_push(my_stack, iptr)
      call assert_equals(i, c_size(my_stack), "Size of stack_type increasing as data pushed")
    end do

    call assert_equals(.false., c_is_empty(my_stack), "Stack is not empty after values pushed")

    do j=1,10
      data => c_pop(my_stack)
      select type(data)
        type is (integer)
          call assert_equals(data, 11-j, "Stack pop gives LIFO value")
          call assert_equals(10-j, c_size(my_stack), "Stack pop removes the LIFO value")
        class default
          call add_fail("Type unknown")
      end select
    end do

    call c_free(my_stack)
  end subroutine test_stack_type_push_pop

  subroutine test_queue_type_pointers
    type(queue_type) :: my_queue
    integer, target :: i
    class(*), pointer :: iptr, data
    i=20
    iptr=>i
    call c_push(my_queue, iptr)
    i=50
    data=>c_pop(my_queue)
    select type(data)
      type is (integer)
      call assert_equals(50, data, "Queue data modified by changing the original variable")
      class default
      call add_fail("Unknown type")
    end select

    call c_free(my_queue)
  end subroutine test_queue_type_pointers

  ! Tests a queue_type push and pop - ensures that it is working in FIFO order
  subroutine test_queue_type_push_pop
    type(queue_type) :: my_queue
    integer :: i,j
    class(*), pointer :: iptr, data

    call assert_equals(.true., c_is_empty(my_queue), "Queue is empty")

    do i=1,10
      allocate(iptr, source=i)
      call c_push(my_queue, iptr)
      call assert_equals(i, c_size(my_queue), "Queue size increases as elements are pushed")
    end do

    call assert_equals(.false., c_is_empty(my_queue), "Queue is not empty after elements pushed")

    do j=1,10
      data => c_pop(my_queue)
      select type(data)
        type is (integer)
          call assert_equals(data, j, "Queue popped element is FIFO")
          call assert_equals(10-j, c_size(my_queue), "Queue pop removes element")
        class default
          call add_fail("Type unknown")
      end select
    end do

    call c_free(my_queue)
  end subroutine test_queue_type_push_pop

  ! Tests adding an element to the list_type and ensures that the list_type sizes up correctly
  subroutine test_list_type_add
    type(list_type) :: my_list
    integer :: i
    class(*), pointer :: iptr

    call assert_equals(.true., c_is_empty(my_list), "List is empty")

    do i=1,10
      allocate(iptr, source=i)
      call c_add(my_list, iptr)
      call assert_equals(i, c_size(my_list), "List add increases list_type size")
    end do

    call assert_equals(.false., c_is_empty(my_list), "List is not empty after element adds")

    call c_free(my_list)
  end subroutine test_list_type_add

  subroutine test_list_type_pointers
    type(list_type) :: my_list
    integer, target :: i
    class(*), pointer :: iptr, data
    i=20
    iptr=>i
    call c_add(my_list, iptr)
    i=50
    data=>c_get(my_list, 0)
    select type(data)
      type is (integer)
        call assert_equals(50, data, "List element modified by changing original value")
      class default
        call add_fail("Unknown type")
    end select

    call c_free(my_list)
  end subroutine test_list_type_pointers

  ! Adds a number of elements to the list_type and then checks each one to ensure that the value
  ! and order has not changed
  subroutine test_list_type_get
    type(list_type) :: my_list
    integer :: i,j
    class(*), pointer :: iptr, data

    do i=1,10
      allocate(iptr, source=i)
      call c_add(my_list, iptr)
    end do

    call assert_equals(10, c_size(my_list), "List size increased after adding elements")

    do j=1,10
      data => c_get(my_list, j)
      select type(data)
        type is(integer)
          call assert_equals(j, data, "Element at location j is consistent with expectations")
        class default
          call add_fail("Type unknown")
      end select
    end do

    call c_free(my_list)
  end subroutine test_list_type_get

  ! Creates a list_type with a number of elements and then inserts new elements into it. Checks all
  ! values at the end to ensure consistency
  subroutine test_list_type_insert
    type(list_type) :: my_list
    integer :: i, element_to_remove, j, k
    class(*), pointer :: iptr=>null(), data=>null()

    do i=1,10
      allocate(iptr, source=i)
      call c_add(my_list, iptr)
    end do

    do i=1,10
      j=i*100
      allocate(iptr, source=j)
      call c_insert(my_list, iptr, i)
    end do

    call assert_equals(20, c_size(my_list), "Post addition and insertion list_type size is correct")

    do k=1,20
      data => c_get(my_list, k)
      select type(data)
        type is (integer)
          call assert_equals(data, merge(k-10, k*100, k .gt. 10), "Element at k is consistent with addition and insertion")
        class default
          call add_fail("Unknown type")
      end select
    end do

    call c_free(my_list)
  end subroutine test_list_type_insert

  ! Will create a list_type and randomly determine an index to remove, will remove this and then check
  ! that it has gone and not affected any other elements in the list_type
  subroutine test_list_type_remove
    type(list_type) :: my_list
    integer :: i, j, element_to_remove
    class(*), pointer :: iptr=>null(), data=>null()
    real :: r

    call init_random_seed
    call random_number(r)
    element_to_remove = int(r*9)+1

    do i=1,10
      allocate(iptr, source=i)
      call c_add(my_list, iptr)
    end do

    data => c_get(my_list, element_to_remove)
    select type(data)
      type is (integer)
        call assert_equals(element_to_remove, data, "Element to remove is consistent")
      class default
        call add_fail("Unknown type")
    end select
    call c_remove(my_list, element_to_remove)

    do j=1,9
      if (j .ne. element_to_remove) then
        data => c_get(my_list, j)
        select type(data)
          type is (integer)
          call assert_equals(merge(j, j+1, j .lt. element_to_remove), data, "After element removal element at j is consistent")
          class default
            call add_fail("Unknown type")
        end select
      end if
    end do

    call c_free(my_list)
  end subroutine test_list_type_remove

  ! Helper subroutine to initialise the random seed (based on the clock)
  subroutine init_random_seed()
    integer :: i, n, clock
    integer, dimension(:), allocatable :: seed

    call random_seed(size = n)
    allocate(seed(n))

    call system_clock(count=clock)

    seed = clock + 37 * (/ (i - 1, i = 1, n) /)
    call random_seed(put=seed)

    deallocate(seed)
  end subroutine init_random_seed

end module test_collections_mod

! Driver for the collections_mod unit tests
program test_collections_driver
  use fruit, only : init_fruit, run_test_case, fruit_summary
  use test_collections_mod, only : test_list_type_add, test_list_type_get, test_list_type_remove, test_list_type_insert, &
       test_list_type_pointers, test_stack_type_push_pop, test_stack_type_pointers, test_queue_type_push_pop, &
       test_queue_type_pointers, test_map_type_put, test_map_type_put_unique, test_map_type_get, test_map_type_remove, &
       test_map_type_pointers

  implicit none

  call init_fruit
  call run_test_case(test_list_type_add, "List addition")
  call run_test_case(test_list_type_get, "List retrieval")
  call run_test_case(test_list_type_remove, "List removal")
  call run_test_case(test_list_type_insert, "List insertion")
  call run_test_case(test_list_type_pointers, "List pointer consistency")
  call run_test_case(test_stack_type_push_pop, "Stack push and pop")
  call run_test_case(test_stack_type_pointers, "Stack pointer consistency")
  call run_test_case(test_queue_type_push_pop, "Queue push and pop")
  call run_test_case(test_queue_type_pointers, "Queue pointer consistency")
  call run_test_case(test_map_type_put, "map_type put")
  call run_test_case(test_map_type_put_unique, "map_type unique key property")
  call run_test_case(test_map_type_get, "map_type retrieval")
  call run_test_case(test_map_type_remove, "map_type removal")
  call run_test_case(test_map_type_pointers, "map_type pointer consistency")
  call fruit_summary
end program test_collections_driver
