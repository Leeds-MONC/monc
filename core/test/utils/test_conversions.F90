! Tests the conversion aspect of the core utilities. This covers all conversion procedures as each test goes to and from
! hence testing two bits of functionality
module test_conversions_mod
  use fruit, only : assert_true, assert_false, assert_not_equals, assert_equals
  use conversions_mod, only : conv_is_integer, conv_is_real, conv_is_logical, conv_to_integer, conv_to_real, &
       conv_to_logical, conv_to_string, conv_to_generic
  implicit none

  contains

  ! Tests the string is an integer functionality
  subroutine test_is_integer()
    call assert_true(conv_is_integer("1"), "Small number")
    call assert_true(conv_is_integer("9874321"), "Large number")
    call assert_true(conv_is_integer("-1024"), "Negative number")
    call assert_true(conv_is_integer("+986"), "Positive number")
    call assert_false(conv_is_integer("1.2"), "Floating point")
    call assert_false(conv_is_integer("ABCD"), "Non number")
    call assert_false(conv_is_integer("87T"), "Numeric and non number characters")
    call assert_false(conv_is_integer("0xAB"), "Hexadecimal")
  end subroutine test_is_integer

  ! Tests the string is a real functionality
  subroutine test_is_real()
    call assert_true(conv_is_real("1"), "Small number")
    call assert_true(conv_is_real("9874321"), "Large number")
    call assert_true(conv_is_real("-1024"), "Negative number")
    call assert_true(conv_is_real("+986"), "Positive number")
    call assert_true(conv_is_real("1.2"), "Floating point")
    call assert_true(conv_is_real("45432.2343"), "Large floating point")
    call assert_true(conv_is_real("1e2"), "Exponent floating point")
    call assert_true(conv_is_real("1.56e+6"), "Large positive exponent floating point")
    call assert_true(conv_is_real("765.98e-6"), "Small negative exponent floating point")
    call assert_false(conv_is_real("ABCD"), "Non number")
    call assert_false(conv_is_real("87T"), "Numeric and non numeric characters")
    call assert_false(conv_is_real("0xAB"), "Hexadecimal")
  end subroutine test_is_real

  ! Tests the string is a logical functionality
  subroutine test_is_logical()
    call assert_true(conv_is_logical("true"), "True value")
    call assert_true(conv_is_logical("false"), "False value")
    call assert_false(conv_is_logical("1"), "One number")
    call assert_false(conv_is_logical("0"), "Zero number")
    call assert_false(conv_is_logical("dsfdsfsd"), "Random characters")
    call assert_false(conv_is_logical("truexyz"), "Append characters to true")
  end subroutine test_is_logical

  ! Tests conversion from real to integer and back again
  subroutine test_real_to_integer()
    real :: test_real, retrieve_real
    integer :: retrieval_int

    test_real=19
    retrieval_int = conv_to_integer(test_real)
    retrieve_real = conv_to_real(retrieval_int)
    call assert_equals(retrieve_real, test_real, "Reals after conversion to and from integer are equal")
  end subroutine test_real_to_integer

  ! Tests conversion from logical to integer and then back again
  subroutine test_logical_to_integer()
    logical :: test_logical, retrieve_logical
    integer :: retrieval_int

    test_logical=.true.
    retrieval_int = conv_to_integer(test_logical)
    retrieve_logical = conv_to_logical(retrieval_int)
    call assert_equals(retrieve_logical, test_logical, "Logicals after conversion to and from integer are equal")
  end subroutine test_logical_to_integer

  ! Tests conversion from logical to real and back again
  subroutine test_logical_to_real()
    logical :: test_logical, retrieve_logical
    real :: retrieval_real

    test_logical=.true.
    retrieval_real = conv_to_real(test_logical)
    retrieve_logical = conv_to_logical(retrieval_real)
    call assert_equals(retrieve_logical, test_logical, "Logicals after conversion to and from real are equal")
  end subroutine test_logical_to_real

  ! Tests conversion from integer to string and bacl
  subroutine test_integer_to_string()
    integer :: test_int, retrieve_int
    character(len=15) :: retrieval_string

    test_int = 92
    retrieval_string = conv_to_string(test_int)
    retrieve_int = conv_to_integer(retrieval_string)
    call assert_equals(retrieve_int, test_int, "Integers after conversion to and from string are equal")
  end subroutine test_integer_to_string

  ! Tests conversion from real to string and back
  subroutine test_real_to_string()
    real :: test_real, retrieve_real, diff
    character(len=30) :: retrieval_string

    test_real = 63.45
    retrieval_string = conv_to_string(test_real)
    retrieve_real = conv_to_real(retrieval_string)
    ! Conversion to string is not exact, therefore take the difference and ensure it is within permissable bounds
    diff = test_real - retrieve_real    
    call assert_true(diff .gt. -0.0001 .and. diff .lt. 0.0001, "Reals after conversion to and from string are equal")
  end subroutine test_real_to_string

  ! Tests conversion from logical to string and back
  subroutine test_logical_to_string()
    logical :: test_logical, retrieve_logical
    character(len=5) :: retrieval_string

    test_logical = .true.
    retrieval_string = conv_to_string(test_logical)
    retrieve_logical = conv_to_logical(retrieval_string)
    call assert_equals(retrieve_logical, test_logical, "Logicals after conversion to and from string are equal")
  end subroutine test_logical_to_string

  ! Tests generic conversion from string and back, ensuring that with false supplied as copy arguments then both point to the same
  ! bit of memory and hence changing one will change the other
  subroutine test_string_to_generic()
    character(len=100) :: test_string
    character(len=100), pointer :: retrieval_string
    class(*), pointer :: generic_data

    integer :: i
    do i=1,100
      test_string(i:i)='C'
    end do

    generic_data => conv_to_generic(test_string, .false.)
    call assert_true(associated(generic_data), "Generic data not null")
    retrieval_string => conv_to_string(generic_data, .false., 100)
    call assert_true(associated(retrieval_string), "Retrieved string not null")
    call assert_equals(test_string, retrieval_string, "To and from generic strings are equal")
    retrieval_string(4:6) = "LKJ"
    call assert_equals(test_string, retrieval_string, "To and from generic strings are equal after one modified")
  end subroutine test_string_to_generic

  ! Tests conversion from string to generic and that specifying copy (in the to or from) will make a separate
  ! copy of the memory so that changing one will not affect the other
  subroutine test_string_to_generic_copy()
    character(len=100) :: test_string
    character(len=100), pointer :: retrieval_string
    class(*), pointer :: generic_data

    integer :: i
    do i=1,100
      test_string(i:i)='C'
    end do

    generic_data => conv_to_generic(test_string, .true.)
    call assert_true(associated(generic_data), "Generic data not null")
    retrieval_string => conv_to_string(generic_data, .false., 100)
    call assert_true(associated(retrieval_string), "Retrieved string not null")
    retrieval_string(4:6) = "LKJ"
    call assert_not_equals(test_string, retrieval_string, "To and from generic strings are different")

    generic_data => conv_to_generic(test_string, .false.)
    retrieval_string => conv_to_string(generic_data, .true., 100)
    call assert_true(associated(generic_data), "Generic data not null")
    retrieval_string(4:6) = "LKJ"
    call assert_not_equals(test_string, retrieval_string, "To and from generic strings are different")
  end subroutine test_string_to_generic_copy

  ! Tests generic conversion from integer and back, ensuring that with false supplied as copy arguments then both point to the same
  ! bit of memory and hence changing one will change the other
  subroutine test_integer_to_generic()
    integer :: test_int
    integer, pointer :: retrieval_int
    class(*), pointer :: generic_data

    test_int = 72

    generic_data => conv_to_generic(test_int, .false.)
    call assert_true(associated(generic_data), "Generic data not null")
    retrieval_int => conv_to_integer(generic_data, .false.)
    call assert_true(associated(retrieval_int), "Retrieved integer not null")
    call assert_equals(test_int, retrieval_int, "To and from generic integers are equal")
    retrieval_int = 13
    call assert_equals(test_int, retrieval_int, "To and from generic integers are equal after modification")
  end subroutine test_integer_to_generic

  ! Tests conversion from integer to generic and that specifying copy (in the to or from) will make a separate
  ! copy of the memory so that changing one will not affect the other
  subroutine test_integer_to_generic_copy()
    integer :: test_int
    integer, pointer :: retrieval_int
    class(*), pointer :: generic_data

    test_int = 72

    generic_data => conv_to_generic(test_int, .true.)
    call assert_true(associated(generic_data), "Generic data not null")
    retrieval_int => conv_to_integer(generic_data, .false.)
    call assert_true(associated(retrieval_int), "Retrieved integer not null")
    retrieval_int = 19
    call assert_not_equals(test_int, retrieval_int, "To and from generic integers are different")

    generic_data => conv_to_generic(test_int, .false.)
    retrieval_int => conv_to_integer(generic_data, .true.)
    call assert_true(associated(generic_data), "Generic data not null")
    retrieval_int = 32
    call assert_not_equals(test_int, retrieval_int, "To and from generic integers are different")
  end subroutine test_integer_to_generic_copy

  ! Tests generic conversion from real and back, ensuring that with false supplied as copy arguments then both point to the same
  ! bit of memory and hence changing one will change the other
  subroutine test_real_to_generic()
    real :: test_real
    real, pointer :: retrieval_real
    class(*), pointer :: generic_data

    test_real = 72.92

    generic_data => conv_to_generic(test_real, .false.)
    call assert_true(associated(generic_data), "Generic data not null")
    retrieval_real => conv_to_real(generic_data, .false.)
    call assert_true(associated(retrieval_real), "Retrieved real not null")
    call assert_equals(test_real, retrieval_real, "To and from generic reals are equal")
    retrieval_real = 13.1
    call assert_equals(test_real, retrieval_real, "To and from generic real are equal after modification")
  end subroutine test_real_to_generic

  ! Tests conversion from real to generic and that specifying copy (in the to or from) will make a separate
  ! copy of the memory so that changing one will not affect the other
  subroutine test_real_to_generic_copy()
    real :: test_real
    real, pointer :: retrieval_real
    class(*), pointer :: generic_data

    test_real = 72.92

    generic_data => conv_to_generic(test_real, .true.)
    call assert_true(associated(generic_data), "Generic data not null")
    retrieval_real => conv_to_real(generic_data, .false.)
    call assert_true(associated(retrieval_real), "Retrieved real not null")
    retrieval_real = 19
    call assert_not_equals(test_real, retrieval_real, "To and from generic reals are different")

    generic_data => conv_to_generic(test_real, .false.)
    call assert_true(associated(generic_data), "Generic data not null")
    retrieval_real => conv_to_real(generic_data, .true.)
    retrieval_real = 32
    call assert_not_equals(test_real, retrieval_real, "To and from generic reals are different")
  end subroutine test_real_to_generic_copy

  ! Tests generic conversion from logical and back, ensuring that with false supplied as copy arguments then both point to the same
  ! bit of memory and hence changing one will change the other
  subroutine test_logical_to_generic()
    logical :: test_logical
    logical, pointer :: retrieval_logical
    class(*), pointer :: generic_data

    test_logical = .true.

    generic_data => conv_to_generic(test_logical, .false.)
    call assert_true(associated(generic_data), "Generic data not null")
    retrieval_logical => conv_to_logical(generic_data, .false.)
    call assert_true(associated(retrieval_logical), "Retrieved logical not null")
    call assert_equals(test_logical, retrieval_logical, "To and from generic logicals are equal")
    retrieval_logical = .false.
    call assert_equals(test_logical, retrieval_logical, "To and from generic logicals are equal after modification")
  end subroutine test_logical_to_generic

  ! Tests conversion from logical to generic and that specifying copy (in the to or from) will make a separate
  ! copy of the memory so that changing one will not affect the other
  subroutine test_logical_to_generic_copy()
    logical :: test_logical
    logical, pointer :: retrieval_logical
    class(*), pointer :: generic_data

    test_logical = .true.

    generic_data => conv_to_generic(test_logical, .true.)
    call assert_true(associated(generic_data), "Generic data not null")
    retrieval_logical => conv_to_logical(generic_data, .false.)
    call assert_true(associated(retrieval_logical), "Retrieved logical not null")
    retrieval_logical = .false.
    call assert_not_equals(test_logical, retrieval_logical, "To and from generic logicals are different")

    generic_data => conv_to_generic(test_logical, .false.)
    call assert_true(associated(generic_data), "Generic data not null")
    retrieval_logical => conv_to_logical(generic_data, .true.)
    retrieval_logical = .false.
    call assert_not_equals(test_logical, retrieval_logical, "To and from generic logicals are different")
  end subroutine test_logical_to_generic_copy

end module test_conversions_mod

! Driver for conversion utility tests
program test_conversion_driver
  use fruit, only : init_fruit, run_test_case, fruit_summary
  use test_conversions_mod, only : test_string_to_generic, test_string_to_generic_copy, test_integer_to_generic, &
       test_integer_to_generic_copy, test_real_to_generic, test_real_to_generic_copy, test_logical_to_generic, &
       test_logical_to_generic_copy, test_integer_to_string, test_real_to_string, test_logical_to_string, &
       test_real_to_integer, test_logical_to_integer, test_logical_to_real, test_is_integer, test_is_real, test_is_logical

  implicit none

  call init_fruit
  call run_test_case(test_string_to_generic, "Test string generic conversion")
  call run_test_case(test_string_to_generic_copy, "Test string generic conversion copy")
  call run_test_case(test_integer_to_generic, "Test integer generic conversion")
  call run_test_case(test_integer_to_generic_copy, "Test integer generic conversion copy")
  call run_test_case(test_real_to_generic, "Test real generic conversion")
  call run_test_case(test_real_to_generic_copy, "Test real generic conversion copy")
  call run_test_case(test_logical_to_generic, "Test logical generic conversion")
  call run_test_case(test_logical_to_generic_copy, "Test logical generic conversion copy")
  call run_test_case(test_integer_to_string, "Test integer to string conversion")
  call run_test_case(test_real_to_string, "Test real to string conversion")
  call run_test_case(test_logical_to_string, "Test logical to string conversion")
  call run_test_case(test_real_to_integer, "Test real to integer conversion")
  call run_test_case(test_logical_to_integer, "Test logical to integer conversion")
  call run_test_case(test_logical_to_real, "Test logical to real conversion")
  call run_test_case(test_is_integer, "Test is string an integer")
  call run_test_case(test_is_real, "Test is string a real")
  call run_test_case(test_is_logical, "Test is string a logical")
  call fruit_summary
end program test_conversion_driver
