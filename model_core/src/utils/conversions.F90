!> Conversion between common inbuilt FORTRAN data types
!!
!! The user will still need to supply conversions between their derived types but this makes it
!! easier when handling common inbuilt type conversions.
module conversions_mod
  use datadefn_mod, only : DEFAULT_PRECISION, SINGLE_PRECISION, DOUBLE_PRECISION
  implicit none

#ifndef TEST_MODE
  private
#endif

  ! This is the rounding applied when going from single to double precision numbers
  integer, parameter :: REAL_ROUNDING_PRECISION=int(1e8)

  !> Converts a data type into the generic (class *) form.
  !!
  !! Will convert structured data into its generic form and return a pointer to this. The caller should also
  !! specify whether to make a copy of the data or not
  !! @param data The structured data to convert into generic
  !! @param copyflag Whether to use a copy of the structured data or not
  !! @returns A pointer to the generic representation of the data
  interface conv_to_generic
    module procedure string_to_generic, integer_to_generic, real_single_to_generic, real_double_to_generic, logical_to_generic
  end interface conv_to_generic

  !> Converts data types to strings.
  !!
  !! For the generic input then a flag indicating whether to make a copy of the underlying generic data and
  !! the length of the resulting string must be specified. For other conversions just the data is needed as
  !! the string length is assumed and the returned string does not point to the data provided
  !! @param data The data to convert into a string
  !! @param copyflag For generic data only: Whether to use a copy of the structured data or not
  !! @param length For generic data only: Length of the resulting string
  !! @returns A string. For generic data a pointer to the string or null if generic conversion not possible
  interface conv_to_string
    module procedure generic_to_string, integer_to_string, real_single_to_string, real_double_to_string, logical_to_string
  end interface conv_to_string

  !> Converts data types to integers
  !!
  !! For the generic input then a flag indicating whether to make a copy of the underlying generic data or not
  !! is required. For all other data this copy is made automatically and no flag is required
  !! @param data The data to convert into an integer
  !! @param copyflag For generic data only: Whether to use a copy of the structured data or not
  !! @returns An integer. For generic data a pointer to the integer or null if generic conversion not possible
  interface conv_to_integer
    module procedure generic_to_integer, string_to_integer, real_to_integer, logical_to_integer
  end interface conv_to_integer

  !> Converts data types to real
  !!
  !! For the generic input then a flag indicating whether to make a copy of the underlying generic data or not
  !! is required. For all other data this copy is made automatically and no flag is required
  !! @param data The data to convert into a real
  !! @param copyflag For generic data only: Whether to use a copy of the structured data or not
  !! @returns A real. For generic data a pointer to the real or null if generic conversion not possible
  interface conv_to_real
    module procedure generic_to_real, string_to_real, integer_to_real, logical_to_real
  end interface conv_to_real

  !> Converts data types to logical
  !!
  !! For the generic input then a flag indicating whether to make a copy of the underlying generic data or not
  !! is required. For all other data this copy is made automatically and no flag is required
  !! @param data The data to convert into a logical
  !! @param copyflag For generic data only: Whether to use a copy of the structured data or not
  !! @returns A logical. For generic data a pointer to the logical or null if generic conversion not possible
  interface conv_to_logical
    module procedure generic_to_logical, string_to_logical, integer_to_logical, real_to_logical
  end interface conv_to_logical

  !> Determines whether a data item can be represented as an integer or not
  !!
  !! Currently works only for strings, reals always can be (strips the fractional component) and logicals
  !! go to one or zero
  !! @param data The data to test
  !! @returns Logical whether or not the data can be represented as an integer
  interface conv_is_integer
    module procedure string_is_integer
  end interface conv_is_integer

  !> Determines whether a data item can be represented as a real or not
  !!
  !! Currently works only for strings, integers always can be and logicals
  !! go to one or zero
  !! @param data The data to test
  !! @returns Logical whether or not the data can be represented as a real
  interface conv_is_real
    module procedure string_is_real
  end interface conv_is_real

  !> Determines whether a data item can be represented as a logical or not
  !!
  !! Currently works only for strings, reals and integers always can be (0=false, otherwise true)
  !! @param data The data to test
  !! @returns Logical whether or not the data can be represented as a logical
  interface conv_is_logical
    module procedure string_is_logical
  end interface conv_is_logical

  public conv_to_generic, conv_to_string, conv_to_integer, conv_to_real, conv_to_logical, &
            conv_is_integer, conv_is_real, conv_is_logical, conv_single_real_to_double, generic_to_double_real

contains

  !> Converts from a single to double precision real. This applies some rounding to a certain number of decimal
  !! places to ignore very small fractions
  !! @param input_real The single precision real to convert
  !! @returns Double precision representation which is smoothed to a specific rounding precision
  real(kind=DOUBLE_PRECISION) function conv_single_real_to_double(input_real)
    real(kind=SINGLE_PRECISION), intent(in) :: input_real

    conv_single_real_to_double=dnint(real(input_real, kind=DEFAULT_PRECISION) * &
         REAL_ROUNDING_PRECISION) / REAL_ROUNDING_PRECISION
  end function conv_single_real_to_double


  !> Determines whether a string is an integer or not
  !! @param string The string to test
  !! @returns Logical whether or not the string can be represented as an integer
  logical function string_is_integer(string)
    character(len=*), intent(in) :: string

    integer :: integer_value, ierr

    if (len(trim(string)) .ne. 0) then
      read(string, '(i10)', iostat=ierr ) integer_value
      string_is_integer = ierr == 0
    else
      string_is_integer=.false.
    end if
  end function string_is_integer

  !> Determines whether a string is a real or not
  !! @param string The string to test
  !! @returns Logical whether or not the string can be represented as a real
  logical function string_is_real(string)
    character(len=*), intent(in) :: string

    integer :: ierr
    real :: real_value

    if (len(trim(string)) .ne. 0) then
      read(string, '(f11.2)', iostat=ierr ) real_value
      string_is_real = ierr == 0
    else
      string_is_real=.false.
    end if
  end function string_is_real

  !> Determines whether a string is a logical or not
  !! @param string The string to test
  !! @returns Logical whether or not the string can be represented as a logical
  logical function string_is_logical(string)
    character(len=*), intent(in) :: string

    string_is_logical = .false.
    if (trim(adjustl(string)) .eq. "true" .or. trim(adjustl(string)) .eq. "false" .or. &
         trim(adjustl(string)) .eq. ".true." .or. trim(adjustl(string)) .eq. ".false." .or. &
         trim(adjustl(string)) .eq. ".true" .or. trim(adjustl(string)) .eq. "true." .or. &
         trim(adjustl(string)) .eq. ".false" .or. trim(adjustl(string)) .eq. "false.") string_is_logical = .true.
  end function string_is_logical

  !> Converts a generic to a string
  !! @param generic The generic to convert into a string
  !! @param makecopy Whether to use a copy of the generic data or not
  !! @param str_length Length of the resulting string
  !! @returns A pointer to the string or null if generic conversion not possible
  function generic_to_string(generic, makecopy, str_length)
    class(*), pointer, intent(in) :: generic
    logical, intent(in) :: makecopy
    integer, intent(in) :: str_length
    character(len=str_length), pointer :: generic_to_string, temporary_generic_ptr

    select type(generic)
      type is (character(len=*))
        if (makecopy) then
          ! Need to do this to enforce string length information
          temporary_generic_ptr=>generic
          allocate(generic_to_string, source=temporary_generic_ptr)
        else
          generic_to_string=>generic
        end if
      class default
        generic_to_string=>null()
    end select
  end function generic_to_string

  !> Converts an integer to a string
  !! @param input The integer to convert into a string
  !! @returns The string of length 15 characters
  function integer_to_string(input)
    integer, intent(in) :: input
    character(len=15) :: integer_to_string

    write(integer_to_string, '(i15)' ) input
    integer_to_string = trim(adjustl(integer_to_string))
  end function integer_to_string

  !> Converts a single precision real to a string
  !! @param input The real to convert into a string
  !! @returns The string of length 30 characters
  function real_single_to_string(input, decimal_places, exponent, exponent_small_numbers)
    real(kind=SINGLE_PRECISION), intent(in) :: input
    character(len=30) :: real_single_to_string
    integer, optional :: decimal_places
    logical, optional :: exponent, exponent_small_numbers

    logical :: transformed
    transformed=.false.

    if (present(exponent)) then
      if (exponent) then
        write(real_single_to_string, '(es30.10)' ) input
        transformed=.true.
      end if
    end if
    if (present(exponent_small_numbers)) then
      if (exponent_small_numbers) then
        write(real_single_to_string, '(g30.10)' ) input
        transformed=.true.
      end if
    end if
    if (.not. transformed) then
      write(real_single_to_string, '(f30.10)' ) input
      if (scan(real_single_to_string, "*") .ne. 0) write(real_single_to_string, '(es30.10)' ) input
    end if
    call trim_trailing_zeros(real_single_to_string, 2)
    if (present(decimal_places)) call limit_to_decimal_places(real_single_to_string, decimal_places)

    real_single_to_string = trim(adjustl(real_single_to_string))
  end function real_single_to_string

  !> Converts a double precision real to a string
  !! @param input The real to convert into a string
  !! @returns The string of length 30 characters
  function real_double_to_string(input, decimal_places, exponent, exponent_small_numbers)
    real(kind=DOUBLE_PRECISION), intent(in) :: input
    character(len=30) :: real_double_to_string
    integer, optional :: decimal_places
    logical, optional :: exponent, exponent_small_numbers

    logical :: transformed
    transformed=.false.

    if (present(exponent)) then
      if (exponent) then
        write(real_double_to_string, '(es30.10)' ) input
        transformed=.true.
      end if
    end if
    if (present(exponent_small_numbers)) then
      if (exponent_small_numbers) then
        write(real_double_to_string, '(g30.10)' ) input
        transformed=.true.
      end if
    end if
    if (.not. transformed) then
      write(real_double_to_string, '(f30.10)' ) input
      if (scan(real_double_to_string, "*") .ne. 0) write(real_double_to_string, '(es30.10)' ) input
    end if
    call trim_trailing_zeros(real_double_to_string, 2)
    if (present(decimal_places)) then
      call limit_to_decimal_places(real_double_to_string, decimal_places)
    end if

    real_double_to_string = trim(adjustl(real_double_to_string))
  end function real_double_to_string

  !> Helper subroutine which trims the string down to an upper limit of decimal places, with all
  !! numbers beyond this point removed
  !! @param stringToParse The raw, uncropped, string to processess which is modified
  !! @param decimalPlaces Number of decimal places to keep
  subroutine limit_to_decimal_places(string_to_parse, decimal_places)
    character(len=*), intent(inout) :: string_to_parse
    integer, intent(in) :: decimal_places

    integer :: decimal_posn, exp_posn

    string_to_parse=adjustl(string_to_parse)
    decimal_posn=index(string_to_parse, ".")
    exp_posn=index(string_to_parse, "E")
    if (decimal_posn .ne. 0 .and. decimal_posn+decimal_places+1 .le. len(string_to_parse)) then
      if (exp_posn .eq. 0) then
        string_to_parse(decimal_posn+decimal_places+1:)=" "
      else
        string_to_parse(decimal_posn+decimal_places+1:)=string_to_parse(exp_posn:)
        string_to_parse(decimal_posn+decimal_places+1+(len(string_to_parse)-exp_posn)+1:)=" "
      end if
    end if    
  end subroutine limit_to_decimal_places

  !> A helper subroutine which trims training zeros from the string after a decimal place
  !! this is to make the string more readable when printed out
  !! @param stringToParse The string to parse which is modified to replace trailing zeros
  !! @param zerosToRetain The number of trailing (after decimal) zeros to retain
  subroutine trim_trailing_zeros(string_to_parse, zeros_to_retain)
    character(len=*), intent(inout) :: string_to_parse
    integer, intent(in) :: zeros_to_retain

    integer :: decimal_posn, i, zero_count, nonzero_hit

    zero_count=0

    decimal_posn=index(string_to_parse, ".")
    if (decimal_posn .ne. 0 .and. decimal_posn .lt. len(string_to_parse)) then
      do i=len(trim(string_to_parse)), decimal_posn, -1
        if (string_to_parse(i:i) .ne. "0") then
          nonzero_hit=i
          exit
        else
          zero_count=zero_count+1
        end if
      end do
      if (zero_count .gt. zeros_to_retain) then
        string_to_parse(nonzero_hit+zeros_to_retain:)=""
      end if
    end if
  end subroutine trim_trailing_zeros  

  !> Converts a logical to a string
  !! @param input The logical to convert into a string
  !! @returns The string of length 5 characters
  function logical_to_string(input)
    logical, intent(in) :: input
    character(len=5) :: logical_to_string

    if (input) then
      logical_to_string = "true"
    else
      logical_to_string = "false"
    end if
  end function logical_to_string

  !> Converts a generic to a logical
  !! @param generic The generic to convert into a logical
  !! @param makecopy Whether to use a copy of the generic data or not
  !! @returns A pointer to the logical or null if generic conversion not possible
  function generic_to_logical(generic, makecopy)
    class(*), pointer, intent(in) :: generic
    logical, intent(in) :: makecopy
    logical, pointer :: generic_to_logical

    select type(generic)
      type is (logical)
        if (makecopy) then
          allocate(generic_to_logical, source=generic)
        else
          generic_to_logical=>generic
        end if
      class default
        generic_to_logical=>null()
    end select
  end function generic_to_logical

  !> Converts a string to a logical
  !! @param string The string to convert into a logical (case sensitive)
  !! @returns The logical
  logical function string_to_logical(string)
    character(len=*), intent(in) :: string

    if (trim(adjustl(string)) .eq. "true" .or. trim(adjustl(string)) .eq. ".true." .or. &
         trim(adjustl(string)) .eq. ".true" .or. trim(adjustl(string)) .eq. "true.") then
      string_to_logical = .true.
    else
      string_to_logical = .false.
    end if
  end function string_to_logical

  !> Converts an integer to a logical
  !! @param input The integer to convert into a logical
  !! @returns The logical
  logical function integer_to_logical(input)
    integer, intent(in) :: input

    if (input .ge. 1) then
      integer_to_logical = .true.
    else
      integer_to_logical = .false.
    end if
  end function integer_to_logical

  !> Converts a real to a logical
  !! @param input The real to convert into a logical
  !! @returns The logical
  logical function real_to_logical(input)
    real, intent(in) :: input

    if (input .ge. 1.0) then
      real_to_logical = .true.
    else
      real_to_logical = .false.
    end if
  end function real_to_logical

  !> Converts a generic to a double real
  !! @param generic The generic to convert into a double real
  !! @param makecopy Whether to use a copy of the generic data or not
  !! @returns A pointer to the double real or null if generic conversion not possible
  function generic_to_double_real(generic, makecopy)
    class(*), pointer, intent(in) :: generic
    logical, intent(in) :: makecopy
    real(kind=DEFAULT_PRECISION), pointer :: generic_to_double_real

    select type(generic)
      type is (real(kind=DEFAULT_PRECISION))
        if (makecopy) then
          allocate(generic_to_double_real, source=generic)
        else
          generic_to_double_real=>generic
        end if
      class default
        generic_to_double_real=>null()
    end select
  end function generic_to_double_real

  !> Converts a generic to a real. If this is infact an integer then will do a conversion and allocate pointed to this
  !! @param generic The generic to convert into a real
  !! @param makecopy Whether to use a copy of the generic data or not
  !! @returns A pointer to the real or null if generic conversion not possible
  function generic_to_real(generic, makecopy)
    class(*), pointer, intent(in) :: generic
    logical, intent(in) :: makecopy
    real, pointer :: generic_to_real

    select type(generic)
      type is (real)
        if (makecopy) then
          allocate(generic_to_real, source=generic)
        else
          generic_to_real=>generic
        end if
        type is (integer)          
          allocate(generic_to_real)        
          generic_to_real=conv_to_real(generic)        
      class default
        generic_to_real=>null()
    end select
  end function generic_to_real

  !> Converts a string to a real
  !! @param string The string to convert into a real
  !! @returns The real
  real function string_to_real(string)
    character(len=*), intent(in) :: string

    if (scan(string, "E") .ne. 0 .or. scan(string, "e") .ne. 0) then
      read(string, '(es30.10)' ) string_to_real
    else
      read(string, '(f10.0)' ) string_to_real
    end if
  end function string_to_real

  !> Converts an integer to a real
  !! @param input The integer to convert into a real
  !! @returns The real
  real function integer_to_real(input)
    integer, intent(in) :: input

    integer_to_real = real(input)
  end function integer_to_real

  !> Converts a logical to a real
  !! @param input The logical to convert into a real
  !! @returns The real
  real function logical_to_real(input)
    logical, intent(in) :: input

    if (input) then
      logical_to_real = 1.0
    else
      logical_to_real = 0.0
    end if
  end function logical_to_real

  !> Converts a generic to an integer
  !! @param generic The generic to convert into an integer
  !! @param makecopy Whether to use a copy of the generic data or not
  !! @returns A pointer to the integer or null if generic conversion not possible
  function generic_to_integer(generic, makecopy)
    class(*), pointer, intent(in) :: generic
    logical, intent(in) :: makecopy
    integer, pointer :: generic_to_integer

    select type(generic)
      type is (integer)
        if (makecopy) then
          allocate(generic_to_integer, source=generic)
        else
          generic_to_integer=>generic
        end if
      class default
        generic_to_integer=>null()
    end select
  end function generic_to_integer

  !> Converts a string to an integer
  !! @param string The string to convert into an integer
  !! @returns The integer
  integer function string_to_integer(string)
    character(len=*), intent(in) :: string

    read(string, '(i15)' ) string_to_integer
  end function string_to_integer

  !> Converts a real to an integer
  !! @param input The real to convert into an integer
  !! @returns The integer
  integer function real_to_integer(input)
    real, intent(in) :: input

    real_to_integer = int(input)
  end function real_to_integer

  !> Converts a logical to an integer
  !! @param input The logical to convert into an integer
  !! @returns The integer
  integer function logical_to_integer(input)
    logical, intent(in) :: input

    if (input) then
      logical_to_integer = 1
    else
      logical_to_integer = 0
    end if
  end function logical_to_integer

  !> Converts a string into its generic data representation
  !! @param string The string to convert into its generic representation
  !! @param makecopy Whether make a copy of the underlying data or just return a simple pointer
  !! @returns A pointer to the generic data
  function string_to_generic(string, makecopy)
    character(len=*), target, intent(in) :: string
    logical, intent(in) :: makecopy
    class(*), pointer :: string_to_generic

    if (makecopy) then
      allocate(string_to_generic, source=string)
    else
      string_to_generic=>string
    end if
  end function string_to_generic

  !> Converts an integer into its generic data representation
  !! @param input The integer to convert into its generic representation
  !! @param makecopy Whether make a copy of the underlying data or just return a simple pointer
  !! @returns A pointer to the generic data
  function integer_to_generic(input, makecopy)
    integer, target , intent(in) :: input
    logical, intent(in) :: makecopy
    class(*), pointer :: integer_to_generic

    if (makecopy) then
      allocate(integer_to_generic, source=input)
    else
      integer_to_generic=>input
    end if
  end function integer_to_generic

  !> Converts a single real into its generic data representation
  !! @param input The real to convert into its generic representation
  !! @param makecopy Whether make a copy of the underlying data or just return a simple pointer
  !! @returns A pointer to the generic data
  function real_single_to_generic(input, makecopy)
    real(kind=SINGLE_PRECISION), target, intent(in) :: input
    logical, intent(in) :: makecopy
    class(*), pointer :: real_single_to_generic

    if (makecopy) then
      allocate(real_single_to_generic, source=input)
    else
      real_single_to_generic=>input
    end if
  end function real_single_to_generic

  !> Converts a double real into its generic data representation
  !! @param input The real to convert into its generic representation
  !! @param makecopy Whether make a copy of the underlying data or just return a simple pointer
  !! @returns A pointer to the generic data
  function real_double_to_generic(input, makecopy)
    real(kind=DOUBLE_PRECISION), target, intent(in) :: input
    logical, intent(in) :: makecopy
    class(*), pointer :: real_double_to_generic

    if (makecopy) then
      allocate(real_double_to_generic, source=input)
    else
      real_double_to_generic=>input
    end if
  end function real_double_to_generic

  !> Converts a logical into its generic data representation
  !! @param input The logical to convert into its generic representation
  !! @param makecopy Whether make a copy of the underlying data or just return a simple pointer
  !! @returns A pointer to the generic data
  function logical_to_generic(input, makecopy)
    logical, target, intent(in) :: input
    logical, intent(in) :: makecopy
    class(*), pointer :: logical_to_generic

    if (makecopy) then
      allocate(logical_to_generic, source=input)
    else
      logical_to_generic=>input
    end if
  end function logical_to_generic
end module conversions_mod
