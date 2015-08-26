! Dummy NetCDF module that can be used instead of the real NetCDF for testing purposes.
! Lots of the functionality is just stubbed out, but the aim is for memory to hold the specific
! values so that testing can be done by changing/inspecting these data structures
module dummy_netcdf_mod
  use collections_mod, only : map_type, c_put, c_contains, c_get, c_free
  use conversions_mod, only : conv_to_generic, conv_to_string, conv_to_integer
  use fruit, only : assert_false
  use datadefn_mod, only : DOUBLE_PRECISION, SINGLE_PRECISION
  implicit none

  ! Whether we are in define mode or not
  logical :: define_mode = .false.

  ! Allows us to wrap integer arrays for generic conversion
  type integer_array_wrapper_type
    integer size
    integer, allocatable, dimension(:) :: data
  end type integer_array_wrapper_type

  ! NetCDF state_mod - note these are global to all currently "open" files
  type(map_type), save :: variable_ids, variable_data, global_attributes, dimension_ids, dimension_lengths, var_data

  ! Parameters that the source code uses
  integer, parameter :: nf90_noerr = 0, dummy_error=100, NF90_GLOBAL=0, nf90_ebaddim=1000, nf90_enotatt=2000, &
    nf90_enotvar = 3000, NF90_REAL=0, NF90_INT=1, NF90_CHAR=2, NF90_DOUBLE=3, NF90_CLOBBER=0, nf90_nowrite=1, &
    NF90_NETCDF4=0, NF90_MPIIO=0

  ! Internal counters to keep track of the latest ids
  integer :: current_ncid=0, current_dim=0, current_var_id = 0

  ! Generic put to support a variety of data types
  interface nf90_put_var
    module procedure nf90_put_var_char, nf90_put_var_integer, nf90_put_var_real, nf90_put_var_real_3d, nf90_put_var_real_scalar, &
         nf90_put_var_double, nf90_put_var_double_3d, nf90_put_var_double_scalar
  end interface nf90_put_var

  ! Generic define to support different sizes of define calls
  interface nf90_def_var
    module procedure nf90_def_var_multiple, nf90_def_var_single, nf90_def_var_atomic
  end interface nf90_def_var

  ! Generic variable get to support different types of variable data
  interface nf90_get_var
    module procedure nf90_get_var_real3d, nf90_get_var_real, nf90_get_var_integer, nf90_get_var_char, nf90_get_var_real_scalar, &
         nf90_get_var_double3d, nf90_get_var_double, nf90_get_var_double_scalar
  end interface nf90_get_var
contains

  ! Resets the NetCDF state_mod (useful between tests)
  subroutine dummy_netcdf_reset()
    call c_free(variable_ids)
    call c_free(global_attributes)
    call c_free(dimension_ids)
    call c_free(dimension_lengths)
    current_ncid=0
    current_dim=0
    current_var_id = 0
    define_mode = .false.
  end subroutine dummy_netcdf_reset

  ! Stub for open a NetCDF file
  integer function nf90_open(path, mode, ncid)
    character(len=*), intent(in) :: path
    integer, intent(in) :: mode
    integer, intent(out) :: ncid

    current_ncid = current_ncid + 1
    ncid = current_ncid
    nf90_open = nf90_noerr
  end function nf90_open

  ! Stub for creating a NetCDF file
  integer function nf90_create(path, mode, ncid, comm, info)
    character(len=*), intent(in) :: path
    integer, intent(in) :: mode
    integer, intent(in), optional :: comm, info
    integer, intent(out) :: ncid

    ncid = 1
    define_mode = .true.
    nf90_create = nf90_noerr
  end function nf90_create

  ! Stub for putting an attribute into the NetCDF
  integer function nf90_put_att(ncid, attribute, key, value)
    integer, intent(in) :: ncid, attribute
    character(len=*), intent(in) :: key, value
    class(*), pointer :: raw_data

    if (attribute == NF90_GLOBAL) then
      raw_data=>conv_to_generic(value, .true.)
      call c_put(global_attributes, key, raw_data)
    else

    end if
    nf90_put_att = nf90_noerr
  end function nf90_put_att

  ! Stub for the enddef call
  integer function nf90_enddef(ncid)
    integer, intent(in) :: ncid

    define_mode = .false.
    nf90_enddef = nf90_noerr
  end function nf90_enddef

  ! Stub for closing the NetCDF file
  integer function nf90_close(ncid)
    integer, intent(in) :: ncid

    nf90_close = nf90_noerr
  end function nf90_close

  ! Stub for inquiring about an attributes size
  integer function nf90_inquire_attribute(ncid, attributeid, key, len)
    integer, intent(in) :: ncid, attributeid
    character(len=*), intent(in) :: key
    integer, intent(out) :: len

    len = 200 ! Dummy 200
    nf90_inquire_attribute = nf90_noerr
  end function nf90_inquire_attribute

  ! Stub for getting an attribute, currently works for global attributes
  integer function nf90_get_att(ncid, attributeid, key, value)
    integer, intent(in) :: ncid, attributeid
    character(len=*), intent(in) :: key
    character(len=100), intent(out) :: value

    class(*), pointer :: raw_data

    if (attributeid == NF90_GLOBAL) then
      raw_data=>c_get(global_attributes, key)
      value = conv_to_string(raw_data, .true., 100)
    end if
    nf90_get_att = nf90_noerr
  end function nf90_get_att

  ! Stub for grabbing the id of a variable's key
  integer function nf90_inq_varid(ncid, key, varid)
    character(len=*), intent(in) :: key
    integer, intent(in) :: ncid
    integer, intent(out) :: varid

    class(*), pointer :: raw_data

    if (c_contains(variable_ids, key)) then
      raw_data=>c_get(variable_ids, key)
      if (associated(raw_data)) then
        varid = conv_to_integer(raw_data, .false.)
        nf90_inq_varid = nf90_noerr
        return
      end if
      nf90_inq_varid = dummy_error
    else
      nf90_inq_varid = nf90_enotvar
    end if
  end function nf90_inq_varid

  ! Stub for getting character variable data
  integer function nf90_get_var_char(ncid, varid, target, indexes, start, count, map)
    integer, intent(in) :: ncid, varid
    character(len=100), intent(out) :: target
    integer, dimension(:), optional, intent(in) :: indexes, start, count, map

    nf90_get_var_char = nf90_noerr
  end function nf90_get_var_char

  ! Stub for getting integer variable data
  integer function nf90_get_var_integer(ncid, varid, target, indexes, start, count, map)
    integer, intent(in) :: ncid, varid
    integer, dimension(*), intent(out) :: target
    integer, dimension(:), optional, intent(in) :: indexes, start, count, map

    class(*), pointer :: raw_data
    raw_data=>c_get(var_data, conv_to_string(varid))
    target(1) = conv_to_integer(raw_data, .false.)

    nf90_get_var_integer = nf90_noerr
  end function nf90_get_var_integer

  ! Stub for getting real variable data
  integer function nf90_get_var_real(ncid, varid, target, indexes, start, count, map)
    integer, intent(in) :: ncid, varid
    real(kind=SINGLE_PRECISION), dimension(*), intent(out) :: target
    integer, dimension(:), optional, intent(in) :: indexes, start, count, map

    nf90_get_var_real = nf90_noerr
  end function nf90_get_var_real

  ! Stub for getting real scalar variable data
  integer function nf90_get_var_real_scalar(ncid, varid, target, indexes, start, count, map)
    integer, intent(in) :: ncid, varid
    real(kind=SINGLE_PRECISION), intent(out) :: target
    integer, dimension(:), optional, intent(in) :: indexes, start, count, map

    nf90_get_var_real_scalar = nf90_noerr
  end function nf90_get_var_real_scalar

  ! Stub for getting read array (rank 3) variable data
  integer function nf90_get_var_real3d(ncid, varid, target, indexes, start, count, map)
    integer, intent(in) :: ncid, varid
    real(kind=SINGLE_PRECISION), dimension(:,:,:), intent(out) :: target
    integer, dimension(:), optional, intent(in) :: indexes, start, count, map

    nf90_get_var_real3d = nf90_noerr
  end function nf90_get_var_real3d

  ! Stub for getting real variable data
  integer function nf90_get_var_double(ncid, varid, target, indexes, start, count, map)
    integer, intent(in) :: ncid, varid
    real(kind=DOUBLE_PRECISION), dimension(*), intent(out) :: target
    integer, dimension(:), optional, intent(in) :: indexes, start, count, map

    nf90_get_var_double = nf90_noerr
  end function nf90_get_var_double

  ! Stub for getting real scalar variable data
  integer function nf90_get_var_double_scalar(ncid, varid, target, indexes, start, count, map)
    integer, intent(in) :: ncid, varid
    real(kind=DOUBLE_PRECISION), intent(out) :: target
    integer, dimension(:), optional, intent(in) :: indexes, start, count, map

    nf90_get_var_double_scalar = nf90_noerr
  end function nf90_get_var_double_scalar

  ! Stub for getting read array (rank 3) variable data
  integer function nf90_get_var_double3d(ncid, varid, target, indexes, start, count, map)
    integer, intent(in) :: ncid, varid
    real(kind=DOUBLE_PRECISION), dimension(:,:,:), intent(out) :: target
    integer, dimension(:), optional, intent(in) :: indexes, start, count, map

    nf90_get_var_double3d = nf90_noerr
  end function nf90_get_var_double3d

  ! Stub for putting character variable data
  integer function nf90_put_var_char(ncid, varid, source, indexes, start, count, map)
    integer, intent(in) :: ncid, varid
    character(len=*), intent(in) :: source
    integer, dimension(:), optional , intent(in) :: indexes, start, count, map

    class(*), pointer :: raw_data

    call assert_false(define_mode, "Switched from define mode")
    raw_data=>conv_to_generic(source, .true.)
    call c_put(var_data, conv_to_string(varid), raw_data)

    nf90_put_var_char = nf90_noerr
  end function nf90_put_var_char

  ! Stub for putting integer variable data
  integer function nf90_put_var_integer(ncid, varid, source, indexes, start, count, map)
    integer, intent(in) :: ncid, varid
    integer, intent(in) :: source
    integer, dimension(:), optional , intent(in) :: indexes, start, count, map

    class(*), pointer :: raw_data

    call assert_false(define_mode, "Switched from define mode")
    raw_data=>conv_to_generic(source, .true.)
    call c_put(var_data, conv_to_string(varid), raw_data)
    nf90_put_var_integer = nf90_noerr
  end function nf90_put_var_integer

  ! Stub for putting real variable data
  integer function nf90_put_var_real(ncid, varid, source, indexes, start, count, map)
    integer, intent(in) :: ncid, varid
    real(kind=SINGLE_PRECISION), dimension(*), intent(in) :: source
    integer, dimension(:), optional , intent(in) :: indexes, start, count, map

    call assert_false(define_mode, "Switched from define mode")
    nf90_put_var_real = nf90_noerr
  end function nf90_put_var_real

  ! Stub for putting real scalar variable data
  integer function nf90_put_var_real_scalar(ncid, varid, source, indexes, start, count, map)
    integer, intent(in) :: ncid, varid
    real(kind=SINGLE_PRECISION), intent(in) :: source
    integer, dimension(:), optional , intent(in) :: indexes, start, count, map

    call assert_false(define_mode, "Switched from define mode")
    nf90_put_var_real_scalar = nf90_noerr
  end function nf90_put_var_real_scalar

  ! Stub for putting real array (rank 3) variable data
  integer function nf90_put_var_real_3d(ncid, varid, source, indexes, start, count, map)
    integer, intent(in) :: ncid, varid
    real(kind=SINGLE_PRECISION), dimension(:,:,:), intent(in) :: source
    integer, dimension(:), optional , intent(in) :: indexes, start, count, map

    call assert_false(define_mode, "Switched from define mode")
    nf90_put_var_real_3d = nf90_noerr
  end function nf90_put_var_real_3d

    ! Stub for putting real variable data
  integer function nf90_put_var_double(ncid, varid, source, indexes, start, count, map)
    integer, intent(in) :: ncid, varid
    real(kind=DOUBLE_PRECISION), dimension(*), intent(in) :: source
    integer, dimension(:), optional , intent(in) :: indexes, start, count, map

    call assert_false(define_mode, "Switched from define mode")
    nf90_put_var_double = nf90_noerr
  end function nf90_put_var_double

  ! Stub for putting real scalar variable data
  integer function nf90_put_var_double_scalar(ncid, varid, source, indexes, start, count, map)
    integer, intent(in) :: ncid, varid
    real(kind=DOUBLE_PRECISION), intent(in) :: source
    integer, dimension(:), optional , intent(in) :: indexes, start, count, map

    call assert_false(define_mode, "Switched from define mode")
    nf90_put_var_double_scalar = nf90_noerr
  end function nf90_put_var_double_scalar

  ! Stub for putting real array (rank 3) variable data
  integer function nf90_put_var_double_3d(ncid, varid, source, indexes, start, count, map)
    integer, intent(in) :: ncid, varid
    real(kind=DOUBLE_PRECISION), dimension(:,:,:), intent(in) :: source
    integer, dimension(:), optional , intent(in) :: indexes, start, count, map

    call assert_false(define_mode, "Switched from define mode")
    nf90_put_var_double_3d = nf90_noerr
  end function nf90_put_var_double_3d

  ! Stub for getting the corresponding id of a dimension name. If none is found returns the appropriate error code
  integer function nf90_inq_dimid(ncid, key, dim_id)
    integer, intent(in) :: ncid
    character(len=*), intent(in) :: key
    integer, intent(out) :: dim_id

    class(*), pointer :: raw_data

    if (c_contains(dimension_ids, key)) then
      raw_data=>c_get(dimension_ids, key)
      dim_id = conv_to_integer(raw_data, .false.)
      nf90_inq_dimid = nf90_noerr
    else
      nf90_inq_dimid = nf90_ebaddim
    end if
  end function nf90_inq_dimid

  ! Gets the length of a dimension from its id
  integer function nf90_inquire_dimension(ncid, id, len)
    integer, intent(in) :: ncid, id
    integer, intent(out) :: len

    class(*), pointer :: raw_data

    if (c_contains(dimension_lengths, conv_to_string(id))) then
      raw_data => c_get(dimension_lengths, conv_to_string(id))
      len = conv_to_integer(raw_data, .false.)
      nf90_inquire_dimension = nf90_noerr
    else
      nf90_inquire_dimension = dummy_error
    end if
  end function nf90_inquire_dimension

  ! Stub for getting the string of an error code
  character(len=10) function nf90_strerror(status)
    integer, intent(in) :: status

    nf90_strerror = "dummy"
  end function nf90_strerror

  ! Defines a dimension of a specific length
  integer function nf90_def_dim(ncid, key, length, dimension_id)
    integer, intent(in) :: ncid, length
    character(len=*), intent(in) :: key
    integer, intent(out) :: dimension_id

    class(*), pointer :: raw_data

    current_dim = current_dim + 1
    dimension_id = current_dim
    raw_data=>conv_to_generic(length, .true.)
    call c_put(dimension_lengths, conv_to_string(dimension_id), raw_data)
    nf90_def_dim = nf90_noerr
  end function nf90_def_dim

  ! Defines a variable with a single dimension
  integer function nf90_def_var_single(ncid, key, type, dim_id, varid)
    integer, intent(in) :: ncid, type, dim_id
    character(len=*), intent(in) :: key
    integer, intent(out) :: varid

    class(*), pointer :: raw_data

    current_var_id = current_var_id + 1
    varid = current_var_id
    raw_data=>conv_to_generic(varid, .true.)
    call c_put(variable_ids, key, raw_data)
    raw_data=>conv_to_generic(dim_id, .true.)
    call c_put(variable_data, key, raw_data)
    nf90_def_var_single = nf90_noerr
  end function nf90_def_var_single

  ! Defines a variable with no dimension (size=1, rank=1)
  integer function nf90_def_var_atomic(ncid, key, type, varid)
    integer, intent(in) :: ncid, type
    character(len=*), intent(in) :: key
    integer, intent(out) :: varid

    class(*), pointer :: raw_data

    current_var_id = current_var_id + 1
    varid = current_var_id
    raw_data=>conv_to_generic(varid, .true.)
    call c_put(variable_ids, key, raw_data)
    raw_data=>conv_to_generic(1, .true.)
    call c_put(variable_data, key, raw_data)
    nf90_def_var_atomic = nf90_noerr
  end function nf90_def_var_atomic

  ! Defines a variable with multiple dimensions
  integer function nf90_def_var_multiple(ncid, key, type, dim_ids, varid)
    integer, intent(in) :: ncid, type
    integer, dimension(:) :: dim_ids
    character(len=*), intent(in) :: key
    integer, intent(out) :: varid

    integer :: i
    type(integer_array_wrapper_type), pointer :: wrapper
    class(*), pointer :: raw_data, raw_id_data

    allocate(wrapper)
    wrapper%size=size(dim_ids)
    allocate(wrapper%data(wrapper%size))

    do i=1,wrapper%size
      wrapper%data(i) = dim_ids(i)
    end do
    raw_data => wrapper
    current_var_id = current_var_id + 1
    varid = current_var_id
    raw_id_data=>conv_to_generic(varid, .true.)
    call c_put(variable_ids, key, raw_id_data)
    call c_put(variable_data, key, raw_data)
    nf90_def_var_multiple = nf90_noerr
  end function nf90_def_var_multiple
end module dummy_netcdf_mod
