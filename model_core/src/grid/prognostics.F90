!> Contains prognostic field definitions and functions
module prognostics_mod
  use grids_mod, only : PRIMAL_GRID
  use datadefn_mod, only : DEFAULT_PRECISION
  use mpi, only : MPI_REQUEST_NULL
  implicit none

#ifndef TEST_MODE
  private
#endif

  !> A prognostic field which is assumed to be 3D
  type, public :: prognostic_field_type
    logical :: active = .false., halo_allocated=.false.
    integer, dimension(3) :: grid  !< The grid that this prognostic lives on per dimension
    real(kind=DEFAULT_PRECISION), dimension(:,:,:), allocatable :: data !< The 3D data associated with this field
    real(kind=DEFAULT_PRECISION), dimension(:), allocatable :: flux_previous_y, flux_y_buffer
    real(kind=DEFAULT_PRECISION), dimension(:,:), allocatable :: flux_previous_x
    integer :: async_flux_handle = MPI_REQUEST_NULL
  end type prognostic_field_type

  !> A pointer to the prognostic field. This is so we can wrap prognostics up in an array
  !! and still refer to the origonal field rather than a copy
  type, public :: prognostic_field_ptr_type
     type(prognostic_field_type), pointer :: ptr
  end type prognostic_field_ptr_type

  public get_field_interpolation_index
contains

  !> Retrieves the index(s) that require interpolation to go from the primal the dual grid
  !! @param field The field to check interpolation for
  function get_field_interpolation_index(field)
    type(prognostic_field_type), intent(inout) :: field
    logical, dimension(3) :: get_field_interpolation_index

    integer :: i
    do i=1,3
      get_field_interpolation_index(i) = merge(.true., .false., field%grid(i) == PRIMAL_GRID)
    end do
  end function get_field_interpolation_index
end module prognostics_mod
