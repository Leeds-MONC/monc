!> Performs the interpolation between the primal and dual grids via a stencil approach. For performance reasons, for each
!! field we store the entirety of the y dimension and the number of x slices required by the stencil. Therefore a new
!! interpolation is simply the calculation of one point and reuse of existing computed points (unless this is the first x or y.)
!! The applicable interpolation stenciled data is copied out, with :,1,1 being the central point, minus and plus in each y and x
!! dimension as determined by the stencil size. This is done so that the stencil size can easily be different for each flow field
!! and this is the case with u (where we need u-2 in the X.)
module stencil_mod
  use prognostics_mod, only : prognostic_field_type, prognostic_field_ptr_type, get_field_interpolation_index
  use grids_mod, only : local_grid_type, X_INDEX, Y_INDEX, Z_INDEX
  use datadefn_mod, only : DEFAULT_PRECISION
  implicit none

#ifndef TEST_MODE
  private
#endif

  !> Configuration for a specific stencil interpolation to perform
  type, public :: grid_stencil_type
    type(prognostic_field_ptr_type), dimension(:), allocatable :: fields
    type(prognostic_field_type), dimension(:,:), allocatable :: interpolated_fields
    integer, dimension(:,:), allocatable :: sizes
    integer :: nfields, max_y_point, max_x_point
  end type grid_stencil_type

  public create_stencil, interpolate_to_dual, free_stencil
contains

  !> Creates a stencil configuration which will then be used for interpolation
  !! @param local_grid The local grid
  !! @param fields The fields to perform the interpolation upon
  !! @param nfields The number of fields included in a specific interpolation
  !! @param interpolations_to_perform The number of distinct interpolations that will be performed
  !! @param sizes The size of the stencil, for each field, in X and Y dimensions
  !! @param xdim Whether the X dimension is active
  !! @param ydim Whether the Y dimension is active
  type(grid_stencil_type) function create_stencil(local_grid, fields, nfields, interpolations_to_perform, sizes, xdim, ydim)
    type(local_grid_type), intent(inout) :: local_grid
    type(prognostic_field_ptr_type), dimension(:), intent(inout) :: fields
    integer, dimension(:,:), intent(in) :: sizes
    integer, intent(in) :: nfields, interpolations_to_perform
    logical, intent(in) :: xdim, ydim

    integer :: j,i

    allocate(create_stencil%fields(0:nfields), create_stencil%interpolated_fields(interpolations_to_perform,nfields))
    allocate(create_stencil%sizes(nfields,2))
    do j=1,nfields
      create_stencil%fields(j) = fields(j)
      create_stencil%sizes(j,1) = merge(sizes(j,1), 0, xdim)
      create_stencil%sizes(j,2) = merge(sizes(j,2), 0, ydim)
      do i=1,interpolations_to_perform
        allocate(create_stencil%interpolated_fields(i, j)%data(local_grid%size(Z_INDEX), &
             local_grid%size(Y_INDEX)+4, create_stencil%sizes(j,1)*2+1))
      end do
    end do
    create_stencil%nfields=nfields
    create_stencil%max_y_point=(local_grid%size(Y_INDEX)+local_grid%halo_size(Y_INDEX) * 2) - 1
    create_stencil%max_x_point=(local_grid%size(X_INDEX)+local_grid%halo_size(X_INDEX) * 2) - 1
  end function create_stencil

  !> Frees up the memory allocated to a stencil
  !! @param stencil The stencil that we will free up
  subroutine free_stencil(stencil)
    type(grid_stencil_type), intent(inout) :: stencil

    integer :: j, i

    if (allocated(stencil%interpolated_fields)) then
      do i=1, size(stencil%interpolated_fields,1)
        do j=1, size(stencil%interpolated_fields,2)
          if (allocated(stencil%interpolated_fields(i,j)%data)) deallocate(stencil%interpolated_fields(i,j)%data)
        end do
      end do
      deallocate(stencil%fields, stencil%interpolated_fields)
    end if
    if (allocated(stencil%sizes)) deallocate(stencil%sizes)
  end subroutine free_stencil  

  !> Interpolates the (vector) flow fields from the primal to dual grid based upon a specific field interpolation
  !! @param local_grid The local grid
  !! @param field The field to base this interpolation upon
  !! @param stencil Stencil configuration
  !! @param x Local column index in X
  !! @param y Local column index in Y
  !! @param interpolated_fields Interpolated prognostics that this will write the output into, central (x,y) point is 1,1
  !! @param interpolation_id Id of the interpolation
  subroutine interpolate_to_dual(local_grid, field, stencil, x, y, interpolated_fields, interpolation_id)
    type(prognostic_field_type), intent(inout) :: field
    type(local_grid_type), intent(inout) :: local_grid
    type(grid_stencil_type), intent(inout) :: stencil
    type(prognostic_field_type), dimension(:), allocatable, intent(inout) :: interpolated_fields
    integer, intent(in) :: x, y, interpolation_id

    integer :: i,j,n, x_top, base_x, y_start, y_top, x_start, abs_x_top
    logical, dimension(3) :: interpolate_in_dimension

    interpolate_in_dimension = get_field_interpolation_index(field)

    do n=1,stencil%nfields
      base_x=x-(stencil%sizes(n,1)+1)
      abs_x_top=x+stencil%sizes(n,1)
      x_start=1
      if (base_x .lt. 0) x_start=x_start+(0-base_x)
      x_top=stencil%sizes(n,1)*2+1
      y_start=y-stencil%sizes(n,2)
      y_top=y+stencil%sizes(n,2)
      if (y_start .lt. 1) y_start=1

      if (x==2) then
        ! If this is the first x slice then compute all points on y dimension up to x top - 1
        do i=x_start,x_top-1
          do j=y_start, y_top
              call calculate_interpolated_cell_value(local_grid, j, base_x+i, interpolate_in_dimension, stencil%fields(n), &
                   stencil%interpolated_fields(interpolation_id,n)%data(:,j,i))            
          end do
        end do
      end if

      if (y==2) then
        ! If this is the first column in y then shift X down (as only store required data) and calculate y to y top -1
        if (x/=2) then
          do i=1, x_top-1
            stencil%interpolated_fields(interpolation_id,n)%data(:,:,i)=&
                 stencil%interpolated_fields(interpolation_id,n)%data(:,:,i+1)
          end do
        end if
        do j=y_start, y_top-1
          call calculate_interpolated_cell_value(local_grid, j, abs_x_top, interpolate_in_dimension, stencil%fields(n), &
               stencil%interpolated_fields(interpolation_id,n)%data(:,j,x_top))
        end do
      end if
      ! Calculate y top, x top which is done for each column
      call calculate_interpolated_cell_value(local_grid, y_top, abs_x_top, interpolate_in_dimension, stencil%fields(n), &
           stencil%interpolated_fields(interpolation_id,n)%data(:,y_top,x_top))
    end do
    call copy_interpolated_to_prognostic(interpolated_fields, stencil, local_grid, x, y, interpolation_id)
  end subroutine interpolate_to_dual

  !> Calculates the interpolated values for each point in a column. The exact method of interpolation depends on whether we
  !! are interpolating in the w, y or x dimension.
  !! @param local_grid The local grid
  !! @param j Current local index in the y dimension
  !! @param i Current local index in the x dimension
  !! @param interpolate_in_dimension Which dimensions to interpolate in
  !! @param field The field to interpolate from
  !! @param column_stencil_values The result of the interpolation for the specific column
  subroutine calculate_interpolated_cell_value(local_grid, j, i, interpolate_in_dimension, field, column_stencil_values)
    type(local_grid_type), intent(inout) :: local_grid
    integer, intent(in) :: j,i
    logical, dimension(3) :: interpolate_in_dimension
    type(prognostic_field_ptr_type), intent(inout) :: field
    real(kind=DEFAULT_PRECISION), dimension(:), intent(inout) :: column_stencil_values

    integer :: k, y_interpol, x_interpol

    y_interpol=merge(j+1, j, interpolate_in_dimension(2))
    x_interpol=merge(i+1, i, interpolate_in_dimension(3))

    do k=1, local_grid%size(Z_INDEX)
      if (y_interpol .lt. 1 .or. y_interpol .gt. local_grid%local_domain_end_index(Y_INDEX)+ local_grid%halo_size(Y_INDEX) &
           .or. x_interpol .lt. 1 .or. x_interpol .gt. local_grid%local_domain_end_index(X_INDEX)+ local_grid%halo_size(X_INDEX) &
           .or. (interpolate_in_dimension(1) .and. k == local_grid%size(Z_INDEX))) then
        column_stencil_values(k)=field%ptr%data(k,j,i)
      else
        column_stencil_values(k)=0.5_DEFAULT_PRECISION*(&
             field%ptr%data(k,j,i) + field%ptr%data(merge(k+1, k, interpolate_in_dimension(1)), y_interpol, x_interpol))
      end if
    end do
  end subroutine calculate_interpolated_cell_value

  !> Copies stencil number of interpolated values from the store used by this functionality to preallocated
  !! prognostics. THe central point is always :,1,1 in the prognostic, with negative and in the x and y
  !! @param interpolated_fields The interpolated field prognostics that will be updated
  !! @param stencil The stencil configuration
  !! @param local_grid Local grid
  !! @param x Local x index of this column
  !! @param y Local y index of this column
  !! @param interpolation_id The ID of this interpolation which is used to retrieve information from the stencil configuration
  subroutine copy_interpolated_to_prognostic(interpolated_fields, stencil, local_grid, x, y, interpolation_id)
    type(prognostic_field_type), dimension(:), allocatable, intent(inout) :: interpolated_fields
    type(local_grid_type), intent(inout) :: local_grid
    type(grid_stencil_type), intent(inout) :: stencil
    integer, intent(in) :: x, y, interpolation_id

    integer :: i, y_min_src, y_max_src, x_max_src, y_min_tgt, y_max_tgt, x_min_tgt, x_max_tgt

    do i=1,stencil%nfields
      y_min_src=y-stencil%sizes(i,2)
      y_min_tgt=1-stencil%sizes(i,2)
      y_max_src=y+stencil%sizes(i,2)
      y_max_tgt=1+stencil%sizes(i,2)
      x_min_tgt=1-stencil%sizes(i,1)
      x_max_src=stencil%sizes(i,1)*2+1
      x_max_tgt=1+stencil%sizes(i,1)
      if (y_min_src .lt. 1) then 
        y_min_tgt=y_min_tgt+(1-y_min_src)
        y_min_src=1
      end if      

      interpolated_fields(i)%data(:,y_min_tgt:y_max_tgt, x_min_tgt:x_max_tgt)=&
           stencil%interpolated_fields(interpolation_id, i)%data(:, y_min_src:y_max_src, 1:x_max_src)
    end do
  end subroutine copy_interpolated_to_prognostic
end module stencil_mod
