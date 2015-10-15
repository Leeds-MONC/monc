!> Diffusion on the TH and Q fields
module diffusion_mod
  use datadefn_mod, only : DEFAULT_PRECISION
  use monc_component_mod, only : COMPONENT_ARRAY_FIELD_TYPE, COMPONENT_DOUBLE_DATA_TYPE, component_descriptor_type, &
       component_field_value_type, component_field_information_type
  use state_mod, only : FORWARD_STEPPING, model_state_type
  use prognostics_mod, only : prognostic_field_type
  use grids_mod, only : Z_INDEX, X_INDEX, Y_INDEX
  use communication_types_mod, only : halo_communication_type, neighbour_description_type, field_data_wrapper_type
  use halo_communication_mod, only : copy_buffer_to_field, perform_local_data_copy_for_field, complete_nonblocking_halo_swap, &
       copy_buffer_to_corner
  use registry_mod, only : is_component_enabled
  use logging_mod, only : LOG_ERROR, log_master_log
  implicit none

#ifndef TEST_MODE
  private
#endif

  real(kind=DEFAULT_PRECISION), dimension(:), allocatable :: th_diffusion
  real(kind=DEFAULT_PRECISION), dimension(:,:), allocatable :: q_diffusion

  public diffusion_get_descriptor

contains

  !> Provides the descriptor back to the caller and is used in component registration
  !! @returns The termination check component descriptor
  type(component_descriptor_type) function diffusion_get_descriptor()
    diffusion_get_descriptor%name="diffusion"
    diffusion_get_descriptor%version=0.1
    diffusion_get_descriptor%initialisation=>initialisation_callback
    diffusion_get_descriptor%timestep=>timestep_callback
    diffusion_get_descriptor%finalisation=>finalisation_callback

    diffusion_get_descriptor%field_value_retrieval=>field_value_retrieval_callback
    diffusion_get_descriptor%field_information_retrieval=>field_information_retrieval_callback
    allocate(diffusion_get_descriptor%published_fields(2))
    diffusion_get_descriptor%published_fields(1)="th_diffusion"
    diffusion_get_descriptor%published_fields(2)="q_diffusion"
  end function diffusion_get_descriptor

  !> Field information retrieval callback, this returns information for a specific components published field
  !! @param current_state Current model state
  !! @param name The name of the field to retrieve information for
  !! @param field_information Populated with information about the field
  subroutine field_information_retrieval_callback(current_state, name, field_information)
    type(model_state_type), target, intent(inout) :: current_state
    character(len=*), intent(in) :: name
    type(component_field_information_type), intent(out) :: field_information

    ! Field description is the same regardless of the specific field being retrieved
    field_information%field_type=COMPONENT_ARRAY_FIELD_TYPE
    field_information%data_type=COMPONENT_DOUBLE_DATA_TYPE
    if (name .eq. "q_diffusion") then
      field_information%number_dimensions=2
    else
      field_information%number_dimensions=1
    end if
    field_information%dimension_sizes(1)=current_state%local_grid%size(Z_INDEX)
    if (name .eq. "q_diffusion") field_information%dimension_sizes(2)=current_state%number_q_fields
    field_information%enabled=.true.
  end subroutine field_information_retrieval_callback

  !> Field value retrieval callback, this returns the value of a specific published field
  !! @param current_state Current model state
  !! @param name The name of the field to retrieve the value for
  !! @param field_value Populated with the value of the field
  subroutine field_value_retrieval_callback(current_state, name, field_value)
    type(model_state_type), target, intent(inout) :: current_state
    character(len=*), intent(in) :: name
    type(component_field_value_type), intent(out) :: field_value
    
    if (name .eq. "th_diffusion") then
      allocate(field_value%real_1d_array(size(th_diffusion)), source=th_diffusion)
    else if (name .eq. "q_diffusion") then
      allocate(field_value%real_2d_array(size(q_diffusion, 1), size(q_diffusion, 2)), source=q_diffusion)
    end if
  end subroutine field_value_retrieval_callback

  !> Sets up the stencil_mod (used in interpolation) and allocates data for the flux fields
  !! @param current_state The current model state_mod
  subroutine initialisation_callback(current_state)
    type(model_state_type), target, intent(inout) :: current_state

    integer :: z_size, y_size, x_size

    if (.not. is_component_enabled(current_state%options_database, "smagorinsky")) then 
      call log_master_log(LOG_ERROR, "Diffusion requires the smagorinsky component to be enabled") 
    end if

    z_size=current_state%global_grid%size(Z_INDEX)
    allocate(th_diffusion(z_size))
    allocate(q_diffusion(z_size, current_state%number_q_fields))

  end subroutine initialisation_callback

  subroutine finalisation_callback(current_state)
    type(model_state_type), target, intent(inout) :: current_state

    if (allocated(th_diffusion)) deallocate(th_diffusion)
    if (allocated(q_diffusion)) deallocate(q_diffusion)
  end subroutine finalisation_callback

  !> At each timestep will compute the diffusion source terms for TH and Q fields per column if these fields are active
  !! @param current_state The current model state
  subroutine timestep_callback(current_state)
    type(model_state_type), target, intent(inout) :: current_state

    integer :: local_y, local_x

    if (.not. current_state%use_viscosity_and_diffusion .or. current_state%halo_column) return
    if (current_state%diffusion_halo_swap_state%swap_in_progress) then
      ! If there is a diffusion halo swap in progress then complete it
      call complete_nonblocking_halo_swap(current_state, current_state%diffusion_halo_swap_state, &
           perform_local_data_copy_for_diff, copy_halo_buffer_to_diff, copy_halo_buffer_to_diff_corners)
    end if

    local_y=current_state%column_local_y
    local_x=current_state%column_local_x

    if (current_state%th%active) call perform_th_diffusion(current_state, local_y, local_x)
    if (current_state%number_q_fields .gt. 0) call perform_q_diffusion(current_state, local_y, local_x)
  end subroutine timestep_callback

  !> Computes the diffusion source terms for each Q field
  !! @param current_state The current model state
  !! @param local_y Local Y index
  !! @param local_x Local X index
  subroutine perform_q_diffusion(current_state, local_y, local_x)
    type(model_state_type), target, intent(inout) :: current_state
    integer, intent(in) :: local_y, local_x

    integer :: n

    do n=1, current_state%number_q_fields
      call general_diffusion(current_state, current_state%zq(n), current_state%sq(n), local_y, local_x, q_diffusion(:,n))
    end do
  end subroutine perform_q_diffusion  

  !> Computes the diffusion source terms for the theta field
  !! @param current_state The current model state
  !! @param local_y Local Y index
  !! @param local_x Local X index
  subroutine perform_th_diffusion(current_state, local_y, local_x)
    type(model_state_type), target, intent(inout) :: current_state
    integer, intent(in) :: local_y, local_x

    integer :: k
    real(kind=DEFAULT_PRECISION) :: adjustment

    call general_diffusion(current_state, current_state%zth, current_state%sth, local_y, local_x, th_diffusion)

    if (current_state%use_anelastic_equations) then
      ! This code only needs to be executed if anelastic, otherwise THREF is constant and the increment calculated here is zero
      do k=2, current_state%local_grid%size(Z_INDEX)
        adjustment=(current_state%global_grid%configuration%vertical%cza(k)*&
             current_state%global_grid%configuration%vertical%dthref(k)*&
             current_state%diff_coefficient%data(k, local_y, local_x) - current_state%global_grid%configuration%vertical%czb(k)*&
             current_state%global_grid%configuration%vertical%dthref(k-1)*&
             current_state%diff_coefficient%data(k-1, local_y, local_x))
        current_state%sth%data(k, local_y, local_x)=current_state%sth%data(k, local_y, local_x)+adjustment
        th_diffusion(k)=th_diffusion(k)+adjustment
      end do
    end if
  end subroutine perform_th_diffusion

  !> General diffusion computation for any field which is provided as arguments. Works in a column
  !! @param current_state The current model state
  !! @param field The field to take values from, typically zth or zq(n)
  !! @param source_field The source target field to update, typically sth or sq(n)
  !! @param local_y Local Y index
  !! @param local_x Local X index
  subroutine general_diffusion(current_state, field, source_field, local_y, local_x, diagnostics)
    type(model_state_type), target, intent(inout) :: current_state
    type(prognostic_field_type), intent(inout) :: field, source_field
    integer, intent(in) :: local_y, local_x
    real(kind=DEFAULT_PRECISION), dimension(:), intent(inout), optional :: diagnostics

    real(kind=DEFAULT_PRECISION) :: term
    integer :: k

    do k=2, current_state%local_grid%size(Z_INDEX)
      term=current_state%global_grid%configuration%horizontal%cx2*0.25_DEFAULT_PRECISION*&
           (((current_state%diff_coefficient%data(k, local_y, local_x)+&
           current_state%diff_coefficient%data(k, local_y, local_x-1))+&
           (current_state%diff_coefficient%data(k-1, local_y, local_x)+&
           current_state%diff_coefficient%data(k-1, local_y, local_x-1)))&
           *(field%data(k, local_y, local_x-1)-field%data(k, local_y, local_x)) -&
           ((current_state%diff_coefficient%data(k, local_y, local_x+1)+&
           current_state%diff_coefficient%data(k, local_y, local_x))+&
           (current_state%diff_coefficient%data(k-1, local_y, local_x+1)+&
           current_state%diff_coefficient%data(k-1, local_y, local_x)))&
           *(field%data(k, local_y, local_x)-field%data(k, local_y, local_x+1)) )&
           +current_state%global_grid%configuration%horizontal%cy2*0.25_DEFAULT_PRECISION*(&
           ((current_state%diff_coefficient%data(k, local_y, local_x)+&
           current_state%diff_coefficient%data(k, local_y-1, local_x))+&
           (current_state%diff_coefficient%data(k-1, local_y, local_x)+&
           current_state%diff_coefficient%data(k-1, local_y-1, local_x)))&
           *(field%data(k, local_y-1, local_x)-field%data(k, local_y, local_x)) -&
           ((current_state%diff_coefficient%data(k, local_y+1, local_x)+&
           current_state%diff_coefficient%data(k, local_y, local_x))+&
           (current_state%diff_coefficient%data(k-1, local_y+1, local_x)+&
           current_state%diff_coefficient%data(k-1, local_y, local_x)))&
           *(field%data(k, local_y, local_x)-field%data(k, local_y+1, local_x)) )&
           +( current_state%global_grid%configuration%vertical%czb(k)*&
           current_state%diff_coefficient%data(k-1, local_y, local_x)*&
           (field%data(k-1, local_y, local_x)-field%data(k, local_y, local_x)))

      if (k .lt. current_state%local_grid%size(Z_INDEX)) then
        term=term - current_state%global_grid%configuration%vertical%cza(k)*&
             current_state%diff_coefficient%data(k, local_y, local_x)*&
             (field%data(k, local_y, local_x)-field%data(k+1, local_y, local_x))
      end if
      source_field%data(k, local_y, local_x)=source_field%data(k, local_y, local_x)+term
      if (present(diagnostics)) diagnostics(k)=term
    end do
  end subroutine general_diffusion

    !> Does local data copying for diffusion coefficient variable halo swap
  !! @param current_state The current model state_mod
  !! @param source_data Optional source data which is written into
  subroutine perform_local_data_copy_for_diff(current_state, halo_depth, involve_corners, source_data)
    type(model_state_type), intent(inout) :: current_state
    integer, intent(in) :: halo_depth
    logical, intent(in) :: involve_corners
    type(field_data_wrapper_type), dimension(:), intent(in), optional :: source_data

    call perform_local_data_copy_for_field(current_state%diff_coefficient%data, current_state%local_grid, &
         current_state%parallel%my_rank, halo_depth, involve_corners)
  end subroutine perform_local_data_copy_for_diff

  !> Copies the halo buffer to halo location for the diffusion coefficient field
  !! @param current_state The current model state
  !! @param neighbour_description The halo swapping description of the neighbour we are accessing the buffer of
  !! @param dim The dimension we receive for
  !! @param target_index The target index for the dimension we are receiving for
  !! @param neighbour_location The location in the local neighbour data stores of this neighbour
  !! @param current_page The current, next, halo swap page to read from (all previous have been read and copied already)
  !! @param source_data Optional source data which is written into
  subroutine copy_halo_buffer_to_diff(current_state, neighbour_description, dim, target_index, &
       neighbour_location, current_page, source_data)
    type(model_state_type), intent(inout) :: current_state
    integer, intent(in) :: dim, target_index, neighbour_location
    integer, intent(inout) :: current_page(:)
    type(neighbour_description_type), intent(inout) :: neighbour_description
    type(field_data_wrapper_type), dimension(:), intent(in), optional :: source_data

    call copy_buffer_to_field(current_state%local_grid, neighbour_description%recv_halo_buffer, &
         current_state%diff_coefficient%data, dim, target_index, current_page(neighbour_location))

    current_page(neighbour_location)=current_page(neighbour_location)+1
  end subroutine copy_halo_buffer_to_diff

  !> Copies the corner halo buffer to the diffusion coefficient field corners
  !! @param current_state The current model state
  !! @param neighbour_description The halo swapping description of the neighbour we are accessing the buffer of
  !! @param corner_loc The corner location
  !! @param x_target_index The X target index for the dimension we are receiving for
  !! @param y_target_index The Y target index for the dimension we are receiving for
  !! @param neighbour_location The location in the local neighbour data stores of this neighbour
  !! @param current_page The current, next, halo swap page to read from (all previous have been read and copied already)
  !! @param source_data Optional source data which is written into
  subroutine copy_halo_buffer_to_diff_corners(current_state, neighbour_description, corner_loc, x_target_index, &
       y_target_index, neighbour_location, current_page, source_data)
    type(model_state_type), intent(inout) :: current_state
    integer, intent(in) :: corner_loc, x_target_index, y_target_index, neighbour_location
    integer, intent(inout) :: current_page(:)
    type(neighbour_description_type), intent(inout) :: neighbour_description
    type(field_data_wrapper_type), dimension(:), intent(in), optional :: source_data

    call copy_buffer_to_corner(current_state%local_grid, neighbour_description%recv_corner_buffer, &
         current_state%diff_coefficient%data, corner_loc, x_target_index, y_target_index, current_page(neighbour_location))

    current_page(neighbour_location)=current_page(neighbour_location)+1
  end subroutine copy_halo_buffer_to_diff_corners  
end module diffusion_mod
