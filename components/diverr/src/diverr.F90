!> Calculates the local divergence error
module diverr_mod
  use monc_component_mod, only : component_descriptor_type
  use grids_mod, only : Z_INDEX, Y_INDEX, X_INDEX
  use state_mod, only : FORWARD_STEPPING, model_state_type
  use prognostics_mod, only : prognostic_field_type
  use datadefn_mod, only : DEFAULT_PRECISION
  implicit none

#ifndef TEST_MODE
  private
#endif

  public diverr_get_descriptor
contains

  !> Descriptor of this component for registration
  !! @returns The diverr component descriptor
  type(component_descriptor_type) function diverr_get_descriptor()
    diverr_get_descriptor%name="diverr"
    diverr_get_descriptor%version=0.1
    diverr_get_descriptor%initialisation=>init_callback
    diverr_get_descriptor%timestep=>timestep_callback
    diverr_get_descriptor%finalisation=>finalisation_callback
  end function diverr_get_descriptor

  !> The initialisation callback will allocate memory for the P field and initialise it
  !! @param current_state The current model state
  subroutine init_callback(current_state)
    type(model_state_type), target, intent(inout) :: current_state

    if (.not. allocated(current_state%p%data)) then
      ! If we are loading from a checkpoint file then this is already present
      allocate(current_state%p%data(current_state%local_grid%size(Z_INDEX) + current_state%local_grid%halo_size(Z_INDEX) * 2, &
           current_state%local_grid%size(Y_INDEX) + current_state%local_grid%halo_size(Y_INDEX) * 2, &
           current_state%local_grid%size(X_INDEX) + current_state%local_grid%halo_size(X_INDEX) * 2))
      current_state%p%data=0.0_DEFAULT_PRECISION
      current_state%p%active=.true.
    end if
  end subroutine init_callback  

  !> Called each timestep this will initialise P and update the local divergence error for the current column
  !! @param current_state The current model state
  subroutine timestep_callback(current_state)
    type(model_state_type), target, intent(inout) :: current_state

    current_state%p%data(:, current_state%column_local_y, current_state%column_local_x) = 0.0_DEFAULT_PRECISION

    if (current_state%halo_column) return
    
    if (current_state%field_stepping == FORWARD_STEPPING) then
      call handle_forward_stepping(current_state)
    else
      call handle_centred_stepping(current_state)
    end if    
    call find_max_divergence_error(current_state)
  end subroutine timestep_callback

  !> Called at finalisation this will deallocate the P field as the model shuts down
  !! @param current_state The current model state
  subroutine finalisation_callback(current_state)
    type(model_state_type), target, intent(inout) :: current_state

    current_state%p%active=.false.
    deallocate(current_state%p%data)
  end subroutine finalisation_callback  

  !> Finds the maximum divergence error locally and writes this into the local variable
  !! @param current_state The current model state
  subroutine find_max_divergence_error(current_state)
    type(model_state_type), target, intent(inout) :: current_state

    integer :: k

    if (current_state%first_timestep_column) current_state%local_divmax=0.0_DEFAULT_PRECISION

    do k=2,current_state%local_grid%size(Z_INDEX)
      if (abs(current_state%p%data(k, current_state%column_local_y, current_state%column_local_x)) .gt. &
           current_state%local_divmax) then
        current_state%local_divmax=abs(current_state%p%data(k,current_state%column_local_y, current_state%column_local_x))
      end if
    end do
  end subroutine find_max_divergence_error  

  !> Handles the calculation of P when the model is forward stepping
  !! @param current_state The current model state
  subroutine handle_forward_stepping(current_state)
    type(model_state_type), target, intent(inout) :: current_state

    real(kind=DEFAULT_PRECISION) :: timec

    timec=1.0_DEFAULT_PRECISION/current_state%dtm
   
    call calculate_p(current_state, current_state%p, current_state%u, current_state%v, current_state%w, timec, &
         current_state%column_local_y, current_state%column_local_x)
  end subroutine handle_forward_stepping

  !> Handles the calculation of P when the model is centred stepping
  !! @param current_state The current model state
  subroutine handle_centred_stepping(current_state)
    type(model_state_type), target, intent(inout) :: current_state

    real(kind=DEFAULT_PRECISION) :: timec

    timec=1.0_DEFAULT_PRECISION/(2.0_DEFAULT_PRECISION*current_state%dtm)

    call calculate_p(current_state, current_state%p, current_state%zu, current_state%zv, current_state%zw, timec, &
         current_state%column_local_y, current_state%column_local_x)
  end subroutine handle_centred_stepping

  !> Calculates P based upon flow fields and the stepping
  !! @param current_state The current model state
  !! @param p The P field to update
  !! @param u The U flow field
  !! @param v The V flow field
  !! @param w The W flow field
  !! @param timec Based upon model dtm this depends on the stepping selected
  !! @param y_local Local Y dimension index
  !! @param x_local Local X dimension index
  subroutine calculate_p(current_state, p, u, v, w, timec, y_local, x_local)
    type(model_state_type), target, intent(inout) :: current_state
    type(prognostic_field_type), intent(inout) :: u, v, w, p
    real(kind=DEFAULT_PRECISION), intent(in) :: timec
    integer, intent(in) :: y_local, x_local
    
    integer :: k
    
    p%data(:, y_local, x_local)=0.0_DEFAULT_PRECISION
    do k=2,current_state%local_grid%size(Z_INDEX)
#ifdef U_ACTIVE
      p%data(k, y_local, x_local)= p%data(k, y_local, x_local)+ &
           current_state%global_grid%configuration%horizontal%cx*(u%data(k, y_local, x_local)-u%data(k, y_local, x_local-1))
#endif
#ifdef V_ACTIVE
      p%data(k, y_local, x_local)= p%data(k, y_local, x_local)+ &
           current_state%global_grid%configuration%horizontal%cy*(v%data(k, y_local, x_local)-v%data(k, y_local-1, x_local))
#endif
#ifdef W_ACTIVE
      p%data(k, y_local, x_local)= p%data(k, y_local, x_local)+ &
           4.0_DEFAULT_PRECISION*(current_state%global_grid%configuration%vertical%tzc2(k)*w%data(k, y_local, x_local)-&
           current_state%global_grid%configuration%vertical%tzc1(k)* w%data(k-1, y_local, x_local))
#endif
      p%data(k, y_local, x_local)=p%data(k, y_local, x_local) * timec
    end do
  end subroutine calculate_p  
end module diverr_mod
