!> Does the field stepping
!! Stepping is called at the end of processing a column and steps the x-2 column
module stepfields_mod
  use monc_component_mod, only : component_descriptor_type
  use state_mod, only : FORWARD_STEPPING, CENTRED_STEPPING, model_state_type
  use prognostics_mod, only : prognostic_field_type
  use grids_mod, only : local_grid_type, X_INDEX, Y_INDEX, Z_INDEX
  use datadefn_mod, only : DEFAULT_PRECISION
  use registry_mod, only : is_component_enabled
  use science_constants_mod, only : rlargep
  implicit none

#ifndef TEST_MODE
  private
#endif

  logical :: determine_flow_minmax=.false., cfl_is_enabled

 
  real(kind=DEFAULT_PRECISION), allocatable :: resetq_min(:)
  logical :: l_nonconservative_positive_q=.true.

  public stepfields_get_descriptor

contains

  !> Provides the descriptor back to the caller and is used in component registration
  !! @returns The KidReader component descriptor
  type(component_descriptor_type) function stepfields_get_descriptor()
    stepfields_get_descriptor%name="stepfields"
    stepfields_get_descriptor%version=0.1
    stepfields_get_descriptor%initialisation=>initialisation_callback
    stepfields_get_descriptor%timestep=>timestep_callback
    stepfields_get_descriptor%finalisation=>finalisation_callback
  end function stepfields_get_descriptor

  !> Initialisation callback
  !! @param current_state The current model state
  subroutine initialisation_callback(current_state)
    type(model_state_type), target, intent(inout) :: current_state

    allocate(resetq_min(current_state%number_q_fields))
    cfl_is_enabled=is_component_enabled(current_state%options_database, "cfltest") 
    if (cfl_is_enabled) call reset_local_minmax_values(current_state)
  end subroutine initialisation_callback


  !> Finalisation callback
  !! @param current_state The current model state
  subroutine finalisation_callback(current_state)
    type(model_state_type), target, intent(inout) :: current_state

    deallocate(resetq_min)
  end subroutine finalisation_callback

  !> Called at each timestep and will perform swapping and smoothing as required
  !! @param current_state The current model state_mod
  subroutine timestep_callback(current_state)
    type(model_state_type), target, intent(inout) :: current_state

    integer :: iq

    if (cfl_is_enabled .and. current_state%first_timestep_column) then
      if (mod(current_state%timestep, current_state%cfl_frequency) == 1 .or. &
         current_state%timestep-current_state%start_timestep .le. current_state%cfl_frequency) then
        determine_flow_minmax=.true.
        call reset_local_minmax_values(current_state)
      else
        determine_flow_minmax=.false.
      end if
    end if

    if (.not. current_state%halo_column) then
      if (determine_flow_minmax .and. cfl_is_enabled) &
         call determine_local_flow_minmax(current_state, current_state%column_local_y,  current_state%column_local_x)
      call step_all_fields(current_state)
    end if

    ! Remove negative rounding errors
    if (l_nonconservative_positive_q)then
      do iq=1,current_state%number_q_fields
        if (current_state%first_timestep_column)then
          resetq_min(iq)=minval(current_state%zq(iq)%data(:,current_state%column_local_y,  current_state%column_local_x))
        else
          resetq_min(iq)=min(resetq_min(iq),&
             minval(current_state%zq(iq)%data(:,current_state%column_local_y,  current_state%column_local_x)))
        end if
        where(current_state%zq(iq)%data(:,current_state%column_local_y,  current_state%column_local_x) < 0.0)
          current_state%zq(iq)%data(:,current_state%column_local_y,  current_state%column_local_x) = 0.0
        end where
      end do
    end if

  end subroutine timestep_callback

  !> Steps all fields
  !! @param current_state The current model state_mod
  subroutine step_all_fields(current_state)
    type(model_state_type), target, intent(inout) :: current_state

    integer :: x_prev, y_prev, i

    x_prev = current_state%column_local_x-2
    y_prev = current_state%column_local_y-1

#ifdef U_ACTIVE   
    call step_single_field(current_state%column_local_x,  current_state%column_local_y, &
         x_prev, y_prev, current_state%u, current_state%zu, current_state%su, current_state%local_grid, .true., &
         current_state%field_stepping, current_state%dtm, current_state%ugal, current_state%savu)
#endif
#ifdef V_ACTIVE
    call step_single_field(current_state%column_local_x,  current_state%column_local_y, &
         x_prev, y_prev, current_state%v, current_state%zv, current_state%sv, current_state%local_grid, .true., &
         current_state%field_stepping, current_state%dtm, current_state%vgal, current_state%savv)
#endif
#ifdef W_ACTIVE
    call step_single_field(current_state%column_local_x,  current_state%column_local_y, &
         x_prev, y_prev, current_state%w, current_state%zw, current_state%sw, current_state%local_grid, .false., &
         current_state%field_stepping, current_state%dtm, real(0., kind=DEFAULT_PRECISION), current_state%savw)
#endif
    if (current_state%th%active) call step_single_field(current_state%column_local_x,  current_state%column_local_y, &
         x_prev, y_prev, current_state%th, current_state%zth, current_state%sth, current_state%local_grid, .false., &
         current_state%field_stepping, current_state%dtm, real(0., kind=DEFAULT_PRECISION))
    do i=1,current_state%number_q_fields
      if (current_state%q(i)%active) then
        call step_single_field(current_state%column_local_x,  current_state%column_local_y, x_prev, y_prev, &
             current_state%q(i), current_state%zq(i), current_state%sq(i), current_state%local_grid, .false., &
             current_state%field_stepping, current_state%dtm, real(0., kind=DEFAULT_PRECISION))
      end if
    end do
  end subroutine step_all_fields

  !> Determines the minimum and maximum values for the local flow field. These are before the stepping, and are all reduced
  !! later on in the cfl test
  !! @param current_state The current model state
  !! @param local_y The local y index
  !! @param local_x The local x index
  subroutine determine_local_flow_minmax(current_state, local_y, local_x)
    type(model_state_type), target, intent(inout) :: current_state
    integer, intent(in) :: local_x, local_y
    
    integer :: k
    
    do k=2, current_state%local_grid%local_domain_end_index(Z_INDEX)
#ifdef U_ACTIVE
      current_state%local_zumax = max(current_state%local_zumax, current_state%zu%data(k,local_y,local_x))
      current_state%local_zumin = min(current_state%local_zumin, current_state%zu%data(k,local_y,local_x))
#endif
#ifdef V_ACTIVE
      current_state%local_zvmax = max(current_state%local_zvmax, current_state%zv%data(k,local_y,local_x))
      current_state%local_zvmin = min(current_state%local_zvmin, current_state%zv%data(k,local_y,local_x))
#endif
#ifdef W_ACTIVE
      if (k .lt. current_state%local_grid%local_domain_end_index(Z_INDEX)) then
        current_state%abswmax(k) = max(current_state%abswmax(k), abs(current_state%zw%data(k,local_y,local_x)))
      end if
#endif
    end do    
  end subroutine determine_local_flow_minmax

  !> Resets the local min and max values for the flow fields
  !! @param current_state The current model state
  subroutine reset_local_minmax_values(current_state)
    type(model_state_type), intent(inout), target :: current_state

    ! Reset the local values for the next timestep
    current_state%local_zumin=rlargep
    current_state%local_zumax=-rlargep
    current_state%local_zvmin=rlargep
    current_state%local_zvmax=-rlargep
    current_state%abswmax=-rlargep
  end subroutine reset_local_minmax_values  

  !> Steps a single specific field. This will step on the yth column of the x-2 slice and x-1 and x if this is the last slice
  !! @param x_local_index The current local x index
  !! @param y_local_index The current local y index
  !! @param x_prev The previous x index to step
  !! @param y_prev The previous y index to step
  !! @param field The prognostic field
  !! @param zfield Z prognostic field
  !! @param sfield Source terms for the prognostic field
  !! @param local_grid Description of the local grid
  !! @param flow_field Whether or not this is a flow field
  !! @param direction The stepping direction (centred or forward)
  !! @param dtm The delta time per timestep
  !! @param gal Galilean transformation
  !! @param sav Optional sav field
  subroutine step_single_field(x_local_index, y_local_index, x_prev, y_prev, field, zfield, sfield, local_grid,&
       flow_field, direction, dtm, gal, sav)
    integer, intent(in) :: x_local_index, y_local_index, x_prev, y_prev, direction
    real(kind=DEFAULT_PRECISION), intent(in) :: dtm, gal
    logical, intent(in) :: flow_field
    type(local_grid_type), intent(inout) :: local_grid
    type(prognostic_field_type), intent(inout) :: field, zfield, sfield
    type(prognostic_field_type), optional, intent(inout) :: sav

    if (x_prev .ge. local_grid%local_domain_start_index(X_INDEX)) then
      if (present(sav)) then
        call step_column_in_slice(y_local_index, x_prev, y_prev, field, zfield, sfield, local_grid, &
             flow_field, direction, dtm, gal, sav)
      else
        call step_column_in_slice(y_local_index, x_prev, y_prev, field, zfield, sfield, local_grid, &
             flow_field, direction, dtm, gal)
      end if
    end if

    if (x_local_index == local_grid%local_domain_end_index(X_INDEX)) then
      ! If this is the last slice then process x-1 (if applicable) and x too
      if (x_local_index .gt. 1) then
        if (present(sav)) then
          call step_column_in_slice(y_local_index, x_local_index-1, y_prev, field, zfield, sfield, local_grid, &
               flow_field, direction, dtm, gal, sav)
        else
          call step_column_in_slice(y_local_index, x_local_index-1, y_prev, field, zfield, sfield, local_grid, &
               flow_field, direction, dtm, gal)
        end if
      end if
      if (present(sav)) then
        call step_column_in_slice(y_local_index, x_local_index, y_prev, field, zfield, sfield, local_grid, &
             flow_field, direction, dtm, gal, sav)     
      else
        call step_column_in_slice(y_local_index, x_local_index, y_prev, field, zfield, sfield, local_grid, &
             flow_field, direction, dtm, gal)     
      end if
    end if
  end subroutine step_single_field

  !> Will step a column in a specific slice. If y_prev is large enough then will step the y-1 column and if this
  !! is the last column of the slice then will also step the current column
  !! @param x_local_index The current local x index
  !! @param y_local_index The current local y index
  !! @param x_prev The previous x index to step
  !! @param y_prev The previous y index to step
  !! @param field The prognostic field
  !! @param zfield Z prognostic field
  !! @param sfield Source terms for the prognostic field
  !! @param local_grid Description of the local grid
  !! @param flow_field Whether or not this is a flow field
  !! @param direction The stepping direction (centred or forward)
  !! @param dtm The delta time per timestep
  !! @param gal Galilean transformation
  !! @param sav Optional sav field
  subroutine step_column_in_slice(y_local_index, x_prev, y_prev, field, zfield, sfield, local_grid,&
       flow_field, direction, dtm, gal, sav)
    integer, intent(in) :: y_local_index, x_prev, y_prev, direction
    real(kind=DEFAULT_PRECISION), intent(in) :: dtm, gal
    logical, intent(in) :: flow_field
    type(local_grid_type), intent(inout) :: local_grid
    type(prognostic_field_type), intent(inout) :: field, zfield, sfield
    type(prognostic_field_type), optional, intent(inout) :: sav

    if (y_prev .ge. local_grid%local_domain_start_index(Y_INDEX)) then
      if (present(sav)) then
        call step_field(x_prev, y_prev, field, zfield, sfield, local_grid, flow_field, direction, dtm, gal, sav)
      else
        call step_field(x_prev, y_prev, field, zfield, sfield, local_grid, flow_field, direction, dtm, gal)
      end if
    end if

    if (y_local_index == local_grid%local_domain_end_index(Y_INDEX)) then
      if (present(sav)) then
        call step_field(x_prev, y_local_index, field, zfield, sfield, local_grid, flow_field, direction, dtm, gal, sav)
      else
        call step_field(x_prev, y_local_index, field, zfield, sfield, local_grid, flow_field, direction, dtm, gal)
      end if
    end if
  end subroutine step_column_in_slice

  !> Will do the actual field stepping
  !! @param flow_field Whether or not we are stepping a flow field
  !! @param direction 1=forward, 0=centred
  !! @param x_index The local X slice index
  !! @param y_index The local Y column index
  !! @param kkp Points in the vertical column
  !! @param dtm The model timestep
  !! @param field The prognostic field
  !! @param zfield The prognostic z field (which goes to timestep t+1)
  !! @param xfield The tendency of the field 
  !! @param gal The galilean transformation
  subroutine step_field(x_local_index, y_local_index, field, zfield, sfield, local_grid, flow_field, direction, dtm, gal, sav)
    integer, intent(in) :: x_local_index, y_local_index, direction
    real(kind=DEFAULT_PRECISION), intent(in) :: dtm, gal
    logical, intent(in) :: flow_field
    type(local_grid_type), intent(inout) :: local_grid
    type(prognostic_field_type), intent(inout) :: field, zfield, sfield
    type(prognostic_field_type), optional, intent(inout) :: sav

    integer :: k
    real(kind=DEFAULT_PRECISION) :: actual_gal, dtm_x2

    dtm_x2 = 2.0_DEFAULT_PRECISION * dtm

    actual_gal = merge(gal, real(0.0_DEFAULT_PRECISION, kind=DEFAULT_PRECISION), flow_field)

    sfield%data(1,y_local_index, x_local_index)=0.0_DEFAULT_PRECISION

    do k=1,local_grid%size(Z_INDEX)
      ! Save the Z field which is used in the Robert filter
      if (present(sav) .and. direction .eq. CENTRED_STEPPING) &
           sav%data(k,y_local_index, x_local_index) = zfield%data(k, y_local_index, x_local_index) + actual_gal
      if (flow_field) field%data(k, y_local_index, x_local_index) = actual_gal + field%data(k, y_local_index, x_local_index)
      if (direction == FORWARD_STEPPING) then
        zfield%data(k, y_local_index, x_local_index) = field%data(k, y_local_index, x_local_index) + dtm * &
             sfield%data(k, y_local_index, x_local_index)
      else
        zfield%data(k, y_local_index, x_local_index) = actual_gal+zfield%data(k, y_local_index, x_local_index)+dtm_x2*&
             sfield%data(k, y_local_index, x_local_index)
      end if
    end do
  end subroutine step_field
end module stepfields_mod
