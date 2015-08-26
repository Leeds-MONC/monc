!> Stepping of the pressure field. Completes the time-stepping of the velocity fields
!! by adding the pressure term (dp/dx_i). In addition, ensures that l_zu and l_zv satisfy the
!! Galilean-transformed boundary condition. This does not do the flow field _p terms which are only
!! needed for diagnostics, nore does it do field halo swapping which is again only needed for diagnostics
module pstep_mod
  use monc_component_mod, only : component_descriptor_type  
  use state_mod, only : model_state_type, CENTRED_STEPPING
  use grids_mod, only : Z_INDEX
  use datadefn_mod, only : DEFAULT_PRECISION
  use registry_mod, only : is_component_enabled
  use logging_mod, only : LOG_ERROR, log_master_log
  implicit none

#ifndef TEST_MODE
  private
#endif

  public pstep_get_descriptor

contains

  !> Descriptor of this component for registration
  !! @returns The pstep component descriptor
  type(component_descriptor_type) function pstep_get_descriptor()
    pstep_get_descriptor%name="pstep"
    pstep_get_descriptor%version=0.1
    pstep_get_descriptor%initialisation=>initialisation_callback
    pstep_get_descriptor%timestep=>timestep_callback
  end function pstep_get_descriptor

  !> Initialisation callback hook which will check the diverr component is enabled (as this allocates p)
  !! @param current_state The current model state
  subroutine initialisation_callback(current_state)
    type(model_state_type), target, intent(inout) :: current_state

    if (.not. is_component_enabled(current_state%options_database, "diverr")) then
      call log_master_log(LOG_ERROR, "The pstep component requires the diverr component to be enabled")
    end if
  end subroutine initialisation_callback  

  !> Called each timestep, this will step the pressure field for the non halo columns
  !! @param current_state The current model state
  subroutine timestep_callback(current_state)
    type(model_state_type), target, intent(inout) :: current_state
    
    if (current_state%galilean_transformation) call perform_galilean_transformation(current_state, &
         current_state%column_local_y, current_state%column_local_x)
    if (.not. current_state%halo_column) call step_pressure_field(current_state)
  end subroutine timestep_callback

  !> Does the actual stepping of the pressure field
  !! @param current_state The current model state
  subroutine step_pressure_field(current_state)
    type(model_state_type), target, intent(inout) :: current_state

    integer :: k, x_index, y_index
    real(kind=DEFAULT_PRECISION) :: dtmtmp

    x_index=current_state%column_local_x
    y_index=current_state%column_local_y

    dtmtmp=merge(current_state%dtm, 0.5_DEFAULT_PRECISION*current_state%dtm, current_state%field_stepping == CENTRED_STEPPING)
    do k=2,current_state%local_grid%size(Z_INDEX)

#ifdef U_ACTIVE
      current_state%zu%data(k, y_index, x_index)= current_state%zu%data(k, y_index, x_index)+ 2.0_DEFAULT_PRECISION*&
           current_state%global_grid%configuration%horizontal%cx*dtmtmp*(current_state%p%data(k, y_index, x_index)-&
           current_state%p%data(k, y_index, x_index+1))      
#endif
#ifdef V_ACTIVE
      current_state%zv%data(k, y_index, x_index)=&
           current_state%zv%data(k, y_index, x_index)+2.0_DEFAULT_PRECISION*&
           current_state%global_grid%configuration%horizontal%cy*dtmtmp*&
           (current_state%p%data(k, y_index, x_index) - current_state%p%data(k, y_index+1, x_index))
#endif
#ifdef W_ACTIVE
      if (k .lt. current_state%local_grid%size(Z_INDEX)) then
        current_state%zw%data(k, y_index, x_index)=current_state%zw%data(k, y_index, x_index)+2.0_DEFAULT_PRECISION*&
             current_state%global_grid%configuration%vertical%rdzn(k+1)*dtmtmp*(current_state%p%data(k, y_index, x_index)-&
             current_state%p%data(k+1, y_index, x_index))
      end if
#endif
    end do
    if (current_state%use_viscosity_and_diffusion .and. current_state%use_surface_boundary_conditions) then
#ifdef U_ACTIVE
      current_state%zu%data(1, y_index, x_index)=-current_state%zu%data(2, y_index, x_index)-&
           2.0_DEFAULT_PRECISION*current_state%ugal
#endif
#ifdef V_ACTIVE
      current_state%zv%data(1, y_index, x_index)=-current_state%zv%data(2, y_index, x_index)-&
           2.0_DEFAULT_PRECISION*current_state%vgal
#endif
    else
#ifdef U_ACTIVE
      current_state%zu%data(1, y_index, x_index)=current_state%zu%data(2, y_index, x_index)
#endif
#ifdef V_ACTIVE
      current_state%zv%data(1, y_index, x_index)=current_state%zv%data(2, y_index, x_index)
#endif
    end if
  end subroutine step_pressure_field  

  !> Performs Galilean transformation of flow current and z fields.
  !! @param current_state The current model state
  !! @param y_index Local y index which we are working with
  !! @param x_index Local x index which we are working with
  subroutine perform_galilean_transformation(current_state, y_index, x_index)
    type(model_state_type), target, intent(inout) :: current_state
    integer, intent(in) :: y_index, x_index

    integer :: k

    do k=1,current_state%local_grid%size(Z_INDEX)
#ifdef U_ACTIVE
      current_state%u%data(k, y_index, x_index)= current_state%u%data(k, y_index, x_index)-current_state%ugal
      current_state%zu%data(k, y_index, x_index)= current_state%zu%data(k, y_index, x_index)-current_state%ugal
#endif
#ifdef V_ACTIVE
      current_state%v%data(k, y_index, x_index)= current_state%v%data(k, y_index, x_index)-current_state%vgal
      current_state%zv%data(k, y_index, x_index)= current_state%zv%data(k, y_index, x_index)-current_state%vgal
#endif
    end do
  end subroutine perform_galilean_transformation
end module pstep_mod
