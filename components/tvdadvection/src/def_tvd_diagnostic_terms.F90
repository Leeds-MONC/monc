module def_tvd_diagnostic_terms

  use state_mod, only : model_state_type
  use datadefn_mod, only : DEFAULT_PRECISION
  use grids_mod, only : local_grid_type, global_grid_type, X_INDEX, Y_INDEX, Z_INDEX 
  
  implicit none

  type str_tvd_diagnostic_terms
    
     ! TVD advective fluxes of u,v, w, q through the bottom face of grid
     ! Needed in profile_diagnostics to calc the turbulent flux diagnostic 
     ! when TVD advection is used
     real(kind=DEFAULT_PRECISION), dimension(:,:,:), allocatable :: adv_u_dgs 
     real(kind=DEFAULT_PRECISION), dimension(:,:,:), allocatable :: adv_v_dgs
     real(kind=DEFAULT_PRECISION), dimension(:,:,:), allocatable :: adv_w_dgs
     real(kind=DEFAULT_PRECISION), dimension(:,:,:), allocatable :: adv_th_dgs
     
     real(kind=DEFAULT_PRECISION), dimension(:,:,:,:), allocatable :: adv_q_dgs

  end type str_tvd_diagnostic_terms

  type(str_tvd_diagnostic_terms) :: tvd_dgs_terms

contains
  
  subroutine allocate_tvd_diagnostic_terms(current_state, tvd_dgs_terms)
    
    implicit none
    
    type(model_state_type), target, intent(in) :: current_state
    type(str_tvd_diagnostic_terms), intent(inout) :: tvd_dgs_terms

    integer :: k_top, x_local, y_local
    
    k_top = current_state%local_grid%size(Z_INDEX) + current_state%local_grid%halo_size(Z_INDEX) * 2
    x_local = current_state%local_grid%size(X_INDEX) + current_state%local_grid%halo_size(X_INDEX) * 2
    y_local = current_state%local_grid%size(Y_INDEX) + current_state%local_grid%halo_size(X_INDEX) * 2
    
    allocate(tvd_dgs_terms%adv_u_dgs(k_top, y_local, x_local), &
         tvd_dgs_terms%adv_v_dgs(k_top, y_local, x_local),     &
         tvd_dgs_terms%adv_w_dgs(k_top, y_local, x_local),     &
         tvd_dgs_terms%adv_th_dgs(k_top, y_local, x_local))
    
    tvd_dgs_terms%adv_u_dgs(:,:,:) = 0.0
    tvd_dgs_terms%adv_v_dgs(:,:,:) = 0.0
    tvd_dgs_terms%adv_w_dgs(:,:,:) = 0.0
    tvd_dgs_terms%adv_th_dgs(:,:,:) = 0.0

    if (current_state%number_q_fields > 0) then  
       allocate(tvd_dgs_terms%adv_q_dgs(k_top, y_local, x_local, current_state%number_q_fields))
       tvd_dgs_terms%adv_q_dgs(:,:,:,:)= 0.0
    endif

  end subroutine allocate_tvd_diagnostic_terms

  subroutine deallocate_tvd_diagnostic_terms(current_state, tvd_dgs_terms)
    type(model_state_type), target, intent(in) :: current_state
    type(str_tvd_diagnostic_terms), intent(inout) :: tvd_dgs_terms

    deallocate(tvd_dgs_terms%adv_u_dgs, &
         tvd_dgs_terms%adv_v_dgs,      &
         tvd_dgs_terms%adv_w_dgs,      &
         tvd_dgs_terms%adv_th_dgs)
    
    if (current_state%number_q_fields > 0) then  
       deallocate(tvd_dgs_terms%adv_q_dgs) 
    endif

  end subroutine deallocate_tvd_diagnostic_terms
    
    
end module def_tvd_diagnostic_terms
