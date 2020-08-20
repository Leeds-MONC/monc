module subgrid_profile_diagnostics_mod
  use monc_component_mod, only : COMPONENT_SCALAR_FIELD_TYPE, COMPONENT_DOUBLE_DATA_TYPE, &
       COMPONENT_ARRAY_FIELD_TYPE, COMPONENT_INTEGER_DATA_TYPE, &
       component_descriptor_type, component_field_value_type, component_field_information_type
  use grids_mod, only : vertical_grid_configuration_type, Z_INDEX, Y_INDEX, X_INDEX
  use optionsdatabase_mod, only : options_get_real, options_get_integer, &
       options_get_logical
  use state_mod, only : model_state_type
  use datadefn_mod, only : DEFAULT_PRECISION
  use q_indices_mod, only: get_q_index, standard_q_names
  use registry_mod, only : is_component_enabled
  use logging_mod, only : LOG_ERROR, LOG_WARN, log_master_log
  use smagorinsky_mod, only: calculate_half_squared_strain_rate, &
                             calculate_richardson_number,        &
                             calculate_thermal_dissipation_rate
  use science_constants_mod, only: G, ratio_mol_wts
  

  implicit none

#ifndef TEST_MODE
  private
#endif

  integer :: total_points, iqv=0, iql=0, iqr=0, iqi=0, iqs=0,    &
       iqg=0

! diagnostic variables
  real(kind=DEFAULT_PRECISION), dimension(:), allocatable ::     &
       uwsg_tot, vwsg_tot, uusg_tot, vvsg_tot, wwsg_tot,         &
       tkesg_tot, wtsg_tot, th2sg_tot, wqsg_tot,                 &
       ! subgrid tke fluxes
       sed_tot,ssub_tot, dissipation_tot,buoysg_tot, wkesg_tot,  &
       theta_dis_tot, vis_coef_tot, diff_coef_tot,               &
       richardson_number_tot, richardson_squared_tot,            &
       ! subgrid moisture fluxes
       wqv_sg_tot, wql_sg_tot, wqr_sg_tot, wqi_sg_tot,           &
       wqs_sg_tot, wqg_sg_tot
! These arrays should in due course be allocated conditionally,
! but let's just get it working first.
  real(kind=DEFAULT_PRECISION), dimension(:), allocatable ::     &
       epsth, ssq, elamr_sq, richardson_number,     &
       subgrid_tke
  real(kind=DEFAULT_PRECISION), dimension(:,:), allocatable ::   &
       subke_2d
! Constants here provisionally
  real(kind=DEFAULT_PRECISION) :: a2_n, ath2_n, pr_n, ri_crit
  real(kind=DEFAULT_PRECISION) :: qlcrit

  type(vertical_grid_configuration_type) :: vertical_grid

  logical :: l_lem_dissipation_rate = .true.
  
  public subgrid_profile_diagnostics_get_descriptor

contains

  !> Provides the component descriptor for the core to register
  !! @returns The descriptor describing this component
  type(component_descriptor_type) function subgrid_profile_diagnostics_get_descriptor()
    subgrid_profile_diagnostics_get_descriptor%name="subgrid_profile_diagnostics"
    subgrid_profile_diagnostics_get_descriptor%version=0.1

    subgrid_profile_diagnostics_get_descriptor%initialisation=>initialisation_callback
    subgrid_profile_diagnostics_get_descriptor%timestep=>timestep_callback

    subgrid_profile_diagnostics_get_descriptor%field_value_retrieval=>field_value_retrieval_callback
    subgrid_profile_diagnostics_get_descriptor%field_information_retrieval=>field_information_retrieval_callback
!   Note: Multiple copies of diagnostics are required for different time processing, so
!   duplicate the variables as often as necessary
    allocate(subgrid_profile_diagnostics_get_descriptor%published_fields(35))

    subgrid_profile_diagnostics_get_descriptor%published_fields(1)="uwsg_total_local"
    subgrid_profile_diagnostics_get_descriptor%published_fields(2)="vwsg_total_local"
    subgrid_profile_diagnostics_get_descriptor%published_fields(3)="uusg_total_local"
    subgrid_profile_diagnostics_get_descriptor%published_fields(4)="vvsg_total_local"
    subgrid_profile_diagnostics_get_descriptor%published_fields(5)="wwsg_total_local"
    subgrid_profile_diagnostics_get_descriptor%published_fields(6)="tkesg_total_local"
    subgrid_profile_diagnostics_get_descriptor%published_fields(7)="wtsg_total_local"
    subgrid_profile_diagnostics_get_descriptor%published_fields(8)="th2sg_total_local"
    subgrid_profile_diagnostics_get_descriptor%published_fields(9)="wqsg_total_local"
    !SF TKE diags
    subgrid_profile_diagnostics_get_descriptor%published_fields(10)="buoysg_total_local"
    subgrid_profile_diagnostics_get_descriptor%published_fields(11)="sed_total_local"
    subgrid_profile_diagnostics_get_descriptor%published_fields(12)="ssub_total_local"
    subgrid_profile_diagnostics_get_descriptor%published_fields(13)="dissipation_total_local"
    !===========
    subgrid_profile_diagnostics_get_descriptor%published_fields(14)="subke"
    subgrid_profile_diagnostics_get_descriptor%published_fields(15)="wqv_sg_total_local"
    subgrid_profile_diagnostics_get_descriptor%published_fields(16)="wql_sg_total_local"
    subgrid_profile_diagnostics_get_descriptor%published_fields(17)="wqr_sg_total_local"
    subgrid_profile_diagnostics_get_descriptor%published_fields(18)="wqi_sg_total_local"
    subgrid_profile_diagnostics_get_descriptor%published_fields(19)="wqs_sg_total_local"
    subgrid_profile_diagnostics_get_descriptor%published_fields(20)="wqg_sg_total_local"

!   =====================================================
!   2nd, provisionally instantaneous, stream

    subgrid_profile_diagnostics_get_descriptor%published_fields(21)="i_uwsg_total_local"
    subgrid_profile_diagnostics_get_descriptor%published_fields(22)="i_vwsg_total_local"
    subgrid_profile_diagnostics_get_descriptor%published_fields(23)="i_uusg_total_local"
    subgrid_profile_diagnostics_get_descriptor%published_fields(24)="i_vvsg_total_local"
    subgrid_profile_diagnostics_get_descriptor%published_fields(25)="i_wwsg_total_local"
    subgrid_profile_diagnostics_get_descriptor%published_fields(26)="i_tkesg_total_local"
    subgrid_profile_diagnostics_get_descriptor%published_fields(27)="i_wtsg_total_local"
    subgrid_profile_diagnostics_get_descriptor%published_fields(28)="i_th2sg_total_local"
    subgrid_profile_diagnostics_get_descriptor%published_fields(29)="i_wqsg_total_local"
!   =====================================================

    subgrid_profile_diagnostics_get_descriptor%published_fields(30)="wkesg_total_local"
    subgrid_profile_diagnostics_get_descriptor%published_fields(31)="theta_dis_total_local"
    subgrid_profile_diagnostics_get_descriptor%published_fields(32)="viscosity_coef_total_local"
    subgrid_profile_diagnostics_get_descriptor%published_fields(33)="diffusion_coef_total_local"
    subgrid_profile_diagnostics_get_descriptor%published_fields(34)="richardson_number_total_local"
    subgrid_profile_diagnostics_get_descriptor%published_fields(35)="richardson_squared_total_local"
    

  end function subgrid_profile_diagnostics_get_descriptor

  subroutine initialisation_callback(current_state)
    type(model_state_type), target, intent(inout) :: current_state

    integer :: k

    if (.not. is_component_enabled(current_state%options_database, "smagorinsky")) then
       call log_master_log(LOG_ERROR, "Subgrid model diags requested but subgrid model not enabled - check config")
    end if

    vertical_grid=current_state%global_grid%configuration%vertical

    total_points=(current_state%local_grid%size(Y_INDEX) * current_state%local_grid%size(X_INDEX))
    
    allocate(uwsg_tot(current_state%local_grid%size(Z_INDEX))   &
         ,   vwsg_tot(current_state%local_grid%size(Z_INDEX))   &
         ,   uusg_tot(current_state%local_grid%size(Z_INDEX))   &
         ,   vvsg_tot(current_state%local_grid%size(Z_INDEX))   &
         ,   wwsg_tot(current_state%local_grid%size(Z_INDEX))   &
         ,   tkesg_tot(current_state%local_grid%size(Z_INDEX))  &
         ,   wkesg_tot(current_state%local_grid%size(Z_INDEX))  &
         ,   vis_coef_tot(current_state%local_grid%size(Z_INDEX))  &
         ,   diff_coef_tot(current_state%local_grid%size(Z_INDEX)) &
         ,   richardson_number_tot(current_state%local_grid%size(Z_INDEX))  &
         ,   richardson_squared_tot(current_state%local_grid%size(Z_INDEX)) &
         ,   dissipation_tot(current_state%local_grid%size(Z_INDEX) ) &
         ,   ssub_tot(current_state%local_grid%size(Z_INDEX))  &
         ,   sed_tot(current_state%local_grid%size(Z_INDEX))  &
         ,   buoysg_tot(current_state%local_grid%size(Z_INDEX)) )

!   Check consistency of any changes with the deallocations later.
    allocate(ssq(current_state%local_grid%size(Z_INDEX)))
    allocate(elamr_sq(current_state%local_grid%size(Z_INDEX)))
    allocate(richardson_number(current_state%local_grid%size(Z_INDEX)))
    allocate(subgrid_tke(current_state%local_grid%size(Z_INDEX)))
    ! subgrid tke 2d scalar field, included in this component to prevent copying
    ! of code. Need to revisit this
    allocate(subke_2d(current_state%local_grid%size(Y_INDEX),current_state%local_grid%size(X_INDEX)))
    if (current_state%th%active) then
       allocate(epsth(current_state%local_grid%size(Z_INDEX))         &
            ,   theta_dis_tot(current_state%local_grid%size(Z_INDEX)) &
            ,   wtsg_tot(current_state%local_grid%size(Z_INDEX))      &
            ,   th2sg_tot(current_state%local_grid%size(Z_INDEX)))    
    endif
        
    if (.not. current_state%passive_q .and. current_state%number_q_fields .gt. 0) then 
       allocate(wqsg_tot(current_state%local_grid%size(Z_INDEX)))
       iqv=get_q_index(standard_q_names%VAPOUR, 'subgrid_profile_diags')
       iql=get_q_index(standard_q_names%CLOUD_LIQUID_MASS, 'subgrid_profile_diags') 
       qlcrit=options_get_real(current_state%options_database, "qlcrit")
       allocate(wqv_sg_tot(current_state%local_grid%size(Z_INDEX)),   &
            wql_sg_tot(current_state%local_grid%size(Z_INDEX)))
    endif
    if (current_state%rain_water_mixing_ratio_index > 0) then
       iqr = current_state%rain_water_mixing_ratio_index
       allocate(wqr_sg_tot(current_state%local_grid%size(Z_INDEX)))
    endif
    if (current_state%ice_water_mixing_ratio_index > 0) then
       iqi = current_state%ice_water_mixing_ratio_index
       allocate(wqi_sg_tot(current_state%local_grid%size(Z_INDEX)))
    endif
    if (current_state%snow_water_mixing_ratio_index > 0) then
       iqs = current_state%snow_water_mixing_ratio_index
       allocate(wqs_sg_tot(current_state%local_grid%size(Z_INDEX)))
    endif
    if (current_state%graupel_water_mixing_ratio_index > 0) then
       iqg = current_state%graupel_water_mixing_ratio_index
       allocate(wqg_sg_tot(current_state%local_grid%size(Z_INDEX)))
    endif
   
    l_lem_dissipation_rate=options_get_logical(current_state%options_database, "l_lem_dissipation_rate")

  end subroutine initialisation_callback  

  subroutine timestep_callback(current_state)
    type(model_state_type), target, intent(inout) :: current_state

    integer :: k, i, j
    integer :: jcol, icol, target_x_index, target_y_index
    integer :: top_index
    real(kind=DEFAULT_PRECISION), dimension(current_state%local_grid%size(Z_INDEX)) :: &
      S11, S22, S33, S12, S23, S13, S33_on_p
    real(kind=DEFAULT_PRECISION), dimension(current_state%local_grid%size(Z_INDEX)) :: &
      tau11,tau22, tau33, tau12, tau23, tau13, tau33_on_p
    real(kind=DEFAULT_PRECISION), dimension(current_state%local_grid%size(Z_INDEX)) :: u_i_prime_tau_i
    real(kind=DEFAULT_PRECISION) :: vistmp, vis12, vis12a, visonp2, visonp2a, vis13, vis23
    real(kind=DEFAULT_PRECISION), dimension(current_state%local_grid%size(Z_INDEX)) :: umean, wu_umean, vmean, wv_vmean, &
         w_pprime_at_w, rke1, w_qvprime_at_w, w_qclprime_at_w, w_thprime_at_w, wq, rho, rec_rho, buoy_cof, &
         uw_tot, vw_tot,qvprime_at_w,qclprime_at_w,thprime_at_w, sed_eq13, sed_eq23, sed33,uwsg,vwsg,sed_swap
    real(kind=DEFAULT_PRECISION) :: C_virtual
    real(kind=DEFAULT_PRECISION) :: upr_at_w,vpr_at_w, vprime_w_local, uprime_w_local, U_1, U_1_bar
    real(kind=DEFAULT_PRECISION) :: buoy_prod, sg_shear_prod, dissipation_rate
    
    logical :: use_Ri_for_buoyant_prod=.TRUE.

    if (.not. current_state%diagnostic_sample_timestep) return

    C_virtual = (ratio_mol_wts-1.0_DEFAULT_PRECISION)
    jcol=current_state%column_local_y
    icol=current_state%column_local_x
    target_y_index=jcol-current_state%local_grid%halo_size(Y_INDEX)
    target_x_index=icol-current_state%local_grid%halo_size(X_INDEX)

      if (current_state%first_timestep_column) then
       sed_tot(:)      = 0.0_DEFAULT_PRECISION 
       ssub_tot(:)     = 0.0_DEFAULT_PRECISION 
       buoysg_tot(:)   = 0.0_DEFAULT_PRECISION 
       sed_eq23(:)     = 0.0_DEFAULT_PRECISION
       sed_eq13(:)     = 0.0_DEFAULT_PRECISION
       sed_swap(:)     = 0.0_DEFAULT_PRECISION
       uwsg_tot(:)     = 0.0_DEFAULT_PRECISION
       vwsg_tot(:)     = 0.0_DEFAULT_PRECISION
       uusg_tot(:)     = 0.0_DEFAULT_PRECISION
       vvsg_tot(:)     = 0.0_DEFAULT_PRECISION
       wwsg_tot(:)     = 0.0_DEFAULT_PRECISION
       tkesg_tot(:)    = 0.0_DEFAULT_PRECISION
       vis_coef_tot(:) = 0.0_DEFAULT_PRECISION
       diff_coef_tot(:)= 0.0_DEFAULT_PRECISION
       dissipation_tot(:)      = 0.0_DEFAULT_PRECISION
       wkesg_tot(:)    = 0.0_DEFAULT_PRECISION
       richardson_number_tot(:) = 0.0_DEFAULT_PRECISION
       richardson_squared_tot(:) = 0.0_DEFAULT_PRECISION
       subke_2d(:,:)   = 0.0_DEFAULT_PRECISION
       if (current_state%th%active) then
         wtsg_tot(:)     = 0.0_DEFAULT_PRECISION
         th2sg_tot(:)    = 0.0_DEFAULT_PRECISION
         theta_dis_tot(:)= 0.0_DEFAULT_PRECISION
       endif
       if (.not. current_state%passive_q .and. &
         current_state%number_q_fields .gt. 0) then
         wqsg_tot(:)     = 0.0_DEFAULT_PRECISION
         wqv_sg_tot(:)=0.0_DEFAULT_PRECISION
         wql_sg_tot(:)=0.0_DEFAULT_PRECISION
         if (iqr > 0) wqr_sg_tot(:) = 0.0_DEFAULT_PRECISION
         if (iqi > 0) wqi_sg_tot(:) = 0.0_DEFAULT_PRECISION
         if (iqs > 0) wqs_sg_tot(:) = 0.0_DEFAULT_PRECISION
         if (iqg > 0) wqg_sg_tot(:) = 0.0_DEFAULT_PRECISION
       endif
    end if

      if (.not. current_state%halo_column) then
    
       do k=2,current_state%local_grid%size(Z_INDEX)-1

! Consistent with ssq calculation from calculate_half_squared_strain_rate
! It could usefully be moved into the Smagorinsky component as a subroutine.
       
#ifdef U_ACTIVE    
      ! Average over 2 points W but 0.5 cancels 2
         S11(k)=current_state%global_grid%configuration%horizontal%cx*(&
           (current_state%u%data(k+1, jcol, icol)-&
           current_state%u%data(k+1, jcol, icol-1))+&
           (current_state%u%data(k, jcol, icol)-&
           current_state%u%data(k, jcol, icol-1)))
#else
         S11(k)=0.0_DEFAULT_PRECISION
#endif
#ifdef V_ACTIVE
      ! Average over 2 points W but 0.5 cancels 2
         S22(k)=current_state%global_grid%configuration%horizontal%cy*(&
           (current_state%v%data(k+1, jcol, icol)-&
           current_state%v%data(k+1, jcol-1, icol))+&
           (current_state%v%data(k, jcol, icol)-&
           current_state%v%data(k, jcol-1, icol)))
#else
         S22(k)=0.0_DEFAULT_PRECISION
#endif
#ifdef W_ACTIVE
      ! Average over 2 points W but 0.5 cancels 2
         S33(k)=((current_state%w%data(k, jcol, icol)-&
           current_state%w%data(k-1, jcol, icol))*&
           current_state%global_grid%configuration%vertical%rdz(k)) +&
           ((current_state%w%data(k+1, jcol, icol)-&
           current_state%w%data(k, jcol, icol))*&
           current_state%global_grid%configuration%vertical%rdz(k+1))
           
      ! No average needed
         S33_on_p(k) = 2.0_DEFAULT_PRECISION * &
         ((current_state%w%data(k,   jcol, icol)-&
           current_state%w%data(k-1, jcol, icol))*&
           current_state%global_grid%configuration%vertical%rdz(k))
#else
         S33(k)=0.0_DEFAULT_PRECISION
         S33_on_p(k) = 0.0_DEFAULT_PRECISION
#endif
#if defined(U_ACTIVE) && defined(W_ACTIVE)
      ! Average over 2 points U and W
         S13(k)=(((current_state%u%data(k+1, jcol, icol)-&
                   current_state%u%data(k, jcol, icol))*&
                  current_state%global_grid%configuration%vertical%rdzn(k+1)+&
                  (current_state%w%data(k, jcol, icol+1)-&
                   current_state%w%data(k, jcol, icol))*&
                  current_state%global_grid%configuration%horizontal%cx)+&
                 ((current_state%u%data(k+1, jcol, icol-1)-&
                   current_state%u%data(k, jcol, icol-1))*&
                  current_state%global_grid%configuration%vertical%rdzn(k+1)+&
                  (current_state%w%data(k, jcol, icol)-&
                   current_state%w%data(k, jcol, icol-1))*&
                  current_state%global_grid%configuration%horizontal%cx))*0.5_DEFAULT_PRECISION
#elif defined(U_ACTIVE) && !defined(W_ACTIVE)
         S13(k)=(((current_state%u%data(k+1, jcol, icol)-&
           current_state%u%data(k, jcol, icol))*&
           current_state%global_grid%configuration%vertical%rdzn(k+1))+&
           ((current_state%u%data(k+1, jcol, icol-1)-&
           current_state%u%data(k, jcol, icol-1))*&
           current_state%global_grid%configuration%vertical%rdzn(k+1)))*0.5_DEFAULT_PRECISION
#elif !defined(U_ACTIVE) && defined(W_ACTIVE)
         S13(k)=(((current_state%w%data(k, jcol, icol+1)-&
           current_state%w%data(k, jcol, icol))*&
           current_state%global_grid%configuration%horizontal%cx)+&
           ((current_state%w%data(k, jcol, icol)-&
           current_state%w%data(k, jcol, icol-1))*&
           current_state%global_grid%configuration%horizontal%cx))*0.5_DEFAULT_PRECISION
#else
         S13(k)=0.0_DEFAULT_PRECISION
#endif
#if defined(W_ACTIVE) && defined(V_ACTIVE)
      ! Average over 2 points W and V
         S23(k)=(((current_state%w%data(k, jcol, icol) - &
                   current_state%w%data(k, jcol-1, icol)) * &
                  current_state%global_grid%configuration%horizontal%cy + &
                  (current_state%v%data(k+1, jcol-1, icol) - &
                   current_state%v%data(k, jcol-1, icol)) * &
                  current_state%global_grid%configuration%vertical%rdzn(k+1)) + &
                 ((current_state%w%data(k, jcol+1, icol) - &
                   current_state%w%data(k, jcol, icol)) * &
                  current_state%global_grid%configuration%horizontal%cy + &
                  (current_state%v%data(k+1, jcol, icol)-&
                   current_state%v%data(k, jcol, icol))*&
                  current_state%global_grid%configuration%vertical%rdzn(k+1)))*0.5_DEFAULT_PRECISION
#elif defined(W_ACTIVE) && !defined(V_ACTIVE)
         S23(k)=(((current_state%w%data(k, jcol, icol)-&
           current_state%w%data(k, jcol-1, icol))*&
           current_state%global_grid%configuration%horizontal%cy)+&
           ((current_state%w%data(k, jcol+1, icol)-&
           current_state%w%data(k, jcol, icol))*&
           current_state%global_grid%configuration%horizontal%cy))*0.5_DEFAULT_PRECISION
#elif !defined(W_ACTIVE) && defined(V_ACTIVE)
         S23(k)=(((current_state%v%data(k+1, jcol-1, icol)-&
           current_state%v%data(k, jcol-1, icol))*&
           current_state%global_grid%configuration%vertical%rdzn(k+1))+&
           ((current_state%v%data(k+1, jcol, icol)-&
           current_state%v%data(k, jcol, icol))*&
           current_state%global_grid%configuration%vertical%rdzn(k+1)))*0.5_DEFAULT_PRECISION
#else
         S23(k)=0.0_DEFAULT_PRECISION
#endif
#if defined(U_ACTIVE) && defined(V_ACTIVE)
      ! Average over 8 points from U and V
         S12(k)=(((((current_state%u%data(k, jcol, icol-1)-&
           current_state%u%data(k, jcol-1, icol-1))*&
           current_state%global_grid%configuration%horizontal%cy+&
           (current_state%v%data(k, jcol-1, icol)-&
           current_state%v%data(k, jcol-1, icol-1))*&
           current_state%global_grid%configuration%horizontal%cx) +&
           ((current_state%u%data(k, jcol+1, icol-1)-&
           current_state%u%data(k, jcol, icol-1))*&
           current_state%global_grid%configuration%horizontal%cy+&
           (current_state%v%data(k, jcol, icol)-&
           current_state%v%data(k, jcol, icol-1))*&
           current_state%global_grid%configuration%horizontal%cx)) +(&
           ((current_state%u%data(k, jcol, icol)-&
           current_state%u%data(k, jcol-1, icol))*&
           current_state%global_grid%configuration%horizontal%cy+&
           (current_state%v%data(k, jcol-1, icol+1)-&
           current_state%v%data(k, jcol-1, icol))*&
           current_state%global_grid%configuration%horizontal%cx) +&
           ((current_state%u%data(k, jcol+1, icol)-&
           current_state%u%data(k, jcol, icol))*&
           current_state%global_grid%configuration%horizontal%cy+&
           (current_state%v%data(k, jcol, icol+1)-&
           current_state%v%data(k, jcol, icol))*&
           current_state%global_grid%configuration%horizontal%cx)))+((&
           ((current_state%u%data(k+1, jcol, icol-1)-&
           current_state%u%data(k+1, jcol-1, icol-1))*&
           current_state%global_grid%configuration%horizontal%cy+&
           (current_state%v%data(k+1, jcol-1, icol)-&
           current_state%v%data(k+1, jcol-1, icol-1))*&
           current_state%global_grid%configuration%horizontal%cx)+&
           ((current_state%u%data(k+1, jcol+1, icol-1)-&
           current_state%u%data(k+1, jcol, icol-1))*&
           current_state%global_grid%configuration%horizontal%cy+&
           (current_state%v%data(k+1, jcol, icol)-&
           current_state%v%data(k+1, jcol, icol-1))*&
           current_state%global_grid%configuration%horizontal%cx))+(&
           ((current_state%u%data(k+1, jcol, icol)-&
           current_state%u%data(k+1, jcol-1, icol))*&
           current_state%global_grid%configuration%horizontal%cy+&
           (current_state%v%data(k+1, jcol-1, icol+1)-&
           current_state%v%data(k+1, jcol-1, icol))*&
           current_state%global_grid%configuration%horizontal%cx)+&
           ((current_state%u%data(k+1, jcol+1, icol)-&
           current_state%u%data(k+1, jcol, icol))*&
           current_state%global_grid%configuration%horizontal%cy+&
           (current_state%v%data(k+1, jcol, icol+1)-&
           current_state%v%data(k+1, jcol, icol))*&
           current_state%global_grid%configuration%horizontal%cx))))*0.125_DEFAULT_PRECISION

#elif defined(U_ACTIVE) && !defined(V_ACTIVE)

         S12(k)=(((((current_state%u%data(k, jcol, icol-1)-&
           current_state%u%data(k, jcol-1, icol-1))*&
           current_state%global_grid%configuration%horizontal%cy) +&
           ((current_state%u%data(k, jcol+1, icol-1)-&
           current_state%u%data(k, jcol, icol-1))*&
           current_state%global_grid%configuration%horizontal%cy)) +(&
           ((current_state%u%data(k, jcol, icol)-&
           current_state%u%data(k, jcol-1, icol))*&
           current_state%global_grid%configuration%horizontal%cy) +&
           ((current_state%u%data(k, jcol+1, icol)-&
           current_state%u%data(k, jcol, icol))*&
           current_state%global_grid%configuration%horizontal%cy)))+((&
           ((current_state%u%data(k+1, jcol, icol-1)-&
           current_state%u%data(k+1, jcol-1, icol-1))*&
           current_state%global_grid%configuration%horizontal%cy)+&
           ((current_state%u%data(k+1, jcol+1, icol-1)-&
           current_state%u%data(k+1, jcol, icol-1))*&
           current_state%global_grid%configuration%horizontal%cy))+(&
           ((current_state%u%data(k+1, jcol, icol)-&
           current_state%u%data(k+1, jcol-1, icol))*&
           current_state%global_grid%configuration%horizontal%cy)+&
           ((current_state%u%data(k+1, jcol+1, icol)-&
           current_state%u%data(k+1, jcol, icol))*&
           current_state%global_grid%configuration%horizontal%cy))))*0.125_DEFAULT_PRECISION

#elif !defined(U_ACTIVE) && defined(V_ACTIVE)

         S12(k)=(((((current_state%v%data(k, jcol-1, icol)-&
           current_state%v%data(k, jcol-1, icol-1))*&
           current_state%global_grid%configuration%horizontal%cx) +&
           ((current_state%v%data(k, jcol, icol)-&
           current_state%v%data(k, jcol, icol-1))*&
           current_state%global_grid%configuration%horizontal%cx)) +&
           (((current_state%v%data(k, jcol-1, icol+1)-&
           current_state%v%data(k, jcol-1, icol))*&
           current_state%global_grid%configuration%horizontal%cx) +&
           ((current_state%v%data(k, jcol, icol+1)-&
           current_state%v%data(k, jcol, icol))*&
           current_state%global_grid%configuration%horizontal%cx)))+((&
           ((current_state%v%data(k+1, jcol-1, icol)-&
           current_state%v%data(k+1, jcol-1, icol-1))*&
           current_state%global_grid%configuration%horizontal%cx)+&
           ((current_state%v%data(k+1, jcol, icol)-&
           current_state%v%data(k+1, jcol, icol-1))*&
           current_state%global_grid%configuration%horizontal%cx))+(&
           ((current_state%v%data(k+1, jcol-1, icol+1)-&
           current_state%v%data(k+1, jcol-1, icol))*&
           current_state%global_grid%configuration%horizontal%cx)+&
           ((current_state%v%data(k+1, jcol, icol+1)-&
           current_state%v%data(k+1, jcol, icol))*&
           current_state%global_grid%configuration%horizontal%cx))))*0.125_DEFAULT_PRECISION
#else
         S12(k)=0.0_DEFAULT_PRECISION
#endif
       enddo   
       
#if defined(U_ACTIVE)
      ! Average over 2 points U and W
       S13(1)=(current_state%u%data(2, jcol, icol) + &
               current_state%u%data(2, jcol, icol-1)) * &
               0.5_DEFAULT_PRECISION / current_state%global_grid%configuration%vertical%zn(2)
#else
       S13(1)=0.0_DEFAULT_PRECISION
#endif
#if defined(V_ACTIVE)
      ! Average over 2 points W and V
       S23(1)=(current_state%v%data(2, jcol, icol) + &
               current_state%v%data(2, jcol-1, icol)) * &
               0.5_DEFAULT_PRECISION / current_state%global_grid%configuration%vertical%zn(2)
#else
       S23(1)=0.0_DEFAULT_PRECISION
#endif
       S12(1)=0.0_DEFAULT_PRECISION
    
       do k=1,current_state%local_grid%size(Z_INDEX)-1
         tau11(k) = current_state%global_grid%configuration%vertical%rho(k) * &
           current_state%vis_coefficient%data(k,jcol,icol) * &
           S11(k)
         tau22(k) = current_state%global_grid%configuration%vertical%rho(k) * &
           current_state%vis_coefficient%data(k,jcol,icol) * &
           S22(k)
         tau33(k) = current_state%global_grid%configuration%vertical%rho(k) * &
           current_state%vis_coefficient%data(k,jcol,icol) * &
           S33(k)
         tau12(k) = current_state%global_grid%configuration%vertical%rho(k) * &
           current_state%vis_coefficient%data(k,jcol,icol) * &
           S12(k)
         tau13(k) = current_state%global_grid%configuration%vertical%rho(k) * &
           current_state%vis_coefficient%data(k,jcol,icol) * &
           S13(k)
         tau23(k) = current_state%global_grid%configuration%vertical%rho(k) * &
           current_state%vis_coefficient%data(k,jcol,icol) * &
           S23(k)
       enddo 
         
       do k=2,current_state%local_grid%size(Z_INDEX)
         tau33_on_p(k) = current_state%global_grid%configuration%vertical%rhon(k) * 0.5 *&
           (current_state%vis_coefficient%data(k-1,jcol,icol) + &
            current_state%vis_coefficient%data(k,  jcol,icol)) * &
           S33_on_p(k)
       enddo   

!      Subgrid diagnostics
       do k=1, current_state%local_grid%size(Z_INDEX)-1
!
!         Field on kth rho-level, so gradient evaluated from k+1st and kth rhon-levels. Horizontally aligned with
!         w-points.
          uwsg_tot(k) = uwsg_tot(k) - 0.5 * &
            current_state%vis_coefficient%data(k,jcol,icol-1) * &
            ( current_state%u%data(k+1,jcol,icol-1) - &
              current_state%u%data(k,jcol,icol-1) ) * &
            vertical_grid%rdzn(k+1)
          uwsg_tot(k) = uwsg_tot(k) - 0.5 * &
            current_state%vis_coefficient%data(k,jcol,icol) * &
            ( current_state%u%data(k+1,jcol,icol) - &
              current_state%u%data(k,jcol,icol) ) * &
            vertical_grid%rdzn(k+1)
!
!         Field on kth rho-level, so gradient evaluated from k+1st and kth rhon-levels. Horizontally aligned with
!         w-points.
          vwsg_tot(k) = vwsg_tot(k) - 0.5 * &
            current_state%vis_coefficient%data(k,jcol-1,icol) * &
            ( current_state%v%data(k+1,jcol-1,icol) - &
              current_state%v%data(k,jcol-1,icol) ) * &
            vertical_grid%rdzn(k+1)
          vwsg_tot(k) = vwsg_tot(k) - 0.5 * &
            current_state%vis_coefficient%data(k,jcol,icol) * &
            ( current_state%v%data(k+1,jcol,icol) - &
              current_state%v%data(k,jcol,icol) ) * &
            vertical_grid%rdzn(k+1)
          !
          ! viscosity and diffusion coefficients
          vis_coef_tot(k) = vis_coef_tot(k) + &
               current_state%vis_coefficient%data(k,jcol,icol)
          diff_coef_tot(k) = diff_coef_tot(k) + &
               current_state%diff_coefficient%data(k,jcol,icol)
       enddo
       !         Field on kth rho-level, so gradient evaluated from k+1st and kth rhon-levels. Horizontally aligned with
       !         w-points.
       if (current_state%th%active) then
          do k=1, current_state%local_grid%size(Z_INDEX)-1
             wtsg_tot(k) = wtsg_tot(k) -  &
                  current_state%diff_coefficient%data(k,jcol,icol) * &
                  ( current_state%th%data(k+1,jcol,icol) + &
                  vertical_grid%thref(k+1) - &
                  current_state%th%data(k,jcol,icol) - &
                  vertical_grid%thref(k) ) * &
                  vertical_grid%rdzn(k+1)
             !
          enddo
       endif

       if (current_state%th%active .and. .not. current_state%passive_q .and. current_state%number_q_fields .gt. 0) then
          
          do k=1, current_state%local_grid%size(Z_INDEX)-1
             wqv_sg_tot(k) = wqv_sg_tot(k) + &
                  current_state%diff_coefficient%data(k, jcol, icol)* &
                  vertical_grid%rdzn(k+1) * &
                  (current_state%q(iqv)%data(k,jcol,icol) - current_state%q(iqv)%data(k+1,jcol,icol))
             wql_sg_tot(k) = wql_sg_tot(k) + &
                  current_state%diff_coefficient%data(k, jcol, icol)* &
                  vertical_grid%rdzn(k+1) * &
                  (current_state%q(iql)%data(k,jcol,icol) - current_state%q(iql)%data(k+1,jcol,icol))
          enddo
          if (iqr > 0) then
             do k=1, current_state%local_grid%size(Z_INDEX)-1
                wqr_sg_tot(k) = wqr_sg_tot(k) + &
                     current_state%diff_coefficient%data(k, jcol, icol)* &
                     vertical_grid%rdzn(k+1) * &
                     (current_state%q(iqr)%data(k,jcol,icol) - current_state%q(iqr)%data(k+1,jcol,icol))
             enddo
          endif
          if (iqi > 0) then
             do k=1, current_state%local_grid%size(Z_INDEX)-1
                wqi_sg_tot(k) = wqi_sg_tot(k) + &
                     current_state%diff_coefficient%data(k, jcol, icol)* &
                     vertical_grid%rdzn(k+1) * &
                     (current_state%q(iqi)%data(k,jcol,icol) - current_state%q(iqi)%data(k+1,jcol,icol))
             enddo
          endif
          if (iqi > 0) then
             do k=1, current_state%local_grid%size(Z_INDEX)-1
                wqs_sg_tot(k) = wqs_sg_tot(k) + &
                     current_state%diff_coefficient%data(k, jcol, icol)* &
                     vertical_grid%rdzn(k+1) * &
                     (current_state%q(iqs)%data(k,jcol,icol) - current_state%q(iqs)%data(k+1,jcol,icol))
             enddo
          endif
          if (iqg > 0) then
             do k=1, current_state%local_grid%size(Z_INDEX)-1
                wqg_sg_tot(k) = wqg_sg_tot(k) + &
                     current_state%diff_coefficient%data(k, jcol, icol)* &
                     vertical_grid%rdzn(k+1) * &
                     (current_state%q(iqg)%data(k,jcol,icol) - current_state%q(iqg)%data(k+1,jcol,icol))
             enddo
          endif
       endif

! Level 1 wind speed - needed in various places

       U_1=SQRT(current_state%u%data(2,jcol,icol) * &
                current_state%u%data(2,jcol,icol) + &
                current_state%v%data(2,jcol,icol) * &
                current_state%v%data(2,jcol,icol))
       U_1_bar=SQRT(current_state%global_grid%configuration%vertical%olubar(2)*&
                    current_state%global_grid%configuration%vertical%olubar(2)+&
                    current_state%global_grid%configuration%vertical%olvbar(2)*&
                    current_state%global_grid%configuration%vertical%olvbar(2))

      ! Do p levels - calculate subgrid TKE diagnostics based on tau, interpolate to w levels
       do k=1, current_state%local_grid%size(Z_INDEX)-1
         rho(k) = current_state%global_grid%configuration%vertical%rho(k)
       end do

! Buoyancy coefficient on u (zn) levels
! note: this is not official MONC elsewhere
       do k=1, current_state%local_grid%size(Z_INDEX)
         buoy_cof(k)=G/current_state%global_grid%configuration%vertical%thref(k) 
       end do


!      =======================================================
!      Calculation of subgrid diagnostics dependent on the dissipation.
!
!      Rationalize this with smagorinsky.
       ri_crit = 0.25_DEFAULT_PRECISION
!
       ssq=calculate_half_squared_strain_rate(current_state, current_state%u, current_state%v, current_state%w)
       richardson_number=calculate_richardson_number(current_state, ssq, current_state%th, current_state%q)
!
       do k=2, current_state%local_grid%size(Z_INDEX)-1
          !        Mimic LEM: I think this is the square of the mixing length (Brown et al. 94, Eq 7)
          elamr_sq(k) = 0.0

          if ( ( &
               current_state%vis_coefficient%data(k, jcol, icol) &
               > 0.0) .and. (richardson_number(k) < ri_crit) ) then
             elamr_sq(k) = &
                  current_state%vis_coefficient%data(k, jcol, icol) / &
                  sqrt( 1.0 - richardson_number(k) * &
                  current_state%diff_coefficient%data(k, jcol, icol) / &
                  current_state%vis_coefficient%data(k, jcol, icol) )
          end if
          if (l_lem_dissipation_rate) then
             dissipation_rate = ssq(k) * ( &
                  current_state%vis_coefficient%data(k, jcol, icol) - &
                  richardson_number(k) * &
                  current_state%diff_coefficient%data(k, jcol, icol) )
          else ! CIRCLE-A code 
             if (use_Ri_for_buoyant_prod) then
      
                buoy_prod = -ssq(k) * richardson_number(k) *  &
                     current_state%diff_coefficient%data(k, jcol, icol)
                
             else
                if (current_state%th%active .and. &
                     .not. current_state%passive_q .and. current_state%number_q_fields .gt. 0) then
                   buoy_prod =  &
                        -current_state%diff_coefficient%data(k,jcol,icol) * &
                        !d/dz - so inside needs to be on p points
                        current_state%global_grid%configuration%vertical%rdzn(k+1) * &
                        !G/th_ref on p
                        (buoy_cof(k+1) * &  ! 1
                        !(th')
                        (current_state%th%data(k+1,jcol,icol) + & ! 2
                        current_state%global_grid%configuration%vertical%thref(k+1) + & ! 2
                        !th_ref * (C_virtual* qv'(k)-qcl'(k))
                        current_state%global_grid%configuration%vertical%thref(k+1) * & ! 2
                        (C_virtual * current_state%q(1)%data(k+1,jcol,icol) -  & ! 3
                        current_state%q(2)%data(k+1,jcol,icol))) - & ! 1
                        !diff for d/dz
                        buoy_cof(k) * & ! 1
                        (current_state%th%data(k,jcol,icol) + & ! 2
                        current_state%global_grid%configuration%vertical%thref(k) + & ! 2
                        current_state%global_grid%configuration%vertical%thref(k) * & ! 2
                        (C_virtual * current_state%q(1)%data(k,jcol,icol) - & ! 3
                        current_state%q(2)%data(k,jcol,icol))) ) ! 0
                else
                   call log_master_log(LOG_ERROR, &
                        "Subgrid diags - buoy_prod calc needs theta and q active, STOP")
                endif
             endif
             ! Circle-A dissipation rate includes buoy_prod term, which is different
             ! to the original LEM version
             dissipation_rate = ssq(k) * &
                  current_state%vis_coefficient%data(k, jcol, icol) + &
                  buoy_prod

             buoysg_tot(k)=buoysg_tot(k) + buoy_prod

          endif
          
          dissipation_tot(k) =  dissipation_tot(k) + dissipation_rate
          
          ! PAC comment: Not clear the code below is correct (where doea a2_n fit in).
          !
          !        This constant needs to go somewhere accessible and where it can be set. The code in the LEM
          !        seems to go round and round in circles here.
          a2_n = 0.23
          subgrid_tke(k) = ( dissipation_rate * dissipation_rate * elamr_sq(k) ) ** (1.0/3.0) / a2_n
          
          !        Assume isotropy.
          tkesg_tot(k) = tkesg_tot(k) + subgrid_tke(k)
          uusg_tot(k) = uusg_tot(k) + (2.0/3.0) * subgrid_tke(k)
          vvsg_tot(k) = vvsg_tot(k) + (2.0/3.0) * subgrid_tke(k)
          wwsg_tot(k) = wwsg_tot(k) + (2.0/3.0) * subgrid_tke(k)
          wkesg_tot(k) = wkesg_tot(k) + subgrid_tke(k) * &
               current_state%w%data(k,jcol,icol)
          richardson_number_tot(k) = richardson_number_tot(k) +  richardson_number(k)
          richardson_squared_tot(k) = richardson_squared_tot(k) +  richardson_number(k)**2.0_DEFAULT_PRECISION
          
          ! Calculate the integrated subgrid TKE for the scalar.
          ! BAD place to put this but convenient
          subke_2d(target_y_index, target_x_index) = subke_2d(target_y_index, target_x_index) + &
               (subgrid_tke(k)*0.5_DEFAULT_PRECISION*                                           &
               vertical_grid%dzn(k+1)*                       &
               vertical_grid%rho(k))
       enddo
       ! *********************** Subgrid buoyant production***************************    
       ! Note - calculating on z levels (i.e. w) 
       ! So need dzn, rdzn.
       ! zn(k)-zn(k-1) is dzn(k), so nothing is stored in k=1
       !------
       if (.not. l_lem_dissipation_rate) then
          if (current_state%th%active .and. &
               .not. current_state%passive_q .and. current_state%number_q_fields .gt. 0) then
             k = 1
             buoy_prod =  &
                  -current_state%diff_coefficient%data(k,jcol,icol) * &
                  !d/dz - so inside needs to be on p points
                  current_state%global_grid%configuration%vertical%rdzn(k+1) * &
                  !G/th_ref on p
                  (buoy_cof(k+1) * &  ! 1
                  !(th')
                  (current_state%th%data(k+1,jcol,icol) + & ! 2
                  current_state%global_grid%configuration%vertical%thref(k+1) + & ! 2
                  !th_ref * (C_virtual* qv'(k)-qcl'(k))
                  current_state%global_grid%configuration%vertical%thref(k+1) * & ! 2
                  (C_virtual * current_state%q(1)%data(k+1,jcol,icol)   & ! 3
                  - current_state%q(2)%data(k+1,jcol,icol))) - & ! 1
                  !diff for d/dz
                  buoy_cof(k) * & ! 1
                  (current_state%th%data(k,jcol,icol) + & ! 2
                  current_state%global_grid%configuration%vertical%thref(k) + & ! 2
                  current_state%global_grid%configuration%vertical%thref(k) * & ! 2
                  (C_virtual * current_state%q(1)%data(k,jcol,icol) & ! 3
                  - current_state%q(2)%data(k,jcol,icol))) ) ! 0
             
             buoysg_tot(k)=buoysg_tot(k) + buoy_prod

          
             dissipation_tot(1) = dissipation_tot(1) + &
                  current_state%vis_coefficient%data(1, jcol, icol) &
                  * U_1 * U_1 / &
                  (current_state%global_grid%configuration%vertical%zn(2) * &
                  current_state%global_grid%configuration%vertical%zn(2)) + buoy_prod
          else
             call log_master_log(LOG_ERROR, &
                "Subgrid diags - buoy_prod calc needs theta and q active, STOP")
          endif
       endif
          
       ! Needs buoyancy part added later
       
       ! divide subke by altitude to make it column mean (as in the LEM)
       subke_2d(target_y_index, target_x_index) =  subke_2d(target_y_index, target_x_index)/ &
            vertical_grid%z(current_state%local_grid%size(Z_INDEX))
       !
       !        Thermal equivalents:
       !        Again, this needs to go somewhere better. pr_n is used in Smagorinsky anyway (as a local variable).
       
       if (current_state%th%active) then
          epsth=calculate_thermal_dissipation_rate(current_state, current_state%th)         
          ath2_n = 0.3
          pr_n   = 0.7
          do k=2, current_state%local_grid%size(Z_INDEX)-1
             ! Theta dissipation rate
             theta_dis_tot(k) = theta_dis_tot(k) + epsth(k)
             if (subgrid_tke(k) > 0.0) &
                  th2sg_tot(k) = th2sg_tot(k) + sqrt( a2_n * elamr_sq(k) / subgrid_tke(k) ) * &
                  epsth(k) / ( ath2_n**2 * pr_n)
          enddo
       endif
       
       
       ! *********************** Subgrid shear production***************************    
       ! Note - calculating on z levels (i.e. w) 
       !
       do k=2,current_state%local_grid%size(Z_INDEX)-2
          
          
          !Subgrid shear-------
          umean(k)=(current_state%global_grid%configuration%vertical%olubar(k+1) - &
               current_state%global_grid%configuration%vertical%olubar(k)) * &
               current_state%global_grid%configuration%vertical%rdzn(k+1)
          vmean(k)=(current_state%global_grid%configuration%vertical%olvbar(k+1) -  &
               current_state%global_grid%configuration%vertical%olvbar(k)) * &
               current_state%global_grid%configuration%vertical%rdzn(k+1)
          
          sg_shear_prod = current_state%vis_coefficient%data(k,jcol,icol)* &
               (S13(k)*umean(k)+S23(k)*vmean(k))
          
          ssub_tot(k)=ssub_tot(k) + sg_shear_prod
          
          !tau is at p points - so convert to w points
          uwsg(k)=0.5*(tau13(k)+tau13(k+1))/rho(k)
          vwsg(k)=0.5*(tau23(k)+tau23(k+1))/rho(k)
          
       end do
       
       k=1
             
       sg_shear_prod = current_state%vis_coefficient%data(1, jcol, icol) * &
            ( 0.5 * & 
            (current_state%u%data(2,jcol,icol-1)  + & 
            current_state%u%data(2,jcol,icol) ) * &
            current_state%global_grid%configuration%vertical%olubar(2) + &
            0.5 * & 
            (current_state%v%data(2,jcol-1,icol)  + & 
            current_state%v%data(2,jcol,icol) ) * &
            current_state%global_grid%configuration%vertical%olvbar(2) ) / &
            (current_state%global_grid%configuration%vertical%zn(2) * &
            current_state%global_grid%configuration%vertical%zn(2)) 
       
       ssub_tot(1)=ssub_tot(1) + sg_shear_prod
       
       ! *********************** Subgrid turbulence transport***************************    
       ! Note - calculating on z levels (i.e. w) 
       ! so need  u_i_prime_tau_i on p levels
       
       do k=2, current_state%local_grid%size(Z_INDEX)-1
          
          u_i_prime_tau_i(k) = ( 0.5_DEFAULT_PRECISION * &
               (current_state%u%data(k,jcol,icol-1)  + &
               current_state%u%data(k,jcol,icol)  ) - &
               current_state%global_grid%configuration%vertical%olubar(k)) * &
               0.5_DEFAULT_PRECISION * ( tau13(k-1)+tau13(k))
          
          u_i_prime_tau_i(k) = u_i_prime_tau_i(k) + ( 0.5_DEFAULT_PRECISION * &
               (current_state%v%data(k,jcol-1,icol)  + &
               current_state%v%data(k,jcol,icol)  ) - &
               current_state%global_grid%configuration%vertical%olvbar(k)) * &
               0.5_DEFAULT_PRECISION * ( tau23(k-1)+tau23(k))
          
          u_i_prime_tau_i(k) = u_i_prime_tau_i(k) + ( 0.5_DEFAULT_PRECISION * &
               (current_state%w%data(k-1,jcol,icol)+&
               current_state%w%data(k,jcol,icol))) *&
               tau33_on_p(k)
          
       end do
       
       u_i_prime_tau_i(1) = &
            ( 0.5_DEFAULT_PRECISION  * & 
            (current_state%u%data(2,jcol,icol-1)  + & 
            current_state%u%data(2,jcol,icol) ) - &
           current_state%global_grid%configuration%vertical%olubar(2) ) * &
           ( 0.5_DEFAULT_PRECISION  * &
           (current_state%u%data(2,jcol,icol-1)  + & 
           current_state%u%data(2,jcol,icol) ) )
       
       u_i_prime_tau_i(1) = u_i_prime_tau_i(1) + &
            ( 0.5_DEFAULT_PRECISION  * & 
            (current_state%v%data(2,jcol-1,icol)  + & 
            current_state%v%data(2,jcol,icol) ) - &
            current_state%global_grid%configuration%vertical%olvbar(2) ) * &
         ( 0.5_DEFAULT_PRECISION  * &
         (current_state%v%data(2,jcol-1,icol)  + & 
         current_state%v%data(2,jcol,icol) ) )
       
       ! Note - here we use the surface vis_coefficient to ensure exact cancellation in
       ! numerical budget
       u_i_prime_tau_i(1) = u_i_prime_tau_i(1) * &        
            current_state%vis_coefficient%data(1,jcol,icol) / & 
            current_state%global_grid%configuration%vertical%zn(2)  
       
       
       do k=2, current_state%local_grid%size(Z_INDEX)-2
          sed_tot(k)=sed_tot(k) + (u_i_prime_tau_i(k+1)-u_i_prime_tau_i(k)) * &
               current_state%global_grid%configuration%vertical%rdzn(k+1) / &
               current_state%global_grid%configuration%vertical%rho(k)
       end do
       
       sed_tot(1)=sed_tot(1) + u_i_prime_tau_i(1) / &
            current_state%global_grid%configuration%vertical%zn(2)
       
       !      =======================================================
    endif ! (.not. current_state%halo_column)
  end subroutine timestep_callback

  !> Field information retrieval callback, this returns information for a specific components published field
  !! @param current_state Current model state
  !! @param name The name of the field to retrieve information for
  !! @param field_information Populated with information about the field
  subroutine field_information_retrieval_callback(current_state, name, field_information)
    type(model_state_type), target, intent(inout) :: current_state
    character(len=*), intent(in) :: name
    type(component_field_information_type), intent(out) :: field_information

    field_information%field_type=COMPONENT_ARRAY_FIELD_TYPE
    field_information%data_type=COMPONENT_DOUBLE_DATA_TYPE
    if (name .eq. 'subke') then
       field_information%number_dimensions=2
       field_information%dimension_sizes(1)=current_state%local_grid%size(Y_INDEX)
       field_information%dimension_sizes(2)=current_state%local_grid%size(X_INDEX)
    else
       field_information%number_dimensions=1
       field_information%dimension_sizes(1)=current_state%local_grid%size(Z_INDEX)
    endif
    
    if (name .eq. "th2sg_total_local" .or. name .eq. "wtsg_total_local" &
         .or. name .eq. "theta_dis_total_local" ) then
      field_information%enabled=current_state%th%active
    else if (name .eq. "wqv_sg_total_local" .or. name .eq. "wql_sg_total_local") then
       field_information%enabled= (.not.current_state%passive_q) .and. &
            (current_state%number_q_fields .gt. 0)
    else if (name .eq. "wqr_sg_total_local") then
       field_information%enabled= current_state%rain_water_mixing_ratio_index .gt. 0
    else if (name .eq. "wqi_sg_total_local") then
       field_information%enabled= current_state%ice_water_mixing_ratio_index .gt. 0
    else if (name .eq. "wqs_sg_total_local") then
       field_information%enabled= current_state%snow_water_mixing_ratio_index .gt. 0
    else if (name .eq. "wqg_sg_total_local") then
       field_information%enabled= current_state%graupel_water_mixing_ratio_index .gt. 0
       
!   ========================================================================
!   2nd stream
    else if (name .eq. "i_th2sg_total_local") then
      field_information%enabled=current_state%th%active
    else if (name .eq. "i_wqsg_total_local") then
      field_information%enabled= (.not.current_state%passive_q) .and. &
        (current_state%liquid_water_mixing_ratio_index > 0)
!   ========================================================================
    else 
      field_information%enabled=.true.
    end if
  end subroutine field_information_retrieval_callback

  !> Field value retrieval callback, this returns the value of a specific published field
  !! @param current_state Current model state
  !! @param name The name of the field to retrieve the value for
  !! @param field_value Populated with the value of the field
  subroutine field_value_retrieval_callback(current_state, name, field_value)
    type(model_state_type), target, intent(inout) :: current_state
    character(len=*), intent(in) :: name
    type(component_field_value_type), intent(out) :: field_value
    
    integer :: k
    
    if (name .eq. "uwsg_total_local") then
       allocate(field_value%real_1d_array(current_state%local_grid%size(Z_INDEX)))
       do k = 1, current_state%local_grid%size(Z_INDEX)
          field_value%real_1d_array(k)=uwsg_tot(k)
       enddo
    else if (name .eq. "sed_total_local") then
       allocate(field_value%real_1d_array(current_state%local_grid%size(Z_INDEX)))
       do k = 1, current_state%local_grid%size(Z_INDEX)
          field_value%real_1d_array(k)=sed_tot(k)
       enddo
    else if (name .eq. "ssub_total_local") then
       allocate(field_value%real_1d_array(current_state%local_grid%size(Z_INDEX)))
       do k = 1, current_state%local_grid%size(Z_INDEX)
          field_value%real_1d_array(k)=ssub_tot(k)
       enddo
    else if (name .eq. "buoysg_total_local") then
       allocate(field_value%real_1d_array(current_state%local_grid%size(Z_INDEX)))
       do k = 1, current_state%local_grid%size(Z_INDEX)
          field_value%real_1d_array(k)=buoysg_tot(k)
       enddo
    else if (name .eq. "vwsg_total_local") then
       allocate(field_value%real_1d_array(current_state%local_grid%size(Z_INDEX)))
       do k = 1, current_state%local_grid%size(Z_INDEX)
          field_value%real_1d_array(k)=vwsg_tot(k)
       enddo
    else if (name .eq. "uusg_total_local") then
       allocate(field_value%real_1d_array(current_state%local_grid%size(Z_INDEX)))
       do k = 1, current_state%local_grid%size(Z_INDEX)
          field_value%real_1d_array(k)=uusg_tot(k)
       enddo
    else if (name .eq. "vvsg_total_local") then
       allocate(field_value%real_1d_array(current_state%local_grid%size(Z_INDEX)))
       do k = 1, current_state%local_grid%size(Z_INDEX)
          field_value%real_1d_array(k)=vvsg_tot(k)
       enddo
    else if (name .eq. "wwsg_total_local") then
       allocate(field_value%real_1d_array(current_state%local_grid%size(Z_INDEX)))
       do k = 1, current_state%local_grid%size(Z_INDEX)
          field_value%real_1d_array(k)=wwsg_tot(k)
       enddo
    else if (name .eq. "tkesg_total_local") then
       allocate(field_value%real_1d_array(current_state%local_grid%size(Z_INDEX)))
       do k = 1, current_state%local_grid%size(Z_INDEX)
          field_value%real_1d_array(k)=tkesg_tot(k)
       enddo
    else if (name .eq. "wkesg_total_local") then
       allocate(field_value%real_1d_array(current_state%local_grid%size(Z_INDEX)))
       do k = 1, current_state%local_grid%size(Z_INDEX)
          field_value%real_1d_array(k)=wkesg_tot(k)
       enddo
    else if (name .eq. "viscosity_coef_total_local") then
       allocate(field_value%real_1d_array(current_state%local_grid%size(Z_INDEX)))
       do k = 1, current_state%local_grid%size(Z_INDEX)
          field_value%real_1d_array(k)=vis_coef_tot(k)
       enddo   
    else if (name .eq. "diffusion_coef_total_local") then
       allocate(field_value%real_1d_array(current_state%local_grid%size(Z_INDEX)))
       do k = 1, current_state%local_grid%size(Z_INDEX)
          field_value%real_1d_array(k)=diff_coef_tot(k)
       enddo  
    else if (name .eq. "richardson_number_total_local") then
       allocate(field_value%real_1d_array(current_state%local_grid%size(Z_INDEX)))
       do k = 1, current_state%local_grid%size(Z_INDEX)
          field_value%real_1d_array(k)=richardson_number_tot(k)
       enddo  
    else if (name .eq. "richardson_squared_total_local") then
       allocate(field_value%real_1d_array(current_state%local_grid%size(Z_INDEX)))
       do k = 1, current_state%local_grid%size(Z_INDEX)
          field_value%real_1d_array(k)=richardson_squared_tot(k)
       enddo
    else if (name .eq. "dissipation_total_local") then
       allocate(field_value%real_1d_array(current_state%local_grid%size(Z_INDEX)))
       do k = 1, current_state%local_grid%size(Z_INDEX)
          field_value%real_1d_array(k)=dissipation_tot(k)
       enddo
    else if (name .eq. "wtsg_total_local") then
       allocate(field_value%real_1d_array(current_state%local_grid%size(Z_INDEX)))
       do k = 1, current_state%local_grid%size(Z_INDEX)
          field_value%real_1d_array(k)=wtsg_tot(k)
       enddo
    else if (name .eq. "th2sg_total_local") then
       allocate(field_value%real_1d_array(current_state%local_grid%size(Z_INDEX)))
       do k = 1, current_state%local_grid%size(Z_INDEX)
          field_value%real_1d_array(k)=th2sg_tot(k)
       enddo
    else if (name .eq. "theta_dis_total_local") then
       allocate(field_value%real_1d_array(current_state%local_grid%size(Z_INDEX)))
       do k = 1, current_state%local_grid%size(Z_INDEX)
          field_value%real_1d_array(k)=theta_dis_tot(k)
       enddo   
    else if (name .eq. "wqv_sg_total_local") then
       allocate(field_value%real_1d_array(current_state%local_grid%size(Z_INDEX)))
       do k = 1, current_state%local_grid%size(Z_INDEX)
          field_value%real_1d_array(k)=wqv_sg_tot(k)
       enddo
    else if (name .eq. "wql_sg_total_local") then
       allocate(field_value%real_1d_array(current_state%local_grid%size(Z_INDEX)))
       do k = 1, current_state%local_grid%size(Z_INDEX)
          field_value%real_1d_array(k)=wql_sg_tot(k)
       enddo
    else if (name .eq. "wqr_sg_total_local") then
       allocate(field_value%real_1d_array(current_state%local_grid%size(Z_INDEX)))
       do k = 1, current_state%local_grid%size(Z_INDEX)
          field_value%real_1d_array(k)=wqr_sg_tot(k)
       enddo
    else if (name .eq. "wqi_sg_total_local") then
       allocate(field_value%real_1d_array(current_state%local_grid%size(Z_INDEX)))
       do k = 1, current_state%local_grid%size(Z_INDEX)
          field_value%real_1d_array(k)=wqi_sg_tot(k)
       enddo
    else if (name .eq. "wqs_sg_total_local") then
       allocate(field_value%real_1d_array(current_state%local_grid%size(Z_INDEX)))
       do k = 1, current_state%local_grid%size(Z_INDEX)
          field_value%real_1d_array(k)=wqs_sg_tot(k)
       enddo
    else if (name .eq. "wqg_sg_total_local") then
       allocate(field_value%real_1d_array(current_state%local_grid%size(Z_INDEX)))
       do k = 1, current_state%local_grid%size(Z_INDEX)
          field_value%real_1d_array(k)=wqg_sg_tot(k)
       enddo
! =====================================================
!   2nd stream
    else if (name .eq. "i_uwsg_total_local") then
       allocate(field_value%real_1d_array(current_state%local_grid%size(Z_INDEX)))
       do k = 1, current_state%local_grid%size(Z_INDEX)
          field_value%real_1d_array(k)=uwsg_tot(k)
       enddo
    else if (name .eq. "i_vwsg_total_local") then
       allocate(field_value%real_1d_array(current_state%local_grid%size(Z_INDEX)))
       do k = 1, current_state%local_grid%size(Z_INDEX)
          field_value%real_1d_array(k)=vwsg_tot(k)
       enddo
    else if (name .eq. "i_uusg_total_local") then
       allocate(field_value%real_1d_array(current_state%local_grid%size(Z_INDEX)))
       do k = 1, current_state%local_grid%size(Z_INDEX)
          field_value%real_1d_array(k)=uusg_tot(k)
       enddo
    else if (name .eq. "i_vvsg_total_local") then
       allocate(field_value%real_1d_array(current_state%local_grid%size(Z_INDEX)))
       do k = 1, current_state%local_grid%size(Z_INDEX)
          field_value%real_1d_array(k)=vvsg_tot(k)
       enddo
    else if (name .eq. "i_wwsg_total_local") then
       allocate(field_value%real_1d_array(current_state%local_grid%size(Z_INDEX)))
       do k = 1, current_state%local_grid%size(Z_INDEX)
          field_value%real_1d_array(k)=wwsg_tot(k)
       enddo
    else if (name .eq. "i_tkesg_total_local") then
       allocate(field_value%real_1d_array(current_state%local_grid%size(Z_INDEX)))
       do k = 1, current_state%local_grid%size(Z_INDEX)
          field_value%real_1d_array(k)=tkesg_tot(k)
       enddo
    else if (name .eq. "i_wtsg_total_local") then
       allocate(field_value%real_1d_array(current_state%local_grid%size(Z_INDEX)))
       do k = 1, current_state%local_grid%size(Z_INDEX)
          field_value%real_1d_array(k)=wtsg_tot(k)
       enddo
    else if (name .eq. "i_th2sg_total_local") then
       allocate(field_value%real_1d_array(current_state%local_grid%size(Z_INDEX)))
       do k = 1, current_state%local_grid%size(Z_INDEX)
          field_value%real_1d_array(k)=th2sg_tot(k)
       enddo
    else if (name .eq. "i_wqsg_total_local") then
       allocate(field_value%real_1d_array(current_state%local_grid%size(Z_INDEX)))
       do k = 1, current_state%local_grid%size(Z_INDEX)
          field_value%real_1d_array(k)=wqsg_tot(k)
       enddo
 ! =====================================================
    else if (name .eq. 'subke') then
       allocate(field_value%real_2d_array(current_state%local_grid%size(Y_INDEX), &
           current_state%local_grid%size(X_INDEX)))
       field_value%real_2d_array(:,:)=subke_2d(:,:)
    end if
  end subroutine field_value_retrieval_callback
  
end module subgrid_profile_diagnostics_mod
