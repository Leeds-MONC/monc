module subgrid_profile_diagnostics_mod
  use monc_component_mod, only : COMPONENT_SCALAR_FIELD_TYPE, COMPONENT_DOUBLE_DATA_TYPE, &
       COMPONENT_ARRAY_FIELD_TYPE, COMPONENT_INTEGER_DATA_TYPE, &
       component_descriptor_type, component_field_value_type, component_field_information_type
  use grids_mod, only : vertical_grid_configuration_type, Z_INDEX, Y_INDEX, X_INDEX
  use optionsdatabase_mod, only : options_get_real
  use state_mod, only : model_state_type
  use datadefn_mod, only : DEFAULT_PRECISION
  use q_indices_mod, only: get_q_index, standard_q_names
  use registry_mod, only : is_component_enabled
  use logging_mod, only : LOG_ERROR, LOG_WARN, log_master_log
  use smagorinsky_mod, only: calculate_half_squared_strain_rate, &
                             calculate_richardson_number,        &
                             calculate_thermal_dissipation_rate

  implicit none

#ifndef TEST_MODE
  private
#endif

  integer :: total_points, iqv=0, iql=0, iqr=0, iqi=0, iqs=0,    &
       iqg=0

! diagnostic variables
  real(kind=DEFAULT_PRECISION), dimension(:), allocatable ::     &
       uwsg_tot, vwsg_tot, uusg_tot, vvsg_tot, wwsg_tot,         &
       tkesg_tot, wtsg_tot, th2sg_tot, wqsg_tot, wkesg_tot,      &
       theta_dis_tot, vis_coef_tot, diff_coef_tot,               &
       richardson_number_tot, richardson_squared_tot,            &
       dis_tot,                                                  &
       ! subgrid moisture fluxes
       wqv_sg_tot, wql_sg_tot, wqr_sg_tot, wqi_sg_tot,           &
       wqs_sg_tot, wqg_sg_tot
! These arrays should in due course be allocated conditionally,
! but let's just get it working first.
  real(kind=DEFAULT_PRECISION), dimension(:), allocatable ::     &
       dissipation, epsth, ssq, elamr_sq, richardson_number,     &
       subgrid_tke
  real(kind=DEFAULT_PRECISION), dimension(:,:), allocatable ::   &
       subke_2d
! Constants here provisionally
  real(kind=DEFAULT_PRECISION) :: a2_n, ath2_n, pr_n, ri_crit
  real(kind=DEFAULT_PRECISION) :: qlcrit

  type(vertical_grid_configuration_type) :: vertical_grid
  
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
    allocate(subgrid_profile_diagnostics_get_descriptor%published_fields(22+9))

    subgrid_profile_diagnostics_get_descriptor%published_fields(1)="uwsg_total_local"
    subgrid_profile_diagnostics_get_descriptor%published_fields(2)="vwsg_total_local"
    subgrid_profile_diagnostics_get_descriptor%published_fields(3)="uusg_total_local"
    subgrid_profile_diagnostics_get_descriptor%published_fields(4)="vvsg_total_local"
    subgrid_profile_diagnostics_get_descriptor%published_fields(5)="wwsg_total_local"
    subgrid_profile_diagnostics_get_descriptor%published_fields(6)="tkesg_total_local"
    subgrid_profile_diagnostics_get_descriptor%published_fields(7)="wtsg_total_local"
    subgrid_profile_diagnostics_get_descriptor%published_fields(8)="th2sg_total_local"
    subgrid_profile_diagnostics_get_descriptor%published_fields(9)="wkesg_total_local"
    subgrid_profile_diagnostics_get_descriptor%published_fields(10)="theta_dis_total_local"
    subgrid_profile_diagnostics_get_descriptor%published_fields(11)="viscosity_coef_total_local"
    subgrid_profile_diagnostics_get_descriptor%published_fields(12)="diffusion_coef_total_local"
    subgrid_profile_diagnostics_get_descriptor%published_fields(13)="richardson_number_total_local"
    subgrid_profile_diagnostics_get_descriptor%published_fields(14)="richardson_squared_total_local"
    subgrid_profile_diagnostics_get_descriptor%published_fields(15)="dissipation_total_local"
    subgrid_profile_diagnostics_get_descriptor%published_fields(16)="subke"
    subgrid_profile_diagnostics_get_descriptor%published_fields(17)="wqv_sg_total_local"
    subgrid_profile_diagnostics_get_descriptor%published_fields(18)="wql_sg_total_local"
    subgrid_profile_diagnostics_get_descriptor%published_fields(19)="wqr_sg_total_local"
    subgrid_profile_diagnostics_get_descriptor%published_fields(20)="wqi_sg_total_local"
    subgrid_profile_diagnostics_get_descriptor%published_fields(21)="wqs_sg_total_local"
    subgrid_profile_diagnostics_get_descriptor%published_fields(22)="wqg_sg_total_local"
    

!   =====================================================
!   2nd, provisionally instantaneous, stream

    subgrid_profile_diagnostics_get_descriptor%published_fields(22+1)="i_uwsg_total_local"
    subgrid_profile_diagnostics_get_descriptor%published_fields(22+2)="i_vwsg_total_local"
    subgrid_profile_diagnostics_get_descriptor%published_fields(22+3)="i_uusg_total_local"
    subgrid_profile_diagnostics_get_descriptor%published_fields(22+4)="i_vvsg_total_local"
    subgrid_profile_diagnostics_get_descriptor%published_fields(22+5)="i_wwsg_total_local"
    subgrid_profile_diagnostics_get_descriptor%published_fields(22+6)="i_tkesg_total_local"
    subgrid_profile_diagnostics_get_descriptor%published_fields(22+7)="i_wtsg_total_local"
    subgrid_profile_diagnostics_get_descriptor%published_fields(22+8)="i_th2sg_total_local"
    subgrid_profile_diagnostics_get_descriptor%published_fields(22+9)="i_wqsg_total_local"
!   =====================================================

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
         ,   dis_tot(current_state%local_grid%size(Z_INDEX)) )

!   Allocation of dissipation should in due course be made conditional.
!   Check sonsistency of any changes with the deallocations later.
    allocate(dissipation(current_state%local_grid%size(Z_INDEX)))
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
    
  end subroutine initialisation_callback  

  subroutine timestep_callback(current_state)
    type(model_state_type), target, intent(inout) :: current_state

    integer :: k, i
    integer :: jcol, icol, target_x_index, target_y_index
    
    jcol=current_state%column_local_y
    icol=current_state%column_local_x
    target_y_index=jcol-current_state%local_grid%halo_size(Y_INDEX)
    target_x_index=icol-current_state%local_grid%halo_size(X_INDEX)
    
    if (current_state%first_timestep_column) then
       uwsg_tot(:)     = 0.0_DEFAULT_PRECISION
       vwsg_tot(:)     = 0.0_DEFAULT_PRECISION
       uusg_tot(:)     = 0.0_DEFAULT_PRECISION
       vvsg_tot(:)     = 0.0_DEFAULT_PRECISION
       wwsg_tot(:)     = 0.0_DEFAULT_PRECISION
       tkesg_tot(:)    = 0.0_DEFAULT_PRECISION
       vis_coef_tot(:) = 0.0_DEFAULT_PRECISION
       diff_coef_tot(:)= 0.0_DEFAULT_PRECISION
       dis_tot(:)      = 0.0_DEFAULT_PRECISION
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
          !
          !        This was called EPS in the LEM
          dissipation(k) = ssq(k) * ( &
               current_state%vis_coefficient%data(k, jcol, icol) - &
               richardson_number(k) * &
               current_state%diff_coefficient%data(k, jcol, icol) )
          !
          !        This constant needs to go somewhere accessible and where it can be set. The code in the LEM
          !        seems to go round and round in circles here.
          a2_n = 0.23
          subgrid_tke(k) = ( dissipation(k) * dissipation(k) * elamr_sq(k) ) ** (1.0/3.0) / a2_n
          !
          !        Assume isotropy.
          tkesg_tot(k) = tkesg_tot(k) + subgrid_tke(k)
          uusg_tot(k) = uusg_tot(k) + (2.0/3.0) * subgrid_tke(k)
          vvsg_tot(k) = vvsg_tot(k) + (2.0/3.0) * subgrid_tke(k)
          wwsg_tot(k) = wwsg_tot(k) + (2.0/3.0) * subgrid_tke(k)
          wkesg_tot(k) = wkesg_tot(k) + subgrid_tke(k) * &
               current_state%w%data(k,jcol,icol)
          richardson_number_tot(k) = richardson_number_tot(k) +  richardson_number(k)
          richardson_squared_tot(k) = richardson_squared_tot(k) +  richardson_number(k)**2.0_DEFAULT_PRECISION
          dis_tot(k) = dis_tot(k) + dissipation(k)
          ! Calculate the integrated subgrid TKE for the scalar.
          ! BAD place to put this but convenient
          subke_2d(target_y_index, target_x_index) = subke_2d(target_y_index, target_x_index) + &
               (subgrid_tke(k)*0.5_DEFAULT_PRECISION*                                           &
               vertical_grid%dzn(k+1)*                       &
               vertical_grid%rho(k))
       enddo
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
!       
!      =======================================================
    endif
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
          field_value%real_1d_array(k)=dis_tot(k)
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
