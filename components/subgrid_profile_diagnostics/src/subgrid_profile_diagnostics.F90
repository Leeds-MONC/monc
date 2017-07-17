module subgrid_profile_diagnostics_mod
  use monc_component_mod, only : COMPONENT_SCALAR_FIELD_TYPE, COMPONENT_DOUBLE_DATA_TYPE, &
       COMPONENT_ARRAY_FIELD_TYPE, COMPONENT_INTEGER_DATA_TYPE, &
       component_descriptor_type, component_field_value_type, component_field_information_type
  use grids_mod, only : Z_INDEX, Y_INDEX, X_INDEX
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

  integer :: total_points, iqv, iql

  real(kind=DEFAULT_PRECISION), dimension(:), allocatable ::     &
       uwsg_tot, vwsg_tot, uusg_tot, vvsg_tot, wwsg_tot,         &
       tkesg_tot, wtsg_tot, th2sg_tot, wqsg_tot
! These arrays should in due course be allocated conditionally,
! but let's just get it working first.
  real(kind=DEFAULT_PRECISION), dimension(:), allocatable ::     &
       dissipation, epsth, ssq, elamr_sq, richardson_number,     &
       subgrid_tke
! Constants here provisionally
  real(kind=DEFAULT_PRECISION) :: a2_n, ath2_n, pr_n, ri_crit
  real(kind=DEFAULT_PRECISION) :: qlcrit

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
    allocate(subgrid_profile_diagnostics_get_descriptor%published_fields(2*(5+2+9)))

    subgrid_profile_diagnostics_get_descriptor%published_fields(1)="uwsg_total_local"
    subgrid_profile_diagnostics_get_descriptor%published_fields(2)="vwsg_total_local"
    subgrid_profile_diagnostics_get_descriptor%published_fields(3)="uusg_total_local"
    subgrid_profile_diagnostics_get_descriptor%published_fields(4)="vvsg_total_local"
    subgrid_profile_diagnostics_get_descriptor%published_fields(5)="wwsg_total_local"
    subgrid_profile_diagnostics_get_descriptor%published_fields(6)="tkesg_total_local"
    subgrid_profile_diagnostics_get_descriptor%published_fields(7)="wtsg_total_local"
    subgrid_profile_diagnostics_get_descriptor%published_fields(8)="th2sg_total_local"
    subgrid_profile_diagnostics_get_descriptor%published_fields(9)="wqsg_total_local"

!   =====================================================
!   2nd, provisionally instantaneous, stream

    subgrid_profile_diagnostics_get_descriptor%published_fields(9+1)="i_uwsg_total_local"
    subgrid_profile_diagnostics_get_descriptor%published_fields(9+2)="i_vwsg_total_local"
    subgrid_profile_diagnostics_get_descriptor%published_fields(9+3)="i_uusg_total_local"
    subgrid_profile_diagnostics_get_descriptor%published_fields(9+4)="i_vvsg_total_local"
    subgrid_profile_diagnostics_get_descriptor%published_fields(9+5)="i_wwsg_total_local"
    subgrid_profile_diagnostics_get_descriptor%published_fields(9+6)="i_tkesg_total_local"
    subgrid_profile_diagnostics_get_descriptor%published_fields(9+7)="i_wtsg_total_local"
    subgrid_profile_diagnostics_get_descriptor%published_fields(9+8)="i_th2sg_total_local"
    subgrid_profile_diagnostics_get_descriptor%published_fields(9+9)="i_wqsg_total_local"
!   =====================================================

  end function subgrid_profile_diagnostics_get_descriptor

  subroutine initialisation_callback(current_state)
    type(model_state_type), target, intent(inout) :: current_state

    integer :: k

    if (.not. is_component_enabled(current_state%options_database, "smagorinsky")) then
      call log_master_log(LOG_ERROR, "Subgrid model diags requested but subgrid model not enabled - check config")
    end if 

    total_points=(current_state%local_grid%size(Y_INDEX) * current_state%local_grid%size(X_INDEX))
    
    allocate(uwsg_tot(current_state%local_grid%size(Z_INDEX))   &
         ,   vwsg_tot(current_state%local_grid%size(Z_INDEX))   &
         ,   uusg_tot(current_state%local_grid%size(Z_INDEX))   &
         ,   vvsg_tot(current_state%local_grid%size(Z_INDEX))   &
         ,   wwsg_tot(current_state%local_grid%size(Z_INDEX))   &
         ,   tkesg_tot(current_state%local_grid%size(Z_INDEX))  &
         ,   wtsg_tot(current_state%local_grid%size(Z_INDEX))   &
         ,   th2sg_tot(current_state%local_grid%size(Z_INDEX))  &
         ,   wqsg_tot(current_state%local_grid%size(Z_INDEX)) )
!   Allocation of dissipation should in due course be made conditional.
!   Check sonsistency of any changes with the deallocations later.
    allocate(dissipation(current_state%local_grid%size(Z_INDEX)))
    allocate(epsth(current_state%local_grid%size(Z_INDEX)))
    allocate(ssq(current_state%local_grid%size(Z_INDEX)))
    allocate(elamr_sq(current_state%local_grid%size(Z_INDEX)))
    allocate(richardson_number(current_state%local_grid%size(Z_INDEX)))
    allocate(subgrid_tke(current_state%local_grid%size(Z_INDEX)))
    
        
    if (.not. current_state%passive_q .and. current_state%number_q_fields .gt. 0) then 
       iqv=get_q_index(standard_q_names%VAPOUR, 'subgrid_profile_diags')
       iql=get_q_index(standard_q_names%CLOUD_LIQUID_MASS, 'subgrid_profile_diags') 
       qlcrit=options_get_real(current_state%options_database, "qlcrit")  
    endif
   
  end subroutine initialisation_callback  

  subroutine timestep_callback(current_state)
    type(model_state_type), target, intent(inout) :: current_state

    integer :: k, i
 
    if (current_state%first_timestep_column) then
       uwsg_tot(:)     = 0.0_DEFAULT_PRECISION
       vwsg_tot(:)     = 0.0_DEFAULT_PRECISION
       uusg_tot(:)     = 0.0_DEFAULT_PRECISION
       vvsg_tot(:)     = 0.0_DEFAULT_PRECISION
       wwsg_tot(:)     = 0.0_DEFAULT_PRECISION
       tkesg_tot(:)    = 0.0_DEFAULT_PRECISION
       wtsg_tot(:)     = 0.0_DEFAULT_PRECISION
       th2sg_tot(:)    = 0.0_DEFAULT_PRECISION
       wqsg_tot(:)     = 0.0_DEFAULT_PRECISION
    end if

    if (.not. current_state%halo_column) then
!      Subgrid diagnostics
       do k=1, current_state%local_grid%size(Z_INDEX)-1
!
!         Field on kth rho-level, so gradient evaluated from k+1st and kth rhon-levels. Horizontally aligned with
!         w-points.
          uwsg_tot(k) = uwsg_tot(k) - 0.5 * &
            current_state%vis_coefficient%data(k,current_state%column_local_y,current_state%column_local_x-1) * &
            ( current_state%u%data(k+1,current_state%column_local_y,current_state%column_local_x-1) - &
              current_state%u%data(k,current_state%column_local_y,current_state%column_local_x-1) ) * &
            current_state%global_grid%configuration%vertical%rdzn(k+1)
          uwsg_tot(k) = uwsg_tot(k) - 0.5 * &
            current_state%vis_coefficient%data(k,current_state%column_local_y,current_state%column_local_x) * &
            ( current_state%u%data(k+1,current_state%column_local_y,current_state%column_local_x) - &
              current_state%u%data(k,current_state%column_local_y,current_state%column_local_x) ) * &
            current_state%global_grid%configuration%vertical%rdzn(k+1)
!
!         Field on kth rho-level, so gradient evaluated from k+1st and kth rhon-levels. Horizontally aligned with
!         w-points.
          vwsg_tot(k) = vwsg_tot(k) - 0.5 * &
            current_state%vis_coefficient%data(k,current_state%column_local_y-1,current_state%column_local_x) * &
            ( current_state%v%data(k+1,current_state%column_local_y-1,current_state%column_local_x) - &
              current_state%v%data(k,current_state%column_local_y-1,current_state%column_local_x) ) * &
            current_state%global_grid%configuration%vertical%rdzn(k+1)
          vwsg_tot(k) = vwsg_tot(k) - 0.5 * &
            current_state%vis_coefficient%data(k,current_state%column_local_y,current_state%column_local_x) * &
            ( current_state%v%data(k+1,current_state%column_local_y,current_state%column_local_x) - &
              current_state%v%data(k,current_state%column_local_y,current_state%column_local_x) ) * &
            current_state%global_grid%configuration%vertical%rdzn(k+1)
          !
       enddo
       !         Field on kth rho-level, so gradient evaluated from k+1st and kth rhon-levels. Horizontally aligned with
       !         w-points.
       if (current_state%th%active) then
          do k=1, current_state%local_grid%size(Z_INDEX)-1
             wtsg_tot(k) = wtsg_tot(k) -  &
                  current_state%diff_coefficient%data(k,current_state%column_local_y,current_state%column_local_x) * &
                  ( current_state%th%data(k+1,current_state%column_local_y,current_state%column_local_x) + &
                  current_state%global_grid%configuration%vertical%thref(k+1) - &
                  current_state%th%data(k,current_state%column_local_y,current_state%column_local_x) - &
                  current_state%global_grid%configuration%vertical%thref(k) ) * &
                  current_state%global_grid%configuration%vertical%rdzn(k+1)
             !
          enddo
       endif

       if (current_state%th%active .and. .not. current_state%passive_q .and. current_state%number_q_fields .gt. 0) then
          
          do k=1, current_state%local_grid%size(Z_INDEX)-1
             ! not coded yet
          enddo
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
               current_state%vis_coefficient%data(k, current_state%column_local_y, current_state%column_local_x) &
               > 0.0) .and. (richardson_number(k) < ri_crit) ) then
             elamr_sq(k) = &
                  current_state%vis_coefficient%data(k, current_state%column_local_y, current_state%column_local_x) / &
                  sqrt( 1.0 - richardson_number(k) * &
                  current_state%diff_coefficient%data(k, current_state%column_local_y, current_state%column_local_x) / &
                  current_state%vis_coefficient%data(k, current_state%column_local_y, current_state%column_local_x) )
          end if
          !
          !        This was called EPS in the LEM
          dissipation(k) = ssq(k) * ( &
               current_state%vis_coefficient%data(k, current_state%column_local_y, current_state%column_local_x) - &
               richardson_number(k) * &
               current_state%diff_coefficient%data(k, current_state%column_local_y, current_state%column_local_x) )
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
       enddo

       !
       !        Thermal equivalents:
       !        Again, this needs to go somewhere better. pr_n is used in Smagorinsky anyway (as a local variable).

       if (current_state%th%active) then 
          epsth=calculate_thermal_dissipation_rate(current_state, current_state%th)          
          ath2_n = 0.3
          pr_n   = 0.7
          do k=2, current_state%local_grid%size(Z_INDEX)-1
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
    field_information%number_dimensions=1
    field_information%dimension_sizes(1)=current_state%local_grid%size(Z_INDEX)
    field_information%data_type=COMPONENT_DOUBLE_DATA_TYPE
    
    if (name .eq. "th2sg_total_local") then
      field_information%enabled=current_state%th%active
    else if (name .eq. "wqsg_total_local") then
      field_information%enabled= (.not.current_state%passive_q) .and. &
        (current_state%liquid_water_mixing_ratio_index > 0)
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
    else if (name .eq. "wqsg_total_local") then
       allocate(field_value%real_1d_array(current_state%local_grid%size(Z_INDEX)))
       do k = 1, current_state%local_grid%size(Z_INDEX)
          field_value%real_1d_array(k)=wqsg_tot(k)
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
    end if
  end subroutine field_value_retrieval_callback
end module subgrid_profile_diagnostics_mod
