module profile_diagnostics_mod
  use monc_component_mod, only : COMPONENT_SCALAR_FIELD_TYPE, COMPONENT_DOUBLE_DATA_TYPE, &
       COMPONENT_ARRAY_FIELD_TYPE, COMPONENT_INTEGER_DATA_TYPE, &
       component_descriptor_type, component_field_value_type, component_field_information_type
  use grids_mod, only : Z_INDEX, Y_INDEX, X_INDEX
  use optionsdatabase_mod, only : options_get_real
  use state_mod, only : model_state_type
  use datadefn_mod, only : DEFAULT_PRECISION
  use q_indices_mod, only: get_q_index, standard_q_names
  use saturation_mod, only: qsaturation

  implicit none

#ifndef TEST_MODE
  private
#endif

  integer :: total_points, iqv, iql
  real(kind=DEFAULT_PRECISION), dimension(:), allocatable ::     &
       tempfac, u_wind_tot, uprime_tot, v_wind_tot, vprime_tot,  &
       ww_tot, theta_tot, qv_tot, ql_tot, w_wind_tot, rh_tot,    &
       wq_tot, wtheta_tot, uw_tot, vw_tot, th2_tot,              &
       thref, prefn, rho, rhon, thinit, uinit, vinit
  real(kind=DEFAULT_PRECISION) :: qlcrit

  public profile_diagnostics_get_descriptor

contains

  !> Provides the component descriptor for the core to register
  !! @returns The descriptor describing this component
  type(component_descriptor_type) function profile_diagnostics_get_descriptor()
    profile_diagnostics_get_descriptor%name="profile_diagnostics"
    profile_diagnostics_get_descriptor%version=0.1

    profile_diagnostics_get_descriptor%initialisation=>initialisation_callback
    profile_diagnostics_get_descriptor%timestep=>timestep_callback

    profile_diagnostics_get_descriptor%field_value_retrieval=>field_value_retrieval_callback
    profile_diagnostics_get_descriptor%field_information_retrieval=>field_information_retrieval_callback
    allocate(profile_diagnostics_get_descriptor%published_fields(2*22))

    profile_diagnostics_get_descriptor%published_fields(1)="theta_total_local"
    profile_diagnostics_get_descriptor%published_fields(2)="vapour_mmr_total_local"
    profile_diagnostics_get_descriptor%published_fields(3)="liquid_mmr_total_local"
    profile_diagnostics_get_descriptor%published_fields(4)="u_wind_total_local"
    profile_diagnostics_get_descriptor%published_fields(5)="uu_total_local"
    profile_diagnostics_get_descriptor%published_fields(6)="v_wind_total_local"
    profile_diagnostics_get_descriptor%published_fields(7)="vv_total_local"
    profile_diagnostics_get_descriptor%published_fields(8)="ww_total_local"
    profile_diagnostics_get_descriptor%published_fields(9)="w_wind_total_local"
    profile_diagnostics_get_descriptor%published_fields(10)="rh_total_local"
    profile_diagnostics_get_descriptor%published_fields(11)="thref_local"
    profile_diagnostics_get_descriptor%published_fields(12)="prefn_local"
    profile_diagnostics_get_descriptor%published_fields(13)="rho_local"
    profile_diagnostics_get_descriptor%published_fields(14)="rhon_local" 
    profile_diagnostics_get_descriptor%published_fields(15)="thinit_local"
    profile_diagnostics_get_descriptor%published_fields(16)="uw_total_local"
    profile_diagnostics_get_descriptor%published_fields(17)="vw_total_local"
    profile_diagnostics_get_descriptor%published_fields(18)="wtheta_total_local"
    profile_diagnostics_get_descriptor%published_fields(19)="th2_total_local"
    profile_diagnostics_get_descriptor%published_fields(20)="wq_total_local"

    profile_diagnostics_get_descriptor%published_fields(21)="uinit_local"
    profile_diagnostics_get_descriptor%published_fields(22)="vinit_local" 

!   =====================================================
!   2nd, provisionally instantaneous, stream
    profile_diagnostics_get_descriptor%published_fields(22+1)="i_theta_total_local"
    profile_diagnostics_get_descriptor%published_fields(22+2)="i_vapour_mmr_total_local"
    profile_diagnostics_get_descriptor%published_fields(22+3)="i_liquid_mmr_total_local"
    profile_diagnostics_get_descriptor%published_fields(22+4)="i_u_wind_total_local"
    profile_diagnostics_get_descriptor%published_fields(22+4)="i_u_wind_total_local"
    profile_diagnostics_get_descriptor%published_fields(22+5)="i_uu_total_local"
    profile_diagnostics_get_descriptor%published_fields(22+6)="i_v_wind_total_local"
    profile_diagnostics_get_descriptor%published_fields(22+7)="i_vv_total_local"
    profile_diagnostics_get_descriptor%published_fields(22+8)="i_ww_total_local"
    profile_diagnostics_get_descriptor%published_fields(22+9)="i_w_wind_total_local"
    profile_diagnostics_get_descriptor%published_fields(22+10)="i_rh_total_local"
    profile_diagnostics_get_descriptor%published_fields(22+11)="i_thref_local"
    profile_diagnostics_get_descriptor%published_fields(22+12)="i_prefn_local"
    profile_diagnostics_get_descriptor%published_fields(22+13)="i_rho_local"
    profile_diagnostics_get_descriptor%published_fields(22+14)="i_rhon_local" 
    profile_diagnostics_get_descriptor%published_fields(22+15)="i_thinit_local"

    profile_diagnostics_get_descriptor%published_fields(22+16)="i_uw_total_local"
    profile_diagnostics_get_descriptor%published_fields(22+17)="i_vw_total_local"
    profile_diagnostics_get_descriptor%published_fields(22+18)="i_wtheta_total_local"
    profile_diagnostics_get_descriptor%published_fields(22+19)="i_th2_total_local"
    profile_diagnostics_get_descriptor%published_fields(22+20)="i_wq_total_local"

    profile_diagnostics_get_descriptor%published_fields(22+21)="i_uinit_local"
    profile_diagnostics_get_descriptor%published_fields(22+22)="i_vinit_local"    

  end function profile_diagnostics_get_descriptor

  subroutine initialisation_callback(current_state)
    type(model_state_type), target, intent(inout) :: current_state

    integer :: k

    total_points=(current_state%local_grid%size(Y_INDEX) * current_state%local_grid%size(X_INDEX))
    
    ! allocate local arrays for the horizontal wind averages
    allocate(u_wind_tot(current_state%local_grid%size(Z_INDEX)) &
         , v_wind_tot(current_state%local_grid%size(Z_INDEX))   &
         , ww_tot(current_state%local_grid%size(Z_INDEX))       &
         , w_wind_tot(current_state%local_grid%size(Z_INDEX))   &
         , prefn(current_state%local_grid%size(Z_INDEX))        &
         , rho(current_state%local_grid%size(Z_INDEX))          &
         , rhon(current_state%local_grid%size(Z_INDEX))         &
         , uinit(current_state%local_grid%size(Z_INDEX))       &
         , vinit(current_state%local_grid%size(Z_INDEX)))
    allocate(wq_tot(current_state%local_grid%size(Z_INDEX))     &
         ,   wtheta_tot(current_state%local_grid%size(Z_INDEX)) &
         ,   uw_tot(current_state%local_grid%size(Z_INDEX))     &
         ,   vw_tot(current_state%local_grid%size(Z_INDEX))     &
         ,   th2_tot(current_state%local_grid%size(Z_INDEX)) )
    
    if (allocated(current_state%global_grid%configuration%vertical%olubar)) then
      allocate(uprime_tot(current_state%local_grid%size(Z_INDEX)))
    end if
    
    if (allocated(current_state%global_grid%configuration%vertical%olvbar)) then
      allocate(vprime_tot(current_state%local_grid%size(Z_INDEX)))
    end if        

    if (current_state%th%active) then
        allocate(theta_tot(current_state%local_grid%size(Z_INDEX)) &
             , thref(current_state%local_grid%size(Z_INDEX))       &
             , thinit(current_state%local_grid%size(Z_INDEX)))
     endif
        
    if (.not. current_state%passive_q .and. current_state%number_q_fields .gt. 0) then                                     
       iqv=get_q_index(standard_q_names%VAPOUR, 'profile_diags')                                                                
       iql=get_q_index(standard_q_names%CLOUD_LIQUID_MASS, 'profile_diags')                                                     
       qlcrit=options_get_real(current_state%options_database, "qlcrit")                                                         
       allocate(qv_tot(current_state%local_grid%size(Z_INDEX))  &
         , ql_tot(current_state%local_grid%size(Z_INDEX)) ) 
       if (current_state%th%active) &
            allocate(rh_tot(current_state%local_grid%size(Z_INDEX)))
    endif

  end subroutine initialisation_callback  

  subroutine timestep_callback(current_state)
    type(model_state_type), target, intent(inout) :: current_state

    integer :: k, i
    real(kind=DEFAULT_PRECISION) :: cltop_col, clbas_col, qv, qc, TdegK, Pmb &
         , qs, exner
    real(kind=DEFAULT_PRECISION) :: uprime_w_local, vprime_w_local &
                                  , thprime_w_local, qprime_w_local

    if (current_state%first_timestep_column) then
       u_wind_tot(:) = 0.0_DEFAULT_PRECISION
       if (allocated(uprime_tot)) uprime_tot(:) = 0.0_DEFAULT_PRECISION
       v_wind_tot(:) = 0.0_DEFAULT_PRECISION
       if (allocated(vprime_tot)) vprime_tot(:) = 0.0_DEFAULT_PRECISION
       w_wind_tot(:) = 0.0_DEFAULT_PRECISION
       ww_tot(:) = 0.0_DEFAULT_PRECISION
       
       wq_tot(:)     = 0.0_DEFAULT_PRECISION
       wtheta_tot(:) = 0.0_DEFAULT_PRECISION
       uw_tot(:)     = 0.0_DEFAULT_PRECISION
       vw_tot(:)     = 0.0_DEFAULT_PRECISION
       th2_tot(:)    = 0.0_DEFAULT_PRECISION
       
       if (current_state%th%active) then 
          theta_tot(:)=0.0_DEFAULT_PRECISION
       endif
       if (.not. current_state%passive_q .and. &
            current_state%number_q_fields .gt. 0) then 
          qv_tot(:)=0.0_DEFAULT_PRECISION
          ql_tot(:)=0.0_DEFAULT_PRECISION
          if (current_state%th%active) &
               rh_tot(:) = 0.0_DEFAULT_PRECISION
       endif
    end if
    if (.not. current_state%halo_column) then
       ! work out the sum of u and v wind over local domaini
       do k=1, current_state%local_grid%size(Z_INDEX)
          u_wind_tot(k) = u_wind_tot(k) + & 
               (current_state%u%data(k,current_state%column_local_y,current_state%column_local_x)  &
                + current_state%ugal)
          if (allocated(uprime_tot)) then
            uprime_tot(k) = uprime_tot(k) + &
                 ((current_state%u%data(k,current_state%column_local_y,current_state%column_local_x) &
                 - (current_state%global_grid%configuration%vertical%olubar(k) - current_state%ugal))**2.)
          end if
          v_wind_tot(k) = v_wind_tot(k) + & 
               (current_state%v%data(k,current_state%column_local_y,current_state%column_local_x)  &
               + current_state%vgal)
          if (allocated(vprime_tot)) then
            vprime_tot(k) = vprime_tot(k) + &
                 ((current_state%v%data(k,current_state%column_local_y,current_state%column_local_x) &
                 - (current_state%global_grid%configuration%vertical%olvbar(k) - current_state%vgal))**2.)
          end if
          ww_tot(k) = ww_tot(k) + &
               (current_state%w%data(k,current_state%column_local_y,current_state%column_local_x)**2.)
          w_wind_tot(k) = w_wind_tot(k) + & 
               (current_state%w%data(k,current_state%column_local_y,current_state%column_local_x))
       enddo

!      <u'w'> and <v'w'> are on w-points, so we interpolate u and v both horizontally and vertically.
       do k=1, current_state%local_grid%size(Z_INDEX)-1
          uprime_w_local =  &
               0.25 * ( current_state%u%data(k,current_state%column_local_y,current_state%column_local_x)   + &
               current_state%u%data(k,current_state%column_local_y,current_state%column_local_x-1) + &
               current_state%u%data(k+1,current_state%column_local_y,current_state%column_local_x) + &
               current_state%u%data(k+1,current_state%column_local_y,current_state%column_local_x-1) ) + &
               current_state%ugal
          if (allocated(current_state%global_grid%configuration%vertical%olubar)) &
               uprime_w_local = uprime_w_local - &
               0.5  * ( current_state%global_grid%configuration%vertical%olubar(k) + &
               current_state%global_grid%configuration%vertical%olubar(k+1) )
          vprime_w_local = &
               0.25 * ( current_state%v%data(k,current_state%column_local_y,current_state%column_local_x)   + &
               current_state%v%data(k,current_state%column_local_y-1,current_state%column_local_x) + &
               current_state%v%data(k+1,current_state%column_local_y,current_state%column_local_x) + &
               current_state%v%data(k+1,current_state%column_local_y-1,current_state%column_local_x) ) + &
               current_state%vgal
          if (allocated(current_state%global_grid%configuration%vertical%olvbar)) &
               vprime_w_local = vprime_w_local - &
               0.5  * ( current_state%global_grid%configuration%vertical%olvbar(k) + &
               current_state%global_grid%configuration%vertical%olvbar(k+1) )
          uw_tot(k) = uw_tot(k) + uprime_w_local * &
               current_state%w%data(k,current_state%column_local_y,current_state%column_local_x)
          vw_tot(k) = vw_tot(k) + vprime_w_local * &
               current_state%w%data(k,current_state%column_local_y,current_state%column_local_x)
       enddo

       if (current_state%th%active) then
          do k=1, current_state%local_grid%size(Z_INDEX)
             theta_tot(k) = theta_tot(k) + & 
                  (current_state%th%data(k,current_state%column_local_y,current_state%column_local_x) &
                  + current_state%global_grid%configuration%vertical%thref(k))
             th2_tot(k) = th2_tot(k) + &
                  (current_state%th%data(k,current_state%column_local_y,current_state%column_local_x) - &
                  current_state%global_grid%configuration%vertical%olthbar(k) )**2
          enddo
!       <w'theta'> is on w-levels, so theta is interpolated to w-levels.
          do k=1, current_state%local_grid%size(Z_INDEX)-1
             thprime_w_local = 0.5 * (  &
                  (current_state%th%data(k,current_state%column_local_y,current_state%column_local_x) - &
                  current_state%global_grid%configuration%vertical%olthbar(k)) + &
                  (current_state%th%data(k+1,current_state%column_local_y,current_state%column_local_x) - &
                  current_state%global_grid%configuration%vertical%olthbar(k+1)) )
             wtheta_tot(k) = wtheta_tot(k) + &
                  current_state%w%data(k,current_state%column_local_y,current_state%column_local_x) * &
                  thprime_w_local
          enddo          
       endif
       if (current_state%th%active .and. .not. current_state%passive_q .and. current_state%number_q_fields .gt. 0) then
          do k=1, current_state%local_grid%size(Z_INDEX)
             qv_tot(k) = qv_tot(k) + (current_state%q(iqv)%data(k,current_state%column_local_y,current_state%column_local_x))   
             ql_tot(k) = ql_tot(k) + (current_state%q(iql)%data(k,current_state%column_local_y,current_state%column_local_x))
             ! temporary code for RH calculation
             exner = current_state%global_grid%configuration%vertical%rprefrcp(k)
             Pmb   = (current_state%global_grid%configuration%vertical%prefn(k)/100.)
             qv    = current_state%q(iqv)%data(k, current_state%column_local_y,current_state%column_local_x) 
             qc    = current_state%q(iql)%data(k, current_state%column_local_y,current_state%column_local_x)
             TdegK = (current_state%th%data(k,current_state%column_local_y,current_state%column_local_x) &
                  + current_state%global_grid%configuration%vertical%thref(k))*exner
             qs = qsaturation(TdegK, Pmb)
             rh_tot(k) = rh_tot(k) + (qv/qs)
          enddo
          do k=1, current_state%local_grid%size(Z_INDEX)-1
             qprime_w_local = 0.5 * (  &
                  (current_state%q(current_state%liquid_water_mixing_ratio_index)% &
                  data(k,current_state%column_local_y,current_state%column_local_x) - &
                  current_state%global_grid%configuration%vertical% &
                  olqbar(k,current_state%liquid_water_mixing_ratio_index)) + &
                  (current_state%q(current_state%liquid_water_mixing_ratio_index)% &
                  data(k+1,current_state%column_local_y,current_state%column_local_x) - &
                  current_state%global_grid%configuration%vertical% &
                  olqbar(k+1,current_state%liquid_water_mixing_ratio_index)) )
             wq_tot(k) = wq_tot(k) + qprime_w_local * &
                  current_state%w%data(k,current_state%column_local_y,current_state%column_local_x)
          end do
       endif
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
    if (name .eq. "theta_total_local") then
      field_information%enabled=current_state%th%active
    else if (name .eq. "vapour_mmr_total_local" .or. name .eq. "liquid_mmr_total_local") then
      field_information%enabled=.not. current_state%passive_q .and. current_state%number_q_fields .gt. 0
    else if (name .eq. "rh_total_local") then
      field_information%enabled=current_state%th%active .and. .not. current_state%passive_q .and. &
           current_state%number_q_fields .gt. 0
    else if (name .eq. "uu_total_local") then
      field_information%enabled=allocated(uprime_tot)
    else if (name .eq. "vv_total_local") then
      field_information%enabled=allocated(vprime_tot)
   else if (name .eq. "wtheta_total_local") then
      field_information%enabled=current_state%th%active
    else if (name .eq. "wq_total_local") then
      field_information%enabled= (.not.current_state%passive_q) .and. &
        (current_state%liquid_water_mixing_ratio_index > 0)
   else if (name .eq. "th2_total_local") then
      field_information%enabled=current_state%th%active
!   ========================================================================
!   2nd stream
    else if (name .eq. "i_theta_total_local") then
      field_information%enabled=current_state%th%active
    else if (name .eq. "i_vapour_mmr_total_local" .or. name .eq. "i_liquid_mmr_total_local") then
      field_information%enabled=.not. current_state%passive_q .and. current_state%number_q_fields .gt. 0
    else if (name .eq. "i_rh_total_local") then
      field_information%enabled=current_state%th%active .and. .not. current_state%passive_q .and. &
           current_state%number_q_fields .gt. 0
    else if (name .eq. "i_uu_total_local") then
      field_information%enabled=allocated(uprime_tot)
    else if (name .eq. "i_vv_total_local") then
      field_information%enabled=allocated(vprime_tot)
    else if (name .eq. "i_wtheta_total_local") then
      field_information%enabled=current_state%th%active
    else if (name .eq. "i_wq_total_local") then
      field_information%enabled= (.not.current_state%passive_q) .and. &
        (current_state%liquid_water_mixing_ratio_index > 0)
   else if (name .eq. "i_th2_total_local") then
      field_information%enabled=current_state%th%active
!   ========================================================================
    else 
      field_information%enabled=.true.
    end if

    if (name .eq. "prefn_local") then
       field_information%units = "Pa"
       field_information%long_name = "reference pressure at cell-faces"
    else if (name .eq. "u_wind_total_local") then
       field_information%units = "m/s"
       field_information%long_name = "per-MONC horizontal sum of u-wind"
    else if (name .eq. "v_wind_total_local") then
       field_information%units = "m/s"
       field_information%long_name = "per-MONC horizontal sum of v-wind"
    else if (name .eq. "w_wind_total_local") then
       field_information%units = "m/s"
       field_information%long_name = "per-MONC horizontal sum of w-wind"
    else if (name .eq. "theta_total_local") then
       field_information%units = "K"
       field_information%long_name = "per-MONC horizontal sum of potential temperature"
    else if (name .eq. "thinit_local") then
       field_information%units = "K"
       field_information%long_name = "initial vertical profile of potential temperature"
    endif
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
    
    if (name .eq. "prefn_local") then
       allocate(field_value%real_1d_array(current_state%local_grid%size(Z_INDEX)))
       do k = 1, current_state%local_grid%size(Z_INDEX)
          field_value%real_1d_array(k)=current_state%global_grid%configuration%vertical%prefn(k)
       enddo
    else if (name .eq. "rho_local") then
       allocate(field_value%real_1d_array(current_state%local_grid%size(Z_INDEX)))
       do k = 1, current_state%local_grid%size(Z_INDEX)
          field_value%real_1d_array(k)=current_state%global_grid%configuration%vertical%rho(k)
       enddo
    else if (name .eq. "rhon_local") then
       allocate(field_value%real_1d_array(current_state%local_grid%size(Z_INDEX)))
       do k = 1, current_state%local_grid%size(Z_INDEX)
          field_value%real_1d_array(k)=current_state%global_grid%configuration%vertical%rhon(k)
       enddo 
    else if (name .eq. "thref_local") then
       allocate(field_value%real_1d_array(current_state%local_grid%size(Z_INDEX)))
       do k = 1, current_state%local_grid%size(Z_INDEX)
          field_value%real_1d_array(k)=current_state%global_grid%configuration%vertical%thref(k)
       enddo   
    else if (name .eq. "thinit_local") then
       allocate(field_value%real_1d_array(current_state%local_grid%size(Z_INDEX)))
       do k = 1, current_state%local_grid%size(Z_INDEX)
          field_value%real_1d_array(k)=current_state%global_grid%configuration%vertical%theta_init(k)
       enddo    
    elseif (name .eq. "u_wind_total_local") then
       allocate(field_value%real_1d_array(current_state%local_grid%size(Z_INDEX)))
       do k = 1, current_state%local_grid%size(Z_INDEX)
          field_value%real_1d_array(k)=u_wind_tot(k)
       enddo
    else if (name .eq. "uu_total_local") then
       allocate(field_value%real_1d_array(current_state%local_grid%size(Z_INDEX)))
       do k = 1, current_state%local_grid%size(Z_INDEX)
          field_value%real_1d_array(k)=uprime_tot(k)
       enddo
    else if (name .eq. "v_wind_total_local") then
       allocate(field_value%real_1d_array(current_state%local_grid%size(Z_INDEX)))
       do k = 1, current_state%local_grid%size(Z_INDEX)
          field_value%real_1d_array(k)=v_wind_tot(k)
       enddo
    else if (name .eq. "vv_total_local") then
       allocate(field_value%real_1d_array(current_state%local_grid%size(Z_INDEX)))
       do k = 1, current_state%local_grid%size(Z_INDEX)
          field_value%real_1d_array(k)=vprime_tot(k)
       enddo
    else if (name .eq. "ww_total_local") then
       allocate(field_value%real_1d_array(current_state%local_grid%size(Z_INDEX)))
       do k = 1, current_state%local_grid%size(Z_INDEX)
          field_value%real_1d_array(k)=ww_tot(k)
       enddo 
    else if (name .eq. "theta_total_local") then
       allocate(field_value%real_1d_array(current_state%local_grid%size(Z_INDEX)))
       do k = 1, current_state%local_grid%size(Z_INDEX)
          field_value%real_1d_array(k)=theta_tot(k)
       enddo
    else if (name .eq. "vapour_mmr_total_local") then
       allocate(field_value%real_1d_array(current_state%local_grid%size(Z_INDEX)))
       do k = 1, current_state%local_grid%size(Z_INDEX)
          field_value%real_1d_array(k)=qv_tot(k)
       enddo
    else if (name .eq. "liquid_mmr_total_local") then
       allocate(field_value%real_1d_array(current_state%local_grid%size(Z_INDEX)))
       do k = 1, current_state%local_grid%size(Z_INDEX)
          field_value%real_1d_array(k)=ql_tot(k)
       enddo
    else if (name .eq. "w_wind_total_local") then
       allocate(field_value%real_1d_array(current_state%local_grid%size(Z_INDEX)))
       do k = 1, current_state%local_grid%size(Z_INDEX)
          field_value%real_1d_array(k)=w_wind_tot(k)
       enddo
    else if (name .eq. "rh_total_local") then
       allocate(field_value%real_1d_array(current_state%local_grid%size(Z_INDEX)))
       do k = 1, current_state%local_grid%size(Z_INDEX)
          field_value%real_1d_array(k)=rh_tot(k)
       enddo
 else if (name .eq. "wq_total_local") then
       allocate(field_value%real_1d_array(current_state%local_grid%size(Z_INDEX)))
       do k = 1, current_state%local_grid%size(Z_INDEX)
          field_value%real_1d_array(k)=wq_tot(k)
       enddo
    else if (name .eq. "wtheta_total_local") then
       allocate(field_value%real_1d_array(current_state%local_grid%size(Z_INDEX)))
       do k = 1, current_state%local_grid%size(Z_INDEX)
          field_value%real_1d_array(k)=wtheta_tot(k)
       enddo
    else if (name .eq. "uw_total_local") then
       allocate(field_value%real_1d_array(current_state%local_grid%size(Z_INDEX)))
       do k = 1, current_state%local_grid%size(Z_INDEX)
          field_value%real_1d_array(k)=uw_tot(k)
       enddo
    else if (name .eq. "vw_total_local") then
       allocate(field_value%real_1d_array(current_state%local_grid%size(Z_INDEX)))
       do k = 1, current_state%local_grid%size(Z_INDEX)
          field_value%real_1d_array(k)=vw_tot(k)
       enddo
    else if (name .eq. "th2_total_local") then
       allocate(field_value%real_1d_array(current_state%local_grid%size(Z_INDEX)))
       do k = 1, current_state%local_grid%size(Z_INDEX)
          field_value%real_1d_array(k)=th2_tot(k)
       enddo
! =====================================================
!   2nd stream
    else if (name .eq. "i_prefn_local") then
       allocate(field_value%real_1d_array(current_state%local_grid%size(Z_INDEX)))
       do k = 1, current_state%local_grid%size(Z_INDEX)
          field_value%real_1d_array(k)=current_state%global_grid%configuration%vertical%prefn(k)
       enddo
    else if (name .eq. "i_rho_local") then
       allocate(field_value%real_1d_array(current_state%local_grid%size(Z_INDEX)))
       do k = 1, current_state%local_grid%size(Z_INDEX)
          field_value%real_1d_array(k)=current_state%global_grid%configuration%vertical%rho(k)
       enddo
    else if (name .eq. "i_rhon_local") then
       allocate(field_value%real_1d_array(current_state%local_grid%size(Z_INDEX)))
       do k = 1, current_state%local_grid%size(Z_INDEX)
          field_value%real_1d_array(k)=current_state%global_grid%configuration%vertical%rhon(k)
       enddo 
    else if (name .eq. "i_thref_local") then
       allocate(field_value%real_1d_array(current_state%local_grid%size(Z_INDEX)))
       do k = 1, current_state%local_grid%size(Z_INDEX)
          field_value%real_1d_array(k)=current_state%global_grid%configuration%vertical%thref(k)
       enddo   
    else if (name .eq. "i_thinit_local") then
       allocate(field_value%real_1d_array(current_state%local_grid%size(Z_INDEX)))
       do k = 1, current_state%local_grid%size(Z_INDEX)
          field_value%real_1d_array(k)=current_state%global_grid%configuration%vertical%theta_init(k)
       enddo    
    else if (name .eq. "i_uinit_local") then
       allocate(field_value%real_1d_array(current_state%local_grid%size(Z_INDEX)))
       do k = 1, current_state%local_grid%size(Z_INDEX)
          field_value%real_1d_array(k)=current_state%global_grid%configuration%vertical%u_init(k)
       enddo    
    else if (name .eq. "i_vinit_local") then
       allocate(field_value%real_1d_array(current_state%local_grid%size(Z_INDEX)))
       do k = 1, current_state%local_grid%size(Z_INDEX)
          field_value%real_1d_array(k)=current_state%global_grid%configuration%vertical%v_init(k)
       enddo    
    elseif (name .eq. "i_u_wind_total_local") then
       allocate(field_value%real_1d_array(current_state%local_grid%size(Z_INDEX)))
       do k = 1, current_state%local_grid%size(Z_INDEX)
          field_value%real_1d_array(k)=u_wind_tot(k)
       enddo
    else if (name .eq. "i_uu_total_local") then
       allocate(field_value%real_1d_array(current_state%local_grid%size(Z_INDEX)))
       do k = 1, current_state%local_grid%size(Z_INDEX)
          field_value%real_1d_array(k)=uprime_tot(k)
       enddo
    else if (name .eq. "i_v_wind_total_local") then
       allocate(field_value%real_1d_array(current_state%local_grid%size(Z_INDEX)))
       do k = 1, current_state%local_grid%size(Z_INDEX)
          field_value%real_1d_array(k)=v_wind_tot(k)
       enddo
    else if (name .eq. "i_vv_total_local") then
       allocate(field_value%real_1d_array(current_state%local_grid%size(Z_INDEX)))
       do k = 1, current_state%local_grid%size(Z_INDEX)
          field_value%real_1d_array(k)=vprime_tot(k)
       enddo
    else if (name .eq. "i_ww_total_local") then
       allocate(field_value%real_1d_array(current_state%local_grid%size(Z_INDEX)))
       do k = 1, current_state%local_grid%size(Z_INDEX)
          field_value%real_1d_array(k)=ww_tot(k)
       enddo 
    else if (name .eq. "i_theta_total_local") then
       allocate(field_value%real_1d_array(current_state%local_grid%size(Z_INDEX)))
       do k = 1, current_state%local_grid%size(Z_INDEX)
          field_value%real_1d_array(k)=theta_tot(k)
       enddo
    else if (name .eq. "i_vapour_mmr_total_local") then
       allocate(field_value%real_1d_array(current_state%local_grid%size(Z_INDEX)))
       do k = 1, current_state%local_grid%size(Z_INDEX)
          field_value%real_1d_array(k)=qv_tot(k)
       enddo
    else if (name .eq. "i_liquid_mmr_total_local") then
       allocate(field_value%real_1d_array(current_state%local_grid%size(Z_INDEX)))
       do k = 1, current_state%local_grid%size(Z_INDEX)
          field_value%real_1d_array(k)=ql_tot(k)
       enddo
    else if (name .eq. "i_w_wind_total_local") then
       allocate(field_value%real_1d_array(current_state%local_grid%size(Z_INDEX)))
       do k = 1, current_state%local_grid%size(Z_INDEX)
          field_value%real_1d_array(k)=w_wind_tot(k)
       enddo
    else if (name .eq. "i_rh_total_local") then
       allocate(field_value%real_1d_array(current_state%local_grid%size(Z_INDEX)))
       do k = 1, current_state%local_grid%size(Z_INDEX)
          field_value%real_1d_array(k)=rh_tot(k)
       enddo
    else if (name .eq. "i_wq_total_local") then
       allocate(field_value%real_1d_array(current_state%local_grid%size(Z_INDEX)))
       do k = 1, current_state%local_grid%size(Z_INDEX)
          field_value%real_1d_array(k)=wq_tot(k)
       enddo
    else if (name .eq. "i_wtheta_total_local") then
       allocate(field_value%real_1d_array(current_state%local_grid%size(Z_INDEX)))
       do k = 1, current_state%local_grid%size(Z_INDEX)
          field_value%real_1d_array(k)=wtheta_tot(k)
       enddo
    else if (name .eq. "i_uw_total_local") then
       allocate(field_value%real_1d_array(current_state%local_grid%size(Z_INDEX)))
       do k = 1, current_state%local_grid%size(Z_INDEX)
          field_value%real_1d_array(k)=uw_tot(k)
       enddo
    else if (name .eq. "i_vw_total_local") then
       allocate(field_value%real_1d_array(current_state%local_grid%size(Z_INDEX)))
       do k = 1, current_state%local_grid%size(Z_INDEX)
          field_value%real_1d_array(k)=vw_tot(k)
       enddo
    else if (name .eq. "i_th2_total_local") then
       allocate(field_value%real_1d_array(current_state%local_grid%size(Z_INDEX)))
       do k = 1, current_state%local_grid%size(Z_INDEX)
          field_value%real_1d_array(k)=th2_tot(k)
       enddo
    end if
  end subroutine field_value_retrieval_callback
end module profile_diagnostics_mod
