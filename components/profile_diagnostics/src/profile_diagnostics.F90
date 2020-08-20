module profile_diagnostics_mod
  use monc_component_mod, only : COMPONENT_SCALAR_FIELD_TYPE, COMPONENT_DOUBLE_DATA_TYPE, &
       COMPONENT_ARRAY_FIELD_TYPE, COMPONENT_INTEGER_DATA_TYPE, &
       component_descriptor_type, component_field_value_type, component_field_information_type
  use grids_mod, only : vertical_grid_configuration_type, Z_INDEX, Y_INDEX, X_INDEX
  use optionsdatabase_mod, only : options_get_real, options_get_string, options_get_logical
  use state_mod, only : model_state_type
  use datadefn_mod, only : DEFAULT_PRECISION, STRING_LENGTH
  use q_indices_mod, only: get_q_index, standard_q_names
  use saturation_mod, only: qsaturation
  use logging_mod, only : LOG_ERROR, log_master_log  
  use def_tvd_diagnostic_terms, only: tvd_dgs_terms, allocate_tvd_diagnostic_terms
  use conversions_mod, only : conv_to_uppercase

  implicit none

#ifndef TEST_MODE
  private
#endif

  integer :: total_points, iqv=0, iql=0, iqr=0, iqi=0, iqs=0,    &
                           iqg=0
  real(kind=DEFAULT_PRECISION), dimension(:), allocatable ::     &
       tempfac, u_wind_tot, uprime_tot, v_wind_tot, vprime_tot,  &
       uprime, vprime, wke_tot, wwww_tot, www_tot, ww_tot,       &
       theta_tot, w_wind_tot, rh_tot, wtheta_ad_tot,             & 
       wtheta_cn_tot, uw_tot, vw_tot, uv_tot, th2_tot,           &
       thref, prefn, rho, rhon, thinit, uinit, vinit,            &
       ! mositure means
       q_temp, qv_tot, ql_tot, qr_tot, qi_tot, qs_tot, qg_tot,   &
       ! moisture flux terms
       wqv_cn_tot, wql_cn_tot, wqr_cn_tot, wqi_cn_tot,           &
       wqs_cn_tot, wqg_cn_tot,                                   &
       wqv_ad_tot, wql_ad_tot, wqr_ad_tot, wqi_ad_tot,           &
       wqs_ad_tot, wqg_ad_tot

  ! 3D and local profile total binary cloud masks
  character(len=STRING_LENGTH)                                :: cloud_mask_method
  real(kind=DEFAULT_PRECISION), dimension(:,:,:), allocatable :: cloud_mask
  real(kind=DEFAULT_PRECISION), dimension(:)    , allocatable :: cloud_mask_tot
  real(kind=DEFAULT_PRECISION), dimension(:)    , allocatable :: cloud_liq_mask_tot
  real(kind=DEFAULT_PRECISION), dimension(:)    , allocatable :: cloud_ice_mask_tot
  logical :: l_partial_liq_ice

  real(kind=DEFAULT_PRECISION) :: qlcrit, qicrit
  ! character string to determine the advection scheme
  character(len=5) :: advect_theta, advect_q, advect_flow

  type(vertical_grid_configuration_type) :: vertical_grid

  public profile_diagnostics_get_descriptor, tvd_dgs_terms

contains

  !> Provides the component descriptor for the core to register
  !! @returns The descriptor describing this component
  type(component_descriptor_type) function profile_diagnostics_get_descriptor()
    profile_diagnostics_get_descriptor%name="profile_diagnostics"
    profile_diagnostics_get_descriptor%version=0.1

    profile_diagnostics_get_descriptor%initialisation=>initialisation_callback
    profile_diagnostics_get_descriptor%timestep=>timestep_callback
   ! profile_diagnostics_get_descriptor%finalisation=>finalisation_callback

    profile_diagnostics_get_descriptor%field_value_retrieval=>field_value_retrieval_callback
    profile_diagnostics_get_descriptor%field_information_retrieval=>field_information_retrieval_callback
    allocate(profile_diagnostics_get_descriptor%published_fields(42+26+4))

    profile_diagnostics_get_descriptor%published_fields(1)="thref_local"
    profile_diagnostics_get_descriptor%published_fields(2)="prefn_local"
    profile_diagnostics_get_descriptor%published_fields(3)="rho_local"
    profile_diagnostics_get_descriptor%published_fields(4)="rhon_local" 
    profile_diagnostics_get_descriptor%published_fields(5)="thinit_local"
    profile_diagnostics_get_descriptor%published_fields(6)="uinit_local"
    profile_diagnostics_get_descriptor%published_fields(7)="vinit_local"
    profile_diagnostics_get_descriptor%published_fields(8)="theta_total_local"
    profile_diagnostics_get_descriptor%published_fields(9)="rh_total_local"
    profile_diagnostics_get_descriptor%published_fields(10)="u_wind_total_local"
    profile_diagnostics_get_descriptor%published_fields(11)="v_wind_total_local"
    profile_diagnostics_get_descriptor%published_fields(12)="w_wind_total_local"

    profile_diagnostics_get_descriptor%published_fields(13)="uu_total_local"
    profile_diagnostics_get_descriptor%published_fields(14)="vv_total_local"
    profile_diagnostics_get_descriptor%published_fields(15)="ww_total_local"
    profile_diagnostics_get_descriptor%published_fields(16)="uw_total_local"
    profile_diagnostics_get_descriptor%published_fields(17)="vw_total_local"
    profile_diagnostics_get_descriptor%published_fields(18)="uv_total_local"
    profile_diagnostics_get_descriptor%published_fields(19)="wke_total_local"
    profile_diagnostics_get_descriptor%published_fields(20)="wwww_total_local"
    profile_diagnostics_get_descriptor%published_fields(21)="www_total_local"
    profile_diagnostics_get_descriptor%published_fields(22)="wtheta_ad_total_local"
    profile_diagnostics_get_descriptor%published_fields(23)="wtheta_cn_total_local"
    profile_diagnostics_get_descriptor%published_fields(24)="th2_total_local"
    
    profile_diagnostics_get_descriptor%published_fields(25)="vapour_mmr_total_local"
    profile_diagnostics_get_descriptor%published_fields(26)="liquid_mmr_total_local"
    profile_diagnostics_get_descriptor%published_fields(27)="rain_mmr_total_local"
    profile_diagnostics_get_descriptor%published_fields(28)="ice_mmr_total_local"
    profile_diagnostics_get_descriptor%published_fields(29)="snow_mmr_total_local"
    profile_diagnostics_get_descriptor%published_fields(30)="graupel_mmr_total_local"

    profile_diagnostics_get_descriptor%published_fields(31)="wqv_ad_total_local"
    profile_diagnostics_get_descriptor%published_fields(32)="wqv_cn_total_local"
    profile_diagnostics_get_descriptor%published_fields(33)="wql_ad_total_local"
    profile_diagnostics_get_descriptor%published_fields(34)="wql_cn_total_local"
    profile_diagnostics_get_descriptor%published_fields(35)="wqr_ad_total_local"
    profile_diagnostics_get_descriptor%published_fields(36)="wqr_cn_total_local"
    profile_diagnostics_get_descriptor%published_fields(37)="wqi_ad_total_local"
    profile_diagnostics_get_descriptor%published_fields(38)="wqi_cn_total_local"
    profile_diagnostics_get_descriptor%published_fields(39)="wqs_ad_total_local"
    profile_diagnostics_get_descriptor%published_fields(40)="wqs_cn_total_local"
    profile_diagnostics_get_descriptor%published_fields(41)="wqg_ad_total_local"
    profile_diagnostics_get_descriptor%published_fields(42)="wqg_cn_total_local"

!   =====================================================
!   2nd, provisionally instantaneous, stream
    profile_diagnostics_get_descriptor%published_fields(42+1)="i_theta_total_local"
    profile_diagnostics_get_descriptor%published_fields(42+2)="i_vapour_mmr_total_local"
    profile_diagnostics_get_descriptor%published_fields(42+3)="i_liquid_mmr_total_local"
    profile_diagnostics_get_descriptor%published_fields(42+4)="i_u_wind_total_local"
    profile_diagnostics_get_descriptor%published_fields(42+5)="i_uu_total_local"
    profile_diagnostics_get_descriptor%published_fields(42+6)="i_v_wind_total_local"
    profile_diagnostics_get_descriptor%published_fields(42+7)="i_vv_total_local"
    profile_diagnostics_get_descriptor%published_fields(42+8)="i_ww_total_local"
    profile_diagnostics_get_descriptor%published_fields(42+9)="i_w_wind_total_local"
    profile_diagnostics_get_descriptor%published_fields(42+10)="i_rh_total_local"
    profile_diagnostics_get_descriptor%published_fields(42+11)="i_thref_local"
    profile_diagnostics_get_descriptor%published_fields(42+12)="i_prefn_local"
    profile_diagnostics_get_descriptor%published_fields(42+13)="i_rho_local"
    profile_diagnostics_get_descriptor%published_fields(42+14)="i_rhon_local" 
    profile_diagnostics_get_descriptor%published_fields(42+15)="i_thinit_local"

    profile_diagnostics_get_descriptor%published_fields(42+16)="i_uw_total_local"
    profile_diagnostics_get_descriptor%published_fields(42+17)="i_vw_total_local"
    profile_diagnostics_get_descriptor%published_fields(42+18)="i_wtheta_cn_total_local"
    profile_diagnostics_get_descriptor%published_fields(42+19)="i_th2_total_local"
   
    profile_diagnostics_get_descriptor%published_fields(42+20)="i_uinit_local"
    profile_diagnostics_get_descriptor%published_fields(42+21)="i_vinit_local"

    profile_diagnostics_get_descriptor%published_fields(42+22)="i_rain_mmr_total_local"
    profile_diagnostics_get_descriptor%published_fields(42+23)="i_ice_mmr_total_local"
    profile_diagnostics_get_descriptor%published_fields(42+24)="i_snow_mmr_total_local"
    profile_diagnostics_get_descriptor%published_fields(42+25)="i_graupel_mmr_total_local"

!   =====================================================
!   cloud mask fields
    profile_diagnostics_get_descriptor%published_fields(42+25+1)="cloud_mask"
    profile_diagnostics_get_descriptor%published_fields(42+25+2)="cloud_mask_total_local"
    profile_diagnostics_get_descriptor%published_fields(42+25+3)="cloud_liq_mask_total_local"
    profile_diagnostics_get_descriptor%published_fields(42+25+4)="cloud_ice_mask_total_local"


  end function profile_diagnostics_get_descriptor

  subroutine initialisation_callback(current_state)
    type(model_state_type), target, intent(inout) :: current_state

    integer :: k

    vertical_grid=current_state%global_grid%configuration%vertical

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
    allocate( uw_tot(current_state%local_grid%size(Z_INDEX))     &
         ,   vw_tot(current_state%local_grid%size(Z_INDEX))     &
         ,   uv_tot(current_state%local_grid%size(Z_INDEX))     &
         ,   wwww_tot(current_state%local_grid%size(Z_INDEX))   &
         ,   www_tot(current_state%local_grid%size(Z_INDEX)) )
    
    if (allocated(current_state%global_grid%configuration%vertical%olubar)) then
       allocate(uprime_tot(current_state%local_grid%size(Z_INDEX)), &
            uprime(current_state%local_grid%size(Z_INDEX)))
    end if
    
    if (allocated(current_state%global_grid%configuration%vertical%olvbar)) then
       allocate(vprime_tot(current_state%local_grid%size(Z_INDEX)), &
            vprime(current_state%local_grid%size(Z_INDEX)))
    end if

    if (allocated(vprime) .and. allocated(uprime)) then
       allocate(wke_tot(current_state%local_grid%size(Z_INDEX)))
    endif

    if (current_state%th%active) then
       allocate(theta_tot(current_state%local_grid%size(Z_INDEX)) &
            , thref(current_state%local_grid%size(Z_INDEX))       &
            , thinit(current_state%local_grid%size(Z_INDEX))      &
            , wtheta_ad_tot(current_state%local_grid%size(Z_INDEX)) &
            , wtheta_cn_tot(current_state%local_grid%size(Z_INDEX)) & 
            , th2_tot(current_state%local_grid%size(Z_INDEX)) )
    endif
    
    ! determine advection scheme used for mom and scalars
    advect_flow = options_get_string(current_state%options_database, "advection_flow_fields")
    advect_theta = options_get_string(current_state%options_database, "advection_theta_field")
    advect_q = options_get_string(current_state%options_database, "advection_q_fields")

    if (.not. current_state%passive_q .and. current_state%number_q_fields .gt. 0) then 
       iqv=get_q_index(standard_q_names%VAPOUR, 'profile_diags')                         
       iql=get_q_index(standard_q_names%CLOUD_LIQUID_MASS, 'profile_diags')  
       qlcrit=options_get_real(current_state%options_database, "qlcrit")  
       qicrit=options_get_real(current_state%options_database, "qicrit")
 
       allocate(qv_tot(current_state%local_grid%size(Z_INDEX))  &
            , ql_tot(current_state%local_grid%size(Z_INDEX)),   &
            q_temp(current_state%local_grid%size(Z_INDEX)),   & 
            wqv_cn_tot(current_state%local_grid%size(Z_INDEX)), &
            wqv_ad_tot(current_state%local_grid%size(Z_INDEX)), &
            wql_cn_tot(current_state%local_grid%size(Z_INDEX)), &
            wql_ad_tot(current_state%local_grid%size(Z_INDEX)))
       if (current_state%th%active) &
            allocate(rh_tot(current_state%local_grid%size(Z_INDEX)))

       ! allocate other hydrometeors. Allocation dependent on index being set in
       ! appropriate microphysics scheme (see casim component from example)
       if (current_state%rain_water_mixing_ratio_index > 0) then
          iqr = current_state%rain_water_mixing_ratio_index
          allocate(qr_tot(current_state%local_grid%size(Z_INDEX)))
          allocate(wqr_cn_tot(current_state%local_grid%size(Z_INDEX)))
          allocate(wqr_ad_tot(current_state%local_grid%size(Z_INDEX)))
       endif
       if (current_state%ice_water_mixing_ratio_index > 0) then
          iqi = current_state%ice_water_mixing_ratio_index 
          allocate(qi_tot(current_state%local_grid%size(Z_INDEX)))
          allocate(wqi_cn_tot(current_state%local_grid%size(Z_INDEX))) 
          allocate(wqi_ad_tot(current_state%local_grid%size(Z_INDEX))) 
       endif
       if (current_state%snow_water_mixing_ratio_index > 0) then
          iqs = current_state%snow_water_mixing_ratio_index
          allocate(qs_tot(current_state%local_grid%size(Z_INDEX)))
          allocate(wqs_cn_tot(current_state%local_grid%size(Z_INDEX)))
          allocate(wqs_ad_tot(current_state%local_grid%size(Z_INDEX)))
       endif
       if (current_state%graupel_water_mixing_ratio_index > 0) then
          iqg = current_state%graupel_water_mixing_ratio_index
          allocate(qg_tot(current_state%local_grid%size(Z_INDEX)))
          allocate(wqg_cn_tot(current_state%local_grid%size(Z_INDEX)))
          allocate(wqg_ad_tot(current_state%local_grid%size(Z_INDEX)))
       endif

       ! arrange and allocate cloud fraction diagnostics...3d mask is optional
       cloud_mask_method = conv_to_uppercase(                                      &
           options_get_string(current_state%options_database, "cloud_mask_method"))
       if (.not. (cloud_mask_method == "DEFAULT"   .or.                      &
                  cloud_mask_method == "SOCRATES"      ) )  then
         call log_master_log(LOG_ERROR,                                            &
          "Requested cloud_mask_method is invalid.  Check profile_diagnostics.F90") 
       end if ! cloud_mask_method validity check
       if (options_get_logical(current_state%options_database, "l_cloud_mask"))    &
         allocate(cloud_mask(current_state%local_grid%size(Z_INDEX),               &
                             current_state%local_grid%size(Y_INDEX),               &
                             current_state%local_grid%size(X_INDEX)))
       allocate(cloud_mask_tot(current_state%local_grid%size(Z_INDEX)))
       allocate(cloud_liq_mask_tot(current_state%local_grid%size(Z_INDEX)))
       allocate(cloud_ice_mask_tot(current_state%local_grid%size(Z_INDEX)))
       l_partial_liq_ice =                                                           &
           options_get_logical(current_state%options_database, "l_partial_liq_ice")
    endif

    call allocate_tvd_diagnostic_terms(current_state, tvd_dgs_terms)

  end subroutine initialisation_callback  

  subroutine timestep_callback(current_state)
    type(model_state_type), target, intent(inout) :: current_state

    integer :: k, i, iq_tmp, icol, jcol
    real(kind=DEFAULT_PRECISION) :: cltop_col, clbas_col, qv, qc, TdegK, Pmb &
         , qs, exner
    real(kind=DEFAULT_PRECISION) :: uprime_w_local, vprime_w_local &
         , thprime_w_local, qprime_w_local 

    
    icol=current_state%column_local_x
    jcol=current_state%column_local_y

    if (current_state%first_timestep_column) then
       u_wind_tot(:) = 0.0_DEFAULT_PRECISION
       if (allocated(uprime_tot)) uprime_tot(:) = 0.0_DEFAULT_PRECISION
       if (allocated(uprime)) uprime(:) = 0.0_DEFAULT_PRECISION
       v_wind_tot(:) = 0.0_DEFAULT_PRECISION
       if (allocated(vprime_tot)) vprime_tot(:) = 0.0_DEFAULT_PRECISION
       if (allocated(vprime)) vprime(:) = 0.0_DEFAULT_PRECISION
       if (allocated(wke_tot)) wke_tot(:) =  0.0_DEFAULT_PRECISION
       
       w_wind_tot(:) = 0.0_DEFAULT_PRECISION
       ww_tot(:) = 0.0_DEFAULT_PRECISION
       wwww_tot(:) = 0.0_DEFAULT_PRECISION
       www_tot(:) = 0.0_DEFAULT_PRECISION
       uw_tot(:)     = 0.0_DEFAULT_PRECISION
       vw_tot(:)     = 0.0_DEFAULT_PRECISION
       uv_tot(:)     = 0.0_DEFAULT_PRECISION
       
       if (current_state%th%active) then
          wtheta_ad_tot(:) = 0.0_DEFAULT_PRECISION
          wtheta_cn_tot(:) = 0.0_DEFAULT_PRECISION
          th2_tot(:)    = 0.0_DEFAULT_PRECISION
          theta_tot(:)=0.0_DEFAULT_PRECISION
       endif
       if (.not. current_state%passive_q .and. &
            current_state%number_q_fields .gt. 0) then 
          q_temp(:)=0.0_DEFAULT_PRECISION
          qv_tot(:)=0.0_DEFAULT_PRECISION
          ql_tot(:)=0.0_DEFAULT_PRECISION
          wqv_cn_tot(:)=0.0_DEFAULT_PRECISION
          wqv_ad_tot(:)=0.0_DEFAULT_PRECISION
          wql_cn_tot(:)=0.0_DEFAULT_PRECISION
          wql_ad_tot(:)=0.0_DEFAULT_PRECISION 
          if (current_state%th%active) &
               rh_tot(:) = 0.0_DEFAULT_PRECISION
          if (iqr > 0) then 
             qr_tot(:) = 0.0_DEFAULT_PRECISION
             wqr_cn_tot(:) = 0.0_DEFAULT_PRECISION
             wqr_ad_tot(:) = 0.0_DEFAULT_PRECISION
          endif
          if (iqi > 0) then
             qi_tot(:) = 0.0_DEFAULT_PRECISION
             wqi_cn_tot(:) = 0.0_DEFAULT_PRECISION
             wqi_ad_tot(:) = 0.0_DEFAULT_PRECISION
          endif
          if (iqs > 0) then
             qs_tot(:) = 0.0_DEFAULT_PRECISION
             wqs_cn_tot(:) = 0.0_DEFAULT_PRECISION
             wqs_ad_tot(:) = 0.0_DEFAULT_PRECISION
          endif
          if (iqg > 0) then 
             qg_tot(:) = 0.0_DEFAULT_PRECISION
             wqg_cn_tot(:) = 0.0_DEFAULT_PRECISION
             wqg_ad_tot(:) = 0.0_DEFAULT_PRECISION
          endif

          if (allocated(cloud_mask)) cloud_mask(:,:,:) = 0.0_DEFAULT_PRECISION
          cloud_mask_tot(:) = 0.0_DEFAULT_PRECISION
          cloud_liq_mask_tot(:) = 0.0_DEFAULT_PRECISION
          cloud_ice_mask_tot(:) = 0.0_DEFAULT_PRECISION
       endif
    end if
    !

    if (.not. current_state%halo_column) then
    ! work out the sum of u and v wind over local domain
    do k=1, current_state%local_grid%size(Z_INDEX)
       u_wind_tot(k) = u_wind_tot(k) + & 
            (current_state%u%data(k,jcol,icol)  &
            + current_state%ugal)
       if (allocated(uprime_tot)) then
          uprime(k) =  &
               (current_state%u%data(k, jcol, icol) &
               - (current_state%global_grid%configuration%vertical%olubar(k) - current_state%ugal))
          !uprime(k) =  &
          !     ((current_state%u%data(k,jcol,icol) &
          !    - (current_state%global_grid%configuration%vertical%olubar(k) - current_state%ugal)))
          uprime_tot(k) = uprime_tot(k) + uprime(k)**2.0_DEFAULT_precision
       end if
       v_wind_tot(k) = v_wind_tot(k) + & 
            (current_state%v%data(k,jcol,icol)  &
            + current_state%vgal)
       if (allocated(vprime_tot)) then
          vprime(k) = &
          (current_state%v%data(k, jcol, icol) &
               - (current_state%global_grid%configuration%vertical%olvbar(k) - current_state%vgal))
          !vprime(k) =  &
          !     ((current_state%v%data(k,jcol,icol) &
          !     - (current_state%global_grid%configuration%vertical%olvbar(k) - current_state%vgal)))
          vprime_tot(k) = vprime_tot(k) + vprime(k)**2.0_DEFAULT_precision
       end if
    enddo
    ! Note: If TVD advection used, a conversion for the grid (from the LEM (RESDGS.653,654))
    !       is used but it is unclear whether that is correct. Also, it is worth noting that 
    !       www and wwww have no offset relating to the advection scheme - AH 21/09/2017
    if (trim(advect_flow) .eq. "pw") then 
       ww_tot(:) = ww_tot(:) + &
            (current_state%w%data(:,jcol,icol)**2.)
    else if (trim(advect_flow) .eq. "tvd") then
       do k=2, current_state%local_grid%size(Z_INDEX)
          ww_tot(k) = ww_tot(k) + 0.5 * &
               ((current_state%w%data(k-1,jcol,icol) + current_state%w%data(k,jcol,icol)) * &
               tvd_dgs_terms%adv_w_dgs(k,jcol,icol))
       enddo
    endif
    
    do k=2, current_state%local_grid%size(Z_INDEX)       
       www_tot(k) = www_tot(k) + &
            (current_state%w%data(k,jcol,icol)**3.)
       wwww_tot(k) = wwww_tot(k) + &
            (current_state%w%data(k,jcol,icol)**4.)
       w_wind_tot(k) = w_wind_tot(k) + & 
            (current_state%w%data(k,jcol,icol))
    enddo

!      <u'w'> and <v'w'> are on w-points, so we interpolate u and v both horizontally and vertically.
!      NOTE: all "prime_w" values are the prognostic interpolated on to the w levels
    if (trim(advect_flow) .eq. "pw") then 
       do k=1, current_state%local_grid%size(Z_INDEX)-1
          uprime_w_local =  &
               0.25 * ( current_state%u%data(k,jcol,icol)   + &
               current_state%u%data(k,jcol,icol-1) + &
               current_state%u%data(k+1,jcol,icol) + &
               current_state%u%data(k+1,jcol,icol-1) ) + &
               current_state%ugal
          if (allocated(current_state%global_grid%configuration%vertical%olubar)) &
               uprime_w_local = uprime_w_local - &
               0.5  * ( current_state%global_grid%configuration%vertical%olubar(k) + &
                  current_state%global_grid%configuration%vertical%olubar(k+1) )
          uw_tot(k) = uw_tot(k) + uprime_w_local * &
               current_state%w%data(k,jcol,icol)
       enddo
    else if (trim(advect_flow) .eq. "tvd") then
       do k=2, current_state%local_grid%size(Z_INDEX)-1
          uw_tot(k) = uw_tot(k) + 0.5 * (      &
               (current_state%w%data(k,jcol,icol) + current_state%w%data(k,jcol,icol+1)) * &
               tvd_dgs_terms%adv_u_dgs(k+1,jcol,icol))
       enddo
    endif
    
    if (trim(advect_flow) .eq. "pw") then 
       do k=1, current_state%local_grid%size(Z_INDEX)-1
          vprime_w_local = &
               0.25 * ( current_state%v%data(k,jcol,icol)   + &
               current_state%v%data(k,jcol-1,icol) + &
               current_state%v%data(k+1,jcol,icol) + &
               current_state%v%data(k+1,jcol-1,icol) ) + &
               current_state%vgal
          if (allocated(current_state%global_grid%configuration%vertical%olvbar)) &
                  vprime_w_local = vprime_w_local - &
                  0.5  * ( current_state%global_grid%configuration%vertical%olvbar(k) + &
                  current_state%global_grid%configuration%vertical%olvbar(k+1) )             
          vw_tot(k) = vw_tot(k) + vprime_w_local * &
               current_state%w%data(k,jcol,icol)
       enddo
    else if (trim(advect_flow) .eq. "tvd") then
       do k=2, current_state%local_grid%size(Z_INDEX)-1
          vw_tot(k) = vw_tot(k) + 0.5 * (      &
               (current_state%w%data(k,jcol,icol) + current_state%w%data(k,jcol+1,icol)) * &
                  tvd_dgs_terms%adv_v_dgs(k+1,jcol,icol))
       enddo
    endif
    
    if (allocated(current_state%global_grid%configuration%vertical%olvbar) .and.               &
         allocated(current_state%global_grid%configuration%vertical%olubar)) then
       ! LEM equivalent code RESDGS.784, 785, comment from LEM
       ! No attempt has been made to provide a calculation of this term which
       ! is consistent with TVD advection 
       do k=1, current_state%local_grid%size(Z_INDEX)-1
          uv_tot(k) = uv_tot(k) + 0.25 * uprime(k) * (vprime(k)          &
               + (current_state%v%data(k,jcol,icol+1))   &
               ! v primed term (not squared) at y - 1 , do whole calc for vprime at y-1
               + (current_state%v%data(k,jcol-1,icol)   &
               - (current_state%global_grid%configuration%vertical%olvbar(k) - current_state%vgal))   &
               + (current_state%v%data(k,jcol-1,icol+1)))
       enddo
    endif
       
    if (allocated(wke_tot)) then
       do k=1, current_state%local_grid%size(Z_INDEX)-1
          uprime_w_local =  &
               0.25_DEFAULT_PRECISION * ( current_state%u%data(k,jcol,icol)   + &
               current_state%u%data(k,jcol,icol-1) + &
               current_state%u%data(k+1,jcol,icol) + &
               current_state%u%data(k+1,jcol,icol-1) ) + &
               current_state%ugal
          if (allocated(current_state%global_grid%configuration%vertical%olubar)) &
               uprime_w_local = uprime_w_local - &
               0.5_DEFAULT_PRECISION * ( current_state%global_grid%configuration%vertical%olubar(k) + &
                  current_state%global_grid%configuration%vertical%olubar(k+1) )
          vprime_w_local = &
               0.25_DEFAULT_PRECISION * ( current_state%v%data(k,jcol,icol)   + &
               current_state%v%data(k,jcol-1,icol) + &
               current_state%v%data(k+1,jcol,icol) + &
               current_state%v%data(k+1,jcol-1,icol) ) + &
               current_state%vgal
          if (allocated(current_state%global_grid%configuration%vertical%olvbar)) &
               vprime_w_local = vprime_w_local - &
               0.5_DEFAULT_PRECISION * ( current_state%global_grid%configuration%vertical%olvbar(k) + &
               current_state%global_grid%configuration%vertical%olvbar(k+1) )

          wke_tot(k) = wke_tot(k) + 0.5_DEFAULT_PRECISION *                     &
               current_state%global_grid%configuration%vertical%rhon(k) *       &
               current_state%w%data(k,jcol,icol) *                              &
               ( uprime_w_local * uprime_w_local +                              &
                 vprime_w_local * vprime_w_local +                              &
                 current_state%w%data(k,jcol,icol) * current_state%w%data(k,jcol,icol) )
       enddo
    endif
       
    if (current_state%th%active) then
       do k=1, current_state%local_grid%size(Z_INDEX)
          theta_tot(k) = theta_tot(k) + & 
               (current_state%th%data(k,jcol,icol) &
               + current_state%global_grid%configuration%vertical%thref(k))
          th2_tot(k) = th2_tot(k) + &
               (current_state%th%data(k,jcol,icol) - &
               current_state%global_grid%configuration%vertical%olthbar(k) )**2
       enddo
       !       <w'theta'> is on w-levels, so theta is interpolated to w-levels.
       !       Set wtheta_ad_tot to the actual model TH flux.
!       Set wtheta_cn_tot slot to flux seen by w times the w-momentum equation. 
       !       This is the flux relevant for the tke balance. 
       !       The two are equal if the centred scheme is used but differ if TVD advection is used on scalars). 
       if (trim(advect_theta) .eq. "pw") then 
          do k=1, current_state%local_grid%size(Z_INDEX)-1
             thprime_w_local = 0.5 * (                       &
                  current_state%th%data(k,jcol,icol) +       &
                  current_state%th%data(k+1,jcol,icol) )
             wtheta_ad_tot(k) = wtheta_ad_tot(k) +           &
                  (current_state%w%data(k,jcol,icol) *        &
                  thprime_w_local)
             wtheta_cn_tot(k) = wtheta_cn_tot(k)
          enddo
       else if (trim(advect_theta) .eq. "tvd") then
          do k=2, current_state%local_grid%size(Z_INDEX)-1
             thprime_w_local = 0.5 * (                        &
                  current_state%th%data(k,jcol,icol) +       &
                  current_state%th%data(k+1,jcol,icol) )
             wtheta_cn_tot(k) = wtheta_cn_tot(k) + &
                  current_state%w%data(k,jcol,icol) * &
                  thprime_w_local
             wtheta_ad_tot(k) = wtheta_ad_tot(k) +           &
                  current_state%w%data(k,jcol,icol) *         &
                  tvd_dgs_terms%adv_th_dgs(k+1,jcol,icol)
          enddo
       endif
    endif
    
    if (.not. current_state%passive_q .and. current_state%number_q_fields .gt. 0) then
       qv_tot(:) = qv_tot(:) + (current_state%q(iqv)%data(:,jcol,icol))   
       ql_tot(:) = ql_tot(:) + (current_state%q(iql)%data(:,jcol,icol))
       if (current_state%th%active) then 
          ! temporary code for RH calculation
          do k=1, current_state%local_grid%size(Z_INDEX)
             exner = current_state%global_grid%configuration%vertical%rprefrcp(k)
             Pmb   = (current_state%global_grid%configuration%vertical%prefn(k)/100.)
             qv    = current_state%q(iqv)%data(k, jcol,icol) 
             qc    = current_state%q(iql)%data(k, jcol,icol)
             TdegK = (current_state%th%data(k,jcol,icol) &
                  + current_state%global_grid%configuration%vertical%thref(k))*exner
             qs = qsaturation(TdegK, Pmb)
             rh_tot(k) = rh_tot(k) + (qv/qs) 
          enddo
       endif
       ! hydrometeor mass profiles
       if (iqr > 0) &
            qr_tot(:) = qr_tot(:) + (current_state%q(iqr)%data(:,jcol,icol))
       if (iqi > 0) &
            qi_tot(:) = qi_tot(:) + (current_state%q(iqi)%data(:,jcol,icol))
       if (iqs > 0) &
            qs_tot(:) = qs_tot(:) + (current_state%q(iqs)%data(:,jcol,icol))
       if (iqg > 0) &
            qg_tot(:) = qg_tot(:) + (current_state%q(iqg)%data(:,jcol,icol))
       !
       ! moisture field fluxes
       ! vapour
       call calculate_wq(current_state, jcol, icol, iqv, wqv_cn_tot, wqv_ad_tot, advect_q)
       ! cloud liquid
       call calculate_wq(current_state, jcol, icol, iql , wql_cn_tot, wql_ad_tot, advect_q) 
       ! rain mass
       if (iqr > 0) &
            call calculate_wq(current_state, jcol, icol, iqr, wqr_cn_tot, wqr_ad_tot, advect_q)
       ! ice mass
       if (iqi > 0) &
            call calculate_wq(current_state, jcol, icol, iqi, wqi_cn_tot, wqi_ad_tot, advect_q)
       ! snow mass
       if (iqs > 0) & 
            call calculate_wq(current_state, jcol, icol, iqs, wqs_cn_tot, wqs_ad_tot, advect_q)
          ! graupel mass
       if (iqg > 0) &
            call calculate_wq(current_state, jcol, icol, iqg, wqg_cn_tot, wqg_ad_tot, advect_q)

       ! cloud mask / cloud fraction
       call calculate_cloud_mask(current_state, jcol, icol)
    endif ! end q_passive and number q field test
   endif
  end subroutine timestep_callback

!  !> Frees up the memory associated with the advection
!  !! @param current_state The current model state
!  subroutine finalisation_callback(current_state)
!    type(model_state_type), target, intent(inout) :: current_state
!
!    call deallocate_tvd_diagnostic_terms(current_state, tvd_dgs_terms)
!    
!  end subroutine finalisation_callback

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
    else if (name .eq. "vapour_mmr_total_local" .or. name .eq. "liquid_mmr_total_local" &
         .or. name .eq. "wqv_ad_total_local" .or. name .eq. "wqv_cn_total_local"        &
         .or. name .eq. "wql_ad_total_local" .or. name .eq. "wql_cn_total_local" ) then
       field_information%enabled=.not. current_state%passive_q .and. current_state%number_q_fields .gt. 0
    else if (name .eq. "rain_mmr_total_local" .or. name .eq. "wqr_ad_total_local" &
         .or. name .eq. "wqr_cn_total_local" ) then
       field_information%enabled= current_state%rain_water_mixing_ratio_index .gt. 0
    else if (name .eq. "ice_mmr_total_local"  .or. name .eq. "wqi_ad_total_local" &
         .or. name .eq. "wqi_cn_total_local" ) then
       field_information%enabled= current_state%ice_water_mixing_ratio_index .gt. 0
    else if (name .eq. "snow_mmr_total_local" .or. name .eq. "wqs_ad_total_local" &
         .or. name .eq. "wqs_cn_total_local"  ) then
       field_information%enabled= current_state%snow_water_mixing_ratio_index .gt. 0
    else if (name .eq. "graupel_mmr_total_local" .or. name .eq. "wqg_ad_total_local" &
         .or. name .eq. "wqg_cn_total_local"  ) then
       field_information%enabled= current_state%graupel_water_mixing_ratio_index .gt. 0
    else if (name .eq. "rh_total_local") then
      field_information%enabled=current_state%th%active .and. .not. current_state%passive_q .and. &
           current_state%number_q_fields .gt. 0
    else if (name .eq. "uu_total_local") then
      field_information%enabled=allocated(uprime_tot)
    else if (name .eq. "vv_total_local") then
      field_information%enabled=allocated(vprime_tot)
   else if (name .eq. "wtheta_ad_total_local" .or. name .eq. "wtheta_cn_total_local") then
      field_information%enabled=current_state%th%active
   else if (name .eq. "th2_total_local") then
      field_information%enabled=current_state%th%active
!   ========================================================================
!   2nd stream
    else if (name .eq. "i_theta_total_local") then
      field_information%enabled=current_state%th%active
    else if (name .eq. "i_vapour_mmr_total_local" .or. name .eq. "i_liquid_mmr_total_local") then
       field_information%enabled=.not. current_state%passive_q .and. current_state%number_q_fields .gt. 0
    else if (name .eq. "i_rain_mmr_total_local" ) then
       field_information%enabled= current_state%rain_water_mixing_ratio_index .gt. 0
    else if (name .eq. "i_ice_mmr_total_local" ) then
       field_information%enabled= current_state%ice_water_mixing_ratio_index .gt. 0
    else if (name .eq. "i_snow_mmr_total_local" ) then
       field_information%enabled= current_state%snow_water_mixing_ratio_index .gt. 0
    else if (name .eq. "i_graupel_mmr_total_local" ) then
       field_information%enabled= current_state%graupel_water_mixing_ratio_index .gt. 0
    else if (name .eq. "i_rh_total_local") then
      field_information%enabled=current_state%th%active .and. .not. current_state%passive_q .and. &
           current_state%number_q_fields .gt. 0
    else if (name .eq. "i_uu_total_local") then
      field_information%enabled=allocated(uprime_tot)
    else if (name .eq. "i_vv_total_local") then
      field_information%enabled=allocated(vprime_tot)
    else if (name .eq. "i_wtheta_cn_total_local") then
      field_information%enabled=current_state%th%active
    else if (name .eq. "i_th2_total_local") then
      field_information%enabled=current_state%th%active
!   ========================================================================
    else if (name .eq. "cloud_mask") then
      field_information%number_dimensions=3
      field_information%dimension_sizes(1)=current_state%local_grid%size(Z_INDEX)
      field_information%dimension_sizes(2)=current_state%local_grid%size(Y_INDEX)
      field_information%dimension_sizes(3)=current_state%local_grid%size(X_INDEX)
      field_information%enabled=allocated(cloud_mask)
    else if (name .eq. "cloud_mask_total_local") then
      field_information%enabled=allocated(cloud_mask_tot)
    else if (name .eq. "cloud_liq_mask_total_local") then
      field_information%enabled=allocated(cloud_liq_mask_tot)
    else if (name .eq. "cloud_ice_mask_total_local") then
      field_information%enabled=allocated(cloud_ice_mask_tot)
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
    else if (name .eq. "www_total_local") then
       allocate(field_value%real_1d_array(current_state%local_grid%size(Z_INDEX)))
       do k = 1, current_state%local_grid%size(Z_INDEX)
          field_value%real_1d_array(k)=www_tot(k)
       enddo
    else if (name .eq. "wwww_total_local") then
       allocate(field_value%real_1d_array(current_state%local_grid%size(Z_INDEX)))
       do k = 1, current_state%local_grid%size(Z_INDEX)
          field_value%real_1d_array(k)=wwww_tot(k)
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
    else if (name .eq. "rain_mmr_total_local" ) then
       allocate(field_value%real_1d_array(current_state%local_grid%size(Z_INDEX)))
       do k = 1, current_state%local_grid%size(Z_INDEX)
          field_value%real_1d_array(k)=qr_tot(k)
       enddo
    else if (name .eq. "ice_mmr_total_local" ) then
       allocate(field_value%real_1d_array(current_state%local_grid%size(Z_INDEX)))
       do k = 1, current_state%local_grid%size(Z_INDEX)
          field_value%real_1d_array(k)=qi_tot(k)
       enddo
    else if (name .eq. "snow_mmr_total_local" ) then
       allocate(field_value%real_1d_array(current_state%local_grid%size(Z_INDEX)))
       do k = 1, current_state%local_grid%size(Z_INDEX)
          field_value%real_1d_array(k)=qs_tot(k)
       enddo
    else if (name .eq. "graupel_mmr_total_local" ) then
       allocate(field_value%real_1d_array(current_state%local_grid%size(Z_INDEX)))
       do k = 1, current_state%local_grid%size(Z_INDEX)
          field_value%real_1d_array(k)=qg_tot(k)
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
    else if (name .eq. "wtheta_cn_total_local") then
       allocate(field_value%real_1d_array(current_state%local_grid%size(Z_INDEX)))
       do k = 1, current_state%local_grid%size(Z_INDEX)
          field_value%real_1d_array(k)=wtheta_cn_tot(k)
       enddo
    else if (name .eq. "wtheta_ad_total_local") then
       allocate(field_value%real_1d_array(current_state%local_grid%size(Z_INDEX)))
       do k = 1, current_state%local_grid%size(Z_INDEX)
          field_value%real_1d_array(k)=wtheta_ad_tot(k)
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
    else if (name .eq. "uv_total_local") then
       allocate(field_value%real_1d_array(current_state%local_grid%size(Z_INDEX)))
       do k = 1, current_state%local_grid%size(Z_INDEX)
          field_value%real_1d_array(k)=uv_tot(k)
       enddo
    else if (name .eq. "wke_total_local") then
       allocate(field_value%real_1d_array(current_state%local_grid%size(Z_INDEX)))
       do k = 1, current_state%local_grid%size(Z_INDEX)
          field_value%real_1d_array(k)=wke_tot(k)
       enddo   
    else if (name .eq. "th2_total_local") then
       allocate(field_value%real_1d_array(current_state%local_grid%size(Z_INDEX)))
       do k = 1, current_state%local_grid%size(Z_INDEX)
          field_value%real_1d_array(k)=th2_tot(k)
       enddo
    else if (name .eq. "wqv_cn_total_local") then
       allocate(field_value%real_1d_array(current_state%local_grid%size(Z_INDEX)))
       do k = 1, current_state%local_grid%size(Z_INDEX)
          field_value%real_1d_array(k)=wqv_cn_tot(k)
       enddo
    else if (name .eq. "wqv_ad_total_local") then
       allocate(field_value%real_1d_array(current_state%local_grid%size(Z_INDEX)))
       do k = 1, current_state%local_grid%size(Z_INDEX)
          field_value%real_1d_array(k)=wqv_ad_tot(k)
       enddo 
    else if (name .eq. "wql_cn_total_local") then
       allocate(field_value%real_1d_array(current_state%local_grid%size(Z_INDEX)))
       do k = 1, current_state%local_grid%size(Z_INDEX)
          field_value%real_1d_array(k)=wql_cn_tot(k)
       enddo
    else if (name .eq. "wql_ad_total_local") then
       allocate(field_value%real_1d_array(current_state%local_grid%size(Z_INDEX)))
       do k = 1, current_state%local_grid%size(Z_INDEX)
          field_value%real_1d_array(k)=wql_ad_tot(k)
       enddo 
    else if (name .eq. "wqr_cn_total_local") then
       allocate(field_value%real_1d_array(current_state%local_grid%size(Z_INDEX)))
       do k = 1, current_state%local_grid%size(Z_INDEX)
          field_value%real_1d_array(k)=wqr_cn_tot(k)
       enddo
    else if (name .eq. "wqr_ad_total_local") then
       allocate(field_value%real_1d_array(current_state%local_grid%size(Z_INDEX)))
       do k = 1, current_state%local_grid%size(Z_INDEX)
          field_value%real_1d_array(k)=wqr_ad_tot(k)
       enddo
    else if (name .eq. "wqi_cn_total_local") then
       allocate(field_value%real_1d_array(current_state%local_grid%size(Z_INDEX)))
       do k = 1, current_state%local_grid%size(Z_INDEX)
          field_value%real_1d_array(k)=wqi_cn_tot(k)
       enddo
    else if (name .eq. "wqi_ad_total_local") then
       allocate(field_value%real_1d_array(current_state%local_grid%size(Z_INDEX)))
       do k = 1, current_state%local_grid%size(Z_INDEX)
          field_value%real_1d_array(k)=wqi_ad_tot(k)
       enddo   
    else if (name .eq. "wqs_cn_total_local") then
       allocate(field_value%real_1d_array(current_state%local_grid%size(Z_INDEX)))
       do k = 1, current_state%local_grid%size(Z_INDEX)
          field_value%real_1d_array(k)=wqs_cn_tot(k)
       enddo
    else if (name .eq. "wqs_ad_total_local") then
       allocate(field_value%real_1d_array(current_state%local_grid%size(Z_INDEX)))
       do k = 1, current_state%local_grid%size(Z_INDEX)
          field_value%real_1d_array(k)=wqs_ad_tot(k)
       enddo 
    else if (name .eq. "wqg_cn_total_local") then
       allocate(field_value%real_1d_array(current_state%local_grid%size(Z_INDEX)))
       do k = 1, current_state%local_grid%size(Z_INDEX)
          field_value%real_1d_array(k)=wqg_cn_tot(k)
       enddo
    else if (name .eq. "wqg_ad_total_local") then
       allocate(field_value%real_1d_array(current_state%local_grid%size(Z_INDEX)))
       do k = 1, current_state%local_grid%size(Z_INDEX)
          field_value%real_1d_array(k)=wqg_ad_tot(k)
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
    else if (name .eq. "rain_mmr_total_local" ) then
       allocate(field_value%real_1d_array(current_state%local_grid%size(Z_INDEX)))
       do k = 1, current_state%local_grid%size(Z_INDEX)
          field_value%real_1d_array(k)=qr_tot(k)
       enddo
    else if (name .eq. "ice_mmr_total_local" ) then
       allocate(field_value%real_1d_array(current_state%local_grid%size(Z_INDEX)))
       do k = 1, current_state%local_grid%size(Z_INDEX)
          field_value%real_1d_array(k)=qi_tot(k)
       enddo
    else if (name .eq. "snow_mmr_total_local" ) then
       allocate(field_value%real_1d_array(current_state%local_grid%size(Z_INDEX)))
       do k = 1, current_state%local_grid%size(Z_INDEX)
          field_value%real_1d_array(k)=qs_tot(k)
       enddo
    else if (name .eq. "graupel_mmr_total_local" ) then
       allocate(field_value%real_1d_array(current_state%local_grid%size(Z_INDEX)))
       do k = 1, current_state%local_grid%size(Z_INDEX)
          field_value%real_1d_array(k)=qg_tot(k)
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
    else if (name .eq. "i_wtheta_ad_total_local") then
       allocate(field_value%real_1d_array(current_state%local_grid%size(Z_INDEX)))
       do k = 1, current_state%local_grid%size(Z_INDEX)
          field_value%real_1d_array(k)=wtheta_ad_tot(k)
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
! =====================================================
    else if (name .eq. "cloud_mask") then
       allocate(field_value%real_3d_array(                                     &
                size(cloud_mask, 1), size(cloud_mask, 2), size(cloud_mask, 3)),&
                source=cloud_mask)
    else if (name .eq. "cloud_mask_total_local") then
       allocate(field_value%real_1d_array(current_state%local_grid%size(Z_INDEX)))
       do k = 1, current_state%local_grid%size(Z_INDEX)
          field_value%real_1d_array(k)=cloud_mask_tot(k)
       enddo
    else if (name .eq. "cloud_liq_mask_total_local") then
       allocate(field_value%real_1d_array(current_state%local_grid%size(Z_INDEX)))
       do k = 1, current_state%local_grid%size(Z_INDEX)
          field_value%real_1d_array(k)=cloud_liq_mask_tot(k)
       enddo
    else if (name .eq. "cloud_ice_mask_total_local") then
       allocate(field_value%real_1d_array(current_state%local_grid%size(Z_INDEX)))
       do k = 1, current_state%local_grid%size(Z_INDEX)
          field_value%real_1d_array(k)=cloud_ice_mask_tot(k)
       enddo
    end if
  end subroutine field_value_retrieval_callback

  subroutine calculate_wq(current_state, jcol, icol, iq, wq_cn, wq_ad, advect_q)
    
    type(model_state_type), target, intent(inout) :: current_state
    
    real(kind=DEFAULT_PRECISION), intent(inout) :: wq_cn(:), wq_ad(:)
    character(len=*), intent(in) :: advect_q
    integer, intent(in) :: jcol, icol, iq
    
    integer :: k
  
    if (trim(advect_q) .eq. "pw") then 
       do k=1, current_state%local_grid%size(Z_INDEX)-1
          wq_ad(k) = wq_ad(k) + ( 0.5 * ( & 
               current_state%q(iq)%data(k,jcol,icol) + &
               current_state%q(iq)%data(k+1,jcol,icol)) * &
               current_state%w%data(k,jcol,icol))
          wq_cn(k) = wq_ad(k)
       enddo
    else if (trim(advect_q) .eq. "tvd") then
       do k=2, current_state%local_grid%size(Z_INDEX)-1
          wq_ad(k) = wq_ad(k) + (current_state%w%data(k,jcol,icol) & 
               * tvd_dgs_terms%adv_q_dgs(k+1, jcol, icol, iq) )
          wq_cn(k) = wq_cn(k) +  ( 0.5 * (  & 
               current_state%q(iq)%data(k,jcol,icol) + &
               current_state%q(iq)%data(k+1,jcol,icol)) * &
               current_state%w%data(k,jcol,icol) )
       enddo
    endif

  end subroutine calculate_wq
    
  !---------------------------------------------------------------------
  !> This routine calculates:
  !    cloud_mask              a binary 3D total cloud mask
  !                            optional: l_cloud_mask
  !    cloud_mask_tot          the total cloud mask profile local sum
  !    cloud_liq_mask_tot      the liquid cloud mask profile local sum
  !    cloud_ice_mask_tot      the ice cloud mask profile local sum
  !  The latter 3 are intended to be transformed into
  !  cloud fraction profile diagnostics via xml processing.
  !
  !
  !> SOCRATES method
  !  The definition used here is consistent with that within 
  !  SOCRATES, where i_cloud_representation = 2, that is, liquid and
  !  ice cloud are treated separately, with their individual fractions
  !  within a cell summing to 1 based their mixing ratios relative to
  !  the cell total.  
  !  Each element of the SOCRATES calculation is reproduced
  !  here because the current formulation within SOCRATES will
  !  not be available here when that component is not enabled.
  !    See: def_merge_atm.F90
  !         merge_atm_data.F90
  !  A better solution might be to have both work from 
  !  the same source for easy consistency in the case where
  !  definitions of cloudy cells were to change.
  !  We do not apply special consideration to the cases where
  !  MONC is run with SOCRATES and i_cloud_representation != 2.
  !  That is, values are calculated here even when cloud is 
  !  off in SOCRATES, and no stratiform/convective distinction
  !  is made. 
  !
  !
  !> DEFAULT method
  !  Cloud fraction is based on exceeding qlcrit and qicrit.
  !  This definition is consistent with that used in the conditional
  !  diagnostics routine (condition "AC").  It does not include
  !  consideration of rain, snow, and graupel fields.  Liquid and
  !  ice cloud are treated separately, with their individual fractions
  !  within a cell summing to 1 based their mixing ratios relative to
  !  the cell total.
  subroutine calculate_cloud_mask(current_state, jcol, icol)
    type(model_state_type), target, intent(inout) :: current_state
    integer, intent(in) :: jcol, icol
    integer :: k, target_y_index, target_x_index
    logical :: l_prepare_3d_mask, cloud_present

    ! The factors below were derived as part of J. Petch PhD
    ! These are used in the SOCRATES method.
    real(kind=DEFAULT_PRECISION) ::                  &
         rainfac  = 0.02,                            &
         snowfac  = 0.40,                            &
         graupfac = 0.05

    ! Local temporary terms
    real(kind=DEFAULT_PRECISION) :: templ, tempi, tempt

    target_y_index = jcol - current_state%local_grid%halo_size(Y_INDEX)
    target_x_index = icol - current_state%local_grid%halo_size(X_INDEX)

    templ = 0.0_DEFAULT_PRECISION
    tempi = 0.0_DEFAULT_PRECISION    
    tempt = 0.0_DEFAULT_PRECISION

    l_prepare_3d_mask = allocated(cloud_mask)

    do k=1, current_state%local_grid%size(Z_INDEX)

      !> Collect available condensate amounts
      if (iql > 0) &
        templ = current_state%q(iql)%data(k, jcol, icol)
      if (iqi > 0) &
        tempi = current_state%q(iqi)%data(k, jcol, icol)

      !> Check cloud_mask_method and modify as needed
      !> The SOCRATES method considers rain, snow, and graupel.
      if (cloud_mask_method == "SOCRATES") then
        if (iqr > 0) &
          templ = templ + rainfac  * current_state%q(iqr)%data(k, jcol, icol)
        if (iqs > 0) &
          tempi = tempi + snowfac  * current_state%q(iqs)%data(k, jcol, icol)
        if (iqg > 0) &
          tempi = tempi + graupfac * current_state%q(iqg)%data(k, jcol, icol)      
      endif ! check cloud_mask_method

      !> Work out cloud fractions
      tempt = templ + tempi
      if (cloud_mask_method == "SOCRATES") then
        cloud_present = (tempt > EPSILON(tempt))
      else ! DEFAULT
        cloud_present = (templ > qlcrit .or. tempi > qicrit .or. (templ+tempi) > qlcrit)
      end if 

      if (cloud_present) then
        if (l_prepare_3d_mask)    &
           cloud_mask(k, target_y_index, target_x_index) = 1.0_DEFAULT_PRECISION
        cloud_mask_tot(k) = cloud_mask_tot(k) + 1.0_DEFAULT_PRECISION

        if (l_partial_liq_ice) then ! separated cloud, partial coverage
          cloud_liq_mask_tot(k) = cloud_liq_mask_tot(k) + (templ / tempt)
          cloud_ice_mask_tot(k) = cloud_ice_mask_tot(k) + (tempi / tempt)
        else ! homogeneous cloud, full coverage for both types of condensate
          if (templ > qlcrit)         &
             cloud_liq_mask_tot(k) = cloud_liq_mask_tot(k) + 1.0_DEFAULT_PRECISION
          if (tempi > qicrit)         &
             cloud_ice_mask_tot(k) = cloud_ice_mask_tot(k) + 1.0_DEFAULT_PRECISION
        end if ! check liquid/ice partition method
      end if ! check cloud_present

    end do ! k loop over vertical model levels

  end subroutine calculate_cloud_mask


end module profile_diagnostics_mod
