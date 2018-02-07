module merge_atm_data

  use datadefn_mod, only : DEFAULT_PRECISION
  use grids_mod, only : Z_INDEX, vertical_grid_configuration_type
  use state_mod, only : model_state_type
  use def_socrates_options, only: str_socrates_options
  use def_merge_atm, only: str_merge_atm
  use def_mcc_profiles, only: str_mcc_profiles
  use def_socrates_derived_fields, only: str_socrates_derived_fields

  implicit none

contains
    
  subroutine merge_data(current_state, socrates_opt, socrates_derived_fields, merge_fields, mcc)
    !
    !------------------------------------------------------------------------
    ! Routine to add ozone, qv and temperature information to the radiation
    ! profiles (including MONC levels). This is called in the timestep
    ! callback and uses irad_levs from calculate_radiation_levels.
    ! Values of t and q above the MONC levels and ozone from all levels are 
    ! taken from a McClattchey profiles, which are defined in the
    ! configuration file.
    !------------------------------------------------------------------------
    !
    type(model_state_type), target, intent(inout) :: current_state
    type (str_socrates_options), intent(in) :: socrates_opt
    type (str_socrates_derived_fields), intent(in) :: socrates_derived_fields
    type (str_mcc_profiles), intent(in)  :: mcc
    type (str_merge_atm), intent(inout) :: merge_fields

    type(vertical_grid_configuration_type) :: vertical_grid

    real(kind=DEFAULT_PRECISION), parameter :: gravity = 9.80665
    ! weights for working out the temperature on grid boundaries
    ! (weights are based on height)
    real(kind=DEFAULT_PRECISION) :: weight_upper, weight_lower

    real(kind=DEFAULT_PRECISION) :: cloud_total
    
    integer :: k_top, icol, jcol ! Shorthand column indices

    integer :: k ! loop counter for z
    integer :: j ! loop counter for mcc_levs

    k_top = current_state%local_grid%size(Z_INDEX)
    icol=current_state%column_local_x
    jcol=current_state%column_local_y
    
    vertical_grid=current_state%global_grid%configuration%vertical
    
    ! set the column vapour and temperature for clear sky
    ! and cloudy calculation by combining the monc profiles and the McClatchey profiles
    ! Remember the McClatchey profiles go from top of atmosphere as does the
    ! radiation calc. 
    ! Also, it is important to note that I have assumed the McClatchey values
    ! exist on the z levels (where w is) not zn.

    merge_fields%pres_level(1:mcc%cut) = mcc%p_level(1:mcc%cut)
    merge_fields%pres_n(1:mcc%cut) = mcc%p_n(1:mcc%cut)
    merge_fields%o3_n(1:mcc%cut) = mcc%o3_n(1:mcc%cut)
    merge_fields%t_level(1:mcc%cut) = mcc%t_level(1:mcc%cut)
    merge_fields%t_n(1:mcc%cut) = mcc%t_n(1:mcc%cut)
    merge_fields%qv_n(1:mcc%cut) = mcc%q_n(1:mcc%cut)
    merge_fields%ql_n(:) = 0.0
    merge_fields%qi_n(:) = 0.0
    merge_fields%total_cloud_fraction(:) = 0.0
    merge_fields%liquid_cloud_fraction(:) = 0.0
    merge_fields%ice_cloud_fraction(:) = 0.0
   
    cloud_total = 0.0
!
! derive absolute temperature
    do k=1, k_top
       merge_fields%t_n_loc(k) = (current_state%th%data(k, jcol, icol)                       &
            + current_state%global_grid%configuration%vertical%thref(k))        &
            * current_state%global_grid%configuration%vertical%rprefrcp(k)
    enddo
! merge the fields centre grid fields that come from MONC first     
    do k = 2, k_top
       ! centre pressure
       merge_fields%pres_n(k+mcc%cut)= &
            current_state%global_grid%configuration%vertical%prefn(k_top+2-k)      
       ! Absolute temperature at centre of the grid
       merge_fields%t_n(k+mcc%cut) = merge_fields%t_n_loc(k_top+2-k)
       ! vapour mixing ratio
       merge_fields%qv_n(k+mcc%cut) = current_state%q(socrates_opt%iqv)%data(k_top+2-k, jcol, icol)
       if ( socrates_opt%cloud_representation == 2) then
          if (socrates_opt%mphys_nq_l > 0) then
             merge_fields%ql_n(k+mcc%cut) =  &
                  current_state%q(socrates_opt%iql)%data(k_top+2-k, jcol, icol)
             if (.not. socrates_opt%l_fix_re) then
                if (socrates_opt%mphys_nd_l > 0) then
                   ! Cloud is defined as double moment so set the cloud number for the
                   ! calculation of effective radius
                   if (socrates_opt%inl > 0) then
                      merge_fields%cloudnumber_n(k+mcc%cut) =  &
                           current_state%q(socrates_opt%inl)%data(k_top+2-k, jcol, icol)
                   else
                      merge_fields%cloudnumber_n(k+mcc%cut) =  &
                           socrates_opt%fixed_cloud_number*1.e6 ! convert to number per m3
                   endif
                endif
             endif      
          endif
          ! NOTE: if rain is on then the mass is added to ql after
          !       multiplication with rainfac. If available, there should
          !       be a seperate optical property for rain.
          if (socrates_opt%mphys_nq_r > 0) then
             merge_fields%ql_n(k+mcc%cut) = merge_fields%ql_n(k+mcc%cut) + &
                  merge_fields%rainfac*        &
                  (current_state%q(socrates_opt%iqr)%data(k_top+2-k, jcol, icol))
          endif
          if (socrates_opt%mphys_nq_i > 0) then
             merge_fields%qi_n(k+mcc%cut) = &
                  current_state%q(socrates_opt%iqi)%data(k_top+2-k, jcol, icol)
          endif
          ! NOTE: if snow and graupel are on then the mass is added to qi after
          !       multiplication with snowfac and graupfac. If available,
          !       there should be a seperate optical property for snow and
          !       graupel.
          if (socrates_opt%mphys_nq_s > 0) then
             merge_fields%qi_n(k+mcc%cut) = merge_fields%qi_n(k+mcc%cut) +         &
                  merge_fields%snowfac*              &
                  (current_state%q(socrates_opt%iqs)%data(k_top+2-k, jcol, icol))
          endif
          if (socrates_opt%mphys_nq_g > 0) then
             merge_fields%qi_n(k+mcc%cut) = merge_fields%qi_n(k+mcc%cut) +         &
                  merge_fields%graupfac*                &
                  (current_state%q(socrates_opt%iqg)%data(k_top+2-k, jcol, icol))
          endif
          ! work out cloud fractions
          cloud_total = merge_fields%ql_n(k+mcc%cut) + merge_fields%qi_n(k+mcc%cut)
          if (cloud_total > EPSILON(cloud_total) ) then
             merge_fields%total_cloud_fraction(k+mcc%cut) = 1.0_DEFAULT_PRECISION
             merge_fields%liquid_cloud_fraction(k+mcc%cut) = &
                  merge_fields%ql_n(k+mcc%cut)/(cloud_total)
             merge_fields%ice_cloud_fraction(k+mcc%cut) = &
                  (1.0_DEFAULT_PRECISION - merge_fields%liquid_cloud_fraction(k+mcc%cut))
          else
             merge_fields%total_cloud_fraction(k+mcc%cut) = 0.0_DEFAULT_PRECISION
             merge_fields%liquid_cloud_fraction(k+mcc%cut) = 0.0_DEFAULT_PRECISION
             merge_fields%ice_cloud_fraction(k+mcc%cut) = 0.0_DEFAULT_PRECISION
          endif
       endif
    enddo

    ! Work out k_top value for the centre fields. This is the mean of the monc centre value and 
    ! McClatchey centre fields
    ! Pressure at the merge cut-off
    merge_fields%pres_n(mcc%cut+1)= &
         (current_state%global_grid%configuration%vertical%prefn(k_top) + &
         mcc%p_n(mcc%cut))/2.0
    ! Temperature at the merge cut-off
    merge_fields%t_n(mcc%cut+1)= &
         (merge_fields%t_n_loc(k_top) + mcc%t_n(mcc%cut))/2.0
    ! Temperature at the merge cut-off
    merge_fields%qv_n(mcc%cut+1)= &
         (current_state%q(socrates_opt%iqv)%data(k_top, jcol, icol) + mcc%q_n(mcc%cut))/2.0

    ! Now sort the Ozone profile which only exists on
    ! McClatchey levels, so needs to be merged on to MONC
    ! levels. Check this code!!

    do k=k_top+mcc%cut,mcc%cut,-1  
       if (merge_fields%pres_level(k).gt.mcc%p_level(mcc%levs)) then
          merge_fields%o3_n(k) = mcc%o3_n(mcc%levs)
       else
          do j=mcc%levs,1,-1  
             if (merge_fields%pres_level(k).gt.mcc%p_level(j)) then
                merge_fields%o3_n(k) = mcc%o3_n(j)
                exit
             endif
          enddo
       endif
    enddo
    
    ! Now work out the temperature and pressure on levels
    ! Set the surface 
    merge_fields%pres_level(mcc%irad_levs) = current_state%surface_pressure
    merge_fields%t_level(mcc%irad_levs) = socrates_derived_fields%srf_temperature
    ! set the zero level (which is TOA) of t_level and pres_level (based on assumptions from LEM)
    merge_fields%pres_level(0) = 0.0
    merge_fields%t_level(0) = merge_fields%t_level(1)

    ! Assume the level is mid-way between n levels. While this is the same as 
    ! the LEM assumption, this is not correct. There should be a weighting factor
    ! that accounts for a stretched MONC grid
    do k = 1, mcc%irad_levs-1
       merge_fields%pres_level(k) = 0.5_DEFAULT_PRECISION * &
            (merge_fields%pres_n(k) + merge_fields%pres_n(k+1))
       
       merge_fields%t_level(k) = 0.5_DEFAULT_PRECISION * &
            (merge_fields%t_n(k) + merge_fields%t_n(k+1))
    enddo

    ! mass of the atmosphere on irad_levs for set_atm
    do k=1, mcc%irad_levs
       merge_fields%mass(k) = &
            (merge_fields%pres_level(k)-merge_fields%pres_level(k-1))/gravity
    enddo

    end subroutine merge_data

end module merge_atm_data
