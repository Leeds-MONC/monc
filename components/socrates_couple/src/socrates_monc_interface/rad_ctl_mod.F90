module rad_ctl_mod

implicit none

contains

subroutine rad_ctl(current_state, sw_spectrum, lw_spectrum, &
     mcc, socrates_opt, merge_fields, socrates_derived_fields)

  USE def_spectrum, ONLY: StrSpecData
  USE sw_control_mod, ONLY: sw_control
  USE lw_control_mod, ONLY: lw_control
  USE lw_rad_input_mod, ONLY: lw_input
  USE def_dimen,    ONLY: StrDim
  USE def_atm,     ONLY: StrAtm, deallocate_atm
  USE def_cld,      ONLY: StrCld, deallocate_cld,          &
                                  deallocate_cld_prsc,     &
                                  deallocate_cld_mcica
  USE def_aer,      ONLY: StrAer,  deallocate_aer,          &
                                   deallocate_aer_prsc
  USE def_bound,    ONLY: StrBound,  deallocate_bound
  USE def_out,      ONLY: StrOut, deallocate_out
  use state_mod, only : model_state_type
  use grids_mod, only : Z_INDEX, Y_INDEX, X_INDEX
  use def_mcc_profiles, only: str_mcc_profiles
  use def_merge_atm, only: str_merge_atm
  use def_socrates_options, only: str_socrates_options
  use def_socrates_derived_fields, only: str_socrates_derived_fields
  use science_constants_mod, only: cp
  use rad_pcf, only: ip_solar, ip_infra_red

  type(model_state_type), target, intent(inout) :: current_state
  TYPE (StrSpecData) :: sw_spectrum
  TYPE (StrSpecData) :: lw_spectrum
  ! Dimensions:
  TYPE(StrDim) :: dimen
  ! Boundary conditions:
  TYPE(StrBound) :: bound
  ! atmosphere profiles
  TYPE(StrAtm) :: atm
  ! Cloud data:
  TYPE(StrCld) :: cld
  ! Aerosol data:
  TYPE(StrAer) :: aer
  ! radiation output
  TYPE(StrOut) :: radout
  ! McClatchey profile data for leves
  type (str_mcc_profiles), intent(in) :: mcc
  type (str_merge_atm), intent(inout) :: merge_fields
  ! MONC  options read from configuration
  type (str_socrates_options), intent(in) :: socrates_opt
  type (str_socrates_derived_fields), intent(inout) :: socrates_derived_fields

  integer :: icol, jcol ! Shorthand column indices from MONC (include halos)
  integer :: target_x_index, target_y_index ! Shorthand column indices for radiation 
                                            !(MONC index less halos)
  integer :: k          ! index for monc k loop to work out heating rates
  integer :: k_top      ! top of the monc domain

  ! work out column indexes
  icol=current_state%column_local_x
  jcol=current_state%column_local_y
  ! Use the radiation interface indexing, which does not include the halo 
  target_y_index=jcol-current_state%local_grid%halo_size(Y_INDEX)
  target_x_index=icol-current_state%local_grid%halo_size(X_INDEX)

  k_top = current_state%local_grid%size(Z_INDEX)

  ! ********************** Shortwave Calculation begin ***************************!
  if (socrates_derived_fields%fraction_lit .ge. 1.0) then 
     ! SOCRATES only call for shortwave when the fraction_lit = 1 
     ! This is set in the interface (socrate_couple.F90). It is important
     ! to note that this method and condition, does not permit a 
     ! partial lit grid
     ! DEPENDS ON: set_control
     call set_control(sw_control, sw_spectrum)
     ! last argument needs to be cloud levels
     ! DEPENDS ON: set_dimen
     call set_dimen(sw_control, dimen, sw_spectrum,                     &
          1, mcc%irad_levs, mcc%irad_levs)
     ! DEPENDS ON: set_atm
     ! n_profile in set_atm is the number of columns being sent to
     ! socrates. Here it is hard-coded to 1, since socrates is called
     ! per column in MONC
     call set_atm(sw_control, atm, dimen, sw_spectrum,                 &
          1 , mcc%irad_levs, socrates_opt, merge_fields)
     ! DEPENDS ON: set_bound
     call set_bound(sw_control, ip_solar, atm, dimen, sw_spectrum, bound, &
          socrates_derived_fields)
     ! DEPENDS ON: set_cld
     call set_cld(sw_control, atm, dimen, sw_spectrum, cld,              &
          1 , mcc%irad_levs, mcc%irad_levs, socrates_opt, merge_fields)
     ! DEPENDS ON: set_aer
     call set_aer(sw_control, atm, dimen, sw_spectrum, aer )

     ! DEPENDS ON: radiance_calc
     call radiance_calc(sw_control, dimen, sw_spectrum, atm, cld, aer, bound, radout)

     ! store the top of atmospher and surface values for diagnostics
     socrates_derived_fields%toa_up_shortwave(target_y_index,target_x_index) = &
          radout%flux_up(1,0,1)
     socrates_derived_fields%toa_down_shortwave(target_y_index,target_x_index) = &
          radout%flux_down(1,0,1)
     socrates_derived_fields%surface_up_shortwave(target_y_index,target_x_index) =  &
          radout%flux_up(1,mcc%irad_levs,1)
     socrates_derived_fields%surface_down_shortwave(target_y_index,target_x_index) = &
          radout%flux_down(1,mcc%irad_levs,1)

     do k = 1, mcc%irad_levs
        merge_fields%sw_heat_rate_radlevs(k) = &
             ! net flux at k-1 
             ((radout%flux_down(1,k-1,1)-radout%flux_up(1,k-1,1)) &
             ! net flux at k
             - (radout%flux_down(1,k,1)-radout%flux_up(1,k,1)))/ &
             (merge_fields%mass(k)*cp)
     enddo
     
     do k = 1, k_top
        socrates_derived_fields%flux_up_sw(k,target_y_index, target_x_index) = &
             radout%flux_up(1,mcc%irad_levs+1-k,1)
        
        socrates_derived_fields%flux_down_sw(k,target_y_index, target_x_index) = &
             radout%flux_down(1,mcc%irad_levs+1-k,1)
     enddo

     current_state%sth_sw%data(1,jcol, icol) = 0.0
     !     NB: the +2 is becasue n=1 is below the surface (comment from LEM)
     current_state%sth_sw%data(2:k_top,jcol, icol) = &
          merge_fields%sw_heat_rate_radlevs(mcc%irad_levs:mcc%irad_levs+2-k_top:-1)
     socrates_derived_fields%swrad_hr(:,target_y_index, target_x_index) = &
          current_state%sth_sw%data(:,jcol, icol)
     ! convert dT/dt to dTH/dt
     current_state%sth_sw%data(:, jcol, icol) = & 
         current_state%sth_sw%data(:, jcol, icol)* &
         current_state%global_grid%configuration%vertical%prefrcp(:)

  else
     do k = 1, k_top
        socrates_derived_fields%flux_up_sw(k,target_y_index, target_x_index) = &
             0.0
        
        socrates_derived_fields%flux_down_sw(k,target_y_index, target_x_index) = &
             0.0
     enddo
     socrates_derived_fields%toa_up_shortwave(target_y_index,target_x_index) = 0.0
     socrates_derived_fields%toa_down_shortwave(target_y_index,target_x_index) = 0.0
     socrates_derived_fields%surface_up_shortwave(target_y_index,target_x_index) = 0.0
     socrates_derived_fields%surface_down_shortwave(target_y_index,target_x_index) = 0.0
     current_state%sth_sw%data(:,jcol, icol) = 0.0
     socrates_derived_fields%swrad_hr(:,target_y_index, target_x_index) = 0.0
     
  endif

  CALL deallocate_out(radout)
  CALL deallocate_aer_prsc(aer)
  CALL deallocate_aer(aer)
  !!$  CALL deallocate_cld_mcica(cld)
  CALL deallocate_bound(bound)
  CALL deallocate_cld_prsc(cld)
  CALL deallocate_cld(cld)
  call deallocate_atm(atm)
  ! ********************** Shortwave Calculation end ******************************!
  
  ! ********************** Longwave Calculation begin *****************************!
  ! DEPENDS ON: set_control
  call set_control(lw_control, lw_spectrum)
  ! DEPENDS ON: set_dimen 
  call set_dimen(lw_control, dimen, lw_spectrum,                                 &
       1, mcc%irad_levs, mcc%irad_levs)
  ! DEPENDS ON: set_atm
  call set_atm(lw_control, atm, dimen, lw_spectrum,                 &
       1 , mcc%irad_levs, socrates_opt, merge_fields)
  ! DEPENDS ON: set_bound
  call set_bound(lw_control, ip_infra_red, atm, dimen, lw_spectrum, bound, &
       socrates_derived_fields)
  ! DEPENDS ON: set_cld
  call set_cld(lw_control, atm, dimen, lw_spectrum, cld ,              &
       1 , mcc%irad_levs, mcc%irad_levs, socrates_opt, merge_fields)
  ! DEPENDS ON: set_aer
  call set_aer(lw_control, atm, dimen, lw_spectrum, aer )
  ! DEPENDS ON: radiance_calc
  call radiance_calc(lw_control, dimen, lw_spectrum, atm, cld, aer, bound, radout)

  socrates_derived_fields%toa_up_longwave(target_y_index,target_x_index) = & 
       radout%flux_up(1,0,1)
  socrates_derived_fields%surface_up_longwave(target_y_index,target_x_index) =  &
       radout%flux_up(1,mcc%irad_levs,1)
  socrates_derived_fields%surface_down_longwave(target_y_index,target_x_index) = &
       radout%flux_down(1,mcc%irad_levs,1)

  do k = 1, mcc%irad_levs
     merge_fields%lw_heat_rate_radlevs(k) = &
          ! net flux at k-1 
          ((radout%flux_down(1,k-1,1)-radout%flux_up(1,k-1,1)) &
          ! net flux at k
          - (radout%flux_down(1,k,1)-radout%flux_up(1,k,1)))/ &
          (merge_fields%mass(k)*cp)
  enddo
     
  current_state%sth_lw%data(1,jcol, icol) = 0.0
     !     NB: the +2 is becasue n=1 is below the surface (comment from LEM)
  current_state%sth_lw%data(2:k_top,jcol, icol) = &
          merge_fields%lw_heat_rate_radlevs(mcc%irad_levs:mcc%irad_levs+2-k_top:-1)
  ! convert dT/dt to dTH/dt
  current_state%sth_lw%data(:, jcol, icol) = & 
         current_state%sth_lw%data(:, jcol, icol)* &
         current_state%global_grid%configuration%vertical%prefrcp(:)

  socrates_derived_fields%lwrad_hr(:,target_y_index, target_x_index) = &
       current_state%sth_lw%data(:,jcol, icol)
  
  do k = 1, k_top
     socrates_derived_fields%flux_up_lw(k,target_y_index, target_x_index) = &
          radout%flux_up(1,mcc%irad_levs+1-k,1)
  
     socrates_derived_fields%flux_down_lw(k,target_y_index, target_x_index) = &
          radout%flux_down(1,mcc%irad_levs+1-k,1)
  enddo

  socrates_derived_fields%totrad_hr(:,target_y_index, target_x_index) = &
       socrates_derived_fields%lwrad_hr(:,target_y_index, target_x_index) + &
       socrates_derived_fields%swrad_hr(:,target_y_index, target_x_index)
  
  CALL deallocate_out(radout)
  CALL deallocate_aer_prsc(aer)
  CALL deallocate_aer(aer)
  CALL deallocate_bound(bound)
  CALL deallocate_cld_prsc(cld)
  CALL deallocate_cld(cld)
  CALL deallocate_atm(atm)
  ! ********************** Longwave Calculation end ********************************!

end subroutine rad_ctl

end module rad_ctl_mod
