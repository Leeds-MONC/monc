!> This module sets up the logicals and parameters for the edward-slingo code
!> from the UM. It also calls the shortwave and longwave ES code and
!> outputs the heating rates and fluxes.
module socrates_couple_mod
  use datadefn_mod, only : DEFAULT_PRECISION
  use monc_component_mod, only : COMPONENT_SCALAR_FIELD_TYPE, COMPONENT_DOUBLE_DATA_TYPE, &
       COMPONENT_ARRAY_FIELD_TYPE, COMPONENT_INTEGER_DATA_TYPE, component_descriptor_type, &
        component_field_value_type, component_field_information_type
  use optionsdatabase_mod, only : options_get_string, options_get_integer, &
       options_get_real, options_get_logical
  use state_mod, only : FORWARD_STEPPING, model_state_type, PRESCRIBED_SURFACE_FLUXES, &
       PRESCRIBED_SURFACE_VALUES
  use grids_mod, only : Z_INDEX, Y_INDEX, X_INDEX
  use science_constants_mod, only : cp
  use conversions_mod, only : conv_to_string
  use logging_mod, only : LOG_INFO, LOG_WARN, LOG_ERROR, &
       LOG_DEBUG, log_master_log, log_log, log_get_logging_level, &
       log_master_log
  use q_indices_mod, only: get_q_index, standard_q_names
  use registry_mod, only : is_component_enabled
  use sw_rad_input_mod, only: sw_input
  use lw_rad_input_mod, only: lw_input
  ! this module exists in the radiance_core dir
  USE def_spectrum, ONLY: StrSpecData
  ! monc socrates couple modules, based on UM
  USE sw_control_mod, ONLY: sw_control
  USE lw_control_mod, ONLY: lw_control
  USE rad_ctl_mod, ONLY: rad_ctl
  ! monc socrates couple modules, only used in MONC
  use science_constants_mod, ONLY: r_over_cp
  use def_merge_atm, only: str_merge_atm,  allocate_merge_data_fields
  use def_mcc_profiles, only: str_mcc_profiles
  use def_socrates_options, only: str_socrates_options
  use def_socrates_derived_fields, only: str_socrates_derived_fields
  use get_and_test_socrates_options_mod, only: set_and_test_socrates_monc_options
  use mcclatchey_profiles, only: read_mcclatchey_profiles, &
       calculate_radiation_levels
  use merge_atm_data, only: merge_data
  use solar_position_angle_mod, only: solar_pos_calculation, solar_angle_calculation

  implicit none

  ! Index for the top of the domain
  integer :: k_top, x_local, y_local, x_nohalos, y_nohalos
  !
  ! local density factor and raidiation factor used in conversions
  real(kind=DEFAULT_PRECISION), allocatable ::  density_factor(:), radiation_factor(:)
  !
  TYPE (StrSpecData), SAVE :: sw_spectrum
  TYPE (StrSpecData), SAVE :: lw_spectrum

  type (str_mcc_profiles) :: mcc
  type (str_merge_atm) :: merge_fields
  type (str_socrates_options) :: socrates_opt
  type (str_socrates_derived_fields) :: socrates_derived_fields

  public socrates_couple_get_descriptor
contains

  !> Provides the descriptor back to the caller and is used in component registration
  !! @returns The termination check component descriptor
  type(component_descriptor_type) function socrates_couple_get_descriptor()
    socrates_couple_get_descriptor%name="socrates_couple"
    socrates_couple_get_descriptor%version=0.1
    socrates_couple_get_descriptor%initialisation=>initialisation_callback
    socrates_couple_get_descriptor%timestep=>timestep_callback
    socrates_couple_get_descriptor%finalisation=>finalisation_callback

    socrates_couple_get_descriptor%field_value_retrieval=>field_value_retrieval_callback
    socrates_couple_get_descriptor%field_information_retrieval=>field_information_retrieval_callback
    
    allocate(socrates_couple_get_descriptor%published_fields(14))
    
    socrates_couple_get_descriptor%published_fields(1)="flux_up_shortwave"
    socrates_couple_get_descriptor%published_fields(2)="flux_down_shortwave"
    socrates_couple_get_descriptor%published_fields(3)="flux_up_longwave"
    socrates_couple_get_descriptor%published_fields(4)="flux_down_longwave"
    socrates_couple_get_descriptor%published_fields(5)="toa_up_longwave"
    socrates_couple_get_descriptor%published_fields(6)="surface_down_longwave"
    socrates_couple_get_descriptor%published_fields(7)="surface_up_longwave"
    socrates_couple_get_descriptor%published_fields(8)="toa_down_shortwave"
    socrates_couple_get_descriptor%published_fields(9)="toa_up_shortwave"
    socrates_couple_get_descriptor%published_fields(10)="surface_down_shortwave"
    socrates_couple_get_descriptor%published_fields(11)="surface_up_shortwave"
    socrates_couple_get_descriptor%published_fields(12)="shortwave_heating_rate"
    socrates_couple_get_descriptor%published_fields(13)="longwave_heating_rate"
    socrates_couple_get_descriptor%published_fields(14)="total_radiative_heating_rate"    
       
  end function socrates_couple_get_descriptor

  !> The initialisation callback sets up the prescribed longwave fluxes and the
  !> exponential decay factor
  !! @param current_state The current model state
  subroutine initialisation_callback(current_state)
    type(model_state_type), target, intent(inout) :: current_state

    integer :: k ! look counter

    k_top=current_state%local_grid%size(Z_INDEX) + current_state%local_grid%halo_size(Z_INDEX) * 2
    y_local=current_state%local_grid%size(Y_INDEX) + current_state%local_grid%halo_size(Y_INDEX) * 2
    x_local=current_state%local_grid%size(X_INDEX) + current_state%local_grid%halo_size(X_INDEX) * 2
    
    y_nohalos=current_state%local_grid%size(Y_INDEX)
    x_nohalos=current_state%local_grid%size(X_INDEX) 
    
    !if (.not. current_state%continuation_run) then
       ! Allocated the longwave and shortwave heating rates.They need to be added to
       ! the dump
    !   allocate(current_state%sth_lw%data(k_top, y_local, x_local))
    !   allocate(current_state%sth_sw%data(k_top, y_local, x_local))
    !   current_state%sth_lw%data(:,:,:) = 0.0
    !   current_state%sth_sw%data(:,:,:) = 0.0
       ! Allocate downward surface fluxes
    !   allocate(current_state%sw_down_surf(y_local, x_local))
    !   allocate(current_state%lw_down_surf(y_local, x_local))
    !   current_state%sw_down_surf(:, :) = 0.0
    !   current_state%lw_down_surf(:, :) = 0.0
    !endif

    ! allocate the density and radiation factor needed for heating rates
    allocate(socrates_derived_fields%density_factor(k_top))
    allocate(socrates_derived_fields%radiation_factor(k_top))
    !
    ! allocate the fluxes, which do not need to be dumped
    allocate(socrates_derived_fields%flux_up_sw(k_top,y_nohalos,x_nohalos))
    allocate(socrates_derived_fields%flux_down_sw(k_top,y_nohalos,x_nohalos))
    allocate(socrates_derived_fields%flux_up_lw(k_top,y_nohalos,x_nohalos))
    allocate(socrates_derived_fields%flux_down_lw(k_top,y_nohalos,x_nohalos))
    allocate(socrates_derived_fields%flux_net_sw(k_top,y_nohalos,x_nohalos))
    allocate(socrates_derived_fields%flux_net_lw(k_top,y_nohalos,x_nohalos))
    !
    ! allocate the surface and TOA fields which are used for diagnostics
    allocate(socrates_derived_fields%toa_up_longwave(y_nohalos,x_nohalos))
    allocate(socrates_derived_fields%toa_down_shortwave(y_nohalos,x_nohalos))
    allocate(socrates_derived_fields%toa_up_shortwave(y_nohalos,x_nohalos))
    allocate(socrates_derived_fields%surface_up_longwave(y_nohalos,x_nohalos))
    allocate(socrates_derived_fields%surface_down_longwave(y_nohalos,x_nohalos))
    allocate(socrates_derived_fields%surface_down_shortwave(y_nohalos,x_nohalos))
    allocate(socrates_derived_fields%surface_up_shortwave(y_nohalos,x_nohalos))
    !
    ! allocate radiative heating rates for diagnostics
    allocate(socrates_derived_fields%swrad_hr(k_top, y_nohalos,x_nohalos))
    allocate(socrates_derived_fields%lwrad_hr(k_top, y_nohalos,x_nohalos))
    allocate(socrates_derived_fields%totrad_hr(k_top, y_nohalos,x_nohalos))

    ! initialise allocated variables to 0 for safety
    socrates_derived_fields%flux_up_sw(:,:,:) = 0.0
    socrates_derived_fields%flux_down_sw(:,:,:) = 0.0
    socrates_derived_fields%flux_up_lw(:,:,:) = 0.0
    socrates_derived_fields%flux_down_lw(:,:,:) = 0.0
    socrates_derived_fields%flux_net_sw(:,:,:) = 0.0
    socrates_derived_fields%flux_net_lw(:,:,:) = 0.0
    socrates_derived_fields%toa_up_longwave(:,:) = 0.0
    socrates_derived_fields%toa_down_shortwave(:,:) = 0.0
    socrates_derived_fields%toa_up_shortwave(:,:) = 0.0
    socrates_derived_fields%surface_up_longwave(:,:) = 0.0
    socrates_derived_fields%surface_down_longwave(:,:) = 0.0
    socrates_derived_fields%surface_down_shortwave(:,:) = 0.0
    socrates_derived_fields%surface_up_shortwave(:,:) = 0.0
    socrates_derived_fields%swrad_hr(:,:,:) = 0.0
    socrates_derived_fields%lwrad_hr(:,:,:) = 0.0
    socrates_derived_fields%totrad_hr(:,:,:) = 0.0

    ! derive density and radiation factor for heating rate calculation
    socrates_derived_fields%density_factor(1) =  0.0 
    do k = 2, k_top
       socrates_derived_fields%density_factor(k) =  & 
            current_state%global_grid%configuration%vertical%rhon(k)* &
            current_state%global_grid%configuration%vertical%dz(k)
    enddo
    socrates_derived_fields%density_factor(1) = socrates_derived_fields%density_factor(2)

    socrates_derived_fields%radiation_factor(2:k_top) =   &
         1.0/(socrates_derived_fields%density_factor(2:k_top)*cp)
    socrates_derived_fields%radiation_factor(1) = socrates_derived_fields%radiation_factor(2)

    call set_and_test_socrates_monc_options(current_state, socrates_opt)

    ! set up the switches and names for Edwards-Slingo code
    call sw_input(current_state)
    ! DEPENDS ON: read_spectrum
    call read_spectrum( sw_control%spectral_file,                  &
         sw_spectrum )
    ! DEPENDS ON: compress_spectrum
    CALL compress_spectrum( sw_control, sw_spectrum )

    call lw_input(current_state, lw_control)
    call read_spectrum( lw_control%spectral_file,                  &
                      lw_spectrum )
    CALL compress_spectrum( lw_control, lw_spectrum )

    ! Read the mcc_profile
    call read_mcclatchey_profiles(current_state, mcc)
    !
    ! Calculate the number of radiation levels with
    ! mcc_levs combined to monc model levels
    call calculate_radiation_levels(current_state, mcc)
    !
    ! allocate fields to pass to socrates
    call allocate_merge_data_fields(current_state, merge_fields, mcc)

  end subroutine initialisation_callback

  !> Called for each column per timestep this will apply a forcing term
  !> to the aerosol fields
  !! @param current_state The current model state
  subroutine timestep_callback(current_state)
    type(model_state_type), target, intent(inout) :: current_state

    real(DEFAULT_PRECISION), parameter :: &
         snotset  = -9999.       ! A value of sec_out meaing darkness
    real(DEFAULT_PRECISION) :: local_dtm  ! Local timestep variable
    integer :: icol, jcol ! Shorthand column indices
    integer :: target_x_index, target_y_index
   
    integer :: k ! look counter

    ! No need to do radiation calculations in the halos or on the first timestep
    !
    if (current_state%halo_column .or. current_state%timestep < 2) return

    local_dtm = current_state%dtm*2.0
    if (current_state%field_stepping == FORWARD_STEPPING) local_dtm=current_state%dtm
    
    ! work out column indexes from MONC (these include the halo)
    icol=current_state%column_local_x 
    jcol=current_state%column_local_y
    ! work out a target index for radiation arrays (no halo)
    target_y_index=jcol-current_state%local_grid%halo_size(Y_INDEX)
    target_x_index=icol-current_state%local_grid%halo_size(X_INDEX)
    
    ! Test whether it is a radiation calc timestep on the first non-halo column
    ! If it is, then calc all year, day, time in hours and timestep.
    ! Note: all socrates time control variables are declared in the socrates_opt
    !       or socrates_derived_fields structure, except rad_last_time, which is 
    !       in current_state. 
    if (current_state%first_nonhalo_timestep_column) then
       !i) 1 call radiation on timestep 2 to initialise the heating rates
       !ii) if rad_int_time less than or equal to 0.0, SOCRATES called on every timestep
       if (socrates_opt%rad_int_time <= 0.0 .or. current_state%timestep == 2 ) then
          socrates_opt%l_rad_calc = .true.
       else
          socrates_opt%l_rad_calc = &
               ((current_state%time - &
               (current_state%rad_last_time + socrates_opt%rad_int_time)) > 0.0)
       endif
    endif

    if (socrates_opt%l_rad_calc) then
       if (current_state%first_nonhalo_timestep_column) then
          call log_master_log &
               (LOG_INFO, "Socrates called, time ="//trim(conv_to_string(current_state%time))//&
               " rad_int_time="//trim(conv_to_string(socrates_opt%rad_int_time))//&
               " local dtm="//trim(conv_to_string(local_dtm)))
          call log_master_log &
               (LOG_INFO, "methane ="//trim(conv_to_string(socrates_opt%ch4_mmr))//&
               " l_ch4="//trim(conv_to_string(lw_control%l_ch4)))
          
          ! Do not really like this but update the rad_time_hours and
          ! rad_day here using time.
          socrates_opt%rad_day = socrates_opt%rad_start_day + &
               int((socrates_opt%rad_start_time + (current_state%time/3600.0))/24.0)
          socrates_opt%rad_time_hours =   &
               (socrates_opt%rad_start_time + ((current_state%time+local_dtm)/3600.0))  &
               - (24.0*(socrates_opt%rad_day-socrates_opt%rad_start_day))
          ! set the radiation timestep
          if (socrates_opt%rad_int_time == 0) then
             socrates_derived_fields%dt_secs = local_dtm
          else
             socrates_derived_fields%dt_secs = socrates_opt%rad_int_time
          endif
          ! Finally,  if we will calculate radiative fluxes on this timestep update the last
          ! radiative timetep to this one.
          current_state%rad_last_time = current_state%time
       endif
       
       ! set surface temperature. Based on LEM, but needs some thought to capture
       ! atmospheric surface and actual surface. Here they are the same
       ! which may lead to too much emission (depending on profile)
       !
       if (current_state%use_surface_boundary_conditions) then
          if (current_state%type_of_surface_boundary_conditions == PRESCRIBED_SURFACE_VALUES) then
             socrates_derived_fields%srf_temperature = &
                  current_state%theta_surf /           &
                  ((current_state%surface_reference_pressure/current_state%surface_pressure)**r_over_cp)
          elseif (current_state%type_of_surface_boundary_conditions == PRESCRIBED_SURFACE_FLUXES) then
             !THIS IS NOT CORRECT - SHOULD WORK OUT TEMP FROM FLUXES
             socrates_derived_fields%srf_temperature = &
                  (current_state%th%data(2, jcol, icol)        &
                  + current_state%global_grid%configuration%vertical%thref(2))   &
                  * current_state%global_grid%configuration%vertical%rprefrcp(2)
          end if
         !else if(DO SOMETHING WITH JULES OUTPUT WHEN AVAILABLE)
       else
          socrates_derived_fields%srf_temperature = &
               (current_state%th%data(2, jcol, icol)        &
               + current_state%global_grid%configuration%vertical%thref(2))   &
               * current_state%global_grid%configuration%vertical%rprefrcp(2)
       end if

       ! merge_data takes t,qv, ql, qi and then returns t_socrates, qv_socrates
       ! ql_socrates, qi_socrates (all declared at head of module so no need
       ! to pass as arguments
       !
       call merge_data(current_state, socrates_opt, socrates_derived_fields, merge_fields, mcc )
       !
       ! intialise variables for rad_calc just to be on safe side
       socrates_derived_fields%sindec = 0.0
       socrates_derived_fields%scs    = 1.0
       socrates_derived_fields%fraction_lit  = 0.0
       socrates_derived_fields%cosz   = snotset

       if (socrates_opt%l_solar_fixed) then
          socrates_derived_fields%sol_const = socrates_opt%solar_fixed
          socrates_derived_fields%sec_out = socrates_opt%sec_fixed
          socrates_derived_fields%fraction_lit = 1.0
       else
          ! this routine uses the rad_start_day, rad_start_year, and l_360
          ! to work out position of sun. Returns sindec and scs to use in
          ! solar angle calculation
          call solar_pos_calculation(socrates_opt, socrates_derived_fields)
          !
          ! set the solar constant by scaling default constant by scs from
          ! solar_pos_calculation
          socrates_derived_fields%sol_const= &
               socrates_opt%default_solar_constant * &
               socrates_derived_fields%scs
          
          ! calulates the solar angle, fraction_lit and cos of zenith angle
          call solar_angle_calculation(socrates_opt, socrates_derived_fields)
          !
          ! set fraction_lit to 1 or 0 depending, i.e. grid is lit or not, there
          ! is no partial lighting
          if (socrates_derived_fields%fraction_lit .ne. 0.0) then
             socrates_derived_fields%fraction_lit = 1.0
             socrates_derived_fields%sec_out =  &
                  1./socrates_derived_fields%cosz
          else
             ! seems a silly value to set sec_out to check as it will break
             ! the variable albedo calc
             socrates_derived_fields%sec_out = snotset
          endif
       endif
       !
       ! Set surface albedo for rad_calc
       !
       if (socrates_derived_fields%fraction_lit .ne. 0.0) then
          if (socrates_opt%l_variable_srf_albedo) then
             stop
!!$              ! albfac1=0.026 albfac2=1.7 albfac3=0.065
!!$              ! albfac4=0.15  albfac5=0.1 albfac6=0.5
!!$              ! albfac7=1.0   alb2_var=0.06
!!$              cosz = 1./sec_out
!!$              albedoin1 = albfac1/(cosz**albfac2 + albfac3)
!!$     &                  + albfac4*(cosz - albfac5)*
!!$     &                   (cosz - albfac6)*(cosz - albfac7)
!!$              albedoin2 = alb2_var
          else
             socrates_derived_fields%albedoin1 = socrates_opt%surface_albedo
             socrates_derived_fields%albedoin2 = socrates_opt%surface_albedo
          endif
       endif
       
       ! AH - after all this testing check whether solar is required. If 
       ! no solar then set fraction_lit = 0.0
       if (socrates_opt%l_no_solar) then 
          socrates_derived_fields%fraction_lit = 0.0
       endif

       call rad_ctl(current_state, sw_spectrum, lw_spectrum,    &
             mcc, socrates_opt, merge_fields, socrates_derived_fields)
       
       ! This is needed for JULES coupling. Including irrespective of JULES enabled
       ! assign downward fluxes at the surface
       !current_state%sw_down_surf(jcol, icol) = &
       !  socrates_derived_fields%flux_down_sw(1, target_y_index, target_x_index)
       !current_state%lw_down_surf(jcol, icol) = &
       !  socrates_derived_fields%flux_down_lw(1, target_y_index, target_x_index)
       
    endif

    ! update the current_state sth
    ! AH - temporary code to check bit comparison of socrates_couple
    !current_state%sth_lw%data(:, jcol, icol) = 0.0
    !current_state%sth_sw%data(:, jcol, icol) = 0.0
    ! AH - end temporary code
    current_state%sth%data(:, jcol, icol) = &
         current_state%sth%data(:, jcol, icol) +  &
         current_state%sth_lw%data(:, jcol, icol) + &
         current_state%sth_sw%data(:, jcol, icol)

  end subroutine timestep_callback

  subroutine finalisation_callback(current_state)
    
    type(model_state_type), target, intent(inout) :: current_state

!    deallocate(merge_fields%t, merge_fields%qv, &
!         ql_socrates, qi_socrates)

  end subroutine finalisation_callback

!> Field information retrieval callback, this returns information for a specific components published field
!! @param current_state Current model state
!! @param name The name of the field to retrieve information for
!! @param field_information Populated with information about the field
  subroutine field_information_retrieval_callback(current_state, name, field_information)
    type(model_state_type), target, intent(inout) :: current_state
   character(len=*), intent(in) :: name
   type(component_field_information_type), intent(out) :: field_information

   if ( name .eq. "flux_up_shortwave" .or. name .eq. "flux_down_shortwave"  .or.  &
        name .eq. "flux_up_longwave" .or. name .eq. "flux_down_longwave" .or.     &
        name .eq. "shortwave_heating_rate" .or. name .eq. "longwave_heating_rate" &
        .or. name .eq. "total_radiative_heating_rate") then
      field_information%field_type=COMPONENT_ARRAY_FIELD_TYPE
      field_information%number_dimensions=3
      field_information%dimension_sizes(1)=current_state%local_grid%size(Z_INDEX)
      field_information%dimension_sizes(2)=current_state%local_grid%size(Y_INDEX)
      field_information%dimension_sizes(3)=current_state%local_grid%size(X_INDEX)
      field_information%data_type=COMPONENT_DOUBLE_DATA_TYPE
   else
      field_information%field_type=COMPONENT_ARRAY_FIELD_TYPE
      field_information%number_dimensions=2
      field_information%dimension_sizes(1)=current_state%local_grid%size(Y_INDEX)
      field_information%dimension_sizes(2)=current_state%local_grid%size(X_INDEX)
      field_information%data_type=COMPONENT_DOUBLE_DATA_TYPE
   endif

   field_information%enabled=.true.

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

   ! 3D radiative flux and heating rates
   if (name .eq. "flux_up_shortwave") then
      allocate(field_value%real_3d_array(current_state%local_grid%size(Z_INDEX),  &
           current_state%local_grid%size(Y_INDEX),                                &
           current_state%local_grid%size(X_INDEX)))
      field_value%real_3d_array(:,:,:) = socrates_derived_fields%flux_up_sw(:,:,:)
   else if (name .eq. "flux_down_shortwave") then
      allocate(field_value%real_3d_array(current_state%local_grid%size(Z_INDEX),  &
           current_state%local_grid%size(Y_INDEX),                                &
           current_state%local_grid%size(X_INDEX)))
      field_value%real_3d_array(:,:,:) = socrates_derived_fields%flux_down_sw(:,:,:)
   else if (name .eq. "flux_up_longwave") then
      allocate(field_value%real_3d_array(current_state%local_grid%size(Z_INDEX),  &
           current_state%local_grid%size(Y_INDEX),                                &
           current_state%local_grid%size(X_INDEX)))
      field_value%real_3d_array(:,:,:) = socrates_derived_fields%flux_up_lw(:,:,:)
   else if (name .eq. "flux_down_longwave") then
      allocate(field_value%real_3d_array(current_state%local_grid%size(Z_INDEX),  &
           current_state%local_grid%size(Y_INDEX),                                &
           current_state%local_grid%size(X_INDEX)))
      field_value%real_3d_array(:,:,:) = socrates_derived_fields%flux_down_lw(:,:,:)
   else if (name .eq. "shortwave_heating_rate") then
      allocate(field_value%real_3d_array(current_state%local_grid%size(Z_INDEX),  &
           current_state%local_grid%size(Y_INDEX),                                &
           current_state%local_grid%size(X_INDEX)))
      field_value%real_3d_array(:,:,:) = socrates_derived_fields%swrad_hr(:,:,:)
   else if (name .eq. "longwave_heating_rate") then
      allocate(field_value%real_3d_array(current_state%local_grid%size(Z_INDEX),  &
           current_state%local_grid%size(Y_INDEX),                                &
           current_state%local_grid%size(X_INDEX)))
      field_value%real_3d_array(:,:,:) = socrates_derived_fields%lwrad_hr(:,:,:)
   else if (name .eq. "total_radiative_heating_rate") then
      allocate(field_value%real_3d_array(current_state%local_grid%size(Z_INDEX),  &
           current_state%local_grid%size(Y_INDEX),                                &
           current_state%local_grid%size(X_INDEX)))
      field_value%real_3d_array(:,:,:) = socrates_derived_fields%totrad_hr(:,:,:)
   !
   ! 2D radiative fluxes   
   else if (name .eq. "toa_up_longwave") then
      allocate(field_value%real_2d_array(current_state%local_grid%size(Y_INDEX),   &
           current_state%local_grid%size(X_INDEX)))
      field_value%real_2d_array(:,:) = socrates_derived_fields%toa_up_longwave(:,:)
   else if (name .eq. "surface_down_longwave") then
      allocate(field_value%real_2d_array(current_state%local_grid%size(Y_INDEX),   &
           current_state%local_grid%size(X_INDEX)))
      field_value%real_2d_array(:,:) = socrates_derived_fields%surface_down_longwave(:,:) 
   else if (name .eq. "surface_up_longwave") then
      allocate(field_value%real_2d_array(current_state%local_grid%size(Y_INDEX),   &
           current_state%local_grid%size(X_INDEX)))
      field_value%real_2d_array(:,:) = socrates_derived_fields%surface_up_longwave(:,:)
   else if (name .eq. "toa_up_shortwave") then
      allocate(field_value%real_2d_array(current_state%local_grid%size(Y_INDEX),   &
           current_state%local_grid%size(X_INDEX)))
      field_value%real_2d_array(:,:) = socrates_derived_fields%toa_up_shortwave(:,:)
   else if (name .eq. "toa_down_shortwave") then
      allocate(field_value%real_2d_array(current_state%local_grid%size(Y_INDEX),   &
           current_state%local_grid%size(X_INDEX)))
      field_value%real_2d_array(:,:) = socrates_derived_fields%toa_down_shortwave(:,:)
   else if (name .eq. "surface_down_shortwave") then
      allocate(field_value%real_2d_array(current_state%local_grid%size(Y_INDEX),   &
           current_state%local_grid%size(X_INDEX)))
      field_value%real_2d_array(:,:) = socrates_derived_fields%surface_down_shortwave(:,:) 
   else if (name .eq. "surface_up_shortwave") then
      allocate(field_value%real_2d_array(current_state%local_grid%size(Y_INDEX),   &
           current_state%local_grid%size(X_INDEX)))
      field_value%real_2d_array(:,:) = socrates_derived_fields%surface_up_shortwave(:,:) 
   end if

 end subroutine field_value_retrieval_callback

end module socrates_couple_mod
