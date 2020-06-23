!> Implimentation of CASIM microphysics
module casim_mod
  use datadefn_mod, only : DEFAULT_PRECISION
  use monc_component_mod, only : component_descriptor_type, &
       COMPONENT_ARRAY_FIELD_TYPE, COMPONENT_DOUBLE_DATA_TYPE, &
       component_field_value_type, component_field_information_type
  use state_mod, only : FORWARD_STEPPING, model_state_type
  use grids_mod, only : Z_INDEX, Y_INDEX, X_INDEX
  use science_constants_mod
  use q_indices_mod, only: get_q_index, standard_q_names
  use optionsdatabase_mod, only : options_get_real, options_get_logical, options_get_integer
  use logging_mod, only : LOG_INFO, LOG_ERROR, log_master_log
  use registry_mod, only : is_component_enabled

  ! casim modules...
  use variable_precision, ONLY: wp
  use initialize, only: mphys_init
  use mphys_switches, only: set_mphys_switches,  &
     l_warm,                                     &
     nq_l, nq_r, nq_i, nq_s, nq_g,               &
     l_2mc, l_2mr, l_2mi, l_2ms, l_2mg,          &
     l_3mr, l_3ms, l_3mg,                        &
     soluble_modes, active_cloud, active_rain,   &
     insoluble_modes, active_ice, active_number, &
     isol, iinsol, option, aerosol_option,        &
     iopt_act, option, aerosol_option        &
     ! All process rates are read in using read_configuration
     ! routine in this component
     , l_aaut, l_aacc, l_aevp, l_ased, l_warm            &
     , l_inuc, iopt_rcrit, iopt_inuc, l_iaut, l_iacw                &
     , l_rain, l_boussinesq, diag_mu_option   &
     , l_sed_3mdiff, l_cons, l_abelshipway, l_sed_icecloud_as_1m         &
     , l_active_inarg2000, process_level, l_separate_rain, l_idep    &
     , max_step_length, max_sed_length, l_sg, l_g, l_passive        &
     , l_passive3m, l_limit_psd, l_override_checks &
     , max_mu, fix_mu, l_raci_g, l_onlycollect, l_inhom_revp &
     , l_tidy_conserve_E , l_tidy_conserve_q & 
     , l_pcond & ! Condensation	 
     , l_praut & ! Autoconversion cloud -> rain
     , l_pracw & ! Accretion  cloud -> rain
     , l_pracr & ! aggregation of rain drops
     , l_prevp & ! evaporation of rain
     , l_psedl & ! sedimentation of cloud
     , l_psedr & ! sedimentation of rain
     , l_ptidy & ! tidying term 1
     , l_ptidy2 & ! tidying term 2
     , l_pinuc & ! ice nucleation
     , l_pidep & ! ice deposition
     , l_piacw & ! ice accreting water
     , l_psaut & ! ice autoconversion ice -> snow
     , l_psdep & ! vapour deposition onto snow
     , l_psacw & ! snow accreting water
     , l_pgdep & ! vapour deposition onto graupel
     , l_pseds & ! snow sedimentation
     , l_psedi & ! ice sedimentation
     , l_psedg & ! graupel sedimentation
     , l_psaci & ! snow accreting ice
     , l_praci & ! rain accreting ice
     , l_psacr & ! snow accreting rain
     , l_pgacr & ! graupel accreting rain
     , l_pgacw & ! graupel accreting cloud water
     , l_pgaci & ! graupel accreting ice
     , l_pgacs & ! graupel accreting snow
     , l_piagg & ! aggregation of ice particles
     , l_psagg & ! aggregation of snow particles
     , l_pgagg & ! aggregation of graupel particles
     , l_psbrk & ! break up of snow flakes
     , l_pgshd & ! shedding of liquid from graupel
     , l_pihal & ! hallet mossop
     , l_psmlt & ! snow melting
     , l_pgmlt & ! graupel melting
     , l_phomr & ! homogeneous freezing of rain
     , l_phomc & ! homogeneous freezing of cloud droplets
     , l_pssub & ! sublimation of snow
     , l_pgsub & ! sublimation of graupel
     , l_pisub & ! sublimation of ice
     , l_pimlt & ! ice melting
     ! New switches for sedimentation (these are sort-of temporary)
     , l_gamma_online & ! when true use standard vn0.3.3 sed, when false use precalced gamma
     , l_subseds_maxv & ! Use a CFL criteria based on max terminal velocity
                    ! and sed_1M_2M 
     , l_sed_eulexp & ! switch for eulexp sed based on UM. Default is false
                      ! so standard casim sed used
     , cfl_vt_max & ! cfl limit for sedimentation (default = 1.0)
     , l_kfsm

  use micro_main, only: shipway_microphysics
  use generic_diagnostic_variables, ONLY: casdiags, allocate_diagnostic_space, &
       deallocate_diagnostic_space
  use casim_monc_dgs_space, only: casim_monc_dgs, allocate_casim_monc_dgs_space, &
       populate_casim_monc_dg

  implicit none

#ifndef TEST_MODE
  private
#endif

  REAL(wp), allocatable :: theta(:,:,:), pressure(:,:,:),  &
     z_half(:,:,:), z_centre(:,:,:), dz(:,:,:), qv(:,:,:),qc(:,:,:) &
     , nc(:,:,:), qr(:,:,:), nr(:,:,:), m3r(:,:,:),rho(:,:,:) &
     , exner(:,:,:), w(:,:,:), tke(:,:,:)                               &
     , qi(:,:,:), ni(:,:,:), qs(:,:,:), ns(:,:,:), m3s(:,:,:) &
     , qg(:,:,:), ng(:,:,:), m3g(:,:,:) 

  REAL(wp), allocatable :: AccumSolMass(:,:,:), AccumSolNumber(:,:,:) ! Accumulation mode aerosol
  REAL(wp), allocatable :: ActiveSolLiquid(:,:,:)                      ! Activated aerosol
  REAL(wp), allocatable :: AitkenSolMass(:,:,:), AitkenSolNumber(:,:,:) ! Aitken mode aerosol
  REAL(wp), allocatable :: CoarseSolMass(:,:,:), CoarseSolNumber(:,:,:) ! Course mode aerosol
  REAL(wp), allocatable :: ActiveSolRain(:,:,:)                      ! Activeated aerosol in rain
  REAL(wp), allocatable :: CoarseDustMass(:,:,:), CoarseDustNumber(:,:,:) ! Coarse Dust
  REAL(wp), allocatable :: ActiveInsolIce(:,:,:)                      ! Activeated dust
  REAL(wp), allocatable :: ActiveSolIce(:,:,:)                      ! Activeated aerosol in ice
  REAL(wp), allocatable :: ActiveInsolLiquid(:,:,:)                      ! Activeated dust in cloud
  REAL(wp), allocatable :: AccumInsolMass(:,:,:)                      ! Accum mode dust mass
  REAL(wp), allocatable :: AccumInsolNumber(:,:,:)                      ! Accum mode dust number
  REAL(wp), allocatable :: ActiveSolNumber(:,:,:)                      ! Activated soluble number (if we need a tracer)
  REAL(wp), allocatable :: ActiveInsolNumber(:,:,:)                      ! Activated insoluble number (if we need a tracer)
  !!AH - these are dummy arrays to accommodate the hygroscopicity from
  !!     UKCA-MODE. At present these are only required so that the shipway_microphysics
  !!     argument list matches between the UM and MONC interface. These may become non-dummy 
  !!     arguments in the future
  REAL(wp), allocatable :: AitkenSolBk(:,:,:)
  REAL(wp), allocatable :: AccumSolBk(:,:,:) 
  REAL(wp), allocatable :: CoarseSolBk(:,:,:)

  ! Tendency from other physics/advection/forcing
  ! NB This is tendency (/s) not an increment over the timestep
  REAL(wp), allocatable :: dqv(:,:,:), dth(:,:,:), dqc(:,:,:), dnc(:,:,:)              &
     , dqr(:,:,:), dnr(:,:,:), dm3r(:,:,:)                                  &
     , dqi(:,:,:), dni(:,:,:), dqs(:,:,:), dns(:,:,:), dm3s(:,:,:)      & 
     , dqg(:,:,:), dng(:,:,:), dm3g(:,:,:) 

  REAL(wp), allocatable :: dAccumSolMass(:,:,:), dAccumSolNumber(:,:,:) ! Accumulation mode aerosol
  REAL(wp), allocatable :: dActiveSolLiquid(:,:,:)                       ! Activated aerosol
  REAL(wp), allocatable :: dAitkenSolMass(:,:,:), dAitkenSolNumber(:,:,:) ! Aitken mode aerosol
  REAL(wp), allocatable :: dCoarseSolMass(:,:,:), dCoarseSolNumber(:,:,:) ! Course mode aerosol
  REAL(wp), allocatable :: dActiveSolRain(:,:,:)                       ! Activeated aerosol in rain
  REAL(wp), allocatable :: dCoarseDustMass(:,:,:), dCoarseDustNumber(:,:,:) ! Dust
  REAL(wp), allocatable :: dActiveInsolIce(:,:,:)                       ! Activeated dust
  REAL(wp), allocatable :: dActiveSolIce(:,:,:)                       ! Activeated aerosol in ice
  REAL(wp), allocatable :: dActiveInsolLiquid(:,:,:)                       ! Activeated dust in cloud
  REAL(wp), allocatable :: dAccumInsolMass(:,:,:)                      ! Accum mode dust mass
  REAL(wp), allocatable :: dAccumInsolNumber(:,:,:)                      ! Accum mode dust number
  REAL(wp), allocatable :: dActiveSolNumber(:,:,:)                      ! Activated soluble number (if we need a tracer)
  REAL(wp), allocatable :: dActiveInsolNumber(:,:,:)                      ! Activated insoluble number (if we need a tracer)

  REAL(wp), allocatable :: surface_precip(:,:)

  INTEGER ::   ils,ile, jls,jle, kls,kle, &
     its,ite, jts,jte, kts,kte

  INTEGER :: iqv=0, iql=0, iqr=0,  iqi=0, iqs=0,  iqg=0
  INTEGER ::        inl=0, inr=0,  ini=0, ins=0,  ing=0
  INTEGER ::               i3mr=0,        i3ms=0, i3mg=0 
  INTEGER ::                &
     i_AccumSolMass=0,      &
     i_AccumSolNumber=0,    &
     i_ActiveSolLiquid=0,   &
     i_AitkenSolMass=0,     &
     i_AitkenSolNumber=0,   &
     i_CoarseSolMass=0,     &
     i_CoarseSolNumber=0,   &
     i_ActiveSolRain=0,     &
     i_CoarseDustMass=0,    &
     i_CoarseDustNumber=0,  &
     i_ActiveInsolIce=0,    &
     i_ActiveSolIce=0,      &
     i_ActiveInsolLiquid=0, &
     i_AccumInsolMass=0,    &
     i_AccumInsolNumber=0,  &
     i_ActiveSolNumber=0,   &
     i_ActiveInsolNumber=0

  real(kind=DEFAULT_PRECISION), dimension(:), allocatable ::     &
       phomc_tot, pinuc_tot, pidep_tot, psdep_tot, piacw_tot, psacw_tot, psacr_tot, pisub_tot,   &
       pssub_tot, pimlt_tot, psmlt_tot, psaut_tot, psaci_tot, praut_tot, pracw_tot, prevp_tot,   &
       pgacw_tot, pgacs_tot, pgmlt_tot, pgsub_tot, psedi_tot, pseds_tot, psedr_tot, psedg_tot,   &
       psedl_tot, pcond_tot
  
  public casim_get_descriptor
contains

  !> Provides the descriptor back to the caller and is used in component registration
  !! @returns The termination check component descriptor
  type(component_descriptor_type) function casim_get_descriptor()
    casim_get_descriptor%name="casim"
    casim_get_descriptor%version=0.1
    casim_get_descriptor%initialisation=>initialisation_callback
    casim_get_descriptor%timestep=>timestep_callback

    ! Set up fields to be published for diagnostics. These are made available to
    ! IO server in the field_value_retrieval_callback at the end of this module
    !
    casim_get_descriptor%field_value_retrieval=>field_value_retrieval_callback
    casim_get_descriptor%field_information_retrieval=>field_information_retrieval_callback
    
    allocate(casim_get_descriptor%published_fields(27))
    
    casim_get_descriptor%published_fields(1)="surface_precip"    
    casim_get_descriptor%published_fields(2)="homogeneous_freezing_rate"
    casim_get_descriptor%published_fields(3)="ice_nucleations_rate"
    casim_get_descriptor%published_fields(4)="ice_deposition_rate"
    casim_get_descriptor%published_fields(5)="snow_deposition_rate"
    casim_get_descriptor%published_fields(6)="ice_acc_cloud_rate"
    casim_get_descriptor%published_fields(7)="snow_acc_cloud_rate"
    casim_get_descriptor%published_fields(8)="snow_acc_rain_rate"
    casim_get_descriptor%published_fields(9)="ice_sublime_rate"
    casim_get_descriptor%published_fields(10)="snow_sublime_rate"
    casim_get_descriptor%published_fields(11)="ice_melt_rate"
    casim_get_descriptor%published_fields(12)="snow_melt_rate"
    casim_get_descriptor%published_fields(13)="snow_autoconversion_rate"
    casim_get_descriptor%published_fields(14)="snow_acc_ice_rate"
    casim_get_descriptor%published_fields(15)="rain_autoconversion_rate"
    casim_get_descriptor%published_fields(16)="rain_acc_cloud_rate"
    casim_get_descriptor%published_fields(17)="rain_evap_rate"
    casim_get_descriptor%published_fields(18)="graup_acc_cloud_rate"
    casim_get_descriptor%published_fields(19)="graup_acc_snow_rate"
    casim_get_descriptor%published_fields(20)="graup_melt_rate"
    casim_get_descriptor%published_fields(21)="graup_sublime_rate"
    casim_get_descriptor%published_fields(22)="ice_sed_rate"
    casim_get_descriptor%published_fields(23)="snow_sed_rate"
    casim_get_descriptor%published_fields(24)="rain_sed_rate"
    casim_get_descriptor%published_fields(25)="graup_sed_rate"
    casim_get_descriptor%published_fields(26)="cloud_sed_rate"
    casim_get_descriptor%published_fields(27)="condensation_rate"
    
  end function casim_get_descriptor

  !> The initialisation callback sets up the microphysics
  !! @param current_state The current model state
  subroutine initialisation_callback(current_state)
    type(model_state_type), target, intent(inout) :: current_state

    integer :: y_size_local, x_size_local

    if (is_component_enabled(current_state%options_database, "simplecloud")) then
      call log_master_log(LOG_ERROR, "Casim and Simplecloud are enabled, this does not work yet. Please disable one")
   end if

   !allocate(psedl_tot(current_state%local_grid%size(Z_INDEX)), &
   !      pcond_tot(current_state%local_grid%size(Z_INDEX)))
    
    y_size_local = current_state%local_grid%size(Y_INDEX)
    x_size_local = current_state%local_grid%size(X_INDEX)
    
    call read_configuration(current_state)   

    ils=1
    ile=1
    jls=1
    jle=1
    kls=2
    kle=current_state%local_grid%size(Z_INDEX)
    its=1
    ite=1
    jts=1
    jte=1
    kts=1
    kte=current_state%local_grid%size(Z_INDEX)

    !> Set up and allocate the local arrays

    allocate(pressure(kte,1,1))
    allocate(z_half(0:kte,1,1))
    allocate(z_centre(kte,1,1))
    allocate(dz(kte,1,1))
    allocate(rho(kte,1,1))
    allocate(exner(kte,1,1))
    allocate(w(kte,1,1))
    allocate(tke(kte,1,1))

    allocate(theta(kte,1,1))
    allocate(qv(kte,1,1))
    allocate(qc(kte,1,1))
    allocate(nc(kte,1,1))
    allocate(qr(kte,1,1))
    allocate(nr(kte,1,1))
    allocate(m3r(kte,1,1))
    allocate(qi(kte,1,1))
    allocate(ni(kte,1,1))
    allocate(qs(kte,1,1))
    allocate(ns(kte,1,1))
    allocate(m3s(kte,1,1))
    allocate(qg(kte,1,1))
    allocate(ng(kte,1,1))
    allocate(m3g(kte,1,1))

    allocate(AccumSolMass(kte,1,1))
    allocate(AccumSolNumber(kte,1,1))
    allocate(ActiveSolLiquid(kte,1,1))
    allocate(AitkenSolMass(kte,1,1))
    allocate(AitkenSolNumber(kte,1,1))
    allocate(CoarseSolMass(kte,1,1))
    allocate(CoarseSolNumber(kte,1,1))
    allocate(ActiveSolRain(kte,1,1))
    allocate(CoarseDustMass(kte,1,1))
    allocate(CoarseDustNumber(kte,1,1))
    allocate(ActiveInsolIce(kte,1,1))
    allocate(ActiveSolIce(kte,1,1))
    allocate(ActiveInsolLiquid(kte,1,1))
    allocate(AccumInsolMass(kte,1,1))
    allocate(AccumInsolNumber(kte,1,1))
    allocate(ActiveSolNumber(kte,1,1))
    allocate(ActiveInsolNumber(kte,1,1))
    ! allocate the hygoscopicity arrays
    allocate(AitkenSolBk(kte,1,1))
    allocate(AccumSolBk(kte,1,1))
    allocate(CoarseSolBk(kte,1,1))

    allocate(dth(kte,1,1))
    allocate(dqv(kte,1,1))
    allocate(dqc(kte,1,1))
    allocate(dnc(kte,1,1))
    allocate(dqr(kte,1,1))
    allocate(dnr(kte,1,1))
    allocate(dm3r(kte,1,1))
    allocate(dqi(kte,1,1))
    allocate(dni(kte,1,1))
    allocate(dqs(kte,1,1))
    allocate(dns(kte,1,1))
    allocate(dm3s(kte,1,1))
    allocate(dqg(kte,1,1))
    allocate(dng(kte,1,1))
    allocate(dm3g(kte,1,1))

    allocate(dAccumSolMass(kte,1,1))
    allocate(dAccumSolNumber(kte,1,1))
    allocate(dActiveSolLiquid(kte,1,1))
    allocate(dAitkenSolMass(kte,1,1))
    allocate(dAitkenSolNumber(kte,1,1))
    allocate(dCoarseSolMass(kte,1,1))
    allocate(dCoarseSolNumber(kte,1,1))
    allocate(dActiveSolRain(kte,1,1))
    allocate(dCoarseDustMass(kte,1,1))
    allocate(dCoarseDustNumber(kte,1,1))
    allocate(dActiveInsolIce(kte,1,1))
    allocate(dActiveSolIce(kte,1,1))
    allocate(dActiveInsolLiquid(kte,1,1))
    allocate(dAccumInsolMass(kte,1,1))
    allocate(dAccumInsolNumber(kte,1,1))
    allocate(dActiveSolNumber(kte,1,1))
    allocate(dActiveInsolNumber(kte,1,1))

    call set_mphys_switches(option,aerosol_option)
    call mphys_init(its, ite, jts, jte, kts, kte, ils, ile, jls, jle, kls, kle, l_tendency=.true.)

    ! Need to allocate the appropriate indices, e.g. iqv, iql...
    ! This needs to be compatible with the rest of the model
    ! This essentially reproduces the switching in the main microphysics
    ! code already done above (set_mphys_switches), so could be combined 
    ! once the MONC method has been finalized.
    ! Note the numbers assigned here may be different from those assigned 
    ! in the microphysics since we share the q array with other components. 

    if (.not. allocated(current_state%cq))then
       allocate(current_state%cq(current_state%number_q_fields))
       current_state%cq=0.0_DEFAULT_PRECISION
    end if

    ! Mass
    iqv = get_q_index(standard_q_names%VAPOUR, 'casim')
    if (nq_l>0)then
       iql = get_q_index(standard_q_names%CLOUD_LIQUID_MASS, 'casim')
       current_state%cq(iql) = -1.0
    end if
    if (nq_r>0)then
       iqr = get_q_index(standard_q_names%RAIN_MASS, 'casim')
       current_state%rain_water_mixing_ratio_index=iqr
      current_state%cq(iqr) = -1.0
    end if
    if (.not. l_warm)then
      if (nq_i>0)then
         iqi = get_q_index(standard_q_names%ICE_MASS, 'casim')
         current_state%ice_water_mixing_ratio_index=iqi
        current_state%cq(iqi) = -1.0
      end if
      if (nq_s>0)then
         iqs = get_q_index(standard_q_names%SNOW_MASS, 'casim')
         current_state%snow_water_mixing_ratio_index=iqs
        current_state%cq(iqs) = -1.0
      end if
      if (nq_g>0)then
         iqg = get_q_index(standard_q_names%GRAUPEL_MASS, 'casim')
         current_state%graupel_water_mixing_ratio_index=iqg
        current_state%cq(iqg) = -1.0
      end if
    end if

    ! Number
    if (l_2mc)inl = get_q_index(standard_q_names%CLOUD_LIQUID_NUMBER, 'casim')
    if (l_2mr)inr = get_q_index(standard_q_names%RAIN_NUMBER, 'casim')
    if (.not. l_warm)then
      if (l_2mi)ini = get_q_index(standard_q_names%ICE_NUMBER, 'casim')
      if (l_2ms)ins = get_q_index(standard_q_names%SNOW_NUMBER, 'casim')
      if (l_2mg)ing = get_q_index(standard_q_names%GRAUPEL_NUMBER, 'casim')
    end if

    ! Third moments
    if (l_3mr)i3mr = get_q_index(standard_q_names%RAIN_THIRD_MOMENT, 'casim')
    if (.not. l_warm)then
      if (l_3ms)i3ms = get_q_index(standard_q_names%SNOW_THIRD_MOMENT, 'casim')
      if (l_3mg)i3mg = get_q_index(standard_q_names%GRAUPEL_THIRD_MOMENT, 'casim')
    end if

    ! Aerosol
    if (soluble_modes(1) > 1)   i_AitkenSolMass     = &
       get_q_index(standard_q_names%AITKEN_SOL_MASS,     'casim') 
    if (soluble_modes(1) > 0)   i_AitkenSolNumber   = &
       get_q_index(standard_q_names%AITKEN_SOL_NUMBER,   'casim') 
    if (soluble_modes(2) > 1)   i_AccumSolMass      = &
       get_q_index(standard_q_names%ACCUM_SOL_MASS,      'casim') 
    if (soluble_modes(2) > 0)   i_AccumSolNumber    = &
       get_q_index(standard_q_names%ACCUM_SOL_NUMBER,    'casim') 
    if (soluble_modes(3) > 1)   i_CoarseSolMass     = &
       get_q_index(standard_q_names%COARSE_SOL_MASS,     'casim') 
    if (soluble_modes(3) > 0)   i_CoarseSolNumber   = &
       get_q_index(standard_q_names%COARSE_SOL_NUMBER,   'casim') 
    if (active_cloud(isol))     i_ActiveSolLiquid   = &
       get_q_index(standard_q_names%ACTIVE_SOL_LIQUID,   'casim') 
    if (active_rain(isol))      i_ActiveSolRain     = &
       get_q_index(standard_q_names%ACTIVE_SOL_RAIN,     'casim') 
    if (insoluble_modes(2) > 1) i_CoarseDustMass    = &
       get_q_index(standard_q_names%COARSE_DUST_MASS,    'casim') 
    if (insoluble_modes(2) > 0) i_CoarseDustnumber  = &
       get_q_index(standard_q_names%COARSE_DUST_NUMBER,  'casim') 
    if (active_ice(iinsol))     i_ActiveInsolIce    = &
       get_q_index(standard_q_names%ACTIVE_INSOL_ICE,    'casim') 
    if (active_ice(isol))       i_ActiveSolIce      = &
       get_q_index(standard_q_names%ACTIVE_SOL_ICE,      'casim') 
    if (active_cloud(iinsol))   i_ActiveInsolLiquid = &
       get_q_index(standard_q_names%ACTIVE_INSOL_LIQUID, 'casim') 
    if (insoluble_modes(1) > 1) i_AccumInsolMass    = &
       get_q_index(standard_q_names%ACCUM_INSOL_MASS,    'casim') 
    if (insoluble_modes(1) > 0) i_AccumInsolNumber  = &
       get_q_index(standard_q_names%ACCUM_INSOL_NUMBER,  'casim') 
    if (active_number(isol))    i_ActiveSolNumber   = &
       get_q_index(standard_q_names%ACTIVE_SOL_NUMBER,   'casim') 
    if (active_number(iinsol))  i_ActiveInsolNumber = &
         get_q_index(standard_q_names%ACTIVE_INSOL_NUMBER, 'casim')

    ! set logicals for the microphysics diagnostics: process rates
    casdiags % l_dth = .TRUE.
    casdiags % l_dqv = .TRUE.
    if ( nq_l>0 ) casdiags % l_dqc = .TRUE.
    if ( nq_r>0 ) casdiags % l_dqr = .TRUE.
    if ( l_pcond ) casdiags % l_pcond = .TRUE.
    if ( l_psedl ) then
       casdiags % l_psedl = .TRUE.
       casdiags % l_surface_rain = .TRUE.
       casdiags % l_precip = .TRUE.
    endif
    if ( l_praut ) casdiags % l_praut = .TRUE.
    if ( l_pracw ) casdiags % l_pracw = .TRUE.
    if ( l_prevp ) casdiags % l_prevp = .TRUE.
    if ( l_psedr ) then
       casdiags % l_psedr = .TRUE.
       casdiags % l_surface_rain = .TRUE.
       casdiags % l_precip = .TRUE.
    endif
    if (.not. l_warm) then
       if ( nq_i>0 ) casdiags % l_dqi = .TRUE.
       if ( nq_s>0 ) casdiags % l_dqs = .TRUE.
       if ( nq_g>0 ) casdiags % l_dqg = .TRUE.
       if ( l_phomc ) casdiags % l_phomc = .TRUE.
       if ( l_pinuc ) casdiags % l_pinuc = .TRUE.
       if ( l_pidep ) casdiags % l_pidep = .TRUE.
       if ( l_piacw ) casdiags % l_piacw = .TRUE.
       if ( l_pisub ) casdiags % l_pisub = .TRUE.
       if ( l_pimlt ) casdiags % l_pimlt = .TRUE.
       if ( l_psedi ) then
          casdiags % l_psedi = .TRUE.
          casdiags % l_surface_snow = .TRUE.
          casdiags % l_precip = .TRUE.
       endif
       if ( l_psmlt ) casdiags % l_psmlt = .TRUE.
       if ( l_psaut ) casdiags % l_psaut = .TRUE.
       if ( l_psaci ) casdiags % l_psaci = .TRUE.
       if ( l_psacw ) casdiags % l_psacw = .TRUE.
       if ( l_psacr ) casdiags % l_psacr = .TRUE.
       if ( l_pssub ) casdiags % l_pssub = .TRUE.
       if ( l_psdep ) casdiags % l_psdep = .TRUE.
       if ( l_pseds ) then
          casdiags % l_pseds = .TRUE.
          casdiags % l_surface_snow = .TRUE.
          casdiags % l_precip = .TRUE.
       endif
       if ( l_pgacw ) casdiags % l_pgacw = .TRUE.
       if ( l_pgacs ) casdiags % l_pgacs = .TRUE.
       if ( l_pgmlt ) casdiags % l_pgmlt = .TRUE.
       if ( l_pgsub ) casdiags % l_pgsub = .TRUE.
       if ( l_psedg ) then
          casdiags % l_psedg = .TRUE.
          casdiags % l_surface_graup = .TRUE.
          casdiags % l_precip = .TRUE.
       endif
    endif

    ! allocate diagnostic space in casdiags depending on the logicals defined above
    CALL allocate_diagnostic_space(its, ite, jts, jte, kts, kte)
    ! this is no longer needed since can use cas_monc_dgs structure but keep for now
    allocate(surface_precip(y_size_local, x_size_local))
    ! allocate diagnostic space for MONC fields to export to IO server
    call allocate_casim_monc_dgs_space(current_state, casdiags)

  end subroutine initialisation_callback

  !> Called for each column per timestep this will calculate the microphysical tendencies
  !! @param current_state The current model state
  subroutine timestep_callback(current_state)
    type(model_state_type), target, intent(inout) :: current_state

    REAL(wp) :: dtwp
    INTEGER :: icol, jcol, iqx, target_x_index, target_y_index

    icol=current_state%column_local_x
    jcol=current_state%column_local_y
    target_y_index=jcol-current_state%local_grid%halo_size(Y_INDEX)
    target_x_index=icol-current_state%local_grid%halo_size(X_INDEX)

    !if (current_state%first_timestep_column) then
    !   psedl_tot(:)= 0.0_DEFAULT_PRECISION
    !   pcond_tot(:)= 0.0_DEFAULT_PRECISION
    !endif
    
    if (current_state%halo_column .or. current_state%timestep < 2) return

    if (current_state%field_stepping == FORWARD_STEPPING)then
      call log_master_log(LOG_ERROR, 'Currently, CASIM assumes CENTERED_STEPPING')
      dtwp = current_state%dtm
    else
      dtwp = 2.0*current_state%dtm
    end if

    ! Initialize aerosol fields to zero...
    AitkenSolMass = 0.0
    dAitkenSolMass = 0.0
    AitkenSolNumber = 0.0
    dAitkenSolNumber = 0.0
    AccumSolMass = 0.0
    dAccumSolMass = 0.0
    AccumSolNumber = 0.0
    dAccumSolNumber = 0.0
    CoarseSolMass = 0.0
    dCoarseSolMass = 0.0
    CoarseSolNumber = 0.0
    dCoarseSolNumber = 0.0
    ActiveSolLiquid = 0.0
    dActiveSolLiquid = 0.0
    CoarseDustMass = 0.0
    dCoarseDustMass = 0.0
    CoarseDustNumber = 0.0
    dCoarseDustNumber = 0.0
    ActiveInsolIce = 0.0
    dActiveInsolIce = 0.0
    ActiveSolIce = 0.0
    dActiveSolIce = 0.0
    ActiveInsolLiquid = 0.0
    dActiveInsolLiquid = 0.0
    AccumInsolMass = 0.0
    dAccumInsolMass = 0.0
    AccumInsolNumber = 0.0
    dAccumInsolNumber = 0.0
    ActiveSolNumber = 0.0
    dActiveSolNumber = 0.0
    ActiveInsolNumber = 0.0
    dActiveInsolNumber = 0.0
    AitkenSolBk = 0.0
    AccumSolBk = 0.0 
    CoarseSolBk = 0.0

    theta(:,1,1) = current_state%zth%data(:, jcol, icol) + current_state%global_grid%configuration%vertical%thref(:)
    dth(:,1,1) = current_state%sth%data(:, jcol, icol)
    exner(:,1,1) = current_state%global_grid%configuration%vertical%rprefrcp(:)
    pressure(:,1,1) = current_state%global_grid%configuration%vertical%prefn(:)
    z_centre(:,1,1) = current_state%global_grid%configuration%vertical%zn(:)
    dz(:,1,1) = current_state%global_grid%configuration%vertical%dz(:)
    z_half(:kte-1,1,1) = current_state%global_grid%configuration%vertical%z(:)
    rho(:,1,1) = current_state%global_grid%configuration%vertical%rhon(:)
    w(:,1,1) = current_state%zw%data(:, jcol, icol)
    tke(:,1,1) = 0.1 ! Test value

    iqx = iqv
    qv(:,1,1) = current_state%zq(iqx)%data(:,jcol,icol)
    dqv(:,1,1) = current_state%sq(iqx)%data(:,jcol,icol)

    ! Warm microphysical fields
    IF (nq_l > 0)then
      iqx = iql
      qc(:,1,1) = current_state%zq(iqx)%data(:,jcol,icol)
      dqc(:,1,1) = current_state%sq(iqx)%data(:,jcol,icol)
    end IF
    IF (nq_r > 0)then
      iqx = iqr
      qr(:,1,1) = current_state%zq(iqx)%data(:,jcol,icol)
      dqr(:,1,1) = current_state%sq(iqx)%data(:,jcol,icol)
    end IF
    IF (nq_l > 1)then
      iqx = inl
      nc(:,1,1) = current_state%zq(iqx)%data(:,jcol,icol)
      dnc(:,1,1) = current_state%sq(iqx)%data(:,jcol,icol)
    end IF
    IF (nq_r > 1)then
      iqx = inr
      nr(:,1,1) = current_state%zq(iqx)%data(:,jcol,icol)
      dnr(:,1,1) = current_state%sq(iqx)%data(:,jcol,icol)
    end IF
    IF (nq_r > 2)then
      iqx = i3mr
      m3r(:,1,1) = current_state%zq(iqx)%data(:,jcol,icol)
      dm3r(:,1,1) = current_state%sq(iqx)%data(:,jcol,icol)
    end IF

    ! Ice microphysical fields
    IF (nq_i > 0)then
      iqx = iqi
      qi(:,1,1) = current_state%zq(iqx)%data(:,jcol,icol)
      dqi(:,1,1) = current_state%sq(iqx)%data(:,jcol,icol)
    end IF
    IF (nq_s > 0)then
      iqx = iqs
      qs(:,1,1) = current_state%zq(iqx)%data(:,jcol,icol)
      dqs(:,1,1) = current_state%sq(iqx)%data(:,jcol,icol)
    end IF
    IF (nq_g > 0)then
      iqx = iqg
      qg(:,1,1) = current_state%zq(iqx)%data(:,jcol,icol)
      dqg(:,1,1) = current_state%sq(iqx)%data(:,jcol,icol)
    end IF
    IF (nq_i > 1)then
      iqx = ini
      ni(:,1,1) = current_state%zq(iqx)%data(:,jcol,icol)
      dni(:,1,1) = current_state%sq(iqx)%data(:,jcol,icol)
    end IF
    IF (nq_s > 1)then
      iqx = ins
      ns(:,1,1) = current_state%zq(iqx)%data(:,jcol,icol)
      dns(:,1,1) = current_state%sq(iqx)%data(:,jcol,icol)
    end IF
    IF (nq_g > 1)then
      iqx = ing
      ng(:,1,1) = current_state%zq(iqx)%data(:,jcol,icol)
      dng(:,1,1) = current_state%sq(iqx)%data(:,jcol,icol)
    end IF
    IF (nq_s > 2)then
      iqx = i3ms
      m3s(:,1,1) = current_state%zq(iqx)%data(:,jcol,icol)
      dm3s(:,1,1) = current_state%sq(iqx)%data(:,jcol,icol)
    end IF
    IF (nq_g > 2)then
      iqx = i3mg
      m3g(:,1,1) = current_state%zq(iqx)%data(:,jcol,icol)
      dm3g(:,1,1) = current_state%sq(iqx)%data(:,jcol,icol)
    end IF

    ! Aerosol fields

    if (i_AitkenSolMass>0)     AitkenSolMass(:,1,1)      = current_state%zq(i_AitkenSolMass)%data(:,jcol,icol)
    if (i_AitkenSolMass>0)    dAitkenSolMass(:,1,1)     = current_state%sq(i_AitkenSolMass)%data(:,jcol,icol)
    if (i_AitkenSolNumber>0)   AitkenSolNumber(:,1,1)    = current_state%zq(i_AitkenSolNumber)%data(:,jcol,icol)
    if (i_AitkenSolNumber>0)  dAitkenSolNumber(:,1,1)   = current_state%sq(i_AitkenSolNumber)%data(:,jcol,icol)
    if (i_AccumSolMass>0)      AccumSolMass(:,1,1)       = current_state%zq(i_AccumSolMass)%data(:,jcol,icol)
    if (i_AccumSolMass>0)     dAccumSolMass(:,1,1)      = current_state%sq(i_AccumSolMass)%data(:,jcol,icol)
    if (i_AccumSolNumber>0)    AccumSolNumber(:,1,1)     = current_state%zq(i_AccumSolNumber)%data(:,jcol,icol)
    if (i_AccumSolNumber>0)   dAccumSolNumber(:,1,1)    = current_state%sq(i_AccumSolNumber)%data(:,jcol,icol)
    if (i_CoarseSolMass>0)     CoarseSolMass(:,1,1)      = current_state%zq(i_CoarseSolMass)%data(:,jcol,icol)
    if (i_CoarseSolMass>0)    dCoarseSolMass(:,1,1)     = current_state%sq(i_CoarseSolMass)%data(:,jcol,icol)
    if (i_CoarseSolNumber>0)   CoarseSolNumber(:,1,1)    = current_state%zq(i_CoarseSolNumber)%data(:,jcol,icol)
    if (i_CoarseSolNumber>0)  dCoarseSolNumber(:,1,1)   = current_state%sq(i_CoarseSolNumber)%data(:,jcol,icol)
    if (i_ActiveSolLiquid>0)   ActiveSolLiquid(:,1,1)    = current_state%zq(i_ActiveSolLiquid)%data(:,jcol,icol)
    if (i_ActiveSolLiquid>0)  dActiveSolLiquid(:,1,1)   = current_state%sq(i_ActiveSolLiquid)%data(:,jcol,icol)
    if (i_CoarseDustMass>0)    CoarseDustMass(:,1,1)     = current_state%zq(i_CoarseDustMass)%data(:,jcol,icol)
    if (i_CoarseDustMass>0)   dCoarseDustMass(:,1,1)    = current_state%sq(i_CoarseDustMass)%data(:,jcol,icol)
    if (i_CoarseDustNumber>0)  CoarseDustNumber(:,1,1)   = current_state%zq(i_CoarseDustNumber)%data(:,jcol,icol)
    if (i_CoarseDustNumber>0) dCoarseDustNumber(:,1,1)  = current_state%sq(i_CoarseDustNumber)%data(:,jcol,icol)
    if (i_ActiveInsolIce>0)    ActiveInsolIce(:,1,1)     = current_state%zq(i_ActiveInsolIce)%data(:,jcol,icol)
    if (i_ActiveInsolIce>0)   dActiveInsolIce(:,1,1)    = current_state%sq(i_ActiveInsolIce)%data(:,jcol,icol)
    if (i_ActiveSolIce>0)      ActiveSolIce(:,1,1)       = current_state%zq(i_ActiveSolIce)%data(:,jcol,icol)
    if (i_ActiveSolIce>0)     dActiveSolIce(:,1,1)      = current_state%sq(i_ActiveSolIce)%data(:,jcol,icol)
    if (i_ActiveInsolLiquid>0) ActiveInsolLiquid(:,1,1)  = current_state%zq(i_ActiveInsolLiquid)%data(:,jcol,icol)
    if (i_ActiveInsolLiquid>0)dActiveInsolLiquid(:,1,1) = current_state%sq(i_ActiveInsolLiquid)%data(:,jcol,icol)
    if (i_AccumInsolMass>0)    AccumInsolMass(:,1,1)     = current_state%zq(i_AccumInsolMass)%data(:,jcol,icol)
    if (i_AccumInsolMass>0)   dAccumInsolMass(:,1,1)    = current_state%sq(i_AccumInsolMass)%data(:,jcol,icol)
    if (i_AccumInsolNumber>0)  AccumInsolNumber(:,1,1)   = current_state%zq(i_AccumInsolNumber)%data(:,jcol,icol)
    if (i_AccumInsolNumber>0) dAccumInsolNumber(:,1,1)  = current_state%sq(i_AccumInsolNumber)%data(:,jcol,icol)
    if (i_ActiveSolNumber>0)   ActiveSolNumber(:,1,1)    = current_state%zq(i_ActiveSolNumber)%data(:,jcol,icol)
    if (i_ActiveSolNumber>0)  dActiveSolNumber(:,1,1)   = current_state%sq(i_ActiveSolNumber)%data(:,jcol,icol)
    if (i_ActiveInsolNumber>0) ActiveInsolNumber(:,1,1)  = current_state%zq(i_ActiveInsolNumber)%data(:,jcol,icol)
    if (i_ActiveInsolNumber>0)dActiveInsolNumber(:,1,1) = current_state%sq(i_ActiveInsolNumber)%data(:,jcol,icol)

    CALL shipway_microphysics(                     &
                                ! in
       its, ite,                                   &
       jts, jte,                                   &
       kts, kte,                                   &
       dtwp,                                       &
       qv, qc, qr,                                 &
       nc, nr, m3r,                                &
       qi, qs, qg,                                 &
       ni, ns, ng,                                 &
       m3s, m3g,                                   &
       theta,                                      &
       AitkenSolMass, AitkenSolNumber,             &
       AccumSolMass, AccumSolNumber,               &
       CoarseSolMass, CoarseSolNumber,             &
       ActiveSolLiquid,                            &
       ActiveSolRain,                              &
       CoarseDustMass, CoarseDustNumber,           &
       ActiveInsolIce,                             &
       ActiveSolIce,                               &
       ActiveInsolLiquid,                          &
       AccumInsolMass,                             &
       AccumInsolNumber,                           &
       ActiveSolNumber,                            &
       ActiveInsolNumber,                          &
       AitkenSolBk,                                &
       AccumSolBk,                                 & 
       CoarseSolBk,                                &
       exner,                                      &
       pressure, rho,                              &
       w, tke,                                     &
       z_half, z_centre,                           &
       dz,                                         &
                                ! in/out
       dqv, dqc, dqr, dnc, dnr, dm3r,              &
       dqi, dqs, dqg, dni, dns, dng, dm3s, dm3g,   &
       dth,                                        &
       dAitkenSolMass, dAitkenSolNumber,           &
       dAccumSolMass, dAccumSolNumber,             &
       dCoarseSolMass, dCoarseSolNumber,           &
       dActiveSolLiquid,                           &
       dActiveSolRain,                             &
       dCoarseDustMass, dCoarseDustNumber,         &
       dActiveInsolIce,                            &
       dActiveSolIce,                              &
       dActiveInsolLiquid,                         &
       dAccumInsolMass,                            &
       dAccumInsolNumber,                          &
       dActiveSolNumber,                           &
       dActiveInsolNumber,                         &
       ils, ile,                                   &
       jls, jle,                                   &
       kls, kle,                                   &
       l_tendency=.TRUE.                           &
       )

    ! write back the tendencies
    current_state%sth%data(:,jcol,icol) = current_state%sth%data(:,jcol,icol) + dth(:,1,1)

    iqx = iqv
    current_state%sq(iqx)%data(:,jcol,icol) = current_state%sq(iqx)%data(:,jcol,icol) + dqv(:,1,1)

    ! Warm microphysical fields
    IF (nq_l > 0)then
      iqx = iql
      current_state%sq(iqx)%data(:,jcol,icol) = current_state%sq(iqx)%data(:,jcol,icol) + dqc(:,1,1)
    end IF

    IF (nq_r > 0)then
      iqx = iqr
      current_state%sq(iqx)%data(:,jcol,icol) = current_state%sq(iqx)%data(:,jcol,icol) + dqr(:,1,1)
    end IF
    IF (nq_l > 1)then
      iqx = inl
      current_state%sq(iqx)%data(:,jcol,icol) = current_state%sq(iqx)%data(:,jcol,icol) + dnc(:,1,1)
    end IF
    IF (nq_r > 1)then
      iqx = inr
      current_state%sq(iqx)%data(:,jcol,icol) = current_state%sq(iqx)%data(:,jcol,icol) + dnr(:,1,1)
    end IF
    IF (nq_r > 2)then
      iqx = i3mr
      current_state%sq(iqx)%data(:,jcol,icol) = current_state%sq(iqx)%data(:,jcol,icol) + dm3r(:,1,1)
    end IF

    ! Ice microphysical fields
    IF (nq_i > 0)then
      iqx = iqi
      current_state%sq(iqx)%data(:,jcol,icol) = current_state%sq(iqx)%data(:,jcol,icol) + dqi(:,1,1)
    end IF
    IF (nq_s > 0)then
      iqx = iqs
      current_state%sq(iqx)%data(:,jcol,icol) = current_state%sq(iqx)%data(:,jcol,icol) + dqs(:,1,1)
    end IF
    IF (nq_g > 0)then
      iqx = iqg
      current_state%sq(iqx)%data(:,jcol,icol) = current_state%sq(iqx)%data(:,jcol,icol) + dqg(:,1,1)
    end IF
    IF (nq_i > 1)then
      iqx = ini
      current_state%sq(iqx)%data(:,jcol,icol) = current_state%sq(iqx)%data(:,jcol,icol) + dni(:,1,1)
    end IF
    IF (nq_s > 1)then
      iqx = ins
      current_state%sq(iqx)%data(:,jcol,icol) = current_state%sq(iqx)%data(:,jcol,icol) + dns(:,1,1)
    end IF
    IF (nq_g > 1)then
      iqx = ing
      current_state%sq(iqx)%data(:,jcol,icol) = current_state%sq(iqx)%data(:,jcol,icol) + dng(:,1,1)
    end IF
    IF (nq_s > 2)then
      iqx = i3ms
      current_state%sq(iqx)%data(:,jcol,icol) = current_state%sq(iqx)%data(:,jcol,icol) + dm3s(:,1,1)
    end IF
    IF (nq_g > 2)then
      iqx = i3mg
      current_state%sq(iqx)%data(:,jcol,icol) = current_state%sq(iqx)%data(:,jcol,icol) + dm3g(:,1,1)
    end IF

    ! Aerosol fields

    if (i_AitkenSolMass>0)     current_state%sq(i_AitkenSolMass)%data(:,jcol,icol) &
       = current_state%sq(i_AitkenSolMass)%data(:,jcol,icol) + dAitkenSolMass(:,1,1)
    if (i_AitkenSolNumber>0)   current_state%sq(i_AitkenSolNumber)%data(:,jcol,icol) &
       =  current_state%sq(i_AitkenSolNumber)%data(:,jcol,icol) + dAitkenSolNumber(:,1,1)
    if (i_AccumSolMass>0)      current_state%sq(i_AccumSolMass)%data(:,jcol,icol) &
       =  current_state%sq(i_AccumSolMass)%data(:,jcol,icol) + dAccumSolMass(:,1,1)
    if (i_AccumSolNumber>0)    current_state%sq(i_AccumSolNumber)%data(:,jcol,icol) &
       =  current_state%sq(i_AccumSolNumber)%data(:,jcol,icol) + dAccumSolNumber(:,1,1)
    if (i_CoarseSolMass>0)     current_state%sq(i_CoarseSolMass)%data(:,jcol,icol) &
       =  current_state%sq(i_CoarseSolMass)%data(:,jcol,icol) + dCoarseSolMass(:,1,1)
    if (i_CoarseSolNumber>0)   current_state%sq(i_CoarseSolNumber)%data(:,jcol,icol) &
       =  current_state%sq(i_CoarseSolNumber)%data(:,jcol,icol) + dCoarseSolNumber(:,1,1)
    if (i_ActiveSolLiquid>0)   current_state%sq(i_ActiveSolLiquid)%data(:,jcol,icol) &
       =  current_state%sq(i_ActiveSolLiquid)%data(:,jcol,icol) + dActiveSolLiquid(:,1,1)
    if (i_CoarseDustMass>0)    current_state%sq(i_CoarseDustMass)%data(:,jcol,icol) &
       =  current_state%sq(i_CoarseDustMass)%data(:,jcol,icol) + dCoarseDustMass(:,1,1)
    if (i_CoarseDustNumber>0)  current_state%sq(i_CoarseDustNumber)%data(:,jcol,icol) &
       =  current_state%sq(i_CoarseDustNumber)%data(:,jcol,icol) + dCoarseDustNumber(:,1,1)
    if (i_ActiveInsolIce>0)    current_state%sq(i_ActiveInsolIce)%data(:,jcol,icol) &
       =  current_state%sq(i_ActiveInsolIce)%data(:,jcol,icol) + dActiveInsolIce(:,1,1)
    if (i_ActiveSolIce>0)      current_state%sq(i_ActiveSolIce)%data(:,jcol,icol) &
       =  current_state%sq(i_ActiveSolIce)%data(:,jcol,icol) + dActiveSolIce(:,1,1)
    if (i_ActiveInsolLiquid>0) current_state%sq(i_ActiveInsolLiquid)%data(:,jcol,icol) &
       =  current_state%sq(i_ActiveInsolLiquid)%data(:,jcol,icol) + dActiveInsolLiquid(:,1,1)
    if (i_AccumInsolMass>0)    current_state%sq(i_AccumInsolMass)%data(:,jcol,icol) &
       =  current_state%sq(i_AccumInsolMass)%data(:,jcol,icol) + dAccumInsolMass(:,1,1)
    if (i_AccumInsolNumber>0)  current_state%sq(i_AccumInsolNumber)%data(:,jcol,icol) &
       =  current_state%sq(i_AccumInsolNumber)%data(:,jcol,icol) + dAccumInsolNumber(:,1,1)
    if (i_ActiveSolNumber>0)   current_state%sq(i_ActiveSolNumber)%data(:,jcol,icol) &
       =  current_state%sq(i_ActiveSolNumber)%data(:,jcol,icol) + dActiveSolNumber(:,1,1)
    if (i_ActiveInsolNumber>0) current_state%sq(i_ActiveInsolNumber)%data(:,jcol,icol) &
       =  current_state%sq(i_ActiveInsolNumber)%data(:,jcol,icol) + dActiveInsolNumber(:,1,1)

    ! for total surface precipitation, sum the surface rain rate (cloud + rain which is precip_r)
    ! and surface
    ! snow rate (precip_s), which is the sum of ice, snow and graupel (See micromain.F90 in casim for
    ! calculation).
    if (l_warm .or. .not. casdiags % l_surface_snow ) then
       surface_precip(target_y_index,target_x_index) = &
            casdiags % SurfaceRainR(1,1)
    else
       surface_precip(target_y_index,target_x_index) = &
            casdiags % SurfaceRainR(1,1) + casdiags % SurfaceSnowR(1,1)
    endif
    call populate_casim_monc_dg(current_state, casdiags)

    
  end subroutine timestep_callback


  !! Reads the casim configuration 
  !! @param current_state The current model state
  subroutine read_configuration(current_state)


    Use mphys_parameters, only: p1, p2, p3, sp1, sp2, sp3   

    type(model_state_type), target, intent(inout) :: current_state

    integer :: ierr

    option          = options_get_integer(current_state%options_database, 'option')
    diag_mu_option  = options_get_integer(current_state%options_database, 'diag_mu_option')
    iopt_act        = options_get_integer(current_state%options_database, 'iopt_act')
    iopt_inuc       = options_get_integer(current_state%options_database, 'iopt_inuc')
    process_level   = options_get_integer(current_state%options_database, 'process_level')
    aerosol_option  = options_get_integer(current_state%options_database, 'aerosol_option')
    max_step_length = options_get_real(current_state%options_database, 'max_step_length')
    max_sed_length  = options_get_real(current_state%options_database, 'max_sed_length')
    p1              = options_get_real(current_state%options_database, 'p1')
    p2              = options_get_real(current_state%options_database, 'p2')
    p3              = options_get_real(current_state%options_database, 'p3')
    sp1             = options_get_real(current_state%options_database, 'sp1')
    sp2             = options_get_real(current_state%options_database, 'sp2')
    sp3             = options_get_real(current_state%options_database, 'sp3')
    max_mu          = options_get_real(current_state%options_database, 'max_mu')
    fix_mu          = options_get_real(current_state%options_database, 'fix_mu')

    l_aaut          = options_get_logical(current_state%options_database, 'l_aaut')
    l_aacc          = options_get_logical(current_state%options_database, 'l_aacc')
    l_aevp          = options_get_logical(current_state%options_database, 'l_aevp')
    l_ased          = options_get_logical(current_state%options_database, 'l_ased')
    l_warm          = options_get_logical(current_state%options_database, 'l_warm')
    l_inuc          = options_get_logical(current_state%options_database, 'l_inuc')
    l_iaut          = options_get_logical(current_state%options_database, 'l_iaut')
    l_idep          = options_get_logical(current_state%options_database, 'l_idep')
    l_iacw          = options_get_logical(current_state%options_database, 'l_iacw')
    l_active_inarg2000 = options_get_logical(current_state%options_database, 'l_active_inarg2000')
    l_separate_rain = options_get_logical(current_state%options_database, 'l_separate_rain')
    l_sg            = options_get_logical(current_state%options_database, 'l_sg')
    l_g             = options_get_logical(current_state%options_database, 'l_g')
    l_passive       = options_get_logical(current_state%options_database, 'l_passive')
    l_passive3m     = options_get_logical(current_state%options_database, 'l_passive3m')
    l_limit_psd     = options_get_logical(current_state%options_database, 'l_limit_psd')
    l_override_checks = options_get_logical(current_state%options_database, 'l_override_checks')
    l_raci_g        = options_get_logical(current_state%options_database, 'l_raci_g')
    l_onlycollect   = options_get_logical(current_state%options_database, 'l_onlycollect')
    l_abelshipway   = options_get_logical(current_state%options_database, 'l_abelshipway')
    l_cons          = options_get_logical(current_state%options_database, 'l_cons')
    l_rain          = options_get_logical(current_state%options_database, 'l_rain')
    l_sed_3mdiff    = options_get_logical(current_state%options_database, 'l_sed_3mdiff')
    l_sed_icecloud_as_1m = options_get_logical(current_state%options_database, 'l_sed_icecloud_as_1m')
    l_tidy_conserve_E = options_get_logical(current_state%options_database, 'l_tidy_conserve_E')
    l_tidy_conserve_q = options_get_logical(current_state%options_database, 'l_tidy_conserve_q')

    l_inhom_revp    = options_get_logical(current_state%options_database, 'l_inhom_revp')
    l_pcond         = options_get_logical(current_state%options_database, 'l_pcond')
    l_praut         = options_get_logical(current_state%options_database, 'l_praut')
    l_pracw         = options_get_logical(current_state%options_database, 'l_pracw')
    l_pracr         = options_get_logical(current_state%options_database, 'l_pracr')
    l_prevp         = options_get_logical(current_state%options_database, 'l_prevp')
    l_psedl         = options_get_logical(current_state%options_database, 'l_psedl')
    l_psedr         = options_get_logical(current_state%options_database, 'l_psedr')
    l_ptidy         = options_get_logical(current_state%options_database, 'l_ptidy')
    l_ptidy2        = options_get_logical(current_state%options_database, 'l_ptidy2')
    l_pinuc         = options_get_logical(current_state%options_database, 'l_pinuc')
    l_pidep         = options_get_logical(current_state%options_database, 'l_pidep')
    l_piacw         = options_get_logical(current_state%options_database, 'l_piacw')
    l_psaut         = options_get_logical(current_state%options_database, 'l_psaut')
    l_psdep         = options_get_logical(current_state%options_database, 'l_psdep')
    l_psacw         = options_get_logical(current_state%options_database, 'l_psacw')
    l_pgdep         = options_get_logical(current_state%options_database, 'l_pgdep')
    l_pseds         = options_get_logical(current_state%options_database, 'l_pseds')
    l_psedi         = options_get_logical(current_state%options_database, 'l_psedi')
    l_psedg         = options_get_logical(current_state%options_database, 'l_psedg')
    l_psaci         = options_get_logical(current_state%options_database, 'l_psaci')
    l_praci         = options_get_logical(current_state%options_database, 'l_praci')
    l_psacr         = options_get_logical(current_state%options_database, 'l_psacr')
    l_pgacr         = options_get_logical(current_state%options_database, 'l_pgacr')
    l_pgacw         = options_get_logical(current_state%options_database, 'l_pgacw')
    l_pgaci         = options_get_logical(current_state%options_database, 'l_pgaci')
    l_pgacs         = options_get_logical(current_state%options_database, 'l_pgacs')
    l_piagg         = options_get_logical(current_state%options_database, 'l_piagg')
    l_psagg         = options_get_logical(current_state%options_database, 'l_psagg')
    l_pgagg         = options_get_logical(current_state%options_database, 'l_pgagg')
    l_psbrk         = options_get_logical(current_state%options_database, 'l_psbrk')
    l_pgshd         = options_get_logical(current_state%options_database, 'l_pgshd')
    l_pihal         = options_get_logical(current_state%options_database, 'l_pihal')
    l_psmlt         = options_get_logical(current_state%options_database, 'l_psmlt')
    l_pgmlt         = options_get_logical(current_state%options_database, 'l_pgmlt')
    l_phomr         = options_get_logical(current_state%options_database, 'l_phomr')
    l_phomc         = options_get_logical(current_state%options_database, 'l_phomc')
    l_pssub         = options_get_logical(current_state%options_database, 'l_pssub')
    l_pgsub         = options_get_logical(current_state%options_database, 'l_pgsub')
    l_pisub         = options_get_logical(current_state%options_database, 'l_pisub')
    l_pimlt         = options_get_logical(current_state%options_database, 'l_pimlt')
    l_gamma_online  = options_get_logical(current_state%options_database, 'l_gamma_online')
    l_subseds_maxv  = options_get_logical(current_state%options_database, 'l_subseds_maxv')
    l_sed_eulexp  = options_get_logical(current_state%options_database, 'l_sed_eulexp')
    cfl_vt_max      = options_get_real(current_state%options_database, 'cfl_vt_max')
    l_kfsm          = options_get_logical(current_state%options_database, 'l_kfsm')

  end subroutine read_configuration

  subroutine field_information_retrieval_callback(current_state, name, field_information)
    type(model_state_type), target, intent(inout) :: current_state
    character(len=*), intent(in) :: name
    type(component_field_information_type), intent(out) :: field_information

    field_information%field_type=COMPONENT_ARRAY_FIELD_TYPE
    field_information%data_type=COMPONENT_DOUBLE_DATA_TYPE
    if (name .eq. "surface_precip") then
       field_information%number_dimensions=2
       field_information%dimension_sizes(1)=current_state%local_grid%size(Y_INDEX)
       field_information%dimension_sizes(2)=current_state%local_grid%size(X_INDEX)
    !else if (name .eq. "pcond_total" .or. name .eq. "psedl_total") then
    !   field_information%number_dimensions=1
    !   field_information%dimension_sizes(1)=current_state%local_grid%size(Z_INDEX)
    else
       field_information%number_dimensions=3
       field_information%dimension_sizes(1)=current_state%local_grid%size(Z_INDEX)
       field_information%dimension_sizes(2)=current_state%local_grid%size(Y_INDEX)
       field_information%dimension_sizes(3)=current_state%local_grid%size(X_INDEX)
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
    
    integer :: i

    if (name .eq. "surface_precip") then
       allocate(field_value%real_2d_array(current_state%local_grid%size(Y_INDEX), &
            current_state%local_grid%size(X_INDEX)))
       field_value%real_2d_array(:,:)= surface_precip(:,:)
    else if (name .eq. "condensation_rate") then
       allocate(field_value%real_3d_array(current_state%local_grid%size(Z_INDEX),  &
            current_state%local_grid%size(Y_INDEX),                                &
            current_state%local_grid%size(X_INDEX)))
       field_value%real_3d_array(:,:,:) = casim_monc_dgs % pcond(:,:,:)
!!$    else if (name .eq. "pcond_total") then
!!$       allocate(field_value%real_1d_array(current_state%local_grid%size(Z_INDEX)))
!!$       field_value%real_1d_array(:)=pcond_tot(:)
    end if
    
  end subroutine field_value_retrieval_callback

end module casim_mod
