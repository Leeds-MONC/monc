module def_socrates_options

  use datadefn_mod, only : DEFAULT_PRECISION
  
  implicit none

  type str_socrates_options

     ! Indices for cloud
     integer :: iqv=0, iql=0, inl=0, iqr=0,  iqi=0, iqs=0,  iqg=0

     integer :: cloud_representation ! This is read from configuration and determines
                                  ! whether socrates should include cloud
                                  ! 2 = clear sky, 5 = cloudy sky rad calc
      ! number of moments for each hydrometeor - read from configuration file
     integer :: mphys_nq_l  ! cloud liquid
     integer :: mphys_nd_l  ! cloud drop number
     integer :: mphys_nq_r  ! rain
     integer :: mphys_nq_i  ! cloud ice
     integer :: mphys_nq_s  ! snow
     integer :: mphys_nq_g  ! graupel
     
     ! Time and location variables for socrates - read from configuration
     logical :: l_solar_fixed   & ! true equals fixed insolation using value in
          ! solar_fixed
          , l_360                 ! 360 days in year as opposed to 365 (a UM thing 
     ! in the LEM, is this still required??)
     
     real(kind=DEFAULT_PRECISION) :: &
          default_solar_constant,    &
          solar_fixed,               &
          sec_fixed,                 &
          latitude,                  &
          longitude,                 &
          rad_start_time,            &
          rad_time_hours,            &
          rad_year,                  &
          rad_start_day,             &
          rad_day,                   &
          rad_int_time
     !
     
     ! Surface albedo variables for socrates -read from configuration
     logical :: l_variable_srf_albedo
     real(kind=DEFAULT_PRECISION) :: surface_albedo
     !
     ! density of water set in config to 997.0 kgm-3
     real(kind=DEFAULT_PRECISION) :: rho_water
     
     ! Fixed effective radius settings for socrates
     ! NOTE: default in global_config is false, so effective radius from
     !       depends on microphysics output 
     logical :: l_fix_re
     ! Use cloud drop number to work out the effective radius
     logical :: l_use_ndrop
     ! Switch for the liu drop spectral param in effective radius
     logical :: l_use_liu_spec

     real(kind=DEFAULT_PRECISION) :: &
          fixed_cloud_re, fixed_ice_re, fixed_cloud_number, kparam 
     
     ! gas mass mixing ratios
     real(kind=DEFAULT_PRECISION) :: & 
          co2_mmr,                   & 
          n2o_mmr,                   &
          ch4_mmr,                   &
          o2_mmr,                    &
          cfc11_mmr,                 &
          cfc12_mmr,                 &
          cfc113_mmr,                &
          cfc114_mmr,                & 
          hcfc22_mmr,                &
          hfc125_mmr,                &
          hfc134a_mmr
     ! logical to decide if it is a radiation timestep. This is true
     ! if (time > ((rad_last_time + rad_int_time)) > 0.0)). Derived
     ! in the timestep_callback in socrates_couple
     logical :: l_rad_calc
      
   end type str_socrates_options
   
 end module def_socrates_options
 
 
