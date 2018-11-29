module get_and_test_socrates_options_mod

 use datadefn_mod, only : DEFAULT_PRECISION
 use state_mod, only : model_state_type
 use registry_mod, only : is_component_enabled
 use logging_mod, only : LOG_INFO, LOG_WARN, LOG_ERROR, &
       LOG_DEBUG, log_master_log, log_log, log_get_logging_level, &
       log_master_log
 use conversions_mod, only : conv_to_string
 use q_indices_mod, only: get_q_index, standard_q_names
 use optionsdatabase_mod, only : options_get_string, options_get_integer, &
       options_get_real, options_get_logical
 
 use def_socrates_options, only: str_socrates_options
 

contains
  
  subroutine set_and_test_socrates_monc_options(current_state, socrates_opt)

    type(model_state_type), target, intent(inout) :: current_state

    type (str_socrates_options), intent(inout) :: socrates_opt
    
    ! First check the MONC switches are sensible to run socrates
    !
    ! check that lwrad_exp is not on, as this will double the cloud top
    ! longwave cooling
    if (is_component_enabled(current_state%options_database, "lwrad_exponential")) then
       call log_master_log &
            (LOG_ERROR, "Socrates and lwrad_exponential both enabled, switch off on in config - STOP")
    endif
    
    ! check whether potential temperature is active, if not stop run
    if (.not. current_state%th%active) then
       call log_master_log &
            (LOG_ERROR, "Socrates called with theta inactive, this is not supported - STOP")
    endif
    ! determine if the radiation is expecting cloud
    socrates_opt%cloud_representation = &
         options_get_integer(current_state%options_database, "i_cloud_representation")
    if (socrates_opt%cloud_representation == 5) then 
       ! only clear sky radiation calculation, initialise vapour index
       if (current_state%number_q_fields < 1) then 
          call log_master_log(LOG_ERROR, "Socrates called for clear sky calc but no vapour field - STOP")
       else
         socrates_opt%iqv=get_q_index(standard_q_names%VAPOUR, 'socrates_couple')
       endif
    else if (socrates_opt%cloud_representation == 2 .or. socrates_opt%cloud_representation == 1) then
       ! Do some tests to make sure config is sensible
       ! (these tests are a bit overkill, better safe than sorry)
       if (.not. is_component_enabled(current_state%options_database, "casim") &
            .and. .not. is_component_enabled(current_state%options_database, "simplecloud")) then
          call log_master_log &
               (LOG_ERROR, "Socrates called for cloudy sky but no microphysics scheme enabled - STOP")
       endif
       if (current_state%passive_q ) then
          call log_master_log &
               (LOG_ERROR, "Socrates called for cloudy sky but q is passive so not cloud or vapour - STOP")
       endif
       if (current_state%number_q_fields < 2) then 
          call log_master_log &
               (LOG_ERROR, "Socrates called for clear and cloud sky calc but no vapour or cloud field - STOP")
       endif
       ! Now monc is happy the set-up seems OK, set-up the microphysics logicals and allocate
       ! arrays
       !
       ! set vapour by default
       socrates_opt%iqv=get_q_index(standard_q_names%VAPOUR, 'socrates_couple')
       ! read from options database to see which hydrometeors are available
       ! for radiation calculation. NOTE: This can differ from the number of hydrometeors
       ! active in the microphysics scheme 
       socrates_opt%mphys_nq_l=options_get_integer(current_state%options_database, "mphys_nq_l")
       socrates_opt%mphys_nd_l=options_get_integer(current_state%options_database, "mphys_nd_l")
       socrates_opt%mphys_nq_r=options_get_integer(current_state%options_database, "mphys_nq_r")
       socrates_opt%mphys_nq_i=options_get_integer(current_state%options_database, "mphys_nq_i")
       socrates_opt%mphys_nq_s=options_get_integer(current_state%options_database, "mphys_nq_s")
       socrates_opt%mphys_nq_g=options_get_integer(current_state%options_database, "mphys_nq_g")

       if (socrates_opt%mphys_nq_l > 0) &
            socrates_opt%iql=get_q_index(standard_q_names%CLOUD_LIQUID_MASS, 'socrates_couple')
       if (socrates_opt%mphys_nq_r > 0) &
            socrates_opt%iqr=get_q_index(standard_q_names%RAIN_MASS, 'socrates_couple')
       if (socrates_opt%mphys_nq_i > 0) &
            socrates_opt%iqi=get_q_index(standard_q_names%ICE_MASS, 'socrates_couple')
       if (socrates_opt%mphys_nq_s > 0) &
            socrates_opt%iqs=get_q_index(standard_q_names%SNOW_MASS, 'socrates_couple')
       if (socrates_opt%mphys_nq_g > 0) &
            socrates_opt%iqg=get_q_index(standard_q_names%GRAUPEL_MASS, 'socrates_couple')
       ! test for fixed effective radius in config
       socrates_opt%l_fix_re = options_get_logical(current_state%options_database, "l_fix_re")
       socrates_opt%l_use_ndrop = options_get_logical(current_state%options_database, "l_use_ndrop")
       if ((socrates_opt%l_fix_re .and. socrates_opt%l_use_ndrop) .or. &
            (.not. socrates_opt%l_fix_re .and. .not. socrates_opt%l_use_ndrop)) then
          call log_master_log &
                     (LOG_ERROR, "Socrates - l_fix_re and l_use_ndrop both true or both false, please pick one - STOP")
       endif
       if (socrates_opt%l_fix_re) then
          socrates_opt%fixed_cloud_re = options_get_real(current_state%options_database, "fixed_cloud_re")
          socrates_opt%fixed_ice_re = options_get_real(current_state%options_database, "fixed_ice_re")
          call log_master_log &
               (LOG_INFO, "Socrates - using fixed cloud effective radius"//trim(conv_to_string(socrates_opt%fixed_cloud_re)))
          call log_master_log &
               (LOG_INFO, "Socrates - using fixed ice effective radius"//trim(conv_to_string(socrates_opt%fixed_ice_re)))
       endif
       if (socrates_opt%l_use_ndrop) then
          socrates_opt%rho_water = options_get_real(current_state%options_database, "rho_water")
          socrates_opt%kparam = options_get_real(current_state%options_database, "kparam")
          socrates_opt%l_use_liu_spec = options_get_logical(current_state%options_database, "l_use_liu_spec")
          if (socrates_opt%mphys_nd_l > 0) then
             if (is_component_enabled(current_state%options_database, "casim")) then
                socrates_opt%inl = get_q_index(standard_q_names%CLOUD_LIQUID_NUMBER, 'socrates_couple')
                socrates_opt%fixed_ice_re = options_get_real(current_state%options_database, "fixed_ice_re")
                call log_master_log &
                     (LOG_INFO, "Socrates using prognostic cloud number from CASIM, fixed ice re="&
                     //trim(conv_to_string(socrates_opt%fixed_ice_re))//" microns")
             else
                call log_master_log &
                     (LOG_ERROR, "Socrates - casim not enabled so no prognostic nd - STOP")
             endif
          else
             socrates_opt%inl = 0
             socrates_opt%fixed_cloud_number = options_get_real(current_state%options_database, "fixed_cloud_number")
             socrates_opt%fixed_ice_re = options_get_real(current_state%options_database, "fixed_ice_re")
             call log_master_log &
                  (LOG_INFO, "Socrates using prescribed fix_cloud_number="&
                  //trim(conv_to_string(socrates_opt%fixed_cloud_number))//" /cm3")
          endif
       endif
    else
       call log_master_log &
            (LOG_ERROR, "Socrates config using unrecognised i_cloud_representation, check config - STOP")    
    endif
    
    ! Get options for time and location variables. These are set in the MONC
    ! configuration. By default they are set to -999.0, which will fail. Hence
    ! user must set them appropriately
    socrates_opt%l_360 = options_get_logical(current_state%options_database, "l_360")
    socrates_opt%l_solar_fixed = &
         options_get_logical(current_state%options_database, "l_solar_fixed")
    socrates_opt%l_no_solar = &
         options_get_logical(current_state%options_database, "l_no_solar")
    socrates_opt%l_rad_calc = .false.

    if (socrates_opt%l_solar_fixed) then
       socrates_opt%solar_fixed = options_get_real(current_state%options_database, "solar_fixed")
       socrates_opt%sec_fixed = options_get_real(current_state%options_database, "sec_fixed")
       if (socrates_opt%solar_fixed .lt. -1.0 .or. socrates_opt%sec_fixed .lt. -1.0) then
          call log_master_log &
               (LOG_ERROR, "Socrates - l_solar_fixed but solar_fixed and/or sec_fixed not set, check config - STOP")
       endif
       socrates_opt%l_variable_srf_albedo = options_get_logical(current_state%options_database, "l_variable_srf_albedo")
       if (socrates_opt%l_variable_srf_albedo) then
          call log_master_log &
               (LOG_INFO, "Socrates - solar_fixed but variable srf albedo, default to fixed srf albedo")
       endif
       socrates_opt%surface_albedo = options_get_real(current_state%options_database, "surface_albedo")
    else ! solar angle will vary so get all the time variables
       ! first get the radiation initial year, day and time variables and add to 
       socrates_opt%rad_year  = options_get_real(current_state%options_database, "rad_start_year")
       if (socrates_opt%rad_year < 0.0) then
          call log_master_log &
               (LOG_ERROR, "Socrates - start year is negative, which is wrong, check config - STOP")
       endif
       !
       socrates_opt%rad_start_day  = options_get_real(current_state%options_database, "rad_start_day")
       if (socrates_opt%rad_start_day < 0.0 .or. socrates_opt%rad_start_day > 360.0) then
          call log_master_log &
               (LOG_ERROR, "Socrates - start day is outside sensible range, check config - STOP")
       endif
       !
       socrates_opt%rad_start_time = options_get_real(current_state%options_database,"rad_start_time") 
       if (socrates_opt%rad_time_hours < 0.0 .or. socrates_opt%rad_time_hours > 24.0) then
          call log_master_log &
               (LOG_ERROR, "Socrates - start time is outside sensible range, check config - STOP")
       endif
       !
       socrates_opt%rad_int_time = options_get_real(current_state%options_database, "rad_int_time")
       if (socrates_opt%rad_int_time <= 0.0) then
          call log_master_log &
               (LOG_WARN, "Socrates - rad_int_time less than 0.0, SOCRATES will be called every timestep")
       endif
       
       ! Now get the surface albedo variables
       socrates_opt%l_variable_srf_albedo = options_get_logical(current_state%options_database, "l_variable_srf_albedo")
       if (socrates_opt%l_variable_srf_albedo) then
          call log_master_log &
            (LOG_ERROR, "Socrates config using variable surface albedo, but this has not been developed. Set to false - STOP")
       endif
       socrates_opt%surface_albedo = options_get_real(current_state%options_database, "surface_albedo")
       
       ! Now get the longitude and latitude variables
       socrates_opt%latitude  = options_get_real(current_state%options_database, "latitude")
       if (socrates_opt%latitude < -90.0 .or. socrates_opt%latitude > 90.0) then
          call log_master_log &
               (LOG_ERROR, "Socrates - latitude is outside sensible range, check config - STOP")
       endif
       socrates_opt%longitude = options_get_real(current_state%options_database, "longitude")
       if (socrates_opt%longitude < -180.0 .or. socrates_opt%longitude > 180.0) then
          call log_master_log &
               (LOG_ERROR, "Socrates - longitude is outside sensible range, check config - STOP")
       endif
    endif ! end l_solar_fixed

    if (socrates_opt%surface_albedo < 0.0 .or. socrates_opt%surface_albedo > 1.0) then
       call log_master_log &
            (LOG_ERROR, "Socrates - surface albedo outside sensible range, check config - STOP")
    endif
     !
     socrates_opt%default_solar_constant =  &
          options_get_real(current_state%options_database, "default_solar_constant")
     
     ! set the well mixed gases. Values are based on UM GA settings
     ! This should be moved to configuration once reading well mixed 
     ! gases from configuration can be shown to bit compare (see #306)  
     ! 
     socrates_opt%co2_mmr = 5.94100e-4
     socrates_opt%n2o_mmr = 4.925e-7
     socrates_opt%ch4_mmr = 9.994e-7
     socrates_opt%o2_mmr = 0.2314
     socrates_opt%cfc12_mmr = 1.129e-9
     socrates_opt%cfc11_mmr = 2.225e-9
     socrates_opt%cfc113_mmr = 0.0
     socrates_opt%cfc114_mmr = 0.0
     socrates_opt%hcfc22_mmr = 0.0
     socrates_opt%hfc125_mmr = 0.0 
     socrates_opt%hfc134a_mmr = 0.0  
    
  end subroutine set_and_test_socrates_monc_options


  end module get_and_test_socrates_options_mod

