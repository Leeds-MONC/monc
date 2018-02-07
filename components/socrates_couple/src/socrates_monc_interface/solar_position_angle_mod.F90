module solar_position_angle_mod
  
  use state_mod, only : model_state_type
  use datadefn_mod, only : DEFAULT_PRECISION
  use def_socrates_options, only: str_socrates_options
  use def_socrates_derived_fields, only: str_socrates_derived_fields
  use science_constants_mod, only : pi
  
  implicit none

contains
  
  subroutine solar_pos_calculation(socrates_opt, socrates_derived_fields)
    ! -------------------------------------------------------------------
    ! Calculations of the earth's orbit described in UMDP 23. Using  
    ! the day of the year and the orbital "constants" (which vary over           
    !  "Milankovitch" timescales) it calculates the sin of the solar     
    !  declination and the inverse-square scaling factor for the solar   
    !  "constant".  
    ! -----------------------------------------------------------------

    type (str_socrates_options), intent(in) :: socrates_opt
    type (str_socrates_derived_fields), intent(inout) :: socrates_derived_fields
                                                                    
     real(kind=DEFAULT_PRECISION) ::      &
           gamma_rad, e, tau0, sinobl, e1, e2, e3, e4, diny_360, & 
           tau1_360, tau1_365  
     real(kind=DEFAULT_PRECISION) :: diny_365              ! number of days in year      
     real(kind=DEFAULT_PRECISION) :: m, v                  ! mean & true anomaly        
                                    
     gamma_rad = 1.352631 
     e=.0167 
     tau0 = 2.5                    ! true date of perihelion  
     sinobl = 0.397789         ! sin (obliquity)    
     e1 = e * (2.0-0.25*e*e)    ! coefficients for 3.1.2
     e2 = 1.25 * e*e          ! coefficients for 3.1.2
     e3 = e*e*e * 13.0/12.0     
     e4 = ((1.0+e*e*0.5)/(1.0-e*e))**2.0    ! constant for 3.1.4 
     diny_360 = 360.                        ! number of days in year 
     tau1_360 = tau0*diny_360/365.25+0.71+.5
     tau1_365 = tau0+.5   
     
      if (.not. socrates_opt%l_360) then                                         
         if (mod(socrates_opt%rad_year,4.0) .eq. 0 .and.      &    ! is this a leap year?  
              (mod(socrates_opt%rad_year,400.0) .eq. 0 .or.   &
              mod(socrates_opt%rad_year,100.0) .ne. 0)) then      
            diny_365 = 366.                                             
         else       
            diny_365 = 365.    
         end if
      end if
                                                                         
      !  tau1 is modified so as to include the conversion of day-ordinal into    
      !  fractional-number-of-days-into-the-year-at-12-z-on-this-day.   
                                        
      if (socrates_opt%l_360) then              
        m = (2*pi) * (socrates_opt%rad_day-tau1_360) / diny_360                ! eq 3.1.1   
      else    
        m = (2*pi) * (socrates_opt%rad_day-tau1_365) / diny_365                ! eq 3.1.1   
      end if                  
      v = m + e1*sin(m) + e2*sin(2.*m) + e3*sin(3.*m)               ! eq 3.1.2   
      
      ! Fields required by monc-socrates coupling
      socrates_derived_fields%scs = e4 * ( 1. + e * cos(v) ) **2    ! eq 3.1.4    
      socrates_derived_fields%sindec = sinobl * sin (v - gamma_rad)     

  end subroutine solar_pos_calculation

  subroutine solar_angle_calculation(socrates_opt, socrates_derived_fields) 

    type (str_socrates_options), intent(in) :: socrates_opt
    type (str_socrates_derived_fields), intent(inout) :: socrates_derived_fields

    real(kind=DEFAULT_PRECISION), parameter ::   &
         degrees_to_radians = 0.017453292,         & ! conversion factor
         seconds_in_hour = 3600.0

     real(kind=DEFAULT_PRECISION) :: &
           sinlat,                    &  ! sin of latitude in radians
           longit                        ! longitude in radians
     
    real(kind=DEFAULT_PRECISION) :: &
         twopi,   &                  ! 2*pi                               
         s2r,  &                        ! seconds-to-radians converter       
         sinsin,   &         ! products of the sines and of the cosines   
         coscos,   &         ! of solar declination and of latitude.      
         hld,   &            ! half-length of the day in radians (equal   
         ! to the hour-angle of sunset, and minus     
         coshld,   &         ! the hour-angle of sunrise)  its cosine.   
         hat,   &            ! local hour angle at the start time.        
         omegab,   &         ! beginning and end of the timestep and      
         omegae,   &         ! of the period over which cosz is           
         omega1,   &         ! integrated, and sunset - all measured in   
         omega2,   &         ! radians after local sunrise, not from      
         omegas,   &         ! local noon as the true hour angle is.      
         difsin,   &         ! a difference-of-sines intermediate value   
         diftim,   &         ! and the corresponding time period
         ! the start-time and length of the timestep (time_radians & dt_radians)
         ! converted to radians after midday gmt, or equivalently, hour       
         ! angle of the sun on the greenwich meridian.
         time_radians,     & 
         dt_radians
    
    twopi = 2. * pi
    s2r = pi / 43200.

    sinlat = sin(socrates_opt%latitude*degrees_to_radians) 
    longit = socrates_opt%longitude*degrees_to_radians

    time_radians = &
         (socrates_opt%rad_time_hours*seconds_in_hour) * s2r - pi                                                  
    dt_radians =  socrates_derived_fields%dt_secs * s2r                                                     
    
    ! original code from LEM looped over k, which is number of points
    ! in the branch for version 1.0 of MONC, number of points will
    ! always be 1, so removed the loop
    hld = 0.                                ! logically unnecessary     
    ! statement without which the cray compiler will not vectorize this code   
    sinsin = socrates_derived_fields%sindec * sinlat                                         
    coscos = sqrt( (1.0- socrates_derived_fields%sindec**2.0) * &
         (1.0- sinlat**2.0) )                 
    coshld = sinsin / coscos                                            
    if (coshld.lt.-1.) then                 ! perpetual night           
       socrates_derived_fields%fraction_lit = 0.                                                      
       socrates_derived_fields%cosz = 0.                                                     
    else                                                               
       hat = longit + time_radians               ! (3.2.2)                   
       if (coshld.gt.1.) then               !   perpetual day - hour    
          omega1 = hat                      ! angles for (3.2.3) are    
          omega2 = hat + dt_radians              ! start & end of timestep   
       else                                !   at this latitude some   
          ! points are sunlit, some not.  different ones need different treatment.   
          hld = acos(-coshld)               ! (3.2.4)                   
          ! the logic seems simplest if one takes all "times" - actually hour        
          ! angles - relative to sunrise (or sunset), but they must be kept in the   
          ! range 0 to 2pi for the tests on their orders to work.                    
          omegab = hat + hld                                            
          if (omegab.lt.0.)   omegab = omegab + twopi                   
          if (omegab.ge.twopi) omegab = omegab - twopi                  
          if (omegab.ge.twopi) omegab = omegab - twopi                  
          !  line repeated - otherwise could have failure if            
          !  longitudes w are > pi rather than < 0.                     
          omegae = omegab + dt_radians                                       
          if (omegae.gt.twopi) omegae = omegae - twopi                  
          omegas = 2. * hld                                             
          ! now that the start-time, end-time and sunset are set in terms of hour    
          ! angle, can set the two hour-angles for (3.2.3).  the simple cases are    
          ! start-to-end-of-timestep, start-to-sunset, sunrise-to-end and sunrise-   
          ! -to-sunset, but two other cases exist and need special treatment.        
          if (omegab.le.omegas .or. omegab.lt.omegae) then              
             omega1 = omegab - hld                                      
          else                                                         
             omega1 = - hld                                             
          endif
          if (omegae.le.omegas) then                                    
             omega2 = omegae - hld                                      
          else                                                         
             omega2 = omegas - hld                                      
          endif
          if (omegae.gt.omegab.and.omegab.gt.omegas) omega2=omega1      
          !  put in an arbitrary marker for the case when the sun does not rise      
          !  during the timestep (though it is up elsewhere at this latitude).       
          !  (cannot set cosz & lit within the else ( coshld < 1 ) block             
          !  because 3.2.3 is done outside this block.)                              
       endif           ! this finishes the else (perpetual day) block   
       difsin = sin(omega2) - sin(omega1)             ! begin (3.2.3)   
       diftim = omega2 - omega1                                         
       ! next, deal with the case where the sun sets and then rises again         
       ! within the timestep.  there the integration has actually been done       
       ! backwards over the night, and the resulting negative difsin and diftim   
       ! must be combined with positive values representing the whole of the      
       ! timestep to get the right answer, which combines contributions from      
       ! the two separate daylit periods.  a simple analytic expression for the   
       ! total sun throughout the day is used.  (this could of course be used     
       ! alone at points where the sun rises and then sets within the timestep)   
       if (diftim.lt.0.) then                                           
          difsin = difsin + 2. * sqrt(1.-coshld**2)                      
          diftim = diftim + 2. * hld                                     
       endif
       if (diftim.eq.0.) then                                           
          ! pick up the arbitrary marker for night points at a partly-lit latitude   
          socrates_derived_fields%cosz = 0.                                                  
          socrates_derived_fields%fraction_lit = 0.                                                   
       else                                                            
          socrates_derived_fields%cosz = difsin*coscos/diftim + sinsin     ! (3.2.3)         
          socrates_derived_fields%fraction_lit = diftim / dt_radians                                       
       endif
    endif            ! this finishes the else (perpetual night) block   
    
  end subroutine solar_angle_calculation

end module solar_position_angle_mod
