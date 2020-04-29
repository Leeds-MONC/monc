module def_socrates_derived_fields

  use datadefn_mod, only : DEFAULT_PRECISION

  implicit none

  type str_socrates_derived_fields

     real(kind=DEFAULT_PRECISION) :: &
          dt_secs   ! radiation timestep in seconds. It will either be
                    ! socrates_opt%rad_int_time, or the MONC timestep
                    ! This is set in the socrates_couple, timestep_callback

     ! declare fields to use in the solar position calculation
     real(kind=DEFAULT_PRECISION) :: &
          default_solar_const,       &  ! solar constant (1365) read from config
          scs,                       &  ! solar constant scaling
          sindec,                    &  ! sin(solar declination)
          cosz
     ! decalre fields derived using output from
     ! solar_position and solar_angle
     real(kind=DEFAULT_PRECISION) :: &
          sol_const,                 &  ! default_solar_constant * scs
          sec_out,                   &  ! 1/cos(solar zenith angle)
          fraction_lit                  ! fraction lit 1 or 0 probably
     real(kind=DEFAULT_PRECISION) :: &
          albedoin1,                 &  ! Albedo between wavelenghs of ??
          albedoin2                     ! Albedo between wavelenghs of ??
     real(kind=DEFAULT_PRECISION) :: &
          srf_temperature

     ! declare density and radiation factor for heating rate calculation
     real(DEFAULT_PRECISION), allocatable :: &
          density_factor(:), radiation_factor(:)

     ! declare radiative heating rates (sw and lw).
     ! This is declared as a 3D array so that the heating rates can be
     ! stored for diagnostics or use when the radiation timestep is
     ! longer than the model timestep.
     ! Allocated in the initialisation callback
     real(DEFAULT_PRECISION), allocatable ::  &
          flux_up_sw(:,:,:),                  & ! shortwave flux up
          flux_down_sw(:,:,:),                & ! shortwave flux down
          flux_net_sw(:,:,:),                 & ! shortwave flux net
          flux_up_lw(:,:,:),                  & ! longwave flux up
          flux_down_lw(:,:,:),                & ! longwave flux down
          flux_net_lw(:,:,:),                 & ! longwave flux net
          swrad_hr(:,:,:),                    & ! shortwave heating rate
          lwrad_hr(:,:,:),                    & ! longwave heating rate
          totrad_hr(:,:,:)                      ! total radiative heating rate
          
     
     ! declare 2-d fields for shortwave and longwave toa and surface
     ! fluxes(sw and lw).
     real(DEFAULT_PRECISION), allocatable ::  &
         toa_up_longwave(:,:),                & ! top-of-atmosphere longwave up
         toa_down_shortwave(:,:),             & ! top-of-atmosphere shortwave down 
         toa_up_shortwave(:,:),               & ! top-of-atmosphere shortwave up
         surface_up_longwave(:,:),            & ! surface longwave up
         surface_down_longwave(:,:),          & ! surface longwave up
         surface_down_shortwave(:,:),         & ! surface shortwave down 
         surface_up_shortwave(:,:)              ! surface shortwave up

  end type str_socrates_derived_fields

end module def_socrates_derived_fields
