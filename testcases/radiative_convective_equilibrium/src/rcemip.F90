module rcemip_mod

  ! This module provides routines to assist in conforming to the RCEMIP specifications of
  ! Wing et al. (2018) Geosci. Model Dev., 11, 793-813, 2018 
  ! https://doi.org/10.5194/gmd-11-793-2018
  ! For use with testcases/radiative_convective_equilibrium/RCEMIP.mcf
       
  use datadefn_mod, only : DEFAULT_PRECISION
  use monc_component_mod, only : component_descriptor_type
  use state_mod, only : model_state_type
  use grids_mod, only : X_INDEX, Y_INDEX, Z_INDEX
  use science_constants_mod, only : r_over_cp, r, g
  use optionsdatabase_mod, only :  options_get_real_array, options_get_real, &
     options_get_logical, options_get_integer, options_get_array_size, options_get_string_array
  use logging_mod, only: log_master_log, LOG_ERROR, LOG_INFO, log_is_master
  use q_indices_mod, only: get_q_index, standard_q_names
  use def_merge_atm, only: str_merge_atm

  implicit none

#ifndef TEST_MODE
  private
#endif
  
  public rcemip_get_descriptor, rcemip_init, rcemip_ozone
contains

  type(component_descriptor_type) function rcemip_get_descriptor()
    rcemip_get_descriptor%name="rcemip"
    rcemip_get_descriptor%version=0.1
    rcemip_get_descriptor%initialisation=>initialisation_callback
    rcemip_get_descriptor%timestep=>timestep_callback
  end function rcemip_get_descriptor

  !! Note that this is not the most efficient way to iterate through theta (j heavy), but it is the same as the LEM set up 
  !! so directly comparable and probably doesn't matter too much as it is just called onec in the initialisation
  subroutine initialisation_callback(current_state)
    type(model_state_type), intent(inout), target :: current_state

  end subroutine initialisation_callback

  subroutine timestep_callback(current_state)
    type(model_state_type), intent(inout), target :: current_state

  end subroutine timestep_callback

  !> Called by gridmanager.F90's calculate_initial_profiles routine for proper order placement
  !! Creates RCEMIP analytic sounding approximating the moist tropical sounding of Dunion (2011)
  !! @param current_state The current model state_mod
  subroutine rcemip_init(current_state)
    type(model_state_type), intent(inout) :: current_state

    real(kind=DEFAULT_PRECISION), allocatable :: zngrid(:), qv(:), tv(:), tabs(:), p(:)
    integer :: i, j, k, nz, iq

    ! Analytic sounding paramters
    real(kind=DEFAULT_PRECISION), parameter :: &
      zt = 15000.0_DEFAULT_PRECISION, &      ! [ m ] (approximate tropopause height)
      q0_295 = 12.00e-3_DEFAULT_PRECISION, & ! [ kg/kg ] (specific humidity)
      q0_300 = 18.65e-3_DEFAULT_PRECISION, & ! [ kg/kg ] (specific humidity)
      q0_305 = 24.00e-3_DEFAULT_PRECISION, & ! [ kg/kg ] (specific humidity)
      qt = 1E-14_DEFAULT_PRECISION, &        ! [ kg/kg ] (specific humidity)
      zq1 = 4000.0_DEFAULT_PRECISION, &      ! [ m ]
      zq2 = 7500.0_DEFAULT_PRECISION, &      ! [ m ]
      gamma = 0.0067_DEFAULT_PRECISION, &    ! [ K/m ]
      p0 = 1014.8_DEFAULT_PRECISION, &       ! [ hPa ]
      const = 0.608_DEFAULT_PRECISION
    real(kind=DEFAULT_PRECISION) :: q0, tv0, tvt, pt
    integer :: sst
    logical :: l_matchthref

    call log_master_log(LOG_INFO, "RCEMIP-specified analytical profiles will be applied.")

    l_matchthref = options_get_logical(current_state%options_database, "l_matchthref")
    sst = options_get_integer(current_state%options_database, "rcemip_sst")

    ! Set surface values (specific humidity, virtual temperature)
    if (sst .eq. 295) q0 = q0_295
    if (sst .eq. 300) q0 = q0_300
    if (sst .eq. 305) q0 = q0_305
    tv0 = real(sst,DEFAULT_PRECISION) * (1.0_DEFAULT_PRECISION + const * q0)
    tvt = tv0 - gamma * zt  ! (virtual temperature at zt)
    pt = p0 * (tvt/tv0) ** (G/(r*gamma))

    nz = current_state%global_grid%size(Z_INDEX)
    allocate(zngrid(nz),qv(nz),tv(nz),tabs(nz),p(nz))
    zngrid(:)=current_state%global_grid%configuration%vertical%zn(:)

    ! Compute profiles
    do k=2,nz
      if (zngrid(k) .ge. 0.0_DEFAULT_PRECISION .and. zngrid(k) .le. zt) then
        qv(k) = q0 * exp(-zngrid(k)/zq1) * exp(-(zngrid(k)/zq2)**2)  ! (specific humidity)
        tv(k) = tv0 - gamma * zngrid(k)
        p(k)  = p0*((tv0 -(gamma  *zngrid(k)))/tv0)**(G/(r *gamma))
      else if (zngrid(k) .gt. zt) then
        qv(k) = qt                                                   ! (specific humidity)
        tv(k) = tvt
        p(k)  = pt * exp(-( G*(zngrid(k)-zt)/(r*tvt) ))
      end if
      tabs(k) = tv(k) / (1.0_DEFAULT_PRECISION + const * qv(k))
    end do
    qv(1)=qv(2)
    tv(1)=tv(2)
    p(1)=p(2)
    tabs(1)=tabs(2)

    ! Conversion: q = w/(1+w)   <===>   w = q/(1-q)
    ! Specific humidity to mixing ratio, store in q_init
    iq = get_q_index("vapour","rcemip_mod")
    current_state%global_grid%configuration%vertical%q_init(:, iq) &
                                                         = qv(:)/(1.0_DEFAULT_PRECISION - qv(:))
    ! Convert absolute temperature to potential temperature
    !   Pressure used in mb [hPa]
    current_state%global_grid%configuration%vertical%theta_init(:) &
                                      = tabs(:)*(1000.0_DEFAULT_PRECISION/p(:))**r_over_cp

    ! Handle mixing ratio initialisation
    if (.not. current_state%continuation_run) then
      do i=current_state%local_grid%local_domain_start_index(X_INDEX), &
           current_state%local_grid%local_domain_end_index(X_INDEX)
        do j=current_state%local_grid%local_domain_start_index(Y_INDEX), &
             current_state%local_grid%local_domain_end_index(Y_INDEX)
          current_state%q(iq)%data(:,j,i) &
                                = current_state%global_grid%configuration%vertical%q_init(:, iq)
        end do
      end do
    end if

    ! Handle theta initialisation
    if (.not. current_state%continuation_run) then
      if (l_matchthref) then 
        if(.not. current_state%use_anelastic_equations) then
          call log_master_log(LOG_ERROR, "Non-anelastic equation set and l_maththref are incompatible")
       end if
       current_state%global_grid%configuration%vertical%thref = current_state%global_grid%configuration%vertical%theta_init
      end if
      
      ! Fill theta data
      do i=current_state%local_grid%local_domain_start_index(X_INDEX), current_state%local_grid%local_domain_end_index(X_INDEX)
        do j=current_state%local_grid%local_domain_start_index(Y_INDEX), current_state%local_grid%local_domain_end_index(Y_INDEX)
            current_state%th%data(:,j,i) = &
                 current_state%global_grid%configuration%vertical%theta_init(:) - &
                 current_state%global_grid%configuration%vertical%thref(:)
        end do
      end do
    end if
  end subroutine rcemip_init


  ! Routine to set the ozone profile by overwriting mcc values.
  ! @param merge_fields: merged model/mcc fields for SOCRATES
  subroutine rcemip_ozone(merge_fields)
    type (str_merge_atm), intent(inout) :: merge_fields

    real(kind=DEFAULT_PRECISION), parameter :: &
      mmr_fac = 1e-6, &     ! conversion factor (ppmv --> mmr)
      hPa_fac = 1e-2, &     ! conversion factor (Pa --> hPa)
      g1 = 3.6478,    &     ! O_3 fit parameter g1, [ppmv hPa**-g2]
      g2 = 0.83209,   &     ! O_3 fit parameter g2, []
      g3 = 11.3515,   &     ! O_3 fit parameter g3, [hPa]
      mma = 28.97,    &     ! mean molar mass of dry air
      mmo = 47.997          ! mean molar mass of ozone
    integer :: k

    do k=1,size(merge_fields%pres_n) !size is mcc%irad_levs
      merge_fields%o3_n(k) = (g1 * ((merge_fields%pres_n(k)*hPa_fac)**g2) &
                                      * exp(-((merge_fields%pres_n(k)*hPa_fac)/g3)) ) &
                               * mmr_fac * (mmo/mma)  ! convert from ppmv to mmr
    end do

  end subroutine rcemip_ozone

end module rcemip_mod
