module def_merge_atm

  use state_mod, only : model_state_type
  use datadefn_mod, only : DEFAULT_PRECISION
  use grids_mod, only : Z_INDEX
  use def_mcc_profiles, only: str_mcc_profiles
  
  implicit none

  type Str_merge_atm
     
     ! atmospheric profiles with dimensions irad_levs, i.e. monc + mcc_cut and
     ! above.
     ! Note: All variables with suffix _n exist in the centre of the grid, while 
     !       all variables with suffix _level exist at the grid boundary
     real(kind=DEFAULT_PRECISION), allocatable ::    &
          qv_n(:), t_n(:)                  & 
          , ql_n(:), qi_n(:)               &
          , pres_n(:), o3_n(:)             &
          , qv_level(:), t_level(:)        & 
          , ql_level(:), qi_level(:)       &
          , pres_level(:), o3_level(:)     &    
          , mass(:)
     ! declare number fields if available
     real(kind=DEFAULT_PRECISION), allocatable ::    &
          cloudnumber_n(:)
     ! Local fields used in merge_atm_fields, with dimension of monc zindex. These are calced
     ! prior to the merge of data
     real(kind=DEFAULT_PRECISION), allocatable ::    &     
          pref_loc(:), t_level_loc(:), t_n_loc(:)
     ! define cloud fractions. This is a assumed to be
     ! 0 or 1, depending on the presence of total water. If using a
     ! fractional cloud scheme, this will need to be changed
     real(kind=DEFAULT_PRECISION), allocatable ::    &     
          liquid_cloud_fraction(:), ice_cloud_fraction(:), total_cloud_fraction(:)
     ! longwave and shortwave heating rates from radiation levels
     real(kind=DEFAULT_PRECISION), allocatable ::    &
          lw_heat_rate_radlevs(:), sw_heat_rate_radlevs(:)
     
     ! The factors below were derived as part of J. Petch PhD
     real(kind=DEFAULT_PRECISION) ::               &
          rainfac=0.02,                            &
          snowfac=0.4,                             &
          graupfac=0.05

  end type Str_merge_atm

contains
  
  subroutine allocate_merge_data_fields(current_state, merge_fields, mcc)

    implicit none
    
    type(model_state_type), target, intent(inout) :: current_state
    type (str_mcc_profiles), intent(in)  :: mcc
    type (str_merge_atm), intent(inout) :: merge_fields

    ! monc fields with the monc vertical grid, derived before merge of McClatchey and monc domains
    allocate(merge_fields%pref_loc(current_state%local_grid%size(Z_INDEX)))
    allocate(merge_fields%t_n_loc(current_state%local_grid%size(Z_INDEX)))
    allocate(merge_fields%t_level_loc(current_state%local_grid%size(Z_INDEX)))

    ! merged fields at the centre of the grid (zn layer)
    allocate(merge_fields%t_n(mcc%irad_levs))         ! absolute temperature
    allocate(merge_fields%qv_n(mcc%irad_levs))        ! vapour mixing ratio
    allocate(merge_fields%ql_n(mcc%irad_levs))        ! liquid water mixing ratio
    allocate(merge_fields%qi_n(mcc%irad_levs))        ! ice mass mixing ratio
    allocate(merge_fields%pres_n(mcc%irad_levs))      ! pressure
    allocate(merge_fields%o3_n(mcc%irad_levs))        ! ozone
    allocate(merge_fields%cloudnumber_n(mcc%irad_levs))    ! cloud drop number concentration
    

    allocate(merge_fields%total_cloud_fraction(mcc%irad_levs))  ! liquid+ice fraction
    allocate(merge_fields%liquid_cloud_fraction(mcc%irad_levs)) ! liquid fraction
    allocate(merge_fields%ice_cloud_fraction(mcc%irad_levs))    ! ice fraction

    ! merged fields at the layer boundary (z layer)
    allocate(merge_fields%t_level(0:mcc%irad_levs))       ! temperature on level
    allocate(merge_fields%qv_level(mcc%irad_levs))        ! vapour mixing ratio
    allocate(merge_fields%ql_level(mcc%irad_levs))        ! liquid water mixing ratio
    allocate(merge_fields%qi_level(mcc%irad_levs))        ! ice mass mixing ratio
    allocate(merge_fields%pres_level(0:mcc%irad_levs))      ! pressure
    allocate(merge_fields%o3_level(mcc%irad_levs))        ! ozone


    allocate(merge_fields%mass(mcc%irad_levs))            ! mass of the atmos at each
                                                          ! level
    
    allocate(merge_fields%lw_heat_rate_radlevs(mcc%irad_levs))
    allocate(merge_fields%sw_heat_rate_radlevs(mcc%irad_levs))
    
    merge_fields%pref_loc(:) = 0.0
    merge_fields%t_n_loc(:) = 0.0
    merge_fields%t_level_loc(:) = 0.0

    merge_fields%t_n(:) = 0.0         ! absolute temperature
    merge_fields%qv_n(:) = 0.0        ! vapour mixing ratio
    merge_fields%ql_n(:) = 0.0        ! liquid water mixing ratio
    merge_fields%qi_n(:) = 0.0        ! rain mass mixing ratio
    merge_fields%pres_n(:) = 0.0      ! pressure
    merge_fields%o3_n(:) = 0.0

    merge_fields%t_level(:) = 0.0         ! absolute temperature
    merge_fields%qv_level(:) = 0.0        ! vapour mixing ratio
    merge_fields%ql_level(:) = 0.0        ! liquid water mixing ratio
    merge_fields%qi_level(:) = 0.0        ! rain mass mixing ratio
    merge_fields%pres_level(:) = 0.0  ! pressure
    merge_fields%o3_level(:) = 0.0

    merge_fields%total_cloud_fraction(:)  = 0.0  ! liquid+ice fraction
    merge_fields%liquid_cloud_fraction(:) = 0.0  ! liquid fraction
    merge_fields%ice_cloud_fraction(:)    = 0.0  ! ice fraction

    merge_fields%mass(:) = 0.0

    merge_fields%lw_heat_rate_radlevs(:) = 0.0
    merge_fields%sw_heat_rate_radlevs(:) = 0.0

  end subroutine allocate_merge_data_fields

end module def_merge_atm
