! Note that this subroutine replaces part of the profile initialisation
! as well as setfluxlook, coriolis, and "forcing from mcf" routines

! NOTES
! - Currently, surface pressure needs to be set in mcf still.
! - GW Damping and grid also need to be set correctly in MCF
! - Module needs to be initialised after gridmanager but before random noise
!   This may need to be checked
! - Ensure PW advection does not clear source terms.
! - Note q in MONC is mixing ratio (rather than specific humidity, as is more usual)

! CURRENTLY TESTING
! - Handle nudging above clouds in prescribed fashion (could be a question for DEPHY community)?
! - Implement lat/lon dependence for radiation (socrates_opt%latitude,socrates_opt%longitude,socrates_opt%surface_albedo)

! TODO
! - Surface pressure initialisation from file?
! - Implement consistency check for use_surface_boundary_conditions flag
! - Implement check that grid manager initialised but random noise hasn't been applied yet
! - Check for possible problematic nature of modifying both current state and vertical grid simultaneously
!   and check what "target" keyword does in this context.
! - Handle (evolving) surface pressure?
! - Deal with surface non-zero height above sea level (better to deal with this in Lagtraj?)
! - Make code self-documenting with Doxygen
! - Add diagnostics
! - Improve interpolation routines? (replace linear interpolation by Steffen interpolation)
! - Code up finalisation callback (deallocation)?
! - Implement a less hacky column mode check?
! - Code restructuring/reformatting/more DRY code
! - Discuss best way to refer to z-coordinates
! - Discuss best way to do lowerbc and z0/z0th reinitialisation
! - Check grid vertical halo size

! TODO: RELATED ISSUES
! - Discuss possible problematic nature of modifying both current_state and current_state%vertical grid simultaneously in general
! - Work on problems with "entire domain" setfluxlook for heterogeneous surface forcings.
! - More systematic implementation of utilities in a separate module.

module dephy_forcings_mod
  use datadefn_mod, only : STRING_LENGTH
  use monc_component_mod, only : component_descriptor_type
  ! Use PRESCRIBED_SURFACE_FLUXES, PRESCRIBED_SURFACE_VALUES to specify boundary conditions for use in other modules
  use state_mod, only : model_state_type, PRESCRIBED_SURFACE_FLUXES, PRESCRIBED_SURFACE_VALUES
  use optionsdatabase_mod, only : options_get_integer, options_get_logical, options_get_real, &
      options_get_string
  use grids_mod, only : vertical_grid_configuration_type, X_INDEX, Y_INDEX, Z_INDEX
  use logging_mod, only : LOG_ERROR, log_master_log, log_log
  ! note z0 and z0th are overwritten during the simulation
  use science_constants_mod, only : cp, rlvap, z0, z0th, G, von_karman_constant, ratio_mol_wts,r_over_cp,&
      alphah,betah,betam,gammah,gammam,rlvap_over_cp
  use q_indices_mod, only: get_q_index, standard_q_names
  use interpolation_mod, only: piecewise_linear_1d, interpolate_point_linear_1d, interpolate_point_linear_2d
  use registry_mod, only : is_component_enabled
  use netcdf, only : nf90_noerr, nf90_global, nf90_nowrite,    &
       nf90_inquire_attribute, nf90_open, nf90_strerror,       &
       nf90_inq_dimid, nf90_inquire_dimension, nf90_inq_varid, &
       nf90_get_var, nf90_inquire, nf90_close, nf90_get_att, &
       nf90_ebaddim, nf90_enotatt, nf90_enotvar, nf90_inquire_attribute
  use configuration_checkpoint_netcdf_parser_mod, only : remove_null_terminator_from_string
  ! use existing fluxlook functionality
  use setfluxlook_mod, only : set_look, change_look
  ! re-initialise lowerbc module as z0 and z0th change over time
  use lowerbc_mod,  only: tstrcona, rhmbc, ddbc, ddbc_x4, eecon, r2ddbc, rcmbc, tstrconb, &
       x4con, xx0con, y2con, yy0con, viscous_courant_coefficient
  use saturation_mod, only : qsaturation
  ! supersede initialisation from mcf file in grid manager
  use gridmanager_mod, only : set_up_vertical_reference_properties,set_anelastic_pressure, &
    setup_reference_state_liquid_water_temperature_and_saturation, &
    calculate_mixing_length_for_neutral_case, set_buoyancy_coefficient
  use datadefn_mod, only : DEFAULT_PRECISION, PRECISION_TYPE
  use mpi, only : MPI_SUM, MPI_IN_PLACE
  use socrates_couple_mod, only: socrates_opt
  use def_socrates_options, only: str_socrates_options

  implicit none

#ifndef TEST_MODE
  private
#endif

  real(kind=DEFAULT_PRECISION), allocatable :: time_dephy(:)
  real(kind=DEFAULT_PRECISION), allocatable :: height_dephy(:)
  real(kind=DEFAULT_PRECISION), allocatable :: module_z(:)
  real(kind=DEFAULT_PRECISION), allocatable :: module_zn(:)
  real(kind=DEFAULT_PRECISION), allocatable :: full_theta(:,:,:)
  real(kind=DEFAULT_PRECISION), parameter :: proper_pi=atan(1.0_DEFAULT_PRECISION) * 4.0_DEFAULT_PRECISION

  character(len=STRING_LENGTH) :: dephy_file
  integer :: ncid_dephy
  integer :: time_len_dephy
  integer :: height_len_dephy
  integer :: kkp !module wide parameter for vertical grid
  logical :: l_verbose=.false. ! Temporary flag for dirty debugging
  ! Three parameters below are meant to catch the model being in column mode.
  integer :: column_check_x !
  integer :: column_check_y !
  integer :: n_dephy_passes=0 ! start checking after the initial two dephy passes

  ! Surface fields (time-dependent)
  real(kind=DEFAULT_PRECISION), allocatable :: lat_traj_dephy(:), &
    lon_traj_dephy(:), &
    ps_forc_dephy(:), &
    ts_dephy(:), &
    sfc_sens_flx_dephy(:), &
    sfc_lat_flx_dephy(:), &
    z0_traj_dephy(:), &
    z0th_traj_dephy(:), &
    ustar_dephy(:), &
    u_traj_dephy(:), &
    v_traj_dephy(:), &
    albedo_traj_dephy(:), &
    q_skin_traj_dephy(:)

  ! Surface fields at time step
  real(kind=DEFAULT_PRECISION) :: lat_traj, &
    lon_traj, &
    ps_forc, &
    ts, &
    sfc_sens_flx, &
    sfc_lat_flx, &
    z0_traj, &
    z0th_traj, &
    ustar, &
    u_traj, &
    v_traj, &
    albedo_traj, &
    q_skin_traj

  ! Initial fields on MONC vertical grid
  real(kind=DEFAULT_PRECISION), allocatable :: u_dephy(:), &
    v_dephy(:), &
    theta_dephy(:), &
    rv_dephy(:), &
    tke_dephy(:)

  ! Forcing fields on MONC vertical grid (time-dependent)
  real(kind=DEFAULT_PRECISION), allocatable :: height_forc_dephy(:,:), &
    pressure_forc_dephy(:,:), &
    ug_dephy(:,:), &
    vg_dephy(:,:), &
    u_adv_dephy(:,:), &
    v_adv_dephy(:,:), &
    theta_adv_dephy(:,:), &
    theta_rad_dephy(:,:), &
    rv_adv_dephy(:,:), &
    w_dephy(:,:), &
    theta_nudging_dephy(:,:), &
    rv_nudging_dephy(:,:), &
    u_nudging_dephy(:,:), &
    v_nudging_dephy(:,:), &
    nudging_inv_u_traj_dephy(:,:), &
    nudging_inv_v_traj_dephy(:,:), &
    nudging_inv_theta_traj_dephy(:,:), &
    nudging_inv_rv_traj_dephy(:,:)

  ! Forcing fields during a time step
  real(kind=DEFAULT_PRECISION), allocatable :: height_forc(:), &
    pressure_forc(:), &
    ug(:), &
    vg(:), &
    u_adv(:), &
    v_adv(:), &
    theta_adv(:), &
    theta_rad(:), &
    rv_adv(:), &
    w(:), &
    theta_nudging(:), &
    rv_nudging(:), &
    u_nudging(:), &
    v_nudging(:), &
    nudging_inv_u_traj(:), &
    nudging_inv_v_traj(:), &
    nudging_inv_theta_traj(:), &
    nudging_inv_rv_traj(:)

  ! Dephy flags: applied only during simulation
  integer :: int_adv_theta, &
        int_adv_rv, &
        int_rad_theta, &
        int_forc_w, &
        int_forc_geo, &
        int_nudging_u, &
        int_nudging_v, &
        int_nudging_theta, &
        int_nudging_rv

  ! EUREC4A nudging procedure
  integer :: int_inversion_nudging=0 ! Flag for inversion nudging (optional DEPHY extension)
  real(kind=DEFAULT_PRECISION), allocatable :: theta_l(:,:,:)
  real(kind=DEFAULT_PRECISION), allocatable :: theta_l_mean(:)
  real(kind=DEFAULT_PRECISION) :: inversion_nudging_height_above
  real(kind=DEFAULT_PRECISION) :: inversion_nudging_transition
  real(kind=DEFAULT_PRECISION) :: inversion_nudging_time

  ! Dephy strings
  character(len=STRING_LENGTH):: str_surfaceType, &
        str_surfaceForcing, &
        str_surfaceForcingWind

  public dephy_forcings_get_descriptor

contains

  !> Provides the descriptor back to the caller and is used in component registration
  !! @returns The dephy_forcings component descriptor
  type(component_descriptor_type) function dephy_forcings_get_descriptor()
    dephy_forcings_get_descriptor%name="dephy_forcings"
    dephy_forcings_get_descriptor%version=0.1
    dephy_forcings_get_descriptor%initialisation=>initialise_callback
    dephy_forcings_get_descriptor%timestep=>timestep_callback
  end function dephy_forcings_get_descriptor


  !> Called during initialisation and will initialise the horizontal and vertical grid configurations
  !! Note that the model state_mod (from a checkpoint or external file) must have been initialised already
  !! @param current_state The current model state_mod
  subroutine initialise_callback(current_state)
    implicit none
    type(model_state_type), target, intent(inout) :: current_state
    type(vertical_grid_configuration_type) :: vertical_grid
    integer :: alloc_z, alloc_y, alloc_x

    !if(current_state%parallel%my_rank==0) then
    !    l_verbose=.true.
    !endif

    if(l_verbose) write(*,*) "initialising dephy"

    alloc_z=current_state%local_grid%size(Z_INDEX) + current_state%local_grid%halo_size(Z_INDEX) * 2
    alloc_y=current_state%local_grid%size(Y_INDEX) + current_state%local_grid%halo_size(Y_INDEX) * 2
    alloc_x=current_state%local_grid%size(X_INDEX) + current_state%local_grid%halo_size(X_INDEX) * 2

    !!! Might need re-thinking: vertical grid passed into subroutines as well as current state
    if (.not. current_state%initialised) then
      call log_log(LOG_ERROR, "Must initialise the model state_mod before constructing the grid properties")
    end if

    vertical_grid=current_state%global_grid%configuration%vertical

    ! get DEPHY filename, which needs to be trimmed
    dephy_file=options_get_string(current_state%options_database, "dephy_file")
    kkp=current_state%local_grid%size(Z_INDEX)

    ! allocate module arrays
    allocate(module_z(kkp))
    allocate(module_zn(kkp))
    allocate(full_theta(alloc_z, alloc_y, alloc_x))
    allocate(height_forc(kkp))
    allocate(pressure_forc(kkp))
    allocate(ug(kkp))
    allocate(vg(kkp))
    allocate(u_adv(kkp))
    allocate(v_adv(kkp))
    allocate(theta_adv(kkp))
    allocate(theta_rad(kkp))
    allocate(rv_adv(kkp))
    allocate(w(kkp))
    allocate(theta_nudging(kkp))
    allocate(rv_nudging(kkp))
    allocate(u_nudging(kkp))
    allocate(v_nudging(kkp))
    allocate(nudging_inv_u_traj(kkp))
    allocate(nudging_inv_v_traj(kkp))
    allocate(nudging_inv_theta_traj(kkp))
    allocate(nudging_inv_rv_traj(kkp))
    allocate(theta_l_mean(kkp))
    allocate(theta_l(alloc_z, alloc_y, alloc_x))

    if(l_verbose) write(*,*) "initialised dephy 1"

    module_z=vertical_grid%z(:)
    module_zn=vertical_grid%zn(:)

    if(l_verbose) write(*,*) "initialised dephy 2"

    call check_status(nf90_open(path = trim(dephy_file), mode = nf90_nowrite, ncid = ncid_dephy))

    if(l_verbose) write(*,*) "initialised dephy 3"

    call dephy_read_dimension_variables() ! reads the forcing time and height variables
    call dephy_read_profile_variables() ! does profile initialisation
    call dephy_read_forcing_variables() ! reads and interpolates forcings
    call dephy_read_surface_variables() ! reads and interpolates forcings
    call dephy_read_integers() ! reads flags
    call dephy_read_strings() ! reads strings
    call dephy_read_inversion_nudging()

    if(l_verbose) write(*,*) "initialised dephy 4"

    call dephy_sanity_checks(current_state) ! checks for incompatible elements

    if(l_verbose) write(*,*) "initialised dephy 5"

    call dephy_time_interpolate(current_state) ! needed to get such things as z0

    if(l_verbose) write(*,*) "initialised dephy 6"

    call dephy_profiles_etc(current_state, vertical_grid) ! initialises proiles/reference profiles

    if(l_verbose) write(*,*) "initialised dephy 7"

    call dephy_setfluxlook_init(current_state) ! initialises surface

    if(l_verbose) write(*,*) "initialised dephy 8"

    call check_status(nf90_close(ncid_dephy))

    if(l_verbose) write(*,*) "initialised dephy 10"

    if(l_verbose) call dephy_bughunting(current_state)

  end subroutine initialise_callback


  subroutine timestep_callback(current_state)
    implicit none
    type(model_state_type), target, intent(inout) :: current_state

    ! update forcing and current fields
    call dephy_time_interpolate(current_state)
    call dephy_column_mode_check(current_state)

    ! Reset z0 and z0th, and related parameters
    z0=z0_traj
    z0th=z0th_traj
    current_state%global_grid%configuration%vertical%zlogm=&
      log(1.0_DEFAULT_PRECISION+current_state%global_grid%configuration%vertical%zn(2)/z0)
    current_state%global_grid%configuration%vertical%zlogth=&
      log((current_state%global_grid%configuration%vertical%zn(2)+z0)/z0th)
    current_state%global_grid%configuration%vertical%vk_on_zlogm=&
      von_karman_constant/current_state%global_grid%configuration%vertical%zlogm

    if(l_verbose) write(*,*) "dephy timestep 1"

    call dephy_setfluxlook_timestep(current_state)

    if(l_verbose) write(*,*) "dephy timestep 2"

    call dephy_apply_forcings(current_state)

    if(l_verbose) write(*,*) "dephy timestep 3"

    ! re-initialisation needed as z0 and z0th means lowerbc needs reinitialisation
    call lowerbc_reset_constants(current_state)

    if(l_verbose) write(*,*) "dephy timestep 4"

    if(l_verbose) call dephy_bughunting(current_state)
    if(l_verbose) call dephy_dirty_diagnostics()

  end subroutine timestep_callback


  subroutine dephy_sanity_checks(current_state)
    implicit none
    type(model_state_type), target, intent(inout) :: current_state

    !! Check for initialisation from mcf file (incompatible)
    logical :: l_init_pl_u     ! if .true. then initialize u field from mcf file
    logical :: l_init_pl_v     ! if .true. then initialize v field
    logical :: l_init_pl_theta ! if .true. then initialize potential temperature field
    logical :: l_init_pl_temp  ! if .true. then initialize temperature field
    logical :: l_init_pl_rh    ! if .true. then initialize relative humidity field
    logical :: l_init_pl_q     ! if .true. then initialize q fields
    real(kind=DEFAULT_PRECISION) :: termination_time, zztop, max_height_cloud

    !! Check if incompatible routines are enabled
    !! Check if incompatible flags are set

    !! Check "standard" coriolis force not activated
    !! Check "other" forcings routine not acticated
    !! Check mean profiles are enabled if nudging is used
    !! Check radiation scheme is compatible (ACTIVATED if 0, DEACTIVATED otherwise)
    !! Check initialisation from profiles not activated
    !! Check not current_state%passive_q .or. current_state%passive_th
    !! Check current_state%number_q_fields > 0
    !! Check for inconsistent surface conditions. Make sure setfluxlook component is not active
    !! Catch ustar-based setups (check earlier work on BOMEX for potential fix).

    if (is_component_enabled(current_state%options_database, "setfluxlook")) then
      call log_master_log(LOG_ERROR, "DEPHY: setfluxlook component incompatible with dephy forcing")
    endif
    if (is_component_enabled(current_state%options_database, "forcing")) then
      call log_master_log(LOG_ERROR, "DEPHY: forcing component incompatible with dephy forcing")
    endif
    if(is_component_enabled(current_state%options_database, "socrates_couple")) then
    if(int_rad_theta==1) then
      call log_master_log(LOG_ERROR, "DEPHY: socrates_couple component incompatible with dephy flag rad_theta==1")
    endif
    endif
    if(.not. is_component_enabled(current_state%options_database, "socrates_couple")) then
    if(int_rad_theta==0) then
      call log_master_log(LOG_ERROR, "DEPHY: absence of socrates_couple component incompatible with dephy flag rad_theta==0")
    endif
    endif
    if (is_component_enabled(current_state%options_database, "lwrad_exponential")) then
      call log_master_log(LOG_ERROR, "DEPHY: lwrad_exponential component incompatible with dephy forcing")
    endif
    if (.not. is_component_enabled(current_state%options_database, "buoyancy")) then
      call log_master_log(LOG_ERROR, "DEPHY: absence of buoyancy component incompatible with dephy forcing")
    endif
    if (.not. is_component_enabled(current_state%options_database, "lower_bc")) then
      call log_master_log(LOG_ERROR, "DEPHY: absence of lower_bc component incompatible with dephy forcing")
    endif
    if (.not. is_component_enabled(current_state%options_database, "set_consistent_lowbc")) then
      call log_master_log(LOG_ERROR, "DEPHY: absence of set_consistent_lowbc component incompatible with dephy forcing")
    endif
    if (.not. is_component_enabled(current_state%options_database, "mean_profiles")) then
      call log_master_log(LOG_ERROR, "DEPHY: absence of mean_profiles component incompatible with dephy forcing")
    endif
    if(.not. (trim(str_surfaceForcing)=="surfaceFlux" .or. trim(str_surfaceForcing)=="ts")) then
      call log_master_log(LOG_ERROR, "DEPHY: surfaceForcing (thermodynamics) not implemented")
    endif
    if(.not. trim(str_surfaceForcingWind)=="z0_traj") then
      call log_master_log(LOG_ERROR, "DEPHY: surfaceForcingWind not implemented")
    endif
    if(.not. current_state%number_q_fields > 0) then
      call log_master_log(LOG_ERROR, "DEPHY: dephy_forcings need current_state%number_q_fields > 0")
    endif
    if(current_state%passive_q) then
      call log_master_log(LOG_ERROR, "DEPHY: dephy_forcings incompatible with passive_q")
    endif
    if(current_state%passive_th) then
      call log_master_log(LOG_ERROR, "DEPHY: dephy_forcings incompatible with passive_th")
    endif
    l_init_pl_theta=options_get_logical(current_state%options_database, "l_init_pl_theta")
    l_init_pl_temp=options_get_logical(current_state%options_database, "l_init_pl_temp")
    l_init_pl_rh=options_get_logical(current_state%options_database, "l_init_pl_rh")
    l_init_pl_q=options_get_logical(current_state%options_database, "l_init_pl_q")
    l_init_pl_u=options_get_logical(current_state%options_database, "l_init_pl_u")
    l_init_pl_v=options_get_logical(current_state%options_database, "l_init_pl_v")
    if(l_init_pl_theta .or. l_init_pl_temp .or. l_init_pl_rh .or. l_init_pl_q .or. &
       l_init_pl_u .or. l_init_pl_v) then
      call log_master_log(LOG_ERROR, &
      "DEPHY: dephy_forcings incompatible with initialisation of profiles using "//&
      "l_init_pl_theta or l_init_pl_temp or l_init_pl_rh or l_init_pl_q or l_init_pl_u or l_init_pl_v")
    endif
    termination_time=options_get_real(current_state%options_database, "termination_time")
    if(termination_time>time_dephy(time_len_dephy)) then
      call log_master_log(LOG_ERROR, "DEPHY: termination time beyond last time in forcing file")
    endif
    zztop=options_get_real(current_state%options_database, "zztop")
    if(zztop>height_dephy(height_len_dephy)) then
      call log_master_log(LOG_ERROR, "DEPHY: zztop beyond highest level in forcing file")
    endif
    max_height_cloud=options_get_real(current_state%options_database, "max_height_cloud")
    if(max_height_cloud<zztop) then
      call log_master_log(LOG_ERROR, "DEPHY: DANGER max_height_cloud below highest level in domain")
    endif

  end subroutine dephy_sanity_checks

  subroutine dephy_column_mode_check(current_state)
    implicit none
    type(model_state_type), target, intent(inout) :: current_state

    ! perform check only after second time-step fully completed
    ! since current column may be different in first time step.
    if(n_dephy_passes<2) then
      column_check_x=current_state%column_local_x
      column_check_y=current_state%column_local_y
      n_dephy_passes=n_dephy_passes+1
    else
      if(.not. (column_check_x==current_state%column_local_x .and. &
      column_check_y==current_state%column_local_y)) then
        call log_master_log(LOG_ERROR, "DEPHY: WAYWARD FRIAR ERROR. MONC SEEMS TO BE RUNNING IN COLUMN MODE WHILE RUNNING DEPHY!")
      end if
      n_dephy_passes=n_dephy_passes+1
    endif

  end subroutine dephy_column_mode_check

  subroutine dephy_read_dimension_variable(field,netcdf_name,dim_length)
    implicit none
    real(kind=DEFAULT_PRECISION), dimension(:), allocatable, intent(inout) :: field
    character(len=*), intent(in) :: netcdf_name
    integer,intent(inout) :: dim_length
    integer :: dim_id

    call check_status(nf90_inq_dimid(ncid_dephy, netcdf_name, dim_id))
    call check_status(nf90_inquire_dimension(ncid_dephy, dim_id, len=dim_length))
    allocate(field(dim_length))
    call check_status(nf90_get_var(ncid_dephy, dim_id, field))

  end subroutine dephy_read_dimension_variable


  subroutine dephy_read_profile_variable(field,netcdf_name)
    implicit none

    real(kind=DEFAULT_PRECISION), dimension(:), allocatable, intent(inout) :: field
    real(kind=DEFAULT_PRECISION), dimension(:,:,:,:), allocatable :: field_in_file
    real(kind=DEFAULT_PRECISION), dimension(:), allocatable :: z_out

    character(len=*), intent(in) :: netcdf_name
    integer :: variable_id
    ! note the DEPHY format includes lat/lon dimensions of length 1
    allocate(z_out(kkp))
    allocate(field(kkp))
    allocate(field_in_file(1,1,height_len_dephy,1))

    z_out = module_zn

    call check_status(nf90_inq_varid(ncid_dephy, netcdf_name, variable_id))
    call check_status(nf90_get_var(ncid_dephy, variable_id, field_in_file))
    call piecewise_linear_1d(height_dephy, field_in_file(1,1,:,1), z_out, field)

  end subroutine dephy_read_profile_variable


  subroutine dephy_read_surface_variable(field,netcdf_name)
    implicit none

    real(kind=DEFAULT_PRECISION), dimension(:), allocatable, intent(inout) :: field
    real(kind=DEFAULT_PRECISION), dimension(:,:,:), allocatable :: field_in_file

    character(len=*), intent(in) :: netcdf_name
    integer :: variable_id
    ! note the DEPHY format includes lat/lon dimensions of length 1
    allocate(field(time_len_dephy))
    allocate(field_in_file(1,1,time_len_dephy))

    call check_status(nf90_inq_varid(ncid_dephy, netcdf_name, variable_id))
    call check_status(nf90_get_var(ncid_dephy, variable_id, field_in_file))
    field=field_in_file(1,1,:)

  end subroutine dephy_read_surface_variable


  subroutine dephy_read_forcing_variable(field,netcdf_name)
    implicit none

    real(kind=DEFAULT_PRECISION), dimension(:,:), allocatable, intent(inout) :: field
    real(kind=DEFAULT_PRECISION), dimension(:,:,:,:), allocatable :: field_in_file
    real(kind=DEFAULT_PRECISION), dimension(:), allocatable :: z_out

    character(len=*), intent(in) :: netcdf_name
    integer :: variable_id
    ! note the DEPHY format includes lat/lon dimensions of length 1

    allocate(z_out(kkp))
    allocate(field(kkp,time_len_dephy))
    allocate(field_in_file(1,1,height_len_dephy,time_len_dephy))

    call check_status(nf90_inq_varid(ncid_dephy, netcdf_name, variable_id))
    call check_status(nf90_get_var(ncid_dephy, variable_id, field_in_file))
    if (netcdf_name .eq. 'w') then
      z_out = module_z
    else
      z_out = module_zn
    end if

    call piecewise_linear_2d_k1(height_dephy, time_dephy, field_in_file(1,1,:,:), z_out, field)

  end subroutine dephy_read_forcing_variable


 subroutine dephy_attribute_exists(does_exist, netcdf_name)
    implicit none
    character(len=*), intent(in) :: netcdf_name
    logical, intent(inout) :: does_exist
    integer :: dephy_integer
    integer :: nc_status

    nc_status=nf90_get_att(ncid_dephy, nf90_global, netcdf_name, dephy_integer)
    call check_status(nc_status, does_exist)
  end subroutine dephy_attribute_exists

  subroutine dephy_read_integer(dephy_integer,netcdf_name)
    implicit none
    integer, intent(inout) :: dephy_integer
    character(len=*), intent(in) :: netcdf_name

    call check_status(nf90_get_att(ncid_dephy, nf90_global, netcdf_name, dephy_integer))

  end subroutine dephy_read_integer


  subroutine dephy_read_real(dephy_real,netcdf_name)
    implicit none
    real(kind=DEFAULT_PRECISION):: dephy_real
    double precision:: dephy_double
    character(len=*), intent(in) :: netcdf_name

    call check_status(nf90_get_att(ncid_dephy, nf90_global, netcdf_name, dephy_double))
    dephy_real=1.0_DEFAULT_PRECISION*dephy_double
  end subroutine dephy_read_real


  subroutine dephy_read_string(dephy_string,netcdf_name)
    character(len=STRING_LENGTH), intent(out) :: dephy_string
    character(len=*), intent(in) :: netcdf_name

    call check_status(nf90_get_att(ncid_dephy, nf90_global, netcdf_name, dephy_string))
    call remove_null_terminator_from_string(dephy_string)

  end subroutine dephy_read_string


  subroutine dephy_read_dimension_variables
    implicit none
    call dephy_read_dimension_variable(time_dephy,'time',time_len_dephy)
    call dephy_read_dimension_variable(height_dephy,'lev',height_len_dephy)

  end subroutine dephy_read_dimension_variables


  subroutine dephy_read_profile_variables
    implicit none

    call dephy_read_profile_variable(u_dephy, 'u')
    call dephy_read_profile_variable(v_dephy, 'v')
    call dephy_read_profile_variable(theta_dephy, 'theta')
    call dephy_read_profile_variable(rv_dephy, 'rv')
    call dephy_read_profile_variable(tke_dephy, 'tke')

  end subroutine dephy_read_profile_variables

  subroutine dephy_read_inversion_nudging
    implicit none
    logical :: l_extended_dephy_format

    call dephy_attribute_exists(l_extended_dephy_format,'inversion_nudging')
    if(l_extended_dephy_format) then
       call dephy_read_integer(int_inversion_nudging,'inversion_nudging')
       if(int_inversion_nudging==1) then
           call dephy_read_real(inversion_nudging_height_above,'inversion_nudging_height_above')
           call dephy_read_real(inversion_nudging_transition,'inversion_nudging_transition')
           call dephy_read_real(inversion_nudging_time,'inversion_nudging_time')
       endif
    end if

  end subroutine dephy_read_inversion_nudging


  subroutine dephy_read_forcing_variables
    implicit none

    call dephy_read_forcing_variable(height_forc_dephy,'height_forc')
    call dephy_read_forcing_variable(pressure_forc_dephy, 'pressure_forc')
    call dephy_read_forcing_variable(ug_dephy,'ug')
    call dephy_read_forcing_variable(vg_dephy,'vg')
    call dephy_read_forcing_variable(u_adv_dephy, 'u_adv')
    call dephy_read_forcing_variable(v_adv_dephy, 'v_adv')
    call dephy_read_forcing_variable(theta_adv_dephy, 'theta_adv')
    call dephy_read_forcing_variable(theta_rad_dephy, 'theta_rad')
    call dephy_read_forcing_variable(rv_adv_dephy, 'rv_adv')
    call dephy_read_forcing_variable(w_dephy, 'w')
    call dephy_read_forcing_variable(theta_nudging_dephy, 'theta_nudging')
    call dephy_read_forcing_variable(rv_nudging_dephy, 'rv_nudging')
    call dephy_read_forcing_variable(u_nudging_dephy, 'u_nudging')
    call dephy_read_forcing_variable(v_nudging_dephy, 'v_nudging')
    call dephy_read_forcing_variable(nudging_inv_u_traj_dephy, 'nudging_inv_u_traj')
    call dephy_read_forcing_variable(nudging_inv_v_traj_dephy, 'nudging_inv_v_traj')
    call dephy_read_forcing_variable(nudging_inv_theta_traj_dephy, 'nudging_inv_theta_traj')
    call dephy_read_forcing_variable(nudging_inv_rv_traj_dephy, 'nudging_inv_rv_traj')

  end subroutine dephy_read_forcing_variables


  subroutine dephy_read_surface_variables
    implicit none

    call dephy_read_surface_variable(lat_traj_dephy, 'lat_traj')
    call dephy_read_surface_variable(lon_traj_dephy, 'lon_traj')
    call dephy_read_surface_variable(ps_forc_dephy, 'ps_forc')
    call dephy_read_surface_variable(ts_dephy, 'ts')
    call dephy_read_surface_variable(sfc_sens_flx_dephy, 'sfc_sens_flx')
    call dephy_read_surface_variable(sfc_lat_flx_dephy, 'sfc_lat_flx')
    call dephy_read_surface_variable(z0_traj_dephy, 'z0_traj')
    call dephy_read_surface_variable(z0th_traj_dephy, 'z0th_traj')
    call dephy_read_surface_variable(ustar_dephy, 'ustar')
    call dephy_read_surface_variable(u_traj_dephy, 'u_traj')
    call dephy_read_surface_variable(v_traj_dephy, 'v_traj')
    call dephy_read_surface_variable(albedo_traj_dephy, 'albedo_traj')
    call dephy_read_surface_variable(q_skin_traj_dephy, 'q_skin_traj')

  end subroutine dephy_read_surface_variables


  subroutine dephy_read_integers
    implicit none

    call dephy_read_integer(int_adv_theta, 'adv_theta')
    call dephy_read_integer(int_adv_rv, 'adv_rv')
    call dephy_read_integer(int_rad_theta, 'rad_theta')
    call dephy_read_integer(int_forc_w, 'forc_w')
    call dephy_read_integer(int_forc_geo, 'forc_geo')
    call dephy_read_integer(int_nudging_u, 'nudging_u')
    call dephy_read_integer(int_nudging_v, 'nudging_v')
    call dephy_read_integer(int_nudging_theta, 'nudging_theta')
    call dephy_read_integer(int_nudging_rv, 'nudging_rv')

  end subroutine dephy_read_integers


  subroutine dephy_read_strings
    implicit none

    call dephy_read_string(str_surfaceType, 'surfaceType')
    call dephy_read_string(str_surfaceForcing, 'surfaceForcing')
    call dephy_read_string(str_surfaceForcingWind, 'surfaceForcingWind')

  end subroutine dephy_read_strings


  subroutine dephy_apply_nudging(prof,prof_targ,inv_nudge_time,tendency)
    implicit none
    real(kind=DEFAULT_PRECISION), intent(in) :: prof(:), prof_targ(:),inv_nudge_time(:)
    real(kind=DEFAULT_PRECISION), intent(inout) :: tendency(:,:,:)
    integer :: ii,jj,kk

    do ii=1,size(tendency,3)
    do jj=1,size(tendency,2)
    do kk=1,size(tendency,1)
      tendency(kk,jj,ii)=tendency(kk,jj,ii)-(prof(kk)-prof_targ(kk))*inv_nudge_time(kk)
    end do
    end do
    end do

  end subroutine


  subroutine dephy_add_profile(field_in,prof,field_out)
    implicit none
    real(kind=DEFAULT_PRECISION), intent(in) :: field_in(:,:,:)
    real(kind=DEFAULT_PRECISION), intent(in) :: prof(:)
    real(kind=DEFAULT_PRECISION), intent(out) :: field_out(:,:,:)
    integer :: ii,jj,kk


    do ii=1,size(field_in,3)
    do jj=1,size(field_in,2)
    do kk=1,size(field_in,1)
      field_out(kk,jj,ii)=field_in(kk,jj,ii)+prof(kk)
    end do
    end do
    end do

  end subroutine dephy_add_profile


  subroutine dephy_apply_tendency(tendency_prof,tendency)
    implicit none
    real(kind=DEFAULT_PRECISION), intent(in) :: tendency_prof(:)
    real(kind=DEFAULT_PRECISION), intent(inout) :: tendency(:,:,:)
    integer :: ii,jj,kk


    do ii=1,size(tendency,3)
    do jj=1,size(tendency,2)
    do kk=1,size(tendency,1)
      tendency(kk,jj,ii)=tendency(kk,jj,ii)+tendency_prof(kk)
    end do
    end do
    end do

  end subroutine


  subroutine dephy_set_profile(profile,field)
    implicit none
    real(kind=DEFAULT_PRECISION), intent(in) :: profile(:)
    real(kind=DEFAULT_PRECISION), intent(inout) :: field(:,:,:)
    integer :: ii,jj,kk


    do ii=1,size(field,3)
    do jj=1,size(field,2)
    do kk=1,size(field,1)
      field(kk,jj,ii)=profile(kk)
    end do
    end do
    end do

  end subroutine dephy_set_profile


  subroutine dephy_apply_subsidence(w_prof,field,tendency)
    implicit none
    real(kind=DEFAULT_PRECISION), intent(in) :: w_prof(:)
    real(kind=DEFAULT_PRECISION), intent(in) :: field(:,:,:)
    real(kind=DEFAULT_PRECISION), intent(inout) :: tendency(:,:,:)
    integer :: ii,jj,kk

    do ii=1,size(tendency,3)
    do jj=1,size(tendency,2)
    tendency(1,jj,ii)=tendency(1,jj,ii)-0.5_DEFAULT_PRECISION*w_prof(2)*&
    (field(2,jj,ii)-field(1,jj,ii))/(module_zn(2)-module_zn(1))
    do kk=2,kkp-1
      if(w_prof(kk+1)<0. .and. w_prof(kk)<0.) then
        ! SUBSIDENCE: GET TENDENCY USING LEVEL ABOVE
        tendency(kk,jj,ii)=tendency(kk,jj,ii)-w_prof(kk+1)*&
        (field(kk+1,jj,ii)-field(kk,jj,ii))/(module_zn(kk+1)-module_zn(kk))
      else if(w_prof(kk+1)>0. .and. w_prof(kk+1)>0.) then
        ! UPSIDENCE: GET TENDENCY USING LEVEL BELOW
        tendency(kk,jj,ii)=tendency(kk,jj,ii)-w_prof(kk)*&
        (field(kk,jj,ii)-field(kk-1,jj,ii))/(module_zn(kk)-module_zn(kk-1))
      else
        ! NO CONSISTENT SIGN OF SUBSIDENCE PROFILE, USE MEAN GRADIENTS AND VELOCITIES?
        tendency(kk,jj,ii)=tendency(kk,jj,ii)-0.5_DEFAULT_PRECISION*&
        (w_prof(kk+1)+w_prof(kk))*(field(kk+1,jj,ii)-field(kk-1,jj,ii))/(module_zn(kk+1)-module_zn(kk-1))
      end if
    end do
    tendency(kkp,jj,ii)=tendency(kkp,jj,ii)-w_prof(kkp)*&
    (field(kkp,jj,ii)-field(kkp-1,jj,ii))/(module_zn(kkp)-module_zn(kkp-1))
    end do
    end do

  end subroutine

  ! Implements energy consistent (non-traditional) coriolis force using time-dependent geostrophic wind
  subroutine dephy_coriolis(u,v,w,u_geo,v_geo,u_gal,v_gal,lat,su,sv,sw)
    implicit none
    real(kind=DEFAULT_PRECISION), intent(in) :: u(:,:,:),v(:,:,:),w(:,:,:)
    real(kind=DEFAULT_PRECISION), intent(in) :: u_geo(:),v_geo(:)
    real(kind=DEFAULT_PRECISION), intent(in) :: u_gal,v_gal
    real(kind=DEFAULT_PRECISION), intent(in) :: lat
    real(kind=DEFAULT_PRECISION), intent(inout) :: su(:,:,:),sv(:,:,:),sw(:,:,:)
    real(kind=DEFAULT_PRECISION) :: fcoriol, fcoriol2
    real(kind=DEFAULT_PRECISION), parameter :: omega_earth=7.2921e-5 ! radial frecuency of earth's rotation
    integer ii,jj,kk

    fcoriol=2.0_DEFAULT_PRECISION*omega_earth*sin(lat*proper_pi/180.0_DEFAULT_PRECISION)
    ! Non-traditional coriolis terms, needed for energy-consistency
    ! See e.g. Igel and Biello 2020
    ! Note geostrophic wind parametrises pressure gradients, and is therefore not in the non-traditional terms.
    fcoriol2=2.0_DEFAULT_PRECISION*omega_earth*cos(lat*proper_pi/180.0_DEFAULT_PRECISION)

#if defined(U_ACTIVE) && defined(V_ACTIVE)
    do ii=2,size(su,3)-1
    do jj=2,size(su,2)-1
    do kk=2,kkp
      su(kk, jj, ii)=su(kk, jj, ii)+fcoriol*&
           (0.25_DEFAULT_PRECISION*(v(kk, jj, ii)+v(kk, jj, ii+1)+&
           v(kk, jj-1, ii)+v(kk, jj-1, ii+1))+&
           v_gal-v_geo(kk))-&
           fcoriol2*&
           (0.25_DEFAULT_PRECISION*(w(kk, jj, ii)+w(kk, jj, ii+1)+&
           w(kk-1, jj, ii)+w(kk-1, jj, ii+1)))

      sv(kk, jj, ii)=sv(kk, jj, ii)-fcoriol*&
           (0.25_DEFAULT_PRECISION*(u(kk, jj, ii)+u(kk, jj, ii-1)+&
           u(kk, jj+1, ii)+u(kk, jj+1, ii-1))+&
           u_gal-u_geo(kk))

    end do
    do kk=2,kkp-1
      sw(kk, jj, ii)=sw(kk,jj,ii)+fcoriol2*&
           (0.25_DEFAULT_PRECISION*(u(kk, jj, ii)+u(kk+1, jj, ii)+&
           u(kk, jj, ii-1)+u(kk+1, jj, ii-1))+u_gal)
    end do
    end do
    end do
#endif

  end subroutine dephy_coriolis


  subroutine dephy_time_interpolate(current_state)
    implicit none
    type(model_state_type), target, intent(inout) :: current_state

    ! interpolate surface variables
    call interpolate_point_linear_1d(time_dephy, lat_traj_dephy, current_state%time, lat_traj)
    call interpolate_point_linear_1d(time_dephy, lon_traj_dephy, current_state%time, lon_traj)
    call interpolate_point_linear_1d(time_dephy, ps_forc_dephy, current_state%time, ps_forc)
    call interpolate_point_linear_1d(time_dephy, ts_dephy, current_state%time, ts)
    call interpolate_point_linear_1d(time_dephy, sfc_sens_flx_dephy, current_state%time, sfc_sens_flx)
    call interpolate_point_linear_1d(time_dephy, sfc_lat_flx_dephy, current_state%time, sfc_lat_flx)
    call interpolate_point_linear_1d(time_dephy, z0_traj_dephy, current_state%time, z0_traj)
    call interpolate_point_linear_1d(time_dephy, z0th_traj_dephy, current_state%time, z0th_traj)
    call interpolate_point_linear_1d(time_dephy, ustar_dephy, current_state%time, ustar)
    call interpolate_point_linear_1d(time_dephy, u_traj_dephy, current_state%time, u_traj)
    call interpolate_point_linear_1d(time_dephy, v_traj_dephy, current_state%time, v_traj)
    call interpolate_point_linear_1d(time_dephy, albedo_traj_dephy, current_state%time, albedo_traj)
    call interpolate_point_linear_1d(time_dephy, q_skin_traj_dephy, current_state%time, q_skin_traj)

    ! interpolate height-dependent variables
    call interpolate_point_linear_2d(time_dephy, height_forc_dephy, current_state%time, height_forc)
    call interpolate_point_linear_2d(time_dephy, pressure_forc_dephy, current_state%time, pressure_forc)
    call interpolate_point_linear_2d(time_dephy, ug_dephy, current_state%time, ug)
    call interpolate_point_linear_2d(time_dephy, vg_dephy, current_state%time, vg)
    call interpolate_point_linear_2d(time_dephy, u_adv_dephy, current_state%time, u_adv)
    call interpolate_point_linear_2d(time_dephy, v_adv_dephy, current_state%time, v_adv)
    call interpolate_point_linear_2d(time_dephy, theta_adv_dephy, current_state%time, theta_adv)
    call interpolate_point_linear_2d(time_dephy, theta_rad_dephy, current_state%time, theta_rad)
    call interpolate_point_linear_2d(time_dephy, rv_adv_dephy, current_state%time, rv_adv)
    call interpolate_point_linear_2d(time_dephy, w_dephy, current_state%time, w)
    call interpolate_point_linear_2d(time_dephy, theta_nudging_dephy, current_state%time, theta_nudging)
    call interpolate_point_linear_2d(time_dephy, rv_nudging_dephy, current_state%time, rv_nudging)
    call interpolate_point_linear_2d(time_dephy, u_nudging_dephy, current_state%time, u_nudging)
    call interpolate_point_linear_2d(time_dephy, v_nudging_dephy, current_state%time, v_nudging)
    if(int_inversion_nudging==0) then
      call interpolate_point_linear_2d(time_dephy, nudging_inv_u_traj_dephy, current_state%time, nudging_inv_u_traj)
      call interpolate_point_linear_2d(time_dephy, nudging_inv_v_traj_dephy, current_state%time, nudging_inv_v_traj)
      call interpolate_point_linear_2d(time_dephy, nudging_inv_theta_traj_dephy, current_state%time, nudging_inv_theta_traj)
      call interpolate_point_linear_2d(time_dephy, nudging_inv_rv_traj_dephy, current_state%time, nudging_inv_rv_traj)
    endif
  end subroutine dephy_time_interpolate

  real(kind=DEFAULT_PRECISION) function cos_transition(absolute_input, transition_start, transition_end)
    real(kind=DEFAULT_PRECISION), intent(in) :: absolute_input
    real(kind=DEFAULT_PRECISION), intent(in) :: transition_start
    real(kind=DEFAULT_PRECISION), intent(in) :: transition_end
    real(kind=DEFAULT_PRECISION) :: normalised_input

    ! function that smoothly transitions from 1 to 0 using a
    ! cosine-shaped transition between start and end
    ! start can be larger than end, in which case is applies in reverse order
    normalised_input = (absolute_input-transition_start)/(transition_end-transition_start)
    if(normalised_input<0.0_DEFAULT_PRECISION) then
      cos_transition=1.0_DEFAULT_PRECISION
    elseif(normalised_input>1.0_DEFAULT_PRECISION) then
      cos_transition=0.0_DEFAULT_PRECISION
    else
      cos_transition=0.5_DEFAULT_PRECISION+0.5_DEFAULT_PRECISION*cos(normalised_input*proper_pi)
    end if
  end function cos_transition


  subroutine dephy_calc_interactive_nudging_profiles(current_state)
    implicit none
    type(model_state_type), target, intent(inout) :: current_state
    real(kind=DEFAULT_PRECISION) :: theta_l_max_grad
    real(kind=DEFAULT_PRECISION) :: z_inversion
    real(kind=DEFAULT_PRECISION) :: this_weight
    integer :: ii,iql,jj,kk,kk_inversion_plus
    iql=get_q_index(standard_q_names%CLOUD_LIQUID_MASS, 'dephy_forcings')

    ! calculate theta_l
    ! average theta_l
    do ii=1,size(theta_l,3)
    do jj=1,size(theta_l,2)
    do kk=1,size(theta_l,1)
      theta_l(kk,jj,ii)=full_theta(kk,jj,ii)-&
      current_state%sq(iql)%data(kk,jj,ii)*rlvap_over_cp/current_state%global_grid%configuration%vertical%rprefrcp(kk)
    end do
    end do
    end do
    call calculate_theta_l_mean(current_state)
    theta_l_max_grad=0.0_DEFAULT_PRECISION
    do kk=2,size(theta_l,1)
      if(module_zn(kk)<6000.0_DEFAULT_PRECISION) then ! clip at 6 km, currently hardcoded
        if((theta_l_mean(kk)-theta_l_mean(kk-1))/(module_zn(kk)-module_zn(kk-1))>theta_l_max_grad) then
          theta_l_max_grad=(theta_l_mean(kk)-theta_l_mean(kk-1))/(module_zn(kk)-module_zn(kk-1))
          kk_inversion_plus=kk
         end if
      else if(kk_inversion_plus==0) then
        kk_inversion_plus=kk
      end if
    end do
    z_inversion=0.5_DEFAULT_PRECISION*(module_zn(kk_inversion_plus)+module_zn(kk_inversion_plus-1))
    do kk=1,size(theta_l,1)
       this_weight=cos_transition(module_zn(kk),z_inversion+inversion_nudging_height_above+inversion_nudging_transition,&
       z_inversion+inversion_nudging_height_above)
       nudging_inv_u_traj(kk)=this_weight/inversion_nudging_time
       nudging_inv_v_traj(kk)=this_weight/inversion_nudging_time
       nudging_inv_theta_traj(kk)=this_weight/inversion_nudging_time
       nudging_inv_rv_traj(kk)=this_weight/inversion_nudging_time
    end do

  end subroutine dephy_calc_interactive_nudging_profiles

  subroutine dephy_apply_forcings(current_state)
    implicit none
    type(model_state_type), target, intent(inout) :: current_state
    ! Surface fields, in time
    integer:: iqv,nn

    iqv=get_q_index(standard_q_names%VAPOUR, 'dephy_forcings')

    ! CALCULATE COMPLETE THETA AND THETA_L
    call dephy_add_profile(current_state%th%data,current_state%global_grid%configuration%vertical%thref,full_theta)
    if(int_inversion_nudging==1) then
      call dephy_calc_interactive_nudging_profiles(current_state)
    endif

    ! APPLY FORCINGS
    ! NUDGINGS NEED MEAN
#ifdef U_ACTIVE
    if(int_nudging_u==1) then
      call dephy_apply_nudging(current_state%global_grid%configuration%vertical%olubar,&
      u_nudging,nudging_inv_u_traj,current_state%su%data)
    end if
#endif
#ifdef V_ACTIVE
    if(int_nudging_v==1) then
      call dephy_apply_nudging(current_state%global_grid%configuration%vertical%olvbar,&
      v_nudging,nudging_inv_v_traj,current_state%sv%data)
    end if
#endif
    if(int_nudging_theta==1) then
      call dephy_apply_nudging(current_state%global_grid%configuration%vertical%olthbar+&
      current_state%global_grid%configuration%vertical%thref,&
      theta_nudging,nudging_inv_theta_traj,current_state%sth%data)
    end if
    if(int_nudging_rv==1) then
      call dephy_apply_nudging(current_state%global_grid%configuration%vertical%olqbar(:,iqv),&
      rv_nudging,nudging_inv_rv_traj,current_state%sq(iqv)%data)
    end if

    ! PROFILES
#ifdef U_ACTIVE
    call dephy_apply_tendency(u_adv,current_state%su%data)
#endif
#ifdef V_ACTIVE
    call dephy_apply_tendency(v_adv,current_state%sv%data)
#endif
    if(int_adv_theta==1) then
      call dephy_apply_tendency(theta_adv,current_state%sth%data)
    end if
    if(int_adv_rv==1) then
      call dephy_apply_tendency(rv_adv,current_state%sq(iqv)%data)
    end if

    ! RADIATION TENDENCIES
    if(int_rad_theta==0) then
      call dephy_update_socrates(socrates_opt,lat_traj,lon_traj,albedo_traj)
    elseif(int_rad_theta==1) then
      call dephy_apply_tendency(theta_rad,current_state%sth%data)
    end if

    ! LARGE-SCALE VERTICAL WIND
    ! USE A DOWNWIND FORMULATION (AS IN DALES),
    ! SO ADVECTION ONLY IN DIRECTION OF WIND
    ! ALWAYS USE LOCAL GRADIENTS, AS NON-LOCAL ONES ARE UNPHYSICAL
    ! APPLY TO ALL Q SPECIES
    if(int_forc_w ==1) then
#ifdef U_ACTIVE
      call dephy_apply_subsidence(w,current_state%u%data,current_state%su%data)
#endif
#ifdef V_ACTIVE
      call dephy_apply_subsidence(w,current_state%v%data,current_state%sv%data)
#endif
    end if
    if((int_forc_w ==1) .OR. (int_forc_w==2)) then
      call dephy_apply_subsidence(w,full_theta,current_state%sth%data)
      DO nn=1,current_state%number_q_fields
        call dephy_apply_subsidence(w,current_state%q(nn)%data,current_state%sq(nn)%data)
      END DO
    end if

    ! IMPLEMENTATION OF FULL CORIOLIS FORCE
    if(int_forc_geo==1) then
      call dephy_coriolis(current_state%u%data,current_state%v%data,current_state%w%data,&
      ug,vg,current_state%ugal,current_state%vgal,lat_traj,&
      current_state%su%data,current_state%sv%data,current_state%sw%data)
    end if

  end subroutine dephy_apply_forcings


  subroutine dephy_initial_profiles(current_state)
    implicit none
    type(model_state_type), target, intent(inout) :: current_state
    integer :: iqv
    logical :: l_matchthref

    iqv=get_q_index(standard_q_names%VAPOUR, 'dephy_forcings')
    l_matchthref=options_get_logical(current_state%options_database, "l_matchthref")
    if(l_matchthref) then
        if(.not. current_state%use_anelastic_equations) then
           call log_master_log(LOG_ERROR, "Non-anelastic equation set and l_maththref are incompatible")
         end if
      current_state%global_grid%configuration%vertical%thref(:)=theta_dephy
    else
      current_state%global_grid%configuration%vertical%thref(:)=current_state%thref0
    endif

    call dephy_set_profile(u_dephy,current_state%u%data)
    call dephy_set_profile(u_dephy,current_state%zu%data)
    call dephy_set_profile(v_dephy,current_state%v%data)
    call dephy_set_profile(v_dephy,current_state%zv%data)
    call dephy_set_profile(theta_dephy-current_state%global_grid%configuration%vertical%thref,current_state%th%data)
    call dephy_set_profile(theta_dephy-current_state%global_grid%configuration%vertical%thref,current_state%zth%data)
    ! Note q in MONC is mixing ratio (rather than specific humidity, as is more usual)
    call dephy_set_profile(rv_dephy,current_state%q(iqv)%data)
    call dephy_set_profile(rv_dephy,current_state%zq(iqv)%data)

  end subroutine dephy_initial_profiles

  !!! THIS SHOULD LIKELY COME FROM A UTILITIES MODULE
  !> Will check a NetCDF status and write to log_log error any decoded statuses. Can be used to decode
  !! whether a dimension or variable exists within the NetCDF data file
  !! @param status The NetCDF status flag
  !! @param foundFlag Whether the field has been found or not
  subroutine check_status(status, found_flag)
    integer, intent(in) :: status
    logical, intent(out), optional :: found_flag

    if (present(found_flag)) then
      found_flag = status /= nf90_ebaddim .and. status /= nf90_enotatt .and. status /= nf90_enotvar
      if (.not. found_flag) return
    end if

    if (status /= nf90_noerr) then
      call log_log(LOG_ERROR, "NetCDF returned error code of "//trim(nf90_strerror(status)))
    end if
  end subroutine check_status

  !!! THIS SUBROUTINE REPLACES set_vertical_reference_profile IN gridmanager.F90
  !> Sets up the vertical grid reference profile at each point
  !! @param current_state The current model state_mod
  !! @param vertical_grid The vertical grid that we are working on
  !! @param kkp Number of grid points in a vertical column
  subroutine dephy_profiles_etc(current_state, vertical_grid)
    implicit none
    type(model_state_type), intent(inout) :: current_state
    type(vertical_grid_configuration_type), intent(inout) :: vertical_grid

    integer :: k

    call dephy_initial_profiles(current_state)
    call set_up_vertical_reference_properties(current_state, vertical_grid, current_state%global_grid%size(Z_INDEX))

    if(l_verbose) write(*,*) "initialised dephy 6.3"

    call set_anelastic_pressure(current_state)

    ! THIS CRUCIAL STATEMENT IS HIDDEN IN set_qv_init_from_rh in gridmanager.F90
    vertical_grid=current_state%global_grid%configuration%vertical

    if(l_verbose) write(*,*) "initialised dephy 6.4"

    do k=2,kkp-1
      ! for diffusion onto p-level from below
      vertical_grid%czb(k)=(vertical_grid%rho(k-1)/vertical_grid%rhon(k))/(vertical_grid%dz(k)*vertical_grid%dzn(k))
      ! for diffusion onto p-level from above
      vertical_grid%cza(k)=(vertical_grid%rho(k)/vertical_grid%rhon(k))/(vertical_grid%dz(k)*vertical_grid%dzn(k+1))
      vertical_grid%czg(k)=-vertical_grid%czb(k)-vertical_grid%cza(k)
      if (k .gt. 2) vertical_grid%czh(k)=vertical_grid%czb(k)*vertical_grid%cza(k-1)
    end do
    do k=2,kkp-1
      ! advection onto p-level from below
      vertical_grid%tzc1(k)=0.25_DEFAULT_PRECISION*vertical_grid%rdz(k)*vertical_grid%rho(k-1)/vertical_grid%rhon(k)
      ! advection onto p-level from above
      vertical_grid%tzc2(k)=0.25_DEFAULT_PRECISION*vertical_grid%rdz(k)*vertical_grid%rho(k)/vertical_grid%rhon(k)
    end do
    do k=2,kkp-1
      ! advection onto w-level (K) from below
      vertical_grid%tzd1(k)=0.25_DEFAULT_PRECISION*vertical_grid%rdzn(k+1)*vertical_grid%rhon(k)/vertical_grid%rho(k)
      ! advection onto w-level (K) from above
      vertical_grid%tzd2(k)=0.25_DEFAULT_PRECISION*vertical_grid%rdzn(k+1)*vertical_grid%rhon(k+1)/vertical_grid%rho(k)
    end do
    k=kkp
    vertical_grid%czb(k)=(vertical_grid%rho(k-1)/vertical_grid%rhon(k))/(vertical_grid%dz(k)*vertical_grid%dzn(k))
    vertical_grid%cza(k)=0.0_DEFAULT_PRECISION
    vertical_grid%czg(k)=-vertical_grid%czb(k)
    vertical_grid%czh(k)=vertical_grid%czb(k)*vertical_grid%cza(k-1)
    vertical_grid%tzc2(k)=0.25_DEFAULT_PRECISION*vertical_grid%rdz(k)*vertical_grid%rho(k)/vertical_grid%rhon(k)
    vertical_grid%tzc1(k)=0.25_DEFAULT_PRECISION*vertical_grid%rdz(k)*vertical_grid%rho(k-1)/vertical_grid%rhon(k)
    vertical_grid%czn=vertical_grid%dzn(2)*0.5_DEFAULT_PRECISION
    vertical_grid%zlogm=log(1.0_DEFAULT_PRECISION+vertical_grid%zn(2)/z0)
    vertical_grid%zlogth=log((vertical_grid%zn(2)+z0)/z0th)
    vertical_grid%vk_on_zlogm=von_karman_constant/vertical_grid%zlogm

    if(l_verbose) write(*,*) "initialised dephy 6.5"

    call setup_reference_state_liquid_water_temperature_and_saturation(&
         current_state, vertical_grid, current_state%global_grid%size(Z_INDEX))

    if(l_verbose) write(*,*) "initialised dephy 6.6"

    call calculate_mixing_length_for_neutral_case(current_state, vertical_grid, current_state%global_grid%size(Z_INDEX))

    if(l_verbose) write(*,*) "initialised dephy 6.7"

    call set_buoyancy_coefficient(current_state, vertical_grid, current_state%global_grid%size(Z_INDEX))

  end subroutine dephy_profiles_etc

  !!! THIS SUBROUTINE REPLACES the initialisation_callback IN setfluxlook.F90
  subroutine dephy_setfluxlook_init(current_state)
    type(model_state_type), intent(inout), target :: current_state

    integer, parameter :: LOOKUP_ENTRIES = 80    !< Number of entries for MO lookup tables
    integer :: iqv

    current_state%lookup_table_entries=LOOKUP_ENTRIES
    current_state%saturated_surface = .true. ! Copied from setfluxlook module
                                             ! We will change this if we find some humidity data
    allocate(current_state%lookup_table_velocity(current_state%lookup_table_entries), &
         current_state%lookup_table_ustr(current_state%lookup_table_entries))

    if(l_verbose) write(*,*) "initialised dephy 7.1"

    if (.not. allocated(current_state%cq))then
     allocate(current_state%cq(current_state%number_q_fields))
     current_state%cq=0.0_DEFAULT_PRECISION
    end if

    iqv = get_q_index(standard_q_names%VAPOUR, 'dephy_forcings')
    current_state%cq(iqv) = ratio_mol_wts-1.0

    if(l_verbose) write(*,*) "initialised dephy 7.2"

    if (trim(str_surfaceForcing) == "ts") then
       current_state%type_of_surface_boundary_conditions = PRESCRIBED_SURFACE_VALUES
       current_state%use_surface_boundary_conditions = .true.
    else if (trim(str_surfaceForcing) == "surfaceFlux") then
       current_state%type_of_surface_boundary_conditions = PRESCRIBED_SURFACE_FLUXES
       current_state%use_surface_boundary_conditions = .true.
    else
       call log_master_log(LOG_ERROR, "Surface condition for dephy not implemented")
    endif

    if(l_verbose) write(*,*) "initialised dephy 7.3"

    if (current_state%type_of_surface_boundary_conditions == PRESCRIBED_SURFACE_FLUXES) then
       call dephy_set_flux(current_state)

       current_state%fbuoy=0.
       current_state%fbuoy=&
       current_state%global_grid%configuration%vertical%buoy_co(1)*current_state%surface_temperature_flux+&
                     current_state%cq(iqv)*current_state%surface_vapour_flux*G
       call set_look(current_state)
       current_state%theta_surf=0.0_DEFAULT_PRECISION
       current_state%surface_vapour_mixing_ratio=0.0_DEFAULT_PRECISION
    else if (current_state%type_of_surface_boundary_conditions == PRESCRIBED_SURFACE_VALUES) then   ! Prescribed surface temperatures
       call dephy_set_flux(current_state)
    end if

    if(l_verbose) write(*,*) "initialised dephy 7.4"

  end subroutine dephy_setfluxlook_init

  !!! THIS SUBROUTINE REPLACES the timestep_callback IN setfluxlook.F90
  subroutine dephy_setfluxlook_timestep(current_state)
    type(model_state_type), intent(inout), target :: current_state

    if(l_verbose) write(*,*) "dephy timestep 1.1"
    call dephy_set_flux(current_state)
    if(l_verbose) write(*,*) "dephy timestep 1.2"
    call change_look(current_state)
    if(l_verbose) write(*,*) "dephy timestep 1.3"

  end subroutine dephy_setfluxlook_timestep

  !!! THIS SUBROUTINE REPLACES the set_flux subroutine IN setfluxlook.F90
  subroutine dephy_set_flux(current_state)
    type(model_state_type), intent(inout), target :: current_state

    integer :: iqv
    iqv = get_q_index(standard_q_names%VAPOUR, 'dephy_forcings')

    if(l_verbose) write(*,*) "initialised dephy 7.3.1"

    if (current_state%type_of_surface_boundary_conditions == PRESCRIBED_SURFACE_FLUXES) then  ! Prescribed surface fluxes

      if(l_verbose) write(*,*) "initialised dephy 7.3.1.0"

      current_state%surface_temperature_flux=sfc_sens_flx/(current_state%global_grid%configuration%vertical%rho(1)*cp)
      current_state%surface_vapour_flux=sfc_lat_flx/(current_state%global_grid%configuration%vertical%rho(1)*rlvap)

      ! Update buoyancy flux...
      current_state%fbuoynew=0.0_DEFAULT_PRECISION
      current_state%fbuoynew=&
         current_state%global_grid%configuration%vertical%buoy_co(1)*current_state%surface_temperature_flux
      current_state%fbuoynew=current_state%fbuoynew+current_state%cq(iqv)*current_state%surface_vapour_flux*G

    else if (current_state%type_of_surface_boundary_conditions == PRESCRIBED_SURFACE_VALUES) then   ! Prescribed surface temperatures

      if(l_verbose) write(*,*) "initialised dephy 7.3.1.1"

      if (current_state%saturated_surface)then
        current_state%surface_vapour_mixing_ratio = qsaturation(ts,current_state%surface_pressure*0.01)
      else
        call log_master_log(LOG_ERROR, "DEPHY: prescribed surface vapout mixing ratio not implemented"//&
        " (note q_skin is reservoir content!)")
      end if

      if(l_verbose) write(*,*) "initialised dephy 7.3.1.2"

      ! Set theta_v
      current_state%theta_surf = ts*&
         (current_state%surface_reference_pressure/current_state%surface_pressure)**r_over_cp
      current_state%theta_virtual_surf = current_state%theta_surf
      current_state%theta_virtual_surf = current_state%theta_surf +  &
             current_state%global_grid%configuration%vertical%thref(2)*  &
             current_state%cq(iqv)*current_state%surface_vapour_mixing_ratio

      if(l_verbose) write(*,*) "initialised dephy 7.3.1.3"

      ! Finally set up new values of THVSURF dependent constants
      current_state%cmbc=betam*current_state%global_grid%configuration%vertical%zn(2)*G*&
         von_karman_constant/current_state%theta_virtual_surf

      if(l_verbose) write(*,*) "initialised dephy 7.3.1.4"

      current_state%rcmbc=1.0_DEFAULT_PRECISION/current_state%cmbc
      current_state%ellmocon=current_state%theta_virtual_surf/(G*von_karman_constant)

      if(l_verbose) write(*,*) "initialised dephy 7.3.1.5"

    end if

  end subroutine dephy_set_flux

  ! lowerbc need to re-initialise due to evolving z0/z0th
  ! It is possible to improve on this by making more changes to the lowerbc code
  subroutine lowerbc_reset_constants(current_state)
    type(model_state_type), target, intent(inout) :: current_state

    real(kind=DEFAULT_PRECISION) :: bhbc

    if ( current_state%use_surface_boundary_conditions .and.  &
         current_state%type_of_surface_boundary_conditions == PRESCRIBED_SURFACE_VALUES) then
      ! variables below are only required when PRESCRIBED_SURFACE_VALUES are used.
       tstrcona=von_karman_constant/alphah*current_state%global_grid%configuration%vertical%zlogth
       bhbc=alphah*current_state%global_grid%configuration%vertical%zlogth
       rhmbc=betah*(current_state%global_grid%configuration%vertical%zn(2)+z0-z0th)/&
            (betam*current_state%global_grid%configuration%vertical%zn(2))
       ddbc=current_state%global_grid%configuration%vertical%zlogm*(bhbc-&
            rhmbc*current_state%global_grid%configuration%vertical%zlogm)
       ddbc_x4=4.*ddbc
       r2ddbc=0.5_DEFAULT_PRECISION/ddbc
       eecon=2.0_DEFAULT_PRECISION*rhmbc*current_state%global_grid%configuration%vertical%zlogm-bhbc
       rcmbc=1.0_DEFAULT_PRECISION/current_state%cmbc
       tstrconb=von_karman_constant/alphah
       x4con=gammam*(current_state%global_grid%configuration%vertical%zn(2)+z0)
       xx0con=gammam*z0
       y2con=gammah*(current_state%global_grid%configuration%vertical%zn(2)+z0)
       yy0con=gammah*z0th
    endif

  end subroutine lowerbc_reset_constants


!!! COPY OF PIECEWISE LINEAR INTERPOLATION ROUTINE, BUT INCLUDING THE K=1 LEVEL
!!! MAINLY TO BE SAFE

  !> Does a simple 1d linear interpolation to a point
  !! @param zvals input z nodes
  !! @param vals  input nodal values
  !! @param z location to interpolate onto
  !! @param f output interpolated value
  subroutine piecewise_linear_2d_k1(zvals, time_vals, vals, z, field)

    ! Assumes input variables (vals) are 2-D, with dims (z, time)

    real(kind=DEFAULT_PRECISION), intent(in) :: zvals(:), time_vals(:)
    real(kind=DEFAULT_PRECISION), intent(in) :: vals(:,:)
    real(kind=DEFAULT_PRECISION), intent(in) :: z(:)
    real(kind=DEFAULT_PRECISION), intent(out) :: field(:,:)

    real(kind=DEFAULT_PRECISION) :: scale_tmp

    integer :: nn, k_monc, k_force                     ! loop counter
    integer :: nz_force, nt_force, nz_monc, nt_monc    ! time and height array sizes for forcing and monc grids

    nz_force = size(zvals)
    nt_force = size(time_vals)
    nz_monc  = size(z)
    nt_monc  = size(time_vals) ! time is intepolated in the timestep callback

    if ( zvals(1) .GT. zvals(nz_force) ) then   ! pressure
       call log_master_log(LOG_ERROR, "Input forcing uses pressure, this has not been coded"// &
            " - please modify your forcing file to using height coordinates or modify the" // &
            " interpolation routine in model_core to work with pressure coords - STOP")
    else
       do k_monc=1,nz_monc
          do k_force=1,nz_force-1
             if( z(k_monc) >= zvals(k_force) .AND. z(k_monc) < zvals(k_force+1) ) then
                scale_tmp = ( z(k_monc) - zvals(k_force) ) /              &
                     ( zvals(k_force+1) - zvals(k_force) )
                do nn=1, nt_force
                   field(k_monc,nn) = vals(k_force,nn) +                  &
                        (  vals(k_force+1,nn) - vals(k_force,nn) )        &
                        * scale_tmp
                enddo
             endif
          enddo
       enddo
       ! now examine the cases below and above forlevs(1) and forlevs(ktmfor
       ! uses the local vertical gradient in the forcing to determine the
       ! new values
       do k_monc=1,nz_monc
          if ( z(k_monc) >= zvals(nz_force) ) then
             scale_tmp = ( z(k_monc) - zvals(nz_force) )                   &
                  / ( zvals(nz_force) - zvals(nz_force-1) )
             do nn=1,nt_force
                field(k_monc,nn) = vals(nz_force,nn) +                  &
                     (  vals(nz_force,nn) - vals(nz_force-1,nn) )        &
                     * scale_tmp
             enddo
          elseif ( z(k_monc) < zvals(1) )THEN
             scale_tmp = ( z(k_monc) - zvals(1) )                        &
                  / ( zvals(1) - zvals(2) )
             do nn=1,nt_force
                field(k_monc,nn) = vals(1,nn) +                  &
                     (  vals(1,nn) - vals(2,nn) )        &
                     * scale_tmp
             enddo
          endif
       enddo
       !
    endif   ! pressure or height

  end subroutine piecewise_linear_2d_k1


integer function maxloc1(field)
   real(kind=DEFAULT_PRECISION), intent(in), dimension(:,:,:) :: field
   integer, dimension(3) :: maxloc_res
   maxloc_res=maxloc(field)
   maxloc1=maxloc_res(1)
end function maxloc1


integer function minloc1(field)
   real(kind=DEFAULT_PRECISION), intent(in), dimension(:,:,:) :: field
   integer, dimension(3) :: minloc_res
   minloc_res=minloc(field)
   minloc1=minloc_res(1)
end function minloc1


subroutine dephy_bughunting(current_state)
   implicit none
   type(model_state_type), intent(in) :: current_state
   integer:: iqv

   iqv=get_q_index(standard_q_names%VAPOUR, 'dephy_forcings')

   ! tendency debugging
   write (*,'(A)') 'DEPHY MANUAL DEBUGGING ROUTINE'
   write (*,*) 'time ',current_state%time

   write (*,'(A)') '                  su          sv          sw         sth         sqv'
   write (*,'(A,5ES12.2)') 'max vals',maxval(current_state%su%data), &
   maxval(current_state%sv%data), maxval(current_state%sw%data), &
   maxval(current_state%sth%data), maxval(current_state%sq(iqv)%data)
   write (*,'(A,5I12)') 'max loc',maxloc1(current_state%su%data), &
   maxloc1(current_state%sv%data), maxloc1(current_state%sw%data), &
   maxloc1(current_state%sth%data), maxloc1(current_state%sq(iqv)%data)
   write (*,'(A)') ' '
   write (*,'(A,5ES12.2)') 'min zvals',minval(current_state%su%data), &
   minval(current_state%sv%data), minval(current_state%sw%data), &
   minval(current_state%sth%data), minval(current_state%sq(iqv)%data)
   write (*,'(A,5I12)') 'min zloc',minloc1(current_state%su%data), &
   minloc1(current_state%sv%data), minloc1(current_state%sw%data), &
   minloc1(current_state%sth%data), minloc1(current_state%sq(iqv)%data)
   write (*,'(A)') ' '
   write (*,'(A)') ' '

   ! value debugging
   write (*,'(A)') '                   u           v           w          th          qv'
   write (*,'(A,5ES12.2)') 'max vals',maxval(current_state%u%data), &
   maxval(current_state%v%data), maxval(current_state%w%data), &
   maxval(current_state%th%data), maxval(current_state%q(iqv)%data)
   write (*,'(A,5I12)') 'max zloc',maxloc1(current_state%u%data), &
   maxloc1(current_state%v%data), maxloc1(current_state%w%data), &
   maxloc1(current_state%th%data), maxloc1(current_state%q(iqv)%data)
   write (*,'(A)') ' '
   write (*,'(A,5ES12.2)') 'min vals',minval(current_state%u%data), &
   minval(current_state%v%data), minval(current_state%w%data), &
   minval(current_state%th%data), minval(current_state%q(iqv)%data)
   write (*,'(A,5I12)') 'min zloc',minloc1(current_state%u%data), &
   minloc1(current_state%v%data), minloc1(current_state%w%data), &
   minloc1(current_state%th%data), minloc1(current_state%q(iqv)%data)
   write (*,'(A)') ' '
   write (*,'(A)') ' '

   !z value debugging
   write (*,'(A)') '                  zu          zv          zw         zth         zqv'
   write (*,'(A,5ES12.2)') 'max vals',maxval(current_state%zu%data), &
   maxval(current_state%zv%data), maxval(current_state%zw%data), &
   maxval(current_state%zth%data), maxval(current_state%zq(iqv)%data)
   write (*,'(A,5I12)') 'max zloc',maxloc1(current_state%zu%data), &
   maxloc1(current_state%zv%data), maxloc1(current_state%zw%data), &
   maxloc1(current_state%zth%data), maxloc1(current_state%zq(iqv)%data)
   write (*,'(A)') ' '
   write (*,'(A,5ES12.2)') 'min vals',minval(current_state%zu%data), &
   minval(current_state%zv%data), minval(current_state%zw%data), &
   minval(current_state%zth%data), minval(current_state%zq(iqv)%data)
   write (*,'(A,5I12)') 'min zloc',minloc1(current_state%zu%data), &
   minloc1(current_state%zv%data), minloc1(current_state%zw%data), &
   minloc1(current_state%zth%data), minloc1(current_state%zq(iqv)%data)

end subroutine dephy_bughunting

subroutine calculate_theta_l_mean(current_state)
  type(model_state_type), intent(inout) :: current_state

  integer :: k, ierr
  real(kind=DEFAULT_PRECISION) :: rnhpts

   rnhpts=1.0_DEFAULT_PRECISION/real(current_state%global_grid%size(X_INDEX)*current_state%global_grid%size(Y_INDEX))

   do k=current_state%local_grid%local_domain_start_index(Z_INDEX), current_state%local_grid%local_domain_end_index(Z_INDEX)
      theta_l_mean(k)=sum(theta_l(k, &
      current_state%local_grid%local_domain_start_index(Y_INDEX):current_state%local_grid%local_domain_end_index(Y_INDEX), &
      current_state%local_grid%local_domain_start_index(X_INDEX):current_state%local_grid%local_domain_end_index(X_INDEX)  &
      ))
    end do

  call mpi_allreduce(MPI_IN_PLACE, theta_l_mean, current_state%local_grid%size(Z_INDEX), PRECISION_TYPE, MPI_SUM, &
       current_state%parallel%monc_communicator, ierr)
  theta_l_mean(:)=theta_l_mean(:)*rnhpts

end subroutine calculate_theta_l_mean

subroutine dephy_update_socrates(socrates_opt,lat_traj,lon_traj,albedo_traj)
  implicit none
  type (str_socrates_options), intent(inout) :: socrates_opt
  real(kind=DEFAULT_PRECISION), intent(in) :: lat_traj
  real(kind=DEFAULT_PRECISION), intent(in) :: lon_traj
  real(kind=DEFAULT_PRECISION), intent(in) :: albedo_traj
  socrates_opt%latitude=lat_traj
  socrates_opt%longitude=lon_traj
  socrates_opt%surface_albedo=albedo_traj
end subroutine dephy_update_socrates


!~ !! SOME MORE DIRTY DIAGNOSTICS JUST ADDED AS COMMENTS

subroutine dephy_dirty_diagnostics()
    implicit none

    write(*,*) 'lat_traj'
    write(*,*) lat_traj
    write(*,*) 'lon_traj'
    write(*,*) lon_traj
    write(*,*) 'ps_forc'
    write(*,*) ps_forc
    write(*,*) 'ts'
    write(*,*) ts
    write(*,*) 'sfc_sens_flx'
    write(*,*) sfc_sens_flx
    write(*,*) 'sfc_lat_flx'
    write(*,*) sfc_lat_flx
    write(*,*) 'z0_traj'
    write(*,*) z0_traj
    write(*,*) 'z0th_traj'
    write(*,*) z0th_traj
    write(*,*) 'ustar'
    write(*,*) ustar
    write(*,*) 'u_traj'
    write(*,*) u_traj
    write(*,*) 'v_traj'
    write(*,*) v_traj
    write(*,*) 'albedo_traj'
    write(*,*) albedo_traj
    write(*,*) 'q_skin_traj'
    write(*,*) q_skin_traj

    write(*,*) 'height_forc'
    write(*,*) height_forc
    write(*,*) 'pressure_forc'
    write(*,*) pressure_forc
    write(*,*) 'ug'
    write(*,*) ug
    write(*,*) 'vg'
    write(*,*) vg
    write(*,*) 'u_adv'
    write(*,*) u_adv
    write(*,*) 'v_adv'
    write(*,*) v_adv
    write(*,*) 'theta_adv'
    write(*,*) theta_adv
    write(*,*) 'theta_rad'
    write(*,*) theta_rad
    write(*,*) 'rv_adv'
    write(*,*) rv_adv
    write(*,*) 'w'
    write(*,*) w
    write(*,*) 'theta_nudging'
    write(*,*) theta_nudging
    write(*,*) 'rv_nudging'
    write(*,*) rv_nudging
    write(*,*) 'u_nudging'
    write(*,*) u_nudging
    write(*,*) 'v_nudging'
    write(*,*) v_nudging
    write(*,*) 'nudging_inv_u_traj'
    write(*,*) nudging_inv_u_traj
    write(*,*) 'nudging_inv_v_traj'
    write(*,*) nudging_inv_v_traj
    write(*,*) 'nudging_inv_theta_traj'
    write(*,*) nudging_inv_theta_traj
    write(*,*) 'nudging_inv_rv_traj'
    write(*,*) nudging_inv_rv_traj

end subroutine dephy_dirty_diagnostics

end module dephy_forcings_mod
