!> The model state which represents the current state of a run
module state_mod
  use collections_mod, only : hashmap_type
  use grids_mod, only : global_grid_type, local_grid_type
  use prognostics_mod, only : prognostic_field_type
  use communication_types_mod, only : halo_communication_type
  use datadefn_mod, only : DEFAULT_PRECISION, STRING_LENGTH
  use logging_mod, only : LOG_WARN, LOG_ERROR, log_log
  implicit none

#ifndef TEST_MODE
  private
#endif

  !> Stepping parameter values which determine centred or forward stepping
  integer, parameter, public :: CENTRED_STEPPING=0, FORWARD_STEPPING=1, PRESCRIBED_SURFACE_FLUXES=0, PRESCRIBED_SURFACE_VALUES=1
  !> The constants defining the reason why the model has terminated
  integer, parameter, public :: TIME_TERMINATION_REASON=0, TIMESTEP_TERMINATION_REASON=1, MESSAGE_TERMINATION_REASON=2, &
       WALLTIME_TERMINATION_REASON=3

  !> Information about the parallel aspects of the system
  type, public :: parallel_state_type
     integer :: processes, & !> Total number of processes
          my_rank,&           !> My process rank in the system
          neighbour_comm,&     !> Neighbour communicator
          monc_communicator=-1, io_communicator=-1, corresponding_io_server_process
     integer, dimension(3) :: &
          my_coords,&         !> My process coordinates in each dimension
          dim_sizes           !> Number of processes in each dimension
     logical, dimension(3,2) :: wrapped_around
     procedure(), nopass, pointer :: decomposition_procedure => null() !> The decomposition procedure to use
  end type parallel_state_type
  
  !> The ModelState which represents the current state of a run
  !!
  !! This state is provided to each callback and may be used and modified as required by
  !! the callbacks. Apart from this state, there should be no other state (global) variables
  !! declared. This allows us to simply persist and retrieve the ModelState when suspending
  !! and reactivating MONC.
  type, public :: model_state_type
    logical :: continue_timestep=.true., initialised=.false., continuation_run=.false.
    logical :: use_viscosity_and_diffusion=.true., &
       use_surface_boundary_conditions=.true., backscatter=.true.

    type(hashmap_type) :: options_database
    type(global_grid_type) :: global_grid
    type(local_grid_type) :: local_grid
    type(parallel_state_type) :: parallel
    type(prognostic_field_type) :: u, w, v, th, p, zu, zw, zv, zth, su, sw, sv, sth, savu, savv, savw, vis_coefficient, &
         diff_coefficient, dis, dis_th
    type(prognostic_field_type), dimension(:), allocatable :: q, zq, sq, disq
    type(halo_communication_type) :: viscosity_halo_swap_state, diffusion_halo_swap_state
    real(kind=DEFAULT_PRECISION) :: time=.0_DEFAULT_PRECISION,& ! Model time in seconds
            dtm,& ! Modeltimestep (s)
            absolute_new_dtm, &
            thref0,&
            rhobous,&
            tsmth=1e-2_DEFAULT_PRECISION,&
            timestep_runtime,&
            local_divmax, global_divmax, cvis=0.0_DEFAULT_PRECISION, surface_temperature_flux, &
            surface_vapour_flux, theta_surf, surface_vapour_mixing_ratio, fbuoy, &
            fbuoynew, theta_virtual_surf, cmbc, rcmbc, ellmocon, velmax, velmin, aloginv, cneut, cfc, &
            surface_pressure=100000.0_DEFAULT_PRECISION, surface_reference_pressure = 100000.0_DEFAULT_PRECISION, &
            cvel, cvel_x, cvel_y, cvel_z, dtm_new, rmlmax, geostrophic_wind_rate_of_change_in_x, &
            geostrophic_wind_rate_of_change_in_y, surface_geostrophic_wind_x, surface_geostrophic_wind_y, &
            local_zumin, local_zumax, local_zvmin, local_zvmax, local_cvel_z
    real(kind=DEFAULT_PRECISION), dimension(:), allocatable :: lookup_table_velocity, &
         lookup_table_ustr, cq, abswmax
    real(kind=DEFAULT_PRECISION), dimension(:,:), allocatable :: psrce_recv_buffer_x, psrce_recv_buffer_y
    integer :: timestep=1, column_global_x, column_global_y, column_local_x, column_local_y, field_stepping, scalar_stepping, &
         momentum_stepping, number_q_fields=0, start_timestep=1, type_of_surface_boundary_conditions, lookup_table_entries, &
         cfl_frequency, termination_reason
    integer :: water_vapour_mixing_ratio_index=0, liquid_water_mixing_ratio_index=0, &
         psrce_x_hs_send_request, psrce_y_hs_send_request, psrce_x_hs_recv_request, psrce_y_hs_recv_request
    logical :: first_timestep_column, last_timestep_column, halo_column, first_nonhalo_timestep_column, &
         passive_q=.false., passive_th=.false., &
         use_time_varying_surface_values, use_anelastic_equations, & ! use_anelastic_equations or use Boussinesq
         saturated_surface, update_dtm=.false., calculate_th_and_q_init, origional_vertical_grid_setup=.true., &
         io_server_enabled
    logical, allocatable :: l_forceq(:)
    double precision :: model_start_wtime

    logical :: galilean_transformation=.true., fix_ugal=.false., fix_vgal=.false.
    real(kind=DEFAULT_PRECISION) :: ugal=0.,vgal=0.
  end type model_state_type

  public get_prognostic_field_units, get_prognostic_field_standard_name, &
         get_prognostic_field_long_name

contains

  !> Provide CF-compliant units for prognostic variables
  !!
  !! @param field_name name of field to look up
  !! @returns units of field
  function get_prognostic_field_units(field_name) result(field_units)
    character(len=*), intent(in) :: field_name
    character(len=STRING_LENGTH) :: field_units

    field_units = ""
    if (trim(field_name) == "p") then
       field_units = "Pa"
    else if (trim(field_name) == "th") then
       field_units = "K"
    else if (trim(field_name) == "q") then
       field_units = "kg/kg"
    else if (trim(field_name) == "thref") then
       field_units = "K"
    else if (trim(field_name) == "u") then
       field_units = "m/s"
    else if (trim(field_name) == "v") then
       field_units = "m/s"
    else if (trim(field_name) == "w") then
       field_units = "m/s"
    else if (trim(field_name) == "x") then
       field_units = "m"
    else if (trim(field_name) == "y") then
       field_units = "m"
    else if (trim(field_name) == "z") then
       field_units = "m"
    else
       call log_log(LOG_WARN, "Units not defined for globally sendable field "//trim(field_name))
    endif

  end function get_prognostic_field_units


  !> Get "long name" (description of a field in terminology of CF conventions)
  !! for prognostic fields globally available for IO
  !! @param field_name name of field to look up
  !! @returns long name of field
  function get_prognostic_field_long_name(field_name) result(field_long_name)
    character(len=*), intent(in) :: field_name
    character(len=STRING_LENGTH) :: field_long_name

    field_long_name = ""
    if (trim(field_name) == "p") then
       field_long_name = "Pressure perturbation"
    else if (trim(field_name) == "th") then
       field_long_name = "Potential temperature perturbation"
    else if (trim(field_name) == "q") then
       field_long_name = "Tracer concentration"
    else if (trim(field_name) == "thref") then
       field_long_name = "Potential temperature vertical profile"
    else if (trim(field_name) == "u") then
       field_long_name = "Zonal wind"
    else if (trim(field_name) == "v") then
       field_long_name = "Meridional wind"
    else if (trim(field_name) == "w") then
       field_long_name = "Vertical velocity"
    else if (trim(field_name) == "x") then
       field_long_name = "East-west displacement of cell centers"
    else if (trim(field_name) == "y") then
       field_long_name = "North-wouth displacement of cell centers"
    else if (trim(field_name) == "z") then
       field_long_name = "Vertical displacement of cell centers"
    else
       call log_log(LOG_WARN, "Description not defined for globally sendable field "//trim(field_name))
    endif

  end function get_prognostic_field_long_name

  !> Provide CF-compliant standard names for prognostic variables
  !! Variables are taken from table at
  !! http://cfconventions.org/Data/cf-standard-names/43/build/cf-standard-name-table.html
  !!
  !! @param field_name name of field to look up
  !! @returns standard name of field
  function get_prognostic_field_standard_name(field_name) result(field_standard_name)
    character(len=*), intent(in) :: field_name
    character(len=STRING_LENGTH) :: field_standard_name

    field_standard_name = ""
    if (trim(field_name) == "p") then
       field_standard_name = "" ! NB: not "air_pressure", since `p` is actually the perturbation in pressure
    else if (trim(field_name) == "th") then
       field_standard_name = "" ! NB: similarly to `p` this variable is the perturbation in potential temperature
    else if (trim(field_name) == "q") then
       field_standard_name = ""  ! no standard name
    else if (trim(field_name) == "thref") then
       field_standard_name = ""
    else if (trim(field_name) == "u") then
       field_standard_name = "eastward_wind"
    else if (trim(field_name) == "v") then
       field_standard_name = "northward_wind"
    else if (trim(field_name) == "w") then
       field_standard_name = "upward_air_velocity"
    else if (trim(field_name) == "x") then
       field_standard_name = "x-coordinate in Cartesian system (cell-centers)"
    else if (trim(field_name) == "y") then
       field_standard_name = "y-coordinate in Cartesian system (cell-centers)"
    else if (trim(field_name) == "z") then
       field_standard_name = "height-coordinate in Cartesian system (cell-edges)"
    else if (trim(field_name) == "zn") then
       field_standard_name = "height-coordinate in Cartesian system (cell-centers)"
    else
       call log_log(LOG_WARN, "Standard name not defined for globally sendable field "//trim(field_name))
    endif

  end function get_prognostic_field_standard_name
end module state_mod
