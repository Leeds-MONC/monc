!> The model state which represents the current state of a run
module state_mod
  use collections_mod, only : hashmap_type
  use grids_mod, only : global_grid_type, local_grid_type
  use prognostics_mod, only : prognostic_field_type
  use communication_types_mod, only : halo_communication_type
  use datadefn_mod, only : DEFAULT_PRECISION
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
          my_rank,&           !> My process rank in the MONC system
          my_global_rank,&    !> My process rank in the global system
          neighbour_comm,&     !> Neighbour communicator
          monc_communicator=-1, io_communicator=-1, corresponding_io_server_process
     integer, dimension(3) :: &
          my_coords,&         !> My process coordinates in each dimension
          dim_sizes           !> Number of processes in each dimension
     logical, dimension(3,2) :: wrapped_around
     procedure(), nopass, pointer :: decomposition_procedure => null() !> The decomposition procedure to use
  end type parallel_state_type


  !> Information about the non-zero sampling intervals
  type, public :: sampling_interval_type
    integer :: interval = 0        ! sampling interval [ts, s if time_basis]
    integer, dimension(:), allocatable :: output   ! output intervals associated with %interval [s]
                                                   ! nint(output_frequency)
    integer :: next_time = 0       ! if time_basis, the next sample time for this %interval [s]
                                   ! if force_output_on_interval, the next output time for 
                                   ! this %interval [s]
    integer :: next_step = 0       ! the next sample timestep for this %interval [ts]
    logical :: active = .false.    ! .true. when sampling for this %interval on the current timestep
  end type sampling_interval_type
  
  !> The ModelState which represents the current state of a run
  !!
  !! This state is provided to each callback and may be used and modified as required by
  !! the callbacks. Apart from this state, there should be no other state (global) variables
  !! declared. This allows us to simply persist and retrieve the ModelState when suspending
  !! and reactivating MONC.
  type, public :: model_state_type
    logical :: continue_timestep=.true., initialised=.false., continuation_run=.false.
    
    logical :: reconfig_run=.false.    ! whether this is the first cycle of a reconfigured 
                                       ! continuation run
    logical :: retain_model_time=.false.    ! by default, reconfigurations have model time 
                                            ! reset to zero 
    logical :: only_compute_on_sample_timestep=.false.    ! by default, diagnostics are available
                                  ! on every timestep.  When .true., certain diagnostics are only
                                  ! computed on specified diagnostic_sample_timesteps
    logical :: diagnostic_sample_timestep=.false.    ! diagnostics should be computed on the 
                                                     ! current timestep
    logical :: normal_step=.true.    ! not a special, shortened timestep due to proximity to a
                                     ! diagnostic_sample_timestep
    logical :: force_output_on_interval=.false.   ! allows the model to adjust the dtm to 
                                                  ! ensure that samples are sent to the IO
                                                  ! server on the output_frequency
                                                  ! time_basis=.true. does this automatically
    logical :: print_debug_data=.false.    ! Prints data for specific variables/points for
                                           ! debugging.  See registry.F90:execute_callbacks
    logical :: use_viscosity_and_diffusion=.true., &
       use_surface_boundary_conditions=.false., backscatter=.true.

    type(hashmap_type) :: options_database
    type(global_grid_type) :: global_grid
    type(local_grid_type) :: local_grid
    type(parallel_state_type) :: parallel
    type(prognostic_field_type) :: u, w, v, th, p, zu, zw, zv, zth, su, sw, sv, sth, savu, savv, & 
         savw, vis_coefficient, &
         diff_coefficient, dis, dis_th, &
         ! Heating rates from socrates
         sth_lw, sth_sw
    type(prognostic_field_type), dimension(:), allocatable :: q, zq, sq, disq
    type(prognostic_field_type), dimension(:), allocatable :: tracer, ztracer, stracer
    ! longwave and shortwave downwelling flux at the surface
    real(kind=DEFAULT_PRECISION), dimension(:,:), allocatable :: sw_down_surf, lw_down_surf
    type(halo_communication_type) :: viscosity_halo_swap_state, diffusion_halo_swap_state
    type(sampling_interval_type), dimension(:), allocatable :: sampling
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
    real(kind=DEFAULT_PRECISION), dimension(:), allocatable :: tracer_decay_rate, tracer_surf_bc
    integer :: timestep=1, column_global_x, column_global_y, column_local_x, column_local_y, field_stepping, scalar_stepping, &
         momentum_stepping, number_q_fields=0, start_timestep=1, type_of_surface_boundary_conditions, lookup_table_entries, &
         cfl_frequency, termination_reason, last_cfl_timestep=0
    integer :: water_vapour_mixing_ratio_index=0, liquid_water_mixing_ratio_index=0, &
         rain_water_mixing_ratio_index=0, ice_water_mixing_ratio_index=0, &
         snow_water_mixing_ratio_index=0, graupel_water_mixing_ratio_index=0, & 
         psrce_x_hs_send_request, psrce_y_hs_send_request, psrce_x_hs_recv_request, psrce_y_hs_recv_request
    integer :: n_tracers=0, n_radioactive_tracers=0
    integer :: traj_tracer_index=0, radioactive_tracer_index=0
    integer, dimension(:), allocatable:: tracer_surf_bc_opt
    logical :: first_timestep_column, last_timestep_column, halo_column, first_nonhalo_timestep_column, &
         passive_q=.false., passive_th=.false., &
         use_time_varying_surface_values, use_anelastic_equations, & ! use_anelastic_equations or use Boussinesq
         saturated_surface, update_dtm=.false., calculate_th_and_q_init, origional_vertical_grid_setup=.true., &
         io_server_enabled, reinit_tracer=.false., time_basis=.false.
    logical, allocatable :: l_forceq(:)
    double precision :: model_start_wtime

    logical :: galilean_transformation=.true., fix_ugal=.false., fix_vgal=.false.
    real(kind=DEFAULT_PRECISION) :: ugal=0.,vgal=0.
    ! SOCRATES time variables are included in state since they need to be dumped
    real(kind=DEFAULT_PRECISION) :: rad_last_time=0.0

  end type model_state_type
end module state_mod
