module simplesetup_mod
  use datadefn_mod, only : DEFAULT_PRECISION
  use monc_component_mod, only : component_descriptor_type
  use state_mod, only : PRESCRIBED_SURFACE_FLUXES, model_state_type
  use conversions_mod, only : conv_to_string
  use logging_mod, only :  LOG_INFO, LOG_ERROR, log_log, log_master_log
  use grids_mod, only : local_grid_type, global_grid_type, X_INDEX, Y_INDEX, Z_INDEX, PRIMAL_GRID, DUAL_GRID
  use prognostics_mod, only : prognostic_field_type
  use optionsdatabase_mod, only : options_get_integer, options_get_logical, options_get_real, &
       options_get_integer_array, options_get_real_array
  use tracers_mod, only : get_tracer_options
  use q_indices_mod, only: get_q_index, standard_q_names
  use registry_mod, only : is_component_enabled

  implicit none

#ifndef TEST_MODE
  private
#endif

  integer :: x_size, y_size, z_size
  real(kind=DEFAULT_PRECISION) :: zztop, dxx, dyy
  logical :: enable_theta=.false.

  public simplesetup_get_descriptor
contains

  type(component_descriptor_type) function simplesetup_get_descriptor()
    simplesetup_get_descriptor%name="simplesetup"
    simplesetup_get_descriptor%version=0.1
    simplesetup_get_descriptor%initialisation=>initialisation_callback
  end function simplesetup_get_descriptor

  subroutine initialisation_callback(current_state)
    type(model_state_type), intent(inout), target :: current_state

    call read_configuration(current_state)
    if (.not. current_state%initialised) then
      current_state%dtm=options_get_real(current_state%options_database, "dtm")
      current_state%dtm_new=current_state%dtm
      call create_grid(current_state, current_state%global_grid)
      call decompose_grid(current_state)
      call allocate_prognostics(current_state)
      current_state%initialised=.true.

    end if

  end subroutine initialisation_callback

  
  subroutine allocate_prognostics(current_state)
    type(model_state_type), intent(inout) :: current_state

    integer :: alloc_z, alloc_y, alloc_x, i

    alloc_z=current_state%local_grid%size(Z_INDEX) + current_state%local_grid%halo_size(Z_INDEX) * 2
    alloc_y=current_state%local_grid%size(Y_INDEX) + current_state%local_grid%halo_size(Y_INDEX) * 2
    alloc_x=current_state%local_grid%size(X_INDEX) + current_state%local_grid%halo_size(X_INDEX) * 2
    
#ifdef U_ACTIVE
    call allocate_prognostic(current_state%u, alloc_z, alloc_y, alloc_x, DUAL_GRID, DUAL_GRID, PRIMAL_GRID)
    call allocate_prognostic(current_state%zu, alloc_z, alloc_y, alloc_x, DUAL_GRID, DUAL_GRID, PRIMAL_GRID)
    call allocate_prognostic(current_state%su, alloc_z, alloc_y, alloc_x, DUAL_GRID, DUAL_GRID, PRIMAL_GRID)
    call allocate_prognostic(current_state%savu, alloc_z, alloc_y, alloc_x, DUAL_GRID, DUAL_GRID, PRIMAL_GRID)
#endif
#ifdef V_ACTIVE
    call allocate_prognostic(current_state%v, alloc_z, alloc_y, alloc_x, DUAL_GRID, PRIMAL_GRID, DUAL_GRID)
    call allocate_prognostic(current_state%zv, alloc_z, alloc_y, alloc_x, DUAL_GRID, PRIMAL_GRID, DUAL_GRID)
    call allocate_prognostic(current_state%sv, alloc_z, alloc_y, alloc_x, DUAL_GRID, PRIMAL_GRID, DUAL_GRID)
    call allocate_prognostic(current_state%savv, alloc_z, alloc_y, alloc_x, DUAL_GRID, PRIMAL_GRID, DUAL_GRID)
#endif
#ifdef W_ACTIVE
    call allocate_prognostic(current_state%w, alloc_z, alloc_y, alloc_x, PRIMAL_GRID, DUAL_GRID, DUAL_GRID)    
    call allocate_prognostic(current_state%zw, alloc_z, alloc_y, alloc_x, PRIMAL_GRID, DUAL_GRID, DUAL_GRID)        
    call allocate_prognostic(current_state%sw, alloc_z, alloc_y, alloc_x, PRIMAL_GRID, DUAL_GRID, DUAL_GRID)        
    call allocate_prognostic(current_state%savw, alloc_z, alloc_y, alloc_x, PRIMAL_GRID, DUAL_GRID, DUAL_GRID)
#endif

    if (enable_theta) then
      call allocate_prognostic(current_state%th, alloc_z, alloc_y, alloc_x, DUAL_GRID, DUAL_GRID, DUAL_GRID)
      call allocate_prognostic(current_state%zth, alloc_z, alloc_y, alloc_x, DUAL_GRID, DUAL_GRID, DUAL_GRID)
      call allocate_prognostic(current_state%sth, alloc_z, alloc_y, alloc_x, DUAL_GRID, DUAL_GRID, DUAL_GRID)
    end if
    
    if (current_state%number_q_fields .gt. 0) then
      allocate(current_state%q(current_state%number_q_fields), &
               current_state%zq(current_state%number_q_fields),&
               current_state%sq(current_state%number_q_fields))
      do i=1, current_state%number_q_fields
        call allocate_prognostic(current_state%q(i), alloc_z, alloc_y, alloc_x, DUAL_GRID, DUAL_GRID, DUAL_GRID)
        call allocate_prognostic(current_state%zq(i), alloc_z, alloc_y, alloc_x, DUAL_GRID, DUAL_GRID, DUAL_GRID)
        call allocate_prognostic(current_state%sq(i), alloc_z, alloc_y, alloc_x, DUAL_GRID, DUAL_GRID, DUAL_GRID)
      end do

      ! Set standard q indices
      current_state%water_vapour_mixing_ratio_index=get_q_index(standard_q_names%VAPOUR, 'simplesetup')
      current_state%liquid_water_mixing_ratio_index=get_q_index(standard_q_names%CLOUD_LIQUID_MASS, 'simplesetup')
    end if
    
    if (current_state%n_tracers .gt. 0) then
      allocate( current_state%tracer(current_state%n_tracers),  &
                current_state%ztracer(current_state%n_tracers), &
                current_state%stracer(current_state%n_tracers))
      do i=1, current_state%n_tracers
        call allocate_prognostic(current_state%tracer(i), alloc_z, alloc_y, alloc_x, DUAL_GRID, DUAL_GRID, DUAL_GRID)
        call allocate_prognostic(current_state%ztracer(i), alloc_z, alloc_y, alloc_x, DUAL_GRID, DUAL_GRID, DUAL_GRID)
        call allocate_prognostic(current_state%stracer(i), alloc_z, alloc_y, alloc_x, DUAL_GRID, DUAL_GRID, DUAL_GRID)
      end do
      
    endif ! allocate tracers

    ! Set arrays for radiative heating rates - Note: this should be protected by a switch
    call allocate_prognostic(current_state%sth_lw, alloc_z, alloc_y, alloc_x, DUAL_GRID, DUAL_GRID, DUAL_GRID)
    call allocate_prognostic(current_state%sth_sw, alloc_z, alloc_y, alloc_x, DUAL_GRID, DUAL_GRID, DUAL_GRID) 

  end subroutine allocate_prognostics

  subroutine allocate_prognostic(field, alloc_z, alloc_y, alloc_x, z_grid, y_grid, x_grid)
    type(prognostic_field_type), intent(inout) :: field
    integer, intent(in) :: alloc_z, alloc_y, alloc_x, z_grid, y_grid, x_grid

    field%active=.true.
    field%grid(Z_INDEX) = z_grid
    field%grid(Y_INDEX) = y_grid
    field%grid(X_INDEX) = x_grid
    allocate(field%data(alloc_z, alloc_y, alloc_x))
    field%data=0.0_DEFAULT_PRECISION
  end subroutine allocate_prognostic  

  subroutine decompose_grid(current_state)
    type(model_state_type), intent(inout) :: current_state

    if (associated(current_state%parallel%decomposition_procedure)) then
      call current_state%parallel%decomposition_procedure(current_state)
    else
      call log_log(LOG_ERROR, "No decomposition specified")
    end if
  end subroutine decompose_grid

  subroutine create_grid(current_state, specific_grid)
    type(model_state_type), intent(inout) :: current_state
    type(global_grid_type), intent(inout) :: specific_grid

    integer, parameter :: KGD_SIZE=200
    integer :: number_kgd, i, kgd(KGD_SIZE)
    real(kind=DEFAULT_PRECISION) :: hgd(KGD_SIZE)

    kgd=-1

    call options_get_integer_array(current_state%options_database, "kgd", kgd)
    call options_get_real_array(current_state%options_database, "hgd", hgd)  

    if (kgd(1)==1)then
      if (hgd(1)/=0.0_DEFAULT_PRECISION)then
        call log_log(LOG_ERROR, "Lowest level is assumed to lie at the surface, check hgd(1)")
      else
        kgd(1:KGD_SIZE-1) = kgd(2:)
        hgd(1:KGD_SIZE-1) = hgd(2:)
      end if
    end if

    do i=1,size(kgd)
      if (kgd(i) == -1) exit      
    end do
    number_kgd=i-1

    if (number_kgd .gt. 0) then
      allocate(current_state%global_grid%configuration%vertical%kgd(number_kgd), &
           current_state%global_grid%configuration%vertical%hgd(number_kgd))
      current_state%global_grid%configuration%vertical%kgd=kgd(1:number_kgd)
      current_state%global_grid%configuration%vertical%hgd=hgd(1:number_kgd)
    end if

    specific_grid%bottom(Z_INDEX) = 0
    specific_grid%bottom(Y_INDEX) = 0
    specific_grid%bottom(X_INDEX) = 0

    specific_grid%top(Z_INDEX) = zztop
    specific_grid%top(Y_INDEX) = dyy * y_size
    specific_grid%top(X_INDEX) = dxx * x_size

    specific_grid%resolution(Z_INDEX) = zztop / z_size
    specific_grid%resolution(Y_INDEX) = dyy
    specific_grid%resolution(X_INDEX) = dxx

    specific_grid%size(Z_INDEX) = z_size
    specific_grid%size(Y_INDEX) = y_size
    specific_grid%size(X_INDEX) = x_size

    specific_grid%active(Z_INDEX) = .true.
    specific_grid%active(Y_INDEX) = .true.
    specific_grid%active(X_INDEX) = .true.

    specific_grid%dimensions = 3
  end subroutine create_grid

  subroutine read_configuration(current_state)
    type(model_state_type), intent(inout), target :: current_state

    current_state%rhobous=options_get_real(current_state%options_database, "rhobous")
    current_state%thref0=options_get_real(current_state%options_database, "thref0")
    current_state%number_q_fields=options_get_integer(current_state%options_database, "number_q_fields")
    current_state%surface_pressure=options_get_real(current_state%options_database, "surface_pressure")
    current_state%surface_reference_pressure=options_get_real(current_state%options_database, "surface_reference_pressure")

    current_state%use_anelastic_equations=options_get_logical(current_state%options_database, "use_anelastic_equations")

    current_state%origional_vertical_grid_setup=options_get_logical(current_state%options_database, &
         "origional_vertical_grid_setup")
    current_state%passive_q=options_get_logical(current_state%options_database, "passive_q")
    current_state%passive_th=options_get_logical(current_state%options_database, "passive_th")
    current_state%rmlmax=options_get_real(current_state%options_database, "rmlmax")
    current_state%calculate_th_and_q_init=options_get_logical(current_state%options_database, "calculate_th_and_q_init")
    current_state%use_viscosity_and_diffusion=options_get_logical(current_state%options_database, "use_viscosity_and_diffusion")
    current_state%backscatter=options_get_logical(current_state%options_database, "backscatter")

    x_size=options_get_integer(current_state%options_database, "x_size")
    y_size=options_get_integer(current_state%options_database, "y_size")
    z_size=options_get_integer(current_state%options_database, "z_size")
    dxx=options_get_real(current_state%options_database, "dxx")
    dyy=options_get_real(current_state%options_database, "dyy")
    zztop=options_get_real(current_state%options_database, "zztop")
    enable_theta=options_get_logical(current_state%options_database, "enable_theta")

    if (current_state%rmlmax<=0.0)current_state%rmlmax=0.23 * max(dxx, dyy)

    if (.not. enable_theta) current_state%passive_th=.true.
    if ( current_state%number_q_fields == 0) current_state%passive_q=.true.

    current_state%galilean_transformation=options_get_logical(current_state%options_database, "galilean_transformation")
    if (current_state%galilean_transformation)then
      current_state%fix_ugal=options_get_logical(current_state%options_database, "fix_ugal")
      current_state%fix_vgal=options_get_logical(current_state%options_database, "fix_vgal")
      if (current_state%fix_ugal)current_state%ugal=options_get_real(current_state%options_database, "ugal")
      if (current_state%fix_vgal)current_state%vgal=options_get_real(current_state%options_database, "vgal")
    end if

    current_state%print_debug_data = options_get_logical(current_state%options_database, "print_debug_data")

    if (.not. current_state%reconfig_run) then        
      call get_tracer_options(current_state)
    end if ! not reconfig

  end subroutine read_configuration
end module simplesetup_mod
