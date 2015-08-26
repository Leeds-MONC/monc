!> Component to set up the model based upon a KiD model configuration
!!
!! These data files are in NetCDF format
module kidreader_mod
  use monc_component_mod, only : component_descriptor_type
  use state_mod, only : model_state_type  
  use logging_mod, only : LOG_INFO, LOG_WARN, LOG_ERROR, LOG_DEBUG, log_master_log, log_log, log_get_logging_level
  use collections_mod, only : map_type
  use grids_mod, only : local_grid_type, global_grid_type, X_INDEX, Y_INDEX, Z_INDEX, PRIMAL_GRID, DUAL_GRID
  use prognostics_mod, only : prognostic_field_type
  use optionsdatabase_mod, only : options_get_string, options_get_real, options_get_integer, options_get_logical, &
       options_has_key, options_get_integer_array
  use datadefn_mod, only : DEFAULT_PRECISION, STRING_LENGTH
  use netcdf, only : nf90_noerr, nf90_global, nf90_nowrite, nf90_inquire_attribute, nf90_open, nf90_strerror, &
       nf90_inq_dimid, nf90_inquire_dimension, nf90_inq_varid, nf90_get_var, nf90_inquire, nf90_close, nf90_get_att
  implicit none

#ifndef TEST_MODE
  private
#endif

  integer, parameter :: NUMBER_Q_COORDS = 100             !< Number of Q field value coords that can be specified
  character(len=*), parameter ::  TIME_KEY = "time",&     !< Corresponding NetCDF data time key
       Z_KEY = "z",&           !< Corresponding NetCDF data z (primal grid) key
       Z_HALF_KEY = "z_half",& !< Corresponding NetCDF data z half (dual grid) key
       X_KEY = "x",&           !< Corresponding NetCDF data x (primal grid) key
       X_HALF_KEY = "x_half",& !< Corresponding NetCDF data x half (primal grid) key
       U_KEY="u",&             !< Corresponding NetCDF data u flow field key
       W_KEY="w"               !< Corresponding NetCDF data w flow field key

  logical :: flood_q = .false., float_q= .false., clone_to_3d = .false., rotate_xy=.false.
  integer :: domain_multiplication=1
  
  !< The Q field coordinates configured by the user
  integer, dimension(NUMBER_Q_COORDS), save :: q_coordinates_x, q_coordinates_y, q_coordinates_z, q_coordinates_value  
  !real :: rhobous, thref0, surface_pressure      !< Boussinesq density, reference potential temperature and surface pressure
  character(len=STRING_LENGTH) :: configuration_file         !< NetCDF model file to load

  public kidreader_get_descriptor

contains

  !> Provides the descriptor back to the caller and is used in component registration
  !! @returns The KidReader component descriptor
  type(component_descriptor_type) function kidreader_get_descriptor()
    kidreader_get_descriptor%name="kidreader"
    kidreader_get_descriptor%version=0.1
    kidreader_get_descriptor%initialisation=>initialise_callback
  end function kidreader_get_descriptor

  !> Initialisation hook which will parse the configuration NetCDF file and set up the model based
  !! upon this
  !! @param current_state The current model state_mod
  subroutine initialise_callback(current_state)
    type(model_state_type), target, intent(inout) :: current_state

    integer :: ncid, time_dim, z_dim, z_half_dim, x_dim, x_half_dim
    real, dimension(:), allocatable :: time, x, z, x_half, z_half
    real, dimension(:,:,:), allocatable :: u, w, v

    if (options_has_key(current_state%options_database, "restart")) then
      call log_master_log(LOG_DEBUG, "Ignoring KiD reader as restart from checkpoint file selected")
      return
    end if

    configuration_file=options_get_string(current_state%options_database, "kid_configuration_file")
    current_state%rhobous=options_get_real(current_state%options_database, "rhobous")
    current_state%thref0=options_get_real(current_state%options_database, "thref0")
    current_state%surface_pressure=options_get_real(current_state%options_database, "surface_pressure")
    flood_q=options_get_logical(current_state%options_database, "flood_q")
    float_q=options_get_logical(current_state%options_database, "float_q")
    clone_to_3d=options_get_logical(current_state%options_database, "clone_to_3d")
    rotate_xy=options_get_logical(current_state%options_database, "rotate_xy")
    domain_multiplication=options_get_integer(current_state%options_database, "domain_multiplication")

    call options_get_integer_array(current_state%options_database, "q_coordinates_x", q_coordinates_x)
    call options_get_integer_array(current_state%options_database, "q_coordinates_y", q_coordinates_y)
    call options_get_integer_array(current_state%options_database, "q_coordinates_z", q_coordinates_z)
    call options_get_integer_array(current_state%options_database, "q_coordinates_value", q_coordinates_value)

    call check_status(nf90_open(path = trim(configuration_file), mode = nf90_nowrite, ncid = ncid))
    call check_kinematics_file(ncid)
    if (log_get_logging_level() .ge. LOG_DEBUG) call read_global_attributes(ncid, current_state%parallel%my_rank)
    call read_dimensions(ncid, time_dim, z_dim, z_half_dim, x_dim, x_half_dim)
    call read_variables(ncid, time_dim, z_dim, z_half_dim, x_dim, x_half_dim, time, x, z, x_half, z_half, u, w, v)
    call check_status(nf90_close(ncid))

    call create_grid(current_state%global_grid, z_half, x_half, z_half_dim, x_half_dim, current_state%parallel%my_rank)
    call define_vertical_levels(current_state, z_half, z_half_dim)

    call decompose_grid(current_state)

    if (rotate_xy) then
      call initialise_velocity_field(current_state%local_grid, current_state%v, DUAL_GRID, PRIMAL_GRID, DUAL_GRID, v)
    else
      call initialise_velocity_field(current_state%local_grid, current_state%u, DUAL_GRID, DUAL_GRID, PRIMAL_GRID, u)
    end if
    call initialise_velocity_field(current_state%local_grid, current_state%w, PRIMAL_GRID, DUAL_GRID, DUAL_GRID, w)
    if (clone_to_3d) call initialise_velocity_field(current_state%local_grid, current_state%v, DUAL_GRID, PRIMAL_GRID, DUAL_GRID, v)
    call initalise_source_and_z_fields(current_state)
    call set_up_q_fields(current_state)
    current_state%initialised=.true.
    call log_master_log(LOG_INFO, "Initialised configuration from KiD model file `"//trim(configuration_file)//"`")
  end subroutine initialise_callback

  subroutine set_up_q_fields(current_state)
    type(model_state_type), intent(inout) :: current_state

    integer :: x_size, y_size, z_size, i

    z_size = current_state%local_grid%size(Z_INDEX) + current_state%local_grid%halo_size(Z_INDEX) * 2
    y_size = current_state%local_grid%size(Y_INDEX) + current_state%local_grid%halo_size(Y_INDEX) * 2
    x_size = current_state%local_grid%size(X_INDEX) + current_state%local_grid%halo_size(X_INDEX) * 2

    if (flood_q) current_state%number_q_fields=current_state%number_q_fields+1
    if (float_q) current_state%number_q_fields=current_state%number_q_fields+1
    allocate(current_state%q(current_state%number_q_fields), current_state%zq(current_state%number_q_fields), &
         current_state%sq(current_state%number_q_fields))
    do i=1,current_state%number_q_fields
      call initialise_single_q_field(current_state, i, z_size, y_size, x_size)
      if (flood_q .and. i == 1) current_state%q(i)%data(:,:,:) = 1.
      if (float_q .and. (.not. flood_q .or. i==2)) then
        call populate_q_tracer(current_state, current_state%q(i))
      end if
    end do
  end subroutine set_up_q_fields

  !> Populates the Q tracer field based upon the configuration that has been read in from the simulation file
  !! @param current_state The current model state_mod
  !! @param q_field The qfield that we are going to fill in
  subroutine populate_q_tracer(current_state, q_field)
    type(model_state_type), intent(inout) :: current_state
    type(prognostic_field_type), intent(inout) :: q_field

    integer :: i, x_local, y_local

    do i=1,NUMBER_Q_COORDS
      if (q_coordinates_x(i) .ne. -1 .and. q_coordinates_y(i) .ne. -1) then
        if (q_coordinates_x(i) .ne. 0 .and. .not. (q_coordinates_x(i) .ge. current_state%local_grid%start(X_INDEX) .and. &
             q_coordinates_x(i) .le. current_state%local_grid%end(X_INDEX))) cycle
        if (q_coordinates_y(i) .ne. 0 .and. .not. (q_coordinates_y(i) .ge. current_state%local_grid%start(Y_INDEX) .and. &
             q_coordinates_y(i) .le. current_state%local_grid%end(Y_INDEX))) cycle
        x_local=(q_coordinates_x(i) - (current_state%local_grid%start(X_INDEX)-1)) + current_state%local_grid%halo_size(X_INDEX)
        y_local=(q_coordinates_y(i) - (current_state%local_grid%start(Y_INDEX)-1)) + current_state%local_grid%halo_size(Y_INDEX)
        if (q_coordinates_z(i) == 0 .and. q_coordinates_y(i) == 0 .and. q_coordinates_x(i) == 0) then
          q_field%data(:,current_state%local_grid%local_domain_start_index(Y_INDEX):&
               current_state%local_grid%local_domain_end_index(Y_INDEX),current_state%local_grid%local_domain_start_index(X_INDEX):&
               current_state%local_grid%local_domain_end_index(X_INDEX)) = q_coordinates_value(i)
        else if (q_coordinates_z(i) == 0 .and. q_coordinates_y(i) == 0 .and. q_coordinates_x(i) .ne. 0) then
          q_field%data(:,current_state%local_grid%local_domain_start_index(Y_INDEX):&
               current_state%local_grid%local_domain_end_index(Y_INDEX),x_local) = q_coordinates_value(i)
        else if (q_coordinates_z(i) == 0 .and. q_coordinates_y(i) .ne. 0 .and. q_coordinates_x(i) == 0) then
          q_field%data(:,y_local,current_state%local_grid%local_domain_start_index(X_INDEX):&
               current_state%local_grid%local_domain_end_index(X_INDEX)) = q_coordinates_value(i)
        else if (q_coordinates_z(i) == 0 .and. q_coordinates_y(i) .ne. 0 .and. q_coordinates_x(i) .ne. 0) then
          q_field%data(:,y_local,x_local) = q_coordinates_value(i)
        else if (q_coordinates_z(i) .ne. 0 .and. q_coordinates_y(i) == 0 .and. q_coordinates_x(i) == 0) then
          q_field%data(q_coordinates_z(i),current_state%local_grid%local_domain_start_index(Y_INDEX):&
               current_state%local_grid%local_domain_end_index(Y_INDEX),current_state%local_grid%local_domain_start_index(X_INDEX):&
               current_state%local_grid%local_domain_end_index(X_INDEX)) = q_coordinates_value(i)
        else if (q_coordinates_z(i) .ne. 0 .and. q_coordinates_y(i) == 0 .and. q_coordinates_x(i) .ne. 0) then
          q_field%data(q_coordinates_z(i),current_state%local_grid%local_domain_start_index(Y_INDEX):&
               current_state%local_grid%local_domain_end_index(Y_INDEX),x_local) = q_coordinates_value(i)
        else if (q_coordinates_z(i) .ne. 0 .and. q_coordinates_y(i) .ne. 0 .and. q_coordinates_x(i) == 0) then
          q_field%data(q_coordinates_z(i),y_local,current_state%local_grid%local_domain_start_index(X_INDEX):&
               current_state%local_grid%local_domain_end_index(X_INDEX)) = q_coordinates_value(i)
        else if (q_coordinates_z(i) .ne. 0 .and. q_coordinates_y(i) .ne. 0 .and. q_coordinates_x(i) .ne. 0) then
          q_field%data(q_coordinates_z(i),y_local,x_local) = q_coordinates_value(i)
        end if
      end if
    end do
  end subroutine populate_q_tracer

  subroutine initialise_single_q_field(current_state, q_id, z_size, y_size, x_size)
    type(model_state_type), intent(inout) :: current_state
    integer, intent(in) :: q_id, x_size, y_size, z_size

    allocate(current_state%q(q_id)%data(z_size, y_size, x_size), current_state%zq(q_id)%data(z_size, y_size, x_size), &
         current_state%sq(q_id)%data(z_size, y_size, x_size))
    current_state%q(q_id)%data(:,:,:) = 0.
    current_state%zq(q_id)%data(:,:,:) = 0.
    current_state%sq(q_id)%data(:,:,:) = 0.
    current_state%q(q_id)%active=.true.
    current_state%zq(q_id)%active=.true.
    current_state%sq(q_id)%active=.true.
  end subroutine initialise_single_q_field

  !> Calls the decomposition procedure to decompose the grid and determine neighbouring processes
  !! If no decomposition procedure is specified then this results in an error
  !! @param current_state The current model state_mod
  subroutine decompose_grid(current_state)
    type(model_state_type), intent(inout) :: current_state

    if (associated(current_state%parallel%decomposition_procedure)) then
      call current_state%parallel%decomposition_procedure(current_state)
    else
      call log_master_log(LOG_ERROR, "No decomposition specified")
    end if
  end subroutine decompose_grid

  !> Based upon the local grid this will initialise the Source, Z and SAV fields for each prognostic
  !! @param current_state
  subroutine initalise_source_and_z_fields(current_state)
    type(model_state_type), intent(inout) :: current_state

    integer :: x_size, y_size, z_size

    z_size = current_state%local_grid%size(Z_INDEX) + current_state%local_grid%halo_size(Z_INDEX) * 2
    y_size = current_state%local_grid%size(Y_INDEX) + current_state%local_grid%halo_size(Y_INDEX) * 2
    x_size = current_state%local_grid%size(X_INDEX) + current_state%local_grid%halo_size(X_INDEX) * 2

#ifdef U_ACTIVE
    allocate(current_state%zu%data(z_size, y_size, x_size))
    allocate(current_state%su%data(z_size, y_size, x_size))
    allocate(current_state%savu%data(z_size, y_size, x_size))
    current_state%zu%data(:,:,:)= 0.0
#endif
#ifdef V_ACTIVE
    allocate(current_state%zv%data(z_size, y_size, x_size))
    allocate(current_state%sv%data(z_size, y_size, x_size))
    allocate(current_state%savv%data(z_size, y_size, x_size))
    current_state%zv%data(:,:,:)= 0.0
#endif
#ifdef W_ACTIVE
    allocate(current_state%zw%data(z_size, y_size, x_size))
    allocate(current_state%sw%data(z_size, y_size, x_size))
    allocate(current_state%savw%data(z_size, y_size, x_size))
    current_state%zw%data(:,:,:)= 0.0
#endif
    if (current_state%th%active) then
      allocate(current_state%zth%data(z_size, y_size, x_size))
      allocate(current_state%sth%data(z_size, y_size, x_size))
      current_state%zth%data(:,:,:)= 0.0
    end if
  end subroutine initalise_source_and_z_fields

  !> Will initialise a velocity field with the loaded data
  !! @param field The velocity field to initialise
  !! @param z_dim The size of the z dimension
  !! @param y_dim The size of the y dimension
  !! @param x_dim The size of the x dimension
  !! @param grid The applicable grid to place this field upon
  !! @param data The NetCDF KiD model data to load into the velocity field
  subroutine initialise_velocity_field(local_grid, field, z_grid, y_grid, x_grid, data)
    type(local_grid_type), intent(inout) :: local_grid
    type(prognostic_field_type), intent(inout) :: field
    integer, intent(in) :: z_grid, y_grid, x_grid
    real, dimension(:,:,:), allocatable, intent(in) :: data

    integer :: i, j, k, preMulSizeY, preMulSizeX

    field%grid(Z_INDEX) = z_grid
    field%grid(Y_INDEX) = y_grid
    field%grid(X_INDEX) = x_grid
    field%active = .true.

    allocate(field%data(local_grid%size(Z_INDEX) + local_grid%halo_size(Z_INDEX) * 2, local_grid%size(Y_INDEX) + &
         local_grid%halo_size(Y_INDEX) * 2, local_grid%size(X_INDEX) + local_grid%halo_size(X_INDEX) * 2))
    field%data=0.0

    ! Divisions here are for the multiplication - we just fill up the original size and then duplicate this across dimensions
    do i=ceiling(local_grid%start(X_INDEX)/real(domain_multiplication)),&
         ceiling(local_grid%end(X_INDEX)/real(domain_multiplication))
      do j=ceiling(local_grid%start(Y_INDEX)/real(domain_multiplication)),&
           ceiling(local_grid%end(Y_INDEX)/real(domain_multiplication))
        do k=local_grid%start(Z_INDEX),local_grid%end(Z_INDEX)          
          field%data(local_grid%halo_size(Z_INDEX)+(k-((local_grid%start(Z_INDEX)-1))), local_grid%halo_size(Y_INDEX)+&
               (j-((ceiling(local_grid%start(Y_INDEX)/ real(domain_multiplication))-1))), local_grid%halo_size(X_INDEX)+&
               (i- ((ceiling(local_grid%start(X_INDEX)/real(domain_multiplication))-1)))) = &
               real(data(1, k, merge(j, i, rotate_xy)), kind=DEFAULT_PRECISION)
        end do
      end do
    end do

    if (domain_multiplication .ge. 2) then
      ! If there is a domain multiplication then duplicate data across dimensions
      preMulSizeY = ceiling(local_grid%size(Y_INDEX) / real(domain_multiplication))
      preMulSizeX = ceiling(local_grid%size(X_INDEX) / real(domain_multiplication))
      if (local_grid%active(Y_INDEX)) then
        do i=1,domain_multiplication-1
          field%data(:, local_grid%halo_size(Y_INDEX) + i*preMulSizeY +1 : local_grid%halo_size(Y_INDEX) + (i+1)*preMulSizeY, &
               local_grid%halo_size(X_INDEX)+1: local_grid%halo_size(X_INDEX)+preMulSizeX) = &
               field%data(:, local_grid%halo_size(Y_INDEX)+1 : local_grid%halo_size(Y_INDEX) + preMulSizeY, &
               local_grid%halo_size(X_INDEX)+1 : local_grid%halo_size(X_INDEX)+preMulSizeX)
        end do
      end if
      if (local_grid%active(X_INDEX)) then
        do i=1,domain_multiplication-1
          field%data(:, local_grid%local_domain_start_index(Y_INDEX) : local_grid%local_domain_end_index(Y_INDEX), &
               local_grid%halo_size(X_INDEX) + i*preMulSizeX + 1 : local_grid%halo_size(X_INDEX) + (i+1)*preMulSizeX) = &
               field%data(:, local_grid%local_domain_start_index(Y_INDEX) : local_grid%local_domain_end_index(Y_INDEX), &
               local_grid%halo_size(X_INDEX)+1 : local_grid%halo_size(X_INDEX)+preMulSizeX)
        end do
      end if
    end if
  end subroutine initialise_velocity_field

  !> Creates a specific grid based upon the data read from the KiD model NetCDF file
  !! @param specific_grid The grid to create
  !! @param z Array of grid points in the z dimension
  !! @param x Array of grid points in the x dimension
  !! @param z_dim The number of points in the z dimension
  !! @param x_dim The number of points in the x dimension
  subroutine create_grid(specific_grid, z, x, z_dim, x_dim, my_rank)
    type(global_grid_type), intent(inout) :: specific_grid
    integer, intent(in) :: z_dim, x_dim, my_rank
    real, dimension(:), intent(in) :: z, x

    specific_grid%bottom(Z_INDEX) = int(z(1))
    specific_grid%bottom(merge(Y_INDEX, X_INDEX, rotate_xy)) = int(x(1))
    if (clone_to_3d) specific_grid%bottom(Y_INDEX) = int(x(1))

    specific_grid%top(Z_INDEX) = int(z(z_dim))
    specific_grid%top(merge(Y_INDEX, X_INDEX, rotate_xy)) = int(x(x_dim)) * domain_multiplication
    if (clone_to_3d) specific_grid%top(Y_INDEX) = int(x(x_dim)) * domain_multiplication

    specific_grid%resolution(Z_INDEX) = int(z(2) - z(1))
    specific_grid%resolution(merge(Y_INDEX, X_INDEX, rotate_xy)) = int(x(2) - x(1))
    if (clone_to_3d) specific_grid%resolution(Y_INDEX) = specific_grid%resolution(X_INDEX)

    specific_grid%size(Z_INDEX) = z_dim
    specific_grid%size(merge(Y_INDEX, X_INDEX, rotate_xy)) = x_dim * domain_multiplication
    if (clone_to_3d) specific_grid%size(Y_INDEX) = x_dim * domain_multiplication

    specific_grid%active(Z_INDEX) = .true.
    specific_grid%active(merge(Y_INDEX, X_INDEX, rotate_xy)) = .true.
    if (clone_to_3d) specific_grid%active(Y_INDEX) = .true.

#ifdef U_ACTIVE
    if (.not. specific_grid%active(X_INDEX)) call log_master_log(LOG_ERROR, &
         "Model compiled with X active but inactive in configuration")
#else
    if (specific_grid%active(X_INDEX)) call log_master_log(LOG_ERROR, &
         "Model compiled with X inactive but active in configuration")
#endif
#ifdef V_ACTIVE
    if (.not. specific_grid%active(Y_INDEX)) call log_master_log(LOG_ERROR, &
         "Model compiled with Y active but inactive in configuration")
#else
    if (specific_grid%active(Y_INDEX)) call log_master_log(LOG_ERROR, &
         "Model compiled with Y inactive but active in configuration")
#endif
#ifdef W_ACTIVE
    if (.not. specific_grid%active(Z_INDEX)) call log_master_log(LOG_ERROR, &
         "Model compiled with Z active but inactive in configuration")
#else
    if (specific_grid%active(Z_INDEX)) call log_master_log(LOG_ERROR, &
         "Model compiled with Z inactive but active in configuration")
#endif
    specific_grid%dimensions = merge(3, 2, clone_to_3d)
  end subroutine create_grid

  !> Defines the vertical levels of the grid. This is both the grid points and corresponding height
  !! for each point in metres
  !! @param current_state The current model state_mod
  !! @param z Array of grid heights in metres
  !! @param z_size Number of grid points
  subroutine define_vertical_levels(current_state, z, z_size)
    type(model_state_type), intent(inout) :: current_state
    integer, intent(in) :: z_size
    real, dimension(:), intent(in) :: z

    integer :: i

    allocate(current_state%global_grid%configuration%vertical%kgd(z_size), &
         current_state%global_grid%configuration%vertical%hgd(z_size))

    do i=1,z_size
      current_state%global_grid%configuration%vertical%kgd(i) = i
      current_state%global_grid%configuration%vertical%hgd(i) = real(z(i))
    end do
  end subroutine define_vertical_levels

  !> Checks that the kinematics file that has been loaded is consistent with what we expect
  !!
  !! This checks the file against a number of limitations of the current KiD reading functionality
  !! to ensure that this component can parse it correctly
  subroutine check_kinematics_file(ncid)
    integer, intent(in) :: ncid

    integer :: ndims_in, nvars_in, ngatts_in, unlimdimid_in

    call check_status(nf90_inquire(ncid, ndims_in, nvars_in, ngatts_in, unlimdimid_in))
    if (ndims_in /= 5) call log_log(LOG_ERROR, "NetCDF KiD number of model dimensions must equal 5")
    if (nvars_in /= 7) call log_log(LOG_ERROR, "NetCDF KiD number of model variables must equal 5")
    if (ngatts_in .le. 0) call log_log(LOG_ERROR, "NetCDF KiD global attributes must be specified")
    if (unlimdimid_in .gt. 0) call log_log(LOG_ERROR, "NetCDF KiD model number of unlimited dimensions must be 0")
  end subroutine check_kinematics_file

  !> Reads the variables from the NetCDF KiD model file
  !! @param ncid The id of the NetCDF file
  !! @param time_dim The number of elements in the time dimension
  !! @param z_dim The number of elements in the z primal grid dimension
  !! @param z_half_dim The number of elements in the z dual grid dimension
  !! @param x_dim The number of elements in the x primal grid dimension
  !! @param x_half_dim The number of elements in the x dual grid dimension
  !! @param time The time data field that is to be read
  !! @param x The points on the primal grid in x dimension
  !! @param z The points on the primal grid z dimension
  !! @param x_half The points on the dual grid x dimension
  !! @param z_half The points on the dual grid z dimension
  !! @param u Velocity field in dimension x to read
  !! @param w Velocity field in dimension z to read
  subroutine read_variables(ncid, time_dim, z_dim, z_half_dim, x_dim, x_half_dim, time, x, z, x_half, z_half, u, w, v)
    integer, intent(in) :: ncid, time_dim, z_dim, z_half_dim, x_dim, x_half_dim
    real, dimension(:), allocatable, intent(inout) :: time, x, z, x_half, z_half
    real, dimension(:,:,:), allocatable, intent(inout) :: u, w, v

    allocate(time(time_dim))
    allocate(z(z_dim))
    allocate(x(x_dim))
    allocate(z_half(z_half_dim))
    allocate(x_half(x_half_dim))
    ! Due to NetCDF being in C, need to reverse the data order for F2003
    if (rotate_xy) then
      allocate(v(time_dim, z_dim, x_half_dim))
    else
      allocate(u(time_dim, z_dim, x_half_dim))
    end if
    allocate(w(time_dim, z_half_dim, x_half_dim))

    call read_single_variable(ncid, TIME_KEY, data1d=time)
    call read_single_variable(ncid, Z_KEY, data1d=z)
    call read_single_variable(ncid, Z_HALF_KEY, data1d=z_half)
    call read_single_variable(ncid, X_KEY, data1d=x)
    call read_single_variable(ncid, X_HALF_KEY, data1d=x_half)
    if (rotate_xy) then
      call read_single_variable(ncid, U_KEY, data3d=v)
    else
      call read_single_variable(ncid, U_KEY, data3d=u)
    end if
    call read_single_variable(ncid, W_KEY, data3d=w)

    if (clone_to_3d) then
      if (.not. allocated(v)) then
        allocate(v(time_dim, z_dim, x_half_dim))
        v=u
      end if
      if (.not. allocated(u)) then
        allocate(u(time_dim, z_dim, x_half_dim))
        u=v
      end if
    end if
  end subroutine read_variables

  !> Reads a single variable out of a NetCDF file
  !! @param ncid The NetCDF file id
  !! @param key The variable key (name) to access
  !! @param data1d Optional one dimensional data to read into
  !! @param data3d Optional three dimensional data to read into
  subroutine read_single_variable(ncid, key, data1d, data3d)
    integer, intent(in) :: ncid
    character(len=*), intent(in) :: key
    real, dimension(:), intent(inout), optional :: data1d
    real, dimension(:,:,:), intent(inout), optional :: data3d

    integer :: variable_id
    real, dimension(:,:,:), allocatable :: sdata

    call check_status(nf90_inq_varid(ncid, key, variable_id))

    if (.not. present(data1d) .and. .not. present(data3d)) return

    if (present(data1d)) then
      call check_status(nf90_get_var(ncid, variable_id, data1d))
    else
      ! 3D will reshape the data to take account of the column-row major C-F transposition
      allocate(sdata(size(data3d,1),size(data3d,3), size(data3d,2)))
      call check_status(nf90_get_var(ncid, variable_id, sdata))
      data3d(:,:,:)=reshape(sdata(:,:,:),(/size(data3d,1),size(data3d,2),size(data3d,3)/))
      deallocate(sdata)
    end if
  end subroutine read_single_variable

  !> Reads the dimensions from the NetCDF file
  !! @param ncid The NetCDF file id
  !! @param time_dim Number of elements in the time dimension
  !! @param z_dim Number of elements in the z dimension of the primal grid
  !! @param z_half_dim Number of elements in the z dimension of the dual grid
  !! @param x_dim Number of elements in the x dimension of the primal grid
  !! @param x_half_dim Number of elements in the x dimension of the dual grid
  subroutine read_dimensions(ncid, time_dim, z_dim, z_half_dim, x_dim, x_half_dim)
    integer, intent(in) :: ncid
    integer, intent(out) ::  time_dim, z_dim, z_half_dim, x_dim, x_half_dim

    integer ::  time_dimid, z_dimid, z_half_dimid, x_dimid, x_half_dimid

    call check_status(nf90_inq_dimid(ncid, TIME_KEY, time_dimid))
    call check_status(nf90_inq_dimid(ncid, Z_KEY, z_dimid))
    call check_status(nf90_inq_dimid(ncid, Z_HALF_KEY, z_half_dimid))
    call check_status(nf90_inq_dimid(ncid, X_KEY, x_dimid))
    call check_status(nf90_inq_dimid(ncid, X_HALF_KEY, x_half_dimid))

    call check_status(nf90_inquire_dimension(ncid, time_dimid, len=time_dim))
    call check_status(nf90_inquire_dimension(ncid, z_dimid, len=z_dim))
    call check_status(nf90_inquire_dimension(ncid, z_half_dimid, len=z_half_dim))
    call check_status(nf90_inquire_dimension(ncid, x_dimid, len=x_dim))
    call check_status(nf90_inquire_dimension(ncid, x_half_dimid, len=x_half_dim))
  end subroutine read_dimensions

  !> Will read the global attributes of the NetCDF KiD model dump and log_log them to debug
  !! @param ncid The NetCDF file id
  subroutine read_global_attributes(ncid, pid)
    integer, intent(in) :: ncid, pid

    character(len=:), allocatable :: attributeValue

    attributeValue=read_specific_global_attribute(ncid, "title")
    call log_master_log(LOG_DEBUG, "KiD file title: "//attributeValue)
    deallocate(attributeValue)

    attributeValue=read_specific_global_attribute(ncid, "creation date")
    call log_master_log(LOG_DEBUG, "KiD file created: "//attributeValue)
    deallocate(attributeValue)
  end subroutine read_global_attributes

  !> Will read a specific global NetCDF attribute
  !! @param ncid The NetCDF file id
  !! @param key The name (key) of the attribute to read
  !! @returns String representation of the attribute value
  function read_specific_global_attribute(ncid, key)
    integer, intent(in) :: ncid
    character(len=*), intent(in) :: key

    integer :: length
    character(len=:),allocatable,target :: read_specific_global_attribute

    call check_status(nf90_inquire_attribute(ncid, nf90_global, key, len = length))
    allocate(character(length) :: read_specific_global_attribute)
    call check_status(nf90_get_att(ncid, nf90_global, key, read_specific_global_attribute))
  end function read_specific_global_attribute

  !> Will check a NetCDF status and write to log_log error any decoded statuses
  !! @param status The NetCDF status flag
  subroutine check_status(status)
    integer, intent(in) :: status

    if (status /= nf90_noerr) then
      call log_log(LOG_ERROR, "NetCDF returned error code of "//trim(nf90_strerror(status)))
    end if
  end subroutine check_status
end module kidreader_mod
