!> Will read in a NetCDF checkpoint file and initialise the model state_mod based upon this
module checkpointer_read_checkpoint_mod
#ifndef TEST_MODE
  use netcdf, only : nf90_global, nf90_nowrite, nf90_inquire_attribute, nf90_open, nf90_inq_dimid, nf90_inquire_dimension, &
       nf90_inq_varid, nf90_get_var, nf90_get_att, nf90_close
#else
  use dummy_netcdf_mod, only : nf90_global, nf90_nowrite, nf90_inquire_attribute, nf90_open, nf90_inq_dimid, &
       nf90_inquire_dimension, nf90_inq_varid, nf90_get_var, nf90_get_att, nf90_close
#endif
  use state_mod, only : model_state_type
  use grids_mod, only : local_grid_type, global_grid_type, X_INDEX, Y_INDEX, Z_INDEX, PRIMAL_GRID, DUAL_GRID
  use prognostics_mod, only : prognostic_field_type
  use logging_mod, only : LOG_INFO, LOG_ERROR, log_log, log_master_log
  use conversions_mod, only : conv_is_integer, conv_to_integer, conv_is_real, conv_to_real, conv_is_logical, conv_to_logical
  use optionsdatabase_mod, only : options_add
  use checkpointer_common_mod, only : EMPTY_DIM_KEY, STRING_DIM_KEY, X_DIM_KEY, Y_DIM_KEY, &
       Z_DIM_KEY, Q_DIM_KEY, Q_KEY, ZQ_KEY, TH_KEY, ZTH_KEY, P_KEY, U_KEY, V_KEY, W_KEY, ZU_KEY, ZV_KEY, ZW_KEY, X_KEY, Y_KEY, &
       Z_KEY, NQFIELDS, UGAL, VGAL, TIME_KEY, TIMESTEP, CREATED_ATTRIBUTE_KEY, TITLE_ATTRIBUTE_KEY, &
       ABSOLUTE_NEW_DTM_KEY, DTM_KEY, DTM_NEW_KEY, Q_INDICES_KEY, check_status, remove_null_terminator_from_string
  use datadefn_mod, only : DEFAULT_PRECISION, STRING_LENGTH
  use q_indices_mod, only : set_q_index
  implicit none

#ifndef TEST_MODE
  private
#endif

  public read_checkpoint_file

contains

  !> Reads in a NetCDF checkpoint file and uses this to initialise the model
  !! @param currentState The current model state_mod
  !! @param filename The filename of the checkpoint file to load
  subroutine read_checkpoint_file(current_state, filename)
    type(model_state_type), intent(inout) :: current_state
    character(len=*), intent(in) :: filename

    integer :: ncid, z_dim, y_dim, x_dim
    logical :: z_found, y_found, x_found
    character(len=:), allocatable :: attribute_value

    call check_status(nf90_open(path = filename, mode = nf90_nowrite, ncid = ncid))
    attribute_value=read_specific_global_attribute(ncid, "created")
    call read_dimensions(ncid, z_dim, y_dim, x_dim, z_found, y_found, x_found)
    call load_global_grid(current_state, ncid, z_dim, y_dim, x_dim, z_found, y_found, x_found)

    call decompose_grid(current_state)

    call load_q_indices(ncid)
    call load_misc(current_state, ncid)
    call load_all_fields(current_state, ncid)
    call initalise_source_and_sav_fields(current_state)  
    call check_status(nf90_close(ncid))

    call log_master_log(LOG_INFO, "Restarted configuration from checkpoint file `"//trim(filename)//"` created at "&
         //attribute_value)
    deallocate(attribute_value)
    current_state%initialised=.true.
  end subroutine read_checkpoint_file

  !> Calls out to the parallel decomposition strategy to decompose the grid based upon the configuration read and number of
  !! processes
  !! @param current_state The current model state
  subroutine decompose_grid(current_state)
    type(model_state_type), intent(inout) :: current_state

    if (associated(current_state%parallel%decomposition_procedure)) then
      call current_state%parallel%decomposition_procedure(current_state)
    else
      call log_log(LOG_ERROR, "No decomposition specified")
    end if
  end subroutine decompose_grid

  !> Initialises the source and sav (for u,v,w) fields for allocated prognostics
  !! @param current_state The model current state
  subroutine initalise_source_and_sav_fields(current_state)
    type(model_state_type), intent(inout) :: current_state

    integer :: x_size, y_size, z_size, i

    z_size = current_state%local_grid%size(Z_INDEX) + current_state%local_grid%halo_size(Z_INDEX) * 2
    y_size = current_state%local_grid%size(Y_INDEX) + current_state%local_grid%halo_size(Y_INDEX) * 2
    x_size = current_state%local_grid%size(X_INDEX) + current_state%local_grid%halo_size(X_INDEX) * 2
#ifdef U_ACTIVE
    allocate(current_state%su%data(z_size, y_size, x_size))
    current_state%su%data(:,:,:) = 0.
    allocate(current_state%savu%data(z_size, y_size, x_size))
    current_state%savu%data(:,:,:) = 0.
#endif
#ifdef V_ACTIVE
    allocate(current_state%sv%data(z_size, y_size, x_size))
    current_state%sv%data(:,:,:) = 0.
    allocate(current_state%savv%data(z_size, y_size, x_size))
    current_state%savv%data(:,:,:) = 0.
#endif
#ifdef W_ACTIVE
    allocate(current_state%sw%data(z_size, y_size, x_size))
    current_state%sw%data(:,:,:) = 0.
    allocate(current_state%savw%data(z_size, y_size, x_size))
    current_state%savw%data(:,:,:) = 0.
#endif
    if (current_state%th%active) then
      allocate(current_state%sth%data(z_size, y_size, x_size))
      current_state%sth%data(:,:,:) = 0.
    end if
    do i=1,current_state%number_q_fields
      allocate(current_state%sq(i)%data(z_size, y_size, x_size))
      current_state%sq(i)%data(:,:,:) = 0.
    end do
  end subroutine initalise_source_and_sav_fields

  !> Loads in misc data from the checkpoint file
  !! @param current_state The current model state_mod
  !! @param ncid The NetCDF file id
  subroutine load_misc(current_state, ncid)
    type(model_state_type), intent(inout) :: current_state
    integer, intent(in) :: ncid

    integer :: i_data(1)  ! Procedure requires a vector rather than scalar
    real(kind=DEFAULT_PRECISION) :: r_data(1)

    call read_single_variable(ncid, TIMESTEP, integer_data_1d=i_data)
    current_state%timestep = i_data(1)+1 ! plus one to increment for next timestep
    !current_state%start_timestep = current_state%timestep+1    
    call read_single_variable(ncid, UGAL, real_data_1d_double=r_data)
    current_state%ugal = r_data(1)
    call read_single_variable(ncid, VGAL, real_data_1d_double=r_data)
    current_state%vgal = r_data(1)
    call read_single_variable(ncid, NQFIELDS, integer_data_1d=i_data)
    current_state%number_q_fields = i_data(1)
    call read_single_variable(ncid, DTM_KEY, real_data_1d_double=r_data)
    current_state%dtm = r_data(1)
    call read_single_variable(ncid, DTM_NEW_KEY, real_data_1d_double=r_data)
    current_state%dtm_new = r_data(1)
    current_state%update_dtm = current_state%dtm .ne. current_state%dtm_new
    call read_single_variable(ncid, ABSOLUTE_NEW_DTM_KEY, real_data_1d_double=r_data)
    current_state%absolute_new_dtm = r_data(1)
    call read_single_variable(ncid, TIME_KEY, real_data_1d_double=r_data)
    ! The time is written into checkpoint as time+dtm, therefore the time as read in has been correctly advanced
    current_state%time = r_data(1)
  end subroutine load_misc

  !> Will read a global attribute from the checkpoint file - note that it allocates string memory
  !! here so the caller should deallocate this to avoid memory leaks
  !! @param ncid The NetCDF file id
  !! @param key The NetCDF global attributes key to look up
  function read_specific_global_attribute(ncid, key)
    integer, intent(in) :: ncid
    character(len=*), intent(in) :: key

    integer :: length
    character(len=:),allocatable,target :: read_specific_global_attribute

    call check_status(nf90_inquire_attribute(ncid, nf90_global, key, len = length))
    allocate(character(length) :: read_specific_global_attribute)
    call check_status(nf90_get_att(ncid, nf90_global, key, read_specific_global_attribute))
  end function read_specific_global_attribute

  !> Reads in and initialises all prognostic fields. It will check the NetCDF checkpoint file
  !! and depending upon what data is in the file then it will set the corresponding fields up
  !! which in essence determines whether it has 1, 2 or 3D fields. This also supports 1, 2 or 3D
  !! grids_mod which will set any missing grid dimensions to be 1
  !! @param current_state The current model state_mod
  !! @param ncid The NetCDF file id
  !! @param z_dim The size of the z dimension
  !! @param y_dim The size of the y dimension
  !! @param x_dim The size of the x dimension
  subroutine load_all_fields(current_state, ncid)
    type(model_state_type), intent(inout) :: current_state
    integer, intent(in) :: ncid

    integer :: i
    logical :: multi_process

    multi_process = current_state%parallel%processes .gt. 1

    if (does_field_exist(ncid, U_KEY)) then
      call load_single_3d_field(ncid, current_state%local_grid, current_state%u, DUAL_GRID, &
           DUAL_GRID, PRIMAL_GRID, U_KEY, multi_process)
      if (current_state%ugal .ne. 0.0) current_state%u%data=current_state%u%data-current_state%ugal
      call load_single_3d_field(ncid, current_state%local_grid, current_state%zu, DUAL_GRID, &
           DUAL_GRID, PRIMAL_GRID, ZU_KEY, multi_process)
    end if
    if (does_field_exist(ncid, V_KEY)) then
      call load_single_3d_field(ncid, current_state%local_grid, current_state%v, DUAL_GRID, &
           PRIMAL_GRID, DUAL_GRID, V_KEY, multi_process)
      if (current_state%vgal .ne. 0.0) current_state%v%data=current_state%v%data-current_state%vgal
      call load_single_3d_field(ncid, current_state%local_grid, current_state%zv, DUAL_GRID, &
           PRIMAL_GRID, DUAL_GRID, ZV_KEY, multi_process)
    end if
    if (does_field_exist(ncid, W_KEY)) then
      call load_single_3d_field(ncid, current_state%local_grid, current_state%w, PRIMAL_GRID, &
           DUAL_GRID, DUAL_GRID, W_KEY, multi_process)
      call load_single_3d_field(ncid, current_state%local_grid, current_state%zw, PRIMAL_GRID, &
           DUAL_GRID, DUAL_GRID, ZW_KEY, multi_process)
    end if
    if (does_field_exist(ncid, TH_KEY)) then
      call load_single_3d_field(ncid, current_state%local_grid, current_state%th, DUAL_GRID, &
           DUAL_GRID, DUAL_GRID, TH_KEY, multi_process)
      call load_single_3d_field(ncid, current_state%local_grid, current_state%zth, DUAL_GRID, &
           DUAL_GRID, DUAL_GRID, ZTH_KEY, multi_process)
    end if
    if (does_field_exist(ncid, P_KEY)) then
      call load_single_3d_field(ncid, current_state%local_grid, current_state%p, DUAL_GRID, &
           DUAL_GRID, DUAL_GRID, P_KEY, multi_process)
    end if    
    allocate(current_state%q(current_state%number_q_fields), current_state%zq(current_state%number_q_fields), &
         current_state%sq(current_state%number_q_fields))
    do i=1,current_state%number_q_fields
      call load_single_3d_field(ncid, current_state%local_grid, current_state%q(i), DUAL_GRID, &
           DUAL_GRID, DUAL_GRID, Q_KEY, multi_process, i)
      call load_single_3d_field(ncid, current_state%local_grid, current_state%zq(i), DUAL_GRID, &
           DUAL_GRID, DUAL_GRID, ZQ_KEY, multi_process, i)
    end do
  end subroutine load_all_fields

  !> Determines whether a variable (field) exists within the NetCDF checkpoint file
  !! @ncid The NetCDF file id
  !! @variable_key The NetCDF variable key we are looking up
  logical function does_field_exist(ncid, variable_key)
    integer, intent(in) :: ncid
    character(len=*), intent(in) :: variable_key

    integer :: variable_id

    call check_status(nf90_inq_varid(ncid, variable_key, variable_id), does_field_exist)
  end function does_field_exist

  !> Loads the Q indices (if they exist) into the model Q index register
  !! @param ncid The NetCDF file id
  subroutine load_q_indices(ncid)
    integer, intent(in) :: ncid

    integer :: number_q_indices, i, q_indices_id
    character(len=STRING_LENGTH) :: key, value

    number_q_indices=get_number_q_indices(ncid)
    if (number_q_indices .gt. 0) then
      call check_status(nf90_inq_varid(ncid, Q_INDICES_KEY, q_indices_id))
      do i=1, number_q_indices
        call check_status(nf90_get_var(ncid, q_indices_id, key, (/ 1, 1, i /)))
        call check_status(nf90_get_var(ncid, q_indices_id, value, (/ 1, 2, i /)))
        call remove_null_terminator_from_string(key)
        call remove_null_terminator_from_string(value)
        call set_q_index(conv_to_integer(trim(value)), key)
      end do
    end if
  end subroutine load_q_indices

  !> Retrieve sthe number of Q indice entries in the NetCDF. This is optional, so may return 0
  !! @param ncid The NetCDF file id
  integer function get_number_q_indices(ncid)
    integer, intent(in) :: ncid

    integer :: q_indices_dimid, q_indices_dim
    logical :: found_flag

    call check_status(nf90_inq_dimid(ncid, Q_INDICES_KEY, q_indices_dimid), found_flag)
    if (found_flag) then
      call check_status(nf90_inquire_dimension(ncid, q_indices_dimid, len=q_indices_dim))
      get_number_q_indices=q_indices_dim
    else
      get_number_q_indices=0
    end if
  end function get_number_q_indices  

  !> Loads in and initialises a single 3D prognostic field. Even if the model is being run in 1 or 2D, the fields
  !! are still stored in 3D with the missing dimension sizes being set to one.
  !! @param ncid The NetCDF file id
  !! @param field The prognostic field that will be initialised from the checkpoint
  !! @param grid The grid that the prognostic field is on
  !! @param variable_key The NetCDF variable name
  !! @param dim_one Size in dimension one
  !! @param dim_two Size in dimension two
  !! @param dim_three Size in dimension three
  subroutine load_single_3d_field(ncid, local_grid, field, z_grid, y_grid, x_grid, variable_key, multi_process, fourth_dim_loc)
    character(len=*), intent(in) :: variable_key
    integer, intent(in) :: z_grid, y_grid, x_grid, ncid
    type(prognostic_field_type), intent(inout) :: field
    type(local_grid_type), intent(inout) :: local_grid
    logical, intent(in) :: multi_process
    integer, optional, intent(in) :: fourth_dim_loc

    integer :: start(4), count(4), i, map(4)

    if (allocated(field%data)) deallocate(field%data)
    allocate(field%data(local_grid%size(Z_INDEX) + local_grid%halo_size(Z_INDEX) * 2, local_grid%size(Y_INDEX) + &
         local_grid%halo_size(Y_INDEX) * 2, local_grid%size(X_INDEX) + local_grid%halo_size(X_INDEX) * 2))
    field%data(:,:,:)=0.

    if (multi_process .or. present(fourth_dim_loc)) then
      do i=1,3
        if (i==1) then
          map(i)=1
        else
          map(i)=map(i-1)*local_grid%size(i-1)
        end if
        start(i) = local_grid%start(i)
        count(i) = local_grid%size(i)
      end do
      if (present(fourth_dim_loc)) then
        start(4) = fourth_dim_loc
        count(4) = 1
        map(4)=map(3)*local_grid%size(3)
      end if

      call read_single_variable(ncid, variable_key, real_data_3d=field%data(local_grid%local_domain_start_index(Z_INDEX):&
           local_grid%local_domain_end_index(Z_INDEX),local_grid%local_domain_start_index(Y_INDEX):&
           local_grid%local_domain_end_index(Y_INDEX), local_grid%local_domain_start_index(X_INDEX):&
           local_grid%local_domain_end_index(X_INDEX)), start=start, count=count, map=map)
    else
      call read_single_variable(ncid, variable_key, real_data_3d=field%data(local_grid%local_domain_start_index(Z_INDEX):&
           local_grid%local_domain_end_index(Z_INDEX),local_grid%local_domain_start_index(Y_INDEX):&
           local_grid%local_domain_end_index(Y_INDEX), local_grid%local_domain_start_index(X_INDEX):&
           local_grid%local_domain_end_index(X_INDEX)))
    end if

    field%grid(Z_INDEX) = z_grid
    field%grid(Y_INDEX) = y_grid
    field%grid(X_INDEX) = x_grid
    field%active = .true.
  end subroutine load_single_3d_field

  !> Reads in and initialises the primal grid from the NetCDF checkpoint and will handle 1, 2 or 3D runs
  !! @param current_state The current model state_mod
  !! @param ncid The NetCDF file id
  !! @param z_dim Size in the z dimension
  !! @param y_dim Size in the y dimension
  !! @param x_sim Size in the x dimension
  !! @param z_found Whether the z grid dimension exists in the checkpoint
  !! @param y_found Whether the y grid dimension exists in the checkpoint
  !! @param x_found Whether the x grid dimension exists in the checkpoint
  subroutine load_global_grid(current_state, ncid, z_dim, y_dim, x_dim, z_found, y_found, x_found)
    type(model_state_type), intent(inout) :: current_state
    integer, intent(in) :: ncid, z_dim, y_dim, x_dim
    logical, intent(in) :: z_found, y_found, x_found

    current_state%global_grid%dimensions = 0

    if (z_found) then
      call read_dimension_of_grid(ncid, current_state%global_grid, Z_KEY, Z_INDEX, z_dim)
      call define_vertical_levels(ncid, current_state, Z_KEY, z_dim)
    end if
    if (y_found) call read_dimension_of_grid(ncid, current_state%global_grid, Y_KEY, Y_INDEX, y_dim)
    if (x_found) call read_dimension_of_grid(ncid, current_state%global_grid, X_KEY, X_INDEX, x_dim)
  end subroutine load_global_grid

  !> Defines the vertical levels of the grid. This is both the grid points and corresponding height
  !! for each point in metres
  !! @param current_state The current model state_mod
  !! @param z_key Key of the Z grid points in the NetCDF
  !! @param z_size Number of grid points
  subroutine define_vertical_levels(ncid, current_state, z_key, z_size)
    type(model_state_type), intent(inout) :: current_state
    character(len=*), intent(in) :: z_key
    integer, intent(in) :: z_size, ncid

    integer :: i
    real, dimension(:), allocatable :: data

    allocate(data(z_size))
    call read_single_variable(ncid, z_key, real_data_1d=data)

    allocate(current_state%global_grid%configuration%vertical%kgd(z_size), &
         current_state%global_grid%configuration%vertical%hgd(z_size))

    do i=1,z_size
      current_state%global_grid%configuration%vertical%kgd(i) = i
      current_state%global_grid%configuration%vertical%hgd(i) = real(data(i))
    end do
    deallocate(data)
  end subroutine define_vertical_levels

  !> Reads a specific dimension of the grid from the NetCDF and sets this up in the state_mod
  !! @param ncid The NetCDF file id
  !! @param grid The model grid that this dimension lies on
  !! @param variable_key The NetCDF variable name that we are reading
  !! @param dimension The dimension corresponding to the definition in the grids_mod module
  !! @param dimension_size Number of grid points in this dimension
  subroutine read_dimension_of_grid(ncid, grid, variable_key, dimension, dimension_size)
    integer, intent(in) :: ncid, dimension_size, dimension
    character(len=*), intent(in) :: variable_key
    type(global_grid_type), intent(inout) :: grid

    real, dimension(:), allocatable :: data

    allocate(data(dimension_size))
    call read_single_variable(ncid, variable_key, real_data_1d=data)
    grid%bottom(dimension) = int(data(1))
    grid%top(dimension) = int(data(dimension_size))
    grid%resolution(dimension) = int(data(2) - data(1))
    grid%size(dimension) = dimension_size
    grid%dimensions = grid%dimensions + 1
    grid%active(dimension) = .true.
    deallocate(data)
  end subroutine read_dimension_of_grid

  !> Reads in a single variable and sets the values of the data based upon this. Handles reading in
  !! and setting vectors of 1 or 3D reals and 1d integers
  !! @param ncid The NetCDF file id
  !! @param key The variable key to read in
  !! @param realData1D Vector of 1D real data to read (optional)
  !! @param realData3D Vector of 3D real data to read (optional)
  !! @param integerData1D Vector of 1D integer data to read (optional)
  subroutine read_single_variable(ncid, key, real_data_1d, real_data_1d_double, real_data_3d, integer_data_1d, start, count, map)
    integer, intent(in) :: ncid
    character(len=*), intent(in) :: key
    real, dimension(:), intent(inout), optional :: real_data_1d
    real(kind=DEFAULT_PRECISION), dimension(:), intent(inout), optional :: real_data_1d_double
    real(kind=DEFAULT_PRECISION), dimension(:,:,:), intent(inout), optional :: real_data_3d
    integer, dimension(:), intent(inout), optional :: integer_data_1d
    integer, dimension(:), intent(in), optional :: start, count, map

    integer :: variable_id

    call check_status(nf90_inq_varid(ncid, key, variable_id))

    if (.not. present(real_data_1d) .and. .not. present(real_data_1d_double) .and. &
         .not. present(real_data_3d) .and. .not. present(integer_data_1d)) return

    if (present(real_data_1d)) then
      if (present(start) .and. present(count) .and. present(map)) then
        call check_status(nf90_get_var(ncid, variable_id, real_data_1d, start=start, count=count, map=map))
      else
        call check_status(nf90_get_var(ncid, variable_id, real_data_1d))
      end if
    else if (present(real_data_1d_double)) then
      if (present(start) .and. present(count) .and. present(map)) then
        call check_status(nf90_get_var(ncid, variable_id, real_data_1d_double, start=start, count=count, map=map))
      else
        call check_status(nf90_get_var(ncid, variable_id, real_data_1d_double))
      end if
    else if (present(real_data_3d)) then
      if (present(start) .and. present(count) .and. present(map)) then
        call check_status(nf90_get_var(ncid, variable_id, real_data_3d, start=start, count=count, map=map))
      else
        call check_status(nf90_get_var(ncid, variable_id, real_data_3d))
      end if
    else if (present(integer_data_1d)) then
      if (present(start) .and. present(count) .and. present(map)) then
        call check_status(nf90_get_var(ncid, variable_id, integer_data_1d, start, count, map))
      else
        call check_status(nf90_get_var(ncid, variable_id, integer_data_1d))
      end if
    end if
  end subroutine read_single_variable

  !> Will read in the dimension sizes from the NetCDF file along with whether the z, y and x dimensions have been found.
  !! If any have not been found then this corresponds to running the model with 2 or 1D grids_mod.
  !! @param ncid The NetCDf file id
  !! @param z_dim Size in the z dimension (set to size of empty dimension if not found)
  !! @param y_dim Size in the y dimension (set to size of empty dimension if not found)
  !! @param x_dim Size in the x dimension (set to size of empty dimension if not found)
  !! @param string_dim Size of the string dimension (size that we use to allocate strings)
  !! @param key_value_pair_dim Size of key-value pair dimension which should be two
  !! @param options_dim Size of the options dimension which corresponds to how many options have been provided
  !! @param z_found Whether the z dimension exists (grid exists in z dimension)
  !! @param y_found Whether the y dimension exists (grid exists in y dimension)
  !! @param x_found Whether the x dimension exists (grid exists in x dimension)
  subroutine read_dimensions(ncid, z_dim, y_dim, x_dim, z_found, y_found, x_found)
    integer, intent(in) :: ncid
    integer, intent(out) :: z_dim, y_dim, x_dim
    logical, intent(out) :: z_found, y_found, x_found

    integer ::  z_dimid, y_dimid, x_dimid, empty_dimid, empty_dim

    call check_status(nf90_inq_dimid(ncid, Z_DIM_KEY, z_dimid), z_found)
    call check_status(nf90_inq_dimid(ncid, Y_DIM_KEY, y_dimid), y_found)
    call check_status(nf90_inq_dimid(ncid, X_DIM_KEY, x_dimid), x_found)
    call check_status(nf90_inq_dimid(ncid, EMPTY_DIM_KEY, empty_dimid))

    call check_status(nf90_inquire_dimension(ncid, empty_dimid, len=empty_dim))

    if (z_found) then
      call check_status(nf90_inquire_dimension(ncid, z_dimid, len=z_dim))
    else
      z_dim = empty_dim
    end if
    if (y_found) then
      call check_status(nf90_inquire_dimension(ncid, y_dimid, len=y_dim))
    else
      y_dim = empty_dim
    end if
    if (x_found) then
      call check_status(nf90_inquire_dimension(ncid, x_dimid, len=x_dim))
    else
      x_dim = empty_dim
    end if
  end subroutine read_dimensions
end module checkpointer_read_checkpoint_mod