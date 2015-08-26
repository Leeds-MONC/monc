!> Writes out model state_mod to a checkpoint NetCDF file
module checkpointer_write_checkpoint_mod
#ifndef TEST_MODE
  use netcdf, only : NF90_DOUBLE, NF90_REAL, NF90_INT, NF90_CHAR, NF90_GLOBAL, NF90_CLOBBER, NF90_NETCDF4, NF90_MPIIO, &
       NF90_COLLECTIVE, nf90_def_var, nf90_var_par_access, nf90_def_var_fill, nf90_put_att, nf90_create, nf90_put_var, &
       nf90_def_dim, nf90_enddef, nf90_close
#else
  use dummy_netcdf_mod, only : NF90_DOUBLE, NF90_REAL, NF90_INT, NF90_CHAR, NF90_GLOBAL, NF90_CLOBBER, NF90_NETCDF4, NF90_MPIIO, &
       nf90_def_var, nf90_put_att, nf90_create, nf90_put_var, nf90_def_dim, nf90_enddef, nf90_close
#endif
  use state_mod, only : model_state_type
  use grids_mod, only : local_grid_type, global_grid_type, vertical_grid_configuration_type, X_INDEX, Y_INDEX, Z_INDEX
  use prognostics_mod, only : prognostic_field_type
  use conversions_mod, only : conv_to_string
  use optionsdatabase_mod, only : options_size, options_value_at, options_key_at
  use checkpointer_common_mod, only : EMPTY_DIM_KEY, KEY_VALUE_PAIR_KEY, OPTIONS_DIM_KEY, OPTIONS_KEY, STRING_DIM_KEY, &
       X_DIM_KEY, Y_DIM_KEY, Z_DIM_KEY, Q_DIM_KEY, Q_KEY, ZQ_KEY, TH_KEY, ZTH_KEY, P_KEY, U_KEY, V_KEY, W_KEY, ZU_KEY, ZV_KEY, &
       ZW_KEY, X_KEY, Y_KEY, Z_KEY, NQFIELDS, UGAL, VGAL, TIME_KEY, TIMESTEP, MAX_STRING_LENGTH, CREATED_ATTRIBUTE_KEY, &
       TITLE_ATTRIBUTE_KEY, ABSOLUTE_NEW_DTM_KEY, DTM_KEY, DTM_NEW_KEY, Q_INDICES_KEY, check_status
  use datadefn_mod, only : DEFAULT_PRECISION, DOUBLE_PRECISION, STRING_LENGTH
  use q_indices_mod, only : q_metadata_type, get_max_number_q_indices, get_indices_descriptor, get_number_active_q_indices
  use mpi, only : MPI_INFO_NULL
  implicit none

#ifndef TEST_MODE
  private
#endif

  character(len=*), parameter :: CHECKPOINT_TITLE = "MONC checkpoint file" !< Title of the NetCDF file

  public write_checkpoint_file

contains

  !> Will write out the current model state_mod into a NetCDF checkpoint file
  !! @param currentState The current model state_mod
  !! @param filename The filename of the NetCDF file that will be written
  subroutine write_checkpoint_file(current_state, filename)
    type(model_state_type), intent(inout) :: current_state
    character(len=*), intent(in) :: filename

    integer :: ncid, z_dim_id, y_dim_id, x_dim_id, q_dim_id, x_id, y_id, z_id, th_id, p_id, time_id,&
         u_id, v_id, w_id, q_id, zu_id, zv_id, zw_id, zth_id, zq_id, timestep_id, ugal_id, &
         vgal_id, number_q_fields_id, string_dim_id, key_value_dim_id, options_id, q_indices_id, &
         dtm_id, dtm_new_id, absolute_new_dtm_id
    logical :: q_indices_declared

    ! If there are multiple processes then open for parallel IO
    if (current_state%parallel%processes .gt. 1) then
      call check_status(nf90_create(filename, ior(NF90_NETCDF4, NF90_MPIIO), ncid, &
           comm = current_state%parallel%monc_communicator, info = MPI_INFO_NULL))
    else
      call check_status(nf90_create(filename, NF90_CLOBBER, ncid))
    end if
    call write_out_global_attributes(ncid)
    call define_grid_dimensions(current_state, ncid, z_dim_id, y_dim_id, x_dim_id)
    if (current_state%number_q_fields .gt. 0) call define_q_field_dimension(current_state, ncid, q_dim_id)

    call define_options_variable(current_state, ncid, string_dim_id, key_value_dim_id, options_id)
    q_indices_declared=define_q_indices_variable(ncid, string_dim_id, key_value_dim_id, q_indices_id)
    call define_grid_variables(current_state, ncid, z_dim_id, y_dim_id, x_dim_id, x_id, y_id, z_id)
    if (current_state%number_q_fields .gt. 0) call define_q_variable(ncid, current_state%parallel%processes .gt. 1, &
         q_dim_id, z_dim_id, y_dim_id, x_dim_id, q_id, zq_id)
    call define_prognostic_variables(current_state, current_state%parallel%processes .gt. 1, ncid, z_dim_id, y_dim_id, &
         x_dim_id, u_id, v_id, w_id, th_id, p_id, zu_id, zv_id, zw_id, zth_id)
    call define_misc_variables(ncid, timestep_id, time_id, ugal_id, vgal_id, number_q_fields_id, &
         dtm_id, dtm_new_id, absolute_new_dtm_id)

    call check_status(nf90_enddef(ncid))

    if (current_state%parallel%my_rank==0) call write_out_grid(ncid, current_state%global_grid, z_id, y_id, x_id)
    call write_out_all_fields(current_state, ncid, u_id, v_id, w_id, zu_id, zv_id, zw_id, th_id, zth_id, q_id, zq_id, p_id)
    if (current_state%parallel%my_rank==0) then
      call write_out_options(current_state, ncid, options_id)
      if (q_indices_declared) call write_out_q_indices(ncid, q_indices_id)
      call write_out_misc_variables(current_state, ncid, timestep_id, time_id, &
           ugal_id, vgal_id, number_q_fields_id, dtm_id, dtm_new_id, absolute_new_dtm_id)
    end if

    call check_status(nf90_close(ncid))
  end subroutine write_checkpoint_file

  !> Writes out global attributes into the checkpoint
  !! @param ncid NetCDF file id
  subroutine write_out_global_attributes(ncid)
    integer, intent(in) :: ncid

    integer :: date_values(8)

    call date_and_time(values=date_values)

    call check_status(nf90_put_att(ncid, NF90_GLOBAL, TITLE_ATTRIBUTE_KEY, CHECKPOINT_TITLE))
    call check_status(nf90_put_att(ncid, NF90_GLOBAL, CREATED_ATTRIBUTE_KEY, trim(conv_to_string(date_values(3)))//"/"//&
         trim(conv_to_string(date_values(2)))//"/"//trim(conv_to_string(date_values(1)))//" "//trim(conv_to_string(&
         date_values(5)))// ":"//trim(conv_to_string(date_values(6)))//":"//trim(conv_to_string(date_values(7)))))
  end subroutine write_out_global_attributes

  !> Writes out the specific Q indicies that are active and need writing
  !! @param ncid The NetCDF file id
  !! @param q_indicies_id The Q indicies NetCDF variable id
  subroutine write_out_q_indices(ncid, q_indices_id)
    integer, intent(in) :: ncid, q_indices_id

    integer :: i, current_index
    type(q_metadata_type) :: specific_q_data

    current_index=1
    do i=1, get_max_number_q_indices()
      specific_q_data=get_indices_descriptor(i)
      if (specific_q_data%l_used) then
        call check_status(nf90_put_var(ncid, q_indices_id, trim(specific_q_data%name), (/ 1, 1, current_index /)))
        call check_status(nf90_put_var(ncid, q_indices_id, trim(conv_to_string(i)), (/ 1, 2, current_index /)))
        current_index=current_index+1
      end if
    end do    
  end subroutine write_out_q_indices

  !> Writes out the options that the model was run with
  !! @param current_state The current model state_mod
  !! @param ncid The NetCDF file id
  !! @param options_id The options NetCDF variable id
  subroutine write_out_options(current_state, ncid, options_id)
    type(model_state_type), intent(inout) :: current_state
    integer, intent(in) :: ncid, options_id

    integer :: i
    character(len=STRING_LENGTH), pointer :: sized_raw_character
    class(*), pointer :: raw_data, raw_to_string

    do i=1,options_size(current_state%options_database)
      raw_data=> options_value_at(current_state%options_database, i)
      raw_to_string=>raw_data
      call check_status(nf90_put_var(ncid, options_id, trim(options_key_at(current_state%options_database, i)), (/ 1, 1, i /)))
      select type (raw_data)
      type is(integer)
        call check_status(nf90_put_var(ncid, options_id, trim(conv_to_string(raw_data)), (/ 1, 2, i /)))
      type is(real)
        call check_status(nf90_put_var(ncid, options_id, trim(conv_to_string(raw_data)), (/ 1, 2, i /)))
      type is(logical)
        call check_status(nf90_put_var(ncid, options_id, trim(conv_to_string(raw_data)), (/ 1, 2, i /)))
      type is(character(len=*))
        ! Done this way to give the character size information and keep the (unsafe) cast in the conversion module
        sized_raw_character=>conv_to_string(raw_to_string, .false., STRING_LENGTH)
        call check_status(nf90_put_var(ncid, options_id, trim(sized_raw_character), (/ 1, 2, i /)))
      end select
    end do
  end subroutine write_out_options

  !> Will write out all prognostic model fields to the checkpoint. It will work in 1, 2 or 3D depending on the model
  !! @param current_state The current model state_mod
  !! @param ncid The NetCDF file id
  !! @param u_id The NetCDF u field dimension id
  !! @param v_id The NetCDF v field dimension id
  !! @param w_id The NetCDF w field dimension id
  subroutine write_out_all_fields(current_state, ncid, u_id, v_id, w_id, zu_id, zv_id, zw_id, th_id, zth_id, q_id, zq_id, p_id)
    type(model_state_type), intent(inout) :: current_state
    integer, intent(in) :: ncid, u_id, v_id, w_id, zu_id, zv_id, zw_id, th_id, zth_id, q_id, zq_id, p_id

    integer :: i
    logical :: multi_process

    multi_process = current_state%parallel%processes .gt. 1
#ifdef U_ACTIVE
    ! Add on the Galilean transformation to the U field for dumping
    if (current_state%ugal .ne. 0.0) current_state%u%data = current_state%u%data + current_state%ugal
    call write_out_velocity_field(ncid, current_state%local_grid, current_state%u, u_id, multi_process)
    if (current_state%ugal .ne. 0.0) current_state%u%data = current_state%u%data - current_state%ugal
    call write_out_velocity_field(ncid, current_state%local_grid, current_state%zu, zu_id, multi_process)
#endif
#ifdef V_ACTIVE
    ! Add on the Galilean transformation to the V field for dumping
    if (current_state%vgal .ne. 0.0) current_state%v%data = current_state%v%data + current_state%vgal
    call write_out_velocity_field(ncid, current_state%local_grid, current_state%v, v_id, multi_process)
    if (current_state%vgal .ne. 0.0) current_state%v%data = current_state%v%data - current_state%vgal
    call write_out_velocity_field(ncid, current_state%local_grid, current_state%zv, zv_id, multi_process)
#endif
#ifdef W_ACTIVE
    call write_out_velocity_field(ncid, current_state%local_grid, current_state%w, w_id, multi_process)
    call write_out_velocity_field(ncid, current_state%local_grid, current_state%zw, zw_id, multi_process)
#endif
    if (current_state%th%active) then
      call write_out_velocity_field(ncid, current_state%local_grid, current_state%th, th_id, multi_process)
      call write_out_velocity_field(ncid, current_state%local_grid, current_state%zth, zth_id, multi_process)
    end if
    if (current_state%p%active) call write_out_velocity_field(ncid, current_state%local_grid, current_state%p, &
         p_id, multi_process)
    do i=1,current_state%number_q_fields
      if (current_state%q(i)%active) then
        call write_out_velocity_field(ncid, current_state%local_grid, current_state%q(i), q_id, multi_process, i)
        call write_out_velocity_field(ncid, current_state%local_grid, current_state%zq(i), zq_id, multi_process, i)
      end if
    end do
  end subroutine write_out_all_fields

  !> Will write out a single velocity field to the checkpoint file. If there are multiple processes then will determine
  !! the bounds, otherwise for serial just dump data field
  !! @param ncid The NetCDF file id
  !! @param field The model prognostic field to write out
  !! @param variable_id The NetCDF variable dimension id
  subroutine write_out_velocity_field(ncid, local_grid, field, variable_id, multi_process, fourth_dim_loc)
    integer, intent(in) :: ncid, variable_id
    type(prognostic_field_type), intent(in) :: field
    type(local_grid_type), intent(inout) :: local_grid
    logical, intent(in) :: multi_process
    integer, optional, intent(in) :: fourth_dim_loc

    integer :: start(4), count(4), i, map(4)

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

      call check_status(nf90_put_var(ncid, variable_id, field%data(local_grid%local_domain_start_index(Z_INDEX):&
           local_grid%local_domain_end_index(Z_INDEX),local_grid%local_domain_start_index(Y_INDEX):&
           local_grid%local_domain_end_index(Y_INDEX), local_grid%local_domain_start_index(X_INDEX):&
           local_grid%local_domain_end_index(X_INDEX)), start=start, count=count))
    else
      call check_status(nf90_put_var(ncid, variable_id, field%data(local_grid%local_domain_start_index(Z_INDEX):&
           local_grid%local_domain_end_index(Z_INDEX),local_grid%local_domain_start_index(Y_INDEX):&
           local_grid%local_domain_end_index(Y_INDEX), local_grid%local_domain_start_index(X_INDEX):&
           local_grid%local_domain_end_index(X_INDEX))))
    end if
  end subroutine write_out_velocity_field

  !> Will write out the grid to the checkpoint, it will work in 1, 2 or 3D depending on what is in the model
  !! @param ncid The NetCDF file id
  !! @param grid The model grid to write out
  !! @param z_id The NetCDF z variable id
  !! @param y_id The NetCDF y variable id
  !! @param x_id The NetCDF x variable id
  subroutine write_out_grid(ncid, grid, z_id, y_id, x_id)
    integer, intent(in) :: ncid, z_id, y_id, x_id
    type(global_grid_type), intent(in) :: grid

    if (grid%active(Z_INDEX)) then
      call write_z_grid_gimension(ncid, grid%configuration%vertical, z_id)
    end if
    if (grid%active(Y_INDEX)) then
      call write_grid_dimension(ncid, grid%size(Y_INDEX), grid%bottom(Y_INDEX), grid%resolution(Y_INDEX), y_id)
    end if
    if (grid%active(X_INDEX)) then
      call write_grid_dimension(ncid, grid%size(X_INDEX), grid%bottom(X_INDEX), grid%resolution(X_INDEX), x_id)
    end if
  end subroutine write_out_grid

  !> Writes out the Z dimension of the grids_mod points which are explicitly calculated
  !! @param ncid The NetCDF file id
  !! @param verticalGrid The vertical grid configuration
  !! @param z_id The NetCDF id of the z field
  subroutine write_z_grid_gimension(ncid, vertical_grid, z_id)
    type(vertical_grid_configuration_type), intent(in) :: vertical_grid
    integer, intent(in) :: ncid, z_id

    call check_status(nf90_put_var(ncid, z_id, vertical_grid%z))
  end subroutine write_z_grid_gimension

  !> Will write out a specific dimension of a grid to the checkpoint file. In the checkpoint we are currently
  !! dumping out the entire grid which is reconstructed from the model parameters (probably want to change this
  !! and write out only the parameters.)
  !! @param ncid The NetCDF file id
  !! @param grids_modize The size of the grid (number of points)
  !! @param gridbottom The smallest point in the grid
  !! @param gridresolution The resolution of the grid
  !! @param dimension_id The NetCDF id of the grid variable to write to
  subroutine write_grid_dimension(ncid, grids_modize, gridbottom, gridresolution, dimension_id)
    integer, intent(in) :: ncid, grids_modize, gridbottom, gridresolution, dimension_id

    real, dimension(:), allocatable :: data
    integer :: i

    allocate(data(grids_modize))
    do i=1,grids_modize
      data(i) = gridbottom + gridresolution * i
    end do
    call check_status(nf90_put_var(ncid, dimension_id, data))
    deallocate(data)
  end subroutine write_grid_dimension

  !> Defines the NetCDF options variable which is basically a 3D character array to form key-value
  !! pair strings for each entry
  !! @param current_state The current model state_mod
  !! @param ncid The NetCDF file id
  !! @param options_id The NetCDF options variable id that is created in this procedure
  subroutine define_options_variable(current_state, ncid, string_dim_id, key_value_dim_id, options_id)
    type(model_state_type), intent(inout) :: current_state
    integer, intent(in) :: ncid
    integer, intent(out) :: string_dim_id, key_value_dim_id, options_id

    integer :: options_dim_id, command_dimensions(3)

    call check_status(nf90_def_dim(ncid, STRING_DIM_KEY, MAX_STRING_LENGTH, string_dim_id))
    call check_status(nf90_def_dim(ncid, KEY_VALUE_PAIR_KEY, 2, key_value_dim_id))
    call check_status(nf90_def_dim(ncid, OPTIONS_DIM_KEY, options_size(current_state%options_database), options_dim_id))

    command_dimensions = (/ string_dim_id, key_value_dim_id, options_dim_id /)

    call check_status(nf90_def_var(ncid, OPTIONS_KEY, NF90_CHAR, command_dimensions, options_id))
  end subroutine define_options_variable

  !> Defines the NetCDF Q indices variable which is, same as the options, stored as key-value pair of strings.
  !! This will only store and create these dimensions if there are any Q indicies to store
  !! @param ncid The NetCDF file id
  !! @param string_dim_id The id of the string dimension
  !! @param key_value_dim_id The id of the key value pair dimension
  !! @param q_indices_id The id that represents this variable in the NetCDF
  !! @returns Whether or not the dimensions and variable has been declared (it won't be if there are no active Q indicies)
  logical function define_q_indices_variable(ncid, string_dim_id, key_value_dim_id, q_indices_id)
    integer, intent(in) :: ncid, string_dim_id, key_value_dim_id
    integer, intent(out) :: q_indices_id

    integer :: q_indices_dim_id, command_dimensions(3), number_active_q

    number_active_q=get_number_active_q_indices()

    if (number_active_q == 0) then
      define_q_indices_variable=.false.
    else
      call check_status(nf90_def_dim(ncid, Q_INDICES_KEY, number_active_q, q_indices_dim_id))

      command_dimensions = (/ string_dim_id, key_value_dim_id, q_indices_dim_id /)

      call check_status(nf90_def_var(ncid, Q_INDICES_KEY, NF90_CHAR, command_dimensions, q_indices_id))
      define_q_indices_variable=.true.
    end if
  end function define_q_indices_variable  

  !> Defines the Q field dimension in the NetCDF
  !! @param current_state The current model state_mod
  !! @param ncid The NetCDF file id
  !! @param q_dim_id Corresponding NetCDF Q dimension id
  subroutine define_q_field_dimension(current_state, ncid, q_dim_id)
    type(model_state_type), intent(inout) :: current_state
    integer, intent(in) :: ncid
    integer, intent(out) :: q_dim_id

    call check_status(nf90_def_dim(ncid, Q_DIM_KEY, current_state%number_q_fields, q_dim_id))
  end subroutine define_q_field_dimension

  !> Will define the grid dimensions and works for 1, 2 or 3D grids_mod
  !! @param current_state The current model state_mod
  !! @param ncid The NetCDF file id
  !! @param z_dim_id The NetCDF z dimension id that is provided by this procedure
  !! @param y_dim_id The NetCDF y dimension id that is provided by this procedure
  !! @param x_dim_id The NetCDF x dimension id that is provided by this procedure
  subroutine define_grid_dimensions(current_state, ncid, z_dim_id, y_dim_id, x_dim_id)
    type(model_state_type), intent(inout) :: current_state
    integer, intent(in) :: ncid
    integer, intent(out) :: z_dim_id, y_dim_id, x_dim_id

    integer :: empty_dim_id

    call check_status(nf90_def_dim(ncid, EMPTY_DIM_KEY, 1, empty_dim_id))

    if (current_state%global_grid%active(Z_INDEX)) then
      call check_status(nf90_def_dim(ncid, Z_DIM_KEY, current_state%global_grid%size(Z_INDEX), z_dim_id))
    else
      z_dim_id = empty_dim_id
    end if
    if (current_state%global_grid%active(Y_INDEX)) then
      call check_status(nf90_def_dim(ncid, Y_DIM_KEY, current_state%global_grid%size(Y_INDEX), y_dim_id))
    else
      y_dim_id = empty_dim_id
    end if
    if (current_state%global_grid%active(X_INDEX)) then
      call check_status(nf90_def_dim(ncid, X_DIM_KEY, current_state%global_grid%size(X_INDEX), x_dim_id))
    else
      x_dim_id = empty_dim_id
    end if
  end subroutine define_grid_dimensions

  !> Defines the NetCDF grid variables. This works for 1, 2 or 3D grids_mod
  !! @param current_state The current model state_mod
  !! @param ncid The NetCDF file id
  !! @param z_dim_id The NetCDF z dimension id
  !! @param y_dim_id The NetCDF y dimension id
  !! @param x_dim_id The NetCDF x dimension id
  !! @param x_id The NetCDF x variable id provided by this procedure
  !! @param y_id The NetCDF y variable id provided by this procedure
  !! @param z_id The NetCDF z variable id provided by this procedure
  subroutine define_grid_variables(current_state, ncid, z_dim_id, y_dim_id, x_dim_id, x_id, y_id, z_id)
    type(model_state_type), intent(inout) :: current_state
    integer, intent(in) :: ncid, z_dim_id, y_dim_id, x_dim_id
    integer, intent(out) :: x_id, y_id, z_id

    if (current_state%global_grid%active(X_INDEX)) then
      call check_status(nf90_def_var(ncid, X_KEY, NF90_REAL, x_dim_id, x_id))
      call check_status(nf90_put_att(ncid, x_id, "units", "m"))
    end if
    if (current_state%global_grid%active(Y_INDEX)) then
      call check_status(nf90_def_var(ncid, Y_KEY, NF90_REAL, y_dim_id, y_id))
      call check_status(nf90_put_att(ncid, y_id, "units", "m"))
    end if
    if (current_state%global_grid%active(Z_INDEX)) then
      call check_status(nf90_def_var(ncid, Z_KEY, NF90_REAL, z_dim_id, z_id))
      call check_status(nf90_put_att(ncid, z_id, "units", "m"))
    end if
  end subroutine define_grid_variables

  !> Defines the Q variable in the checkpoint file
  !! @param ncid The NetCDF file id
  !! @param multi_process Whether to support parallel IO operations or not
  !! @param q_dim_id The NetCDF q dimension id
  !! @param z_dim_id The NetCDF z dimension id
  !! @param y_dim_id The NetCDF y dimension id
  !! @param x_dim_id The NetCDF x dimension id
  !! @param q_id The NetCDF q variable id provided by this procedure
  !! @param zq_id The NetCDF zq variable id provided by this procedure
  subroutine define_q_variable(ncid, multi_process, q_dim_id, z_dim_id, y_dim_id, x_dim_id, q_id, zq_id)
    logical, intent(in) :: multi_process
    integer, intent(in) :: ncid, z_dim_id, y_dim_id, x_dim_id, q_dim_id
    integer, intent(out) :: q_id, zq_id

    integer, dimension(:), allocatable :: dimids

    allocate(dimids(4))
    dimids = (/ z_dim_id, y_dim_id, x_dim_id, q_dim_id /)

    call check_status(nf90_def_var(ncid, Q_KEY, merge(NF90_DOUBLE, NF90_REAL, &
         DEFAULT_PRECISION == DOUBLE_PRECISION), dimids, q_id))
    call check_status(nf90_def_var(ncid, ZQ_KEY, merge(NF90_DOUBLE, NF90_REAL, &
         DEFAULT_PRECISION == DOUBLE_PRECISION), dimids, zq_id))
    
    if (multi_process) then
      call check_status(nf90_def_var_fill(ncid, q_id, 1, 1))
      call check_status(nf90_def_var_fill(ncid, zq_id, 1, 1))
      call check_status(nf90_var_par_access(ncid, q_id, NF90_COLLECTIVE))
      call check_status(nf90_var_par_access(ncid, zq_id, NF90_COLLECTIVE))
    end if
  end subroutine define_q_variable

  !> Defines prognostic variables in the NetCDF. This handles 1, 2 and 3D grids_mod and 1, 2 and 3D fields, which
  !! most likely have the same dimensions but this is not mandatory here. All prognostic fields are 3D, if
  !! the grid is not 3D then the empty (size 1) dimension is used in that dimension
  !! @param current_state The current model state_mod
  !! @param ncid The NetCDF file id
  !! @param z_dim_id The NetCDF z dimension id
  !! @param y_dim_id The NetCDF y dimension id
  !! @param x_dim_id The NetCDF x dimension id
  !! @param u_id The u prognostic field id provided by this procedure
  !! @param v_id The v prognostic field id provided by this procedure
  !! @param w_id The w prognostic field id provided by this procedure
  !! @param th_id The theta prognostic field id provided by this procedure
  subroutine define_prognostic_variables(current_state, multi_process, ncid, z_dim_id, &
       y_dim_id, x_dim_id, u_id, v_id, w_id, th_id, p_id, zu_id, zv_id, zw_id, zth_id)
    type(model_state_type), intent(inout) :: current_state
    logical, intent(in) :: multi_process
    integer, intent(in) :: ncid, z_dim_id, y_dim_id, x_dim_id
    integer, intent(out) :: u_id, v_id, w_id, th_id, p_id, zu_id, zv_id, zw_id, zth_id

#ifdef U_ACTIVE
    call define_velocity_variable(ncid, multi_process, z_dim_id, y_dim_id, x_dim_id, field_name=U_KEY, field_id=u_id)
    call define_velocity_variable(ncid, multi_process, z_dim_id, y_dim_id, x_dim_id, field_name=ZU_KEY, field_id=zu_id)
#endif
#ifdef V_ACTIVE
    call define_velocity_variable(ncid, multi_process, z_dim_id, y_dim_id, x_dim_id, field_name=V_KEY, field_id=v_id)
    call define_velocity_variable(ncid, multi_process, z_dim_id, y_dim_id, x_dim_id, field_name=ZV_KEY, field_id=zv_id)
#endif
#ifdef W_ACTIVE
    call define_velocity_variable(ncid, multi_process, z_dim_id, y_dim_id, x_dim_id, field_name=W_KEY, field_id=w_id)
    call define_velocity_variable(ncid, multi_process, z_dim_id, y_dim_id, x_dim_id, field_name=ZW_KEY, field_id=zw_id)
#endif
    if (current_state%th%active) then
      call define_velocity_variable(ncid, multi_process, z_dim_id, y_dim_id, x_dim_id, field_name=TH_KEY, field_id=th_id)
      call define_velocity_variable(ncid, multi_process, z_dim_id, y_dim_id, x_dim_id, field_name=ZTH_KEY, field_id=zth_id)
    end if
    if (current_state%p%active) then 
      call define_velocity_variable(ncid, multi_process, z_dim_id, y_dim_id, x_dim_id, field_name=P_KEY, field_id=p_id)
    end if
  end subroutine define_prognostic_variables

  !> Defines misc variables in the NetCDF file
  !! @param ncid The NetCDF file id
  !! @param timestep_id The NetCDF timestep variable
  subroutine define_misc_variables(ncid, timestep_id, time_id, ugal_id, vgal_id, number_q_fields_id, &
       dtm_id, dtm_new_id, absolute_new_dtm_id)
    integer, intent(in) :: ncid
    integer, intent(out) :: timestep_id, time_id, ugal_id, vgal_id, number_q_fields_id, dtm_id, dtm_new_id, absolute_new_dtm_id

    call check_status(nf90_def_var(ncid, TIMESTEP, NF90_INT, timestep_id))
    call check_status(nf90_def_var(ncid, TIME_KEY, NF90_DOUBLE, time_id))
    call check_status(nf90_def_var(ncid, UGAL, NF90_DOUBLE, ugal_id))
    call check_status(nf90_def_var(ncid, VGAL, NF90_DOUBLE, vgal_id))
    call check_status(nf90_def_var(ncid, NQFIELDS, NF90_INT, number_q_fields_id))
    call check_status(nf90_def_var(ncid, DTM_KEY, NF90_DOUBLE, dtm_id))
    call check_status(nf90_def_var(ncid, DTM_NEW_KEY, NF90_DOUBLE, dtm_new_id))
    call check_status(nf90_def_var(ncid, ABSOLUTE_NEW_DTM_KEY, NF90_DOUBLE, absolute_new_dtm_id))
  end subroutine define_misc_variables

  !> Will dump out (write) misc model data to the checkpoint
  !! @param current_state The current model state_mod
  !! @param ncid The NetCDF file id
  !! @param timestep_id The NetCDF timestep variable id
  subroutine write_out_misc_variables(current_state, ncid, timestep_id, time_id, ugal_id, &
       vgal_id, number_q_fields_id, dtm_id, dtm_new_id, absolute_new_dtm_id)
    type(model_state_type), intent(inout) :: current_state
    integer, intent(in) :: ncid, timestep_id, time_id, ugal_id, vgal_id, number_q_fields_id, &
         dtm_id, dtm_new_id, absolute_new_dtm_id

    call check_status(nf90_put_var(ncid, timestep_id, current_state%timestep))
    ! The time is incremented with dtm as the model was about to increment for the next step and this is needed for diagnostics
    call check_status(nf90_put_var(ncid, time_id, current_state%time+current_state%dtm))
    call check_status(nf90_put_var(ncid, ugal_id, current_state%ugal))
    call check_status(nf90_put_var(ncid, vgal_id, current_state%vgal))
    call check_status(nf90_put_var(ncid, number_q_fields_id, current_state%number_q_fields))
    call check_status(nf90_put_var(ncid, dtm_id, current_state%dtm))
    call check_status(nf90_put_var(ncid, dtm_new_id, current_state%dtm_new))
    call check_status(nf90_put_var(ncid, absolute_new_dtm_id, current_state%absolute_new_dtm))
  end subroutine write_out_misc_variables

  !> Will define a single velocity variable in the NetCDF file
  !! @param ncid The NetCDF file id
  !! @param dimone Size in first dimension
  !! @param dimtwo Size in second dimension (optional)
  !! @param dimthree Size in third dimension (optional)
  !! @param field_name The NetCDF name of the variable
  !! @param field_id The NetCDF variable id produced by this procedure
  subroutine define_velocity_variable(ncid, multi_process, dimone, dimtwo, dimthree, field_name, field_id)
    integer, intent(in) :: ncid, dimone
    integer, intent(in), optional :: dimtwo, dimthree
    integer, intent(out) :: field_id
    character(len=*), intent(in) :: field_name
    logical, intent(in) :: multi_process

    integer, dimension(:), allocatable :: dimids

    if (present(dimtwo) .and. present(dimthree)) then
      allocate(dimids(3))
      dimids = (/ dimone, dimtwo, dimthree /)
    else if (present(dimtwo) .or. present(dimthree)) then
      allocate(dimids(2))
      dimids = (/ dimone, merge(dimtwo, dimthree, present(dimtwo)) /)
    else
      allocate(dimids(1))
      dimids = (/ dimone /)
    end if

    call check_status(nf90_def_var(ncid, field_name, merge(NF90_DOUBLE, NF90_REAL, DEFAULT_PRECISION == DOUBLE_PRECISION), &
         dimids, field_id))
    if (multi_process) then
      call check_status(nf90_def_var_fill(ncid, field_id, 1, 1))
      call check_status(nf90_var_par_access(ncid, field_id, NF90_COLLECTIVE))
    end if
    call check_status(nf90_put_att(ncid, field_id, "units", "m/s"))
  end subroutine define_velocity_variable
end module checkpointer_write_checkpoint_mod
