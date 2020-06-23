!> The NetCDF file type writer which performs actual writing of NetCDF files to the parallel filesystem. These are opened by
!! all IO servers and all IO servers can participate as variables might be located across the different IO processes
module netcdf_filetype_writer_mod
  use datadefn_mod, only : DEFAULT_PRECISION, STRING_LENGTH, SINGLE_PRECISION, DOUBLE_PRECISION
  use configuration_parser_mod, only : INSTANTANEOUS_TYPE, TIME_AVERAGED_TYPE, io_configuration_type, &
       data_values_type, get_data_value_by_field_name, &
       cond_request, diag_request, cond_long, diag_long
  use collections_mod, only : hashmap_type, hashset_type, list_type, map_type, iterator_type, mapentry_type, &
       c_get_generic, c_get_integer, c_get_string, c_contains, c_generic_at, c_real_at, c_integer_at, c_put_generic, &
       c_put_integer, c_remove, c_free, c_has_next, c_get_iterator, c_next_mapentry, c_next_generic, c_get_real, c_size, &
       c_next_string, c_is_empty, c_add_string
  use conversions_mod, only : conv_to_integer, conv_to_string, conv_to_real
  use logging_mod, only : LOG_ERROR, LOG_WARN, LOG_DEBUG, log_log, log_master_log, log_get_logging_level, log_is_master
  use writer_types_mod, only : writer_type, writer_field_type, write_field_collective_values_type, &
       netcdf_diagnostics_timeseries_type, netcdf_diagnostics_type, write_field_collective_descriptor_type, &
       write_field_collective_monc_info_type
  use forthread_mod, only : forthread_mutex_init, forthread_mutex_lock, forthread_mutex_unlock, forthread_mutex_destroy, &
       forthread_rwlock_rdlock, forthread_rwlock_wrlock, forthread_rwlock_unlock, forthread_rwlock_init, forthread_rwlock_destroy
  use threadpool_mod, only : check_thread_status
  use netcdf, only : NF90_DOUBLE, NF90_REAL, NF90_INT, NF90_CHAR, NF90_GLOBAL, NF90_CLOBBER, NF90_NETCDF4, NF90_MPIIO, &
       NF90_COLLECTIVE, NF90_UNLIMITED, nf90_def_var, nf90_var_par_access, nf90_def_var_fill, nf90_put_att, &
       nf90_create, nf90_put_var, nf90_def_dim, nf90_enddef, nf90_close, nf90_ebaddim, nf90_enotatt, nf90_enotvar, &
       nf90_noerr, nf90_strerror, nf90_redef, nf90_inq_varid
  use io_server_client_mod, only : ARRAY_FIELD_TYPE, SCALAR_FIELD_TYPE, MAP_FIELD_TYPE, DOUBLE_DATA_TYPE, &
       INTEGER_DATA_TYPE, STRING_DATA_TYPE
  use io_server_state_writer_mod, only : define_io_server_state_contributions, write_io_server_state
  use mpi, only : MPI_INFO_NULL
  use mpi_communication_mod, only : lock_mpi, unlock_mpi, wait_for_mpi_request
  use grids_mod, only : Z_INDEX, Y_INDEX, X_INDEX
  use optionsdatabase_mod, only : options_size, options_value_at, options_key_at
  use netcdf_misc_mod, only : check_netcdf_status
  use mpi, only : MPI_STATUS_IGNORE, MPI_REQUEST_NULL, MPI_INT

  implicit none

#ifndef TEST_MODE
  private
#endif

  type(hashmap_type), volatile :: file_states
  integer, volatile :: file_states_rwlock, netcdf_mutex

  logical :: l_nc_dim, l_nd_dim
  integer :: nc_dim_id, nd_dim_id, nopt_dim_id
  integer :: nc_var_id_s, nd_var_id_s, nc_var_id_l, nd_var_id_l, nopt_var_id

  public initialise_netcdf_filetype, finalise_netcdf_filetype, define_netcdf_file, write_variable, close_netcdf_file, &
       store_io_server_state, get_writer_entry_from_netcdf
contains

  !> Initialises the NetCDF writing functionality
  subroutine initialise_netcdf_filetype()
    call check_thread_status(forthread_rwlock_init(file_states_rwlock, -1))
    call check_thread_status(forthread_mutex_init(netcdf_mutex, -1))
  end subroutine initialise_netcdf_filetype

  !> Finalises  the NetCDF writing functionality
  subroutine finalise_netcdf_filetype()
    call check_thread_status(forthread_rwlock_destroy(file_states_rwlock))
    call check_thread_status(forthread_mutex_destroy(netcdf_mutex))
  end subroutine finalise_netcdf_filetype  
  
  !> Defines a NetCDF file - which creates it, defines all dimensions and variables. This must be called by all IO server
  !! processes as the NetCDF operations here are collective
  !! @param io_configuration The IO server configuration
  !! @param file_writer_information The writer entry that is being written
  !! @param timestep The write timestep
  !! @param time The write time
  subroutine define_netcdf_file(io_configuration, file_writer_information, timestep, time, time_points, termination_write)
    type(io_configuration_type), intent(inout) :: io_configuration
    type(writer_type), intent(inout), target :: file_writer_information
    type(map_type), intent(inout) :: time_points
    integer, intent(in) :: timestep
    real, intent(in) :: time
    logical, intent(in) :: termination_write
    
    character(len=STRING_LENGTH) :: unique_filename
    type(netcdf_diagnostics_type), pointer :: ncdf_writer_state
    class(*), pointer :: generic
    integer :: zn_var_id
    integer :: z_var_id

    ncdf_writer_state=>get_file_state(file_writer_information%filename, timestep, .true.)
    if (.not. associated(ncdf_writer_state)) then
      call check_thread_status(forthread_rwlock_wrlock(file_states_rwlock))
      ncdf_writer_state=>get_file_state(file_writer_information%filename, timestep, .false.)
      if (.not. associated(ncdf_writer_state)) then
        allocate(ncdf_writer_state)
        ncdf_writer_state%corresponding_writer_entry=>file_writer_information
        ncdf_writer_state%termination_write=termination_write
        call check_thread_status(forthread_mutex_init(ncdf_writer_state%mutex, -1))
        call check_thread_status(forthread_mutex_lock(ncdf_writer_state%mutex))
        generic=>ncdf_writer_state
        call c_put_generic(file_states, trim(file_writer_information%filename)//"#"//trim(conv_to_string(timestep)), &
             generic, .false.)
        call check_thread_status(forthread_rwlock_unlock(file_states_rwlock))
        
        if (file_writer_information%write_on_model_time) then
          call generate_unique_filename(file_writer_information%filename, unique_filename, &
               file_writer_information%defined_write_time)
        else
          call generate_unique_filename(file_writer_information%filename, unique_filename, timestep=timestep)
        end if
        call check_thread_status(forthread_mutex_lock(netcdf_mutex))
        call lock_mpi()
        call check_netcdf_status(nf90_create(unique_filename, ior(NF90_NETCDF4, NF90_MPIIO), ncdf_writer_state%ncid, &
             comm = io_configuration%io_communicator, info = MPI_INFO_NULL))
        call unlock_mpi()
        call write_out_global_attributes(io_configuration, ncdf_writer_state%ncid, file_writer_information, timestep, time)
        call define_dimensions(ncdf_writer_state, io_configuration%dimension_sizing)
        call define_time_series_dimensions(ncdf_writer_state, file_writer_information, time, time_points, termination_write)
        call define_variables(io_configuration, ncdf_writer_state, file_writer_information)
        zn_var_id = define_coordinate_variable(ncdf_writer_state,"zn")  
        z_var_id  = define_coordinate_variable(ncdf_writer_state,"z")
        nopt_var_id = define_options_database_variable(ncdf_writer_state)
        call lock_mpi()
        call check_netcdf_status(nf90_enddef(ncdf_writer_state%ncid))
        call unlock_mpi()

        !> Isolate these writes, as they may trip over one another
        if (io_configuration%my_io_rank == 0) then
          !> Write conditional diagnostics descriptors to file if file contains corresponding dimension.
          if (l_nc_dim) then
            call write_condition_variable(ncdf_writer_state, nc_var_id_s, cond_request)
            call write_condition_variable(ncdf_writer_state, nc_var_id_l, cond_long)
          end if
          if (l_nd_dim) then
            call write_condition_variable(ncdf_writer_state, nd_var_id_s, diag_request)
            call write_condition_variable(ncdf_writer_state, nd_var_id_l, diag_long)
          end if

          !> Write options_database and height coordinates to all files.
          call write_out_options(io_configuration, ncdf_writer_state)
          call write_coordinate_variable(ncdf_writer_state, zn_var_id, io_configuration%zn_field)
          call write_coordinate_variable(ncdf_writer_state, z_var_id,  io_configuration%z_field)
        end if ! end write isolation

        call check_thread_status(forthread_mutex_unlock(ncdf_writer_state%mutex))        
        call check_thread_status(forthread_mutex_unlock(netcdf_mutex))
 
      else
        call check_thread_status(forthread_rwlock_unlock(file_states_rwlock))
      end if
    end if
  end subroutine define_netcdf_file

  !> Stores the IO server state in the NetCDF file
  !! @param io_configuration The IO server configuration
  !! @param file_writer_information The file writer information
  !! @param timestep The write timestep
  subroutine store_io_server_state(io_configuration, writer_entries, time_points, file_writer_information, timestep)
    type(io_configuration_type), intent(inout) :: io_configuration
    type(writer_type), volatile, dimension(:), intent(inout) :: writer_entries
    type(hashmap_type), volatile, intent(inout) :: time_points
    type(writer_type), intent(inout), target :: file_writer_information
    integer, intent(in) :: timestep

    type(netcdf_diagnostics_type), pointer :: ncdf_writer_state

    ncdf_writer_state=>get_file_state(file_writer_information%filename, timestep, .true.)
    call lock_mpi()
    call check_netcdf_status(nf90_redef(ncdf_writer_state%ncid))
    call unlock_mpi()
    call define_io_server_state_contributions(io_configuration, writer_entries, time_points, ncdf_writer_state)
    call lock_mpi()
    call check_netcdf_status(nf90_enddef(ncdf_writer_state%ncid))
    call unlock_mpi()
    call write_io_server_state(io_configuration, writer_entries, time_points, ncdf_writer_state)
  end subroutine store_io_server_state

  !> Looks up and retrieves the writer entry that corresponds to this NetCDF file state
  !! @param field_name The field name that is being communicated
  !! @param timestep The write timestep
  !! @returns The writer entry
  function get_writer_entry_from_netcdf(field_name, timestep, terminated)
    character(len=STRING_LENGTH) :: field_name
    integer :: timestep
    logical, intent(out), optional :: terminated
    type(writer_type), pointer :: get_writer_entry_from_netcdf

    type(netcdf_diagnostics_type), pointer :: file_state

    file_state=>get_file_state(field_name, timestep, .true.)
    if (present(terminated)) terminated=file_state%termination_write
    get_writer_entry_from_netcdf=>file_state%corresponding_writer_entry
  end function get_writer_entry_from_netcdf  

  !> Call back for the inter IO reduction which actually does the NetCDF file closing which is a 
  !! collective (synchronous) operation. This also cleans up the file state as it is no longer required
  !! @param io_configuration The IO server configuration
  !! @param field_name The field name that is being communicated
  !! @param timestep The write timestep
  subroutine close_netcdf_file(io_configuration, field_name, timestep)
    type(io_configuration_type), intent(inout) :: io_configuration
    character(len=STRING_LENGTH) :: field_name
    integer :: timestep

    type(iterator_type) :: iterator
    type(netcdf_diagnostics_timeseries_type), pointer :: ptr
    class(*), pointer :: generic

    type(netcdf_diagnostics_type), pointer :: file_state

    file_state=>get_file_state(field_name, timestep, .true.)
    call check_thread_status(forthread_mutex_lock(file_state%mutex))
    call check_thread_status(forthread_mutex_lock(netcdf_mutex))
    call lock_mpi()
    call check_netcdf_status(nf90_close(file_state%ncid))
    call unlock_mpi()
    call check_thread_status(forthread_mutex_unlock(netcdf_mutex))
    call check_thread_status(forthread_mutex_unlock(file_state%mutex))
    call check_thread_status(forthread_mutex_destroy(file_state%mutex))
    call c_free(file_state%dimension_to_id)
    call c_free(file_state%variable_to_id)
    iterator=c_get_iterator(file_state%timeseries_dimension)
    do while (c_has_next(iterator))
      generic=>c_get_generic(c_next_mapentry(iterator))
      select type(generic)
      type is(netcdf_diagnostics_timeseries_type)
        ptr => generic
        deallocate(ptr)
      end select      
    end do  
    call c_free(file_state%timeseries_dimension)
    call check_thread_status(forthread_rwlock_wrlock(file_states_rwlock))
    call c_remove(file_states, trim(field_name)//"#"//trim(conv_to_string(timestep)))
    call check_thread_status(forthread_rwlock_unlock(file_states_rwlock))    
    if (log_get_logging_level() .ge. LOG_DEBUG .and. log_is_master()) then
      call log_master_log(LOG_DEBUG, "Done physical close for NetCDF file at timestep "//trim(conv_to_string(timestep)))
    end if    
  end subroutine close_netcdf_file 

  !> Writes the contents of a variable to the NetCDF file. This also removes the written entries from the field information
  !! type in order to conserve memory
  !! @param field_to_write_information The field that is going to be written
  !! @param filename The name of the diagnostics file
  !! @param timestep The write timestep
  !! @param time The write time
  subroutine write_variable(io_configuration, field_to_write_information, filename, timestep, time)
    type(io_configuration_type), intent(inout) :: io_configuration
    type(writer_field_type), intent(inout) :: field_to_write_information
    character(len=*), intent(in) :: filename
    integer, intent(in) :: timestep
    real, intent(in) :: time

    type(netcdf_diagnostics_type), pointer :: file_state
            
    file_state=>get_file_state(filename, timestep, .true.)
    if (field_to_write_information%collective_write) then
      if (field_to_write_information%collective_contiguous_optimisation) then
        call write_contiguous_collective_variable_to_diagnostics(io_configuration, field_to_write_information, &
             timestep, time, file_state)
      else
        call write_collective_variable_to_diagnostics(io_configuration, field_to_write_information, timestep, time, file_state)
      end if
    else
      call write_independent_variable_to_diagnostics(field_to_write_information, timestep, time, file_state)
    end if
  end subroutine write_variable 

  !> Translates a dimension name to its numeric corresponding identifier
  !! @param dim_name The name of the dimension to look up
  !! @param is_auto_dimension Optional parameter determining whether dimension is auto or not
  !! @returns Corresponding identifier
  integer function get_dimension_identifier(dim_name, is_auto_dimension)
    character(len=*), intent(in) :: dim_name
    logical, intent(out), optional :: is_auto_dimension

    integer :: dash_idx
    logical :: is_modified_size

    dash_idx=index(dim_name, "_")
    dash_idx=dash_idx-1
    is_modified_size=dash_idx .ne. -1
    if (.not. is_modified_size) dash_idx=len_trim(dim_name)

    if (dim_name(:dash_idx) .eq. "z" .or. dim_name(:dash_idx) .eq. "zn") then
      get_dimension_identifier=Z_INDEX
    else if (dim_name(:dash_idx) .eq. "y") then
      get_dimension_identifier=Y_INDEX
    else if (dim_name(:dash_idx) .eq. "x") then
      get_dimension_identifier=X_INDEX
    else
      get_dimension_identifier=-1
    end if

    if (present(is_auto_dimension)) is_auto_dimension=is_modified_size
  end function get_dimension_identifier

  !> Retrieves the original size of a specific dimension (which is auto)
  !! @param dim_name The auto dimension (full) name
  !! @param dimension_store The map of dimensions we are looking up
  !! @returns The original size of this auto dimension
  integer function get_dimension_original_size(dim_name, dimension_store)
    character(len=*), intent(in) :: dim_name
    type(map_type), intent(inout) :: dimension_store

    integer :: dash_idx

    dash_idx=index(dim_name, "_")
    dash_idx=dash_idx-1
    if (dash_idx .eq. -1) dash_idx=len_trim(dim_name)

    get_dimension_original_size=c_get_integer(dimension_store, dim_name(:dash_idx))
  end function get_dimension_original_size  

  !> Writes a coordinate variable into the NetCDF file
  !! @param file_state The NetCDF file state
  !! @param coord_var_id The coordinate variable id in the file
  !! @param field_values The field values to write
  subroutine write_coordinate_variable(file_state, coord_var_id, field_values)
    type(netcdf_diagnostics_type), intent(inout) :: file_state
    integer, intent(in) :: coord_var_id
    real(kind=DEFAULT_PRECISION), dimension(:), intent(in) :: field_values

    integer :: count_to_write(1)

    count_to_write(1)=size(field_values)
    call lock_mpi()
    call check_netcdf_status(nf90_put_var(file_state%ncid, coord_var_id, field_values, count=count_to_write))
    call unlock_mpi()
  end subroutine write_coordinate_variable

  !> Writes the conditional diagnostic variable names into the NetCDF file
  !! @param file_state The NetCDF file state
  !! @param c_var_id The variable id in the file
  !! @param field_values The field values to write
  subroutine write_condition_variable(file_state, c_var_id, field_values)
    type(netcdf_diagnostics_type), intent(inout) :: file_state
    integer, intent(in) :: c_var_id 
    character(len=STRING_LENGTH), dimension(:), intent(in) :: field_values

    integer :: count_to_write(2), start_pos(2)
    integer :: pos, string_size
    character(len=STRING_LENGTH) :: dum_string

    count_to_write(2)=1    ! element count
    start_pos(1)=1         ! string character start position
    call lock_mpi()
    do pos=1,size(field_values)
      dum_string = trim(field_values(pos))
      count_to_write(1) = len(trim(field_values(pos)))
      start_pos(2)=pos
      call check_netcdf_status(nf90_put_var(file_state%ncid, c_var_id, dum_string, &
                               start=start_pos, count=count_to_write))
    end do
    call unlock_mpi()
  end subroutine write_condition_variable

  !> Writes contiguous collective variable blocks into the NetCDF files. These are blocks of data spanning multiple variables
  !! and multiple time points, the idea being to minimise the number of overall writes into the file. These variables are
  !! defined collective, so it might be that other IO servers have less contiguous data and hence more writes, therefore
  !! empty writes might need to be issued (at the end) to match up against these other servers
  !! @param io_configuration The IO server configuration
  !! @param field_to_write_information The field to write
  !! @param timestep The current model timestep
  !! @param time The current model time
  !! @param file_state Specific file state descriptor
  subroutine write_contiguous_collective_variable_to_diagnostics(io_configuration, field_to_write_information, timestep, &
       time, file_state)
    type(io_configuration_type), intent(inout) :: io_configuration
    type(writer_field_type), intent(inout) :: field_to_write_information
    integer, intent(in) :: timestep
    real, intent(in) :: time
    type(netcdf_diagnostics_type), intent(inout) :: file_state

    real :: value_to_test
    real(kind=DEFAULT_PRECISION), dimension(:,:,:,:), allocatable :: contiguous_values
    type(hashset_type) :: items_to_remove
    real(kind=DEFAULT_PRECISION), dimension(:), allocatable :: timeseries_time_to_write
    type(iterator_type) :: value_to_write_iterator, collective_descriptor_iterator, monc_iterator, value_to_remove_iterator
    type(mapentry_type) :: value_to_write_map_entry
    class(*), pointer :: generic
    type(write_field_collective_descriptor_type), pointer :: collective_descriptor
    type(netcdf_diagnostics_timeseries_type), pointer :: timeseries_diag
    type(write_field_collective_monc_info_type), pointer :: monc_descriptor
    type(write_field_collective_values_type), pointer :: multi_monc_entries
    type(data_values_type), pointer :: data_value
    integer :: number_time_entries, i, start(4), count(4), field_id, ierr
    character(len=STRING_LENGTH) :: removal_key

    timeseries_diag=>get_specific_timeseries_dimension(file_state, field_to_write_information%output_frequency, &
         field_to_write_information%timestep_frequency)
    field_id=c_get_integer(file_state%variable_to_id, get_field_key(field_to_write_information))
    if (.not. timeseries_diag%variable_written) allocate(timeseries_time_to_write(timeseries_diag%num_entries))

    number_time_entries=0
    value_to_write_iterator=c_get_iterator(field_to_write_information%values_to_write)
    do while (c_has_next(value_to_write_iterator))
      value_to_write_map_entry=c_next_mapentry(value_to_write_iterator)
      value_to_test=conv_to_real(value_to_write_map_entry%key)
      if (value_to_test .le. time .and. value_to_test .gt. field_to_write_information%previous_write_time) then
        number_time_entries=number_time_entries+1
        call c_add_string(items_to_remove, value_to_write_map_entry%key)
      end if
    end do
    
    if (number_time_entries .ne. timeseries_diag%num_entries) then
      call log_log(LOG_WARN, "Expected "//trim(conv_to_string(timeseries_diag%num_entries))//&
           " but have "//trim(conv_to_string(number_time_entries)))
      if (number_time_entries .gt. timeseries_diag%num_entries) number_time_entries=timeseries_diag%num_entries
    end if

    collective_descriptor_iterator=c_get_iterator(field_to_write_information%collective_descriptors)
    do while (c_has_next(collective_descriptor_iterator))
      collective_descriptor=>get_next_collective_descriptor(collective_descriptor_iterator) 
      allocate(contiguous_values(collective_descriptor%count(1), collective_descriptor%count(2), &
           collective_descriptor%count(3), number_time_entries))
      monc_iterator=c_get_iterator(collective_descriptor%specific_monc_info)
      do while (c_has_next(monc_iterator))
        monc_descriptor=>get_next_specific_monc_info(monc_iterator)
        value_to_write_iterator=c_get_iterator(field_to_write_information%values_to_write)
        i=0
        do while (c_has_next(value_to_write_iterator))
          value_to_write_map_entry=c_next_mapentry(value_to_write_iterator)
          value_to_test=conv_to_real(value_to_write_map_entry%key)
          if (value_to_test .le. time .and. value_to_test .gt. field_to_write_information%previous_write_time) then
            i=i+1
            if (allocated(timeseries_time_to_write)) timeseries_time_to_write(i)=value_to_test
            generic=>c_get_generic(value_to_write_map_entry)
            select type(generic)
            type is(write_field_collective_values_type)
              multi_monc_entries=>generic
            end select
            data_value=>get_data_value_by_field_name(multi_monc_entries%monc_values, &
                 trim(conv_to_string(monc_descriptor%monc_source)))
            if (collective_descriptor%split_dim == Y_INDEX) then
              contiguous_values(:,monc_descriptor%relative_dimension_start:monc_descriptor%relative_dimension_start+&
                   monc_descriptor%counts(Y_INDEX)-1,:,i)=reshape(data_value%values, (/ monc_descriptor%counts(Z_INDEX), &
                   monc_descriptor%counts(Y_INDEX), monc_descriptor%counts(X_INDEX)/))
            else
              contiguous_values(:,:,monc_descriptor%relative_dimension_start:monc_descriptor%relative_dimension_start+&
                   monc_descriptor%counts(X_INDEX)-1,i)=reshape(data_value%values, (/ monc_descriptor%counts(Z_INDEX), &
                   monc_descriptor%counts(Y_INDEX), monc_descriptor%counts(X_INDEX)/))
            end if
            deallocate(data_value%values)
            deallocate(data_value)            
          end if
        end do
      end do
      count(1:3)=collective_descriptor%count
      count(4)=number_time_entries
      start(1:3)=collective_descriptor%absolute_start
      start(4)=1
      call check_thread_status(forthread_mutex_lock(file_state%mutex))
      call check_thread_status(forthread_mutex_lock(netcdf_mutex))
      call lock_mpi()
      call check_netcdf_status(nf90_put_var(file_state%ncid, field_id, contiguous_values, start=start, count=count))
      call unlock_mpi()
      call check_thread_status(forthread_mutex_unlock(netcdf_mutex))
      call check_thread_status(forthread_mutex_unlock(file_state%mutex))
      deallocate(contiguous_values)
      if (allocated(timeseries_time_to_write)) then
        call check_thread_status(forthread_mutex_lock(netcdf_mutex))
        call lock_mpi()
        call check_netcdf_status(nf90_put_var(file_state%ncid, timeseries_diag%netcdf_var_id, &
             timeseries_time_to_write, count=(/ timeseries_diag%num_entries /)))
        call unlock_mpi()
        call check_thread_status(forthread_mutex_unlock(netcdf_mutex))
        timeseries_diag%variable_written=.true.
        deallocate(timeseries_time_to_write)
      end if
    end do
    if (.not. c_is_empty(items_to_remove)) then
      value_to_remove_iterator=c_get_iterator(items_to_remove)
      do while (c_has_next(value_to_remove_iterator))
        removal_key=c_next_string(value_to_remove_iterator)
        generic=>c_get_generic(field_to_write_information%values_to_write, removal_key)
        select type(generic)
        type is(write_field_collective_values_type)
          multi_monc_entries=>generic
        end select
        call c_free(multi_monc_entries%monc_values)
        deallocate(multi_monc_entries)
        call c_remove(field_to_write_information%values_to_write, removal_key)
      end do
      call c_free(items_to_remove)
    end if
    if (field_to_write_information%max_num_collective_writes_request_handle .ne. MPI_REQUEST_NULL) then
      call wait_for_mpi_request(field_to_write_information%max_num_collective_writes_request_handle)
    end if
    if (c_size(field_to_write_information%collective_descriptors) .lt. field_to_write_information%max_num_collective_writes) then      
      call check_thread_status(forthread_mutex_lock(file_state%mutex))
      call check_thread_status(forthread_mutex_lock(netcdf_mutex))
      call lock_mpi()
      do i=c_size(field_to_write_information%collective_descriptors), field_to_write_information%max_num_collective_writes-1        
        call check_netcdf_status(nf90_put_var(file_state%ncid, field_id, (/1.0/), start=(/1/), count=(/0/)))
      end do
      call unlock_mpi()
      call check_thread_status(forthread_mutex_unlock(netcdf_mutex))
      call check_thread_status(forthread_mutex_unlock(file_state%mutex))
    end if
  end subroutine write_contiguous_collective_variable_to_diagnostics
  
  !> Retrieves the next collective descriptor based upon the iterator
  !! @param iterator The iterator to retrieve from, this is updated to reference the proceeding entry
  !! @returns The next collective descriptor
  function get_next_collective_descriptor(iterator)
    type(iterator_type), intent(inout) :: iterator
    type(write_field_collective_descriptor_type), pointer :: get_next_collective_descriptor

    class(*), pointer :: generic

    generic=>c_next_generic(iterator)
    select type(generic)
    type is (write_field_collective_descriptor_type)
      get_next_collective_descriptor=>generic
    end select
  end function get_next_collective_descriptor

  !> Retrieves the next specific monc information item from the iterator
  !! @param iterator The iterator to retrieve from, this is updated to reference the proceeding entry
  !! @returns The next specific monc information item
  function get_next_specific_monc_info(iterator)
    type(iterator_type) :: iterator
    type(write_field_collective_monc_info_type), pointer :: get_next_specific_monc_info

    class(*), pointer :: generic

    generic=>c_next_generic(iterator)
    select type(generic)
    type is (write_field_collective_monc_info_type)
      get_next_specific_monc_info=>generic
    end select
  end function get_next_specific_monc_info

  !> Writes collective variables, where we are working with the values from multiple MONCs and storing these in their own
  !! specific relative location in the diagnostics file. 
  !! @param field_to_write_information The field that is going to be written
  !! @param timestep The write timestep
  !! @param time The write time
  !! @param file_state File storate state
  subroutine write_collective_variable_to_diagnostics(io_configuration, field_to_write_information, timestep, time, file_state)
    type(io_configuration_type), intent(inout) :: io_configuration
    type(writer_field_type), intent(inout) :: field_to_write_information
    integer, intent(in) :: timestep
    real, intent(in) :: time
    type(netcdf_diagnostics_type), intent(inout) :: file_state

    real :: value_to_test
    integer :: i, k, included_num, field_id, start(field_to_write_information%dimensions+1), &
         count(field_to_write_information%dimensions+1), monc_location, dim_identifier, auto_period, dim_start
    class(*), pointer :: generic
    type(write_field_collective_values_type), pointer :: multi_monc_entries
    logical :: is_auto_dimension
    real(kind=DEFAULT_PRECISION), dimension(:), allocatable :: timeseries_time_to_write
    character(len=STRING_LENGTH), dimension(:), allocatable :: items_to_remove
    type(data_values_type), pointer :: data_value
    type(netcdf_diagnostics_timeseries_type), pointer :: timeseries_diag
    type(iterator_type) :: value_to_write_iterator, monc_entries_iterator
    type(mapentry_type) :: value_to_write_map_entry, monc_entries_map_entry

    timeseries_diag=>get_specific_timeseries_dimension(file_state, field_to_write_information%output_frequency, &
         field_to_write_information%timestep_frequency)
    if (.not. timeseries_diag%variable_written) allocate(timeseries_time_to_write(timeseries_diag%num_entries))

    allocate(items_to_remove(timeseries_diag%num_entries))
    included_num=1
    field_id=c_get_integer(file_state%variable_to_id, get_field_key(field_to_write_information))
    value_to_write_iterator=c_get_iterator(field_to_write_information%values_to_write)
    do while (c_has_next(value_to_write_iterator))
      value_to_write_map_entry=c_next_mapentry(value_to_write_iterator)
      value_to_test=conv_to_real(value_to_write_map_entry%key)
      if (value_to_test .le. time .and. value_to_test .gt. field_to_write_information%previous_write_time) then
        if (included_num .le. timeseries_diag%num_entries) then
          if (allocated(timeseries_time_to_write)) timeseries_time_to_write(included_num)=value_to_test
          generic=>c_get_generic(value_to_write_map_entry)
          select type(generic)
          type is(write_field_collective_values_type)
            multi_monc_entries=>generic
          end select
          monc_entries_iterator=c_get_iterator(multi_monc_entries%monc_values)
          do while(c_has_next(monc_entries_iterator))
            monc_entries_map_entry=c_next_mapentry(monc_entries_iterator)
            data_value=>get_data_value_by_field_name(multi_monc_entries%monc_values, monc_entries_map_entry%key)           
            monc_location=c_get_integer(io_configuration%monc_to_index, monc_entries_map_entry%key)
            do k=1, field_to_write_information%dimensions
              dim_identifier=get_dimension_identifier(field_to_write_information%dim_size_defns(k), is_auto_dimension)
              if (dim_identifier .gt. -1) then
                start(k)=io_configuration%registered_moncs(monc_location)%local_dim_starts(dim_identifier)
                count(k)=io_configuration%registered_moncs(monc_location)%local_dim_sizes(dim_identifier)
                if (is_auto_dimension) then              
                  auto_period=ceiling(real(get_dimension_original_size(field_to_write_information%dim_size_defns(k), &
                       io_configuration%dimension_sizing))/field_to_write_information%actual_dim_size(k))
                  start(k)=(start(k)/auto_period)+1
                  if (io_configuration%registered_moncs(monc_location)%local_dim_starts(dim_identifier)==1) then
                    dim_start=1
                  else
                    dim_start=auto_period - &
                         mod(io_configuration%registered_moncs(monc_location)%local_dim_starts(dim_identifier)-2, auto_period)
                  end if
                  count(k)=ceiling(real(io_configuration%registered_moncs(monc_location)%local_dim_sizes(dim_identifier) - &
                       (dim_start-1))/auto_period)
                end if
              else
                call log_log(LOG_ERROR, "Can not locate dimension "//trim(field_to_write_information%dim_size_defns(k)))
              end if
            end do
            start(field_to_write_information%dimensions+1) = included_num
            count(field_to_write_information%dimensions+1) = 1
            call check_thread_status(forthread_mutex_lock(file_state%mutex))
            call check_thread_status(forthread_mutex_lock(netcdf_mutex))
            call lock_mpi()
            call check_netcdf_status(nf90_put_var(file_state%ncid, field_id, data_value%values, start=start, count=count))
            call unlock_mpi()
            call check_thread_status(forthread_mutex_unlock(netcdf_mutex))
            call check_thread_status(forthread_mutex_unlock(file_state%mutex))
            deallocate(data_value%values)
            deallocate(data_value)
          end do
          items_to_remove(included_num)=value_to_write_map_entry%key
          included_num=included_num+1
          call c_free(multi_monc_entries%monc_values)
          deallocate(multi_monc_entries)
        else
          call log_log(LOG_WARN, "Omitted time entry of field '"//trim(field_to_write_information%field_name)//&
               "' as past dimension length at time "//conv_to_string(value_to_test))
        end if
      end if
    end do
    if (included_num-1 .ne. timeseries_diag%num_entries) then
      call log_log(LOG_WARN, "Miss match of time entries for field '"//trim(field_to_write_information%field_name)//&
           "', included entries="//trim(conv_to_string(included_num-1))//" but expected entries="//&
           trim(conv_to_string(timeseries_diag%num_entries)))
    end if
    if (allocated(timeseries_time_to_write)) then
      call check_thread_status(forthread_mutex_lock(netcdf_mutex))
      call lock_mpi()
      call check_netcdf_status(nf90_put_var(file_state%ncid, timeseries_diag%netcdf_var_id, &
           timeseries_time_to_write, count=(/ timeseries_diag%num_entries /)))
      call unlock_mpi()
      call check_thread_status(forthread_mutex_unlock(netcdf_mutex))
      timeseries_diag%variable_written=.true.
    end if
    if (included_num .gt. 1) then
      do i=1, included_num-1
        call c_remove(field_to_write_information%values_to_write, items_to_remove(i))
      end do
    end if
    deallocate(items_to_remove)
    if (allocated(timeseries_time_to_write)) deallocate(timeseries_time_to_write)
  end subroutine write_collective_variable_to_diagnostics  

  !> Writes independent variables to the diagnostics file. This writes the entire variable and works by writing all the 
  !! memory up and then doing the write in one call
  !! @param field_to_write_information The field that is going to be written
  !! @param timestep The write timestep
  !! @param time The write time
  !! @param file_state File storate state
  subroutine write_independent_variable_to_diagnostics(field_to_write_information, timestep, time, file_state)
    type(writer_field_type), intent(inout) :: field_to_write_information
    integer, intent(in) :: timestep
    real, intent(in) :: time
    type(netcdf_diagnostics_type), intent(inout) :: file_state

    if (field_to_write_information%field_type == ARRAY_FIELD_TYPE .or. &
         field_to_write_information%field_type == SCALAR_FIELD_TYPE) then
      if (field_to_write_information%data_type == DOUBLE_DATA_TYPE .or. &
           field_to_write_information%data_type == INTEGER_DATA_TYPE) then
        call write_out_number_values(field_to_write_information, timestep, time, file_state)
      end if
    else if (field_to_write_information%field_type == MAP_FIELD_TYPE) then
      call write_out_map(field_to_write_information, timestep, time, file_state)
    end if
  end subroutine write_independent_variable_to_diagnostics

  subroutine write_out_number_values(field_to_write_information, timestep, time, file_state)
    type(writer_field_type), intent(inout) :: field_to_write_information
    integer, intent(in) :: timestep
    real, intent(in) :: time
    type(netcdf_diagnostics_type), intent(inout) :: file_state

    real(kind=DEFAULT_PRECISION), dimension(:), allocatable :: values_to_write, timeseries_time_to_write
    real :: value_to_test
    integer :: i, field_id, next_entry_index, array_size, included_num
    integer, dimension(:), allocatable :: count_to_write
    type(iterator_type) :: iterator
    type(mapentry_type) :: map_entry

    type(data_values_type), pointer :: data_value
    character(len=STRING_LENGTH), dimension(:), allocatable :: items_to_remove
    type(netcdf_diagnostics_timeseries_type), pointer :: timeseries_diag

    timeseries_diag=>get_specific_timeseries_dimension(file_state, field_to_write_information%output_frequency, &
         field_to_write_information%timestep_frequency)

    next_entry_index=1
    included_num=1
    array_size=1

    if (.not. timeseries_diag%variable_written) allocate(timeseries_time_to_write(timeseries_diag%num_entries))

    allocate(count_to_write(field_to_write_information%dimensions+1))
    if (field_to_write_information%dimensions .gt. 0) then
      do i=1, field_to_write_information%dimensions
        array_size=array_size*field_to_write_information%actual_dim_size(i)
        count_to_write(i)=field_to_write_information%actual_dim_size(i)
      end do
    end if
    allocate(values_to_write(array_size*timeseries_diag%num_entries))
    allocate(items_to_remove(timeseries_diag%num_entries))
    iterator=c_get_iterator(field_to_write_information%values_to_write)
    do while (c_has_next(iterator))
      map_entry=c_next_mapentry(iterator)
      value_to_test=conv_to_real(map_entry%key)
      if (value_to_test .le. time .and. value_to_test .gt. field_to_write_information%previous_write_time) then
        data_value=>get_data_value_by_field_name(field_to_write_information%values_to_write, map_entry%key)
        if (size(values_to_write) .ge. next_entry_index+size(data_value%values)-1) then
          values_to_write(next_entry_index: next_entry_index+size(data_value%values)-1)=data_value%values(:)
          next_entry_index=next_entry_index+size(data_value%values)
          deallocate(data_value%values)
          deallocate(data_value)
          items_to_remove(included_num)=map_entry%key
          if (allocated(timeseries_time_to_write)) timeseries_time_to_write(included_num)=value_to_test
          included_num=included_num+1
        else
          call log_log(LOG_WARN, "Omitted time entry of field '"//trim(field_to_write_information%field_name)//&
               "' as past time dimension length")
        end if
      end if
    end do
    count_to_write(size(count_to_write))=included_num-1
    field_id=c_get_integer(file_state%variable_to_id, get_field_key(field_to_write_information))
    if (included_num-1 .ne. timeseries_diag%num_entries) then
      call log_log(LOG_WARN, "Miss match of time entries for field '"//trim(field_to_write_information%field_name)//&
           "', included entries="//trim(conv_to_string(included_num-1))//" but expected entries="//&
           trim(conv_to_string(timeseries_diag%num_entries)))
    end if
    call check_thread_status(forthread_mutex_lock(file_state%mutex))
    call check_thread_status(forthread_mutex_lock(netcdf_mutex))
    call lock_mpi()
    call check_netcdf_status(nf90_put_var(file_state%ncid, field_id, values_to_write, count=count_to_write))
    if (allocated(timeseries_time_to_write)) then
      call check_netcdf_status(nf90_put_var(file_state%ncid, timeseries_diag%netcdf_var_id, &
           timeseries_time_to_write, count=(/ timeseries_diag%num_entries /)))
      timeseries_diag%variable_written=.true.
    end if
    call unlock_mpi()
    call check_thread_status(forthread_mutex_unlock(netcdf_mutex))
    call check_thread_status(forthread_mutex_unlock(file_state%mutex))
    deallocate(values_to_write)
    if (included_num .gt. 1) then
      do i=1, included_num-1
        call c_remove(field_to_write_information%values_to_write, items_to_remove(i))
      end do
    end if
    deallocate(items_to_remove)
    if (allocated(timeseries_time_to_write)) deallocate(timeseries_time_to_write)
  end subroutine write_out_number_values

  subroutine write_out_map(field_to_write_information, timestep, time, file_state)
    type(writer_field_type), intent(inout) :: field_to_write_information
    integer, intent(in) :: timestep
    real, intent(in) :: time
    type(netcdf_diagnostics_type), intent(inout) :: file_state

    integer :: i, j, field_id, included_num
    real :: value_to_test
    type(iterator_type) :: iterator, map_data_iterator
    type(mapentry_type) :: map_entry, map_data_entry
    type(data_values_type), pointer :: data_value
    character(len=STRING_LENGTH), dimension(:), allocatable :: items_to_remove

    field_id=c_get_integer(file_state%variable_to_id, get_field_key(field_to_write_information))
    included_num=1
    j=1
    allocate(items_to_remove(c_size(field_to_write_information%values_to_write)))
    iterator=c_get_iterator(field_to_write_information%values_to_write)
    do while (c_has_next(iterator))
      map_entry=c_next_mapentry(iterator)
      value_to_test=conv_to_real(map_entry%key)
      if (value_to_test .le. time .and. value_to_test .gt. field_to_write_information%previous_write_time) then
        items_to_remove(included_num)=map_entry%key
        included_num=included_num+1
        data_value=>get_data_value_by_field_name(field_to_write_information%values_to_write, map_entry%key)
        map_data_iterator=c_get_iterator(data_value%map_values)
        i=1
        call check_thread_status(forthread_mutex_lock(file_state%mutex))
        call check_thread_status(forthread_mutex_lock(netcdf_mutex))
        call lock_mpi()
        do while (c_has_next(map_data_iterator))
          map_data_entry=c_next_mapentry(map_data_iterator)          
          call check_netcdf_status(nf90_put_var(file_state%ncid, field_id, trim(map_data_entry%key), (/ 1, 1, i, j /)))
          call check_netcdf_status(nf90_put_var(file_state%ncid, field_id, trim(c_get_string(map_data_entry)), (/ 1, 2, i, j /)))
          call c_remove(data_value%map_values, map_data_entry%key)
          i=i+1
        end do
        call unlock_mpi()
        call check_thread_status(forthread_mutex_unlock(netcdf_mutex))
        call check_thread_status(forthread_mutex_unlock(file_state%mutex))
        call c_free(data_value%map_values)
        j=j+1
      end if
    end do
    if (included_num .gt. 1) then
      do i=1, included_num-1
        call c_remove(field_to_write_information%values_to_write, items_to_remove(i))
      end do      
    end if
    deallocate(items_to_remove)
  end subroutine write_out_map

  !> Writes out the options_database defining this model run.
  !! @param io_configuration The configuration of the IO server
  subroutine write_out_options(io_configuration, file_state)
    type(io_configuration_type), intent(inout) :: io_configuration
    type(netcdf_diagnostics_type), intent(inout) :: file_state

    integer :: i
    character(len=STRING_LENGTH), pointer :: sized_raw_character
    class(*), pointer :: raw_data, raw_to_string

    call lock_mpi()
    do i=1,options_size(io_configuration%options_database)
      raw_data=> options_value_at(io_configuration%options_database, i)
      raw_to_string=>raw_data
      call check_netcdf_status(nf90_put_var(file_state%ncid, nopt_var_id, &
           trim(options_key_at(io_configuration%options_database, i)), (/ 1, 1, i /)))
      select type (raw_data)
      type is(integer)
        call check_netcdf_status(nf90_put_var(file_state%ncid, nopt_var_id, trim(conv_to_string(raw_data)), (/ 1, 2, i /)))
      type is(real(kind=SINGLE_PRECISION))
        call check_netcdf_status(nf90_put_var(file_state%ncid, nopt_var_id, trim(conv_to_string(raw_data)), (/ 1, 2, i /)))
      type is(real(kind=DOUBLE_PRECISION))
        call check_netcdf_status(nf90_put_var(file_state%ncid, nopt_var_id, trim(conv_to_string(raw_data)), (/ 1, 2, i /)))
      type is(logical)
        call check_netcdf_status(nf90_put_var(file_state%ncid, nopt_var_id, trim(conv_to_string(raw_data)), (/ 1, 2, i /)))
      type is(character(len=*))
        ! Done this way to give the character size information and keep the (unsafe) cast in the conversion module
        sized_raw_character=>conv_to_string(raw_to_string, .false., STRING_LENGTH)
        call check_netcdf_status(nf90_put_var(file_state%ncid, nopt_var_id, trim(sized_raw_character), (/ 1, 2, i /)))
      end select
    end do
    call unlock_mpi()

  end subroutine write_out_options

  !> Defines dimensions for all required dimensions. This is usually the number required plus one, but in some cases is entirely
  !! required depending how the output frequency and diagnostics write times match up
  !! @param file_state The state of the NetCDF file
  !! @param file_writer_information Writer information
  !! @param time The model write time
  subroutine define_time_series_dimensions(file_state, file_writer_information, time, time_points, termination_write)
    type(netcdf_diagnostics_type), intent(inout) :: file_state
    type(writer_type), intent(inout) :: file_writer_information
    real, intent(in) :: time
    type(map_type), intent(inout) :: time_points
    logical, intent(in) :: termination_write

    integer :: i
    character(len=STRING_LENGTH) :: dim_key
    type(netcdf_diagnostics_timeseries_type), pointer :: timeseries_diag
    class(*), pointer :: generic

    do i=1, size(file_writer_information%contents)
      if (file_writer_information%contents(i)%output_frequency .lt. 0.0) then
        dim_key="time_series_"//trim(conv_to_string(file_writer_information%contents(i)%timestep_frequency))
      else
        dim_key="time_series_"//trim(conv_to_string(file_writer_information%contents(i)%timestep_frequency))//"_"//&
             trim(conv_to_string(file_writer_information%contents(i)%output_frequency))
      end if
      if (.not. c_contains(file_state%timeseries_dimension, dim_key)) then
        allocate(timeseries_diag)
        timeseries_diag%variable_written=.false.
        timeseries_diag%num_entries=get_number_timeseries_entries(time_points, &
             file_writer_information%contents(i)%previous_tracked_write_point, &
             file_writer_information%contents(i)%output_frequency, file_writer_information%contents(i)%timestep_frequency, &
             termination_write, timeseries_diag%last_write_point)
        call lock_mpi()
        call check_netcdf_status(nf90_def_dim(file_state%ncid, dim_key, timeseries_diag%num_entries, &
             timeseries_diag%netcdf_dim_id))
        call unlock_mpi()
        generic=>timeseries_diag
        call c_put_generic(file_state%timeseries_dimension, dim_key, generic, .false.)
      end if
      file_writer_information%contents(i)%previous_tracked_write_point=timeseries_diag%last_write_point
    end do
  end subroutine define_time_series_dimensions

  !> Retrieves the number of timeseries entries for a specific frequency and previous write time. This is based on
  !! the range of time points that are provided to the call
  !! @param time_points The list of times that data has been sent over that is applicable to this write
  !! @param previous_write_time When the field was previously written
  !! @param frequency The frequency of outputs of the field
  integer function get_number_timeseries_entries(time_points, previous_write_time, output_frequency, timestep_frequency, &
       termination_write, last_write_entry)
    type(map_type), intent(inout) :: time_points
    real, intent(in) :: output_frequency, previous_write_time
    integer, intent(in) :: timestep_frequency
    logical, intent(in) :: termination_write
    real, intent(out) :: last_write_entry

    integer :: ts
    real :: tp_entry, write_point
    type(iterator_type) :: iterator
    type(mapentry_type) :: map_entry
    logical :: include_item

    get_number_timeseries_entries=0
    write_point=previous_write_time
    iterator=c_get_iterator(time_points)
    do while (c_has_next(iterator))
      map_entry=c_next_mapentry(iterator)
      ts=conv_to_integer(map_entry%key)
      if (timestep_frequency .gt. 0) then
        include_item=mod(ts, timestep_frequency) == 0
      else
        include_item=.false.
      end if
      if (include_item .or. (.not. c_has_next(iterator) .and. termination_write)) then
        tp_entry=c_get_real(map_entry)
        if (tp_entry .ge. write_point+output_frequency) then
          get_number_timeseries_entries=get_number_timeseries_entries+1
          write_point=tp_entry
          last_write_entry=tp_entry
        end if
      end if
    end do
  end function get_number_timeseries_entries

  !> Defines a coordinate variable in the NetCDF file
  !! @param file_state The NetCDF file state
  !! @param coord_name The name of the coordinate
  !! @returns The resulting coordinate variable id
  integer function define_coordinate_variable(file_state, coord_name)
    type(netcdf_diagnostics_type), intent(inout) :: file_state
    character(len=*), intent(in) :: coord_name

    integer :: field_id, dimension_ids(1)

    dimension_ids(1)=c_get_integer(file_state%dimension_to_id, trim(coord_name))
    call lock_mpi()
    call check_netcdf_status(nf90_def_var(file_state%ncid, trim(coord_name), NF90_DOUBLE, dimension_ids, field_id))
    call unlock_mpi()
    define_coordinate_variable=field_id
  end function define_coordinate_variable  

  !> Defines the options_database variable in the NetCDF file
  !! @param file_state The NetCDF file state
  !! @returns The resulting options_database variable id
  integer function define_options_database_variable(file_state)
    type(netcdf_diagnostics_type), intent(inout) :: file_state

    integer :: field_id, dimension_ids(3)

    dimension_ids(1)=file_state%string_dim_id
    dimension_ids(2)=file_state%key_value_dim_id
    dimension_ids(3)=nopt_dim_id
    call lock_mpi()
    call check_netcdf_status(nf90_def_var(file_state%ncid, "options_database", NF90_CHAR, dimension_ids, field_id))
    call unlock_mpi()
    define_options_database_variable=field_id
  end function define_options_database_variable

  !> Defines all variables in the file writer state
  !! @param file_state The NetCDF file state
  !! @param file_writer_information The file writer information
  subroutine define_variables(io_configuration, file_state, file_writer_information)
    type(io_configuration_type), intent(inout) :: io_configuration
    type(netcdf_diagnostics_type), intent(inout) :: file_state
    type(writer_type), intent(in) :: file_writer_information

    integer :: i, j, data_type, field_id, map_dim_id
    integer, dimension(:), allocatable :: dimension_ids
    character(len=STRING_LENGTH) :: variable_key
    type(netcdf_diagnostics_timeseries_type), pointer :: timeseries_diag
    class(*), pointer :: generic
    type(iterator_type) :: iterator
    type(mapentry_type) :: map_entry

    l_nc_dim = .false.
    l_nd_dim = .false.

    iterator=c_get_iterator(file_state%timeseries_dimension)
    do while (c_has_next(iterator))
      map_entry=c_next_mapentry(iterator)
      generic=>c_get_generic(map_entry)
      select type(generic)
        type is(netcdf_diagnostics_timeseries_type)
          timeseries_diag=>generic
      end select
      call lock_mpi()
      call check_netcdf_status(nf90_def_var(file_state%ncid, map_entry%key, &
             NF90_DOUBLE, timeseries_diag%netcdf_dim_id, timeseries_diag%netcdf_var_id))
      call unlock_mpi()
    end do

    do i=1, size(file_writer_information%contents)
      if (.not. file_writer_information%contents(i)%enabled) cycle
      if (file_writer_information%contents(i)%data_type == DOUBLE_DATA_TYPE) then
        data_type=NF90_DOUBLE
      else if (file_writer_information%contents(i)%data_type == INTEGER_DATA_TYPE) then
        data_type=NF90_INT
      else if (file_writer_information%contents(i)%data_type == STRING_DATA_TYPE) then
        data_type=NF90_CHAR
      end if
      variable_key=get_field_key(file_writer_information%contents(i))      
      if (file_writer_information%contents(i)%field_type == ARRAY_FIELD_TYPE) then
        allocate(dimension_ids(file_writer_information%contents(i)%dimensions+1))
        do j=1, file_writer_information%contents(i)%dimensions
          if (c_contains(file_state%dimension_to_id, file_writer_information%contents(i)%dim_size_defns(j))) then
            dimension_ids(j)=c_get_integer(file_state%dimension_to_id, file_writer_information%contents(i)%dim_size_defns(j))
          else
            call log_log(LOG_ERROR, "Can not find information for dimension named '"//&
                 trim(file_writer_information%contents(i)%dim_size_defns(j))//"'")
          end if
        end do
        dimension_ids(j)=retrieve_time_series_dimension_id_for_field(file_state, file_writer_information, i)

        !> Only enable writing of conditional diagnostics descriptor fields to this file if there is
        !  a data field using its dimensions.
        if (any(dimension_ids .eq. nc_dim_id) .or. any(dimension_ids .eq. nd_dim_id)) then
          l_nc_dim = .true.
          l_nd_dim = .true.
        end if 

        call lock_mpi()
        call check_netcdf_status(nf90_def_var(file_state%ncid, variable_key, &
             data_type, dimension_ids, field_id))

        if (file_writer_information%contents(i)%collective_write .and. &
             file_writer_information%contents(i)%collective_contiguous_optimisation .and. &
             io_configuration%number_of_io_servers .gt. 1) then
          call check_netcdf_status(nf90_def_var_fill(file_state%ncid, field_id, 1, 1))
          call check_netcdf_status(nf90_var_par_access(file_state%ncid, field_id, NF90_COLLECTIVE))
        end if        
        call unlock_mpi()
        deallocate(dimension_ids)
      else if (file_writer_information%contents(i)%field_type == SCALAR_FIELD_TYPE) then
        call lock_mpi()
        call check_netcdf_status(nf90_def_var(file_state%ncid, variable_key, &
             data_type, retrieve_time_series_dimension_id_for_field(file_state, file_writer_information, i), field_id))
        call unlock_mpi()
      else if (file_writer_information%contents(i)%field_type == MAP_FIELD_TYPE) then        
        allocate(dimension_ids(4))
        dimension_ids(1)=file_state%string_dim_id
        dimension_ids(2)=file_state%key_value_dim_id
        dimension_ids(3)=c_get_integer(file_state%dimension_to_id, file_writer_information%contents(i)%dim_size_defns(1))
        dimension_ids(4)=retrieve_time_series_dimension_id_for_field(file_state, file_writer_information, i)
        call lock_mpi()
        call check_netcdf_status(nf90_def_var(file_state%ncid, variable_key, data_type, dimension_ids, field_id))
        call unlock_mpi()
        deallocate(dimension_ids)
      end if
      call c_put_integer(file_state%variable_to_id, variable_key, field_id)
      if (len_trim(file_writer_information%contents(i)%units) .gt. 0) then
        call lock_mpi()
        call check_netcdf_status(nf90_put_att(file_state%ncid, field_id, "units", file_writer_information%contents(i)%units))
        call unlock_mpi()
      end if
    end do

    !> Provide conditional diagnostic descriptions if data being used in this file.
    if (l_nc_dim) then
      allocate(dimension_ids(2))
      dimension_ids(1)=file_state%string_dim_id
      dimension_ids(2)=nc_dim_id
      call lock_mpi()
      call check_netcdf_status(nf90_def_var(file_state%ncid, "conditions_fields_short", &
             NF90_CHAR, dimension_ids, nc_var_id_s))
      call check_netcdf_status(nf90_def_var(file_state%ncid, "conditions_fields_long", &
             NF90_CHAR, dimension_ids, nc_var_id_l))
      call unlock_mpi()
      deallocate(dimension_ids)
    end if
    if (l_nd_dim) then
      allocate(dimension_ids(2))
      dimension_ids(1)=file_state%string_dim_id
      dimension_ids(2)=nd_dim_id
      call lock_mpi()
      call check_netcdf_status(nf90_def_var(file_state%ncid, "diagnostics_fields_short", &
             NF90_CHAR, dimension_ids, nd_var_id_s))
      call check_netcdf_status(nf90_def_var(file_state%ncid, "diagnostics_fields_long", &
             NF90_CHAR, dimension_ids, nd_var_id_l))
      call unlock_mpi()
      deallocate(dimension_ids)
    end if

  end subroutine define_variables


  !> For a specific field will retrieve the NetCDF id of the time series dimension most appropriate for this field. If
  !! a dimension can not be located then an error is raised
  !! @param file_state NetCDF file state
  !! @param file_writer_information The writer information
  !! @param field_index Index of the field in the writer information contents that we are writing
  !! @returns The NetCDF dimension id
  integer function retrieve_time_series_dimension_id_for_field(file_state, file_writer_information, field_index)
    type(netcdf_diagnostics_type), intent(inout) :: file_state
    type(writer_type), intent(in) :: file_writer_information
    integer, intent(in) :: field_index

    type(netcdf_diagnostics_timeseries_type), pointer :: timeseries_diag

    timeseries_diag=>get_specific_timeseries_dimension(file_state, &
         file_writer_information%contents(field_index)%output_frequency, &
         file_writer_information%contents(field_index)%timestep_frequency)
    if (associated(timeseries_diag)) then
      retrieve_time_series_dimension_id_for_field=timeseries_diag%netcdf_dim_id
    else
      call log_log(LOG_ERROR, "Can not find time series dimension with output frequency "//&
           trim(conv_to_string(file_writer_information%contents(field_index)%output_frequency)))
    end if
  end function retrieve_time_series_dimension_id_for_field

  !> Given the file state and the output frequency of a field will retrive the appropriate time series dimension entry
  !! that corresponds to this or null if none can be found
  !! @param file_state The NetCDF file state that is being written
  !! @param output_frequency Time frequency of writes
  !! @param timestep_frequency The timestep frequency of data arrival/generation
  !! @returns The corresponding time series dimension entry or null if none is found
  function get_specific_timeseries_dimension(file_state, output_frequency, timestep_frequency)
    type(netcdf_diagnostics_type), intent(inout) :: file_state
    real, intent(in) :: output_frequency
    integer, intent(in) :: timestep_frequency
    type(netcdf_diagnostics_timeseries_type), pointer :: get_specific_timeseries_dimension

    character(len=STRING_LENGTH) :: dim_key
    class(*), pointer :: generic

    if (output_frequency .lt. 0.0) then
      dim_key="time_series_"//trim(conv_to_string(timestep_frequency))
    else
      dim_key="time_series_"//trim(conv_to_string(timestep_frequency))//"_"// trim(conv_to_string(output_frequency))
    end if
    generic=>c_get_generic(file_state%timeseries_dimension, dim_key)

    if (associated(generic)) then
      select type(generic)
      type is(netcdf_diagnostics_timeseries_type)
        get_specific_timeseries_dimension=>generic
      end select
    else
      get_specific_timeseries_dimension=>null()
    end if
  end function get_specific_timeseries_dimension  

  !> Defines spatial dimensions in the diagnostics file
  !! @param file_state The NetCDF file state
  !! @param dimension_sizing All dimension name to size key value pairs
  subroutine define_dimensions(file_state, dimension_sizing)
    type(netcdf_diagnostics_type), intent(inout) :: file_state
    type(map_type), intent(inout) :: dimension_sizing

    integer :: ncdf_dimid, dim_length
    type(iterator_type) :: iterator
    type(mapentry_type) :: map_entry


    iterator=c_get_iterator(dimension_sizing)
    call lock_mpi()
    do while (c_has_next(iterator))
      map_entry=c_next_mapentry(iterator)
      dim_length=c_get_integer(map_entry)
      if (dim_length .gt. 0) then
        call check_netcdf_status(nf90_def_dim(file_state%ncid, map_entry%key, dim_length, ncdf_dimid))
        call c_put_integer(file_state%dimension_to_id, map_entry%key, ncdf_dimid)
        if (map_entry%key == "nc") then
          nc_dim_id=ncdf_dimid
        end if
        if (map_entry%key == "nd") then
          nd_dim_id=ncdf_dimid
        end if
        if (map_entry%key == "number_options") then
          nopt_dim_id=ncdf_dimid
        end if
      end if
    end do
    call check_netcdf_status(nf90_def_dim(file_state%ncid, "string", STRING_LENGTH, file_state%string_dim_id))
    call check_netcdf_status(nf90_def_dim(file_state%ncid, "kvp", 2, file_state%key_value_dim_id))
    call unlock_mpi()
  end subroutine define_dimensions

  !> Retrieves a file state based upon its timestep or null if none is found
  !! @param filename The filename to look up
  !! @param timestep The timestep to look up
  !! @param dolock Whether to issue a read lock or not
  !! @returns The corresponding file state entry or null if none is found  
  function get_file_state(filename, timestep, dolock)
    character(len=*), intent(in) :: filename
    integer, intent(in) :: timestep
    logical, intent(in) :: dolock
    type(netcdf_diagnostics_type), pointer :: get_file_state

    class(*), pointer :: generic

    if (dolock) call check_thread_status(forthread_rwlock_rdlock(file_states_rwlock))
    generic=>c_get_generic(file_states, trim(filename)//"#"//trim(conv_to_string(timestep)))
    if (dolock) call check_thread_status(forthread_rwlock_unlock(file_states_rwlock))

    if (associated(generic)) then
      select type(generic)
      type is (netcdf_diagnostics_type)      
        get_file_state=>generic
      end select
    else
      get_file_state=>null()
    end if
  end function get_file_state

  !> Retrieves the field key, corresponding to the field name in the NetCDF file and what we store against the id in
  !! the NetCDF file state internally. This is either the field name itself or the manipulation type appended if there is
  !! duplication
  !! @param field_to_write_information The field to get the key from
  !! @returns The field name key
  character(len=STRING_LENGTH) function get_field_key(field_to_write_information)
    type(writer_field_type), intent(in) :: field_to_write_information    
    
    get_field_key=field_to_write_information%field_name
    if (field_to_write_information%duplicate_field_name) then
      if (field_to_write_information%time_manipulation_type == INSTANTANEOUS_TYPE) then
        get_field_key=trim(get_field_key)//"_instantaneous"
      else if (field_to_write_information%time_manipulation_type == TIME_AVERAGED_TYPE) then
        get_field_key=trim(get_field_key)//"_timeaveraged"
      end if
    end if
  end function get_field_key
  
  !> Generates a unique filename based upon the base one specified and the number 
  !! of completed timesteps
  !! @param old_name The existing name that is used as a base
  !! @param timestep The current model timestep
  !! @param new_name The new name that is produced by this subroutine
  subroutine generate_unique_filename(old_name, new_name, configured_write_time, timestep)
    character(len=STRING_LENGTH), intent(in) :: old_name
    real, intent(in), optional :: configured_write_time
    integer, intent(in), optional :: timestep
    character(len=STRING_LENGTH), intent(out) :: new_name

    integer :: dot_posn
    
    dot_posn=index(old_name, ".")
    if (dot_posn .gt. 0) then
      new_name = old_name(1:dot_posn-1)
    else
      new_name=old_name
    end if
    if (present(configured_write_time)) then
      new_name=trim(new_name)//"_"//trim(conv_to_string(configured_write_time))
    else if (present(timestep)) then
      new_name=trim(new_name)//"_"//trim(conv_to_string(timestep))
    end if
    if (dot_posn .gt. 0) then
      new_name=trim(new_name)//old_name(dot_posn:len(old_name))
    end if    
  end subroutine generate_unique_filename

  !> Writes out global attributes into the checkpoint
  !! @param ncid NetCDF file id
  subroutine write_out_global_attributes(io_configuration, ncid, file_writer_information, timestep, time)
    type(io_configuration_type), intent(inout) :: io_configuration
    integer, intent(in) :: ncid, timestep
    type(writer_type), intent(inout) :: file_writer_information
    real, intent(in) :: time

    integer :: date_values(8), ierr
    character(len=50) :: date_time

    call date_and_time(values=date_values)
    call lock_mpi()
    call mpi_bcast(date_values, 8, MPI_INT, 0, io_configuration%io_communicator, ierr)
    call unlock_mpi()
    date_time=trim(conv_to_string(date_values(3)))//"/"//&
         trim(conv_to_string(date_values(2)))//"/"//trim(conv_to_string(date_values(1)))//" "//trim(conv_to_string(&
         date_values(5)))// ":"//trim(conv_to_string(date_values(6)))//":"//trim(conv_to_string(date_values(7)))

    call lock_mpi()
    call check_netcdf_status(nf90_put_att(ncid, NF90_GLOBAL, "title", file_writer_information%title))
    call check_netcdf_status(nf90_put_att(ncid, NF90_GLOBAL, "created", date_time))
    call check_netcdf_status(nf90_put_att(ncid, NF90_GLOBAL, "MONC time", trim(conv_to_string(time))))
    call check_netcdf_status(nf90_put_att(ncid, NF90_GLOBAL, "MONC timestep", trim(conv_to_string(timestep))))
    call check_netcdf_status(nf90_put_att(ncid, NF90_GLOBAL, "Diagnostic write frequency", &
         trim(conv_to_string(file_writer_information%write_time_frequency))))
    call check_netcdf_status(nf90_put_att(ncid, NF90_GLOBAL, "Previous diagnostic write at", &
         trim(conv_to_string(file_writer_information%previous_write_time))))
    call unlock_mpi()
  end subroutine write_out_global_attributes
end module netcdf_filetype_writer_mod
