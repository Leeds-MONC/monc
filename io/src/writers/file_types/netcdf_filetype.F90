!> The NetCDF file type writer which performs actual writing of NetCDF files to the parallel filesystem. These are opened by
!! all IO servers and all IO servers can participate as variables might be located across the different IO processes
module netcdf_filetype_writer_mod
  use datadefn_mod, only : DEFAULT_PRECISION, STRING_LENGTH
  use configuration_parser_mod, only : INSTANTANEOUS_TYPE, TIME_AVERAGED_TYPE, io_configuration_type, &
       data_values_type, get_data_value_by_field_name
  use collections_mod, only : hashmap_type, hashset_type, list_type, map_type, c_get, c_contains, &
       c_key_at, c_value_at, c_put, c_size, c_add, c_remove, c_free
  use conversions_mod, only : conv_to_generic, conv_to_integer, conv_to_string, conv_to_real
  use logging_mod, only : LOG_ERROR, LOG_WARN, LOG_DEBUG, log_log, log_master_log, log_get_logging_level, log_is_master
  use writer_types_mod, only : writer_type, writer_field_type, write_field_collective_values_type
  use forthread_mod, only : forthread_mutex_init, forthread_mutex_lock, forthread_mutex_unlock, forthread_mutex_destroy, &
       forthread_rwlock_rdlock, forthread_rwlock_wrlock, forthread_rwlock_unlock, forthread_rwlock_init, forthread_rwlock_destroy
  use threadpool_mod, only : check_thread_status
  use netcdf, only : NF90_DOUBLE, NF90_REAL, NF90_INT, NF90_CHAR, NF90_GLOBAL, NF90_CLOBBER, NF90_NETCDF4, NF90_MPIIO, &
       NF90_COLLECTIVE, NF90_UNLIMITED, nf90_def_var, nf90_var_par_access, nf90_def_var_fill, nf90_put_att, &
       nf90_create, nf90_put_var, nf90_def_dim, nf90_enddef, nf90_close, nf90_ebaddim, nf90_enotatt, nf90_enotvar, &
       nf90_noerr, nf90_strerror, nf90_redef, nf90_inq_varid
  use io_server_client_mod, only : ARRAY_FIELD_TYPE, SCALAR_FIELD_TYPE, DOUBLE_DATA_TYPE, INTEGER_DATA_TYPE  
  use mpi, only : MPI_INFO_NULL
  use grids_mod, only : Z_INDEX, Y_INDEX, X_INDEX
  implicit none

#ifndef TEST_MODE
  private
#endif

  type netcdf_diagnostics_timeseries_type
     integer :: netcdf_dim_id, netcdf_var_id, num_entries
     real :: last_write_point
     logical :: variable_written
  end type netcdf_diagnostics_timeseries_type  

  !< Keeps track of a specific diagnostics NetCDF file
  type netcdf_diagnostics_type
     integer :: ncid, mutex
     type(map_type) :: dimension_to_id
     type(hashmap_type) :: variable_to_id, timeseries_dimension
     type(writer_type), pointer :: corresponding_writer_entry
  end type netcdf_diagnostics_type

  type(hashmap_type), volatile :: file_states
  integer, volatile :: file_states_rwlock

  public initialise_netcdf_filetype, finalise_netcdf_filetype, define_netcdf_file, write_variable, close_netcdf_file
contains

  !> Initialises the NetCDF writing functionality
  subroutine initialise_netcdf_filetype()
    call check_thread_status(forthread_rwlock_init(file_states_rwlock, -1))
  end subroutine initialise_netcdf_filetype

  !> Finalises  the NetCDF writing functionality
  subroutine finalise_netcdf_filetype()
    call check_thread_status(forthread_rwlock_destroy(file_states_rwlock))
  end subroutine finalise_netcdf_filetype  
  
  !> Defines a NetCDF file - which creates it, defines all dimensions and variables. This must be called by all IO server
  !! processes as the NetCDF operations here are collective
  !! @param io_configuration The IO server configuration
  !! @param file_writer_information The writer entry that is being written
  !! @param timestep The write timestep
  !! @param time The write time
  subroutine define_netcdf_file(io_configuration, file_writer_information, timestep, time, time_points)
    type(io_configuration_type), intent(inout) :: io_configuration
    type(writer_type), intent(inout), target :: file_writer_information
    type(map_type), intent(inout) :: time_points
    integer :: timestep
    real :: time
    
    character(len=STRING_LENGTH) :: unique_filename
    type(netcdf_diagnostics_type), pointer :: ncdf_writer_state
    class(*), pointer :: generic

    ncdf_writer_state=>get_file_state(file_writer_information%filename, timestep, .true.)
    if (.not. associated(ncdf_writer_state)) then
      call check_thread_status(forthread_rwlock_wrlock(file_states_rwlock))
      ncdf_writer_state=>get_file_state(file_writer_information%filename, timestep, .false.)
      if (.not. associated(ncdf_writer_state)) then
        allocate(ncdf_writer_state)
        ncdf_writer_state%corresponding_writer_entry=>file_writer_information
        call check_thread_status(forthread_mutex_init(ncdf_writer_state%mutex, -1))
        call check_thread_status(forthread_mutex_lock(ncdf_writer_state%mutex))
        generic=>ncdf_writer_state
        call c_put(file_states, trim(file_writer_information%filename)//"#"//trim(conv_to_string(timestep)), generic)
        call check_thread_status(forthread_rwlock_unlock(file_states_rwlock))
        
        call generate_unique_filename(file_writer_information%filename, timestep, unique_filename)
        if (io_configuration%number_of_io_servers .gt. 1) then
          call check_status(nf90_create(unique_filename, ior(NF90_NETCDF4, NF90_MPIIO), ncdf_writer_state%ncid, &
               comm = io_configuration%io_communicator, info = MPI_INFO_NULL))
        else
          call check_status(nf90_create(unique_filename, NF90_CLOBBER, ncdf_writer_state%ncid))
        end if
        call write_out_global_attributes(ncdf_writer_state%ncid, file_writer_information, timestep, time)
        call define_dimensions(ncdf_writer_state, io_configuration%dimension_sizing)
        call define_time_series_dimensions(ncdf_writer_state, file_writer_information, time, time_points)
        call define_variables(ncdf_writer_state, file_writer_information)
        call check_status(nf90_enddef(ncdf_writer_state%ncid))
        call check_thread_status(forthread_mutex_unlock(ncdf_writer_state%mutex))
      else
        call check_thread_status(forthread_rwlock_unlock(file_states_rwlock))
      end if
    end if
  end subroutine define_netcdf_file

  !> Call back for the inter IO reduction which actually does the NetCDF file closing which is a 
  !! collective (synchronous) operation. This also cleans up the file state as it is no longer required
  !! @param io_configuration The IO server configuration
  !! @param field_name The field name that is being communicated
  !! @param timestep The write timestep
  function close_netcdf_file(io_configuration, field_name, timestep)
    type(io_configuration_type), intent(inout) :: io_configuration
    character(len=STRING_LENGTH) :: field_name
    integer :: timestep
    type(writer_type), pointer :: close_netcdf_file

    integer :: i, entries
    class(*), pointer :: generic

    type(netcdf_diagnostics_type), pointer :: file_state

    file_state=>get_file_state(field_name, timestep, .true.)
    call check_thread_status(forthread_mutex_lock(file_state%mutex))
    call check_status(nf90_close(file_state%ncid))
    call check_thread_status(forthread_mutex_unlock(file_state%mutex))
    call check_thread_status(forthread_mutex_destroy(file_state%mutex))
    call c_free(file_state%dimension_to_id)
    call c_free(file_state%variable_to_id)
    entries=c_size(file_state%timeseries_dimension)
    do i=1, entries
      generic=>c_value_at(file_state%timeseries_dimension, i)
      deallocate(generic)
    end do    
    call c_free(file_state%timeseries_dimension)
    call check_thread_status(forthread_rwlock_wrlock(file_states_rwlock))
    call c_remove(file_states, trim(field_name)//"#"//trim(conv_to_string(timestep)))
    call check_thread_status(forthread_rwlock_unlock(file_states_rwlock))    
    if (log_get_logging_level() .ge. LOG_DEBUG .and. log_is_master()) then
      call log_master_log(LOG_DEBUG, "Done physical close for NetCDF file at timestep "//trim(conv_to_string(timestep)))
    end if    
    close_netcdf_file=>file_state%corresponding_writer_entry
  end function close_netcdf_file 

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
      call write_collective_variable_to_diagnostics(io_configuration, field_to_write_information, timestep, time, file_state)
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

    get_dimension_original_size=conv_to_integer(c_get(dimension_store, dim_name(:dash_idx)), .false.)
  end function get_dimension_original_size  

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
    integer :: i, j, k, entries, included_num, monc_entries, source, field_id, start(field_to_write_information%dimensions+1), &
         count(field_to_write_information%dimensions+1), monc_location, dim_identifier, auto_period, dim_start
    class(*), pointer :: generic
    type(write_field_collective_values_type), pointer :: multi_monc_entries
    logical :: is_auto_dimension
    real(kind=DEFAULT_PRECISION), dimension(:), allocatable :: timeseries_time_to_write
    character(len=STRING_LENGTH), dimension(:), allocatable :: items_to_remove
    type(data_values_type), pointer :: data_value
    type(netcdf_diagnostics_timeseries_type), pointer :: timeseries_diag

    timeseries_diag=>get_specific_timeseries_dimension(file_state, field_to_write_information%output_frequency, &
         field_to_write_information%timestep_frequency)
    if (.not. timeseries_diag%variable_written) allocate(timeseries_time_to_write(timeseries_diag%num_entries))

    allocate(items_to_remove(timeseries_diag%num_entries))
    included_num=1
    field_id=conv_to_integer(c_get(file_state%variable_to_id, get_field_key(field_to_write_information)), .false.)
    entries=c_size(field_to_write_information%values_to_write)
    do i=1, entries
      value_to_test=conv_to_real(c_key_at(field_to_write_information%values_to_write, i))
      if (value_to_test .le. time .and. value_to_test .gt. field_to_write_information%previous_write_time) then
        if (included_num .le. timeseries_diag%num_entries) then
          if (allocated(timeseries_time_to_write)) timeseries_time_to_write(included_num)=value_to_test
          generic=>c_get(field_to_write_information%values_to_write, c_key_at(field_to_write_information%values_to_write, i))
          select type(generic)
          type is(write_field_collective_values_type)
            multi_monc_entries=>generic
          end select
          monc_entries=c_size(multi_monc_entries%monc_values)
          do j=1, monc_entries
            source=conv_to_integer(c_key_at(multi_monc_entries%monc_values, j))
            data_value=>get_data_value_by_field_name(multi_monc_entries%monc_values, conv_to_string(source))           
            generic=>c_get(io_configuration%monc_to_index, conv_to_string(source))
            monc_location=conv_to_integer(generic, .false.)
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
            call check_status(nf90_put_var(file_state%ncid, field_id, data_value%values, start=start, count=count))
            call check_thread_status(forthread_mutex_unlock(file_state%mutex))
            deallocate(data_value%values)
            deallocate(data_value)
          end do
          items_to_remove(included_num)=c_key_at(field_to_write_information%values_to_write, i)
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
      call check_status(nf90_put_var(file_state%ncid, timeseries_diag%netcdf_var_id, &
           timeseries_time_to_write, count=(/ timeseries_diag%num_entries /)))
      timeseries_diag%variable_written=.true.
    end if
    if (included_num .gt. 1) then
      do i=1, included_num-1
        call c_remove(field_to_write_information%values_to_write, items_to_remove(i))
      end do
      deallocate(items_to_remove)
    end if
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

    real(kind=DEFAULT_PRECISION), dimension(:), allocatable :: values_to_write, timeseries_time_to_write
    real :: value_to_test
    integer :: i, field_id, entries, next_entry_index, array_size, included_num
    integer, dimension(:), allocatable :: count_to_write

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
    entries=c_size(field_to_write_information%values_to_write)
    do i=1, entries
      value_to_test=conv_to_real(c_key_at(field_to_write_information%values_to_write, i))
      if (value_to_test .le. time .and. value_to_test .gt. field_to_write_information%previous_write_time) then
        data_value=>get_data_value_by_field_name(field_to_write_information%values_to_write, &
             c_key_at(field_to_write_information%values_to_write, i))
        if (size(values_to_write) .ge. next_entry_index+size(data_value%values)-1) then
          values_to_write(next_entry_index: next_entry_index+size(data_value%values)-1)=data_value%values(:)
          next_entry_index=next_entry_index+size(data_value%values)
          deallocate(data_value%values)
          deallocate(data_value)
          items_to_remove(included_num)=c_key_at(field_to_write_information%values_to_write, i)
          if (allocated(timeseries_time_to_write)) timeseries_time_to_write(included_num)=value_to_test
          included_num=included_num+1
        else
          call log_log(LOG_WARN, "Omitted time entry of field '"//trim(field_to_write_information%field_name)//&
               "' as past time dimension length")
        end if
      end if
    end do
    count_to_write(size(count_to_write))=included_num-1
    field_id=conv_to_integer(c_get(file_state%variable_to_id, get_field_key(field_to_write_information)), .false.)
    if (included_num-1 .ne. timeseries_diag%num_entries) then
      call log_log(LOG_WARN, "Miss match of time entries for field '"//trim(field_to_write_information%field_name)//&
           "', included entries="//trim(conv_to_string(included_num-1))//" but expected entries="//&
           trim(conv_to_string(timeseries_diag%num_entries)))
    end if
    call check_thread_status(forthread_mutex_lock(file_state%mutex))
    call check_status(nf90_put_var(file_state%ncid, field_id, values_to_write, count=count_to_write))
    if (allocated(timeseries_time_to_write)) then
      call check_status(nf90_put_var(file_state%ncid, timeseries_diag%netcdf_var_id, &
           timeseries_time_to_write, count=(/ timeseries_diag%num_entries /)))
      timeseries_diag%variable_written=.true.
    end if    
    call check_thread_status(forthread_mutex_unlock(file_state%mutex))
    deallocate(values_to_write)
    if (included_num .gt. 1) then
      do i=1, included_num-1
        call c_remove(field_to_write_information%values_to_write, items_to_remove(i))
      end do
      deallocate(items_to_remove)
    end if
    if (allocated(timeseries_time_to_write)) deallocate(timeseries_time_to_write)
  end subroutine write_independent_variable_to_diagnostics

  !> Defines dimensions for all required dimensions. This is usually the number required plus one, but in some cases is entirely
  !! required depending how the output frequency and diagnostics write times match up
  !! @param file_state The state of the NetCDF file
  !! @param file_writer_information Writer information
  !! @param time The model write time
  subroutine define_time_series_dimensions(file_state, file_writer_information, time, time_points)
    type(netcdf_diagnostics_type), intent(inout) :: file_state
    type(writer_type), intent(inout) :: file_writer_information
    real, intent(in) :: time
    type(map_type), intent(inout) :: time_points

    integer :: i
    character(len=STRING_LENGTH) :: dim_key
    type(netcdf_diagnostics_timeseries_type), pointer :: timeseries_diag
    class(*), pointer :: generic

    do i=1, size(file_writer_information%contents)
      dim_key="time_series_"//trim(conv_to_string(file_writer_information%contents(i)%timestep_frequency))//"_"//&
           trim(conv_to_string(file_writer_information%contents(i)%output_frequency))
      if (.not. c_contains(file_state%timeseries_dimension, dim_key)) then
        allocate(timeseries_diag)
        timeseries_diag%variable_written=.false.
        timeseries_diag%num_entries=get_number_timeseries_entries(time_points, &
             file_writer_information%contents(i)%previous_tracked_write_point, &
             file_writer_information%contents(i)%output_frequency, file_writer_information%contents(i)%timestep_frequency, &
             timeseries_diag%last_write_point)
        call check_status(nf90_def_dim(file_state%ncid, dim_key, timeseries_diag%num_entries, timeseries_diag%netcdf_dim_id))        
        generic=>timeseries_diag
        call c_put(file_state%timeseries_dimension, dim_key, generic)
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
       last_write_entry)
    type(map_type), intent(inout) :: time_points
    real, intent(in) :: output_frequency, previous_write_time
    integer, intent(in) :: timestep_frequency
    real, intent(out) :: last_write_entry

    integer :: i, entries, ts
    real :: tp_entry, write_point

    get_number_timeseries_entries=0
    write_point=previous_write_time
    entries=c_size(time_points)
    do i=1, entries
      ts=conv_to_integer(c_key_at(time_points, i))
      if (mod(ts, timestep_frequency) == 0) then
        tp_entry=conv_to_real(c_value_at(time_points, i), .false.)
        if (tp_entry .ge. write_point+output_frequency) then
          get_number_timeseries_entries=get_number_timeseries_entries+1
          write_point=tp_entry
          last_write_entry=tp_entry
        end if
      end if
    end do
  end function get_number_timeseries_entries

  !> Defines all variables in the file writer state
  !! @param file_state The NetCDF file state
  !! @param file_writer_information The file writer information
  subroutine define_variables(file_state, file_writer_information)
    type(netcdf_diagnostics_type), intent(inout) :: file_state
    type(writer_type), intent(in) :: file_writer_information

    integer :: i, j, data_type, field_id, entries
    integer, dimension(:), allocatable :: dimension_ids
    character(len=STRING_LENGTH) :: variable_key
    type(netcdf_diagnostics_timeseries_type), pointer :: timeseries_diag
    class(*), pointer :: generic

    entries=c_size(file_state%timeseries_dimension)
    do i=1, entries
      generic=>c_value_at(file_state%timeseries_dimension, i)
      select type(generic)
        type is(netcdf_diagnostics_timeseries_type)
          timeseries_diag=>generic
      end select   
      variable_key=c_key_at(file_state%timeseries_dimension, i)
      call check_status(nf90_def_var(file_state%ncid, variable_key, &
             NF90_DOUBLE, timeseries_diag%netcdf_dim_id, timeseries_diag%netcdf_var_id))
    end do

    do i=1, size(file_writer_information%contents)
      if (.not. file_writer_information%contents(i)%enabled) cycle
      if (file_writer_information%contents(i)%data_type == DOUBLE_DATA_TYPE) then
        data_type=NF90_DOUBLE
      else if (file_writer_information%contents(i)%data_type == INTEGER_DATA_TYPE) then
        data_type=NF90_INT
      end if
      variable_key=get_field_key(file_writer_information%contents(i))      
      if (file_writer_information%contents(i)%field_type == ARRAY_FIELD_TYPE) then
        allocate(dimension_ids(file_writer_information%contents(i)%dimensions+1))        
        do j=1, file_writer_information%contents(i)%dimensions
          if (c_contains(file_state%dimension_to_id, file_writer_information%contents(i)%dim_size_defns(j))) then
          dimension_ids(j)=conv_to_integer(c_get(file_state%dimension_to_id, &
               file_writer_information%contents(i)%dim_size_defns(j)), .false.)
          else
            call log_log(LOG_ERROR, "Can not find information for dimension named '"//&
                 trim(file_writer_information%contents(i)%dim_size_defns(j))//"'")
          end if          
        end do
        dimension_ids(j)=retrieve_time_series_dimension_id_for_field(file_state, file_writer_information, i)
        call check_status(nf90_def_var(file_state%ncid, variable_key, &
             data_type, dimension_ids, field_id))
        deallocate(dimension_ids)
      else if (file_writer_information%contents(i)%field_type == SCALAR_FIELD_TYPE) then
        call check_status(nf90_def_var(file_state%ncid, variable_key, &
             data_type, retrieve_time_series_dimension_id_for_field(file_state, file_writer_information, i), field_id))
      end if
      call c_put(file_state%variable_to_id, variable_key, conv_to_generic(field_id, .true.))
      if (len_trim(file_writer_information%contents(i)%units) .gt. 0) then
        call check_status(nf90_put_att(file_state%ncid, field_id, "units", file_writer_information%contents(i)%units))
      end if
    end do
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

    dim_key="time_series_"//trim(conv_to_string(timestep_frequency))//"_"// trim(conv_to_string(output_frequency))
    generic=>c_get(file_state%timeseries_dimension, dim_key)

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

    integer :: ncdf_dimid, i, entries
    character(len=STRING_LENGTH) :: dim_str

    entries=c_size(dimension_sizing)

    do i=1, entries
      dim_str=c_key_at(dimension_sizing, i)
      call check_status(nf90_def_dim(file_state%ncid, dim_str, &
           conv_to_integer(c_value_at(dimension_sizing, i), .false.), ncdf_dimid))
      call c_put(file_state%dimension_to_id, dim_str, conv_to_generic(ncdf_dimid, .true.))
    end do
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
    generic=>c_get(file_states, trim(filename)//"#"//trim(conv_to_string(timestep)))
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
  subroutine generate_unique_filename(old_name, timestep, new_name)
    character(len=STRING_LENGTH), intent(in) :: old_name
    integer, intent(in) :: timestep
    character(len=STRING_LENGTH), intent(out) :: new_name

    integer :: dot_posn
    
    dot_posn=index(old_name, ".")
    if (dot_posn .gt. 0) then
      new_name = old_name(1:dot_posn-1)
    else
      new_name=old_name
    end if
    new_name=trim(new_name)//"_"//trim(conv_to_string(timestep))
    if (dot_posn .gt. 0) then
      new_name=trim(new_name)//old_name(dot_posn:len(old_name))
    end if    
  end subroutine generate_unique_filename

  !> Writes out global attributes into the checkpoint
  !! @param ncid NetCDF file id
  subroutine write_out_global_attributes(ncid, file_writer_information, timestep, time)
    integer, intent(in) :: ncid, timestep
    type(writer_type), intent(inout) :: file_writer_information
    real, intent(in) :: time

    integer :: date_values(8)

    call date_and_time(values=date_values)

    call check_status(nf90_put_att(ncid, NF90_GLOBAL, "title", "MONC diagnostics"))
    call check_status(nf90_put_att(ncid, NF90_GLOBAL, "created", trim(conv_to_string(date_values(3)))//"/"//&
         trim(conv_to_string(date_values(2)))//"/"//trim(conv_to_string(date_values(1)))//" "//trim(conv_to_string(&
         date_values(5)))// ":"//trim(conv_to_string(date_values(6)))//":"//trim(conv_to_string(date_values(7)))))
    call check_status(nf90_put_att(ncid, NF90_GLOBAL, "MONC time", trim(conv_to_string(time))))
    call check_status(nf90_put_att(ncid, NF90_GLOBAL, "MONC timestep", trim(conv_to_string(timestep))))
    call check_status(nf90_put_att(ncid, NF90_GLOBAL, "Diagnostic write frequency", &
         trim(conv_to_string(file_writer_information%write_time_frequency))))
    call check_status(nf90_put_att(ncid, NF90_GLOBAL, "Previous diagnostic write at", &
         trim(conv_to_string(file_writer_information%previous_write_time))))
  end subroutine write_out_global_attributes

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
end module netcdf_filetype_writer_mod
