!> The IO server state module which will write out the current state of the IO server to a NetCDF file
module io_server_state_writer_mod
  use datadefn_mod, only : DEFAULT_PRECISION
  use netcdf, only : NF90_CHAR, NF90_BYTE, NF90_INT, nf90_def_var, nf90_put_var, nf90_def_dim
  use configuration_parser_mod, only : io_configuration_type
  use writer_types_mod, only : netcdf_diagnostics_type, serialise_writer_type, writer_type
  use collections_mod, only : hashmap_type, mapentry_type, iterator_type, c_put_integer, c_get_integer, c_get_iterator, &
       c_next_mapentry, c_has_next, c_get_real, c_size
  use conversions_mod, only : conv_to_integer
  use logging_mod, only : LOG_ERROR, log_log
  use io_server_client_mod, only : pack_scalar_field
  use mpi, only : MPI_IN_PLACE, MPI_INT, MPI_SUM, MPI_STATUS_IGNORE, mpi_wait
  use mpi_communication_mod, only : lock_mpi, unlock_mpi, wait_for_mpi_request
  use netcdf_misc_mod, only : check_netcdf_status
  use timeaveraged_time_manipulation_mod, only : serialise_time_averaged_state
  use instantaneous_time_manipulation_mod, only : serialise_instantaneous_state
  implicit none

#ifndef TEST_MODE
  private
#endif

  abstract interface
     subroutine writer_field_manager_serialise_state(byte_data)
       character, dimension(:), allocatable, intent(out) :: byte_data
     end subroutine writer_field_manager_serialise_state

     logical function write_field_manager_determine_if_up_to_date(timestep)
       integer, intent(in) :: timestep
     end function write_field_manager_determine_if_up_to_date
  end interface

  character, dimension(:), allocatable :: serialised_writer_entries, serialised_timeaveraged_manipulation_state, &
       serialised_instantaneous_manipulation_state, serialised_writer_field_manager_state, serialised_writer_entries_time_points
  integer :: global_writer_entry_byte_size(5), global_writer_entry_byte_size_request, my_writer_entry_start_point(5), &
       my_writer_entry_start_request
  procedure(writer_field_manager_serialise_state), pointer :: serialise_writer_field_manager_state
  procedure(write_field_manager_determine_if_up_to_date), pointer :: is_write_field_manager_up_to_date
  
  public define_io_server_state_contributions, write_io_server_state, package_io_server_state, &
       set_serialise_write_field_manager_state, is_io_server_state_writer_ready
contains

  !> Sets the procedure to call for serialises the field manager state, this is handled in this manner due to the ordering
  !! dependencies between these modules
  !! @param serialise_writer_field_manager_state_arg The procedure to call for serialising the field manager state
  subroutine set_serialise_write_field_manager_state(serialise_writer_field_manager_state_arg, &
       is_write_field_manager_up_to_date_arg)
    procedure(writer_field_manager_serialise_state) :: serialise_writer_field_manager_state_arg
    procedure(write_field_manager_determine_if_up_to_date) :: is_write_field_manager_up_to_date_arg
    
    serialise_writer_field_manager_state=>serialise_writer_field_manager_state_arg
    is_write_field_manager_up_to_date=>is_write_field_manager_up_to_date_arg
  end subroutine set_serialise_write_field_manager_state

  !> Determines whether the IO server state writer is ready (i.e. state is at a specific level for the timestep)
  !! @timestep The timestep to check the IO server is up to date with
  !! @returns Whether the IO server is ready or not
  logical function is_io_server_state_writer_ready(timestep)
    integer, intent(in) :: timestep

    is_io_server_state_writer_ready=is_write_field_manager_up_to_date(timestep)
  end function is_io_server_state_writer_ready  

  !> Will initially package the IO server state and kick off non-blocking collective calls to aggregate the sizes on each
  !! IO server process which is needed for the NetCDF dimension sizing and location in this.
  !! @param io_configuration IO server configuration
  !! @param writer_entries The writer entries that we are going to store
  subroutine package_io_server_state(io_configuration, writer_entries, time_points)
    type(io_configuration_type), intent(inout) :: io_configuration
    type(writer_type), volatile, dimension(:), intent(inout) :: writer_entries
    type(hashmap_type), volatile, intent(inout) :: time_points

    integer :: i, current_data_point, ierr
    character, dimension(:), allocatable :: dvt_byte_data, temp

    allocate(serialised_writer_entries(kind(i)))
    current_data_point=1
    current_data_point=pack_scalar_field(serialised_writer_entries, current_data_point, size(writer_entries))
    
    do i=1, size(writer_entries)
      call serialise_writer_type(writer_entries(i), dvt_byte_data)
      allocate(temp(size(serialised_writer_entries) + size(dvt_byte_data) + kind(i)))
      temp(1:size(serialised_writer_entries)) = serialised_writer_entries
      current_data_point=size(serialised_writer_entries)+1
      call move_alloc(from=temp, to=serialised_writer_entries)
      current_data_point=pack_scalar_field(serialised_writer_entries, current_data_point, size(dvt_byte_data))
      serialised_writer_entries(current_data_point:current_data_point+size(dvt_byte_data)-1)=dvt_byte_data
      deallocate(dvt_byte_data)
    end do
    global_writer_entry_byte_size(1)=size(serialised_writer_entries)
    my_writer_entry_start_point(1)=size(serialised_writer_entries)

    call serialise_time_averaged_state(serialised_timeaveraged_manipulation_state)
    global_writer_entry_byte_size(2)=size(serialised_timeaveraged_manipulation_state)
    my_writer_entry_start_point(2)=size(serialised_timeaveraged_manipulation_state)

    call serialise_instantaneous_state(serialised_instantaneous_manipulation_state)
    global_writer_entry_byte_size(3)=size(serialised_instantaneous_manipulation_state)
    my_writer_entry_start_point(3)=size(serialised_instantaneous_manipulation_state)

    call serialise_writer_field_manager_state(serialised_writer_field_manager_state)
    global_writer_entry_byte_size(4)=size(serialised_writer_field_manager_state)
    my_writer_entry_start_point(4)=size(serialised_writer_field_manager_state)

    call serialise_writer_entries_time_points(time_points)
    global_writer_entry_byte_size(5)=size(serialised_writer_entries_time_points)
    my_writer_entry_start_point(5)=size(serialised_writer_entries_time_points)

    call lock_mpi()
    call mpi_iallreduce(MPI_IN_PLACE, global_writer_entry_byte_size, size(global_writer_entry_byte_size), MPI_INT, MPI_SUM, &
         io_configuration%io_communicator, global_writer_entry_byte_size_request, ierr)
    call mpi_iscan(MPI_IN_PLACE, my_writer_entry_start_point, size(my_writer_entry_start_point), MPI_INT, MPI_SUM, &
         io_configuration%io_communicator, my_writer_entry_start_request, ierr)
    call unlock_mpi()
  end subroutine package_io_server_state

  !> Serialises the writer entry time points which are held in a hashmap
  !! @param time_points The input time points which are to be serialised
  subroutine serialise_writer_entries_time_points(time_points)
    type(hashmap_type), volatile, intent(inout) :: time_points

    integer :: current_data_point, key
    type(mapentry_type) :: map_entry
    type(iterator_type) :: iterator
    real(kind=DEFAULT_PRECISION) :: r_value
    
    allocate(serialised_writer_entries_time_points(kind(current_data_point) + &
         ((kind(key) + kind(r_value)) * c_size(time_points))))
    current_data_point=1
    current_data_point=pack_scalar_field(serialised_writer_entries_time_points, current_data_point, c_size(time_points))
    iterator=c_get_iterator(time_points)
    do while (c_has_next(iterator))
      map_entry=c_next_mapentry(iterator)
      key=conv_to_integer(map_entry%key)
      r_value=c_get_real(map_entry)
      current_data_point=pack_scalar_field(serialised_writer_entries_time_points, current_data_point, key)
      current_data_point=pack_scalar_field(serialised_writer_entries_time_points, current_data_point, double_real_value=r_value)
    end do    
  end subroutine serialise_writer_entries_time_points  

  !> Defines the dimensions and variables in a NetCDF file that consitute the IO server current state
  !! @param io_configuration IO server configuration
  !! @param netcdf_file The NetCDF file state
  subroutine define_io_server_state_contributions(io_configuration, netcdf_file)
    type(io_configuration_type), intent(inout) :: io_configuration
    type(netcdf_diagnostics_type), intent(inout) :: netcdf_file
    
    integer :: ncdf_dimid, ncdf_varid
    character, dimension(:), allocatable :: byte

    if (.not. allocated(serialised_writer_entries)) then
      call log_log(LOG_ERROR, "Must package the IO server state before commencing file writer operations")
    end if    

    call lock_mpi()
    call check_netcdf_status(nf90_def_dim(netcdf_file%ncid, "io_configuration_dim", &
         size(io_configuration%text_configuration), ncdf_dimid))
    call c_put_integer(netcdf_file%dimension_to_id, "io_configuration_dim", ncdf_dimid)

    call check_netcdf_status(nf90_def_var(netcdf_file%ncid, "io_configuration", NF90_CHAR, ncdf_dimid, ncdf_varid))
    call c_put_integer(netcdf_file%variable_to_id, "io_configuration", ncdf_varid)

    call check_netcdf_status(nf90_def_dim(netcdf_file%ncid, "entries_directory_dim", io_configuration%number_of_io_servers, &
         ncdf_dimid))
    call unlock_mpi()
    call c_put_integer(netcdf_file%dimension_to_id, "entries_directory_dim", ncdf_dimid)

    call wait_for_mpi_request(global_writer_entry_byte_size_request)

    call define_state_storage(netcdf_file, ncdf_dimid, "serialised_writer_entry", global_writer_entry_byte_size(1))
    call define_state_storage(netcdf_file, ncdf_dimid, "serialised_timeaveraged_manipulation", &
         global_writer_entry_byte_size(2))
    call define_state_storage(netcdf_file, ncdf_dimid, "serialised_instantaneous_manipulation", &
         global_writer_entry_byte_size(3))
    call define_state_storage(netcdf_file, ncdf_dimid, "serialised_writer_manager", global_writer_entry_byte_size(4))
    call define_state_storage(netcdf_file, ncdf_dimid, "serialised_timepoints", global_writer_entry_byte_size(5))
  end subroutine define_io_server_state_contributions

  !> Defines some state storate for a specific facet of the IO server. This creates the directory (location for each IO server
  !! where to load their data), the dimension and the variable
  !! @param netcdf_file The NetCDF file state
  !! @param entries_directory_dim_id NetCDF dimension ID for the directory entries
  !! @param base_key The base key for this storage
  !! @param expected_global_entries The number of expected global entries that will be stored
  subroutine define_state_storage(netcdf_file, entries_directory_dim_id, base_key, expected_global_entries)
    type(netcdf_diagnostics_type), intent(inout) :: netcdf_file
    integer, intent(in) :: entries_directory_dim_id, expected_global_entries
    character(len=*), intent(in) :: base_key

    integer :: ncdf_varid, ncdf_dimid

    call lock_mpi()
    call check_netcdf_status(nf90_def_var(netcdf_file%ncid, trim(base_key)//"_directory", &
         NF90_INT, entries_directory_dim_id, ncdf_varid))
    call c_put_integer(netcdf_file%variable_to_id, trim(base_key)//"_directory", ncdf_varid)

    call check_netcdf_status(nf90_def_dim(netcdf_file%ncid, trim(base_key)//"_dim", expected_global_entries, ncdf_dimid))
    call c_put_integer(netcdf_file%dimension_to_id,  trim(base_key)//"_dim", ncdf_dimid)

    call check_netcdf_status(nf90_def_var(netcdf_file%ncid, trim(base_key), NF90_CHAR, ncdf_dimid, ncdf_varid))
    call c_put_integer(netcdf_file%variable_to_id, trim(base_key), ncdf_varid)
    call unlock_mpi()
  end subroutine define_state_storage  

  !> Writes the actual IO server state into the NetCDF file
  !! @param io_configuration IO server configuration
  !! @param netcdf_file The NetCDF file state
  subroutine write_io_server_state(io_configuration, netcdf_file)
    type(io_configuration_type), intent(inout) :: io_configuration
    type(netcdf_diagnostics_type), intent(inout) :: netcdf_file

    integer :: field_id

    if (.not. allocated(serialised_writer_entries)) then
      call log_log(LOG_ERROR, "Must package the IO server state before commencing file writer operations")
    end if

    field_id=c_get_integer(netcdf_file%variable_to_id, "io_configuration")

    if (io_configuration%my_io_rank == 0) then
      call lock_mpi()
      call check_netcdf_status(nf90_put_var(netcdf_file%ncid, field_id, io_configuration%text_configuration, &
           count=(/size(io_configuration%text_configuration)/)))
      call unlock_mpi()
    end if
    call wait_for_mpi_request(my_writer_entry_start_request)

    call write_state_storage(netcdf_file, my_writer_entry_start_point(1), io_configuration%my_io_rank, &
         "serialised_writer_entry", serialised_writer_entries)
    deallocate(serialised_writer_entries)

    call write_state_storage(netcdf_file, my_writer_entry_start_point(2), io_configuration%my_io_rank, &
         "serialised_timeaveraged_manipulation", serialised_timeaveraged_manipulation_state)
    deallocate(serialised_timeaveraged_manipulation_state)

    call write_state_storage(netcdf_file, my_writer_entry_start_point(3), io_configuration%my_io_rank, &
         "serialised_instantaneous_manipulation", serialised_instantaneous_manipulation_state)
    deallocate(serialised_instantaneous_manipulation_state)

    call write_state_storage(netcdf_file, my_writer_entry_start_point(4), io_configuration%my_io_rank, &
         "serialised_writer_manager", serialised_writer_field_manager_state)
    deallocate(serialised_writer_field_manager_state)

    call write_state_storage(netcdf_file, my_writer_entry_start_point(5), io_configuration%my_io_rank, &
         "serialised_timepoints", serialised_writer_entries_time_points)
    deallocate(serialised_writer_entries_time_points)
  end subroutine write_io_server_state

  !> Writes out the state for a specific facet into the NetCDF file
  !! @param netcdf_file The NetCDF file state
  !! @param writer_entry_start_point The start point, this is uncorrected so it is actually the end point
  !! @param my_io_rank My IO server rank
  !! @param base_key The base key to use for look up
  !! @param raw_byte_code The raw byte code to write
  subroutine write_state_storage(netcdf_file, writer_entry_start_point, my_io_rank, base_key, raw_byte_code)
    type(netcdf_diagnostics_type), intent(inout) :: netcdf_file
    integer, intent(in) :: writer_entry_start_point, my_io_rank
    character(len=*), intent(in) :: base_key
    character, dimension(:), intent(in) :: raw_byte_code
    
    integer :: field_id, corrected_writer_entry_start

    corrected_writer_entry_start=(writer_entry_start_point-size(raw_byte_code))+1
    
    field_id=c_get_integer(netcdf_file%variable_to_id, trim(base_key)//"_directory")
    call lock_mpi()
    call check_netcdf_status(nf90_put_var(netcdf_file%ncid, field_id, corrected_writer_entry_start, (/ my_io_rank+1 /)))
    
    field_id=c_get_integer(netcdf_file%variable_to_id, trim(base_key))
    call check_netcdf_status(nf90_put_var(netcdf_file%ncid, field_id, raw_byte_code, &
         start=(/corrected_writer_entry_start/), count=(/size(raw_byte_code)/)))
    call unlock_mpi()
  end subroutine write_state_storage  
end module io_server_state_writer_mod
