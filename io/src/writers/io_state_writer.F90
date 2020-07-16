!> The IO server state module which will write out the current state of the IO server to a NetCDF file
module io_server_state_writer_mod
  use iso_c_binding, only: c_int, c_char, c_null_char, c_size_t, c_intptr_t, c_ptr, c_loc, c_sizeof, c_long
  use datadefn_mod, only : DEFAULT_PRECISION
  use netcdf, only : NF90_CHAR, NF90_BYTE, NF90_INT, nf90_def_var, nf90_put_var, nf90_def_dim
  use configuration_parser_mod, only : io_configuration_type
  use writer_types_mod, only : netcdf_diagnostics_type, serialise_writer_type, prepare_to_serialise_writer_type, writer_type
  use collections_mod, only : hashmap_type, mapentry_type, iterator_type, c_put_integer, c_get_integer, c_get_iterator, &
       c_next_mapentry, c_has_next, c_get_real, c_size
  use conversions_mod, only : conv_to_integer, conv_to_string
  use logging_mod, only : LOG_ERROR, LOG_WARN, log_log, log_master_log
  use io_server_client_mod, only : pack_scalar_field
  use mpi, only : MPI_IN_PLACE, MPI_INT, MPI_LONG_LONG_INT, MPI_SUM, MPI_STATUS_IGNORE, mpi_wait
  use mpi_communication_mod, only : lock_mpi, unlock_mpi, wait_for_mpi_request
  use netcdf_misc_mod, only : check_netcdf_status
  use timeaveraged_time_manipulation_mod, only : prepare_to_serialise_time_averaged_state, serialise_time_averaged_state
  use instantaneous_time_manipulation_mod, only : prepare_to_serialise_instantaneous_state, serialise_instantaneous_state
  implicit none

#ifndef TEST_MODE
  private
#endif

  interface
     !> ISO C binding interface for NetCDF dimension definition, needed to support 64 bit lengths
     function nc_def_dim(ncid, name, nlen, idp) bind(C)
       use iso_c_binding, only: c_int, c_size_t, c_char
       
       integer(kind=c_int), value :: ncid
       character(kind=c_char), intent(in) :: name(*)
       integer(kind=c_size_t), value :: nlen
       integer(kind=c_int), intent(inout) :: idp
       integer(kind=c_int) :: nc_def_dim
     end function nc_def_dim

     !> ISO C binding interface for NetCDF put variable text, needed to support 64 bit starts, counts and strides
     function nc_put_vars_text(ncid, varid, startp, countp, stridep, op)  bind(C)
       use iso_c_binding, only: c_int, c_ptr, c_char

       integer(kind=c_int), value :: ncid, varid
       type(c_ptr), value :: startp, countp, stridep
       character(kind=c_char), intent(in) :: op(*)
       integer(kind=c_int) :: nc_put_vars_text
     end function nc_put_vars_text

     !> ISO C binding interface for NetCDF put long scalar, needed to support putting longs into file
     function nc_put_var1_long(ncid, varid, indexp, op) bind(C)
       use iso_c_binding, only: c_int, c_ptr, c_long

       integer(kind=c_int), value :: ncid, varid
       type(c_ptr), value :: indexp
       integer(kind=c_long), intent(in) :: op
       integer(kind=c_int) :: nc_put_var1_long
     end function nc_put_var1_long

     !> ISO C binding interface for NetCDF define variable, needed to support defining a long scalar variable
     function nc_def_var(ncid, name, xtype, ndims, dimidsp, varidp) bind(C)
       use iso_c_binding, only: c_int, c_char

       integer(kind=c_int), value :: ncid
       character(kind=c_char), intent(in) :: name(*)
       integer(kind=c_int), value :: xtype
       integer(kind=c_int), value :: ndims
       integer(kind=c_int), intent(in) :: dimidsp(*)
       integer(kind=c_int), intent(out) :: varidp
       integer(kind=c_int) :: nc_def_var
     end function nc_def_var
  end interface

  abstract interface
     integer(kind=8) function writer_field_manager_prepare_to_serialise_state()
     end function writer_field_manager_prepare_to_serialise_state     
     
     subroutine writer_field_manager_serialise_state(byte_data)
       character, dimension(:), allocatable, intent(inout) :: byte_data
     end subroutine writer_field_manager_serialise_state

     logical function write_field_manager_determine_if_up_to_date(timestep)
       integer, intent(in) :: timestep
     end function write_field_manager_determine_if_up_to_date
  end interface

  character, dimension(:), allocatable :: serialised_writer_entries, serialised_timeaveraged_manipulation_state, &
       serialised_instantaneous_manipulation_state, serialised_writer_field_manager_state, serialised_writer_entries_time_points
  integer(kind=8), dimension(:), allocatable :: global_writer_entry_byte_size, my_writer_entry_start_point, &
       local_writer_entry_byte_size
  integer :: global_writer_entry_byte_size_request, my_writer_entry_start_request
  procedure(writer_field_manager_serialise_state), pointer :: serialise_writer_field_manager_state
  procedure(writer_field_manager_prepare_to_serialise_state), pointer :: prepare_to_serialise_writer_field_manager_state
  procedure(write_field_manager_determine_if_up_to_date), pointer :: is_write_field_manager_up_to_date
  
  public define_io_server_state_contributions, write_io_server_state, set_serialise_write_field_manager_state, &
       is_io_server_state_writer_ready
contains

  !> Sets the procedure to call for serialises the field manager state, this is handled in this manner due to the ordering
  !! dependencies between these modules
  !! @param serialise_writer_field_manager_state_arg The procedure to call for serialising the field manager state
  !! @param prepare_to_serialise_writer_field_manager_state_arg Preparation of field manager state procedure
  !! @param is_write_field_manager_up_to_date_arg Procedure to determine whether field manager is up to date
  subroutine set_serialise_write_field_manager_state(serialise_writer_field_manager_state_arg, &
       prepare_to_serialise_writer_field_manager_state_arg, is_write_field_manager_up_to_date_arg)    
    procedure(writer_field_manager_serialise_state) :: serialise_writer_field_manager_state_arg
    procedure(writer_field_manager_prepare_to_serialise_state) :: prepare_to_serialise_writer_field_manager_state_arg
    procedure(write_field_manager_determine_if_up_to_date) :: is_write_field_manager_up_to_date_arg

    integer(kind=c_size_t) :: size_t_test
    
    serialise_writer_field_manager_state=>serialise_writer_field_manager_state_arg
    prepare_to_serialise_writer_field_manager_state=>prepare_to_serialise_writer_field_manager_state_arg
    is_write_field_manager_up_to_date=>is_write_field_manager_up_to_date_arg

    if (c_sizeof(size_t_test) .lt. 8) then
      call log_master_log(LOG_WARN, &
           "Your system's size_t is not 64 bit, this will limit the size of IO server state storage to 4GB")
    end if    
  end subroutine set_serialise_write_field_manager_state

  !> Determines whether the IO server state writer is ready (i.e. state is at a specific level for the timestep)
  !! @timestep The timestep to check the IO server is up to date with
  !! @returns Whether the IO server is ready or not
  logical function is_io_server_state_writer_ready(timestep)
    integer, intent(in) :: timestep

    is_io_server_state_writer_ready=is_write_field_manager_up_to_date(timestep)
  end function is_io_server_state_writer_ready  

  !> Will determine the size of the package for different facets of the IO server state and kick off non-blocking collective
  !! calls to aggregate the sizes on each IO server process which is needed for the NetCDF dimension sizing and location in this.
  !! Note that during this certain locks are issued to ensure sizes don't change between this and physical packaging
  !! @param io_configuration IO server configuration
  !! @param writer_entries The writer entries that we are going to store
  subroutine prepare_io_server_state_storage(io_configuration, writer_entries, time_points)
    type(io_configuration_type), intent(inout) :: io_configuration
    type(writer_type), volatile, dimension(:), intent(inout) :: writer_entries
    type(hashmap_type), volatile, intent(inout) :: time_points

    integer :: i, current_data_point, number_writer_entries_included, ierr, current_index
    character, dimension(:), allocatable :: dvt_byte_data, temp

    number_writer_entries_included=0
    do i=1, size(writer_entries)
      if (writer_entries(i)%include_in_io_state_write) number_writer_entries_included=number_writer_entries_included+1
    end do

    allocate(local_writer_entry_byte_size(number_writer_entries_included+4), &
         global_writer_entry_byte_size(number_writer_entries_included+4), &
         my_writer_entry_start_point(number_writer_entries_included+4))

    current_index=1
    do i=1, size(writer_entries)
      if (writer_entries(i)%include_in_io_state_write) then
        local_writer_entry_byte_size(current_index)=prepare_to_serialise_writer_type(writer_entries(i))
        current_index=current_index+1
      end if
    end do

    local_writer_entry_byte_size(current_index)=prepare_to_serialise_time_averaged_state()
    local_writer_entry_byte_size(current_index+1)=prepare_to_serialise_instantaneous_state()
    local_writer_entry_byte_size(current_index+2)=prepare_to_serialise_writer_field_manager_state()
    local_writer_entry_byte_size(current_index+3)=prepare_to_serialise_writer_entries_time_points(time_points)
    
    global_writer_entry_byte_size=local_writer_entry_byte_size
    my_writer_entry_start_point=local_writer_entry_byte_size
    call lock_mpi()
    call mpi_iallreduce(MPI_IN_PLACE, global_writer_entry_byte_size, size(global_writer_entry_byte_size), MPI_LONG_LONG_INT, &
         MPI_SUM, io_configuration%io_communicator, global_writer_entry_byte_size_request, ierr)
    call mpi_iscan(MPI_IN_PLACE, my_writer_entry_start_point, size(my_writer_entry_start_point), MPI_LONG_LONG_INT, &
         MPI_SUM, io_configuration%io_communicator, my_writer_entry_start_request, ierr)
    call unlock_mpi()
  end subroutine prepare_io_server_state_storage

  !> Prepares to serialise the writer entry time points
  !! @param time_points Hashmap of the timepoints that will need to be serialised
  !! @returns The byte size required for this serialisation
  integer(kind=8) function prepare_to_serialise_writer_entries_time_points(time_points)
    type(hashmap_type), volatile, intent(inout) :: time_points

    real(kind=DEFAULT_PRECISION) :: r_value

    prepare_to_serialise_writer_entries_time_points=kind(prepare_to_serialise_writer_entries_time_points) + &
         ((kind(prepare_to_serialise_writer_entries_time_points) + kind(r_value)) * c_size(time_points))
  end function prepare_to_serialise_writer_entries_time_points

  !> Serialises the writer entry time points which are held in a hashmap
  !! @param time_points The input time points which are to be serialised
  !! @param byte_data Allocated byte data which will hold the serialised values
  subroutine serialise_writer_entries_time_points(time_points, byte_data)
    type(hashmap_type), volatile, intent(inout) :: time_points
    character, dimension(:), allocatable, intent(inout) :: byte_data

    integer :: key, current_data_point
    type(mapentry_type) :: map_entry
    type(iterator_type) :: iterator
    real(kind=DEFAULT_PRECISION) :: r_value
        
    current_data_point=1
    current_data_point=pack_scalar_field(byte_data, current_data_point, c_size(time_points))
    iterator=c_get_iterator(time_points)
    do while (c_has_next(iterator))
      map_entry=c_next_mapentry(iterator)
      key=conv_to_integer(map_entry%key)
      r_value=c_get_real(map_entry)
      current_data_point=pack_scalar_field(byte_data, current_data_point, key)
      current_data_point=pack_scalar_field(byte_data, current_data_point, double_real_value=r_value)
    end do    
  end subroutine serialise_writer_entries_time_points  

  !> Defines the dimensions and variables in a NetCDF file that consitute the IO server current state. This will call out
  !! to prepare all IO state for storage (determines the size of each byte code and issues locks for consistency.)
  !! @param io_configuration IO server configuration
  !! @param netcdf_file The NetCDF file state
  subroutine define_io_server_state_contributions(io_configuration, writer_entries, time_points, netcdf_file)
    type(io_configuration_type), intent(inout) :: io_configuration
    type(writer_type), volatile, dimension(:), intent(inout) :: writer_entries
    type(hashmap_type), volatile, intent(inout) :: time_points
    type(netcdf_diagnostics_type), intent(inout) :: netcdf_file
    
    integer :: ncdf_dimid, ncdf_varid, current_index, i
    character, dimension(:), allocatable :: byte

    call prepare_io_server_state_storage(io_configuration, writer_entries, time_points)
    
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

    current_index=1
    do i=1, size(writer_entries)
      if (writer_entries(i)%include_in_io_state_write) then
        call define_state_storage(netcdf_file, ncdf_dimid, "serialised_writer_entry_"//trim(conv_to_string(i)), &
             global_writer_entry_byte_size(current_index))
        current_index=current_index+1
      end if
    end do        
    
    call define_state_storage(netcdf_file, ncdf_dimid, "serialised_timeaveraged_manipulation", &
         global_writer_entry_byte_size(current_index))
    call define_state_storage(netcdf_file, ncdf_dimid, "serialised_instantaneous_manipulation", &
         global_writer_entry_byte_size(current_index+1))
    call define_state_storage(netcdf_file, ncdf_dimid, "serialised_writer_manager", global_writer_entry_byte_size(current_index+2))
    call define_state_storage(netcdf_file, ncdf_dimid, "serialised_timepoints", global_writer_entry_byte_size(current_index+3))
  end subroutine define_io_server_state_contributions

  !> Defines some state storate for a specific facet of the IO server. This creates the directory (location for each IO server
  !! where to load their data), the dimension and the variable. Note this uses the ISO C bindings to NetCDF to support
  !! 64 bit dimension and field lengths
  !! @param netcdf_file The NetCDF file state
  !! @param entries_directory_dim_id NetCDF dimension ID for the directory entries
  !! @param base_key The base key for this storage
  !! @param expected_global_entries The number of expected global entries that will be stored
  subroutine define_state_storage(netcdf_file, entries_directory_dim_id, base_key, expected_global_entries)
    type(netcdf_diagnostics_type), intent(inout) :: netcdf_file
    integer, intent(in) :: entries_directory_dim_id
    integer(kind=8), intent(in) :: expected_global_entries
    character(len=*), intent(in) :: base_key

    integer :: ncdf_varid
    integer :: ncdf_dimid
    integer(kind=c_int) :: cncid, cdimid, cvarid, cxtype, cnvdims
    integer(kind=c_int) :: cvdims(1)

    cncid = netcdf_file%ncid
    cnvdims=1
    cvdims(1)=entries_directory_dim_id-1
    
    cxtype=10 ! NC_INT64 type, note that this is signed as Fortran doesn't really like unsigned integers
    call lock_mpi()
    call check_netcdf_status(nc_def_var(cncid, trim(base_key)//"_directory"//C_NULL_CHAR, cxtype, &
         cnvdims, cvdims, cvarid))
    ncdf_varid=cvarid+1

    call c_put_integer(netcdf_file%variable_to_id, trim(base_key)//"_directory", ncdf_varid)
  
    call check_netcdf_status(nc_def_dim(cncid, trim(base_key)//"_dim"//C_NULL_CHAR, &
         int(expected_global_entries, c_size_t), cdimid))
    ncdf_dimid=cdimid+1
    call c_put_integer(netcdf_file%dimension_to_id,  trim(base_key)//"_dim", ncdf_dimid)
    call check_netcdf_status(nf90_def_var(netcdf_file%ncid, trim(base_key), NF90_CHAR, ncdf_dimid, ncdf_varid))
    call c_put_integer(netcdf_file%variable_to_id, trim(base_key), ncdf_varid)
    call unlock_mpi()
  end subroutine define_state_storage

  !> Packags up and writes the actual IO server state into the NetCDF file. The act of serialisation will effectively unlock
  !! the IO server state locks that were issued in the preparation call (this is to ensure consistency between reported size
  !! and actual size when it comes time to package up.)
  !! @param io_configuration IO server configuration
  !! @param writer_entries The writer types which need to be serialised
  !! @param time_points The time points that need to be serialised
  !! @param netcdf_file The NetCDF file state
  subroutine write_io_server_state(io_configuration, writer_entries, time_points, netcdf_file)
    type(io_configuration_type), intent(inout) :: io_configuration
    type(writer_type), volatile, dimension(:), intent(inout) :: writer_entries
    type(hashmap_type), volatile, intent(inout) :: time_points
    type(netcdf_diagnostics_type), intent(inout) :: netcdf_file

    integer :: field_id, var_id, current_index, i
    character, dimension(:), allocatable :: serialised_bytes

    field_id=c_get_integer(netcdf_file%variable_to_id, "io_configuration")

    if (io_configuration%my_io_rank == 0) then
      call lock_mpi()
      call check_netcdf_status(nf90_put_var(netcdf_file%ncid, field_id, io_configuration%text_configuration, &
           count=(/size(io_configuration%text_configuration)/)))
      call unlock_mpi()
    end if
    call wait_for_mpi_request(my_writer_entry_start_request)

    current_index=1
    do i=1, size(writer_entries)
      if (writer_entries(i)%include_in_io_state_write) then
        allocate(serialised_bytes(local_writer_entry_byte_size(current_index)))
        call serialise_writer_type(writer_entries(i), serialised_bytes)
        call write_state_storage(netcdf_file, my_writer_entry_start_point(current_index), io_configuration%my_io_rank, &
             "serialised_writer_entry_"//trim(conv_to_string(i)), serialised_bytes)
        deallocate(serialised_bytes)
        current_index=current_index+1
      end if
    end do

    allocate(serialised_bytes(local_writer_entry_byte_size(current_index)))
    call serialise_time_averaged_state(serialised_bytes)
    call write_state_storage(netcdf_file, my_writer_entry_start_point(current_index), io_configuration%my_io_rank, &
         "serialised_timeaveraged_manipulation", serialised_bytes)
    deallocate(serialised_bytes)

    allocate(serialised_bytes(local_writer_entry_byte_size(current_index+1)))
    call serialise_instantaneous_state(serialised_bytes)
    call write_state_storage(netcdf_file, my_writer_entry_start_point(current_index+1), io_configuration%my_io_rank, &
         "serialised_instantaneous_manipulation", serialised_bytes)
    deallocate(serialised_bytes)

    allocate(serialised_bytes(local_writer_entry_byte_size(current_index+2)))
    call serialise_writer_field_manager_state(serialised_bytes)
    call write_state_storage(netcdf_file, my_writer_entry_start_point(current_index+2), io_configuration%my_io_rank, &
         "serialised_writer_manager", serialised_bytes)
    deallocate(serialised_bytes)

    allocate(serialised_bytes(local_writer_entry_byte_size(current_index+3)))
    call serialise_writer_entries_time_points(time_points, serialised_bytes)
    call write_state_storage(netcdf_file, my_writer_entry_start_point(current_index+3), io_configuration%my_io_rank, &
         "serialised_timepoints", serialised_bytes)
    deallocate(serialised_bytes)

    deallocate(local_writer_entry_byte_size, global_writer_entry_byte_size, my_writer_entry_start_point)
  end subroutine write_io_server_state

  !> Writes out the state for a specific facet into the NetCDF file. Note that this uses the ISO C bindings into NetCDF
  !! to support 64 bit counts, starts and strides
  !! @param netcdf_file The NetCDF file state
  !! @param writer_entry_start_point The start point, this is uncorrected so it is actually the end point
  !! @param my_io_rank My IO server rank
  !! @param base_key The base key to use for look up
  !! @param raw_byte_code The raw byte code to write
  subroutine write_state_storage(netcdf_file, writer_entry_start_point, my_io_rank, base_key, raw_byte_code)
    type(netcdf_diagnostics_type), intent(inout) :: netcdf_file
    integer, intent(in) :: my_io_rank
    integer(kind=8), intent(in) :: writer_entry_start_point 
    character(len=*), intent(in) :: base_key
    character, dimension(:), intent(in) :: raw_byte_code

    integer(kind=c_int) :: cncid, cvarid, cndims, cstat1, cstatus
    integer(kind=c_size_t), target :: cstart(1), ccounts(1)
    integer(kind=c_intptr_t), target :: cstrides(1)
    integer(kind=c_long) :: c_writer_corrected_start_point
    type(c_ptr) :: cstartptr, ccountsptr, cstridesptr

    cncid=netcdf_file%ncid
    cstart(1)=my_io_rank
    cvarid=c_get_integer(netcdf_file%variable_to_id, trim(base_key)//"_directory")-1
    c_writer_corrected_start_point=(writer_entry_start_point-size(raw_byte_code))+1

    ccounts(1)=size(raw_byte_code)
    cstrides(1)=1

    cstartptr=c_loc(cstart)
    ccountsptr=c_loc(ccounts)
    cstridesptr=c_loc(cstrides)
    
    call lock_mpi()
    call check_netcdf_status(nc_put_var1_long(cncid, cvarid, cstartptr, c_writer_corrected_start_point))
    cvarid=c_get_integer(netcdf_file%variable_to_id, trim(base_key))-1
    cstart(1)=c_writer_corrected_start_point-1    
    call check_netcdf_status(nc_put_vars_text(cncid, cvarid, cstartptr, ccountsptr, cstridesptr, raw_byte_code))
    call unlock_mpi()
  end subroutine write_state_storage
end module io_server_state_writer_mod
