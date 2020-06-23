!> Reads the IO server state that was stored in a NetCDF checkpoint file
module io_server_state_reader_mod
  use iso_c_binding, only: c_int, c_char, c_null_char, c_size_t, c_intptr_t, c_ptr, c_loc, c_sizeof, c_long
  use datadefn_mod, only : DEFAULT_PRECISION, STRING_LENGTH
  use netcdf, only : nf90_global, nf90_nowrite, nf90_inquire_attribute, nf90_open, nf90_inq_dimid, nf90_inquire_dimension, &
       nf90_inq_varid, nf90_get_var, nf90_get_att, nf90_close
  use netcdf_misc_mod, only : check_netcdf_status
  use collections_mod, only : hashmap_type, c_put_real
  use conversions_mod, only : conv_to_string 
  use writer_types_mod, only : writer_type, unserialise_writer_type
  use logging_mod, only : LOG_ERROR, LOG_WARN, log_log, log_master_log
  use mpi, only : mpi_comm_rank, mpi_comm_size
  use data_utils_mod, only : unpack_scalar_integer_from_bytedata, unpack_scalar_dp_real_from_bytedata
  use timeaveraged_time_manipulation_mod, only : unserialise_time_averaged_state
  use instantaneous_time_manipulation_mod, only : unserialise_instantaneous_state
  use configuration_parser_mod, only : io_configuration_type
  use optionsdatabase_mod, only : options_get_string
  implicit none

#ifndef TEST_MODE
  private
#endif

  interface
     !> ISO C binding for NetCDF inquire dimension, required for 64 bit dimension length
     function nc_inq_dim(ncid, dimid, name, lenp) bind(C)
       use iso_c_binding, only: c_int, c_size_t, c_char

       integer(kind=c_int), value :: ncid
       integer(kind=c_int), value :: dimid
       character(kind=c_char), intent(inout) :: name(*)
       integer(kind=c_size_t), intent(out) :: lenp
       integer(kind=c_int) :: nc_inq_dim
     end function nc_inq_dim

     !> ISO C binding for NetCDF get long scalar variable, required for retrieving long variables
     function nc_get_vara_long(ncid, varid, startp, countp, ip) bind(C)
       use iso_c_binding, only: c_int, c_long, c_ptr

       integer(kind=c_int), value :: ncid, varid
       type(c_ptr), value :: startp, countp
       integer(kind=c_long), intent(out) :: ip(*)
       integer(kind=c_int) :: nc_get_vara_long
     end function nc_get_vara_long

     !> ISO C binding for NetCDF get text vars, required for 64 bit start, count & stride
     function nc_get_vars_text(ncid, varid, startp, countp, stridep, ip) bind(C)
       use iso_c_binding, only: c_int, c_ptr, c_char

       integer(kind=c_int), value :: ncid, varid
       type(c_ptr), value :: startp, countp, stridep
       character(kind=c_char), intent(out) :: ip(*)
       integer(kind=c_int) :: nc_get_vars_text
     end function nc_get_vars_text
  end interface

  abstract interface
     subroutine writer_field_manager_unserialise_state(byte_data)
       character, dimension(:), intent(in) :: byte_data
     end subroutine writer_field_manager_unserialise_state
  end interface

  public read_io_server_configuration, reactivate_writer_federator_state, reactivate_writer_field_manager_state
contains

  !> Reads the IO server configuration, which is the XML configuration initially run with and stored in the checkpoint. Note
  !! that this will open, read the XML in and then close the file.
  !! @param checkpoint_filename The checkpoint filename to open and read from
  !! @param io_xml_configuration XML configuration is read from the checkpoint and placed into here
  !! @param io_communicator_arg The MPI IO server communicator
  subroutine read_io_server_configuration(checkpoint_filename, io_xml_configuration, io_communicator_arg)
    character(len=STRING_LENGTH), intent(in) :: checkpoint_filename
    character, dimension(:), allocatable, intent(inout) :: io_xml_configuration
    integer, intent(in) :: io_communicator_arg

    integer :: ncid, number_io_server, my_io_server_rank, ierr
    integer :: dim_id, dim_size
    logical :: found

    call mpi_comm_rank(io_communicator_arg, my_io_server_rank, ierr)
    call mpi_comm_size(io_communicator_arg, number_io_server, ierr)
    call check_netcdf_status(nf90_open(path = checkpoint_filename, mode = nf90_nowrite, ncid = ncid))
    call check_netcdf_status(nf90_inq_dimid(ncid, "entries_directory_dim", dim_id), found)
    if (.not. found) then
      if (my_io_server_rank==0) then
        call log_log(LOG_WARN, "Restarting the IO server fresh as the checkpoint file does not contain IO state")
      end if      
      return
    end if    
    call check_netcdf_status(nf90_inquire_dimension(ncid, dim_id, len=dim_size))    

    if (dim_size .ne. number_io_server) then
      call log_log(LOG_ERROR, "Can not restart IO server with a different number of IO servers")
    end if
    call get_io_server_configuration(ncid, io_xml_configuration)
    call check_netcdf_status(nf90_close(ncid))
  end subroutine read_io_server_configuration

  !> Reactivates the writer federator and everything beneath it (i.e. just not the writer field manager.) For memory reasons
  !! this explicitly reopens the checkpoint file, will read each individual byte code entry in & repackage before deallocating
  !! memory and moving onto the next entry. The file is then closed
  !! @param io_configuration The IO server configuration
  !! @param writer_entries The configured but empty writer entries to unpack the state into
  !! @param time_points Time points to unpack the state into
  subroutine reactivate_writer_federator_state(io_configuration, writer_entries, time_points)
    type(io_configuration_type), intent(inout) :: io_configuration
    type(writer_type), dimension(:) :: writer_entries
    type(hashmap_type), volatile, intent(inout) :: time_points

    integer :: ncid, ierr, i
    character, dimension(:), allocatable :: raw_bytes

    call check_netcdf_status(nf90_open(path = options_get_string(io_configuration%options_database, "checkpoint"), &
         mode = nf90_nowrite, ncid = ncid))

    do i=1, size(writer_entries)
      ! Note that the different writer entries are dealt with separately for memory reasons
      if (writer_entries(i)%include_in_io_state_write) then
        call get_io_server_serialised_bytes(ncid, io_configuration%number_of_io_servers, io_configuration%my_io_rank, &
             "serialised_writer_entry_"//trim(conv_to_string(i)), raw_bytes)
        call unserialise_writer_type(writer_entries(i), raw_bytes)
        deallocate(raw_bytes)
      end if
    end do

    call get_io_server_serialised_bytes(ncid, io_configuration%number_of_io_servers, io_configuration%my_io_rank, &
         "serialised_timeaveraged_manipulation", raw_bytes)
    call restart_timeaveraged_state_from_checkpoint(raw_bytes)
    deallocate(raw_bytes)

    call get_io_server_serialised_bytes(ncid, io_configuration%number_of_io_servers, io_configuration%my_io_rank, &
         "serialised_instantaneous_manipulation", raw_bytes)
    call restart_instantaneous_state_from_checkpoint(raw_bytes)
    deallocate(raw_bytes)

    call get_io_server_serialised_bytes(ncid, io_configuration%number_of_io_servers, io_configuration%my_io_rank, &
         "serialised_timepoints", raw_bytes)
    call restart_writer_state_timepoints(time_points, raw_bytes)
    deallocate(raw_bytes)
    
    call check_netcdf_status(nf90_close(ncid))
  end subroutine reactivate_writer_federator_state

  !> Reactivates the writer field manager state from the checkpoint file, for memory reasons this will open the file,
  !! read in and deserialise the byte code before closing it
  !! @param io_configuration IO server configuration
  !! @param unserialise_writer_field_manager Procedure pointer to the unserialisation for the writer field manager
  subroutine reactivate_writer_field_manager_state(io_configuration, unserialise_writer_field_manager)
    type(io_configuration_type), intent(inout) :: io_configuration
    procedure(writer_field_manager_unserialise_state) :: unserialise_writer_field_manager

    integer :: ncid, ierr
    character, dimension(:), allocatable :: raw_bytes

    call check_netcdf_status(nf90_open(path = options_get_string(io_configuration%options_database, "checkpoint"), &
         mode = nf90_nowrite, ncid = ncid))

    call get_io_server_serialised_bytes(ncid, io_configuration%number_of_io_servers, io_configuration%my_io_rank, &
         "serialised_writer_manager", raw_bytes)
    call restart_writer_field_manager_from_checkpoint(unserialise_writer_field_manager, raw_bytes)
    deallocate(raw_bytes)

    call check_netcdf_status(nf90_close(ncid))
  end subroutine reactivate_writer_field_manager_state  

  !> Retrieves the IO server XML configuration from the checkpoint file
  !! @param ncid The NetCDF checkpoint file ID
  !! @param io_xml_configuration XML configuration is read from the checkpoint and placed into here
  subroutine get_io_server_configuration(ncid, io_xml_configuration)
    integer, intent(in) :: ncid
    character, dimension(:), allocatable, intent(inout) :: io_xml_configuration

    integer :: dim_id, var_id, dim_size
    logical :: found

    call check_netcdf_status(nf90_inq_dimid(ncid, "io_configuration_dim", dim_id), found)
    if (.not. found) return
    call check_netcdf_status(nf90_inquire_dimension(ncid, dim_id, len=dim_size))

    call check_netcdf_status(nf90_inq_varid(ncid, "io_configuration", var_id), found)
    if (.not. found) return
    allocate(io_xml_configuration(dim_size))
    call check_netcdf_status(nf90_get_var(ncid, var_id, io_xml_configuration, count=(/dim_size/)))
  end subroutine get_io_server_configuration

  !> Retrieves some IO server serialised bytes which will make up the state of a specific facet. Note that this uses the
  !! ISO C bindings in order to support 64 bit starts, counts & strides along with 64 bit scalar long variable fields.
  !! @param ncid The NetCDF checkpoint file ID
  !! @param number_io_server The number of IO servers
  !! @param my_io_server_rank My IO server rank
  !! @param base_key The base look up key to retrieve the state which is stored in the file
  !! @param raw_bytes The retrieved raw bytes as they sit in the file, this is allocated to hold the required data
  subroutine get_io_server_serialised_bytes(ncid, number_io_server, my_io_server_rank, base_key, raw_bytes)
    integer, intent(in) :: ncid, number_io_server, my_io_server_rank
    character(len=*), intent(in) :: base_key
    character, dimension(:), allocatable, intent(out) :: raw_bytes
    
    integer :: dim_id, var_id
    logical :: found
    integer(kind=8) :: dim_size, serialised_range(2), number_serialised_entries

    integer(kind=c_int)        :: cncid, cdimid, cstatus, cvarid
    integer(kind=c_size_t)     :: cdlen
    character(len=256) :: tmpname
    integer(KIND=c_size_t), target :: cstart(1), ccounts(1)
    Integer(KIND=c_intptr_t), target :: cstrides(1)
    type(c_ptr) :: cstartptr, ccountsptr, cstridesptr

    cncid=ncid
    cstartptr=c_loc(cstart)
    ccountsptr=c_loc(ccounts)
    cstridesptr=c_loc(cstrides)

    call check_netcdf_status(nf90_inq_dimid(ncid, trim(base_key)//"_dim", dim_id), found)
    if (.not. found) return
    
    cdimid=dim_id-1
    call check_netcdf_status(nc_inq_dim(cncid, cdimid, tmpname, cdlen))
    dim_size=cdlen

    call check_netcdf_status(nf90_inq_varid(ncid, trim(base_key)//"_directory", var_id), found)
    if (.not. found) return
    cvarid=var_id-1
    cstart(1)=my_io_server_rank
    if (my_io_server_rank .lt. number_io_server-1) then
      ccounts(1)=2
      call check_netcdf_status(nc_get_vara_long(cncid, cvarid, cstartptr, ccountsptr, serialised_range))
      if (serialised_range(2) .gt. dim_size) then
        call log_log(LOG_ERROR, "Serialised entry beyond size in the file")
      end if
    else
      ccounts(1)=1
      call check_netcdf_status(nc_get_vara_long(cncid, cvarid, cstartptr, ccountsptr, serialised_range))
      serialised_range(2)=dim_size
    end if    
    number_serialised_entries=(serialised_range(2)-serialised_range(1)) + 1    
    call check_netcdf_status(nf90_inq_varid(ncid, trim(base_key), var_id), found)
    if (.not. found) return
    allocate(raw_bytes(number_serialised_entries))

    cvarid=var_id-1
    cstart=serialised_range(1)-1
    ccounts=number_serialised_entries
    cstrides(1)=1
    call check_netcdf_status(nc_get_vars_text(cncid, cvarid, cstartptr, ccountsptr, cstridesptr, raw_bytes))
  end subroutine get_io_server_serialised_bytes

  !> Restarts the writer state from a specific checkpoint byte data chunk of memory
  !! @param writer_entries The array of writer entries to fill in
  !! @param raw_byte The serialised byte state to unpackage and restart from
  subroutine restart_writer_state_from_checkpoint(writer_entries, raw_bytes)
    type(writer_type), dimension(:) :: writer_entries
    character, dimension(:), allocatable :: raw_bytes

    integer :: i, number_entries, current_point, byte_size
    
    if (.not. allocated(raw_bytes)) then
      call log_master_log(LOG_WARN, "On restart no writer state in checkpoint file")
      return
    end if
    current_point=1
    number_entries=unpack_scalar_integer_from_bytedata(raw_bytes, current_point)
    if (number_entries .ne. size(writer_entries)) then
      call log_log(LOG_ERROR, "On restart have a different number of configured entries than those in the checkpoint file")
    end if
    do i=1, size(writer_entries)
      if (writer_entries(i)%include_in_io_state_write) then
        byte_size=unpack_scalar_integer_from_bytedata(raw_bytes, current_point)
        call unserialise_writer_type(writer_entries(i), raw_bytes(current_point:current_point+byte_size-1))
        current_point=current_point+byte_size
      end if
    end do
  end subroutine restart_writer_state_from_checkpoint

  !> Restarts the writer state timepoints held in the writer federator
  !! @param time_points The timepoints hashmap which is filled in from the serialised version
  !! @param raw_byte The serialised byte state to unpackage and restart from
  subroutine restart_writer_state_timepoints(time_points, raw_bytes)
    type(hashmap_type), volatile, intent(inout) :: time_points
    character, dimension(:), allocatable :: raw_bytes

    integer :: i, number_entries, current_point, byte_size, timestep_key
    real(kind=DEFAULT_PRECISION) :: r_value

    if (.not. allocated(raw_bytes)) then
      call log_master_log(LOG_WARN, "On restart no writer state timepoints in checkpoint file")
      return
    end if
    current_point=1
    number_entries=unpack_scalar_integer_from_bytedata(raw_bytes, current_point)
    do i=1, number_entries
      timestep_key=unpack_scalar_integer_from_bytedata(raw_bytes, current_point)
      r_value=unpack_scalar_dp_real_from_bytedata(raw_bytes, current_point)
      call c_put_real(time_points, conv_to_string(timestep_key), r_value)
    end do
  end subroutine restart_writer_state_timepoints  

  !> Will restart the time averaged manipulation state from the checkpoint file
  !! @param raw_byte The serialised byte state to unpackage and restart from
  subroutine restart_timeaveraged_state_from_checkpoint(raw_bytes)
    character, dimension(:), allocatable :: raw_bytes
    
    if (.not. allocated(raw_bytes)) then
      call log_master_log(LOG_WARN, "On restart no time averaged state in checkpoint file")
      return
    end if
    call unserialise_time_averaged_state(raw_bytes)
  end subroutine restart_timeaveraged_state_from_checkpoint

  !> Will restart the instantaneous manipulation state from the checkpoint file
  !! @param raw_byte The serialised byte state to unpackage and restart from
  subroutine restart_instantaneous_state_from_checkpoint(raw_bytes)
    character, dimension(:), allocatable :: raw_bytes
    
    if (.not. allocated(raw_bytes)) then
      call log_master_log(LOG_WARN, "On restart no instantaneous state in checkpoint file")
      return
    end if
    call unserialise_instantaneous_state(raw_bytes)
  end subroutine restart_instantaneous_state_from_checkpoint

  !> Will restart the field manager state from the checkpoint file
  !! @param unserialise_writer_field_manager The unserialise field manager procedure, done this way due to module ordering
  !! @param raw_byte The serialised byte state to unpackage and restart from
  subroutine restart_writer_field_manager_from_checkpoint(unserialise_writer_field_manager, raw_bytes)
    procedure(writer_field_manager_unserialise_state) :: unserialise_writer_field_manager
    character, dimension(:), allocatable :: raw_bytes

    if (.not. allocated(raw_bytes)) then
      call log_master_log(LOG_WARN, "On restart no writer field manager state in checkpoint file")
      return
    end if
    call unserialise_writer_field_manager(raw_bytes)
  end subroutine restart_writer_field_manager_from_checkpoint  
end module io_server_state_reader_mod
