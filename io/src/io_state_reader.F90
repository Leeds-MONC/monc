!> Reads the IO server state that was stored in a NetCDF checkpoint file
module io_server_state_reader_mod
  use datadefn_mod, only : DEFAULT_PRECISION, STRING_LENGTH
  use netcdf, only : nf90_global, nf90_nowrite, nf90_inquire_attribute, nf90_open, nf90_inq_dimid, nf90_inquire_dimension, &
       nf90_inq_varid, nf90_get_var, nf90_get_att, nf90_close
  use netcdf_misc_mod, only : check_netcdf_status
  use collections_mod, only : hashmap_type, c_put_real
  use conversions_mod, only : conv_to_string 
  use writer_types_mod, only : writer_type, unserialise_writer_type
  use logging_mod, only : LOG_ERROR, LOG_WARN, log_log
  use mpi, only : mpi_comm_rank, mpi_comm_size
  use data_utils_mod, only : unpack_scalar_integer_from_bytedata, unpack_scalar_dp_real_from_bytedata
  use timeaveraged_time_manipulation_mod, only : unserialise_time_averaged_state
  use instantaneous_time_manipulation_mod, only : unserialise_instantaneous_state
  implicit none

#ifndef TEST_MODE
  private
#endif

  abstract interface
     subroutine writer_field_manager_unserialise_state(byte_data)
       character, dimension(:), intent(in) :: byte_data
     end subroutine writer_field_manager_unserialise_state
  end interface

  character, dimension(:), allocatable :: serialised_writer_bytes, serialised_timeaveraged_bytes, &
       serialised_instantaneous_bytes, serialised_writer_field_manager_bytes, serialised_writer_entries_time_points
  integer :: my_io_server_rank

  public read_io_server_state, restart_writer_state_from_checkpoint, restart_timeaveraged_state_from_checkpoint, &
       restart_instantaneous_state_from_checkpoint, restart_writer_field_manager_from_checkpoint, &
       restart_writer_state_timepoints
contains

  !> Reads the IO server state
  !! @param checkpoint_filename The checkpoint filename to open and read from
  !! @param io_xml_configuration XML configuration is read from the checkpoint and placed into here
  subroutine read_io_server_state(checkpoint_filename, io_xml_configuration, io_communicator_arg)
    character(len=STRING_LENGTH), intent(in) :: checkpoint_filename
    character, dimension(:), allocatable, intent(inout) :: io_xml_configuration
    integer, intent(in) :: io_communicator_arg

    integer :: ncid, number_io_server, ierr
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
    call get_io_server_serialised_bytes(ncid, number_io_server, my_io_server_rank, &
         "serialised_writer_entry", serialised_writer_bytes)
    call get_io_server_serialised_bytes(ncid, number_io_server, my_io_server_rank, &
         "serialised_timeaveraged_manipulation", serialised_timeaveraged_bytes)
    call get_io_server_serialised_bytes(ncid, number_io_server, my_io_server_rank, &
         "serialised_instantaneous_manipulation", serialised_instantaneous_bytes)
    call get_io_server_serialised_bytes(ncid, number_io_server, my_io_server_rank, &
         "serialised_writer_manager", serialised_writer_field_manager_bytes)
    call get_io_server_serialised_bytes(ncid, number_io_server, my_io_server_rank, &
         "serialised_timepoints", serialised_writer_entries_time_points)
    call check_netcdf_status(nf90_close(ncid))
  end subroutine read_io_server_state

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

  !> Retrieves some IO server serialised bytes which will make up the state of a specific facet
  !! @param ncid The NetCDF checkpoint file ID
  !! @param number_io_server The number of IO servers
  !! @param my_io_server_rank My IO server rank
  !! @param base_key The base look up key to retrieve the state which is stored in the file
  !! @param raw_bytes The retrieved raw bytes as they sit in the file, this is allocated to hold the required data
  subroutine get_io_server_serialised_bytes(ncid, number_io_server, my_io_server_rank, base_key, raw_bytes)
    integer, intent(in) :: ncid, number_io_server, my_io_server_rank
    character(len=*), intent(in) :: base_key
    character, dimension(:), allocatable, intent(out) :: raw_bytes
    
    integer :: dim_id, var_id, dim_size, serialised_range(2), number_serialised_entries
    logical :: found

    call check_netcdf_status(nf90_inq_dimid(ncid, trim(base_key)//"_dim", dim_id), found)
    if (.not. found) return
    call check_netcdf_status(nf90_inquire_dimension(ncid, dim_id, len=dim_size))

    call check_netcdf_status(nf90_inq_varid(ncid, trim(base_key)//"_directory", var_id), found)
    if (.not. found) return
    if (my_io_server_rank .lt. number_io_server-1) then
      call check_netcdf_status(nf90_get_var(ncid, var_id, serialised_range, start=(/my_io_server_rank+1/), count=(/2/)))
    else
      call check_netcdf_status(nf90_get_var(ncid, var_id, serialised_range(1), start=(/my_io_server_rank+1/)))
      serialised_range(2)=dim_size
    end if
    if (serialised_range(2) .gt. dim_size) then
      call log_log(LOG_ERROR, "Serialised entry beyond size in the file")
    end if
    number_serialised_entries=(serialised_range(2)-serialised_range(1)) + 1    
    call check_netcdf_status(nf90_inq_varid(ncid, trim(base_key), var_id), found)
    if (.not. found) return
    allocate(raw_bytes(number_serialised_entries))
    call check_netcdf_status(nf90_get_var(ncid, var_id, raw_bytes, start=(/serialised_range(1)/),&
         count=(/number_serialised_entries/)))
  end subroutine get_io_server_serialised_bytes

  !> Restarts the writer state from a specific checkpoint byte data chunk of memory
  !! @param writer_entries The array of writer entries to fill in
  subroutine restart_writer_state_from_checkpoint(writer_entries)
    type(writer_type), dimension(:) :: writer_entries

    integer :: i, number_entries, current_point, byte_size
    
    if (.not. allocated(serialised_writer_bytes)) then
      if (my_io_server_rank == 0) then
        call log_log(LOG_WARN, "On restart no writer state in checkpoint file")
      end if
      return
    end if
    current_point=1
    number_entries=unpack_scalar_integer_from_bytedata(serialised_writer_bytes, current_point)
    if (number_entries .ne. size(writer_entries)) then
      call log_log(LOG_ERROR, "On restart have a different number of configured entries than those in the checkpoint file")
    end if
    do i=1, size(writer_entries)
      byte_size=unpack_scalar_integer_from_bytedata(serialised_writer_bytes, current_point)
      call unserialise_writer_type(writer_entries(i), serialised_writer_bytes(current_point:current_point+byte_size-1))
      current_point=current_point+byte_size
    end do    
    deallocate(serialised_writer_bytes)
  end subroutine restart_writer_state_from_checkpoint

  !> Restarts the writer state timepoints held in the writer federator
  !! @param time_points The timepoints hashmap which is filled in from the serialised version
  subroutine restart_writer_state_timepoints(time_points)
    type(hashmap_type), volatile, intent(inout) :: time_points

    integer :: i, number_entries, current_point, byte_size, timestep_key
    real(kind=DEFAULT_PRECISION) :: r_value

    if (.not. allocated(serialised_writer_entries_time_points)) then
      if (my_io_server_rank == 0) then
        call log_log(LOG_WARN, "On restart no writer state timepoints in checkpoint file")
      end if
      return
    end if
    current_point=1
    number_entries=unpack_scalar_integer_from_bytedata(serialised_writer_entries_time_points, current_point)
    do i=1, number_entries
      timestep_key=unpack_scalar_integer_from_bytedata(serialised_writer_entries_time_points, current_point)
      r_value=unpack_scalar_dp_real_from_bytedata(serialised_writer_entries_time_points, current_point)
      call c_put_real(time_points, conv_to_string(timestep_key), r_value)
    end do
    deallocate(serialised_writer_entries_time_points)
  end subroutine restart_writer_state_timepoints  

  !> Will restart the time averaged manipulation state from the checkpoint file
  subroutine restart_timeaveraged_state_from_checkpoint()
    if (.not. allocated(serialised_timeaveraged_bytes)) then
      if (my_io_server_rank == 0) then
        call log_log(LOG_WARN, "On restart no time averaged state in checkpoint file")
      end if
      return
    end if
    call unserialise_time_averaged_state(serialised_timeaveraged_bytes)
    deallocate(serialised_timeaveraged_bytes)
  end subroutine restart_timeaveraged_state_from_checkpoint

  !> Will restart the instantaneous manipulation state from the checkpoint file
  subroutine restart_instantaneous_state_from_checkpoint()
    if (.not. allocated(serialised_instantaneous_bytes)) then
      if (my_io_server_rank == 0) then
        call log_log(LOG_WARN, "On restart no instantaneous state in checkpoint file")
      end if
      return
    end if
    call unserialise_instantaneous_state(serialised_instantaneous_bytes)
    deallocate(serialised_instantaneous_bytes)
  end subroutine restart_instantaneous_state_from_checkpoint

  !> Will restart the field manager state from the checkpoint file
  !! @param unserialise_writer_field_manager The unserialise field manager procedure, done this way due to module ordering
  subroutine restart_writer_field_manager_from_checkpoint(unserialise_writer_field_manager)
    procedure(writer_field_manager_unserialise_state) :: unserialise_writer_field_manager

    if (.not. allocated(serialised_writer_field_manager_bytes)) then
      if (my_io_server_rank == 0) then
        call log_log(LOG_WARN, "On restart no writer field manager state in checkpoint file")
      end if
      return
    end if
    call unserialise_writer_field_manager(serialised_writer_field_manager_bytes)
    deallocate(serialised_writer_field_manager_bytes)
  end subroutine restart_writer_field_manager_from_checkpoint  
end module io_server_state_reader_mod
