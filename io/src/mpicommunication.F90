!> Abstraction layer around MPI, this issues and marshals the lower level communication details
module mpi_communication_mod
  use datadefn_mod, only : STRING_LENGTH
  use collections_mod, only : map_type, c_put
  use conversions_mod, only : conv_to_generic
  use configuration_parser_mod, only : io_configuration_data_definition_type, io_configuration_inter_communication_description
  use logging_mod, only : LOG_ERROR, log_log
  use io_server_client_mod, only : ARRAY_FIELD_TYPE, MAP_FIELD_TYPE, INTEGER_DATA_TYPE, BOOLEAN_DATA_TYPE, STRING_DATA_TYPE, &
       FLOAT_DATA_TYPE, DOUBLE_DATA_TYPE, COMMAND_TAG, DATA_TAG, INTER_IO_COMMUNICATION, get_data_description_from_name, &
       data_sizing_description_type, populate_mpi_type_extents, append_mpi_datatype, get_mpi_datatype_from_internal_representation
  use forthread_mod, only : forthread_mutex_lock, forthread_mutex_unlock
  use threadpool_mod, only : check_thread_status
  use mpi, only : MPI_COMM_WORLD, MPI_SOURCE, MPI_INT, MPI_BYTE, MPI_STATUS_SIZE, MPI_REQUEST_NULL, &
       MPI_STATUS_IGNORE, MPI_ANY_SOURCE
  implicit none

#ifndef TEST_MODE
  private
#endif

  integer :: command_buffer,& !< Buffer used to receive the command data into when it arrives on that channel
       command_request_handle !< Request handle representing the asynchronous P2P command request

  public build_mpi_datatype, data_receive, test_for_command, register_command_receive, &
       cancel_requests, free_mpi_type, get_number_io_servers, get_my_io_rank, test_for_inter_io
contains

  !> Retrieves the number of IO servers that are running in total
  !! @param io_comm The IO server communicator
  !! @returns The number of running IO servers
  integer function get_number_io_servers(io_comm)
    integer, intent(in) :: io_comm

    integer :: number, ierr

    call mpi_comm_size(io_comm, number, ierr)
    get_number_io_servers=number
  end function get_number_io_servers

  !> Retrieves my IO server rank out of the number of IO servers that are running
  !! @param io_comm The IO server communicator
  !! @returns My IO server rank
  integer function get_my_io_rank(io_comm)
    integer, intent(in) :: io_comm

    integer :: number, ierr

    call mpi_comm_rank(io_comm, number, ierr)
    get_my_io_rank=number
  end function get_my_io_rank  

  !> Registers a request for receiving a command from any MONC process on the command channel
  subroutine register_command_receive()
    integer :: ierr

    call mpi_irecv(command_buffer, 1, MPI_INT, MPI_ANY_SOURCE, COMMAND_TAG, &
         MPI_COMM_WORLD, command_request_handle, ierr)
  end subroutine register_command_receive

  !> Awaits some data on the data channel. This is of the type, size from the source provided and can either be written into
  !! a byte buffer or integer buffer depending upon the arguments provided
  !! @param mpi_datatype The MPI type of the data we are receiving
  !! @param size Number of elements to receive
  !! @param source The PID of the MONC process to receieve this data from
  !! @param dump_data (Optional) byte data buffer, for the data dump
  !! @param description_data (Optional) integer data buffer, for data description
  integer function data_receive(mpi_datatype, num_elements, source, dump_data, data_dump_id, description_data)
    integer, intent(in) :: mpi_datatype, num_elements, source
    integer, intent(in), optional :: data_dump_id
    character, dimension(:), allocatable, intent(inout), optional :: dump_data
    type(data_sizing_description_type), dimension(:), intent(inout), optional :: description_data
    integer :: ierr, status(MPI_STATUS_SIZE), recv_count, tag_to_use

    if (present(dump_data)) then
      tag_to_use=DATA_TAG
      if (present(data_dump_id)) tag_to_use=tag_to_use+data_dump_id
      call mpi_recv(dump_data, num_elements, mpi_datatype, source, tag_to_use, MPI_COMM_WORLD, status, ierr)
      call mpi_get_count(status, mpi_datatype, recv_count, ierr)
      data_receive=recv_count
    else if (present(description_data)) then
      call mpi_recv(description_data, num_elements, mpi_datatype, source, DATA_TAG, MPI_COMM_WORLD, status, ierr)
      call mpi_get_count(status, mpi_datatype, recv_count, ierr)
      data_receive=recv_count
    end if    
  end function data_receive

  !> Cancels all outstanding communication requests
  subroutine cancel_requests()
    call cancel_request(command_request_handle)
  end subroutine cancel_requests  

  !> Cancels a specific communication request
  !! @param req Handle of the request to cancel
  subroutine cancel_request(req)
    integer, intent(in) :: req

    integer :: ierr

    if (req .ne. MPI_REQUEST_NULL) call mpi_cancel(req, ierr)
  end subroutine cancel_request  

  !> Tests for a command message based upon the request already registered
  !! @param command The command which is received is returned to the caller
  !! @param source The PID of the source MONC process is returned to the caller
  !! @returns Whether a message is received or not
  logical function test_for_command(command, source)
    integer, intent(out) :: command, source

    integer :: ierr, status(MPI_STATUS_SIZE), complete

    call mpi_test(command_request_handle, complete, status, ierr)

    if (complete .eq. 1) then
      command = command_buffer
      source = status(MPI_SOURCE)
      call register_command_receive()
      test_for_command=.true.
    else
      test_for_command=.false.
    end if    
  end function test_for_command

  !> Tests for inter IO server communication
  !! @param inter_io_communications Data structures representing the possible inter IO communication sources
  !! @param number_of_inter_io Number of inter IO communication descriptions registered
  !! @param command The command which is received is returned to the caller
  !! @param source The source of the inter IO communication, which is set to the index of the inter descriptor
  !! @returns Whether a message is received or not
  logical function test_for_inter_io(inter_io_communications, number_of_inter_io, io_communicator, command, source, data_buffer)
    integer, intent(in) :: number_of_inter_io, io_communicator
    integer, intent(out) :: command, source
    type(io_configuration_inter_communication_description), dimension(:), intent(inout) :: inter_io_communications
    character, dimension(:), allocatable, intent(inout) :: data_buffer

    integer :: i, completed, ierr, status(MPI_STATUS_SIZE), message_size
    logical :: message_pending

    do i=1, number_of_inter_io
      call mpi_iprobe(MPI_ANY_SOURCE, inter_io_communications(i)%message_tag, io_communicator, message_pending, status, ierr)
      if (message_pending) then
        call mpi_get_count(status, MPI_BYTE, message_size, ierr)
        allocate(data_buffer(message_size))
        call mpi_recv(data_buffer, message_size, MPI_BYTE, MPI_ANY_SOURCE, inter_io_communications(i)%message_tag, &
             io_communicator, MPI_STATUS_IGNORE, ierr)
        command=INTER_IO_COMMUNICATION
        source=i
        test_for_inter_io=.true.
        return
      end if
    end do
    test_for_inter_io=.false.
  end function test_for_inter_io
  
  !> Frees an MPI type, used in clean up
  !! @param the_type The MPI type to free up
  subroutine free_mpi_type(the_type)
    integer, intent(in) :: the_type

    integer :: ierr

    call mpi_type_free(the_type, ierr)
  end subroutine free_mpi_type  

  !> Builds the MPI type that corresponds to the data which will be received from a specific MONC process. Two factors
  !! determine the structure and size of this - the XML configuration which has been parsed and also specific details
  !! of array sizes sent by each process as part of its registration process
  !! @param io_configuration IO server representation of the configuration which contains data structure layout
  !! @param array_sizes Sizes of each data array which has been received from a MONC process when it registeres
  !! @param data_size The data size corresponding to this type is returned (i.e. the buffer size in bytes required to hold it)
  !! @param field_start_locations For the MONC process the start location for each field, keyed on field name
  !! @param field_end_locations For the MONC process the end location for each field, keyed on field name
  !! @param field_dimensions Optional map of dimensions, if provided will store the number of dimensions for each field
  !! @returns The MPI data type representation
  integer function build_mpi_datatype(data_definition, data_size_info, data_size, field_start_locations, &
       field_end_locations, field_dimensions)
    type(io_configuration_data_definition_type), intent(in) :: data_definition
    type(data_sizing_description_type), dimension(:), intent(in) :: data_size_info
    integer, intent(out) :: data_size
    type(map_type), intent(out) :: field_start_locations, field_end_locations
    type(map_type), intent(out), optional :: field_dimensions

    integer :: type_extents(5), type_counts, i, j, field_start, data_type, field_array_sizes, &
         temp_size, prev_data_type, old_types(20), offsets(20), block_counts(20), new_type, current_location, ierr, field_ignores
    logical :: field_found
    type(data_sizing_description_type) :: field_size_info
    class(*), pointer :: generic

    type_extents=populate_mpi_type_extents()

    field_start=1
    data_type=0
    type_counts=0
    field_array_sizes=0
    field_ignores=0
    current_location=1
    do i=1,data_definition%number_of_data_fields
      if (data_type == 0) then        
        prev_data_type=data_type        
        data_type=data_definition%fields(i)%data_type
      else
        if (data_type .ne. data_definition%fields(i)%data_type) then
          ! For efficient type packing, combine multiple fields with the same type - therefore when the type changes work the previous one pack
          call append_mpi_datatype(field_start, i-1-field_ignores, field_array_sizes, data_type, &
               type_extents, prev_data_type, type_counts+1, old_types, offsets, block_counts)          
          field_start=i
          field_array_sizes=0
          field_ignores=0
          prev_data_type=data_type                   
          data_type=data_definition%fields(i)%data_type
          type_counts=type_counts+1
        end if
      end if
      generic=>conv_to_generic(current_location, .true.)
      call c_put(field_start_locations, data_definition%fields(i)%name, generic)
      if (data_definition%fields(i)%field_type .eq. ARRAY_FIELD_TYPE .or. &
           data_definition%fields(i)%field_type .eq. MAP_FIELD_TYPE) then
        ! Grab the field info based upon the field name
        field_found=get_data_description_from_name(data_size_info, data_definition%fields(i)%name, field_size_info)
        if (.not. field_found .or. field_size_info%dimensions == 0) then
          ! If no field info, or the dimension is 0 then this MONC process is not sending that field - check it is optional
          if (.not. data_definition%fields(i)%optional) then
            call log_log(LOG_ERROR, "Non optional field `"//trim(data_definition%fields(i)%name)//&
                 "' omitted from MONC IO server registration")
          end if
          field_ignores=field_ignores+1
        else
          ! If the field is specified then use the size data to assemble the field size and append to current size
          temp_size=1
          do j=1, field_size_info%dimensions
            temp_size=temp_size*field_size_info%dim_sizes(j)
          end do          
          if (data_definition%fields(i)%field_type .eq. MAP_FIELD_TYPE) then
            field_array_sizes=(field_array_sizes+temp_size*STRING_LENGTH*2)-1
            current_location=current_location+temp_size*STRING_LENGTH*2
          else
            field_array_sizes=(field_array_sizes+temp_size)-1
            current_location=current_location+temp_size*type_extents(data_type)
          end if
        end if
      else
        if (data_definition%fields(i)%optional) then
          if (get_data_description_from_name(data_size_info, data_definition%fields(i)%name)) then
            if (data_type==STRING_DATA_TYPE) then
              field_array_sizes=(field_array_sizes+STRING_LENGTH)-1
              current_location=current_location+type_extents(data_type)*STRING_LENGTH
            else
              current_location=current_location+type_extents(data_type)
            end if
          else
            field_ignores=field_ignores+1
          end if
        else
          if (data_type==STRING_DATA_TYPE) then
            field_array_sizes=(field_array_sizes+STRING_LENGTH)-1
            current_location=current_location+type_extents(data_type)*STRING_LENGTH
          else
            current_location=current_location+type_extents(data_type)
          end if
        end if
      end if
      generic=>conv_to_generic(current_location-1, .true.)
      call c_put(field_end_locations, data_definition%fields(i)%name, generic)
      if (present(field_dimensions)) then
        generic=>conv_to_generic(field_size_info%dimensions, .true.)
        call c_put(field_dimensions, data_definition%fields(i)%name, generic)
      end if
    end do
    if (field_start .le. i-1) then
      ! If there are outstanding fields to process then we do this here
      call append_mpi_datatype(field_start, i-1, field_array_sizes, data_type, &
               type_extents, prev_data_type, type_counts+1, old_types, offsets, block_counts)
      type_counts=type_counts+1
    end if
    call mpi_type_struct(type_counts, block_counts, offsets, old_types, new_type, ierr) 
    call mpi_type_commit(new_type, ierr)
    call mpi_type_size(new_type, data_size, ierr)
    build_mpi_datatype=new_type
  end function build_mpi_datatype
end module mpi_communication_mod
