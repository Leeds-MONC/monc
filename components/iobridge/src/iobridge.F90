!> Bridge between MONC and the IO server, this registers the current MONC process, will issue data dumps
!! and deregister MONCs when the model run is completed.
module iobridge_mod
  use monc_component_mod, only : COMPONENT_DOUBLE_DATA_TYPE, COMPONENT_INTEGER_DATA_TYPE, component_descriptor_type, &
       component_field_value_type, component_field_information_type
  use collections_mod, only : hashmap_type, map_type, list_type, c_contains, c_get_generic, c_get_string, c_put_generic, &
       c_put_integer, c_size, c_key_at, c_free
  use conversions_mod, only : conv_to_string
  use state_mod, only : model_state_type
  use grids_mod, only : X_INDEX, Y_INDEX, Z_INDEX, local_grid_type
  use optionsdatabase_mod, only : options_size, options_get_logical
  use prognostics_mod, only : prognostic_field_type
  use datadefn_mod, only : DEFAULT_PRECISION, SINGLE_PRECISION, DOUBLE_PRECISION, STRING_LENGTH
  use logging_mod, only : LOG_ERROR, LOG_WARN, log_log, log_master_log
  use optionsdatabase_mod, only : options_get_integer
  use q_indices_mod, only : q_metadata_type, get_indices_descriptor
  use registry_mod, only : get_all_component_published_fields, get_component_field_value, &
       get_component_field_information, is_component_enabled
  use io_server_client_mod, only : COMMAND_TAG, DATA_TAG, REGISTER_COMMAND, DEREGISTER_COMMAND, DATA_COMMAND_START, &
       ARRAY_FIELD_TYPE, SCALAR_FIELD_TYPE, MAP_FIELD_TYPE, INTEGER_DATA_TYPE, BOOLEAN_DATA_TYPE, STRING_DATA_TYPE, &
       FLOAT_DATA_TYPE, DOUBLE_DATA_TYPE, LOCAL_SIZES_KEY, LOCAL_START_POINTS_KEY, LOCAL_END_POINTS_KEY, NUMBER_Q_INDICIES_KEY, &
       data_sizing_description_type, definition_description_type, field_description_type, build_mpi_type_data_sizing_description,&
       build_mpi_type_field_description, build_mpi_type_definition_description, populate_mpi_type_extents, append_mpi_datatype, &
       get_mpi_datatype_from_internal_representation, pack_scalar_field, pack_array_field, pack_map_field
  use mpi, only : MPI_COMM_WORLD, MPI_INT, MPI_BYTE, MPI_REQUEST_NULL, MPI_STATUSES_IGNORE, MPI_STATUS_IGNORE, MPI_STATUS_SIZE
  use q_indices_mod, only : q_metadata_type, get_max_number_q_indices, get_indices_descriptor, get_number_active_q_indices

  use conditional_diagnostics_column_mod, only : ncond, ndiag, cond_request, diag_request, cond_long, diag_long


  implicit none

#ifndef TEST_MODE
  private
#endif

  type io_server_sendable_field_sizing
     integer :: number_dimensions, dimensions(4)
  end type io_server_sendable_field_sizing

  type io_configuration_field_type
     character(len=STRING_LENGTH) :: name
     integer :: field_type, data_type
     logical :: optional, enabled
  end type io_configuration_field_type

  type io_configuration_data_definition_type
     character(len=STRING_LENGTH) :: name
     logical :: send_on_terminate
     integer :: number_of_data_fields, frequency, mpi_datatype
     type(io_configuration_field_type), dimension(:), allocatable :: fields
     integer :: dump_requests(2) !< Dump non blocking send request handles
     character, dimension(:), allocatable :: send_buffer !< Send buffer which holds the model during a dump
  end type io_configuration_data_definition_type

  type(io_configuration_data_definition_type), dimension(:), allocatable :: data_definitions
  type(map_type) :: unique_field_names, sendable_fields, component_field_descriptions
  logical :: io_server_enabled, in_finalisation_callback

  public iobridge_get_descriptor

contains 

  type(component_descriptor_type) function iobridge_get_descriptor()
    iobridge_get_descriptor%name="iobridge"
    iobridge_get_descriptor%version=0.1
    iobridge_get_descriptor%initialisation=>init_callback
    iobridge_get_descriptor%timestep=>timestep_callback
    iobridge_get_descriptor%finalisation=>finalisation_callback
  end function iobridge_get_descriptor

  !> Initialisation call back, called at the start of the model run
  !! @param current_state The current model state
  subroutine init_callback(current_state)
    type(model_state_type), target, intent(inout) :: current_state

    integer :: mpi_type_data_sizing_description, mpi_type_definition_description, mpi_type_field_description, ierr

    if (.not. options_get_logical(current_state%options_database, "enable_io_server")) then
      io_server_enabled=.false.
      call log_master_log(LOG_WARN, "Enabled IO bridge but missing IO server compilation, therefore ignoring IO bridge component")
      return
    end if

    io_server_enabled=.true.
    in_finalisation_callback=.false.

    call populate_sendable_fields(current_state)

    mpi_type_data_sizing_description=build_mpi_type_data_sizing_description()
    mpi_type_definition_description=build_mpi_type_definition_description()
    mpi_type_field_description=build_mpi_type_field_description()

    call register_with_io_server(current_state, mpi_type_definition_description, mpi_type_field_description)
    call send_monc_specific_data_to_server(current_state, mpi_type_data_sizing_description)

    call mpi_type_free(mpi_type_data_sizing_description, ierr)
    call mpi_type_free(mpi_type_definition_description, ierr)
    call mpi_type_free(mpi_type_field_description, ierr)

    call build_mpi_data_types()
    
  end subroutine init_callback

  !> Model dump call back, called for each model dump phase
  !! @param current_state The current model state
  subroutine timestep_callback(current_state)
    type(model_state_type), target, intent(inout) :: current_state

    integer :: i

    if (.not. io_server_enabled) return

    do i=1, size(data_definitions)
      if (data_definitions(i)%frequency .gt. 0) then
        if (mod(current_state%timestep, data_definitions(i)%frequency) == 0) then
          call send_data_to_io_server(current_state, i)
        end if
      end if
    end do
  end subroutine timestep_callback

  !> Sends data to the IO server
  !! @param current_state The current model state
  !! @param data_index The specific data index to send over
  subroutine send_data_to_io_server(current_state, data_index)
    type(model_state_type), target, intent(inout) :: current_state
    integer, intent(in) :: data_index

    integer :: command_to_send, ierr

    if (data_definitions(data_index)%dump_requests(1) .ne. MPI_REQUEST_NULL .or. &
         data_definitions(data_index)%dump_requests(2) .ne. MPI_REQUEST_NULL) then
      ! Here wait for previous data dump to complete (consider extending to using buffers for performance)
      call mpi_waitall(2, data_definitions(data_index)%dump_requests, MPI_STATUSES_IGNORE, ierr)
    end if

    ! Pack the send buffer and send it to the IO server
    call pack_send_buffer(current_state, data_definitions(data_index))

    command_to_send=DATA_COMMAND_START+data_index
    call mpi_issend(command_to_send, 1, MPI_INT, current_state%parallel%corresponding_io_server_process, &
         COMMAND_TAG, MPI_COMM_WORLD, data_definitions(data_index)%dump_requests(1), ierr)
    call mpi_issend(data_definitions(data_index)%send_buffer, 1, data_definitions(data_index)%mpi_datatype, &
         current_state%parallel%corresponding_io_server_process, DATA_TAG+data_index, MPI_COMM_WORLD, &
         data_definitions(data_index)%dump_requests(2), ierr)
  end subroutine send_data_to_io_server  

  !> Finalisation call back, called at the end of the model run
  !! @param current_state The current model state
  subroutine finalisation_callback(current_state)
    type(model_state_type), target, intent(inout) :: current_state

    integer :: ierr, i

    if (.not. io_server_enabled) return
    in_finalisation_callback=.true.

    do i=1, size(data_definitions)      
      if (data_definitions(i)%send_on_terminate) then
        call send_data_to_io_server(current_state, i)
      end if
      if (data_definitions(i)%dump_requests(1) .ne. MPI_REQUEST_NULL .or. &
           data_definitions(i)%dump_requests(2) .ne. MPI_REQUEST_NULL) then
        call mpi_waitall(2, data_definitions(i)%dump_requests, MPI_STATUSES_IGNORE, ierr)
      end if
      if (allocated(data_definitions(i)%send_buffer)) deallocate(data_definitions(i)%send_buffer)
      call mpi_type_free(data_definitions(i)%mpi_datatype, ierr)
    end do
    call mpi_send(DEREGISTER_COMMAND, 1, MPI_INT, current_state%parallel%corresponding_io_server_process, &
         COMMAND_TAG, MPI_COMM_WORLD, ierr)
  end subroutine finalisation_callback

  !> Builds the MPI data types that correspond to the field descriptions and sizings
  subroutine build_mpi_data_types()
    integer :: i, dump_send_buffer_size

    do i=1, size(data_definitions)
      dump_send_buffer_size=build_mpi_data_type_for_definition(data_definitions(i))
      data_definitions(i)%dump_requests=(/MPI_REQUEST_NULL, MPI_REQUEST_NULL/)
      allocate(data_definitions(i)%send_buffer(dump_send_buffer_size))
    end do    
  end subroutine build_mpi_data_types

  !> Builds the MPI data type for a specific definition with sizing information
  !! @param specific_data_definition The data definition to build the type for
  !! @returns The size (in bytes) that the send buffer needs to be to store the data for the MPI operation
  integer function build_mpi_data_type_for_definition(specific_data_definition)
    type(io_configuration_data_definition_type), intent(inout) :: specific_data_definition

    integer :: type_extents(5), type_counts, i, j, tempsize, field_start, data_type, field_array_sizes, &
         temp_size, prev_data_type, old_types(20), offsets(20), block_counts(20), ierr, field_ignores
    logical :: field_found
    type(io_server_sendable_field_sizing) :: field_size_info

    type_extents=populate_mpi_type_extents()

    field_start=1
    data_type=0
    type_counts=0
    field_array_sizes=0
    field_ignores=0
    do i=1, specific_data_definition%number_of_data_fields
      if (data_type == 0) then
        prev_data_type=data_type        
        data_type=specific_data_definition%fields(i)%data_type
      else
        if (data_type .ne. specific_data_definition%fields(i)%data_type) then
          ! For efficient type packing, combine multiple fields with the same type - therefore when the type changes work the previous one pack
          call append_mpi_datatype(field_start, i-1-field_ignores, field_array_sizes, data_type, &
               type_extents, prev_data_type, type_counts+1, old_types, offsets, block_counts)
          field_start=i
          field_array_sizes=0
          field_ignores=0
          prev_data_type=data_type                   
          data_type=specific_data_definition%fields(i)%data_type
          type_counts=type_counts+1
        end if
      end if

      if (specific_data_definition%fields(i)%field_type .eq. ARRAY_FIELD_TYPE .or. &
           specific_data_definition%fields(i)%field_type .eq. MAP_FIELD_TYPE) then
        ! Grab the field info based upon the field name
        field_size_info=get_sendable_field_sizing(specific_data_definition%fields(i)%name, field_found)
        specific_data_definition%fields(i)%enabled=field_found
        if (.not. field_found .or. field_size_info%number_dimensions == 0) then
          ! If no field info, or the dimension is 0 then this MONC process is not sending that field - check it is optional
          if (.not. specific_data_definition%fields(i)%optional) then
            call log_log(LOG_ERROR, "Non optional field `"//trim(specific_data_definition%fields(i)%name)//&
                 "' omitted from MONC IO server registration")
          end if
          field_ignores=field_ignores+1
        else
          ! If the field is specified then use the size data to assemble the field size and append to current size
          temp_size=1
          do j=1, field_size_info%number_dimensions
            temp_size=temp_size*field_size_info%dimensions(j)
          end do
          if (specific_data_definition%fields(i)%field_type .eq. MAP_FIELD_TYPE) then
            field_array_sizes=(field_array_sizes+temp_size*STRING_LENGTH*2)-1
          else
            field_array_sizes=(field_array_sizes+temp_size)-1
          end if
        end if
      else
        if (specific_data_definition%fields(i)%optional) then
          field_size_info=get_sendable_field_sizing(specific_data_definition%fields(i)%name, field_found)
          specific_data_definition%fields(i)%enabled=field_found
          if (.not. field_found) field_ignores=field_ignores+1
        end if        
      end if
    end do
    if (field_start .le. i-1) then
      ! If there are outstanding fields to process then we do this here
      call append_mpi_datatype(field_start, i-1, field_array_sizes, data_type, &
               type_extents, prev_data_type, type_counts+1, old_types, offsets, block_counts)
      type_counts=type_counts+1
    end if

    call mpi_type_struct(type_counts, block_counts, offsets, old_types, specific_data_definition%mpi_datatype, ierr) 
    call mpi_type_commit(specific_data_definition%mpi_datatype, ierr)
    call mpi_type_size(specific_data_definition%mpi_datatype, tempsize, ierr)
    build_mpi_data_type_for_definition=tempsize
  end function build_mpi_data_type_for_definition

  !> Populates the sendable field definitions with the field sizing information
  !! @param current_state The current model state
  subroutine populate_sendable_fields(current_state)
    type(model_state_type), target, intent(inout) :: current_state

    type(list_type) :: published_field_descriptors
    integer :: i

    call populate_globally_visible_sendable_fields(current_state)
    published_field_descriptors=get_all_component_published_fields()
    do i=1, c_size(published_field_descriptors)
      call populate_component_public_field(current_state, c_get_string(published_field_descriptors, i))
    end do    
  end subroutine populate_sendable_fields

  !> Populates the field information for a specific publically available field offered by one of the components
  !! @param field_visibility_descriptor The field descriptor which contains sizing information
  subroutine populate_component_public_field(current_state, field_name)
    type(model_state_type), target, intent(inout) :: current_state
    character(len=*), intent(in) :: field_name

    class(*), pointer :: generic_data
    type(io_server_sendable_field_sizing), pointer :: field_sizing
    type(component_field_information_type) :: field_information
    type(component_field_information_type), pointer :: field_information_info_alloc


    field_information=get_component_field_information(current_state, field_name)
    if (field_information%enabled) then
      allocate(field_information_info_alloc, source=field_information)
      generic_data=>field_information_info_alloc
      call c_put_generic(component_field_descriptions, field_name, generic_data, .false.)

      allocate(field_sizing)
      field_sizing%number_dimensions=field_information%number_dimensions
      field_sizing%dimensions=field_information%dimension_sizes
      generic_data=>field_sizing
      call c_put_generic(sendable_fields, field_name, generic_data, .false.)
    end if
  end subroutine populate_component_public_field  

  !> Populates the globally visible sendable fields which is a key value pair mapping between name and description of that field
  !! @param current_state The current model state
  subroutine populate_globally_visible_sendable_fields(current_state)
    type(model_state_type), target, intent(inout) :: current_state

    integer :: x_size, y_size, z_size
    class(*), pointer :: raw_generic

    z_size=current_state%local_grid%size(Z_INDEX)
    y_size=current_state%local_grid%size(Y_INDEX)
    x_size=current_state%local_grid%size(X_INDEX)

    raw_generic=>generate_sendable_description(options_size(current_state%options_database))
    call c_put_generic(sendable_fields, "options_database", raw_generic, .false.)
    if (get_number_active_q_indices() .gt. 0) then
      raw_generic=>generate_sendable_description(get_number_active_q_indices())
      call c_put_generic(sendable_fields, "q_indicies", raw_generic, .false.)
    end if
    raw_generic=>generate_sendable_description(3)
    call c_put_generic(sendable_fields, "local_grid_size", raw_generic, .false.)
    call c_put_generic(sendable_fields, "local_grid_start", raw_generic, .false.)
#ifdef U_ACTIVE
    raw_generic=>generate_sendable_description()
    call c_put_generic(sendable_fields, "x_size", raw_generic, .false.)
    call c_put_generic(sendable_fields, "x_bottom", raw_generic, .false.)
    call c_put_generic(sendable_fields, "x_resolution", raw_generic, .false.)
    raw_generic=>generate_sendable_description(z_size, y_size, x_size)
    call c_put_generic(sendable_fields, "u", raw_generic, .false.)
    call c_put_generic(sendable_fields, "u_nogal", raw_generic, .false.)
    call c_put_generic(sendable_fields, "zu", raw_generic, .false.)
    if (allocated(current_state%global_grid%configuration%vertical%olubar)) then
      raw_generic=>generate_sendable_description(z_size)
      call c_put_generic(sendable_fields, "olubar", raw_generic, .false.)
      call c_put_generic(sendable_fields, "olzubar", raw_generic, .false.)
    end if
#endif
#ifdef V_ACTIVE
    raw_generic=>generate_sendable_description()
    call c_put_generic(sendable_fields, "y_size", raw_generic, .false.)
    call c_put_generic(sendable_fields, "y_bottom", raw_generic, .false.)
    call c_put_generic(sendable_fields, "y_resolution", raw_generic, .false.)
    raw_generic=>generate_sendable_description(z_size, y_size, x_size)
    call c_put_generic(sendable_fields, "v", raw_generic, .false.)
    call c_put_generic(sendable_fields, "v_nogal", raw_generic, .false.)
    call c_put_generic(sendable_fields, "zv", raw_generic, .false.)
    if (allocated(current_state%global_grid%configuration%vertical%olvbar)) then
      raw_generic=>generate_sendable_description(z_size)
      call c_put_generic(sendable_fields, "olvbar", raw_generic, .false.)
      call c_put_generic(sendable_fields, "olzvbar", raw_generic, .false.)
    end if
#endif
#ifdef W_ACTIVE
    raw_generic=>generate_sendable_description(z_size)
    call c_put_generic(sendable_fields, "z", raw_generic, .false.)
    call c_put_generic(sendable_fields, "thref", raw_generic, .false.)
    call c_put_generic(sendable_fields, "prefn", raw_generic, .false.)
    call c_put_generic(sendable_fields, "rhon", raw_generic, .false.)
    call c_put_generic(sendable_fields, "rho", raw_generic, .false.)
    raw_generic=>generate_sendable_description(z_size, y_size, x_size)
    call c_put_generic(sendable_fields, "w", raw_generic, .false.)
    call c_put_generic(sendable_fields, "zw", raw_generic, .false.)
#endif          
    if (current_state%number_q_fields .gt. 0) then
      raw_generic=>generate_sendable_description(z_size, y_size, x_size, current_state%number_q_fields)
      call c_put_generic(sendable_fields, "q", raw_generic, .false.)
      call c_put_generic(sendable_fields, "zq", raw_generic, .false.)
      if (allocated(current_state%global_grid%configuration%vertical%olqbar)) then
        raw_generic=>generate_sendable_description(z_size, current_state%number_q_fields)
        call c_put_generic(sendable_fields, "olqbar", raw_generic, .false.)
        call c_put_generic(sendable_fields, "olzqbar", raw_generic, .false.)
      end if
    end if
    if (current_state%th%active) then
      raw_generic=>generate_sendable_description(z_size, y_size, x_size)
      call c_put_generic(sendable_fields, "th", raw_generic, .false.)
      call c_put_generic(sendable_fields, "zth", raw_generic, .false.)
      if (allocated(current_state%global_grid%configuration%vertical%olthbar)) then
        raw_generic=>generate_sendable_description(z_size)
        call c_put_generic(sendable_fields, "olthbar", raw_generic, .false.)
        call c_put_generic(sendable_fields, "olzthbar", raw_generic, .false.)
      end if
    end if
    if (current_state%p%active) then
      raw_generic=>generate_sendable_description(z_size, y_size, x_size)
      call c_put_generic(sendable_fields, "p", raw_generic, .false.)   
    end if
    ! need to dump heating rate tendency from socrates radiation
    if (is_component_enabled(current_state%options_database, "socrates_couple")) then
       raw_generic=>generate_sendable_description(z_size, y_size, x_size)
       call c_put_generic(sendable_fields, "sth_lw", raw_generic, .false.)
       call c_put_generic(sendable_fields, "sth_sw", raw_generic, .false.)
    endif
       
    if (is_component_enabled(current_state%options_database, "pdf_analysis")) then
      if (allocated(current_state%global_grid%configuration%vertical%w_up)) then
        raw_generic=>generate_sendable_description(z_size)
        call c_put_generic(sendable_fields, "w_up", raw_generic, .false.)
      end if

      if (allocated(current_state%global_grid%configuration%vertical%w_dwn)) then
        raw_generic=>generate_sendable_description(z_size)
        call c_put_generic(sendable_fields, "w_dwn", raw_generic, .false.)
      end if
    end if 

  end subroutine populate_globally_visible_sendable_fields  

  !> Generates a sendable description based upon the dimension information supplied, missing arguments means that dimension
  !! does not exist
  !! @param dim1 Optional size of dimension one
  !! @param dim2 Optional size of dimension two
  !! @param dim3 Optional size of dimension three
  !! @param dim4 Optional size of dimension four
  !! @returns The corresponding sendable description for the field, with the number of dimensions and sizing of each of these
  function generate_sendable_description(dim1, dim2, dim3, dim4)
    integer, intent(in), optional :: dim1, dim2, dim3, dim4
    class(*), pointer :: generate_sendable_description

    type(io_server_sendable_field_sizing), pointer :: field
    integer :: number_dimensions

    allocate(field)
    number_dimensions=0
    field%dimensions=0
    if (present(dim1)) then
      field%dimensions(1)=dim1
      number_dimensions=number_dimensions+1
    end if
    if (present(dim2)) then
      field%dimensions(2)=dim2
      number_dimensions=number_dimensions+1
    end if
    if (present(dim3)) then
      field%dimensions(3)=dim3
      number_dimensions=number_dimensions+1
    end if
    if (present(dim4)) then
      field%dimensions(4)=dim4
      number_dimensions=number_dimensions+1
    end if
    field%number_dimensions=number_dimensions
    generate_sendable_description=>field
  end function generate_sendable_description

  !> Sends this MONC specific information to the IO server, which is field info (sizing & availability) as well as meta data
  !! such as ZN field and Q field names
  !! @param current_state The current model state
  !! @param mpi_type_data_sizing_description MPI data type representing the sizing message
  subroutine send_monc_specific_data_to_server(current_state, mpi_type_data_sizing_description)
    type(model_state_type), target, intent(inout) :: current_state
    integer, intent(in) :: mpi_type_data_sizing_description

    type(data_sizing_description_type), dimension(:), allocatable :: data_description
    character, dimension(:), allocatable :: buffer
    integer :: number_unique_fields, buffer_size, request_handles(2), ierr
    real(kind=DEFAULT_PRECISION) :: dreal

    number_unique_fields=c_size(unique_field_names)
    allocate(data_description(number_unique_fields+4))
    request_handles(1)=send_data_field_sizes_to_server(current_state, mpi_type_data_sizing_description, &
       data_description, number_unique_fields)
    buffer_size=(kind(dreal)*current_state%local_grid%size(Z_INDEX))*2 + (STRING_LENGTH * current_state%number_q_fields &
                 + 4*ncond*STRING_LENGTH + 2*ndiag*STRING_LENGTH )
    allocate(buffer(buffer_size))
    request_handles(2)=send_general_monc_information_to_server(current_state, buffer)
    call mpi_waitall(2, request_handles, MPI_STATUSES_IGNORE, ierr)
    deallocate(data_description)
    deallocate(buffer)
  end subroutine send_monc_specific_data_to_server  

  !> Assembles all the data field sizing information and sends this to the IO server
  !! @param current_state The current model state
  !! @param mpi_type_data_sizing_description MPI data type representing the sizing message
  !! @param data_description Data descriptions which will be populated and then sent
  !! @param number_unique_fields The number of unique fields that we are sending over
  integer function send_data_field_sizes_to_server(current_state, mpi_type_data_sizing_description, &
       data_description, number_unique_fields)
    type(model_state_type), target, intent(inout) :: current_state
    integer, intent(in) :: mpi_type_data_sizing_description, number_unique_fields
    type(data_sizing_description_type), dimension(:), intent(inout) :: data_description

    integer :: ierr, i, next_index, request_handle
    character(len=STRING_LENGTH) :: field_name
    
    call package_local_monc_decomposition_into_descriptions(current_state, data_description)
    next_index=5
    do i=1, number_unique_fields
      field_name=c_key_at(unique_field_names, i)
      if (c_contains(sendable_fields, field_name)) then
        call assemble_individual_description(data_description, next_index, field_name, get_sendable_field_sizing(field_name))
        next_index=next_index+1
      end if      
    end do    

    call mpi_isend(data_description, next_index-1, mpi_type_data_sizing_description, &
         current_state%parallel%corresponding_io_server_process, DATA_TAG, MPI_COMM_WORLD, request_handle, ierr)
    send_data_field_sizes_to_server=request_handle
  end function send_data_field_sizes_to_server

  !> Sends the general MONC information (ZN field and Q field names) to the IO server
  !! @param current_state The current model state
  !! @param buffer The communication buffer to use
  !! @returns Handle to nonblocking send
  integer function send_general_monc_information_to_server(current_state, buffer)
    type(model_state_type), target, intent(inout) :: current_state
    character, dimension(:), intent(inout) :: buffer
    
    character(len=STRING_LENGTH) :: q_field_name, cd_field_name
    type(q_metadata_type) :: q_meta_data
    integer :: current_loc, n, ierr, request_handle    
    
    current_loc=1
    current_loc=pack_array_field(buffer, current_loc, real_array_1d=current_state%global_grid%configuration%vertical%zn)
    if (current_state%number_q_fields .gt. 0) then
      do n=1, current_state%number_q_fields
        q_meta_data=get_indices_descriptor(n)       
        if (q_meta_data%l_used) then
          q_field_name=q_meta_data%name
        else
          q_field_name="qfield_"//trim(conv_to_string(n))
        end if        
        current_loc=pack_scalar_field(buffer, current_loc, string_value=q_field_name)
      end do
    end if
    current_loc=pack_array_field(buffer, current_loc, real_array_1d=current_state%global_grid%configuration%vertical%z)

    do n=1,ncond*2 
      if (n .le. ncond) then
        cd_field_name = cond_request(n)
        current_loc=pack_scalar_field(buffer, current_loc, string_value=cd_field_name)
       cd_field_name = cond_long(n)
        current_loc=pack_scalar_field(buffer, current_loc, string_value=cd_field_name)
      else
        cd_field_name = ".not. "//trim(cond_request(n-ncond))
        current_loc=pack_scalar_field(buffer, current_loc, string_value=cd_field_name)
        cd_field_name = ".not. "//trim(cond_long(n-ncond))
        current_loc=pack_scalar_field(buffer, current_loc, string_value=cd_field_name)
      end if
    end do    
    do n=1,ndiag
      cd_field_name = diag_request(n)
      current_loc=pack_scalar_field(buffer, current_loc, string_value=cd_field_name)
      cd_field_name = diag_long(n)
      current_loc=pack_scalar_field(buffer, current_loc, string_value=cd_field_name)
    end do


    call mpi_isend(buffer, current_loc-1, MPI_BYTE, current_state%parallel%corresponding_io_server_process, &
         DATA_TAG, MPI_COMM_WORLD, request_handle, ierr)    
    send_general_monc_information_to_server=request_handle
  end function send_general_monc_information_to_server  

  !> Packages the local MONC decomposition information into descriptions for communication
  !! @param current_state The current model state
  !! @param data_description THe data description to pack into
  subroutine package_local_monc_decomposition_into_descriptions(current_state, data_description)
    type(model_state_type), target, intent(inout) :: current_state
    type(data_sizing_description_type), dimension(:), intent(inout) :: data_description

    type(io_server_sendable_field_sizing) :: sizing_info

    sizing_info%number_dimensions=3
    sizing_info%dimensions(Z_INDEX)=current_state%local_grid%size(Z_INDEX)
    sizing_info%dimensions(Y_INDEX)=current_state%local_grid%size(Y_INDEX)
    sizing_info%dimensions(X_INDEX)=current_state%local_grid%size(X_INDEX)
    call assemble_individual_description(data_description, 1, LOCAL_SIZES_KEY, sizing_info)
    sizing_info%dimensions(Z_INDEX)=current_state%local_grid%start(Z_INDEX)
    sizing_info%dimensions(Y_INDEX)=current_state%local_grid%start(Y_INDEX)
    sizing_info%dimensions(X_INDEX)=current_state%local_grid%start(X_INDEX)
    call assemble_individual_description(data_description, 2, LOCAL_START_POINTS_KEY, sizing_info)
    sizing_info%dimensions(Z_INDEX)=current_state%local_grid%end(Z_INDEX)
    sizing_info%dimensions(Y_INDEX)=current_state%local_grid%end(Y_INDEX)
    sizing_info%dimensions(X_INDEX)=current_state%local_grid%end(X_INDEX)
    call assemble_individual_description(data_description, 3, LOCAL_END_POINTS_KEY, sizing_info)
    sizing_info%number_dimensions=1
    sizing_info%dimensions(1)=get_number_active_q_indices()
    call assemble_individual_description(data_description, 4, NUMBER_Q_INDICIES_KEY, sizing_info)
  end subroutine package_local_monc_decomposition_into_descriptions  

  !> Retrieves the sizing information associated with a specific field
  !! @param field_name The field name to look up
  !! @param field_found Optional flag depicting whether the field was found or not
  type(io_server_sendable_field_sizing) function get_sendable_field_sizing(field_name, field_found)
    character(len=*), intent(in) :: field_name
    logical, intent(out), optional :: field_found

    class(*), pointer :: generic

    if (present(field_found)) field_found=.false.
    if (c_contains(sendable_fields, field_name)) then
      generic=>c_get_generic(sendable_fields, field_name)
      select type(generic)
      type is (io_server_sendable_field_sizing)
        get_sendable_field_sizing=generic
        if (present(field_found)) field_found=.true.
      end select
    end if
  end function get_sendable_field_sizing

  !> Retrieves the descriptor associated with some component's field based upon the field name
  !! @param field_name The field name
  type(component_field_information_type) function get_component_field_descriptor(field_name)
    character(len=*), intent(in) :: field_name

    class(*), pointer :: generic
    if (c_contains(component_field_descriptions, field_name)) then
      generic=>c_get_generic(component_field_descriptions, field_name)
      select type(generic)
      type is (component_field_information_type)
        get_component_field_descriptor=generic
      end select
    end if
  end function get_component_field_descriptor   

  !> Will assemble an individual description of an array data field
  !! @param data_description The data structure used to describe the size of arrays
  !! @param index The index of this current field
  !! @param field_name The corresponding field name that we are describing
  !! @param dimensions The number of dimensions (zero means the field is inactive)
  !! @param dim1 Optional size of dimension 1
  !! @param dim2 Optional size of dimension 2
  !! @param dim3 Optional size of dimension 3
  !! @param dim4 Optional size of dimension 4
  subroutine assemble_individual_description(data_description, index, field_name, field_sizing_description)
    integer, intent(in) :: index
    character(len=*), intent(in) :: field_name
    type(io_server_sendable_field_sizing), intent(in) :: field_sizing_description
    type(data_sizing_description_type), dimension(:), intent(inout) :: data_description

    data_description(index)%field_name=field_name
    data_description(index)%dimensions=field_sizing_description%number_dimensions
    data_description(index)%dim_sizes=field_sizing_description%dimensions
  end subroutine assemble_individual_description 
  
  !> Registers this MONC with the corresponding IO server. This will encapsulate the entire protocol, which is sending the
  !! registration command, receiving the data and field definitions from the IO server and then sending back the sizing
  !! for the fields that this MONC will contribute.
  !! @param current_state The current model state
  !! @param mpi_type_definition_description MPI data type for data definition message
  !! @param mpi_type_field_description MPI data type for field definition message
  subroutine register_with_io_server(current_state, mpi_type_definition_description, mpi_type_field_description)
    type(model_state_type), target, intent(inout) :: current_state
    integer, intent(in) :: mpi_type_definition_description, mpi_type_field_description

    type(definition_description_type), dimension(:), allocatable :: definition_descriptions
    type(field_description_type), dimension(:), allocatable :: field_descriptions
    integer :: number_defns, number_fields, status(MPI_STATUS_SIZE), ierr

    call mpi_send(REGISTER_COMMAND, 1, MPI_INT, current_state%parallel%corresponding_io_server_process, &
         COMMAND_TAG, MPI_COMM_WORLD, ierr)

    call mpi_probe(current_state%parallel%corresponding_io_server_process, DATA_TAG, MPI_COMM_WORLD, status, ierr)

    call mpi_get_count(status, mpi_type_definition_description, number_defns, ierr)

    allocate(definition_descriptions(number_defns))

    call mpi_recv(definition_descriptions, number_defns, mpi_type_definition_description, &
         current_state%parallel%corresponding_io_server_process, DATA_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
    number_fields=get_total_number_of_fields(definition_descriptions, number_defns)

    allocate(field_descriptions(number_fields))
    call mpi_recv(field_descriptions, number_fields, mpi_type_field_description, &
         current_state%parallel%corresponding_io_server_process, DATA_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
    call populate_data_definition_configuration(definition_descriptions, number_defns, field_descriptions, number_fields)
    deallocate(definition_descriptions)
    
  end subroutine register_with_io_server

  !> Retrieve the total number of fields, which is all the fields in all the data definitions
  !! @param definition_descriptions Data definition descriptions
  !! @param number_defns The number of data definitions
  !! @returns The total number of fields
  integer function get_total_number_of_fields(definition_descriptions, number_defns)
    type(definition_description_type), dimension(:), intent(inout) :: definition_descriptions
    integer, intent(in) :: number_defns

    integer :: i

    get_total_number_of_fields=0
    do i=1, number_defns
      get_total_number_of_fields=get_total_number_of_fields+definition_descriptions(i)%number_fields
    end do    
  end function get_total_number_of_fields

  !> Based upon the received data and field definitions this will configure the IO bridge internal representation of these
  !! facets, which is a structured tree, data defintions holding their own fields rather than the unstructured data
  !! we get from the IO server
  !! @param definition_descriptions data definitions received from the IO server
  !! @param number_defns The number of data definitions
  !! @param field_descriptions The field descriptions that we received from the IO server
  !! @param number_fields The total number of field descriptions received
  subroutine populate_data_definition_configuration(definition_descriptions, number_defns, field_descriptions, number_fields)
    type(definition_description_type), dimension(:), intent(inout) :: definition_descriptions
    type(field_description_type), dimension(:), intent(inout) :: field_descriptions
    integer, intent(in) :: number_defns, number_fields

    integer :: i, definition_index, field_index

    allocate(data_definitions(number_defns))
    do i=1, number_defns
      data_definitions(i)%name=definition_descriptions(i)%definition_name
      data_definitions(i)%send_on_terminate=definition_descriptions(i)%send_on_terminate
      data_definitions(i)%number_of_data_fields=0 ! Will increment this for each field
      data_definitions(i)%frequency=definition_descriptions(i)%frequency
      allocate(data_definitions(i)%fields(definition_descriptions(i)%number_fields))
    end do
    do i=1, number_fields
      definition_index=get_definition_index(field_descriptions(i)%definition_name)
      field_index=data_definitions(definition_index)%number_of_data_fields+1
      data_definitions(definition_index)%number_of_data_fields=field_index
      data_definitions(definition_index)%fields(field_index)%name=field_descriptions(i)%field_name
      data_definitions(definition_index)%fields(field_index)%field_type=field_descriptions(i)%field_type      
      data_definitions(definition_index)%fields(field_index)%data_type=field_descriptions(i)%data_type
      data_definitions(definition_index)%fields(field_index)%optional=field_descriptions(i)%optional
      if (field_descriptions(i)%optional .or. field_descriptions(i)%field_type == ARRAY_FIELD_TYPE .or. &
           field_descriptions(i)%field_type == MAP_FIELD_TYPE) then
        call c_put_integer(unique_field_names, field_descriptions(i)%field_name, 1)
      end if
      if (.not. field_descriptions(i)%optional) data_definitions(definition_index)%fields(field_index)%enabled=.true.
    end do    
  end subroutine populate_data_definition_configuration

  !> Looks up a specific definition based upon its name and returns the index
  !! @param name The defintion name to search for
  !! @returns The index where the corresponding definition can be found or -1 if no definition is located with that name
  integer function get_definition_index(name)
    character(len=*), intent(in) :: name

    integer :: i
    do i=1, size(data_definitions)
      if (data_definitions(i)%name .eq. name) then
        get_definition_index=i
        return
      end if
    end do
    get_definition_index=-1
  end function get_definition_index  

  !> Packs the current state into the send buffer. This iterates through each field in the data description and adds it to the
  !! send buffer
  !! @param current_state The current model state
  !! @param data_definition The definition of the data which hold the send buffer and the fields
  subroutine pack_send_buffer(current_state, data_definition)
    type(model_state_type), target, intent(inout) :: current_state
    type(io_configuration_data_definition_type), intent(inout) :: data_definition

    integer :: current_buffer_point, i

    current_buffer_point=1
    do i=1, data_definition%number_of_data_fields
      if (data_definition%fields(i)%enabled) then
        if (data_definition%fields(i)%field_type == ARRAY_FIELD_TYPE) then
          current_buffer_point=pack_array_into_send_buffer(current_state, data_definition, data_definition%fields(i), &
               current_buffer_point)
        else if (data_definition%fields(i)%field_type == MAP_FIELD_TYPE) then
          current_buffer_point=pack_map_into_send_buffer(current_state, data_definition, data_definition%fields(i), &
               current_buffer_point)
        else if (data_definition%fields(i)%field_type == SCALAR_FIELD_TYPE) then
          current_buffer_point=pack_scalar_into_send_buffer(current_state, data_definition, data_definition%fields(i), &
               current_buffer_point)
        end if
      end if
    end do
  end subroutine pack_send_buffer

  !> Packs scalar fields into the send bufer
  !! @param current_state The current model state
  !! @param data_definition The data definition description
  !! @param field The specific field we are looking up
  !! @param current_buffer_point The current point in the buffer where this data will be entered
  !! @returns The new current buffer point which is after the data addition has taken place
  integer function pack_scalar_into_send_buffer(current_state, data_definition, field, current_buffer_point)
    type(model_state_type), target, intent(inout) :: current_state
    type(io_configuration_data_definition_type), intent(inout) :: data_definition
    type(io_configuration_field_type), intent(in) :: field
    integer, intent(in) :: current_buffer_point

    if (field%name .eq. "timestep") then
      pack_scalar_into_send_buffer=pack_scalar_field(data_definition%send_buffer, current_buffer_point, &
           int_value=current_state%timestep)
    else if (field%name .eq. "terminated") then
      pack_scalar_into_send_buffer=pack_scalar_field(data_definition%send_buffer, current_buffer_point, &
           logical_value=.not. current_state%continue_timestep .and. in_finalisation_callback)
    else if (field%name .eq. "z_size") then
      pack_scalar_into_send_buffer=pack_scalar_field(data_definition%send_buffer, current_buffer_point, &
           int_value=current_state%global_grid%size(Z_INDEX))
    else if (field%name .eq. "y_size") then
      pack_scalar_into_send_buffer=pack_scalar_field(data_definition%send_buffer, current_buffer_point, &
           int_value=current_state%global_grid%size(Y_INDEX))
    else if (field%name .eq. "y_bottom") then
      pack_scalar_into_send_buffer=pack_scalar_field(data_definition%send_buffer, current_buffer_point, &
           real_value=current_state%global_grid%bottom(Y_INDEX))
    else if (field%name .eq. "y_top") then
      pack_scalar_into_send_buffer=pack_scalar_field(data_definition%send_buffer, current_buffer_point, &
           real_value=current_state%global_grid%top(Y_INDEX))
    else if (field%name .eq. "y_resolution") then
      pack_scalar_into_send_buffer=pack_scalar_field(data_definition%send_buffer, current_buffer_point, &
           real_value=current_state%global_grid%resolution(Y_INDEX))
    else if (field%name .eq. "x_size") then
      pack_scalar_into_send_buffer=pack_scalar_field(data_definition%send_buffer, current_buffer_point, &
           int_value=current_state%global_grid%size(X_INDEX))
    else if (field%name .eq. "x_bottom") then
      pack_scalar_into_send_buffer=pack_scalar_field(data_definition%send_buffer, current_buffer_point, &
           real_value=current_state%global_grid%bottom(X_INDEX))
    else if (field%name .eq. "x_top") then
      pack_scalar_into_send_buffer=pack_scalar_field(data_definition%send_buffer, current_buffer_point, &
           real_value=current_state%global_grid%top(X_INDEX))
    else if (field%name .eq. "x_resolution") then
      pack_scalar_into_send_buffer=pack_scalar_field(data_definition%send_buffer, current_buffer_point, &
           real_value=current_state%global_grid%resolution(X_INDEX))
    else if (field%name .eq. "time") then
      ! The time is incremented with dtm as the model was about to increment for the next step and this is needed for diagnostics
      pack_scalar_into_send_buffer=pack_scalar_field(data_definition%send_buffer, current_buffer_point, &
           real_value=current_state%time+current_state%dtm)
    else if (field%name .eq. "ugal") then
      pack_scalar_into_send_buffer=pack_scalar_field(data_definition%send_buffer, current_buffer_point, &
           real_value=current_state%ugal)
    else if (field%name .eq. "vgal") then
      pack_scalar_into_send_buffer=pack_scalar_field(data_definition%send_buffer, current_buffer_point, &
           real_value=current_state%vgal)
    else if (field%name .eq. "nqfields") then
      pack_scalar_into_send_buffer=pack_scalar_field(data_definition%send_buffer, current_buffer_point, &
           int_value=current_state%number_q_fields)
    else if (field%name .eq. "dtm") then
      pack_scalar_into_send_buffer=pack_scalar_field(data_definition%send_buffer, current_buffer_point, &
           real_value=current_state%dtm)
    else if (field%name .eq. "dtm_new") then
      pack_scalar_into_send_buffer=pack_scalar_field(data_definition%send_buffer, current_buffer_point, &
           real_value=current_state%dtm_new)
    else if (field%name .eq. "absolute_new_dtm") then
      pack_scalar_into_send_buffer=pack_scalar_field(data_definition%send_buffer, current_buffer_point, &
           real_value=current_state%absolute_new_dtm)
    else if (field%name .eq. "rad_last_time") then
       pack_scalar_into_send_buffer=pack_scalar_field(data_definition%send_buffer, current_buffer_point, &
            real_value=current_state%rad_last_time)
    else
      ! Handle component field here
      pack_scalar_into_send_buffer=handle_component_field_scalar_packing_into_send_buffer(current_state, &
           data_definition, field, current_buffer_point)
    end if
  end function pack_scalar_into_send_buffer

  !> Packs a components field scalar into the send buffer, these are fields that are served up by components rather than
  !! explicitly available
  !! @param current_state The current model state
  !! @param data_definition The data definition description
  !! @param field The specific field we are looking up
  !! @param current_buffer_point The current point in the buffer where this data will be entered
  !! @returns The new current buffer point which is after the data addition has taken place
  integer function handle_component_field_scalar_packing_into_send_buffer(current_state, data_definition, &
       field, current_buffer_point)
    type(model_state_type), target, intent(inout) :: current_state
    type(io_configuration_data_definition_type), intent(inout) :: data_definition
    type(io_configuration_field_type), intent(in) :: field
    integer, intent(in) :: current_buffer_point

    type(component_field_information_type) :: field_descriptor
    type(component_field_value_type) :: published_value

    field_descriptor=get_component_field_descriptor(field%name)
    published_value=get_component_field_value(current_state, field%name)
    if (field_descriptor%data_type == COMPONENT_DOUBLE_DATA_TYPE) then
      handle_component_field_scalar_packing_into_send_buffer=pack_scalar_field(data_definition%send_buffer, &
           current_buffer_point, real_value=published_value%scalar_real)
    else if (field_descriptor%data_type == COMPONENT_INTEGER_DATA_TYPE) then
      handle_component_field_scalar_packing_into_send_buffer=pack_scalar_field(data_definition%send_buffer, &
           current_buffer_point, int_value=published_value%scalar_int)         
    end if
  end function handle_component_field_scalar_packing_into_send_buffer 

  !> Packs map fields into the send buffer
  !! @param current_state The current model state
  !! @param data_definition The data definition description
  !! @param field The specific field we are looking up
  !! @param current_buffer_point The current point in the buffer where this data will be entered
  !! @returns The new current buffer point which is after the data addition has taken place
  integer function pack_map_into_send_buffer(current_state, data_definition, field, current_buffer_point)
    type(model_state_type), target, intent(inout) :: current_state
    type(io_configuration_data_definition_type), intent(inout) :: data_definition
    type(io_configuration_field_type), intent(in) :: field
    integer, intent(in) :: current_buffer_point

    integer :: i
    type(q_metadata_type) :: specific_q_data
    type(hashmap_type) :: q_indicies_map

    if (field%name .eq. "options_database") then
      pack_map_into_send_buffer=pack_map_field(data_definition%send_buffer, current_buffer_point, current_state%options_database)
    else if (field%name .eq. "q_indicies") then
      do i=1, get_max_number_q_indices()
        specific_q_data=get_indices_descriptor(i)
        if (specific_q_data%l_used) then
          call c_put_integer(q_indicies_map, specific_q_data%name, i)
        end if
      end do
      pack_map_into_send_buffer=pack_map_field(data_definition%send_buffer, current_buffer_point, q_indicies_map)
      call c_free(q_indicies_map)
    end if    
  end function pack_map_into_send_buffer  

  !> Packs array fields into the send bufer
  !! @param current_state The current model state
  !! @param data_definition The data definition description
  !! @param field The specific field we are looking up
  !! @param current_buffer_point The current point in the buffer where this data will be entered
  !! @returns The new current buffer point which is after the data addition has taken place
  integer function pack_array_into_send_buffer(current_state, data_definition, field, current_buffer_point)
    type(model_state_type), target, intent(inout) :: current_state
    type(io_configuration_data_definition_type), intent(inout) :: data_definition
    type(io_configuration_field_type), intent(in) :: field
    integer, intent(in) :: current_buffer_point

    if (field%name .eq. "local_grid_size") then
      pack_array_into_send_buffer=pack_array_field(data_definition%send_buffer, current_buffer_point, &
           int_array=current_state%local_grid%size)
    else if (field%name .eq. "local_grid_start") then
      pack_array_into_send_buffer=pack_array_field(data_definition%send_buffer, current_buffer_point, &
           int_array=current_state%local_grid%start)
    else if (field%name .eq. "z") then
      pack_array_into_send_buffer=pack_array_field(data_definition%send_buffer, current_buffer_point, &
           real_array_1d=current_state%global_grid%configuration%vertical%z)
    else if (field%name .eq. "olubar") then
      pack_array_into_send_buffer=pack_array_field(data_definition%send_buffer, current_buffer_point, &
           real_array_1d=current_state%global_grid%configuration%vertical%olubar)
    else if (field%name .eq. "olzubar") then
      pack_array_into_send_buffer=pack_array_field(data_definition%send_buffer, current_buffer_point, &
           real_array_1d=current_state%global_grid%configuration%vertical%olzubar)
    else if (field%name .eq. "olvbar") then
      pack_array_into_send_buffer=pack_array_field(data_definition%send_buffer, current_buffer_point, &
           real_array_1d=current_state%global_grid%configuration%vertical%olvbar)
    else if (field%name .eq. "olzvbar") then
      pack_array_into_send_buffer=pack_array_field(data_definition%send_buffer, current_buffer_point, &
           real_array_1d=current_state%global_grid%configuration%vertical%olzvbar)
    else if (field%name .eq. "olthbar") then
      pack_array_into_send_buffer=pack_array_field(data_definition%send_buffer, current_buffer_point, &
           real_array_1d=current_state%global_grid%configuration%vertical%olthbar)
    else if (field%name .eq. "olzthbar") then
      pack_array_into_send_buffer=pack_array_field(data_definition%send_buffer, current_buffer_point, &
           real_array_1d=current_state%global_grid%configuration%vertical%olzthbar)
    else if (field%name .eq. "olqbar") then
      pack_array_into_send_buffer=pack_array_field(data_definition%send_buffer, current_buffer_point, &
           real_array_2d=current_state%global_grid%configuration%vertical%olqbar)
    else if (field%name .eq. "olzqbar") then
      pack_array_into_send_buffer=pack_array_field(data_definition%send_buffer, current_buffer_point, &
           real_array_2d=current_state%global_grid%configuration%vertical%olzqbar)
    else if (field%name .eq. "thref") then
      pack_array_into_send_buffer=pack_array_field(data_definition%send_buffer, current_buffer_point, &
           real_array_1d=current_state%global_grid%configuration%vertical%thref)
    else if (field%name .eq. "prefn") then
      pack_array_into_send_buffer=pack_array_field(data_definition%send_buffer, current_buffer_point, &
           real_array_1d=current_state%global_grid%configuration%vertical%prefn)
    else if (field%name .eq. "rhon") then
      pack_array_into_send_buffer=pack_array_field(data_definition%send_buffer, current_buffer_point, &
           real_array_1d=current_state%global_grid%configuration%vertical%rhon)
    else if (field%name .eq. "rho") then
      pack_array_into_send_buffer=pack_array_field(data_definition%send_buffer, current_buffer_point, &
           real_array_1d=current_state%global_grid%configuration%vertical%rho)
    else if (field%name .eq. "u") then
      current_state%u%data=current_state%u%data+current_state%ugal
      pack_array_into_send_buffer=pack_prognostic_flow_field(data_definition%send_buffer, current_state%u, &
           current_buffer_point, current_state%local_grid)
      current_state%u%data=current_state%u%data-current_state%ugal
    else if (field%name .eq. "u_nogal") then
      pack_array_into_send_buffer=pack_prognostic_flow_field(data_definition%send_buffer, current_state%u, current_buffer_point, &
           current_state%local_grid)
    else if (field%name .eq. "zu") then
      pack_array_into_send_buffer=pack_prognostic_flow_field(data_definition%send_buffer, current_state%zu, current_buffer_point,&
           current_state%local_grid)
    else if (field%name .eq. "v") then
      current_state%v%data=current_state%v%data+current_state%vgal
      pack_array_into_send_buffer=pack_prognostic_flow_field(data_definition%send_buffer, current_state%v, current_buffer_point, &
           current_state%local_grid)
      current_state%v%data=current_state%v%data-current_state%vgal
    else if (field%name .eq. "v_nogal") then
      pack_array_into_send_buffer=pack_prognostic_flow_field(data_definition%send_buffer, current_state%v, current_buffer_point, &
           current_state%local_grid)
    else if (field%name .eq. "zv") then
      pack_array_into_send_buffer=pack_prognostic_flow_field(data_definition%send_buffer, current_state%zv, current_buffer_point,&
           current_state%local_grid)
    else if (field%name .eq. "w") then
      pack_array_into_send_buffer=pack_prognostic_flow_field(data_definition%send_buffer, current_state%w, current_buffer_point, &
           current_state%local_grid)
    else if (field%name .eq. "zw") then
      pack_array_into_send_buffer=pack_prognostic_flow_field(data_definition%send_buffer, current_state%zw, current_buffer_point,&
           current_state%local_grid)
    else if (field%name .eq. "q") then
      pack_array_into_send_buffer=pack_q_fields(data_definition%send_buffer, current_state%q, current_state%number_q_fields, &
           current_buffer_point, current_state%local_grid)
    else if (field%name .eq. "zq") then
      pack_array_into_send_buffer=pack_q_fields(data_definition%send_buffer, current_state%zq, current_state%number_q_fields, &
           current_buffer_point, current_state%local_grid)
    else if (field%name .eq. "th") then
      pack_array_into_send_buffer=pack_prognostic_flow_field(data_definition%send_buffer, current_state%th, current_buffer_point,&
           current_state%local_grid)
    else if (field%name .eq. "zth") then
      pack_array_into_send_buffer=pack_prognostic_flow_field(data_definition%send_buffer, current_state%zth, current_buffer_point,&
           current_state%local_grid)
    else if (field%name .eq. "p") then
      pack_array_into_send_buffer=pack_prognostic_flow_field(data_definition%send_buffer, current_state%p, current_buffer_point, &
           current_state%local_grid)
    else if (field%name .eq. "sth_lw") then
      pack_array_into_send_buffer=pack_prognostic_flow_field(data_definition%send_buffer,  &
            current_state%sth_lw, current_buffer_point, current_state%local_grid) 
    else if (field%name .eq. "sth_sw") then
      pack_array_into_send_buffer=pack_prognostic_flow_field(data_definition%send_buffer,  &
            current_state%sth_sw, current_buffer_point, current_state%local_grid)
    else if (field%name .eq. "w_up") then
      pack_array_into_send_buffer=pack_array_field(data_definition%send_buffer, current_buffer_point, &
           real_array_1d=current_state%global_grid%configuration%vertical%w_up)
    else if (field%name .eq. "w_dwn") then
      pack_array_into_send_buffer=pack_array_field(data_definition%send_buffer, current_buffer_point, &
           real_array_1d=current_state%global_grid%configuration%vertical%w_dwn)
!!$    else if (field%name .eq. "sw_down_surf") then
!!$      pack_array_into_send_buffer=pack_prognostic_flow_field(data_definition%send_buffer,  &
!!$            current_state%sth_sw, current_buffer_point, current_state%local_grid)
!!$    else if (field%name .eq. "lww_down_surf") then
!!$      pack_array_into_send_buffer=pack_prognostic_flow_field(data_definition%send_buffer,  &
!!$            current_state%sth_sw, current_buffer_point, current_state%local_grid)  
    else
       ! Handle component field here
       pack_array_into_send_buffer=handle_component_field_array_packing_into_send_buffer(current_state, &
            data_definition, field, current_buffer_point)
    end if
  end function pack_array_into_send_buffer

  !> Packs a components field array into the send buffer, these are fields that are served up by components rather than
  !! explicitly available
  !! @param current_state The current model state
  !! @param data_definition The data definition description
  !! @param field The specific field we are looking up
  !! @param current_buffer_point The current point in the buffer where this data will be entered
  !! @returns The new current buffer point which is after the data addition has taken place
  integer function handle_component_field_array_packing_into_send_buffer(current_state, data_definition, &
       field, current_buffer_point)
    type(model_state_type), target, intent(inout) :: current_state
    type(io_configuration_data_definition_type), intent(inout) :: data_definition
    type(io_configuration_field_type), intent(in) :: field
    integer, intent(in) :: current_buffer_point

    type(component_field_information_type) :: field_descriptor
    type(component_field_value_type) :: published_value

    field_descriptor=get_component_field_descriptor(field%name)
    published_value=get_component_field_value(current_state, field%name)
    if (field_descriptor%data_type == COMPONENT_DOUBLE_DATA_TYPE) then
      if (field_descriptor%number_dimensions == 1) then
        handle_component_field_array_packing_into_send_buffer=pack_array_field(data_definition%send_buffer, &
             current_buffer_point, real_array_1d=published_value%real_1d_array)
        deallocate(published_value%real_1d_array)
      else if (field_descriptor%number_dimensions == 2) then
        handle_component_field_array_packing_into_send_buffer=pack_array_field(data_definition%send_buffer, &
             current_buffer_point, real_array_2d=published_value%real_2d_array)
        deallocate(published_value%real_2d_array)
      else if (field_descriptor%number_dimensions == 3) then
        handle_component_field_array_packing_into_send_buffer=pack_array_field(data_definition%send_buffer, &
             current_buffer_point, real_array_3d=published_value%real_3d_array)
        deallocate(published_value%real_3d_array)
      else if (field_descriptor%number_dimensions == 4) then
        handle_component_field_array_packing_into_send_buffer=pack_array_field(data_definition%send_buffer, &
             current_buffer_point, real_array_4d=published_value%real_4d_array)
        deallocate(published_value%real_4d_array)
      end if
    end if
  end function handle_component_field_array_packing_into_send_buffer

  !> Packs the data of a specific prognostic field into a buffer
  !! @param buffer The buffer to pack the field into
  !! @param prognostic The prognostic field  
  !! @param start_offset The starting offset to write into the buffer
  !! @param local_grid Description of the local grid
  !! @returns The next location in the buffer to write to (next start offset)
  integer function pack_prognostic_flow_field(buffer, prognostic, start_offset, local_grid)
    character, dimension(:), allocatable, intent(inout) :: buffer
    type(prognostic_field_type), intent(inout) :: prognostic
    integer, intent(in) :: start_offset
    type(local_grid_type), intent(inout) :: local_grid    

    integer :: target_end

    target_end=start_offset + (local_grid%size(Z_INDEX)*local_grid%size(Y_INDEX)*local_grid%size(X_INDEX)*kind(prognostic%data)-1)

    buffer(start_offset : target_end) = transfer(prognostic%data(&
         local_grid%local_domain_start_index(Z_INDEX): local_grid%local_domain_end_index(Z_INDEX),&
         local_grid%local_domain_start_index(Y_INDEX): local_grid%local_domain_end_index(Y_INDEX), &
         local_grid%local_domain_start_index(X_INDEX): local_grid%local_domain_end_index(X_INDEX)), &
         buffer(start_offset : target_end))
    pack_prognostic_flow_field=target_end+1
  end function pack_prognostic_flow_field  

  !> Packs the Q fields into the send buffer
  !! @param buffer The send buffer to pack into
  !! @param q_fields Q prognostic fields
  !! @param number_q_fields The number of Q fields
  !! @param start_offset Starting offset in the buffer to pack into
  !! @param local_grid Local grid description
  !! @returns Updated write location, which is the next location in the buffer to write to
  integer function pack_q_fields(buffer, q_fields, number_q_fields, start_offset, local_grid)
    character, dimension(:), allocatable, intent(inout) :: buffer
    type(prognostic_field_type), dimension(:), intent(inout) :: q_fields
    integer, intent(in) :: start_offset, number_q_fields
    type(local_grid_type), intent(inout) :: local_grid    

    integer :: target_end, i, current_starting_index

    current_starting_index=start_offset

    do i=1,number_q_fields
      target_end=current_starting_index + (local_grid%size(Z_INDEX)*local_grid%size(Y_INDEX)*&
           local_grid%size(X_INDEX)*kind(q_fields(i)%data)-1)
      buffer(current_starting_index : target_end) = transfer(q_fields(i)%data(&
           local_grid%local_domain_start_index(Z_INDEX): local_grid%local_domain_end_index(Z_INDEX),&
           local_grid%local_domain_start_index(Y_INDEX): local_grid%local_domain_end_index(Y_INDEX), &
           local_grid%local_domain_start_index(X_INDEX): local_grid%local_domain_end_index(X_INDEX)), &
           buffer(current_starting_index : target_end))
      current_starting_index=target_end+1      
    end do
    pack_q_fields=target_end+1
  end function pack_q_fields   
end module iobridge_mod
