!> The main IO server functionality which handles waiting for commands and data both of which are delt with.
!! The lower level details of the communication, configuration parsing etc are all held elsewhere. The server
!! can be thought of similar to a bus, with command and data channels. The command gives context to what is on
!! the data channel and not all commands require data (such as deregistration of MONC process)
module io_server_mod
  use datadefn_mod, only : DEFAULT_PRECISION, STRING_LENGTH, LONG_STRING_LENGTH
  use configuration_parser_mod, only : DATA_SIZE_STRIDE, io_configuration_type, io_configuration_data_definition_type, &
       io_configuration_registered_monc_type, configuration_parse, extend_registered_moncs_array, retrieve_data_definition, &
       build_definition_description_type_from_configuration, build_field_description_type_from_configuration, get_monc_location, &
       get_io_xml, cond_request, diag_request, cond_long, diag_long, ncond, ndiag, l_thoff

  use mpi_communication_mod, only : build_mpi_datatype, data_receive, test_for_command, register_command_receive, &
       cancel_requests, free_mpi_type, get_number_io_servers, get_my_io_rank, test_for_inter_io, lock_mpi, unlock_mpi, &
       waitall_for_mpi_requests, initialise_mpi_communication, pause_for_mpi_interleaving
  use diagnostic_federator_mod, only : initialise_diagnostic_federator, finalise_diagnostic_federator, &
       check_diagnostic_federator_for_completion, pass_fields_to_diagnostics_federator, determine_diagnostics_fields_available
  use writer_federator_mod, only : initialise_writer_federator, finalise_writer_federator, check_writer_for_trigger, &
       inform_writer_federator_fields_present, inform_writer_federator_time_point, provide_q_field_names_to_writer_federator, &
       provide_tracer_names_to_writer_federator, any_pending
  use writer_field_manager_mod, only : initialise_writer_field_manager, finalise_writer_field_manager, &
       provide_monc_data_to_writer_federator
  use collections_mod, only : hashset_type, hashmap_type, map_type, iterator_type, c_get_integer, c_put_integer, c_is_empty, &
       c_remove, c_add_string, c_integer_at, c_free, c_get_iterator, c_has_next, c_next_mapentry
  use conversions_mod, only : conv_to_string
  use string_utils_mod, only : replace_character
  use io_server_client_mod, only : REGISTER_COMMAND, DEREGISTER_COMMAND, INTER_IO_COMMUNICATION, DATA_COMMAND_START, DATA_TAG, &
       LOCAL_SIZES_KEY, LOCAL_START_POINTS_KEY, LOCAL_END_POINTS_KEY, NUMBER_Q_INDICIES_KEY, SCALAR_FIELD_TYPE, &
       data_sizing_description_type, definition_description_type, field_description_type, build_mpi_type_data_sizing_description,&
       get_data_description_from_name, build_mpi_type_field_description, build_mpi_type_definition_description
  use forthread_mod, only : forthread_rwlock_rdlock, forthread_rwlock_wrlock, forthread_rwlock_tryrdlock, &
       forthread_rwlock_unlock, forthread_rwlock_init, forthread_rwlock_destroy, forthread_mutex_init, forthread_mutex_lock, &
       forthread_mutex_unlock, forthread_cond_wait, forthread_cond_signal, forthread_cond_init
  use threadpool_mod, only : threadpool_init, threadpool_finalise, threadpool_start_thread, check_thread_status, &
       threadpool_deactivate, threadpool_is_idle
  use global_callback_inter_io_mod, only : perform_global_callback
  use logging_mod, only : LOG_ERROR, LOG_WARN, log_log, initialise_logging
  use optionsdatabase_mod, only : options_get_logical, options_get_string
  use mpi, only : MPI_COMM_WORLD, MPI_STATUSES_IGNORE, MPI_BYTE, MPI_INT
  use io_server_state_reader_mod, only : read_io_server_configuration
  implicit none

#ifndef TEST_MODE
  private
#endif  

  integer :: mpi_type_data_sizing_description, & !< The MPI type for field sizing (i.e. array size etc send when MONCs register)
       mpi_type_definition_description, & !< The MPI data type for data descriptions sent to MONCs
       mpi_type_field_description !< The MPI data type for field descriptions sent to MONCs
  type(io_configuration_type), volatile, save :: io_configuration !< Internal representation of the IO configuration
  logical, volatile :: continue_poll_messages, & !< Whether to continue waiting command messages from any MONC processes
       initialised_present_data
  logical, volatile :: continue_poll_interio_messages, already_registered_finishing_call
  type(field_description_type), dimension(:), allocatable :: registree_field_descriptions
  type(definition_description_type), dimension(:), allocatable :: registree_definition_descriptions
  integer, dimension(:,:), allocatable :: sample_output_pairs

  integer, volatile :: monc_registration_lock

  public io_server_run
contains

  !> Called to start the IO server and once this subroutine returns then it indicates that the IO server has finished.
  !! The runtime is spent in here awaiting commands and then dealing with them. Termination occurs when all MONC processes
  !! have deregistered, note that to trigger this then at least one MONC process must first register
  !! @param io_communicator_arg The IO communicator containing just the IO servers
  !! @param io_xml_configuration Textual XML configuration that is used to set up the IO server
  subroutine io_server_run(options_database, io_communicator_arg, &
       provided_threading, total_global_processes, continuation_run, reconfig_initial_time, io_configuration_file, &
       my_global_rank)
    type(hashmap_type), intent(inout) :: options_database
    integer, intent(in) :: io_communicator_arg, provided_threading, total_global_processes, my_global_rank
    logical, intent(in) :: continuation_run
    real(kind=DEFAULT_PRECISION), intent(in) :: reconfig_initial_time
    character(len=LONG_STRING_LENGTH), intent(in) :: io_configuration_file

    integer :: command, source, my_rank, ierr
    character, dimension(:), allocatable :: data_buffer, io_xml_configuration
    type(hashmap_type) :: diagnostic_generation_frequency


    if (continuation_run) then
      ! Handle case where we need to allocate this due to no IO server config
      call read_io_server_configuration(options_get_string(options_database, "checkpoint"), &
           io_xml_configuration, io_communicator_arg)
    end if

    if (.not. allocated(io_xml_configuration)) then
      io_xml_configuration=get_io_xml(io_configuration_file)
      if (continuation_run) then
        call mpi_comm_rank(io_communicator_arg, my_rank, ierr)
        if (my_rank == 0) then
          call log_log(LOG_WARN, "No IO server configuration in checkpoint file - starting from XML provided file instead")
        end if
      end if      
    end if    

    call check_for_condi_conflict(io_xml_configuration, options_database)
    call configuration_parse(options_database, io_xml_configuration, io_configuration)
    deallocate(io_xml_configuration)
    call threadpool_init(io_configuration)
    call initialise_mpi_communication(provided_threading)
    call check_thread_status(forthread_rwlock_init(monc_registration_lock, -1))
    call check_thread_status(forthread_mutex_init(io_configuration%general_info_mutex, -1))
    initialised_present_data=.false.
    continue_poll_messages=.true.
    continue_poll_interio_messages=.true.
    already_registered_finishing_call=.false.
    io_configuration%io_communicator=io_communicator_arg
    io_configuration%number_of_io_servers=get_number_io_servers(io_communicator_arg)
    io_configuration%number_of_global_moncs=total_global_processes-io_configuration%number_of_io_servers
    io_configuration%my_io_rank=get_my_io_rank(io_communicator_arg)
    io_configuration%my_global_rank=my_global_rank
    call initialise_logging(io_configuration%my_io_rank)
    registree_definition_descriptions=build_definition_description_type_from_configuration(io_configuration)
    registree_field_descriptions=build_field_description_type_from_configuration(io_configuration)
    diagnostic_generation_frequency=initialise_diagnostic_federator(io_configuration)
    call initialise_writer_federator(io_configuration, diagnostic_generation_frequency, continuation_run, &
                                     reconfig_initial_time, sample_output_pairs)
    call c_free(diagnostic_generation_frequency)
    call initialise_writer_field_manager(io_configuration, continuation_run, reconfig_initial_time)
    mpi_type_data_sizing_description=build_mpi_type_data_sizing_description()
    mpi_type_definition_description=build_mpi_type_definition_description()
    mpi_type_field_description=build_mpi_type_field_description()

    call register_command_receive()

    do while (await_command(command, source, data_buffer))      
      call handle_command_message(command, source, data_buffer)
    end do
    call threadpool_deactivate()
    call finalise_writer_field_manager()
    call finalise_writer_federator()
    call finalise_diagnostic_federator(io_configuration)
    call check_thread_status(forthread_rwlock_destroy(monc_registration_lock))
    call free_individual_registered_monc_aspects()
    call cancel_requests()
    call free_mpi_type(mpi_type_data_sizing_description)
    call free_mpi_type(mpi_type_definition_description)
    call free_mpi_type(mpi_type_field_description)    
    call threadpool_finalise()
  end subroutine io_server_run

  !> Handle potential conditional diagnostics conflict
  !! Provides a more helpful error in the case where conditional diagnostics are requested as output,
  !! but their components are not enabled.
  !! We check this by searching the io_xml_configuration.
  !! @param raw_contents, intended to be the io_xml_configuration character array
  !! @param options_database
  subroutine check_for_condi_conflict(raw_contents, options_database)
    character, dimension(:), intent(in) :: raw_contents
    type(hashmap_type), intent(inout) :: options_database
    character(len=size(raw_contents)) :: string_to_process
    integer :: i

    if (.not. options_get_logical(options_database, "conditional_diagnostics_column_enabled")) then
      do i=1, size(raw_contents)
        string_to_process(i:i)=raw_contents(i)
      end do
      if (index(string_to_process,"CondDiags_") .ne. 0) then
        call log_log(LOG_ERROR, &
            "Conditional diagnostics are DISABLED but requested via xml.  Enable or remove request to resolve.")
      end if
    end if
  end subroutine check_for_condi_conflict

  !> Awaits a command or shutdown from MONC processes and other IO servers
  !! @param command The command received is output
  !! @param source The source process received is output
  !! @returns Whether to continue polling for commands (and whether to process the current output)
  logical function await_command(command, source, data_buffer)
    integer, intent(out) :: command, source
    character, dimension(:), allocatable :: data_buffer

    logical :: completed, inter_io_complete

    completed=.false.
    await_command=.false.
    do while(.not. completed)
      if (.not. continue_poll_messages .and. .not. continue_poll_interio_messages) return
      if (continue_poll_messages) then
        if (test_for_command(command, source)) then
          await_command=.true.
          return
        end if
      end if
      if (continue_poll_interio_messages .and. allocated(io_configuration%inter_io_communications)) then       
        inter_io_complete=test_for_inter_io(io_configuration%inter_io_communications, &
             io_configuration%number_inter_io_communications, io_configuration%io_communicator, command, source, data_buffer) 
        if (inter_io_complete) then
          await_command=.true.
          return
        end if
      end if
      if (.not. continue_poll_messages .and. .not. already_registered_finishing_call) then
        if (check_diagnostic_federator_for_completion(io_configuration) .and. &
            (.not. any_pending()) .and. threadpool_is_idle()) then
          already_registered_finishing_call=.true.          
          call perform_global_callback(io_configuration, "termination", 1, termination_callback)          
        end if
      end if  
      if (.not. completed) call pause_for_mpi_interleaving()
    end do    
  end function await_command

  !> This is the termination callback which is called once all MONCs have deregistered, no sends are active by inter IO
  !! communications and all threads are idle. This shuts down the inter IO listening and kickstarts finalisation and closure
  !! @param io_configuration The IO server configuration
  !! @param values Values (ignored)
  !! @param field_name Field name identifier
  !! @param timestep Timestep identifier
  subroutine termination_callback(io_configuration, values, field_name, timestep)
    type(io_configuration_type), intent(inout) :: io_configuration
    real(DEFAULT_PRECISION), dimension(:) :: values
    character(len=STRING_LENGTH) :: field_name
    integer :: timestep

    continue_poll_interio_messages=.false.
  end subroutine termination_callback  

  !> Called to handle a specific command that has been recieved
  !! @param command The command which has been received from some process
  !! @param source The PID of the source (MONC) process
  subroutine handle_command_message(command, source, data_buffer)
    integer, intent(in) :: command, source
    character, dimension(:), allocatable, intent(inout) :: data_buffer

    if (command == REGISTER_COMMAND) then
      if (l_thoff) then
        call handle_monc_registration((/ source /))
      else
        call threadpool_start_thread(handle_monc_registration, (/ source /))
      end if
    else if (command == DEREGISTER_COMMAND) then
      if (l_thoff) then
        call handle_deregistration_command((/ source /))
      else
        call threadpool_start_thread(handle_deregistration_command, (/ source /))
      end if
    else if (command == INTER_IO_COMMUNICATION) then
      if (l_thoff) then
        call handle_inter_io_communication_command((/ source /), data_buffer=data_buffer)
      else
        call threadpool_start_thread(handle_inter_io_communication_command, (/ source /), data_buffer=data_buffer)      
      end if
      deallocate(data_buffer)
    else if (command .ge. DATA_COMMAND_START) then      
      call pull_back_data_message_and_handle(source, command-DATA_COMMAND_START)
    end if    
  end subroutine handle_command_message

  !> Handles inter IO server communications
  !! @param arguments The thread based arguments, this is the index of the inter IO server description
  subroutine handle_inter_io_communication_command(arguments, data_buffer)
    integer, dimension(:), intent(in) :: arguments
    character, dimension(:), allocatable, intent(inout), optional :: data_buffer

    integer :: source

    source=arguments(1)

    call io_configuration%inter_io_communications(source)%handling_procedure(io_configuration, data_buffer, source)
  end subroutine handle_inter_io_communication_command

  !> Frees up the memory associated with individual registered MONCs. This is done at the end for all MONCs as we can't
  !! deallocate dynamically in a threaded environment without excessive ordering and locking in case some data processing
  !! is queued or in progress
  subroutine free_individual_registered_monc_aspects()
    integer :: i, specific_monc_data_type
    type(iterator_type) :: types_iterator

    do i=1, size(io_configuration%registered_moncs)
      types_iterator=c_get_iterator(io_configuration%registered_moncs(i)%registered_monc_types)
      do while (c_has_next(types_iterator))
        specific_monc_data_type=c_get_integer(c_next_mapentry(types_iterator))
        call free_mpi_type(specific_monc_data_type)
      end do      
      if (allocated(io_configuration%registered_moncs(i)%field_start_locations)) &
           deallocate(io_configuration%registered_moncs(i)%field_start_locations)
      if (allocated(io_configuration%registered_moncs(i)%field_end_locations)) &
           deallocate(io_configuration%registered_moncs(i)%field_end_locations)
      if (allocated(io_configuration%registered_moncs(i)%definition_names)) &
           deallocate(io_configuration%registered_moncs(i)%definition_names)
      if (allocated(io_configuration%registered_moncs(i)%dimensions)) deallocate(io_configuration%registered_moncs(i)%dimensions)
    end do
  end subroutine free_individual_registered_monc_aspects  

  !> Deregisteres a specific MONC source process
  !! @param source The MONC process PID that we are deregistering
  subroutine handle_deregistration_command(arguments, data_buffer)
    integer, dimension(:), intent(in) :: arguments
    character, dimension(:), allocatable, intent(inout), optional :: data_buffer

    integer :: monc_location, source

    source=arguments(1)
    monc_location=get_monc_location(io_configuration, source)
    call check_thread_status(forthread_mutex_lock(io_configuration%registered_moncs(monc_location)%active_mutex))
    do while (io_configuration%registered_moncs(monc_location)%active_threads .gt. 0)
      call check_thread_status(forthread_cond_wait(io_configuration%registered_moncs(monc_location)%deactivate_condition_variable,&
             io_configuration%registered_moncs(monc_location)%active_mutex))
    end do
    call check_thread_status(forthread_mutex_unlock(io_configuration%registered_moncs(monc_location)%active_mutex))
    call check_thread_status(forthread_rwlock_wrlock(monc_registration_lock))
    io_configuration%active_moncs=io_configuration%active_moncs-1
    if (io_configuration%active_moncs==0) continue_poll_messages=.false.
    call check_thread_status(forthread_rwlock_unlock(monc_registration_lock))
  end subroutine handle_deregistration_command

  !> Retrieves the message from MONC off the data channel and throws this to a thread in the thread pool to actually process
  !! We do it this way to enforce ordering between the command (including the data set ID) and the raw data itself
  !! @param source Source PID of the MONC process
  !! @param data_set ID of the data set being communicated
  subroutine pull_back_data_message_and_handle(source, data_set)
    integer, intent(in) :: source, data_set

    integer :: specific_monc_data_type, specific_monc_buffer_size, recv_count, monc_location, matched_datadefn_index
    character, dimension(:), allocatable :: data_buffer

    call check_thread_status(forthread_rwlock_rdlock(monc_registration_lock))
    monc_location=get_monc_location(io_configuration, source)

    specific_monc_data_type=c_get_integer(io_configuration%registered_moncs(monc_location)%registered_monc_types, &
         conv_to_string(data_set))
    specific_monc_buffer_size=c_get_integer(io_configuration%registered_moncs(monc_location)%registered_monc_buffer_sizes, &
         conv_to_string(data_set))

    allocate(data_buffer(specific_monc_buffer_size))
    recv_count=data_receive(specific_monc_data_type, 1, source, dump_data=data_buffer, data_dump_id=data_set)


    ! This call is not handled by threading...should aid in ensuring that all time points are listed appropriately
    matched_datadefn_index=retrieve_data_definition(io_configuration, &
         io_configuration%registered_moncs(monc_location)%definition_names(data_set))
    if (matched_datadefn_index .gt. 0) then
      call inform_writer_federator_time_point(io_configuration, source, data_set, data_buffer)
    end if

    call check_thread_status(forthread_rwlock_unlock(monc_registration_lock))

    if (l_thoff) then
      call handle_data_message((/ source,  data_set /), data_buffer=data_buffer)
    else
      call threadpool_start_thread(handle_data_message, (/ source,  data_set /), data_buffer=data_buffer)
    end if

    deallocate(data_buffer)
  end subroutine pull_back_data_message_and_handle  

  !> Handles the command for data download from a specific process. This will allocate the receive buffer
  !! and then call to get the data. Once it has been received then the data is run against handling rules
  !! @param arguments, element 1 is the source & element 2 is the data_set 
  !! @param data_buffer The actual data from MONC read from the data channel
  subroutine handle_data_message(arguments, data_buffer)
    integer, dimension(:), intent(in) :: arguments
    character, dimension(:), allocatable, intent(inout), optional :: data_buffer

    integer :: monc_location, data_set, source, matched_datadefn_index

    source=arguments(1)
    data_set=arguments(2)

    call check_thread_status(forthread_rwlock_rdlock(monc_registration_lock))
    monc_location=get_monc_location(io_configuration, source)

    call check_thread_status(forthread_mutex_lock(io_configuration%registered_moncs(monc_location)%active_mutex))
    io_configuration%registered_moncs(monc_location)%active_threads=&
         io_configuration%registered_moncs(monc_location)%active_threads+1
    call check_thread_status(forthread_mutex_unlock(io_configuration%registered_moncs(monc_location)%active_mutex))
    
    matched_datadefn_index=retrieve_data_definition(io_configuration, &
         io_configuration%registered_moncs(monc_location)%definition_names(data_set))
    if (matched_datadefn_index .gt. 0) then
      call pass_fields_to_diagnostics_federator(io_configuration, source, data_set, data_buffer)
      call provide_monc_data_to_writer_federator(io_configuration, source, data_set, data_buffer)
      call check_writer_for_trigger(io_configuration, source, data_set, data_buffer)
    else
      call log_log(LOG_WARN, "IO server can not find matching data definition with name "&
           //io_configuration%registered_moncs(monc_location)%definition_names(data_set))
    end if    

    call check_thread_status(forthread_mutex_lock(io_configuration%registered_moncs(monc_location)%active_mutex))
    io_configuration%registered_moncs(monc_location)%active_threads=&
         io_configuration%registered_moncs(monc_location)%active_threads-1
    call check_thread_status(forthread_cond_signal(io_configuration%registered_moncs(monc_location)%deactivate_condition_variable))
    call check_thread_status(forthread_mutex_unlock(io_configuration%registered_moncs(monc_location)%active_mutex))
    call check_thread_status(forthread_rwlock_unlock(monc_registration_lock))
  end subroutine handle_data_message

  !> Handles registration from some MONC process. The source process sends some data description to this IO server which
  !! basically tells the IO server the size of the array datas (which might be different on different processes in the case
  !! of uneven decomposition.) Based upon this a communication (MPI) data type is constructed and the data size in bytes determined
  !! @param source The PID of the MONC process that is registering itself
  subroutine handle_monc_registration(arguments, data_buffer)
    integer, dimension(:), intent(in) :: arguments
    character, dimension(:), allocatable, intent(inout), optional :: data_buffer

    integer :: configuration_send_request(3), ierr, number_data_definitions, this_monc_index, source

    source=arguments(1)
    configuration_send_request=send_configuration_to_registree(source)
    number_data_definitions=io_configuration%number_of_data_definitions

    call check_thread_status(forthread_rwlock_wrlock(monc_registration_lock))

    io_configuration%number_of_moncs=io_configuration%number_of_moncs+1
    this_monc_index=io_configuration%number_of_moncs
    if (io_configuration%number_of_moncs .gt. size(io_configuration%registered_moncs)) then
      call log_log(LOG_ERROR, "You have a high ratio of computational cores to IO servers, the limit is currently 100")
      ! The extension of the MONC registration array is broken as the pointers involved in the map does not get copied across
      ! we could manually do this, but that is for another day! If you need to extend these limits either increase the constants
      ! or fix the extension, I don't think it will be too hard to fix the extension bit (copy the maps manually)
      call extend_registered_moncs_array(io_configuration)      
    end if

    io_configuration%active_moncs=io_configuration%active_moncs+1
    call check_thread_status(forthread_rwlock_unlock(monc_registration_lock))

    call c_put_integer(io_configuration%monc_to_index, conv_to_string(source), this_monc_index)

    call check_thread_status(forthread_mutex_init(io_configuration%registered_moncs(this_monc_index)%active_mutex, -1))
    call check_thread_status(forthread_cond_init(&
         io_configuration%registered_moncs(this_monc_index)%deactivate_condition_variable, -1))
    io_configuration%registered_moncs(this_monc_index)%active_threads=0
    io_configuration%registered_moncs(this_monc_index)%source_id=source

    allocate(io_configuration%registered_moncs(this_monc_index)%field_start_locations(number_data_definitions), &
         io_configuration%registered_moncs(this_monc_index)%field_end_locations(number_data_definitions), &
         io_configuration%registered_moncs(this_monc_index)%definition_names(number_data_definitions), &
         io_configuration%registered_moncs(this_monc_index)%dimensions(number_data_definitions))

    ! Wait for configuration to have been sent to registree
    call waitall_for_mpi_requests(configuration_send_request, 3)
    call init_data_definition(source, io_configuration%registered_moncs(this_monc_index))
  end subroutine handle_monc_registration

  !> Sends the data and field descriptions to the MONC process that just registered with the IO server
  !! @param source The MPI rank (MPI_COMM_WORLD) of the registree
  !! @returns The nonblocking send request handles which can be waited for completion later (overlap compute and communication)
  function send_configuration_to_registree(source)
    integer, intent(in) :: source
    integer :: send_configuration_to_registree(3)
    
    integer :: ierr, srequest(3)

    call lock_mpi()
    call mpi_isend(registree_definition_descriptions, size(registree_definition_descriptions), mpi_type_definition_description, &
         source, DATA_TAG, MPI_COMM_WORLD, srequest(1), ierr)
    call mpi_isend(registree_field_descriptions, size(registree_field_descriptions), mpi_type_field_description, &
         source, DATA_TAG, MPI_COMM_WORLD, srequest(2), ierr)
    call mpi_isend(sample_output_pairs, size(sample_output_pairs), MPI_INT, &
         source, DATA_TAG, MPI_COMM_WORLD, srequest(3), ierr)
    call unlock_mpi()

    send_configuration_to_registree=srequest    
  end function send_configuration_to_registree  

  !> Initialise the sizing of data definitions from a MONC process. The IO server determines, from configuration, the
  !! structure of each data definition but the size of the arrays depends upon the MONC process (due to uneven distribution
  !! of data etc...) This receives the sizing message and then builds the MPI datatype for each data definition that the IO 
  !! server will receive from that specific MONC process. The field sizings are for all fields in every data definition, and
  !! these are applied to each data definition which will simply ignore non matching fields
  !! @param source The source MONC PID
  !! @param monc_defn The corresponding MONC definition data structure
  subroutine init_data_definition(source, monc_defn)
    integer, intent(in) :: source
    type(io_configuration_registered_monc_type), intent(inout) :: monc_defn

    type(data_sizing_description_type) :: data_description(io_configuration%number_of_distinct_data_fields+4)
    integer :: created_mpi_type, data_size, recv_count, i
    type(data_sizing_description_type) :: field_description
    logical :: field_found
    
    recv_count=data_receive(mpi_type_data_sizing_description, io_configuration%number_of_distinct_data_fields+4, &
         source, description_data=data_description)

    call handle_monc_dimension_information(data_description, monc_defn)
     
    do i=1, io_configuration%number_of_data_definitions
      created_mpi_type=build_mpi_datatype(io_configuration%data_definitions(i), data_description, data_size, &
           monc_defn%field_start_locations(i), monc_defn%field_end_locations(i), monc_defn%dimensions(i))            

      call c_put_integer(monc_defn%registered_monc_types, conv_to_string(i), created_mpi_type)
      call c_put_integer(monc_defn%registered_monc_buffer_sizes, conv_to_string(i), data_size)

      monc_defn%definition_names(i)=io_configuration%data_definitions(i)%name
    end do
    if (.not. initialised_present_data) then
      initialised_present_data=.true.
      field_found=get_data_description_from_name(data_description, NUMBER_Q_INDICIES_KEY, field_description)
      call c_put_integer(io_configuration%dimension_sizing, "active_q_indicies", field_description%dim_sizes(1))
      call register_present_field_names_to_federators(data_description, recv_count)
    end if
    call get_monc_information_data(source)
  end subroutine init_data_definition

  !> Retrieves MONC information data, this is sent by MONC (and received) regardless, but only actioned if the data has not
  !! already been set
  !! @param source MONC source process
  subroutine get_monc_information_data(source)
    integer, intent(in) :: source

    character, dimension(:), allocatable :: buffer
    character(len=STRING_LENGTH) :: q_field_name, tracer_name, cd_field_name
    integer :: buffer_size, z_size, num_q_fields, num_tracers, n, current_point, recv_count
    type(data_sizing_description_type) :: field_description
    real(kind=DEFAULT_PRECISION) :: dreal
    logical :: field_found


    z_size=c_get_integer(io_configuration%dimension_sizing, "z")
    num_q_fields=c_get_integer(io_configuration%dimension_sizing, "qfields")
    num_tracers=c_get_integer(io_configuration%dimension_sizing, "tfields")

    buffer_size=(kind(dreal)*z_size)*2 + (STRING_LENGTH * num_q_fields) + STRING_LENGTH * num_tracers &
                 + 2*ncond*STRING_LENGTH + 2*ndiag*STRING_LENGTH 
    allocate(buffer(buffer_size))
    recv_count=data_receive(MPI_BYTE, buffer_size, source, buffer)
    if (.not. io_configuration%general_info_set) then
      call check_thread_status(forthread_mutex_lock(io_configuration%general_info_mutex))
      if (.not. io_configuration%general_info_set) then
        io_configuration%general_info_set=.true.
        allocate(io_configuration%zn_field(z_size))
        allocate(io_configuration%z_field(z_size))
        io_configuration%zn_field=transfer(buffer(1:kind(dreal)*z_size), io_configuration%zn_field)
        current_point=(kind(dreal)*z_size)
        if (num_q_fields .gt. 0) then
          do n=1, num_q_fields
            q_field_name=transfer(buffer(current_point+1:current_point+STRING_LENGTH), q_field_name)
            current_point=current_point+STRING_LENGTH
            call replace_character(q_field_name, " ", "_")
            call c_add_string(io_configuration%q_field_names, q_field_name)
          end do
        end if
    
        if (num_tracers .gt. 0) then
          do n=1, num_tracers
            tracer_name=transfer(buffer(current_point+1:current_point+STRING_LENGTH), tracer_name)
            current_point=current_point+STRING_LENGTH
            call c_add_string(io_configuration%tracer_names, tracer_name)
          end do
        end if

        io_configuration%z_field=transfer(buffer(current_point+1:current_point+kind(dreal)*z_size), &
                                          io_configuration%z_field)
        current_point=current_point+(kind(dreal)*z_size)

        do n=1,ncond
          cond_request(n)=transfer(buffer(current_point+1:current_point+STRING_LENGTH), cd_field_name)
          current_point=current_point+STRING_LENGTH
          cond_long(n)=transfer(buffer(current_point+1:current_point+STRING_LENGTH), cd_field_name)
          current_point=current_point+STRING_LENGTH
        end do

        do n=1,ndiag
          diag_request(n)=transfer(buffer(current_point+1:current_point+STRING_LENGTH), cd_field_name)
          current_point=current_point+STRING_LENGTH
          diag_long(n)=transfer(buffer(current_point+1:current_point+STRING_LENGTH), cd_field_name)
          current_point=current_point+STRING_LENGTH
        end do

      end if
      call provide_q_field_names_to_writer_federator(io_configuration%q_field_names)
      call provide_tracer_names_to_writer_federator(io_configuration%tracer_names)
      call check_thread_status(forthread_mutex_unlock(io_configuration%general_info_mutex))
    end if
    deallocate(buffer)
  end subroutine get_monc_information_data  

  !> Registers with the writer federator the set of fields (prognostic and diagnostic) that are available, this is based on
  !! the array/optional fields present from MONC and the non-optional scalars. This is quite an expensive operation, so only
  !! done once
  !! @param data_description Array of data descriptions from MONC
  !! @param recv_count Number of data descriptions
  subroutine register_present_field_names_to_federators(data_description, recv_count)
    type(data_sizing_description_type), dimension(:), intent(in) :: data_description
    integer, intent(in) :: recv_count

    type(hashset_type) :: present_field_names
    type(hashmap_type) :: diagnostics_field_names_and_roots
    integer :: i, j

    do i=1, recv_count
      call c_add_string(present_field_names, data_description(i)%field_name)
    end do
    do i=1, io_configuration%number_of_data_definitions
      do j=1, io_configuration%data_definitions(i)%number_of_data_fields
        if (io_configuration%data_definitions(i)%fields(j)%field_type == SCALAR_FIELD_TYPE .and. .not. &
             io_configuration%data_definitions(i)%fields(j)%optional) then
          call c_add_string(present_field_names, io_configuration%data_definitions(i)%fields(j)%name)
        end if        
      end do      
    end do
    call c_add_string(present_field_names, "time")
    call c_add_string(present_field_names, "timestep")
    call inform_writer_federator_fields_present(io_configuration, present_field_names)
    diagnostics_field_names_and_roots=determine_diagnostics_fields_available(present_field_names)
    call inform_writer_federator_fields_present(io_configuration, diag_field_names_and_roots=diagnostics_field_names_and_roots)
    call c_free(present_field_names)
    call c_free(diagnostics_field_names_and_roots)
  end subroutine register_present_field_names_to_federators  

  !> Handles the provided local MONC dimension and data layout information
  !! @param data_description The data descriptions sent over from MONC
  !! @param monc_defn The corresponding MONC definition data structure
  subroutine handle_monc_dimension_information(data_description, monc_defn)
    type(io_configuration_registered_monc_type), intent(inout) :: monc_defn
    type(data_sizing_description_type), dimension(:) :: data_description

    type(data_sizing_description_type) :: field_description
    integer :: i
    logical :: field_found

    field_found=get_data_description_from_name(data_description, LOCAL_SIZES_KEY, field_description)
    if (.not. field_found) call log_log(LOG_ERROR, "Malformed MONC registration, no local size information")
    do i=1,3
      monc_defn%local_dim_sizes(i)=field_description%dim_sizes(i)
    end do
    field_found=get_data_description_from_name(data_description, LOCAL_START_POINTS_KEY, field_description)
    if (.not. field_found) call log_log(LOG_ERROR, "Malformed MONC registration, no local start point information")
    do i=1,3
      monc_defn%local_dim_starts(i)=field_description%dim_sizes(i)
    end do
    field_found=get_data_description_from_name(data_description, LOCAL_END_POINTS_KEY, field_description)
    if (.not. field_found) call log_log(LOG_ERROR, "Malformed MONC registration, no local end point information")
    do i=1,3
      monc_defn%local_dim_ends(i)=field_description%dim_sizes(i)
    end do
  end subroutine handle_monc_dimension_information
end module io_server_mod
