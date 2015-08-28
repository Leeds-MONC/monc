!> Main core entry point to the rest of the model, this is called by the program main
module monc_mod
  use datadefn_mod, only : LONG_STRING_LENGTH, STRING_LENGTH
  use monc_component_mod, only : component_descriptor_type
  use collections_mod, only : list_type, map_type, c_size, c_key_at, c_value_at, c_get
  use conversions_mod, only : conv_to_string, conv_to_real
  use state_mod, only : model_state_type
  use registry_mod, only : get_all_registered_components, get_component_info, get_missing_registry_components, &
       register_component, execute_initialisation_callbacks, &
       execute_finalisation_callbacks, init_registry, order_all_callbacks, display_callbacks_in_order_at_each_stage
  use timestepper_mod, only : init_timestepper, timestep, finalise_timestepper
  use logging_mod, only : LOG_INFO, LOG_WARN, LOG_ERROR, LOG_DEBUG, log_log, log_get_logging_level, log_set_logging_level, &
       log_master_log, initialise_logging
  use optionsdatabase_mod, only : load_command_line_into_options_database, options_get_integer, options_has_key, &
       options_get_string, options_get_logical
  use configuration_file_parser_mod, only : parse_configuration_file
  use configuration_checkpoint_netcdf_parser_mod, only : parse_configuration_checkpoint_netcdf
  use science_constants_mod, only : initialise_science_constants
  use mpi, only : MPI_COMM_WORLD, MPI_THREAD_MULTIPLE
  use datadefn_mod, only : DEFAULT_PRECISION, init_data_defn
  implicit none

#ifndef TEST_MODE
  private
#endif
  abstract interface
     !> IO server entry procedure which may be passed to the core entry point (if IO server is enabled)
     !! @param io_communicator_arg The IO communicator
     !! @param io_xml_configuration The IO server textual configuration
     subroutine io_server_run_procedure(options_database, io_communicator_arg, io_xml_configuration, provided_threading, &
          total_global_processes)
       import map_type
       type(map_type), intent(inout) :: options_database
       integer, intent(in) :: io_communicator_arg, provided_threading, total_global_processes
       character, dimension(:), allocatable, intent(inout) :: io_xml_configuration       
     end subroutine io_server_run_procedure
  end interface

  !< For reading the IO XML configuration, these are string length constants which can be increased if required
  integer, parameter :: FILE_STR_STRIDE=10000, FILE_LINE_LEN=2000

  public monc_core_bootstrap

contains

  !> Main core entry point to bootstrap running the model
  !!
  !! Reads in command line arguments, sets up the model state and registers components. Then runs through
  !! and calls the execution of each model stage.
  !! @param componentDescriptions Descriptions of existing components which should be registered
  !! @param io_server_run Optional IO server entry procedure
  subroutine monc_core_bootstrap(component_descriptions, io_server_run)
    type(list_type), intent(inout) :: component_descriptions
    procedure(io_server_run_procedure) :: io_server_run

    type(model_state_type) :: state
    integer :: ierr, myrank, size, io_server_placement_period, provided_threading
    logical :: i_am_monc_process, enable_io_server
    character(len=LONG_STRING_LENGTH) :: io_server_config_file
    character, dimension(:), allocatable :: io_server_configuration_contents

    call load_model_configuration(state%options_database)

    enable_io_server=determine_if_io_server_enabled(state%options_database)

    if (enable_io_server) then
      call mpi_init_thread(MPI_THREAD_MULTIPLE, provided_threading, ierr)
    else
      call mpi_init(ierr)
    end if
    
    call init_data_defn()
    ! Set up the logging with comm world PIDs initially for logging from the configuration parsing
    call mpi_comm_rank(MPI_COMM_WORLD, myrank, ierr)
    call initialise_logging(myrank)
    
    call log_set_logging_level(options_get_integer(state%options_database, "logging_modlevel"))

    if (enable_io_server) then
      call mpi_comm_size(MPI_COMM_WORLD, size, ierr)
      if (size==1) call log_log(LOG_ERROR, &
           "Run with 1 process, With IO server enabled then the minimum process size is 2 (1 for IO, 1 for MONC)")
      call get_io_configuration(state%options_database, io_server_config_file, io_server_placement_period)
      call split_communicator_into_monc_and_io(io_server_placement_period, state%parallel%monc_communicator, &
           state%parallel%io_communicator, i_am_monc_process, state%parallel%corresponding_io_server_process)
      if (.not. i_am_monc_process) then        
        io_server_configuration_contents=get_io_xml(io_server_config_file)
        call io_server_run(state%options_database, state%parallel%io_communicator, &
             io_server_configuration_contents, provided_threading, size)       
      else  
        call monc_run(component_descriptions, state)
      end if
    else
      state%parallel%monc_communicator=MPI_COMM_WORLD
      call monc_run(component_descriptions, state)
    end if

    call mpi_finalize(ierr)
  end subroutine monc_core_bootstrap

  !> Determines whether the IO server should be enabled or not
  !! @param options_database The options database
  !! @returns Whether to enable the IO server or not
  logical function determine_if_io_server_enabled(options_database)
    type(map_type), intent(inout) :: options_database

    determine_if_io_server_enabled=options_get_logical(options_database, "enable_io_server")
    if (determine_if_io_server_enabled) then
      determine_if_io_server_enabled=options_get_logical(options_database, "iobridge_enabled")      
    end if    
  end function determine_if_io_server_enabled  

  !> Loads the configuration into the options database, either from a file or checkpoint
  !! @param options_database The options database
  subroutine load_model_configuration(options_database)
    type(map_type), intent(inout) :: options_database

    call load_command_line_into_options_database(options_database)
    if (options_has_key(options_database, "config")) then
      call parse_configuration_file(options_database, options_get_string(options_database, "config"))
    else if (options_has_key(options_database, "checkpoint")) then
      call parse_configuration_checkpoint_netcdf(options_database, &
           options_get_string(options_database, "checkpoint"), MPI_COMM_WORLD)
    else
      call log_master_log(LOG_ERROR, "You must either provide a configuration file or checkpoint to restart from")
      call mpi_barrier(MPI_COMM_WORLD) ! All other processes barrier here to ensure 0 displays the message before quit
      stop
    end if
    ! Reload command line arguments to override any stuff in the configuration files
    call load_command_line_into_options_database(options_database)
  end subroutine load_model_configuration 

  !> Called by MONC processes to run the MONC model
  !! @param componentDescriptions Descriptions of existing components which should be registered
  !! @param state The current model state
  subroutine monc_run(component_descriptions, state)
    type(list_type), intent(inout) :: component_descriptions
    type(model_state_type), intent(inout) :: state

    integer :: ierr
    real(kind=DEFAULT_PRECISION) :: start_time, end_time, timestepping_time, modeldump_time

    call cpu_time(start_time)
    call mpi_comm_rank(state%parallel%monc_communicator, state%parallel%my_rank, ierr)
    call mpi_comm_size(state%parallel%monc_communicator, state%parallel%processes, ierr)

    call initialise_logging(state%parallel%my_rank)
    
    call log_master_log(LOG_INFO,"MONC running with "//trim(conv_to_string(state%parallel%processes))//" processes")

#ifdef DEBUG_MODE
    call log_master_log(LOG_WARN,"MONC compiled with debug options, you probably want to recompile without for production runs")
#endif    

    call init_registry(state%options_database) ! Initialise the registry

    call fill_registry_with_components(state%options_database, component_descriptions)
    call initialise_science_constants(state)
    call order_all_callbacks()
    ! If the option has been provided then display the registered component information
    if (is_present_and_true(state%options_database, "registered") .and. state%parallel%my_rank==0) &
         call display_registed_components()
    if (is_present_and_true(state%options_database, "warnmissing") .and. state%parallel%my_rank==0) &
         call issue_missing_registry_component_warnings()
    if (is_present_and_true(state%options_database, "showcallbacks") .and. state%parallel%my_rank==0) &
         call display_callbacks_in_order_at_each_stage()

    if (.not. is_present_and_true(state%options_database, "norun")) then
      ! Unless configured otherwise then run through the different stages of execution
      call perform_model_steps(state, timestepping_time, modeldump_time)
    end if
    call mpi_barrier(state%parallel%monc_communicator, ierr)
    call cpu_time(end_time)
    if (state%parallel%my_rank==0) then
      call log_log(LOG_INFO, "Entire MONC run completed in "//trim(conv_to_string(int((end_time-start_time) * 1000)))//&
           "ms (timestepping="//trim(conv_to_string(int(timestepping_time * 1000)))//"ms, modeldump="//&
           trim(conv_to_string(int(modeldump_time * 1000)))//"ms, misc="//trim(conv_to_string((&
           int((end_time-start_time) * 1000)) - (int(timestepping_time * 1000) + int(modeldump_time * 1000))))//"ms)")
    end if    
  end subroutine monc_run

  !> Will run through the actual model stages and call the appropriate callbacks at each stage
  !! @param state The model state
  !! @timestepping_time The time spent in doing actual timestepping (computation)
  !! @modeldump_time The time spent in doing the model dump
  subroutine perform_model_steps(state, timestepping_time, modeldump_time)
    type(model_state_type), intent(inout) :: state
    real(kind=DEFAULT_PRECISION), intent(out) :: timestepping_time, modeldump_time

    integer :: logging_mod_level
    real(kind=DEFAULT_PRECISION) :: start_time, end_time, start_iteration_time

    timestepping_time=0.0_DEFAULT_PRECISION
    modeldump_time=0.0_DEFAULT_PRECISION

    call init_timestepper()

    logging_mod_level = log_get_logging_level()
    call execute_initialisation_callbacks(state)
    state%continue_timestep=.true.
    call cpu_time(start_time)
    do while (state%continue_timestep)
      if (state%update_dtm) state%dtm=state%dtm_new
      ! The start of a timestep
      if (logging_mod_level .ge. LOG_DEBUG) call cpu_time(start_iteration_time)
      call timestep(state) ! Call out to the timestepper to do the actual timestepping per component
      if (logging_mod_level .ge. LOG_DEBUG .and. state%parallel%my_rank==0) &
           call display_timestep_information(state%timestep, start_iteration_time)
      state%timestep = state%timestep+1
      state%time = state%time + state%dtm
    end do
    call cpu_time(end_time)
    state%timestep_runtime=end_time-start_time
    timestepping_time=timestepping_time+state%timestep_runtime
    call execute_finalisation_callbacks(state)

    call finalise_timestepper()
  end subroutine perform_model_steps

  !> Provides timestepping information about the current step and performance
  !! @param timestep The current timestep which has been completed
  !! @param startTime The F95 CPU time that the current timestep was started at
  subroutine display_timestep_information(timestep, start_time)
    integer, intent(in) :: timestep
    real(kind=DEFAULT_PRECISION), intent(in) :: start_time

    real(kind=DEFAULT_PRECISION) :: end_time

    call cpu_time(end_time)
    call log_log(LOG_DEBUG, "Timestep "//trim(conv_to_string(timestep))//" completed in "//&
         trim(conv_to_string(int((end_time-start_time) * 1000)))//"ms")
  end subroutine display_timestep_information  

  !> Registers each supplied component description
  subroutine fill_registry_with_components(options_database, component_descriptions)
    type(map_type), intent(inout) :: options_database
    type(list_type), intent(inout) :: component_descriptions

    integer :: i, number_of_components
    class(*), pointer :: raw_data

    number_of_components = c_size(component_descriptions)
    do i=1,number_of_components
      raw_data=>c_get(component_descriptions, i)
      select type(raw_data)
        type is (component_descriptor_type)
          call register_component(options_database, raw_data)
        class default
          call log_log(LOG_WARN, "Can not register component due to corrupted data")
      end select
    end do
  end subroutine fill_registry_with_components

  !> Determines whether an option is present in the database and true. This combines the key check and getting
  !! the value. Just calling to get the value directly will error if it does not exist, we don't nescesarily
  !! want for checking optional command line flags
  !! @param options_database The options database
  !! @param key The key to test for
  logical function is_present_and_true(options_database, key)
    type(map_type), intent(inout) :: options_database
    character(len=*), intent(in) :: key
    
    if (options_has_key(options_database, key)) then
      is_present_and_true=options_get_logical(options_database, key)
      return
    end if
    is_present_and_true=.false.
  end function is_present_and_true 

  !> Issues warning log_log messages for each component which has not been registered
  subroutine issue_missing_registry_component_warnings()
    type(list_type) :: missing_components
    integer :: i, missing_size
    class(*), pointer :: raw_data

    missing_components = get_missing_registry_components()
    missing_size = c_size(missing_components)
    do i=1,missing_size
      raw_data=>c_get(missing_components, i)
      call log_log(LOG_WARN, "Missing registry component: "//conv_to_string(raw_data, .false., STRING_LENGTH))
    end do
  end subroutine issue_missing_registry_component_warnings

  !> Displays the registered components and their version numbers
  subroutine display_registed_components()
    type(map_type) :: registered_components
    integer :: i
    character(len=STRING_LENGTH) :: component_name
    type(component_descriptor_type), pointer :: detailed_component_info
    real :: component_version
    class(*), pointer :: raw_data

    registered_components = get_all_registered_components()
    call log_log(LOG_INFO, "Registered components: "//conv_to_string(c_size(registered_components)))
    do i=1,c_size(registered_components)
      component_name = c_key_at(registered_components, i)
      raw_data=>c_value_at(registered_components, i)
      component_version = conv_to_real(raw_data, .false.)
      detailed_component_info => get_component_info(component_name)
      call log_log(LOG_INFO, trim(component_name)//" "//trim(conv_to_string(component_version)))
    end do
  end subroutine display_registed_components

  !> Splits the MPI_COMM_WORLD communicator into MONC and IO separate communicators. The size of each depends
  !! on the stride supplied.
  !! @param io_stride The absolute process id stride for IO processes
  !! @param monc_communicator The communicator associated with MONC processes
  !! @param io_communicator The communicator associated with IO processes
  subroutine split_communicator_into_monc_and_io(io_stride, monc_communicator, io_communicator, &
       am_i_monc_process, corresponding_io_server_process)
    integer, intent(in) :: io_stride
    integer, intent(out) :: monc_communicator, io_communicator, corresponding_io_server_process
    logical, intent(out) :: am_i_monc_process

    integer, dimension(:), allocatable :: members_monc_group, members_io_group
    integer :: total_ranks, monc_group, io_group, io_processes, monc_processes, i, io_index, &
         monc_index, my_rank, ierr, global_group

    call mpi_comm_size(MPI_COMM_WORLD, total_ranks, ierr)
    call mpi_comm_rank(MPI_COMM_WORLD, my_rank, ierr)

    io_processes=get_number_io_processes(total_ranks, io_stride)
    monc_processes=total_ranks-io_processes
    allocate(members_io_group(io_processes), members_monc_group(monc_processes))
    io_index=1
    monc_index=1
    am_i_monc_process=.true.

    do i=0, total_ranks-1
      if (mod(i, io_stride) == 0 .and. i .lt. total_ranks-1) then
        members_io_group(io_index)=i
        io_index=io_index+1
        if (my_rank == i) am_i_monc_process=.false.
        if (my_rank .gt. i .and. my_rank .lt. i+io_stride) corresponding_io_server_process=i
      else
        members_monc_group(monc_index)=i
        monc_index=monc_index+1
      end if
    end do

    call mpi_comm_group(MPI_COMM_WORLD, global_group, ierr)
    call mpi_group_incl(global_group, monc_processes, members_monc_group, monc_group, ierr)
    call mpi_group_incl(global_group, io_processes, members_io_group, io_group, ierr)
    call mpi_comm_create(MPI_COMM_WORLD, monc_group, monc_communicator, ierr)
    call mpi_comm_create(MPI_COMM_WORLD, io_group, io_communicator, ierr)
    deallocate(members_io_group, members_monc_group)
  end subroutine split_communicator_into_monc_and_io

  !> Based upon the total number of processes and the IO process id stride determines the number of 
  !! processes that will be used for the IO server. The MONC processes is total processes - io processes
  !! @param total_ranks Total number of processes in use
  !! @param io_stride The absolute process id stride for IO processes
  !! @returns The number of processes used for running the IO server
  integer function get_number_io_processes(total_ranks, io_stride)
    integer, intent(in) :: total_ranks, io_stride

    get_number_io_processes=total_ranks/io_stride
    if (get_number_io_processes * io_stride .lt. total_ranks-1) get_number_io_processes=get_number_io_processes+1
  end function get_number_io_processes

  !> Reads the IO server configuration and populates the required variables of the configuration
  !! file name and the placement period
  !! @param options_database The options database
  subroutine get_io_configuration(options_database, ioserver_configuration_file, ioserver_placement_period)
    type(map_type), intent(inout) :: options_database
    character(len=LONG_STRING_LENGTH), intent(out) :: ioserver_configuration_file
    integer, intent(out) :: ioserver_placement_period
   
    integer :: myrank, ierr

    ioserver_configuration_file=options_get_string(options_database, "ioserver_configuration_file")
    ioserver_placement_period=options_get_integer(options_database, "ioserver_placement_period")

    if (ioserver_placement_period == -1 .or. ioserver_configuration_file == "") then
      call mpi_comm_rank(MPI_COMM_WORLD, myrank, ierr)
      if (myrank == 0) call log_log(LOG_ERROR, "To run an IO server you must provide the placement period and configuration file")
      call mpi_barrier(MPI_COMM_WORLD) ! All other processes barrier here to ensure 0 displays the message before quit
      stop
    end if
  end subroutine get_io_configuration

  !> Reads in textual data from a file and returns this, used to read the IO server XML configuration file. Returned is a
  !! character array of exactly the correct size filled with all the configuration
  !! @param filename Name of the file to read
  !! @returns The contents of the XML file
  recursive function get_io_xml(filename, funit_num) result(io_xml)
    character(len=*), intent(in) :: filename
    integer, intent(in), optional :: funit_num
    character, dimension(:), allocatable :: io_xml, temp_io_xml

    character(len=FILE_LINE_LEN) :: temp_line, adjusted_io_line
    character(len=FILE_STR_STRIDE) :: reading_buffer
    integer :: ierr, first_quote, last_quote, chosen_unit

    if (present(funit_num)) then
      chosen_unit=funit_num
    else
      chosen_unit=2
    end if

    reading_buffer=""
    open (unit=chosen_unit, file=filename, status='OLD', iostat=ierr)
    if (ierr .ne. 0) call log_log(LOG_ERROR, "Error opening file '"//trim(filename)//"'")
    do while (ierr == 0)
      read(chosen_unit,"(A)",iostat=ierr) temp_line
      adjusted_io_line=adjustl(temp_line)
      if (ierr == 0 .and. adjusted_io_line(1:1) .ne. "!" .and. adjusted_io_line(1:2) .ne. "//") then
        if (index(temp_line, "#include") .ne. 0) then
          first_quote=index(temp_line, """")
          last_quote=index(temp_line, """", back=.true.)
          if (first_quote .ne. 0 .and. last_quote .ne. 0) then
            call add_in_specific_line(io_xml, reading_buffer)
            temp_io_xml=get_io_xml(temp_line(first_quote+1:last_quote-1), chosen_unit+1)
            call combine_xml_arrays(io_xml, temp_io_xml)
            deallocate(temp_io_xml)
            reading_buffer=new_line("A")
          else
            call log_log(LOG_ERROR, "Malformed IO XML, include directives must have filename in quotes")
          end if          
        else
          if (len_trim(reading_buffer) + len_trim(temp_line) .ge. FILE_STR_STRIDE) then
            call add_in_specific_line(io_xml, reading_buffer)
            reading_buffer=""
          end if
          reading_buffer=trim(reading_buffer)//trim(temp_line)//new_line("A")
        end if
      end if
    end do
    if (len_trim(reading_buffer) .gt. 0) call add_in_specific_line(io_xml, reading_buffer)
    close(chosen_unit)
  end function get_io_xml

  !> Adds a specific line into the io xml. The IO XML is always exactly the correct size, so here is either allocated or
  !! resized to match what the read buffer requires
  !! @param io_xml The IO XML which holds all the configuration and is exactly the correct size
  !! @param reading_buffer A buffer which will be copied into the resized/allocated IO XML
  subroutine add_in_specific_line(io_xml, reading_buffer)
    character, dimension(:), allocatable, intent(inout) :: io_xml
    character(len=*), intent(in) :: reading_buffer

    character, dimension(:), allocatable :: temp_io_xml
    integer :: i

    if (.not. allocated(io_xml)) then
      allocate(io_xml(len_trim(reading_buffer)))
      do i=1, len_trim(reading_buffer)
        io_xml(i)=reading_buffer(i:i)
      end do
    else
      allocate(temp_io_xml(size(io_xml)+len_trim(reading_buffer)))
      temp_io_xml(:size(io_xml)) = io_xml
      do i=1, len_trim(reading_buffer)
        temp_io_xml(size(io_xml)+i) = reading_buffer(i:i)
      end do
      call move_alloc(from=temp_io_xml,to=io_xml)
    end if
  end subroutine add_in_specific_line

  !> Combines two IO XML arrays together (for instance one returned from a recursive include)
  !! @param io_xml The IO XML is a source and target, this is allocated or resized to hold its contents + other arrays contents
  !! @param other_xml_array The other XML array which will be copied into the IO XML
  subroutine combine_xml_arrays(io_xml, other_xml_array)
    character, dimension(:), allocatable, intent(inout) :: io_xml, other_xml_array

    character, dimension(:), allocatable :: temp_io_xml

    if (.not. allocated(other_xml_array)) return

    if (.not. allocated(io_xml)) then
      allocate(io_xml(size(other_xml_array)), source=other_xml_array)
    else
      allocate(temp_io_xml(size(io_xml)+size(other_xml_array)))
      temp_io_xml(:size(io_xml)) = io_xml
      temp_io_xml(size(io_xml)+1:) = other_xml_array
      call move_alloc(from=temp_io_xml,to=io_xml)
    end if
  end subroutine combine_xml_arrays
end module monc_mod
