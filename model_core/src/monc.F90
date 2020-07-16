!> Main core entry point to the rest of the model, this is called by the program main
module monc_mod
  use datadefn_mod, only : LONG_STRING_LENGTH, STRING_LENGTH
  use monc_component_mod, only : component_descriptor_type
  use collections_mod, only : list_type, hashmap_type, map_type, iterator_type, mapentry_type, c_size, c_next_generic, &
       c_get_real, c_has_next, c_get_iterator, c_next_mapentry
  use conversions_mod, only : conv_to_string, conv_to_real, conv_is_logical, conv_to_logical
  use state_mod, only : model_state_type
  use registry_mod, only : get_all_registered_components, get_component_info, register_component, &
       execute_initialisation_callbacks, execute_finalisation_callbacks, init_registry, order_all_callbacks, &
       display_callbacks_in_order_at_each_stage
  use timestepper_mod, only : init_timestepper, timestep, finalise_timestepper
  use logging_mod, only : LOG_INFO, LOG_WARN, LOG_ERROR, LOG_DEBUG, log_log, log_get_logging_level, log_set_logging_level, &
       log_master_log, initialise_logging, log_master_newline
  use optionsdatabase_mod, only : load_command_line_into_options_database, options_get_integer, options_has_key, &
       options_get_string, options_get_logical, options_add, options_remove_key
  use configuration_file_parser_mod, only : parse_configuration_file
  use configuration_checkpoint_netcdf_parser_mod, only : parse_configuration_checkpoint_netcdf
  use science_constants_mod, only : initialise_science_constants
  use mpi, only : MPI_COMM_WORLD, MPI_THREAD_MULTIPLE, MPI_THREAD_SERIALIZED, MPI_THREAD_SINGLE, MPI_THREAD_FUNNELED, mpi_wtime
  use datadefn_mod, only : DEFAULT_PRECISION, init_data_defn
#ifndef TEST_MODE
  use netcdf, only : nf90_nowrite, nf90_open, nf90_inq_varid, nf90_get_var, nf90_close
#else
  use dummy_netcdf_mod, only : nf90_nowrite, nf90_open, nf90_inq_varid, nf90_get_var, nf90_close
#endif
  use checkpointer_common_mod, only : TIME_KEY, check_status

  implicit none

#ifndef TEST_MODE
  private
#endif
  abstract interface
     !> IO server entry procedure which may be passed to the core entry point (if IO server is enabled)
     !! @param io_communicator_arg The IO communicator
     !! @param io_xml_configuration The IO server textual configuration
     subroutine io_server_run_procedure(options_database, io_communicator_arg, provided_threading, &
          total_global_processes, continuation_run, reconfig_initial_time, io_configuration_file,  &
          my_global_rank)
       import hashmap_type, LONG_STRING_LENGTH, DEFAULT_PRECISION
       type(hashmap_type), intent(inout) :: options_database
       integer, intent(in) :: io_communicator_arg, provided_threading, total_global_processes, my_global_rank
       logical, intent(in) :: continuation_run
       real(kind=DEFAULT_PRECISION), intent(in) :: reconfig_initial_time
       character(len=LONG_STRING_LENGTH), intent(in) :: io_configuration_file
     end subroutine io_server_run_procedure
  end interface

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
    integer :: ierr, myrank, size, io_server_placement_period, provided_threading, selected_threading_mode
    logical :: i_am_monc_process
    logical :: io_continuation
    real(kind=DEFAULT_PRECISION) :: reconfig_initial_time
    character(len=LONG_STRING_LENGTH) :: io_server_config_file

    ! Initialise MPI
    selected_threading_mode=get_mpi_threading_mode()
    call mpi_init_thread(selected_threading_mode, provided_threading, ierr)

    call init_data_defn()

    ! Set up the logging with comm world PIDs initially for logging from the configuration parsing
    call mpi_comm_rank(MPI_COMM_WORLD, myrank, ierr)
    state%parallel%my_global_rank = myrank
    call initialise_logging(myrank)
    call log_master_log(LOG_INFO, "Starting MONC...")
    call log_master_newline()

    if (selected_threading_mode .gt. provided_threading) then
      call log_master_log(LOG_ERROR, "You have selected to thread at level '"//&
           trim(mpi_threading_level_to_string(selected_threading_mode))//&
           "' but the maximum level your MPI implementation can provide is '"//&
           trim(mpi_threading_level_to_string(provided_threading))//"'")
    end if

    ! Load model configuration
    reconfig_initial_time = 0.0_DEFAULT_PRECISION ! set locally, should precede load_model_configuration
    call load_model_configuration(state, state%options_database, io_continuation, reconfig_initial_time)
    call log_set_logging_level(options_get_integer(state%options_database, "logging"))
    call perform_options_compatibility_checks(state%options_database)

    ! Check on io_server settings and start MONC
    state%io_server_enabled=determine_if_io_server_enabled(state%options_database)
    if (state%io_server_enabled) then
      call mpi_comm_size(MPI_COMM_WORLD, size, ierr)
      if (size==1) call log_log(LOG_ERROR, &
           "Run with 1 process, With IO server enabled then the minimum process size is 2 (1 for IO, 1 for MONC)")
      call get_io_configuration(state%options_database, io_server_config_file, io_server_placement_period)
      call split_communicator_into_monc_and_io(io_server_placement_period, state%parallel%monc_communicator, &
           state%parallel%io_communicator, i_am_monc_process, state%parallel%corresponding_io_server_process)
      if (.not. i_am_monc_process) then        
        call io_server_run(state%options_database, state%parallel%io_communicator, provided_threading, &
             size, io_continuation, reconfig_initial_time, io_server_config_file, myrank)
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
    type(hashmap_type), intent(inout) :: options_database

    determine_if_io_server_enabled=options_get_logical(options_database, "enable_io_server")
    if (determine_if_io_server_enabled) then
      determine_if_io_server_enabled=options_get_logical(options_database, "iobridge_enabled")      
    end if    
  end function determine_if_io_server_enabled  

  !> Loads the configuration into the options database, either from a file or checkpoint
  !! @param options_database The options database
  !! @param io_continuation  Whether the io_server should be reinitialized from the checkpoint.
  subroutine load_model_configuration(state, options_database, io_continuation, reconfig_initial_time)
    type(model_state_type), intent(inout) :: state
    type(hashmap_type), intent(inout) :: options_database
    logical, intent(out) :: io_continuation
    real(kind=DEFAULT_PRECISION), intent(inout) :: reconfig_initial_time

    call load_command_line_into_options_database(options_database)

    ! Cold start from mcf config file.
    if (options_has_key(options_database, "config")) then
      state%continuation_run=.false.
      io_continuation=.false.
      call log_master_log(LOG_INFO, "This cycle is a cold start using config: '"//&
                          trim(options_get_string(options_database, "config"))//"'")
      call log_master_newline()
      call parse_configuration_file(options_database, options_get_string(options_database, "config"))

    ! Reconfiguration reads configuration from mcf and data from netcdf checkpoint.  
    ! Calling it reconfig allows this startup option.
    ! This is a continuation run for MONCs, but not for the IOserver.
    else if (options_has_key(options_database, "reconfig") .and. &
             options_has_key(options_database, "checkpoint")) then
      state%reconfig_run=.true.    ! this is specific to the initial cycle of this kind of run
      state%continuation_run=.true.
      io_continuation=.false.
      call log_master_log(LOG_INFO, "This cycle is a reconfigured start using config: '"//&
             trim(options_get_string(options_database, "reconfig"))//&
             "' from checkpoint: '"//trim(options_get_string(options_database, "checkpoint"))//"'")
                  
      if (options_get_logical(options_database, "retain_model_time")) then
        call extract_time_from_checkpoint_file(options_get_string(options_database, "checkpoint"),&
                                               reconfig_initial_time)
        state%retain_model_time = .true.
      end if

      call log_master_log(LOG_INFO, "Reconfiguration starting from time: "//trim(conv_to_string(reconfig_initial_time)))
      call log_master_newline()

      call parse_configuration_file(options_database, &
               options_get_string(options_database, "reconfig"))

    ! Continuation
    else if (options_has_key(options_database, "checkpoint")) then
      state%continuation_run=.true.
      io_continuation=.true.
      call log_master_log(LOG_INFO, "This cycle is a continuation from checkpoint: '"//&
                          trim(options_get_string(options_database, "checkpoint"))//"'")
      call parse_configuration_checkpoint_netcdf(options_database, &
               options_get_string(options_database, "checkpoint"), MPI_COMM_WORLD)
      call log_master_log(LOG_INFO, "Continuation uses config: '"//&
                          trim(options_get_string(options_database, "config"))//"'")
      call log_master_newline()

    ! Error - appropriate start conditions not met.
    else
      call log_master_log(LOG_ERROR, "You must provide a configuration file for a cold start,"//&
       " a checkpoint to restart from, or a reconfig file and a checkpoint to reconfigure from.")
      call mpi_barrier(MPI_COMM_WORLD) ! All other processes barrier here to ensure 0 displays the message before quit
      stop
    end if

    ! Reload command line arguments to override any stuff in the configuration files
    call load_command_line_into_options_database(options_database, .true.)

    ! In the case of reconfig, we won't want it to do this again on a later continuation cycle, so we'll remove
    ! the reconfig key, after recording the source file as config.
    if (options_has_key(options_database, "reconfig")) then
      call options_add(options_database, "config", options_get_string(options_database, "reconfig"))
      call options_remove_key(options_database, "reconfig")
      call options_remove_key(options_database, "retain_model_time")
    end if

  end subroutine load_model_configuration 

  !> Performs options_database compatibility checks.
  !! @param options_database The options database
  subroutine perform_options_compatibility_checks(options_database)
    type(hashmap_type), intent(inout) :: options_database

    !> If conditional diagnostics are operating, both components should be enabled.
    if (is_present_and_true(options_database, "conditional_diagnostics_column_enabled") .NEQV. &
        is_present_and_true(options_database, "conditional_diagnostics_whole_enabled")  ) then

      if ( .not. is_present_and_true(options_database, "conditional_diagnostics_column_enabled")) &
        call options_add(options_database, "conditional_diagnostics_column_enabled", .true.)
      if ( .not. is_present_and_true(options_database, "conditional_diagnostics_whole_enabled")) &
        call options_add(options_database, "conditional_diagnostics_whole_enabled", .true.)

      call log_master_log(LOG_INFO, "Only one conditional_diagnostics component is enabled, but both are required to function.")
      call log_master_log(LOG_INFO, "We assume you would like conditional_diagnostics enabled so have enabled the other, too.")
    end if

    !> In order to use time_basis or force_output_on_interval as intended, the cfltest and iobridge
    !  components need to be enabled, and the io_server must be enabled.
    if (is_present_and_true(options_database, "time_basis") .or.   &
        is_present_and_true(options_database, "force_output_on_interval") ) then
      if ( .not. (is_present_and_true(options_database, "iobridge_enabled") .and.  &
                  is_present_and_true(options_database, "cfltest_enabled")  .and.  &
                  is_present_and_true(options_database, "enable_io_server")          )) &
        call log_master_log(LOG_ERROR, "In order to use time_basis or force_output_on_interval"//&
         " as intended, the cfltest and iobridge components need to be enabled, and the"//&
         " io_server must be enabled (to permit full function of iobridge).")
    end if 

  end subroutine perform_options_compatibility_checks

  !> Called by MONC processes to run the MONC model
  !! @param componentDescriptions Descriptions of existing components which should be registered
  !! @param state The current model state
  subroutine monc_run(component_descriptions, state)
    type(list_type), intent(inout) :: component_descriptions
    type(model_state_type), intent(inout) :: state

    integer :: ierr, total_size
    double precision :: end_time, timestepping_time, modeldump_time

    state%model_start_wtime=mpi_wtime()
    call mpi_comm_rank(state%parallel%monc_communicator, state%parallel%my_rank, ierr)
    call mpi_comm_size(state%parallel%monc_communicator, state%parallel%processes, ierr)
    call mpi_comm_size(MPI_COMM_WORLD, total_size, ierr)

    call initialise_logging(state%parallel%my_rank)
    
    call log_master_log(LOG_INFO,"MONC running with "//trim(conv_to_string(state%parallel%processes))//" processes, "// &
         trim(conv_to_string(total_size-state%parallel%processes))// " IO server(s)")

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
    if (is_present_and_true(state%options_database, "showcallbacks") .and. state%parallel%my_rank==0) &
         call display_callbacks_in_order_at_each_stage()

    if (.not. is_present_and_true(state%options_database, "norun")) then
      ! Unless configured otherwise then run through the different stages of execution
      call perform_model_steps(state, timestepping_time, modeldump_time)
    end if
    call mpi_barrier(state%parallel%monc_communicator, ierr)
    end_time=mpi_wtime()
    if (state%parallel%my_rank==0) then
      call log_log(LOG_INFO, "Entire MONC run completed in "//trim(conv_to_string(int((end_time-state%model_start_wtime)*1000)))//&
           "ms (timestepping="//trim(conv_to_string(int(timestepping_time * 1000)))//"ms, modeldump="//&
           trim(conv_to_string(int(modeldump_time * 1000)))//"ms, misc="//trim(conv_to_string((&
           int((end_time-state%model_start_wtime) * 1000)) - (int(timestepping_time * 1000) + int(modeldump_time * 1000))))//"ms)")
    end if    
  end subroutine monc_run

  !> Will run through the actual model stages and call the appropriate callbacks at each stage
  !! @param state The model state
  !! @timestepping_time The time spent in doing actual timestepping (computation)
  !! @modeldump_time The time spent in doing the model dump
  subroutine perform_model_steps(state, timestepping_time, modeldump_time)
    type(model_state_type), intent(inout) :: state
    double precision, intent(out) :: timestepping_time, modeldump_time

    integer :: logging_mod_level
    double precision :: start_time, end_time, start_iteration_time

    timestepping_time=0.0_DEFAULT_PRECISION
    modeldump_time=0.0_DEFAULT_PRECISION

    call init_timestepper()

    logging_mod_level = log_get_logging_level()
    call execute_initialisation_callbacks(state)
    state%continue_timestep=.true.
    start_time=mpi_wtime()
    do while (state%continue_timestep)
      if (state%update_dtm) state%dtm=state%dtm_new
      ! The start of a timestep
      if (logging_mod_level .ge. LOG_DEBUG) start_iteration_time=mpi_wtime()
      call timestep(state) ! Call out to the timestepper to do the actual timestepping per component
      if (logging_mod_level .ge. LOG_DEBUG .and. state%parallel%my_rank==0) &
           call display_timestep_information(state%timestep, start_iteration_time)
      if (state%continue_timestep) then
        state%timestep = state%timestep+1
        state%time = state%time + state%dtm
      end if
    end do
    end_time=mpi_wtime()
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
    double precision, intent(in) :: start_time

    double precision :: end_time

    end_time=mpi_wtime()
    call log_log(LOG_DEBUG, "Timestep "//trim(conv_to_string(timestep))//" completed in "//&
         trim(conv_to_string(int((end_time-start_time) * 1000)))//"ms")
  end subroutine display_timestep_information  

  !> Registers each supplied component description
  subroutine fill_registry_with_components(options_database, component_descriptions)
    type(hashmap_type), intent(inout) :: options_database
    type(list_type), intent(inout) :: component_descriptions

    class(*), pointer :: raw_data
    type(iterator_type) :: iterator

    iterator=c_get_iterator(component_descriptions)
    do while (c_has_next(iterator))
      raw_data=>c_next_generic(iterator)
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
    type(hashmap_type), intent(inout) :: options_database
    character(len=*), intent(in) :: key
    
    if (options_has_key(options_database, key)) then
      is_present_and_true=options_get_logical(options_database, key)
      return
    end if
    is_present_and_true=.false.
  end function is_present_and_true

  !> Displays the registered components and their version numbers
  subroutine display_registed_components()
    type(map_type) :: registered_components
    type(iterator_type) :: iterator
    type(mapentry_type):: map_entry

    registered_components = get_all_registered_components()
    call log_log(LOG_INFO, "Registered components: "//conv_to_string(c_size(registered_components)))
    iterator=c_get_iterator(registered_components)
    do while (c_has_next(iterator))
      map_entry=c_next_mapentry(iterator)
      call log_log(LOG_INFO, trim(map_entry%key)//" "//trim(conv_to_string(c_get_real(map_entry))))
    end do    
  end subroutine display_registed_components

  !> Splits the MPI_COMM_WORLD communicator into MONC and IO separate communicators. The size of each depends
  !! on the stride supplied. This will deal with the case where you only have 1 extra process, for instance 3 MONCs to an
  !! IO server with 5 processes. 0=IO server, 1-3 are MONCS but by rights 4 would be an IO server. However we dont want to
  !! waste a process as an IO server which is not serving anything, hence in this edge case it will be used as a MONC instead
  !! @param io_stride The absolute process id stride for IO processes
  !! @param monc_communicator The communicator associated with MONC processes
  !! @param io_communicator The communicator associated with IO processes
  subroutine split_communicator_into_monc_and_io(moncs_per_io, monc_communicator, io_communicator, &
       am_i_monc_process, corresponding_io_server_process)
    integer, intent(in) :: moncs_per_io
    integer, intent(out) :: monc_communicator, io_communicator, corresponding_io_server_process
    logical, intent(out) :: am_i_monc_process

    integer, dimension(:), allocatable :: members_monc_group, members_io_group
    integer :: total_ranks, monc_group, io_group, io_processes, monc_processes, i, io_index, &
         monc_index, my_rank, ierr, global_group, io_stride

    call mpi_comm_size(MPI_COMM_WORLD, total_ranks, ierr)
    call mpi_comm_rank(MPI_COMM_WORLD, my_rank, ierr)

    io_stride=moncs_per_io+1
    io_processes=get_number_io_processes(total_ranks, io_stride)
    monc_processes=total_ranks-io_processes
    allocate(members_io_group(io_processes), members_monc_group(monc_processes))
    io_index=1
    monc_index=1
    corresponding_io_server_process=-1
    am_i_monc_process=.true.

    do i=0, total_ranks-1
      if (mod(i, io_stride) == 0 .and. i .lt. total_ranks) then
        if (io_index .le. io_processes) then
          members_io_group(io_index)=i
        else
          members_monc_group(monc_index)=i
          monc_index=monc_index+1
        end if        
        io_index=io_index+1
        if (my_rank == i) am_i_monc_process=.false.
        if (my_rank .gt. i .and. my_rank .lt. i+io_stride) then
          corresponding_io_server_process=i
        end if
      else
        members_monc_group(monc_index)=i
        monc_index=monc_index+1
      end if
    end do

    if (.not. am_i_monc_process .and. my_rank .eq. total_ranks-1) then
      am_i_monc_process=.true.
      corresponding_io_server_process=my_rank-io_stride
    end if

    if (am_i_monc_process .and. corresponding_io_server_process .lt. 0) then
      call log_log(LOG_ERROR, "MONC can not deduce its IO server rank, try with a different number of IO to MONC setting")
    end if    

    if (log_get_logging_level() .ge. LOG_DEBUG) then
      call log_log(LOG_DEBUG, "IO server assignment, rank="//conv_to_string(my_rank)//" IO server="//&
           conv_to_string(corresponding_io_server_process)//" am I a MONC="//conv_to_string(am_i_monc_process))
    end if    

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
  integer function get_number_io_processes(total_ranks, moncs_per_io)
    integer, intent(in) :: total_ranks, moncs_per_io

    integer :: io_stride

    io_stride=moncs_per_io
    get_number_io_processes=total_ranks/io_stride
    if (get_number_io_processes * io_stride .lt. total_ranks-1) get_number_io_processes=get_number_io_processes+1
  end function get_number_io_processes

  !> Reads the IO server configuration and populates the required variables of the configuration
  !! file name and the placement period
  !! @param options_database The options database
  subroutine get_io_configuration(options_database, ioserver_configuration_file, moncs_per_io_server)
    type(hashmap_type), intent(inout) :: options_database
    character(len=LONG_STRING_LENGTH), intent(out) :: ioserver_configuration_file
    integer, intent(out) :: moncs_per_io_server
   
    integer :: myrank, ierr

    ioserver_configuration_file=options_get_string(options_database, "ioserver_configuration_file")
    moncs_per_io_server=options_get_integer(options_database, "moncs_per_io_server")

    if (moncs_per_io_server == -1 .or. ioserver_configuration_file == "") then
      call mpi_comm_rank(MPI_COMM_WORLD, myrank, ierr)
      if (myrank == 0) call log_log(LOG_ERROR, "To run an IO server you must provide the placement period and configuration file")
      call mpi_barrier(MPI_COMM_WORLD) ! All other processes barrier here to ensure 0 displays the message before quit
      stop
    end if
  end subroutine get_io_configuration  

  !> Retrives the configured MPI threading mode, this is serialized by default but can be overridden via environment variable
  !! @returns The MONC MPI threading mode
  integer function get_mpi_threading_mode()
    character(len=STRING_LENGTH) :: thread_multiple_config_value
    integer :: status

    call get_environment_variable("MONC_THREAD_MULTIPLE", thread_multiple_config_value, status=status)

    if (status == 0 .and. conv_is_logical(trim(thread_multiple_config_value))) then
      if (conv_to_logical(trim(thread_multiple_config_value))) then
        get_mpi_threading_mode=MPI_THREAD_MULTIPLE
      else
        get_mpi_threading_mode=MPI_THREAD_SERIALIZED
      end if
    else
      get_mpi_threading_mode=MPI_THREAD_SERIALIZED
    end if
  end function get_mpi_threading_mode

  !> Converts an MPI threading level to the string representation of it
  !! @param lvl The integer MPI level
  !! @returns The string representation of the level
  character(len=STRING_LENGTH) function mpi_threading_level_to_string(lvl)
    integer, intent(in) :: lvl

    if (lvl == MPI_THREAD_SINGLE) then
      mpi_threading_level_to_string="single"
    else if (lvl == MPI_THREAD_FUNNELED) then
      mpi_threading_level_to_string="funneled"
    else if (lvl == MPI_THREAD_SERIALIZED) then
      mpi_threading_level_to_string="serialized"
    else if (lvl == MPI_THREAD_MULTIPLE) then
      mpi_threading_level_to_string="multiple"
    else
      mpi_threading_level_to_string="unknown"
    end if
  end function mpi_threading_level_to_string  

  !> Reads the NetCDF checkpoint file to obtain model time
  !! @param filename The filename of the checkpoint file to load
  !! @param checkpoint_time The model time from the the checkpoint file
  subroutine extract_time_from_checkpoint_file(filename, checkpoint_time)
    character(len=*), intent(in) :: filename
    real(kind=DEFAULT_PRECISION), intent(out) :: checkpoint_time

    integer :: ncid, variable_id
    real(kind=DEFAULT_PRECISION) :: r_data(1)

    call check_status(nf90_open(path = filename, mode = nf90_nowrite, ncid = ncid))
    call check_status(nf90_inq_varid(ncid, TIME_KEY, variable_id))
    call check_status(nf90_get_var(ncid, variable_id, r_data))

    checkpoint_time = r_data(1)

    call check_status(nf90_close(ncid))
  end subroutine extract_time_from_checkpoint_file


end module monc_mod
