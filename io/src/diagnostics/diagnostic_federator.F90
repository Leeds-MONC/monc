!> This diagnostics federator will take in data fields sent from a MONC, perform operators on these as required by the
!! diagnostics definition to produce diagnostic instantaneous fields which are then sent along to anything interested in them
module diagnostic_federator_mod
  use datadefn_mod, only : DEFAULT_PRECISION, STRING_LENGTH
  use configuration_parser_mod, only : io_configuration_type, io_configuration_misc_item_type, data_values_type, &
       io_configuration_field_type, io_configuration_diagnostic_field_type, get_number_field_dimensions, &
       get_data_value_by_field_name, get_monc_location
  use collections_mod, only : hashmap_type, hashset_type, map_type, list_type, c_get, c_size, c_put, c_contains, c_is_empty, &
       c_add, c_remove, c_value_at, c_key_at, c_free
  use conversions_mod, only : conv_to_generic, conv_to_string, conv_to_real, conv_to_integer, generic_to_double_real
  use reduction_inter_io_mod, only : init_reduction_inter_io, check_reduction_inter_io_for_completion, &
       finalise_reduction_inter_io, perform_inter_io_reduction, get_reduction_operator
  use broadcast_inter_io_mod, only : init_broadcast_inter_io, perform_inter_io_broadcast, finalise_broadcast_inter_io, &
       check_broadcast_inter_io_for_completion
  use allreduction_inter_io_mod, only : init_allreduction_inter_io, finalise_allreduction_inter_io, &
       perform_inter_io_allreduction, check_allreduction_inter_io_for_completion
  use global_callback_inter_io_mod, only : init_global_callback_inter_io, finalise_global_callback_inter_io
  use forthread_mod, only : forthread_rwlock_rdlock, forthread_rwlock_wrlock, forthread_rwlock_unlock, &
       forthread_rwlock_init, forthread_rwlock_destroy, forthread_mutex_init, forthread_mutex_lock, forthread_mutex_trylock, &
       forthread_mutex_unlock, forthread_mutex_destroy
  use threadpool_mod, only : check_thread_status
  use data_utils_mod, only : get_scalar_integer_from_monc, get_scalar_real_from_monc, is_field_present, &
       get_action_attribute_logical, get_action_attribute_integer, get_action_attribute_string, &
       get_array_double_from_monc, get_array_integer_from_monc
  use logging_mod, only : LOG_WARN, LOG_ERROR, log_log
  use operator_mod, only : perform_activity, initialise_operators, finalise_operators, get_operator_required_fields, &
       get_operator_perform_procedure, get_operator_auto_size
  use io_server_client_mod, only : DOUBLE_DATA_TYPE, INTEGER_DATA_TYPE
  use writer_field_manager_mod, only : provide_field_to_writer_federator
  implicit none

#ifndef TEST_MODE
  private
#endif

  !< The type of activity
  integer, parameter :: OPERATOR_TYPE=1, REDUCTION_TYPE=2, BROADCAST_TYPE=3, ALLREDUCTION_TYPE=4, PERFORM_CLEAN_EVERY=100

  !< A wrapper type containing all the diagnostics for MONC source processes at a specific timestep
  type all_diagnostics_at_timestep_type
     integer :: communication_corresponding_activities_rwlock, completed_diagnostics_rwlock
     type(list_type) :: diagnostic_entries
     type(hashset_type) :: completed_diagnostics
     type(hashmap_type) :: communication_corresponding_activities
  end type all_diagnostics_at_timestep_type  

  !< The diagnostics at a timestep for a specific MONC - this holds all the state required for executing the diagnostics
  type diagnostics_at_timestep_type
     integer :: timestep, completed_fields_rwlock, outstanding_fields_rwlock, activity_completion_mutex, source, &
          source_location, number_diags_outstanding
     real(kind=DEFAULT_PRECISION) :: time
     type(hashset_type) :: outstanding_fields, completed_activities
     type(hashmap_type) :: completed_fields
  end type diagnostics_at_timestep_type  

  !< A diagnostic which is a name and then the list of activities require to be executed
  type diagnostics_type
     character(len=STRING_LENGTH) :: diagnostic_name
     type(list_type) :: activities
     integer :: generation_timestep_frequency
     logical :: collective
  end type diagnostics_type  

  !< A diagnostic activity which is executed at some point with an input and returns an output
  type diagnostics_activity_type
     integer :: activity_type, communication_operator, root
     real(kind=DEFAULT_PRECISION) :: result
     type(list_type) :: required_fields
     type(map_type) :: activity_attributes
     character(len=STRING_LENGTH) :: result_name, activity_name
     procedure(perform_activity), pointer, nopass :: operator_procedure
  end type diagnostics_activity_type

  type(hashmap_type), volatile :: diagnostics_per_monc_at_timestep, all_diagnostics_at_timestep
  type(hashset_type), volatile :: all_outstanding_fields, available_fields
  type(diagnostics_type), volatile, dimension(:), allocatable :: diagnostic_definitions
  integer, volatile :: timestep_entries_rwlock, all_diagnostics_per_timestep_rwlock, clean_progress_mutex, &
       previous_clean_point, previous_viewed_timestep, current_point

 public initialise_diagnostic_federator, finalise_diagnostic_federator, check_diagnostic_federator_for_completion, &
      pass_fields_to_diagnostics_federator, determine_diagnostics_fields_available
contains  

  !> Initialises the diagnostics action and sets up the diagnostics master definitions
  !! @param io_configuration The IO server configuration
  !! @returns The map of diagnostic fields to the frequency (in timesteps) of generation
  type(hashmap_type) function initialise_diagnostic_federator(io_configuration)
    type(io_configuration_type), intent(inout) :: io_configuration

    call initialise_operators()
    call check_thread_status(forthread_rwlock_init(timestep_entries_rwlock, -1))
    call check_thread_status(forthread_rwlock_init(all_diagnostics_per_timestep_rwlock, -1))
    call check_thread_status(forthread_mutex_init(clean_progress_mutex, -1))    
    call init_reduction_inter_io(io_configuration)
    call init_broadcast_inter_io(io_configuration)
    call init_allreduction_inter_io(io_configuration)
    call init_global_callback_inter_io(io_configuration)
    call define_diagnostics(io_configuration, initialise_diagnostic_federator)
    previous_clean_point=0
    previous_viewed_timestep=0
    current_point=0
  end function initialise_diagnostic_federator
  
  !> Checks whether the diagnostics federator has completed or not, this is really checking all the underlying operations
  !! and communications to ensure that the data has been sent all the way through
  !! @param io_configuration Configuration state of the IO server
  !! @returns Whether the federator has completed
  logical function check_diagnostic_federator_for_completion(io_configuration)
    type(io_configuration_type), intent(inout) :: io_configuration

    check_diagnostic_federator_for_completion=check_reduction_inter_io_for_completion(io_configuration)
    if (check_diagnostic_federator_for_completion) then
      check_diagnostic_federator_for_completion=check_broadcast_inter_io_for_completion(io_configuration)
      if (check_diagnostic_federator_for_completion) then
        check_diagnostic_federator_for_completion=check_allreduction_inter_io_for_completion(io_configuration)
      end if      
    end if
  end function check_diagnostic_federator_for_completion

  !> Finalises the diagnostics federator, waiting for all outstanding requests and then freeing data
  !! @param io_configuration Configuration state of the IO server
  subroutine finalise_diagnostic_federator(io_configuration)
    type(io_configuration_type), intent(inout) :: io_configuration

    call finalise_broadcast_inter_io()
    call finalise_reduction_inter_io(io_configuration)
    call finalise_allreduction_inter_io(io_configuration)
    call finalise_global_callback_inter_io(io_configuration)
    call check_thread_status(forthread_rwlock_destroy(timestep_entries_rwlock))
    call check_thread_status(forthread_rwlock_destroy(all_diagnostics_per_timestep_rwlock))
    call check_thread_status(forthread_mutex_destroy(clean_progress_mutex))
    call finalise_operators()   
  end subroutine finalise_diagnostic_federator

  !> Determines the diagnostics fields that are available based upon the input MONC fields on registration
  !! that will be sent from MONC to the IO server
  !! @param monc_field_names Set of field names that are made available from MONC to the IO server
  !! @returns The set of available diagnostics
  type(hashset_type) function determine_diagnostics_fields_available(monc_field_names)
    type(hashset_type), intent(inout) :: monc_field_names

    integer :: i, j, k, num_activities, num_fields
    type(diagnostics_activity_type) :: specific_activity
    type(hashset_type) :: result_names_for_activities
    character(len=STRING_LENGTH) :: specific_field_name
    logical :: diagnostic_provided

    do i=1, size(diagnostic_definitions)
      diagnostic_provided=.true.
      num_activities=c_size(diagnostic_definitions(i)%activities)
      do j=1, num_activities
        specific_activity=get_activity_at_index(diagnostic_definitions(i)%activities, j)
        call c_add(result_names_for_activities, specific_activity%result_name)
      end do
      do j=1, num_activities
        specific_activity=get_activity_at_index(diagnostic_definitions(i)%activities, j)
        num_fields=c_size(specific_activity%required_fields)
        do k=1, num_fields
          specific_field_name=conv_to_string(c_get(specific_activity%required_fields, k), .false., STRING_LENGTH)          
          if (.not. c_contains(result_names_for_activities, specific_field_name) .and. &
               .not. c_contains(monc_field_names, specific_field_name)) then            
            diagnostic_provided=.false.
            exit
          end if          
        end do
        if (.not. diagnostic_provided) exit
      end do
      if (diagnostic_provided) then
        call c_add(determine_diagnostics_fields_available, diagnostic_definitions(i)%diagnostic_name)
        call c_add(available_fields, diagnostic_definitions(i)%diagnostic_name)
      end if
      call c_free(result_names_for_activities)
    end do    
  end function determine_diagnostics_fields_available  

  !> Entry point into the diagnostics federator this runs the diagnostics, executing the defined rules based upon the input
  !! data
  !! @param io_configuration Configuration of the IO server
  !! @param source The source PID of the MONC process
  !! @param data_id The ID of the data definition that is represented by the dump
  !! @param data_dump The data we have received from MONC
  subroutine pass_fields_to_diagnostics_federator(io_configuration, source, data_id, data_dump)
    type(io_configuration_type), intent(inout) :: io_configuration
    integer, intent(in) :: source, data_id
    character, dimension(:), allocatable :: data_dump

    integer :: timestep
    real(kind=DEFAULT_PRECISION) :: time
    type(diagnostics_at_timestep_type), pointer :: timestep_entry
    type(all_diagnostics_at_timestep_type), pointer :: diagnostics_by_timestep

    if (.not. allocated(diagnostic_definitions)) return

    if (is_field_present(io_configuration, source, data_id, "timestep") .and. &
         is_field_present(io_configuration, source, data_id, "time")) then
      timestep=get_scalar_integer_from_monc(io_configuration, source, data_id, data_dump, "timestep")
      time=get_scalar_real_from_monc(io_configuration, source, data_id, data_dump, "time")
      timestep_entry=>find_or_register_timestep_entry(io_configuration, timestep, source, time)
      diagnostics_by_timestep=>get_diagnostics_by_timestep(timestep, .true.)
      call clean_diagnostic_states(timestep)
      call check_diagnostics_entries_against_data(io_configuration, source, data_id, data_dump, timestep_entry)
      call issue_communication_calls(io_configuration, timestep_entry, diagnostics_by_timestep, source, data_id, data_dump)      
      call check_all_activities_against_completed_fields(io_configuration, timestep_entry, diagnostics_by_timestep)
    else
      call log_log(LOG_WARN, "Can not run the diagnostics federator without a timestep and time field in the MONC data")
    end if
  end subroutine pass_fields_to_diagnostics_federator

  !> Checks all pending activities against the completed fields and runs them if the required fields are now available
  !! @param io_configuration The IO server configuration
  !! @param timestep_entry The timestep entry that we are checking against
  subroutine check_all_activities_against_completed_fields(io_configuration, timestep_entry, diagnostics_by_timestep)
    type(io_configuration_type), intent(inout) :: io_configuration
    type(diagnostics_at_timestep_type), intent(inout) :: timestep_entry
    type(all_diagnostics_at_timestep_type), intent(inout) :: diagnostics_by_timestep

    integer :: i, j, entries, num_diags
    type(diagnostics_activity_type), pointer :: activity
    type(list_type) :: activities_to_remove
    character(len=STRING_LENGTH) :: field_name
    logical :: updated_entry, entry_in_completed_diagnostics, operator_produced_values
    type(data_values_type) :: value_to_send

    updated_entry=.true.
    
    call check_thread_status(forthread_mutex_lock(timestep_entry%activity_completion_mutex))
    do while (updated_entry)
      updated_entry=.false.      
      num_diags=size(diagnostic_definitions)
      do j=1, num_diags
        call check_thread_status(forthread_rwlock_rdlock(diagnostics_by_timestep%completed_diagnostics_rwlock))
        entry_in_completed_diagnostics=c_contains(diagnostics_by_timestep%completed_diagnostics, &
             diagnostic_definitions(j)%diagnostic_name)
        call check_thread_status(forthread_rwlock_unlock(diagnostics_by_timestep%completed_diagnostics_rwlock))
        if (diagnostic_definitions(j)%collective .or. .not. entry_in_completed_diagnostics) then
          entries=c_size(diagnostic_definitions(j)%activities)
          do i=1, entries
            if (.not. c_contains(timestep_entry%completed_activities, trim(diagnostic_definitions(j)%diagnostic_name)//"#"//&
                 trim(conv_to_string(i)))) then
              activity=>get_activity_at_index(diagnostic_definitions(j)%activities, i)      
              if (are_fields_available_for_activity(timestep_entry, activity)) then                
                call c_add(timestep_entry%completed_activities, trim(diagnostic_definitions(j)%diagnostic_name)//"#"//&
                     trim(conv_to_string(i)))
                updated_entry=.true.              
                if (activity%activity_type == OPERATOR_TYPE) then                        
                  operator_produced_values=handle_operator_completion(io_configuration, timestep_entry, activity)
                  call c_add(activities_to_remove, conv_to_generic(i, .true.))
                  if (operator_produced_values .and. activity%result_name == diagnostic_definitions(j)%diagnostic_name) then
                    call handle_diagnostic_calculation_completed(io_configuration, j, timestep_entry, diagnostics_by_timestep)
                  end if
                else if (activity%activity_type == REDUCTION_TYPE .or. activity%activity_type == BROADCAST_TYPE &
                     .or. activity%activity_type == ALLREDUCTION_TYPE) then
                  field_name=conv_to_string(c_get(activity%required_fields, 1), .false., STRING_LENGTH)
                  call check_thread_status(forthread_rwlock_rdlock(timestep_entry%completed_fields_rwlock))
                  value_to_send=get_data_value_by_field_name(timestep_entry%completed_fields, field_name)
                  call check_thread_status(forthread_rwlock_unlock(timestep_entry%completed_fields_rwlock))
                  call check_thread_status(forthread_mutex_unlock(timestep_entry%activity_completion_mutex))
                  call perform_inter_io_communication(io_configuration, timestep_entry, diagnostics_by_timestep, &
                       activity, value_to_send, field_name, diagnostic_definitions(j)%collective)
                  call check_thread_status(forthread_mutex_lock(timestep_entry%activity_completion_mutex))
                end if
              end if
            end if
          end do
        end if
      end do  
    end do
    call check_thread_status(forthread_mutex_unlock(timestep_entry%activity_completion_mutex))
  end subroutine check_all_activities_against_completed_fields

  !> Handles the completion of the operator
  !! @param timestep_entry The specific timestep entry
  !! @param specific_activity The specific activity that just completed
  logical function handle_operator_completion(io_configuration, timestep_entry, specific_activity)
    type(io_configuration_type), intent(inout) :: io_configuration
    type(diagnostics_at_timestep_type), intent(inout) :: timestep_entry
    type(diagnostics_activity_type), intent(inout) :: specific_activity

    type(data_values_type), pointer :: operator_result
    real(kind=DEFAULT_PRECISION), dimension(:), allocatable :: operator_result_values
    class(*), pointer :: generic

    call specific_activity%operator_procedure(io_configuration, timestep_entry%completed_fields, &
         specific_activity%activity_attributes, timestep_entry%source_location, timestep_entry%source, operator_result_values)    

    if (allocated(operator_result_values)) then
      allocate(operator_result)
      operator_result%values=operator_result_values
      deallocate(operator_result_values)
      generic=>operator_result
      call check_thread_status(forthread_rwlock_wrlock(timestep_entry%completed_fields_rwlock))
      call c_put(timestep_entry%completed_fields, specific_activity%result_name, generic)
      call check_thread_status(forthread_rwlock_unlock(timestep_entry%completed_fields_rwlock))
      call check_thread_status(forthread_rwlock_wrlock(timestep_entry%outstanding_fields_rwlock))       
      call c_remove(timestep_entry%outstanding_fields, specific_activity%result_name)        
      call check_thread_status(forthread_rwlock_unlock(timestep_entry%outstanding_fields_rwlock)) 
      handle_operator_completion=.true.
    else
      handle_operator_completion=.false.
    end if
  end function handle_operator_completion

  !> Handles inter io reduction completion, it adds the resulting value to the appropriate completion lists
  !! and then checks the pending activities and runs them if we can execute any of these based upon this value
  !! @param io_configuration Configuration of the IO server
  !! @param values Array of values resulting from the communication
  !! @param field_name Corresponding field name
  !! @param timestep Corresponding timestep
  subroutine handle_completion(io_configuration, values, field_name, timestep)
    type(io_configuration_type), intent(inout) :: io_configuration
    real(DEFAULT_PRECISION), dimension(:) :: values
    character(len=STRING_LENGTH) :: field_name
    integer :: timestep

    integer :: i, entries
    
    type(all_diagnostics_at_timestep_type), pointer :: diagnostics_by_timestep
    type(data_values_type), pointer :: result_to_add
    class(*), pointer :: generic

    diagnostics_by_timestep=>get_diagnostics_by_timestep(timestep, .true.)
    if (.not. c_is_empty(diagnostics_by_timestep%diagnostic_entries)) then
      entries=c_size(diagnostics_by_timestep%diagnostic_entries)
      do i=1, entries
        allocate(result_to_add)
        result_to_add%values=values
        call handle_completion_for_specific_monc_timestep_entry(io_configuration, result_to_add, &
             field_name, timestep, get_specific_monc_timestep_entry_at_index(diagnostics_by_timestep, i), &
             diagnostics_by_timestep)
      end do
    end if
    call check_thread_status(forthread_rwlock_wrlock(diagnostics_by_timestep%communication_corresponding_activities_rwlock))
    call c_remove(diagnostics_by_timestep%communication_corresponding_activities, trim(field_name))
    call check_thread_status(forthread_rwlock_unlock(diagnostics_by_timestep%communication_corresponding_activities_rwlock))
  end subroutine handle_completion

  !> This handles inter IO completion for a specific timestep entry. This is required as at a timestep there can be
  !! multiple entires based on each MONC which communicates with the IO server so we handle an individual one here
  !! @param io_configuration Configuration of the IO server
  !! @param values Array of values resulting from the communication
  !! @param field_name Corresponding field name
  !! @param timestep Corresponding timestep
  !! @param timestep_entry The specific timestep entry which relates to an individual MONC communicating with the IO server
  subroutine handle_completion_for_specific_monc_timestep_entry(io_configuration, result_to_add, &
       field_name, timestep, timestep_entry, diagnostics_by_timestep)
    type(io_configuration_type), intent(inout) :: io_configuration
    type(data_values_type), pointer, intent(in) :: result_to_add
    character(len=STRING_LENGTH) :: field_name
    integer :: timestep, i
    type(diagnostics_at_timestep_type), pointer :: timestep_entry
    type(all_diagnostics_at_timestep_type), intent(inout) :: diagnostics_by_timestep

    logical :: entry_in_completed_diagnostics
    type(diagnostics_activity_type), pointer :: activity
    class(*), pointer :: generic
    
    activity=>get_comm_activity_from_fieldname(diagnostics_by_timestep, field_name)
        
    generic=>result_to_add
    call check_thread_status(forthread_rwlock_wrlock(timestep_entry%completed_fields_rwlock))
    call c_put(timestep_entry%completed_fields, trim(activity%result_name), generic)
    call check_thread_status(forthread_rwlock_unlock(timestep_entry%completed_fields_rwlock))

    call check_thread_status(forthread_rwlock_wrlock(timestep_entry%completed_fields_rwlock))
    ! is field name here correct?
    call c_remove(timestep_entry%outstanding_fields, trim(field_name))
    call check_thread_status(forthread_rwlock_unlock(timestep_entry%completed_fields_rwlock))
    
    do i=1, size(diagnostic_definitions)
      if (activity%result_name == diagnostic_definitions(i)%diagnostic_name) then                      
        if (diagnostic_definitions(i)%collective) then          
          call handle_diagnostic_calculation_completed(io_configuration, i, timestep_entry, diagnostics_by_timestep)
        else
          call check_thread_status(forthread_rwlock_rdlock(diagnostics_by_timestep%completed_diagnostics_rwlock))
          entry_in_completed_diagnostics=c_contains(diagnostics_by_timestep%completed_diagnostics, &
               diagnostic_definitions(i)%diagnostic_name)
          call check_thread_status(forthread_rwlock_unlock(diagnostics_by_timestep%completed_diagnostics_rwlock))
          if (.not. entry_in_completed_diagnostics) then
            call handle_diagnostic_calculation_completed(io_configuration, i, timestep_entry, diagnostics_by_timestep)
          end if
        end if
        exit
      end if
    end do

    call check_all_activities_against_completed_fields(io_configuration, timestep_entry, diagnostics_by_timestep)
  end subroutine handle_completion_for_specific_monc_timestep_entry

  !> Handles completion of a diagnostic calculation and will then pass this onto interested parties
  !! @param diagnostic_index The diagnostics which has completed
  !! @param timestep_entry The timestep entry which has just completed this diagnostic
  subroutine handle_diagnostic_calculation_completed(io_configuration, diagnostic_index, timestep_entry, diagnostics_by_timestep)
    type(io_configuration_type), intent(inout) :: io_configuration
    integer, intent(in) :: diagnostic_index
    type(diagnostics_at_timestep_type), intent(inout) :: timestep_entry
    type(all_diagnostics_at_timestep_type), intent(inout) :: diagnostics_by_timestep

    type(data_values_type), pointer :: diagnostics_value_entry
    integer :: i, entries
    type(diagnostics_at_timestep_type), pointer :: activity_at_index

    call check_thread_status(forthread_rwlock_rdlock(timestep_entry%completed_fields_rwlock))
    diagnostics_value_entry=>get_data_value_by_field_name(timestep_entry%completed_fields, &
         trim(diagnostic_definitions(diagnostic_index)%diagnostic_name))
    call check_thread_status(forthread_rwlock_unlock(timestep_entry%completed_fields_rwlock))

    if (diagnostic_definitions(diagnostic_index)%collective) then
      timestep_entry%number_diags_outstanding=timestep_entry%number_diags_outstanding-1
      call provide_field_to_writer_federator(io_configuration, diagnostic_definitions(diagnostic_index)%diagnostic_name, &
               diagnostics_value_entry%values, timestep_entry%timestep, timestep_entry%time, &
               diagnostic_definitions(diagnostic_index)%generation_timestep_frequency, timestep_entry%source)
      if (allocated(diagnostics_value_entry%values)) deallocate(diagnostics_value_entry%values)
    else
      call check_thread_status(forthread_rwlock_rdlock(diagnostics_by_timestep%completed_diagnostics_rwlock))
      if (.not. c_contains(diagnostics_by_timestep%completed_diagnostics, &
           diagnostic_definitions(diagnostic_index)%diagnostic_name)) then
        call check_thread_status(forthread_rwlock_unlock(diagnostics_by_timestep%completed_diagnostics_rwlock))
        call check_thread_status(forthread_rwlock_wrlock(diagnostics_by_timestep%completed_diagnostics_rwlock))
        if (.not. c_contains(diagnostics_by_timestep%completed_diagnostics,&
             diagnostic_definitions(diagnostic_index)%diagnostic_name)) then
          call c_add(diagnostics_by_timestep%completed_diagnostics, diagnostic_definitions(diagnostic_index)%diagnostic_name)
          call check_thread_status(forthread_rwlock_unlock(diagnostics_by_timestep%completed_diagnostics_rwlock))
          entries=c_size(diagnostics_by_timestep%diagnostic_entries)
          do i=1, entries
            activity_at_index=>get_specific_monc_timestep_entry_at_index(diagnostics_by_timestep, i)
            activity_at_index%number_diags_outstanding=activity_at_index%number_diags_outstanding-1
          end do
          call provide_field_to_writer_federator(io_configuration, diagnostic_definitions(diagnostic_index)%diagnostic_name, &
               diagnostics_value_entry%values, timestep_entry%timestep, timestep_entry%time, &
               diagnostic_definitions(diagnostic_index)%generation_timestep_frequency)
          if (allocated(diagnostics_value_entry%values)) deallocate(diagnostics_value_entry%values)
        else
          call check_thread_status(forthread_rwlock_unlock(diagnostics_by_timestep%completed_diagnostics_rwlock))
        end if
      else
        call check_thread_status(forthread_rwlock_unlock(diagnostics_by_timestep%completed_diagnostics_rwlock))
      end if
    end if        
  end subroutine handle_diagnostic_calculation_completed  

  !> Determines whether the fields required for an activity are available so that activity can be run
  !! @param timestep_entry The timestep entry that we are checking the completed fields for
  !! @param activity The activity that we are checking
  !! @returns Whether the required fields for the activity are available
  logical function are_fields_available_for_activity(timestep_entry, activity)
    type(diagnostics_at_timestep_type), intent(inout) :: timestep_entry
    type(diagnostics_activity_type) :: activity

    integer :: i, entries
    character(len=STRING_LENGTH) :: field_name

    are_fields_available_for_activity=.true.
    entries=c_size(activity%required_fields)
    if (entries .gt. 0) then
      call check_thread_status(forthread_rwlock_rdlock(timestep_entry%completed_fields_rwlock))
      do i=1, entries
        field_name=conv_to_string(c_get(activity%required_fields, i), .true., STRING_LENGTH)
        if (.not. c_contains(timestep_entry%completed_fields, field_name)) then          
          are_fields_available_for_activity=.false.
          exit
        else
        end if
      end do
      call check_thread_status(forthread_rwlock_unlock(timestep_entry%completed_fields_rwlock))
    end if
  end function are_fields_available_for_activity

  !> Performs the actual inter IO communication by calling out to the appropriate inter IO module
  !! @param io_configuration Configuration of the IO server
  !! @param timestep_entry The timestep entry
  !! @param activity The activity this is executing
  !! @param value_to_send The value to communicate
  !! @param communication_field_name Name of the field that we are communicating
  subroutine perform_inter_io_communication(io_configuration, timestep_entry, all_entries_at_timestep, &
       activity, value_to_send, communication_field_name, collective_diagnostic)
    type(io_configuration_type), intent(inout) :: io_configuration
    type(diagnostics_at_timestep_type), intent(inout) :: timestep_entry
    type(all_diagnostics_at_timestep_type), intent(inout) :: all_entries_at_timestep
    type(diagnostics_activity_type), pointer, intent(in) :: activity
    type(data_values_type), intent(in) :: value_to_send
    character(len=STRING_LENGTH), intent(in) :: communication_field_name
    logical, intent(in) :: collective_diagnostic

    class(*), pointer :: generic

    generic=>activity
    call check_thread_status(forthread_rwlock_wrlock(all_entries_at_timestep%communication_corresponding_activities_rwlock))
    call c_put(all_entries_at_timestep%communication_corresponding_activities, trim(communication_field_name), generic)
    call check_thread_status(forthread_rwlock_unlock(all_entries_at_timestep%communication_corresponding_activities_rwlock))

    if (activity%activity_type == REDUCTION_TYPE) then
      call perform_inter_io_reduction(io_configuration, value_to_send%values, size(value_to_send%values), &
           communication_field_name, activity%communication_operator, activity%root, timestep_entry%timestep, handle_completion)
      if (activity%root .ne. io_configuration%my_io_rank) then
          timestep_entry%number_diags_outstanding=timestep_entry%number_diags_outstanding-1
      end if
    else if (activity%activity_type == ALLREDUCTION_TYPE) then
      call perform_inter_io_allreduction(io_configuration, value_to_send%values, size(value_to_send%values), &
           communication_field_name, activity%communication_operator, activity%root, timestep_entry%timestep, handle_completion)
    else if (activity%activity_type == BROADCAST_TYPE) then
      call perform_inter_io_broadcast(io_configuration, value_to_send%values, size(value_to_send%values), &
           communication_field_name, activity%root, timestep_entry%timestep, handle_completion)
    end if
  end subroutine perform_inter_io_communication
  
  !> Issues any inter io communucation calls that are appropriate based upon the data recieved from MONC
  !! @param io_configuration Configuration of the IO server
  !! @param timestep_entry The timestep entry
  !! @param source The source PID of the MONC process
  !! @param data_id The ID of the data definition that is represented by the dump
  !! @param data_dump The data we have received from MONC
  subroutine issue_communication_calls(io_configuration, timestep_entry, diagnostics_by_timestep, source, data_id, data_dump)
    type(io_configuration_type), intent(inout) :: io_configuration
    type(diagnostics_at_timestep_type), intent(inout) :: timestep_entry
    integer, intent(in) :: source, data_id
    character, dimension(:), allocatable :: data_dump
    type(all_diagnostics_at_timestep_type), intent(inout) :: diagnostics_by_timestep

    logical :: completed_diagnostics_entry
    integer :: i, j, num_diags, entries
    type(diagnostics_activity_type), pointer :: activity
    character(len=STRING_LENGTH) :: communication_field_name
    
    num_diags=size(diagnostic_definitions)
    do j=1, num_diags
      call check_thread_status(forthread_rwlock_rdlock(diagnostics_by_timestep%completed_diagnostics_rwlock))
      completed_diagnostics_entry=c_contains(diagnostics_by_timestep%completed_diagnostics, &
           diagnostic_definitions(j)%diagnostic_name)
      call check_thread_status(forthread_rwlock_unlock(diagnostics_by_timestep%completed_diagnostics_rwlock))
      if (.not. completed_diagnostics_entry) then
        entries=c_size(diagnostic_definitions(j)%activities)
        do i=1, entries
          if (.not. c_contains(timestep_entry%completed_activities, trim(diagnostic_definitions(j)%diagnostic_name)//"#"//&
               trim(conv_to_string(i)))) then
            activity=>get_activity_at_index(diagnostic_definitions(j)%activities, i)
            if ((activity%activity_type == REDUCTION_TYPE .or. activity%activity_type == BROADCAST_TYPE &
                 .or. activity%activity_type == ALLREDUCTION_TYPE)) then
              if (.not. c_is_empty(activity%required_fields)) then
                communication_field_name=conv_to_string(c_get(activity%required_fields, 1), .false., STRING_LENGTH)
                if (is_field_present(io_configuration, source, data_id, communication_field_name)) then
                  call c_add(timestep_entry%completed_activities, trim(diagnostic_definitions(j)%diagnostic_name)//"#"//&
                       trim(conv_to_string(i)))                  
                  call perform_inter_io_communication(io_configuration, timestep_entry, diagnostics_by_timestep, activity, &
                       get_value_from_monc_data(io_configuration, source, data_id, data_dump, communication_field_name), &
                       communication_field_name, diagnostic_definitions(j)%collective)            
                end if
              end if
            end if
          end if
        end do
      end if
    end do
  end subroutine issue_communication_calls

  !> Checks the outstanding fields of a time step entry against the data recieved from MONC and moves any acquired data
  !! from the outstanding set to the completed, along with its value
  !! @param io_configuration Configuration of the IO server  
  !! @param source The source PID of the MONC process
  !! @param data_id The ID of the data definition that is represented by the dump
  !! @param data_dump The data we have received from MONC
  !! @param timestep_diagnostics_entry The timestep entry
  subroutine check_diagnostics_entries_against_data(io_configuration, source, data_id, data_dump, &
       timestep_diagnostics_entry)
    type(io_configuration_type), intent(inout) :: io_configuration
    integer, intent(in) :: source, data_id
    character, dimension(:), allocatable :: data_dump
    type(diagnostics_at_timestep_type), intent(inout) :: timestep_diagnostics_entry

    integer :: i, entries
    character(len=STRING_LENGTH) :: field_name
    type(data_values_type), pointer :: field_value
    type(hashset_type) :: removed_entries
    class(*), pointer :: generic
        
    if (.not. c_is_empty(timestep_diagnostics_entry%outstanding_fields)) then
      call check_thread_status(forthread_rwlock_rdlock(timestep_diagnostics_entry%outstanding_fields_rwlock))
      entries=c_size(timestep_diagnostics_entry%outstanding_fields)
      do i=1, entries
        field_name=c_get(timestep_diagnostics_entry%outstanding_fields, i)
        if (is_field_present(io_configuration, source, data_id, field_name)) then
          field_value=>get_value_from_monc_data(io_configuration, source, data_id, data_dump, field_name)
          generic=>field_value
          call check_thread_status(forthread_rwlock_wrlock(timestep_diagnostics_entry%completed_fields_rwlock))
          call c_put(timestep_diagnostics_entry%completed_fields, trim(field_name), generic)
          call check_thread_status(forthread_rwlock_unlock(timestep_diagnostics_entry%completed_fields_rwlock))
          call c_add(removed_entries, field_name)
        end if
      end do
      call check_thread_status(forthread_rwlock_unlock(timestep_diagnostics_entry%outstanding_fields_rwlock))
      if (.not. c_is_empty(removed_entries)) then
        entries=c_size(removed_entries)
        call check_thread_status(forthread_rwlock_wrlock(timestep_diagnostics_entry%outstanding_fields_rwlock))
        do i=1, entries
          call c_remove(timestep_diagnostics_entry%outstanding_fields, c_get(removed_entries, i))
        end do
        call check_thread_status(forthread_rwlock_unlock(timestep_diagnostics_entry%outstanding_fields_rwlock))
      end if
    end if
  end subroutine check_diagnostics_entries_against_data

  !> Cleans the diagnostic states if required (based on the timestep period)
  !! @param current_timestep The current timestep which is checked against the previous run timestep
  subroutine clean_diagnostic_states(current_timestep)
    integer, intent(in) :: current_timestep

    integer :: i, j, entries, ts_entries, have_lock
    type(list_type) :: entries_to_remove
    type(all_diagnostics_at_timestep_type), pointer :: specific_all_diagnostics_for_ts
    type(diagnostics_at_timestep_type), pointer :: specific_monc_timestep_entry
    logical :: all_completed
    character(len=STRING_LENGTH) :: entry_key

    have_lock=forthread_mutex_trylock(clean_progress_mutex)
    if (have_lock == 0) then
      if (previous_viewed_timestep .ne. current_timestep) then
        current_point=current_point+1
        previous_viewed_timestep=current_timestep
      end if
      if (previous_clean_point + PERFORM_CLEAN_EVERY .lt. current_point) then
        previous_clean_point=current_point
        call check_thread_status(forthread_rwlock_rdlock(all_diagnostics_per_timestep_rwlock))
        entries=c_size(all_diagnostics_at_timestep)
        do i=1, entries
          specific_all_diagnostics_for_ts=>get_diagnostics_at_index(i, .false., entry_key)
          ts_entries=c_size(specific_all_diagnostics_for_ts%diagnostic_entries)
          all_completed=.true.
          do j=1, ts_entries
            specific_monc_timestep_entry=>get_specific_monc_timestep_entry_at_index(specific_all_diagnostics_for_ts, j)
            if (specific_monc_timestep_entry%number_diags_outstanding .gt. 0) then
              all_completed=.false.
              exit
            end if
          end do
          if (all_completed) then
            call c_add(entries_to_remove, conv_to_generic(entry_key, .true.))
          end if
        end do
        call check_thread_status(forthread_rwlock_unlock(all_diagnostics_per_timestep_rwlock))

        if (.not. c_is_empty(entries_to_remove)) then
          call check_thread_status(forthread_rwlock_wrlock(all_diagnostics_per_timestep_rwlock))
          entries=c_size(entries_to_remove)
          do i=1, entries
            entry_key=conv_to_string(c_get(entries_to_remove, i), .false., STRING_LENGTH)
            call deallocate_diagnostics_at_timestep(entry_key)
            call c_remove(all_diagnostics_at_timestep, entry_key)
          end do
          call check_thread_status(forthread_rwlock_unlock(all_diagnostics_per_timestep_rwlock))
          call c_free(entries_to_remove)
        end if
      end if
      call check_thread_status(forthread_mutex_unlock(clean_progress_mutex))
    end if
  end subroutine clean_diagnostic_states 

  !> Deallocates all the diagnostics at a specific timestep, this removes all the individual MONC timestep entries and deallocates
  !! internal all diagnostic data, but keeps the all diagnostic entry in the list (which is removed by the caller)
  !! @param key The look up key which corresponds to the all diagnostics entry
  subroutine deallocate_diagnostics_at_timestep(key)
    character(len=*), intent(in) :: key

    type(all_diagnostics_at_timestep_type), pointer :: all_diagnostics_at_ts
    type(diagnostics_at_timestep_type), pointer :: specific_monc_timestep_entry
    type(data_values_type), pointer :: field_data_value
    integer :: i, entries, cfentries, j
    class(*), pointer :: generic
    all_diagnostics_at_ts=>get_diagnostic_by_key(key)
    if (associated(all_diagnostics_at_ts)) then
      entries=c_size(all_diagnostics_at_ts%diagnostic_entries)
      do i=1, entries
        specific_monc_timestep_entry=>get_specific_monc_timestep_entry_at_index(all_diagnostics_at_ts, i)
        call check_thread_status(forthread_rwlock_destroy(specific_monc_timestep_entry%completed_fields_rwlock))
        call check_thread_status(forthread_rwlock_destroy(specific_monc_timestep_entry%outstanding_fields_rwlock))
        call check_thread_status(forthread_mutex_destroy(specific_monc_timestep_entry%activity_completion_mutex))
        call c_free(specific_monc_timestep_entry%outstanding_fields)
        call c_free(specific_monc_timestep_entry%completed_activities)
        cfentries=c_size(specific_monc_timestep_entry%completed_fields)
        do j=1, cfentries
          field_data_value=>get_data_value_by_field_name(specific_monc_timestep_entry%completed_fields, &
               c_key_at(specific_monc_timestep_entry%completed_fields, j))          
          if (allocated(field_data_value%values)) deallocate(field_data_value%values)
          deallocate(field_data_value)
        end do                
        call c_free(specific_monc_timestep_entry%completed_fields)
      end do
      call c_free(all_diagnostics_at_ts%completed_diagnostics)
      call c_free(all_diagnostics_at_ts%communication_corresponding_activities)
      entries=c_size(all_diagnostics_at_ts%diagnostic_entries)
      do i=1, entries
        generic=>c_get(all_diagnostics_at_ts%diagnostic_entries, i)
        if (associated(generic)) deallocate(generic)
      end do
      call c_free(all_diagnostics_at_ts%diagnostic_entries)
      call check_thread_status(forthread_rwlock_destroy(all_diagnostics_at_ts%communication_corresponding_activities_rwlock))
      call check_thread_status(forthread_rwlock_destroy(all_diagnostics_at_ts%completed_diagnostics_rwlock))
    end if
  end subroutine deallocate_diagnostics_at_timestep

  !> Retrieves a value from the communicated MONC data. If this was an integer then converts to a real
  !! @param io_configuration Configuration of the IO server  
  !! @param source The source PID of the MONC process
  !! @param data_id The ID of the data definition that is represented by the dump
  !! @param data_dump The data we have received from MONC
  !! @param field_name The field to retrieve
  !! @returns Data value wrapper containing the values and meta-data
  function get_value_from_monc_data(io_configuration, source, data_id, data_dump, field_name)
    type(io_configuration_type), intent(inout) :: io_configuration
    integer, intent(in) :: source, data_id
    character, dimension(:), allocatable :: data_dump
    character(len=*), intent(in) :: field_name
    type(data_values_type), pointer :: get_value_from_monc_data

    integer :: field_data_type, i
    integer, dimension(:), allocatable :: int_values
    
    allocate(get_value_from_monc_data)
    field_data_type=get_datatype_of_field(io_configuration%data_definitions(data_id)%fields, field_name)
    if (field_data_type == 0) then
      call log_log(LOG_ERROR, "No data type for field '"//trim(field_name)//"'")
    end if
    get_value_from_monc_data%dimensions=get_number_field_dimensions(io_configuration, field_name, source, data_id)
    if (field_data_type == DOUBLE_DATA_TYPE) then
      get_value_from_monc_data%values=get_array_double_from_monc(io_configuration, source, data_id, data_dump, field_name)
    else if (field_data_type == INTEGER_DATA_TYPE) then
      int_values=get_array_integer_from_monc(io_configuration, source, data_id, data_dump, field_name)
      allocate(get_value_from_monc_data%values(size(int_values)))
      do i=1, size(int_values)
        get_value_from_monc_data%values(i)=conv_to_real(int_values(i))
      end do
      deallocate(int_values)
    end if    
  end function get_value_from_monc_data

  !> Retrieves the data type of a field or 0 if the field was not found
  !! @param fields Array of fields to search
  !! @param field_name The name of the field to locate
  !! @returns The data type of this field or 0 if the field was not found
  integer function get_datatype_of_field(fields, field_name)
    type(io_configuration_field_type), dimension(:), intent(in) :: fields
    character(len=*), intent(in) :: field_name

    integer :: i

    do i=1, size(fields)
      if (fields(i)%name .eq. field_name) then
        get_datatype_of_field=fields(i)%data_type
        return
      end if      
    end do
    get_datatype_of_field=0
  end function get_datatype_of_field

  !> Locates or registers a new (if it does not exist) time step entry based upon the timestep and source MONC process. This is
  !! a timestep/source specific state which represents progress through the overall diagnostic configuration
  !! @param timestep The timestep that we are locating for
  !! @param source The source MONC process ID
  !! @returns The specific timestep entry
  function find_or_register_timestep_entry(io_configuration, timestep, source, time)
    type(io_configuration_type), intent(inout) :: io_configuration
    integer, intent(in) :: timestep, source
    real(kind=DEFAULT_PRECISION), intent(in) :: time
    type(diagnostics_at_timestep_type), pointer :: find_or_register_timestep_entry

    class(*), pointer :: generic
    type(all_diagnostics_at_timestep_type), pointer :: all_diags_by_timestep

    find_or_register_timestep_entry=>get_timestep_entry(timestep, source, .true.)
    if (.not. associated(find_or_register_timestep_entry)) then
      call check_thread_status(forthread_rwlock_wrlock(timestep_entries_rwlock))
      find_or_register_timestep_entry=>get_timestep_entry(timestep, source, .false.)
      if (.not. associated(find_or_register_timestep_entry)) then
        find_or_register_timestep_entry=>create_timestep_entry(io_configuration, timestep, time, source)
        generic=>find_or_register_timestep_entry
        call c_put(diagnostics_per_monc_at_timestep, conv_to_string(timestep)//"#"//conv_to_string(source), generic)
        all_diags_by_timestep=>find_or_add_diagnostics_by_timestep(timestep)
        call check_thread_status(forthread_rwlock_wrlock(all_diagnostics_per_timestep_rwlock))
        call c_add(all_diags_by_timestep%diagnostic_entries, generic)
        call check_thread_status(forthread_rwlock_unlock(all_diagnostics_per_timestep_rwlock))
      end if
      call check_thread_status(forthread_rwlock_unlock(timestep_entries_rwlock))
    end if    
  end function find_or_register_timestep_entry  

  !> Creates a timestep entry and processes all members, determining activities their required fields etc...
  !! @param timestep The timestep that we are locating for
  !! @returns The new timestep entry
  function create_timestep_entry(io_configuration, timestep, time, source)
    type(io_configuration_type), intent(inout) :: io_configuration
    integer, intent(in) :: timestep, source
    real(kind=DEFAULT_PRECISION), intent(in) :: time
    type(diagnostics_at_timestep_type), pointer :: create_timestep_entry

    integer :: i, entries

    allocate(create_timestep_entry)
    entries=c_size(all_outstanding_fields)
    do i=1, entries
      create_timestep_entry%timestep=timestep
      create_timestep_entry%time=time
      create_timestep_entry%number_diags_outstanding=c_size(available_fields)
      create_timestep_entry%source=source
      create_timestep_entry%source_location=get_monc_location(io_configuration, source)
      call c_add(create_timestep_entry%outstanding_fields, c_get(all_outstanding_fields, i))
      call check_thread_status(forthread_mutex_init(create_timestep_entry%activity_completion_mutex, -1))
      call check_thread_status(forthread_rwlock_init(create_timestep_entry%completed_fields_rwlock, -1))
      call check_thread_status(forthread_rwlock_init(create_timestep_entry%outstanding_fields_rwlock, -1))
    end do
  end function create_timestep_entry

  !> Adds the required fields of an activity to the overall required fields which are cloned for each new timestep entry
  !! @param required_fields The required fields of an activity
  subroutine add_required_fields_if_needed(required_fields)
    type(list_type), intent(inout) :: required_fields

    integer :: i, entries
    character(len=STRING_LENGTH) :: field_name

    entries=c_size(required_fields)
    if (entries .gt. 0) then
      do i=1, entries
        field_name=conv_to_string(c_get(required_fields, i), .false., STRING_LENGTH)
        call c_add(all_outstanding_fields, field_name)
      end do     
    end if
  end subroutine add_required_fields_if_needed

  !> Retrieves a specific activity at an index
  !! @param activities The list of activities to look up
  !! @param index The index to retrieve the activity at
  !! @returns The activity at the specific index or null if none is found
  function get_activity_at_index(activities, index)
    integer, intent(in) :: index
    type(list_type), intent(inout) :: activities
    type(diagnostics_activity_type), pointer :: get_activity_at_index

    class(*), pointer :: generic

    generic=>c_get(activities, index)
    if (associated(generic)) then
      select type(generic)
        type is(diagnostics_activity_type)
          get_activity_at_index=>generic
      end select      
    else
      get_activity_at_index=>null()
    end if
  end function get_activity_at_index  

  !> Retrieves the timestep at a specific timestep and source MONC
  !! @param timestep The timestep to look up
  !! @param source The source MONC process id
  !! @param do_lock Whether to issue a read lock over the timeseries collection
  !! @returns The timestep entry or null if none is found
  function get_timestep_entry(timestep, source, do_lock)
    integer, intent(in) :: timestep, source
    logical, intent(in) :: do_lock
    type(diagnostics_at_timestep_type), pointer :: get_timestep_entry
    
    class(*), pointer :: generic

    if (do_lock) call check_thread_status(forthread_rwlock_rdlock(timestep_entries_rwlock))
    generic=>c_get(diagnostics_per_monc_at_timestep, conv_to_string(timestep)//"#"//conv_to_string(source))
    if (do_lock) call check_thread_status(forthread_rwlock_unlock(timestep_entries_rwlock))
    if (associated(generic)) then
      select type(generic)
        type is(diagnostics_at_timestep_type)
          get_timestep_entry=>generic
      end select      
    else
      get_timestep_entry=>null()
    end if
  end function get_timestep_entry

  !> Finds or adds diagnostics by timestep. This is used to maintain a list of all diagnostic entries for a specific timestep
  !! for every MONC source process
  !! @param timestep The corresponding timestep
  !! @returns The diagnostics at the timestep or a new diagnostics if none was found
  function find_or_add_diagnostics_by_timestep(timestep)
    integer, intent(in) :: timestep
    type(all_diagnostics_at_timestep_type), pointer :: find_or_add_diagnostics_by_timestep

    class(*), pointer :: generic

    find_or_add_diagnostics_by_timestep=>get_diagnostics_by_timestep(timestep, .true.)
    if (.not. associated(find_or_add_diagnostics_by_timestep)) then
      call check_thread_status(forthread_rwlock_wrlock(all_diagnostics_per_timestep_rwlock))
      find_or_add_diagnostics_by_timestep=>get_diagnostics_by_timestep(timestep, .false.)
      if (.not. associated(find_or_add_diagnostics_by_timestep)) then
        allocate(find_or_add_diagnostics_by_timestep)
        call check_thread_status(forthread_rwlock_init(&
             find_or_add_diagnostics_by_timestep%communication_corresponding_activities_rwlock, -1))
        call check_thread_status(forthread_rwlock_init(find_or_add_diagnostics_by_timestep%completed_diagnostics_rwlock, -1))
        generic=>find_or_add_diagnostics_by_timestep
        call c_put(all_diagnostics_at_timestep, conv_to_string(timestep), generic)
      end if
      call check_thread_status(forthread_rwlock_unlock(all_diagnostics_per_timestep_rwlock))
    end if
  end function find_or_add_diagnostics_by_timestep

  !> Retrieves the diagnostics list (each MONC source) at a specific timestep
  !! @param timestep The timestep to look up
  !! @param do_lock Whether to issue read locks
  !! @returns The diagnostics container at this timestep or null if none is found
  function get_diagnostics_by_timestep(timestep, do_lock)
    integer, intent(in) :: timestep
    logical, intent(in) :: do_lock
    type(all_diagnostics_at_timestep_type), pointer :: get_diagnostics_by_timestep

    class(*), pointer :: generic
    
    if (do_lock) call check_thread_status(forthread_rwlock_rdlock(all_diagnostics_per_timestep_rwlock))
    generic=>c_get(all_diagnostics_at_timestep, conv_to_string(timestep))
    if (do_lock) call check_thread_status(forthread_rwlock_unlock(all_diagnostics_per_timestep_rwlock))
    if (associated(generic)) then
      select type(generic)
      type is(all_diagnostics_at_timestep_type)
        get_diagnostics_by_timestep=>generic
      end select
    else
      get_diagnostics_by_timestep=>null()
    end if
  end function get_diagnostics_by_timestep

  !> Retrieves a specific MONC timestep entry from the all diagnostics at timestep based upon its index
  !! @param all_diags_by_timestep All diagnostics by timestep
  !! @param index The index to look up
  !! @returns The specific MONC timestep entry or null if none is found
  function get_specific_monc_timestep_entry_at_index(all_diags_by_timestep, index)
    type(all_diagnostics_at_timestep_type), intent(inout) :: all_diags_by_timestep
    integer, intent(in) :: index
    type(diagnostics_at_timestep_type), pointer :: get_specific_monc_timestep_entry_at_index

    class(*), pointer :: generic

    generic=>c_get(all_diags_by_timestep%diagnostic_entries, index)

    if (associated(generic)) then
      select type(generic)
      type is (diagnostics_at_timestep_type)
        get_specific_monc_timestep_entry_at_index=>generic
      end select
    else
      get_specific_monc_timestep_entry_at_index=>null()
    end if
  end function get_specific_monc_timestep_entry_at_index

  !> Retrieves all diagnostics at a timestep by its key
  !! @param key The key to look up
  !! @returns The corresponding all diagnostics at this timestep or null if none is found
  function get_diagnostic_by_key(key)
    character(len=*), intent(in) :: key

    type(all_diagnostics_at_timestep_type), pointer :: get_diagnostic_by_key

    class(*), pointer :: generic
    
    generic=>c_get(all_diagnostics_at_timestep, key)    
    if (associated(generic)) then
      select type(generic)
      type is(all_diagnostics_at_timestep_type)
        get_diagnostic_by_key=>generic
      end select
    else
      get_diagnostic_by_key=>null()
    end if
  end function get_diagnostic_by_key  

  !> Retrieves the all diagnostics at a specific timestep
  !! @param index The index to look up
  !! @param do_lock Whether to issue a read lock or not
  !! @param entry_key is set to be the corresponding key of this diagnostics entry
  !! @returns The corresponding all diagnostics or null if none is found
  function get_diagnostics_at_index(index, do_lock, entry_key)
    integer, intent(in) :: index
    logical, intent(in) :: do_lock
    character(len=*), intent(out) :: entry_key
    type(all_diagnostics_at_timestep_type), pointer :: get_diagnostics_at_index

    class(*), pointer :: generic
    
    if (do_lock) call check_thread_status(forthread_rwlock_rdlock(all_diagnostics_per_timestep_rwlock))
    generic=>c_value_at(all_diagnostics_at_timestep, index)
    entry_key=c_key_at(all_diagnostics_at_timestep, index)
    if (do_lock) call check_thread_status(forthread_rwlock_unlock(all_diagnostics_per_timestep_rwlock))
    if (associated(generic)) then
      select type(generic)
      type is(all_diagnostics_at_timestep_type)
        get_diagnostics_at_index=>generic
      end select
    else
      get_diagnostics_at_index=>null()
    end if
  end function get_diagnostics_at_index  

  !> Retrieves a communication activity from its field name
  !! @param timestep_entry The timestep entry to base the look up upon
  !! @param field_name The field name to look up
  !! @returns The activity or null if none is found
  function get_comm_activity_from_fieldname(diagnostics_by_timestep, field_name)
    type(all_diagnostics_at_timestep_type), intent(inout) :: diagnostics_by_timestep
    character(len=*), intent(in) :: field_name
    type(diagnostics_activity_type), pointer :: get_comm_activity_from_fieldname

    class(*), pointer :: generic
    
    call check_thread_status(forthread_rwlock_rdlock(diagnostics_by_timestep%communication_corresponding_activities_rwlock))
    generic=>c_get(diagnostics_by_timestep%communication_corresponding_activities, field_name)
    call check_thread_status(forthread_rwlock_unlock(diagnostics_by_timestep%communication_corresponding_activities_rwlock))
    if (associated(generic)) then
      select type(generic)
        type is(diagnostics_activity_type)
          get_comm_activity_from_fieldname=>generic
      end select      
    else
      get_comm_activity_from_fieldname=>null()
    end if    
  end function get_comm_activity_from_fieldname  

  !> Retrieves a misc action from the parsed user XML configuration at a specific index
  !! @param action_members The members to extract from
  !! @param index The index to look up
  !! @returns The misc item at this index or null if none is found
  function get_misc_action_at_index(action_members, index)
    type(list_type), intent(inout) :: action_members
    integer, intent(in) :: index
    type(io_configuration_misc_item_type), pointer :: get_misc_action_at_index

    class(*), pointer :: generic

    generic=>c_get(action_members, index)
    if (associated(generic)) then
      select type(generic)
      type is(io_configuration_misc_item_type)
        get_misc_action_at_index=>generic
      end select
    else
      get_misc_action_at_index=>null()
    end if
  end function get_misc_action_at_index

  !> Based upon the IO configuration this will define the diagnostics structure. It is done once at initialisation and then this
  !! same information is used for execution at each data arrival point.
  !! @param io_configuration The IO server configuration
  !! @param diagnostic_generation_frequency Map of diagnostic name to the frequency (in timesteps) of generation
  subroutine define_diagnostics(io_configuration, diagnostic_generation_frequency)
    type(io_configuration_type), intent(inout) :: io_configuration
    type(hashmap_type), intent(out) :: diagnostic_generation_frequency

    integer :: i, j, entries, action_entities, activity_freq
    type(io_configuration_misc_item_type), pointer :: misc_action
    type(diagnostics_activity_type), pointer :: item
    character(len=STRING_LENGTH) :: activity_name, timestep_string

    class(*), pointer :: generic

    entries=io_configuration%number_of_diagnostics
    if (entries .gt. 0) then
      allocate(diagnostic_definitions(entries))

      do i=1, entries
        diagnostic_definitions(i)%generation_timestep_frequency=0
        diagnostic_definitions(i)%diagnostic_name=io_configuration%diagnostics(i)%name
        diagnostic_definitions(i)%collective=io_configuration%diagnostics(i)%collective
        action_entities=c_size(io_configuration%diagnostics(i)%members)
        if (action_entities .gt. 0) then
          do j=1, action_entities
            misc_action=>get_misc_action_at_index(io_configuration%diagnostics(i)%members, j)
            allocate(item)          
            item%result_name=conv_to_string(c_get(misc_action%embellishments, "result"), .false., STRING_LENGTH)
            if (c_contains(misc_action%embellishments, "root")) then
              if (get_action_attribute_string(misc_action%embellishments, "root") .eq. "auto") then
                item%root=mod(i, io_configuration%number_of_io_servers)
              else
                item%root=get_action_attribute_integer(misc_action%embellishments, "root")
              end if
            else
              item%root=-1
            end if
            if (misc_action%type .eq. "operator") then
              activity_name=conv_to_string(c_get(misc_action%embellishments, "name"), .false., STRING_LENGTH)
              item%activity_name=activity_name
              item%required_fields=get_operator_required_fields(activity_name, misc_action%embellishments)
              item%activity_attributes=misc_action%embellishments
              item%operator_procedure=>get_operator_perform_procedure(activity_name)
              item%activity_type=OPERATOR_TYPE          
            else if (misc_action%type .eq. "communication") then
              if (item%root .lt. 0) call log_log(LOG_ERROR, "Root must be supplied and 0 or greater for communication actions")
              activity_name=conv_to_string(c_get(misc_action%embellishments, "name"), .false., STRING_LENGTH)
              if (activity_name .eq. "reduction" .or. activity_name .eq. "allreduction") then            
                call c_add(item%required_fields, c_get(misc_action%embellishments, "field"))
                item%activity_type=merge(REDUCTION_TYPE, ALLREDUCTION_TYPE, activity_name .eq. "reduction")
                item%communication_operator=get_reduction_operator(conv_to_string(c_get(&
                     misc_action%embellishments, "operator"), .false., STRING_LENGTH))
              else if (activity_name .eq. "broadcast") then            
                call c_add(item%required_fields, c_get(misc_action%embellishments, "field"))
                item%activity_type=BROADCAST_TYPE
              end if
            end if
            call add_required_fields_if_needed(item%required_fields)
            activity_freq=get_diagnostic_generation_frequency(io_configuration, item%required_fields)
            if (diagnostic_definitions(i)%generation_timestep_frequency .lt. activity_freq) then
              diagnostic_definitions(i)%generation_timestep_frequency=activity_freq
            end if
            generic=>item
            call c_add(diagnostic_definitions(i)%activities, generic)
          end do
        end if
        call c_put(diagnostic_generation_frequency, diagnostic_definitions(i)%diagnostic_name, &
             conv_to_generic(diagnostic_definitions(i)%generation_timestep_frequency, .true.))
        call process_auto_dimensions(io_configuration, io_configuration%diagnostics(i), i)
      end do
    end if
  end subroutine define_diagnostics

  !> Retrives a diagnostic activity based upon its result name or null if none is found
  !! @param result_name The name of the result we are looking up
  !! @param diagnostic_entry_index The diagnostic index that we are concerned with
  !! @returns The corresponding activity or null if none is found
  function get_diagnostic_activity_by_result_name(result_name, diagnostic_entry_index)
    character(len=STRING_LENGTH), intent(inout) :: result_name
    integer, intent(in) :: diagnostic_entry_index
    type(diagnostics_activity_type), pointer :: get_diagnostic_activity_by_result_name

    integer :: i, entries

    entries=c_size(diagnostic_definitions(diagnostic_entry_index)%activities)
    do i=1, entries
      get_diagnostic_activity_by_result_name=>get_activity_at_index(diagnostic_definitions(diagnostic_entry_index)%activities, i)
      if (get_diagnostic_activity_by_result_name%result_name == result_name) then
        return
      end if      
    end do
    get_diagnostic_activity_by_result_name=>null()
  end function get_diagnostic_activity_by_result_name

  !> Processes all auto dimensions by looking them up and resolving them based upon the operators
  !! @param io_configuration Configuration of the IO server
  !! @param diagnostic_configuration Configuration of the diagnostic field
  !! @param entry_index The specific diagnostic entry that we care about
  subroutine process_auto_dimensions(io_configuration, diagnostic_configuration, entry_index)
    type(io_configuration_type), intent(inout) :: io_configuration
    type(io_configuration_diagnostic_field_type), intent(inout) :: diagnostic_configuration
    integer, intent(in) :: entry_index

    integer :: i, auto_index, number_auto_dims, diag_modified_dim_size
    character(len=STRING_LENGTH) :: specific_dimension
    type(diagnostics_activity_type), pointer :: diagnostic_activity

    diagnostic_activity=>get_diagnostic_activity_by_result_name(diagnostic_definitions(entry_index)%diagnostic_name, entry_index)
    if (associated(diagnostic_activity)) then
      if (diagnostic_activity%activity_type==OPERATOR_TYPE) then
        do i=1, diagnostic_configuration%dimensions
          auto_index=index(diagnostic_configuration%dim_size_defns(i), "-auto")
          if (auto_index .ne. 0) then
            specific_dimension=diagnostic_configuration%dim_size_defns(i)(1:auto_index-1)
            diag_modified_dim_size=get_operator_auto_size(io_configuration, diagnostic_activity%activity_name, &
                 specific_dimension, diagnostic_activity%activity_attributes)
            if (diag_modified_dim_size .ge. 0) then
              specific_dimension=trim(specific_dimension)//"_"//trim(conv_to_string(diag_modified_dim_size))
              diagnostic_configuration%dim_size_defns(i)=specific_dimension
              call c_put(io_configuration%dimension_sizing, specific_dimension, conv_to_generic(diag_modified_dim_size, .true.))
            else
              diagnostic_configuration%dim_size_defns(i)=specific_dimension
            end if
          end if
        end do
      end if
    end if
  end subroutine process_auto_dimensions

  !> Retrieves the max diagnostic generation frequency for a set of fields
  !! @param io_configuration The IO server configuration
  !! @param required_fields List of required fields that we are looking up
  !! @returns The max generation frequency
  integer function get_diagnostic_generation_frequency(io_configuration, required_fields)
    type(io_configuration_type), intent(inout) :: io_configuration
    type(list_type), intent(inout) :: required_fields

    integer :: i, entries, field_freq
    character(len=STRING_LENGTH) :: required_field_name

    entries=c_size(required_fields)
    get_diagnostic_generation_frequency=0
    do i=1, entries
      required_field_name=conv_to_string(c_get(required_fields, i), .false., STRING_LENGTH)
      field_freq=get_field_frequency(io_configuration, required_field_name)
      if (get_diagnostic_generation_frequency .lt. field_freq) get_diagnostic_generation_frequency=field_freq
    end do    
  end function get_diagnostic_generation_frequency  

  !> Retrieves the generation frequency for a specific field
  !! @param io_configuration IO server configuration
  !! @param field_name The name of the field we are looking up
  !! @returns The generation frequency for this specific field
  integer function get_field_frequency(io_configuration, field_name)
    type(io_configuration_type), intent(inout) :: io_configuration
    character(len=*), intent(in) :: field_name

    integer :: i, j
    do i=1, io_configuration%number_of_data_definitions
      do j=1, io_configuration%data_definitions(i)%number_of_data_fields
        if (io_configuration%data_definitions(i)%fields(j)%name == field_name) then
          get_field_frequency=io_configuration%data_definitions(i)%frequency
          return
        end if
      end do
    end do
    get_field_frequency=0
  end function get_field_frequency
end module diagnostic_federator_mod
