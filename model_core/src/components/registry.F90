!> MONC component registry
!!
!! Supports management of components. Each stage is called via the registry which will execute the
!! installed callback hooks for that stage in order.
module registry_mod
  use datadefn_mod, only : STRING_LENGTH
  use collections_mod, only : list_type, hashmap_type, map_type, iterator_type, c_size, c_generic_at, c_key_at, c_get_integer, &
       c_get_string, c_get_generic, c_remove, c_put_generic, c_put_string, c_put_integer, c_put_real, c_is_empty, &
       c_contains, c_add_generic, c_add_string, c_free, c_get_iterator, c_has_next, c_next_mapentry, mapentry_type
  use monc_component_mod, only : component_descriptor_type, component_field_value_type, component_field_information_type, &
       FINALISATION_PRIORITY_INDEX, INIT_PRIORITY_INDEX, TIMESTEP_PRIORITY_INDEX, &
       pointer_wrapper_value_type, pointer_wrapper_info_type
  use conversions_mod, only : conv_to_string
  use state_mod, only : model_state_type
  use optionsdatabase_mod, only : options_has_key, options_get_string, options_get_logical, options_get_array_size
  use logging_mod, only : LOG_INFO, LOG_ERROR, LOG_WARN, log_master_log, log_is_master
  use grids_mod, only : X_INDEX, Y_INDEX, Z_INDEX

  implicit none

#ifndef TEST_MODE
  private
#endif

  integer, parameter :: GROUP_TYPE_WHOLE=0, & !< Execute the callbacks in this group once per timestep
       GROUP_TYPE_COLUMN=1 !< Execute the callbacks in this group for each column per timestep

  type(list_type), save :: field_information

  integer, dimension(:), allocatable :: group_types !< Group types
  character(len=STRING_LENGTH), dimension(:), allocatable :: enabled_component_input_keys,& !< Temporary read array of component enable names
       group_locations !< Provides an id to each group

  !> Private helper type which wraps a procedure pointer. This is needed for storage in our collections_mod
  !! and allows us to associate additional information, such as number of times called, performance etc
  !! in future if we wanted to
  type pointer_wrapper_type
    procedure(), nopass, pointer :: ptr !< Procedure pointer which is used to point to the callback hook
  end type pointer_wrapper_type

  !> Descriptor of a group
  type group_descriptor_type
     character(len=STRING_LENGTH) :: name !< Name of the group
     character(len=STRING_LENGTH), dimension(30) :: group_members
     integer :: type, & !< Type (execute once per timestep or per column)
          id, & !< Id number of the group which is also the order executed
          number_of_members
  end type group_descriptor_type  

  type(map_type), save :: init_callbacks,&     !< Callback hooks for the initialisation stage       
                    finalisation_callbacks,&   !< Callback hooks for the finalisation stage
                    component_descriptions,&   !< Copies of component descriptors
                    group_descriptors,&        !< Group descriptors for each group, name->descriptor     
                    component_groups, init_orderings, finalisation_orderings
  type(hashmap_type), save :: field_procedure_retrievals, field_procedure_sizings

  type(map_type), dimension(:), allocatable :: timestep_callbacks !< Callback hooks for the timestep stage

  public GROUP_TYPE_WHOLE, GROUP_TYPE_COLUMN, group_descriptor_type, register_component, deregister_component, &
       execute_initialisation_callbacks, execute_timestep_callbacks, &
       execute_finalisation_callbacks, get_component_info, get_all_registered_components, free_registry, init_registry, &
       order_all_callbacks, display_callbacks_in_order_at_each_stage, get_ordered_groups, &
       is_component_enabled, get_all_component_published_fields, get_component_field_value, get_component_field_information, &
       is_component_field_available
contains

  !> Initialises the registry with the provided configuration file
  !! @param configurationFileName The filename of the configuration file to parse
  subroutine init_registry(options_database)
    type(hashmap_type), intent(inout) :: options_database

    call read_group_configurations(options_database)
    call read_initialisation_and_finalisation_orders(options_database)
    allocate(timestep_callbacks(c_size(group_descriptors)))
  end subroutine init_registry

  !> Will deregister all components and free up the registry data structures. This can either be called
  !! at the end of execution to clean memory up or used to clear the registry
  subroutine free_registry()
    integer :: i, entries
   
    entries = c_size(component_descriptions)
    do i=1, entries
      ! Key de registering key at element one as each deregister removes elements from map_type
      call deregister_component(c_key_at(component_descriptions, 1))
    end do

    call c_free(init_callbacks)
    do i=1,size(timestep_callbacks)
      call c_free(timestep_callbacks(i))
    end do
    deallocate(timestep_callbacks)
    call c_free(finalisation_callbacks)
    call c_free(component_descriptions)
    call c_free(component_groups)
    call c_free(init_orderings)
    call c_free(finalisation_orderings)
  end subroutine free_registry

  !> Will register a component and install the nescesary callback hooks
  !! @param descriptor The component descriptor and a separate copy of this it stored as reference
  subroutine register_component(options_database, descriptor)
    type(hashmap_type), intent(inout) :: options_database
    type(component_descriptor_type), intent(in) :: descriptor

    type(component_descriptor_type), pointer :: registry_descriptor
    class(*), pointer :: description_data
    character(len=STRING_LENGTH) :: group_name
    logical :: component_enabled

    allocate(registry_descriptor, source=descriptor) ! Make copy of the descriptor (which registry might modify)

    if (options_has_key(options_database, trim(descriptor%name)//"_enabled")) then
      component_enabled=options_get_logical(options_database, trim(descriptor%name)//"_enabled")
    else
      component_enabled=.false.
      call log_master_log(LOG_WARN, "No enabled configuration for component "//trim(descriptor%name)//" therefore disabling this")
    end if

    if (component_enabled) then
      if (c_contains(component_groups, trim(descriptor%name))) then
        group_name=c_get_string(component_groups, trim(descriptor%name))
        call load_callback_hooks(registry_descriptor, group_name)
      else
        call load_callback_hooks(registry_descriptor)
      end if
      call load_published_fields(descriptor)
    end if    

    description_data => registry_descriptor
    call c_put_generic(component_descriptions, descriptor%name, description_data, .false.)
  end subroutine register_component

  !> Determines whether a specific published field is available or not
  !! @param name The name of the field to search for
  !! @returns Whether this field is found
  logical function is_component_field_available(name)
    character(len=*), intent(in) :: name

    is_component_field_available=c_contains(field_procedure_retrievals, name)
  end function is_component_field_available  

  !> Retrieves the value wrapper of a components published field
  !! @param name The name of the field to look up
  !! @returns The value wrapper that is associated with this field
  function get_component_field_value(current_state, name)
    type(model_state_type), target, intent(inout) :: current_state
    character(len=*), intent(in) :: name
    type(component_field_value_type) :: get_component_field_value

    class(*), pointer :: data

    if (c_contains(field_procedure_retrievals, name)) then
      data=>c_get_generic(field_procedure_retrievals, name)
      select type(data)
      type is (pointer_wrapper_value_type)
        call data%ptr(current_state, name, get_component_field_value)
      end select
    else
      call log_master_log(LOG_ERROR, "Published field '"//trim(name)//"' is not found in any enabled components")
    end if
  end function get_component_field_value
  
  !> Retrieves information about a components published field which includes its type and size
  !! @param name The name of the field to look up
  !! @returns The value wrapper that is associated with this field
  function get_component_field_information(current_state, name)
    type(model_state_type), target, intent(inout) :: current_state
    character(len=*), intent(in) :: name
    type(component_field_information_type) :: get_component_field_information

    class(*), pointer :: data

    if (c_contains(field_procedure_sizings, name)) then
      data=>c_get_generic(field_procedure_sizings, name)
      select type(data)
      type is (pointer_wrapper_info_type)
        call data%ptr(current_state, name, get_component_field_information)
      end select
    else
      call log_master_log(LOG_ERROR, "Published field '"//trim(name)//"' is not found in any enabled components")
    end if
  end function get_component_field_information  

  !> Retrieves all of the published field information
  !! @returns The list of published fields
  function get_all_component_published_fields()
    type(list_type) :: get_all_component_published_fields

    get_all_component_published_fields=field_information
  end function get_all_component_published_fields  

  !> Loads the published fields information for an entire component into the registry's definition list
  !! @param descriptor The field descriptor to load in
  subroutine load_published_fields(descriptor)
    type(component_descriptor_type), intent(in) :: descriptor

    integer :: i
    class(*), pointer :: field_generic_description
    type(pointer_wrapper_info_type), pointer :: wrapper_info
    type(pointer_wrapper_value_type), pointer :: wrapper_value
    class(*), pointer :: genericwrapper

    if (associated(descriptor%published_fields) .and. associated(descriptor%field_value_retrieval) .and. &
         associated(descriptor%field_information_retrieval)) then
      do i=1, size(descriptor%published_fields)
        field_generic_description=>descriptor%published_fields(i)
        call c_add_generic(field_information, field_generic_description, .false.)
        
        allocate(wrapper_value) ! We allocate our own copy of the descriptor here to ensure the consistency of registry information
        wrapper_value%ptr => descriptor%field_value_retrieval
        genericwrapper=>wrapper_value
        call c_put_generic(field_procedure_retrievals, descriptor%published_fields(i), genericwrapper, .false.)
        
        allocate(wrapper_info)
        wrapper_info%ptr => descriptor%field_information_retrieval
        genericwrapper=>wrapper_info
        call c_put_generic(field_procedure_sizings, descriptor%published_fields(i), genericwrapper, .false.)
      end do

      else if (associated(descriptor%published_fields) .or. associated(descriptor%field_value_retrieval) .or. &
         associated(descriptor%field_information_retrieval)) then
        call log_master_log(LOG_WARN, "Component "//trim(descriptor%name)//&
             " has provided incomplete configuration for published fields, therefore ignoring these")
    end if
  end subroutine load_published_fields

  !> Will deregister a component, remove all callback hooks and free registry specific memory allocated
  !! to the component
  !! @param name The name of the component to de-register
  subroutine deregister_component(name)
    character(len=*), intent(in) :: name
    class(*), pointer :: description_data

    description_data=>c_get_generic(component_descriptions, name)
    select type(description_data)
      type is (component_descriptor_type)
        call remove_descriptor(description_data)
    end select
    deallocate(description_data)  ! Is added as a clone of the original so free the memory
  end subroutine deregister_component

  !> Retrieves detailed information about a specific component
  !! @param name The name of the component to retrieve information for
  !! @returns The registry's copy of the component descriptor or null if it does not exist
  function get_component_info(name)
    character(len=*), intent(in) :: name
    type(component_descriptor_type), pointer :: get_component_info
    class(*), pointer :: description_data

    get_component_info => null()
    description_data => c_get_generic(component_descriptions, name)
    if (associated(description_data)) then
      select type(description_data)
        type is (component_descriptor_type)
        get_component_info => description_data
      end select
    end if
  end function get_component_info

  !> Returns a brief summary of all registered components
  !! @returns A map_type where the keys are the component names and value the corresponding version number
  type(map_type) function get_all_registered_components()
    integer :: i, number_of_components
    class(*), pointer :: description_data
    type(iterator_type) :: iterator

    iterator=c_get_iterator(component_descriptions)
    do while (c_has_next(iterator))
      description_data=>c_get_generic(c_next_mapentry(iterator))
      select type(description_data)
        type is(component_descriptor_type)
          call c_put_real(get_all_registered_components, description_data%name, description_data%version)
      end select
    end do    
  end function get_all_registered_components

  !> Calls all initialisation callbacks with the specified state
  !! @param currentState The current model state which may (and often is) modified
  subroutine execute_initialisation_callbacks(current_state)
    type(model_state_type), intent(inout) :: current_state

    call execute_callbacks(init_callbacks, current_state)
  end subroutine execute_initialisation_callbacks

  !> Calls all timestep callbacks with the specified state
  !! @param currentState The current model state which may (and often is) modified
  subroutine execute_timestep_callbacks(current_state, group_id)
    type(model_state_type), intent(inout) :: current_state
    integer :: group_id

    if (.not. c_is_empty(timestep_callbacks(group_id))) then
      call execute_callbacks(timestep_callbacks(group_id), current_state)
    end if
  end subroutine execute_timestep_callbacks

  !> Calls all finalisation callbacks with the specified state
  !! @param currentState The current model state which may (and often is) modified
  subroutine execute_finalisation_callbacks(current_state)
    type(model_state_type), intent(inout) :: current_state

    call execute_callbacks(finalisation_callbacks, current_state)
  end subroutine execute_finalisation_callbacks

  !> Orders all callbacks in the prospective stages based upon the priorities of each descriptor.
  !!
  !! Note that this is an expensive operation as it involves multiple searches of the component
  !! list_types so should be called sparingly. When two priorities are equal then the order depends
  !! upon which one was registered first.
  subroutine order_all_callbacks()
    integer :: i
    type(map_type) :: specific_ts

    call rebalance_callbacks(init_callbacks, init_orderings, "initialisation")
    do i=1,size(timestep_callbacks)
      specific_ts=order_grouped_timstep_callbacks(i)
      call rebalance_callbacks(timestep_callbacks(i), specific_ts, "time stepping")
      call c_free(specific_ts)
    end do    
    call rebalance_callbacks(finalisation_callbacks, finalisation_orderings, "finalisation")
  end subroutine order_all_callbacks

  type(map_type) function order_grouped_timstep_callbacks(group_id)
    integer, intent(in) :: group_id
    
    integer :: group_size, i

    type(group_descriptor_type) :: descriptor

    descriptor=get_group_descriptor_from_id(group_id)
    group_size=descriptor%number_of_members
    do i=1, group_size
      call c_put_integer(order_grouped_timstep_callbacks, trim(descriptor%group_members(i)), i)
    end do    
  end function order_grouped_timstep_callbacks

  !> Determines whether or not a specific component is registered and enabled
  !! @param component_name The name of the component to check for
  logical function is_component_enabled(options_database, component_name)
    type(hashmap_type), intent(inout) :: options_database
    character(len=*), intent(in) :: component_name
    
    if (c_contains(component_descriptions, component_name)) then
      if (options_has_key(options_database, trim(component_name)//"_enabled")) then
        is_component_enabled=options_get_logical(options_database, trim(component_name)//"_enabled")
        return
      end if
    end if
    is_component_enabled=.false.
  end function is_component_enabled  

  !> Displays the registered callbacks of each stage in the order that they will be called
  subroutine display_callbacks_in_order_at_each_stage()
    integer :: i
    type(group_descriptor_type) :: group_descriptor

    call display_callbacks_in_order(init_callbacks, "init")
    do i=1,size(timestep_callbacks)
      group_descriptor=get_group_descriptor_from_id(i)
      call display_callbacks_in_order(timestep_callbacks(i), "timestep group '"//trim(group_descriptor%name)//"'")
    end do    
    call display_callbacks_in_order(finalisation_callbacks, "finalisation")
  end subroutine display_callbacks_in_order_at_each_stage

  !> Orders all the groups (in the order that they will be called in) and returns an array with these in order. 
  !! This is useful for prefetching the groups in order so that per timestep we can just iterate through the array which is O(n)
  !! @returns An array of groups in the order that they will be executed in
  subroutine get_ordered_groups(ordered_groups)
    type(group_descriptor_type), dimension(:), allocatable :: ordered_groups

    integer :: i

    allocate(ordered_groups(c_size(group_descriptors)))
    do i=1,c_size(group_descriptors)
      ordered_groups(i)=get_group_descriptor_from_id(i)
    end do    
  end subroutine get_ordered_groups

  !--------------------------------------------------------------------------
  ! Private procedures acting as helpers to the functionality of the registry
  !--------------------------------------------------------------------------

  !> Given a group name this returns the id (i.e. order) of that group
  !! @param group_name The group name to look up
  !! @returns The id (also order) of the group
  integer function get_group_id(group_name)
    character(len=*), intent(in) :: group_name

    type(group_descriptor_type) :: descriptor

    descriptor=get_group_descriptor_from_name(group_name)
    get_group_id=descriptor%id
  end function get_group_id  

  !> Given a group name this returns the group descriptor corresponding to that or an error if none is found
  !! @param group_name Name of the group to look up
  !! @returns Descriptor of the corresponding group
  type(group_descriptor_type) function get_group_descriptor_from_name(group_name)
    character(len=*), intent(in) :: group_name

    class(*), pointer :: generic       

    generic=>c_get_generic(group_descriptors, group_name)
    if (associated(generic)) then
      select type(generic)
        type is(group_descriptor_type)
          get_group_descriptor_from_name=generic
      end select      
    else
      call log_master_log(LOG_ERROR, "No configuration specified for group "//group_name)
    end if    
  end function get_group_descriptor_from_name

  !> Given the id of a group this will return the corresponding descriptor
  !! @param group_id Id of the group to find
  !! @returns The descriptor that has the correct id number
  type(group_descriptor_type) function get_group_descriptor_from_id(group_id)
    integer, intent(in) :: group_id

    type(iterator_type) :: iterator
    class(*), pointer :: generic

    iterator=c_get_iterator(group_descriptors)
    do while (c_has_next(iterator))
      generic=>c_get_generic(c_next_mapentry(iterator))
       if (associated(generic)) then
        select type(generic)
        type is(group_descriptor_type)
          if (generic%id == group_id) then
            get_group_descriptor_from_id=generic
            return
          end if
        end select
      end if
    end do     
  end function get_group_descriptor_from_id
  
  !> Displays the registered callbacks of a specific stage in the order that they will be called
  !! @param stageCallbacks The registered callbacks for a stage
  !! @param stagetitle The title of the stage - used for printing out information
  subroutine display_callbacks_in_order(stage_callbacks, stagetitle)
    type(map_type), intent(inout) :: stage_callbacks
    character(len=*), intent(in) :: stagetitle

    integer :: i, entries

    entries = c_size(stage_callbacks)
    do i=1,entries
      call log_master_log(LOG_INFO, "Stage: "//stagetitle//" at: "//trim(conv_to_string(i))//&
                          "  "//trim(c_key_at(stage_callbacks, i)) )
    end do
  end subroutine display_callbacks_in_order

  subroutine read_initialisation_and_finalisation_orders(options_database)
    type(hashmap_type), intent(inout) :: options_database

    call read_specific_orders(options_database, "initialisation_stage_ordering", init_orderings)
    call read_specific_orders(options_database, "finalisation_stage_ordering", finalisation_orderings)
  end subroutine read_initialisation_and_finalisation_orders

  subroutine read_specific_orders(options_database, key, data_structure)
    type(hashmap_type), intent(inout) :: options_database
    type(map_type), intent(inout) :: data_structure
    character(len=*) :: key

    integer :: number_of_elements, i
    character(len=STRING_LENGTH) :: component_name

    number_of_elements=options_get_array_size(options_database, key)
    do i=1, number_of_elements
      component_name=options_get_string(options_database, key, i)
      call c_put_integer(data_structure, trim(component_name), i)
    end do
  end subroutine read_specific_orders  

  subroutine read_group_configurations(options_database)
    type(hashmap_type), intent(inout) :: options_database

    integer :: group_elements, i, j
    class(*), pointer :: generic_to_add
    type(group_descriptor_type) :: group_description
    character(len=STRING_LENGTH) :: group_type

    group_elements=options_get_array_size(options_database, "group_names")
    if (group_elements .lt. 1) call log_master_log(LOG_ERROR, "You must provide some group definitions")
    do i=1, group_elements
      group_description%name=trim(options_get_string(options_database, "group_names", i))
      if (options_has_key(options_database, trim(group_description%name)//"_group_type")) then
        group_type=trim(options_get_string(options_database, trim(group_description%name)//"_group_type"))
        if (trim(group_type) .eq. "entire") then
          group_description%type=0
        else if (trim(group_type) .eq. "column") then
          group_description%type=1
        else if (trim(group_type) .eq. "slice") then
          group_description%type=2
        else
          call log_master_log(LOG_ERROR, "Group type "//trim(group_type)//" for group "&
               //trim(group_description%name)//" not understood")
        end if
      else
        call log_master_log(LOG_ERROR, "No group type for group "//trim(group_description%name))
      end if
      group_description%id=i
      if (.not. options_has_key(options_database, trim(group_description%name)//"_group_contents")) then
        call log_master_log(LOG_ERROR, "No component contents specified for group "//trim(group_description%name))
      end if
      group_description%number_of_members=options_get_array_size(options_database, &
           trim(group_description%name)//"_group_contents")    
      if (group_description%number_of_members == 0 .or. .not. options_has_key(options_database, &
           trim(group_description%name)//"_group_contentsa_size")) then
        if (options_has_key(options_database, trim(group_description%name)//"_group_contents")) then
          group_description%number_of_members=1
          group_description%group_members(1)=trim(options_get_string(options_database, &
                 trim(group_description%name)//"_group_contents"))
          call c_put_string(component_groups, trim(group_description%group_members(1)), group_description%name)
        else
          call log_master_log(LOG_ERROR, "No contents specified for group "//trim(group_description%name))
        end if
      else
        do j=1, group_description%number_of_members
          group_description%group_members(j)=trim(options_get_string(options_database, &
               trim(group_description%name)//"_group_contents", j))
          call c_put_string(component_groups, trim(group_description%group_members(j)), group_description%name)
        end do
      end if
      allocate(generic_to_add, source=group_description)
      call c_put_generic(group_descriptors, group_description%name, generic_to_add, .false.)
    end do    
  end subroutine read_group_configurations  

  !> Will remove a specific descriptor from the registry table and uninstall the corresponding
  !! callback hooks for each state
  !! @param descriptor The component descriptor which is to be removed
  subroutine remove_descriptor(descriptor)
    type(component_descriptor_type), intent(inout) :: descriptor

    character(len=STRING_LENGTH) :: group_name

    if (c_contains(component_groups, trim(descriptor%name))) then
      group_name=c_get_string(component_groups, trim(descriptor%name))
      call unload_callback_hooks(descriptor, group_name)
    else
      call unload_callback_hooks(descriptor)
    end if
    call c_remove(component_descriptions, descriptor%name)
  end subroutine remove_descriptor

  !> Will unload the callback hooks that have been installed for each state
  !! @param descriptor The component descriptor which is to be unloaded
  subroutine unload_callback_hooks(descriptor, group_name)
    type(component_descriptor_type), intent(in) :: descriptor
    character(len=*), intent(in), optional :: group_name

    if (associated(descriptor%initialisation)) call c_remove(init_callbacks, descriptor%name)
    if (associated(descriptor%timestep) .and. present(group_name)) &
         call c_remove(timestep_callbacks(get_group_id(group_name)), descriptor%name)    
    if (associated(descriptor%finalisation)) call c_remove(finalisation_callbacks, descriptor%name)
  end subroutine unload_callback_hooks

  !> Will install the callback hooks for each state
  !! @param descriptor The component descriptor which is to be installed
  subroutine load_callback_hooks(descriptor, group_name)
    type(component_descriptor_type), intent(in) :: descriptor
    character(len=*), intent(in), optional :: group_name

    if (associated(descriptor%initialisation)) call add_callback(init_callbacks, descriptor%name, descriptor%initialisation)
    if (associated(descriptor%timestep)) then
      if (present(group_name)) then
        call add_callback(timestep_callbacks(get_group_id(group_name)), descriptor%name, descriptor%timestep)
      else
        call log_master_log(LOG_ERROR, "In the configuration you must provide a group for component "&
             //trim(descriptor%name)//" which has a timestep callback")
      end if
    end if
    if (associated(descriptor%finalisation)) call add_callback(finalisation_callbacks, descriptor%name, descriptor%finalisation)
  end subroutine load_callback_hooks

  subroutine rebalance_callbacks(callbacks, priorities, stage_name)
    type(map_type), intent(inout) :: callbacks, priorities
    character(len=*), intent(in) :: stage_name

    type(map_type) :: ordered_callbacks
    integer :: i, entries_in_list, current_item
    class(*), pointer :: generic
    
    entries_in_list=c_size(callbacks)
    do i=1, entries_in_list
      current_item=get_highest_callback_priority(callbacks, priorities)
      if (current_item .ge. 1) then
        generic=>c_generic_at(callbacks, current_item)
        call c_put_generic(ordered_callbacks, c_key_at(callbacks, current_item), generic, .false.)
        call c_remove(callbacks, c_key_at(callbacks, current_item))
      end if
    end do

    entries_in_list=c_size(callbacks)
    do i=1, entries_in_list
      generic=>c_generic_at(callbacks, i)
      call c_put_generic(ordered_callbacks, c_key_at(callbacks, i), generic, .false.)
      call log_master_log(LOG_WARN, "Run order callback for component "//trim(c_key_at(callbacks, i))//&
           " at stage "//stage_name//" not specified")
    end do
    callbacks=ordered_callbacks
  end subroutine rebalance_callbacks  

  integer function get_highest_callback_priority(callbacks, priorities)
    type(map_type), intent(inout) :: callbacks, priorities

    integer :: i, entries_in_list, min_priority, min_location, priority
    character(len=STRING_LENGTH) :: key

    min_location=0
    min_priority=1000

    entries_in_list=c_size(callbacks)
    do i=1, entries_in_list
      key=c_key_at(callbacks, i)
      if (c_contains(priorities, key)) then
        priority=c_get_integer(priorities, key)
        if (priority .lt. min_priority) then
          min_priority=priority
          min_location=i
        end if      
      end if
    end do
    get_highest_callback_priority=min_location
  end function get_highest_callback_priority

  !> Will execute the appropriate callbacks in a specific map_type given the current state
  !! @param callbackmap_type The map_type of callback hooks to execute
  !! @param currentState The model state which may be (and likely is) modified in callbacks
  subroutine execute_callbacks(callback_map, current_state)
    type(map_type), intent(inout) :: callback_map
    type(model_state_type), intent(inout) :: current_state

    class(*), pointer :: data
    type(iterator_type) :: iterator
    type(mapentry_type) :: map_entry
    integer :: k,j,i

    iterator=c_get_iterator(callback_map)
    do while (c_has_next(iterator))
      map_entry=c_next_mapentry(iterator)
      data=>c_get_generic(map_entry)
      select type(data)
      type is (pointer_wrapper_type)
      call data%ptr(current_state)

      if (current_state%print_debug_data) then
        if (log_is_master()) then
          k=current_state%local_grid%size(Z_INDEX)/2
          j=current_state%local_grid%local_domain_start_index(Y_INDEX)
          i=current_state%local_grid%local_domain_start_index(X_INDEX)
          if (allocated(current_state%u%data) .and. allocated(current_state%sth%data) &
              .and. allocated(current_state%zu%data) .and. allocated(current_state%sw%data) &
              .and. current_state%column_local_x == i .and. current_state%column_local_y == j) then
            print *, trim(map_entry%key),' ', k,j,i, &
              current_state%zu%data(k,j,i), current_state%u%data(k,j,i), &
              current_state%sth%data(k,j,i), current_state%sw%data(k,j,i)
          end if
        end if
      end if

!        type is (pointer_wrapper_init_type)
!          call data%ptr(current_state)
!        type is (pointer_wrapper_timestep_type)
!          call data%ptr(current_state)
!        type is (pointer_wrapper_finalisation_type)
!          call data%ptr(current_state)
      end select
    end do
  end subroutine execute_callbacks

  !> Will install a specific callback hook into the specified map_type of existing hooks
  !! @param callbackmap_type The map_type of existing callbacks which we are going to install this one into
  !! @param name The name of the callback that we are installing
  !! @param procedurePointer Pointer to the procedure which implements the callback
  subroutine add_callback(callback_map, name, procedure_pointer)
    type(map_type), intent(inout) :: callback_map
    procedure(), pointer :: procedure_pointer
    character(len=*), intent(in) :: name

    type(pointer_wrapper_type), pointer :: wrapper
    class(*), pointer :: genericwrapper

    allocate(wrapper) ! We allocate our own copy of the descriptor here to ensure the consistency of registry information
    wrapper%ptr => procedure_pointer
    genericwrapper=>wrapper
    call c_put_generic(callback_map, name, genericwrapper, .false.)
  end subroutine add_callback
  
!  subroutine add_callback_init(callback_map, name, procedure_pointer)
!    type(map_type), intent(inout) :: callback_map
!    procedure(component_initialisation), pointer :: procedure_pointer
!    character(len=*), intent(in) :: name

!    type(pointer_wrapper_init_type), pointer :: wrapper
!    class(*), pointer :: genericwrapper

!    allocate(wrapper) ! We allocate our own copy of the descriptor here to ensure the consistency of registry information
!    wrapper%ptr => procedure_pointer
!    genericwrapper=>wrapper
!    call c_put_generic(callback_map, name, genericwrapper, .false.)
!  end subroutine add_callback_init
  
!  subroutine add_callback_timestep(callback_map, name, procedure_pointer)
!    type(map_type), intent(inout) :: callback_map
!    procedure(component_timestep), pointer :: procedure_pointer
!    character(len=*), intent(in) :: name

!    type(pointer_wrapper_timestep_type), pointer :: wrapper
!    class(*), pointer :: genericwrapper

!    allocate(wrapper)
!    wrapper%ptr => procedure_pointer
!    genericwrapper=>wrapper
!    call c_put_generic(callback_map, name, genericwrapper, .false.)
!  end subroutine add_callback_timestep
  
!  subroutine add_callback_finalisation(callback_map, name, procedure_pointer)
!    type(map_type), intent(inout) :: callback_map
!    procedure(component_finalisation), pointer :: procedure_pointer
!    character(len=*), intent(in) :: name

!    type(pointer_wrapper_finalisation_type), pointer :: wrapper
!    class(*), pointer :: genericwrapper

!    allocate(wrapper)
!    wrapper%ptr => procedure_pointer
!    genericwrapper=>wrapper
!    call c_put_generic(callback_map, name, genericwrapper, .false.)
!  end subroutine add_callback_finalisation
  
end module registry_mod
