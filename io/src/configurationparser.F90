!> Parses the XML configuration file to produce the io configuration description which contains the data layout
!! specification and rules for handling received data
module configuration_parser_mod
  use datadefn_mod, only : DEFAULT_PRECISION, STRING_LENGTH
  use sax_xml_parser_mod, only : xml_parse
  use conversions_mod, only : conv_to_string, conv_to_integer, conv_to_real, conv_to_uppercase
  use collections_mod, only : hashmap_type, hashset_type, map_type, list_type, mapentry_type, c_get_generic, c_get_integer, &
       c_free, c_size, c_put_integer, c_put_string, c_add_generic, c_add_string
  use logging_mod, only : LOG_WARN, LOG_INFO, LOG_ERROR, log_log, log_master_log, log_master_newline
  use optionsdatabase_mod, only : options_has_key, options_get_logical, options_get_integer, options_get_string, options_get_real, &
                                  options_get_array_size
  use io_server_client_mod, only : ARRAY_FIELD_TYPE, SCALAR_FIELD_TYPE, MAP_FIELD_TYPE, INTEGER_DATA_TYPE, BOOLEAN_DATA_TYPE, &
       STRING_DATA_TYPE, FLOAT_DATA_TYPE, DOUBLE_DATA_TYPE, definition_description_type, field_description_type
  use q_indices_mod, only : get_number_active_q_indices
  use netcdf, only : NF90_DOUBLE, NF90_REAL
  implicit none

#ifndef TEST_MODE
  private
#endif

  character(len=STRING_LENGTH), parameter :: DEFAULT_FILE_TITLE = "MONC diagnostics"
  integer, parameter :: EQ_OPERATOR_TYPE=1, LT_OPERATOR_TYPE=2, GT_OPERATOR_TYPE=3, LTE_OPERATOR_TYPE=4, &
       GTE_OPERATOR_TYPE=5, ADD_OPERATOR_TYPE=6, SUBTRACT_OPERATOR_TYPE=7, MULTIPLY_OPERATOR_TYPE=8, DIV_OPERATOR_TYPE=9, &
       MOD_OPERATOR_TYPE=10, AND_OPERATOR_TYPE=11, OR_OPERATOR_TYPE=12
  integer, parameter :: MONC_SIZE_STRIDE=100, DATA_SIZE_STRIDE=10

  integer, parameter :: TIME_AVERAGED_TYPE=1, INSTANTANEOUS_TYPE=2, NONE_TYPE=3, GROUP_TYPE=1, FIELD_TYPE=2, IO_STATE_TYPE=3

  interface get_data_value_by_field_name
     module procedure get_data_value_from_hashmap_by_field_name, get_data_value_from_map_by_field_name
  end interface get_data_value_by_field_name  

  type data_values_type
     integer :: data_type, dimensions, dim_sizes(4)
     real(kind=DEFAULT_PRECISION), dimension(:), allocatable :: values
     character(len=STRING_LENGTH), dimension(:), allocatable :: string_values
     type(map_type) :: map_values
  end type data_values_type  

  !> Configuration that representes the state of a registered MONC process
  type io_configuration_registered_monc_type
     type(map_type) :: registered_monc_types, registered_monc_buffer_sizes
     type(map_type), dimension(:), allocatable :: field_start_locations, field_end_locations, dimensions
     character(len=STRING_LENGTH), dimension(:), allocatable :: definition_names
     integer :: active_threads, active_mutex, deactivate_condition_variable, local_dim_sizes(3), local_dim_starts(3), &
          local_dim_ends(3), source_id
  end type io_configuration_registered_monc_type  

  !> Configuration associated with the representation of a specific data field
  type io_configuration_field_type
     character(len=STRING_LENGTH) :: name, namespace, dim_size_defns(4), units
     integer :: field_type, data_type, dimensions
     logical :: optional, collective
  end type io_configuration_field_type

  !> Configuration of a specific data definition
  type io_configuration_data_definition_type
     character(len=STRING_LENGTH) :: name, namespace
     logical :: send_on_terminate
     integer :: number_of_data_fields, frequency
     type(map_type) :: compiled_fields, trigger_field_types
     type(io_configuration_field_type), dimension(:), allocatable :: fields
  end type io_configuration_data_definition_type

  type io_configuration_inter_communication_description
     character(len=STRING_LENGTH) :: name
     integer :: message_tag
     procedure(handle_recv_data_from_io_server), pointer, nopass :: handling_procedure
  end type io_configuration_inter_communication_description

  type io_configuration_misc_item_type
     character(len=STRING_LENGTH) :: type, namespace
     type(map_type) :: embellishments
  end type io_configuration_misc_item_type

  type io_configuration_diagnostic_field_type
     character(len=STRING_LENGTH) :: name, dim_size_defns(4), units, namespace
     integer :: field_type, data_type, dimensions
     logical :: collective
     type(list_type) :: members
  end type io_configuration_diagnostic_field_type

  type io_configuration_group_type
     character(len=STRING_LENGTH) :: name, namespace
     type(list_type) :: members
  end type io_configuration_group_type

  type io_configuration_file_writer_facet_type
     integer :: facet_type, time_manipulation_type
     real :: output_time_frequency
     character(len=STRING_LENGTH) :: facet_name     
  end type io_configuration_file_writer_facet_type  

  type io_configuration_file_writer_type
     character(len=STRING_LENGTH) :: file_name, title
     integer :: number_of_contents, write_timestep_frequency, write_precision
     real :: write_time_frequency
     logical :: write_on_model_time, write_on_terminate, include_in_io_state_write
     type(io_configuration_file_writer_facet_type), dimension(:), allocatable :: contents     
  end type io_configuration_file_writer_type

  !> Overall IO configuration
  type io_configuration_type
     integer :: number_of_data_definitions, number_of_diagnostics, io_communicator, number_of_moncs, &
          number_of_io_servers, my_io_rank, active_moncs, number_inter_io_communications, number_of_threads, number_of_groups, &
          number_of_writers, number_of_distinct_data_fields, number_of_global_moncs, general_info_mutex, my_global_rank
     type(io_configuration_data_definition_type), dimension(:), allocatable :: data_definitions
     type(io_configuration_diagnostic_field_type), dimension(:), allocatable :: diagnostics
     type(io_configuration_group_type), dimension(:), allocatable :: groups
     type(io_configuration_file_writer_type), dimension(:), allocatable :: file_writers
     type(io_configuration_registered_monc_type), dimension(:), allocatable :: registered_moncs
     type(io_configuration_inter_communication_description), dimension(:), allocatable :: inter_io_communications
     type(map_type) :: monc_to_index, dimension_sizing
     type(hashmap_type) :: options_database
     real(kind=DEFAULT_PRECISION), dimension(:), allocatable :: zn_field
     real(kind=DEFAULT_PRECISION), dimension(:), allocatable :: z_field
     type(list_type) :: q_field_names
     type(list_type) :: tracer_names
     logical :: general_info_set
     character, dimension(:), allocatable :: text_configuration
  end type io_configuration_type

  abstract interface
     subroutine handle_recv_data_from_io_server(io_configuration, data_buffer, inter_io_index)
       import io_configuration_type
       type(io_configuration_type), intent(inout) :: io_configuration
       character, dimension(:), intent(inout) :: data_buffer
       integer, intent(in) :: inter_io_index
     end subroutine handle_recv_data_from_io_server
  end interface

  !< For reading the IO XML configuration, these are string length constants which can be increased if required
  integer, parameter :: FILE_STR_STRIDE=10000, FILE_LINE_LEN=2000

  logical :: inside_data_definition, inside_handling_definition, inside_server_config, inside_action_config, &
       inside_diagnostic_config, inside_group_config, inside_generic_writing, inside_specific_file_writing
  integer :: current_building_definition, current_building_field, current_building_diagnostic, &
       current_trigger_index, current_building_group, current_building_file_writer
  type(io_configuration_type), save :: building_config !< IO configuration that is build built up from XML parsing
  character(len=STRING_LENGTH) :: data_handling_namespace

  type(hashmap_type) :: options_database
  type(hashset_type) :: data_field_names

  character(len=STRING_LENGTH), dimension(:), allocatable :: cond_request, diag_request, cond_long, diag_long
  integer :: ncond, ndiag
  
  logical :: l_thoff=.false.

  public EQ_OPERATOR_TYPE, LT_OPERATOR_TYPE, GT_OPERATOR_TYPE, LTE_OPERATOR_TYPE, GTE_OPERATOR_TYPE, ADD_OPERATOR_TYPE, &
       SUBTRACT_OPERATOR_TYPE, MULTIPLY_OPERATOR_TYPE, DIV_OPERATOR_TYPE, MOD_OPERATOR_TYPE, DATA_SIZE_STRIDE, &
       TIME_AVERAGED_TYPE, INSTANTANEOUS_TYPE, NONE_TYPE, GROUP_TYPE, FIELD_TYPE, IO_STATE_TYPE, handle_recv_data_from_io_server, &
       io_configuration_type, io_configuration_field_type, io_configuration_data_definition_type, &
       io_configuration_diagnostic_field_type, io_configuration_registered_monc_type, data_values_type, configuration_parse, &
       io_configuration_misc_item_type, io_configuration_inter_communication_description, &
       extend_registered_moncs_array, retrieve_data_definition, retrieve_monc_definition, extend_inter_io_comm_array, &
       build_definition_description_type_from_configuration, build_field_description_type_from_configuration, &
       get_number_field_dimensions, get_data_value_by_field_name, get_data_value_from_map_entry, get_monc_location, &
       get_diagnostic_field_configuration, get_prognostic_field_configuration, get_io_xml, &
       cond_request, diag_request, cond_long, diag_long, ncond, ndiag, l_thoff

contains

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
    logical :: comment_block

    comment_block = .false.

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
      if (adjusted_io_line(1:4) .eq. "<!--" .or. comment_block) then
        comment_block = .true.
        if (index(temp_line, "-->") .ne. 0) comment_block = .false.
        cycle
      end if
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

  !> This will parse an XML string into the IO configuration
  !! @param raw_configuration The raw (unparsed) XML string to process
  !! @param parsed_configuration Configuration determining the layout and handling of data
  subroutine configuration_parse(provided_options_database, raw_configuration, parsed_configuration)
    type(hashmap_type), intent(inout) :: provided_options_database
    character, dimension(:), intent(in) :: raw_configuration
    type(io_configuration_type), intent(out) :: parsed_configuration

    options_database=provided_options_database

    l_thoff = options_get_logical(provided_options_database, "l_thoff")

    inside_data_definition=.false.
    inside_handling_definition=.false.
    inside_server_config=.false.
    inside_diagnostic_config=.false.
    inside_group_config=.false.
    inside_generic_writing=.false.
    inside_specific_file_writing=.false.
    current_building_field=1
    current_building_definition=1
    current_building_diagnostic=1
    current_building_file_writer=1
    current_building_group=1
    current_trigger_index=1
    building_config%number_of_writers=0
    building_config%number_of_groups=0
    building_config%number_of_threads=-1
    allocate(building_config%data_definitions(DATA_SIZE_STRIDE))
    allocate(building_config%diagnostics(DATA_SIZE_STRIDE))
    allocate(building_config%groups(DATA_SIZE_STRIDE))
    allocate(building_config%file_writers(DATA_SIZE_STRIDE))
    allocate(building_config%inter_io_communications(DATA_SIZE_STRIDE))
    call xml_parse(raw_configuration, start_element_callback, end_element_callback)
    call add_in_dimensions(provided_options_database)
    building_config%options_database=options_database
    building_config%number_of_distinct_data_fields=c_size(data_field_names)
    building_config%number_of_moncs=0
    building_config%active_moncs=0
    allocate(building_config%registered_moncs(MONC_SIZE_STRIDE))
    parsed_configuration=building_config
    parsed_configuration%number_inter_io_communications=0
    parsed_configuration%general_info_set=.false.    
    call c_free(data_field_names)
    allocate(parsed_configuration%text_configuration(size(raw_configuration)), source=raw_configuration)
  end subroutine configuration_parse

  subroutine add_in_dimensions(provided_options_database)
    type(hashmap_type), intent(inout) :: provided_options_database

    integer :: dim_size, n_tracers=0

    call c_put_integer(building_config%dimension_sizing, "x", options_get_integer(provided_options_database, "x_size"))
    call c_put_integer(building_config%dimension_sizing, "y", options_get_integer(provided_options_database, "y_size"))

    dim_size=options_get_integer(provided_options_database, "z_size")
    call c_put_integer(building_config%dimension_sizing, "z", dim_size)
    call c_put_integer(building_config%dimension_sizing, "zn", dim_size)
    call c_put_integer(building_config%dimension_sizing, "qfields", &
         options_get_integer(provided_options_database, "number_q_fields"))
    
    call c_put_integer(building_config%dimension_sizing, "number_options", c_size(provided_options_database))
    call c_put_integer(building_config%dimension_sizing, "active_q_indicies", get_number_active_q_indices())
         
    if (options_get_logical(options_database, "tracers_enabled")) then
      if (options_get_logical(options_database, "trajectories_enabled")) then
        n_tracers = 5
      end if    
      if (options_get_logical(options_database, "radioactive_tracers_enabled")) then
        n_tracers = n_tracers + options_get_integer(provided_options_database, "n_radioactive_tracers")
      end if
    end if  ! tracers_enabled
    call c_put_integer(building_config%dimension_sizing, "tfields", n_tracers)

    !> Since the model appears to write each of the items in dimension_sizing as dimensions in every
    !  netcdf file, only let these exist when the corresponding code is enabled.
    if (options_get_logical(options_database, "conditional_diagnostics_column_enabled")) then
      ncond = options_get_array_size(options_database, "cond_request")*2 
      allocate(cond_request(ncond))
      allocate(cond_long(ncond))
      call c_put_integer(building_config%dimension_sizing, "nc", ncond)
      ndiag = options_get_array_size(options_database, "diag_request")
      allocate(diag_request(ndiag))
      allocate(diag_long(ndiag))
      call c_put_integer(building_config%dimension_sizing, "nd", ndiag)
    end if 
    if (options_get_logical(options_database, "pdf_analysis_enabled" )) then
      dim_size = options_get_integer(options_database, "n_w_bins")
      call c_put_integer(building_config%dimension_sizing, "n_w_bins", dim_size)
    end if 

  end subroutine add_in_dimensions  
  
  !> XML element start (opening) call back. This handles most of the configuration parsing
  !! @param element_name Name of the XML element
  !! @param number_of_attributes The number of attributes associated with this XML element
  !! @param attribute_names Each attribute name (same location as attribute value)
  !! @param attribute_values Each attribute value (same location as attribute name)
  subroutine start_element_callback(element_name, number_of_attributes, attribute_names, attribute_values)
    character(len=*), intent(in) :: element_name       
    character(len=*), dimension(:), intent(in) :: attribute_names, attribute_values
    integer, intent(in) :: number_of_attributes

    integer :: namespace_index

    if (inside_data_definition) then
      if (element_name == "field") then
        call process_xml_into_field_description(attribute_names, attribute_values)
      end if
    else if (inside_handling_definition) then
      if (inside_diagnostic_config) then
        call add_misc_member_to_diagnostic(element_name, attribute_names, attribute_values)
      else
        if (element_name == "diagnostic") then
          call define_diagnostic(attribute_names, attribute_values)
          inside_diagnostic_config=.true.        
        end if
      end if
    else if (inside_server_config) then
      if (element_name == "thread_pool") then
        call handle_thread_pool_configuration(attribute_names, attribute_values)
      end if
    else if (inside_group_config) then
      call add_diagnostic_field_to_group(element_name, attribute_names, attribute_values)
    else if (inside_generic_writing) then
      if (inside_specific_file_writing) then
        if (element_name == "include") then
          call add_include_to_file_writer(attribute_names, attribute_values)
        end if
      else if (element_name == "file") then
        inside_specific_file_writing=.true.
        call define_file_writer(attribute_names, attribute_values)
      end if
    else if (element_name == "data-writing") then
      inside_generic_writing=.true.
    else if (element_name == "data-definition") then
      inside_data_definition=.true.
      call handle_new_data_definition(attribute_names, attribute_values)
    else if (element_name == "data-handling") then
      namespace_index=get_field_index_from_name(attribute_names, "namespace")
      if (namespace_index == 0) then
        data_handling_namespace=""
      else
        data_handling_namespace=retrieve_string_value(attribute_values(namespace_index), STRING_DATA_TYPE)
      end if
      inside_handling_definition=.true.
    else if (element_name == "group") then
      call define_group(attribute_names, attribute_values)
      inside_group_config=.true.
    else if (element_name == "server-configuration") then
      inside_server_config=.true.
    end if
  end subroutine start_element_callback

  !> XML element end (closing) call back.
  !! @param element_name Name of the XML element
  subroutine end_element_callback(element_name)
    character(len=*), intent(in) :: element_name

    integer :: i

    if (element_name == "data-definition") then
      building_config%data_definitions(current_building_definition)%number_of_data_fields=current_building_field-1
      do i=1, current_building_field-1
        building_config%data_definitions(current_building_definition)%fields(i)%namespace=&
             building_config%data_definitions(current_building_definition)%namespace        
      end do
      current_building_field=1
      building_config%number_of_data_definitions=current_building_definition
      current_building_definition=current_building_definition+1
      inside_data_definition=.false.
    else if (element_name == "data-handling") then
       building_config%number_of_diagnostics=current_building_diagnostic-1
       inside_handling_definition=.false.
    else if (element_name == "server-configuration") then
      inside_server_config=.false.
    else if (element_name == "action") then
      inside_action_config=.false.
    else if (element_name == "diagnostic") then
      current_building_diagnostic=current_building_diagnostic+1
      inside_diagnostic_config=.false.
    else if (element_name == "data-writing") then
      inside_generic_writing=.false.
      building_config%number_of_writers=current_building_file_writer-1
    else if (element_name == "group") then
      current_building_group=current_building_group+1
      building_config%number_of_groups=building_config%number_of_groups+1
      inside_group_config=.false.
    else if (element_name == "file") then
      current_building_file_writer=current_building_file_writer+1
      inside_specific_file_writing=.false.
    end if    
  end subroutine end_element_callback

  subroutine handle_thread_pool_configuration(attribute_names, attribute_values)
    character(len=*), dimension(:), intent(in) :: attribute_names, attribute_values

    integer :: number_index

    number_index=get_field_index_from_name(attribute_names, "number")
    
    if (number_index /= 0) then
      building_config%number_of_threads=conv_to_integer(attribute_values(number_index))
    end if    
  end subroutine handle_thread_pool_configuration  

  !> Creates a new data definition configuration item based upon the attributes supplied
  !! @param number_of_attributes Number of XML attributes associated with this element
  !! @param attribute_names Each attribute name (same location as attribute value)
  !! @param attribute_values Each attribute value (same location as attribute name)
  subroutine handle_new_data_definition(attribute_names, attribute_values)
    character(len=*), dimension(:), intent(in) :: attribute_names, attribute_values

    integer :: name_index, frequency_index, namespace_index, send_on_termination_index
    character(len=STRING_LENGTH) :: namespace

    name_index=get_field_index_from_name(attribute_names, "name")
    frequency_index=get_field_index_from_name(attribute_names, "frequency")
    namespace_index=get_field_index_from_name(attribute_names, "namespace")
    send_on_termination_index=get_field_index_from_name(attribute_names, "send_on_terminate")
    if (name_index /= 0 .and. frequency_index /=0) then
      if (current_building_definition .gt. size(building_config%data_definitions)) call extend_data_definition_array()
      building_config%data_definitions(current_building_definition)%name=&
           retrieve_string_value(attribute_values(name_index), STRING_DATA_TYPE)
      building_config%data_definitions(current_building_definition)%frequency=&
           conv_to_integer(retrieve_string_value(attribute_values(frequency_index), INTEGER_DATA_TYPE))
      if (namespace_index /= 0) then
        namespace=retrieve_string_value(attribute_values(namespace_index), STRING_DATA_TYPE)
        building_config%data_definitions(current_building_definition)%namespace=namespace             
      else
        namespace=""
        building_config%data_definitions(current_building_definition)%namespace=""
      end if
      if (send_on_termination_index /= 0) then
        building_config%data_definitions(current_building_definition)%send_on_terminate=&
             retrieve_string_value(attribute_values(send_on_termination_index), STRING_DATA_TYPE) == "true"
      else
        building_config%data_definitions(current_building_definition)%send_on_terminate=.false.
      end if      
    else
      call log_log(LOG_ERROR, "A data definition requires a name and frequency")
    end if    
    allocate(building_config%data_definitions(current_building_definition)%fields(DATA_SIZE_STRIDE))

    building_config%data_definitions(current_building_definition)%fields(1)%name="timestep"
    building_config%data_definitions(current_building_definition)%fields(1)%namespace=namespace
    building_config%data_definitions(current_building_definition)%fields(1)%field_type=SCALAR_FIELD_TYPE
    building_config%data_definitions(current_building_definition)%fields(1)%dimensions=0
    building_config%data_definitions(current_building_definition)%fields(1)%data_type=INTEGER_DATA_TYPE
    building_config%data_definitions(current_building_definition)%fields(1)%optional=.false.
    building_config%data_definitions(current_building_definition)%fields(1)%collective=.false.
    building_config%data_definitions(current_building_definition)%fields(1)%units=""

    building_config%data_definitions(current_building_definition)%fields(2)%name="time"
    building_config%data_definitions(current_building_definition)%fields(2)%namespace=namespace
    building_config%data_definitions(current_building_definition)%fields(2)%field_type=SCALAR_FIELD_TYPE
    building_config%data_definitions(current_building_definition)%fields(2)%dimensions=0
    building_config%data_definitions(current_building_definition)%fields(2)%data_type=DOUBLE_DATA_TYPE
    building_config%data_definitions(current_building_definition)%fields(2)%optional=.false.
    building_config%data_definitions(current_building_definition)%fields(2)%collective=.false.
    building_config%data_definitions(current_building_definition)%fields(2)%units=""

    building_config%data_definitions(current_building_definition)%fields(3)%name="terminated"
    building_config%data_definitions(current_building_definition)%fields(3)%namespace=namespace
    building_config%data_definitions(current_building_definition)%fields(3)%field_type=SCALAR_FIELD_TYPE
    building_config%data_definitions(current_building_definition)%fields(3)%dimensions=0
    building_config%data_definitions(current_building_definition)%fields(3)%data_type=BOOLEAN_DATA_TYPE
    building_config%data_definitions(current_building_definition)%fields(3)%optional=.false.
    building_config%data_definitions(current_building_definition)%fields(3)%collective=.false.
    building_config%data_definitions(current_building_definition)%fields(3)%units=""

    current_building_field=4
  end subroutine handle_new_data_definition

  subroutine add_misc_member_to_diagnostic(element_name, attribute_names, attribute_values)
    character(len=*), intent(in) :: element_name       
    character(len=*), dimension(:), intent(in) :: attribute_names, attribute_values

    type(io_configuration_misc_item_type), pointer :: misc_member
    class(*), pointer :: generic
    integer :: i

    allocate(misc_member)

    misc_member%type=retrieve_string_value(element_name, STRING_DATA_TYPE)
    misc_member%namespace=data_handling_namespace
    do i=1, size(attribute_names)
      call c_put_string(misc_member%embellishments, retrieve_string_value(attribute_names(i), STRING_DATA_TYPE), &
           retrieve_string_value(attribute_values(i), STRING_DATA_TYPE))
    end do
    generic=>misc_member
    call c_add_generic(building_config%diagnostics(current_building_diagnostic)%members, generic, .false.)
  end subroutine add_misc_member_to_diagnostic

  subroutine add_diagnostic_field_to_group(element_name, attribute_names, attribute_values)
    character(len=*), intent(in) :: element_name       
    character(len=*), dimension(:), intent(in) :: attribute_names, attribute_values

    integer :: field_index

    if (element_name == "member") then
      field_index=get_field_index_from_name(attribute_names, "name")
      if (field_index .gt. 0) then
        call c_add_string(building_config%groups(current_building_group)%members, &
             retrieve_string_value(attribute_values(field_index), STRING_DATA_TYPE))
      else
        call log_log(LOG_ERROR, "A diagnostics group member requires a name")
      end if
    else      
        call log_log(LOG_ERROR, "Unrecognised diagnostics group participant, name is '"//trim(element_name)//"'")
    end if
  end subroutine add_diagnostic_field_to_group  

  subroutine define_group(attribute_names, attribute_values)      
    character(len=*), dimension(:), intent(in) :: attribute_names, attribute_values
   
    integer :: field_index

    if (current_building_group .gt. size(building_config%groups)) call extend_groups_array()

    field_index=get_field_index_from_name(attribute_names, "name")
    if (field_index .gt. 0) then
      building_config%groups(current_building_group)%name=retrieve_string_value(attribute_values(field_index), STRING_DATA_TYPE)
      field_index=get_field_index_from_name(attribute_names, "namespace")
      if (field_index .gt. 0) then
        building_config%groups(current_building_group)%namespace=&
             retrieve_string_value(attribute_values(field_index), STRING_DATA_TYPE)
      else
        building_config%groups(current_building_group)%namespace=""
      end if
    else
      call log_log(LOG_ERROR, "A diagnostics group requires a name")
    end if
  end subroutine define_group

  subroutine add_include_to_file_writer(attribute_names, attribute_values)
    character(len=*), dimension(:), intent(in) :: attribute_names, attribute_values

    integer :: field_index, number_of_contents

    number_of_contents=building_config%file_writers(current_building_file_writer)%number_of_contents+1

    if (number_of_contents .gt. &
         size(building_config%file_writers(current_building_file_writer)%contents)) call extend_file_writer_contents_array()

    building_config%file_writers(current_building_file_writer)%contents(number_of_contents)%facet_type=0
    field_index=get_field_index_from_name(attribute_names, "group")
    if (field_index .gt. 0) then
      building_config%file_writers(current_building_file_writer)%contents(number_of_contents)%facet_name=&
           retrieve_string_value(attribute_values(field_index), STRING_DATA_TYPE)
      building_config%file_writers(current_building_file_writer)%contents(number_of_contents)%facet_type=GROUP_TYPE
      call add_include_group_or_field_to_file_writer(attribute_names, attribute_values, number_of_contents)
    else
      field_index=get_field_index_from_name(attribute_names, "field")
      if (field_index .gt. 0) then
        building_config%file_writers(current_building_file_writer)%contents(number_of_contents)%facet_name=&
             retrieve_string_value(attribute_values(field_index), STRING_DATA_TYPE)
        building_config%file_writers(current_building_file_writer)%contents(number_of_contents)%facet_type=FIELD_TYPE
        call add_include_group_or_field_to_file_writer(attribute_names, attribute_values, number_of_contents)
      else
        field_index=get_field_index_from_name(attribute_names, "state")
        if (field_index .gt. 0) then
          if (trim(retrieve_string_value(attribute_values(field_index), STRING_DATA_TYPE)) .eq. "io") then
            building_config%file_writers(current_building_file_writer)%contents(number_of_contents)%facet_type=IO_STATE_TYPE
          end if
        end if        
      end if      
    end if

    if (building_config%file_writers(current_building_file_writer)%contents(number_of_contents)%facet_type == 0) then
      call log_log(LOG_ERROR, "Inclusion to file writer requires a field or group to include")
    end if

    building_config%file_writers(current_building_file_writer)%number_of_contents=&
         building_config%file_writers(current_building_file_writer)%number_of_contents+1
  end subroutine add_include_to_file_writer  
  
  subroutine add_include_group_or_field_to_file_writer(attribute_names, attribute_values, number_of_contents)
    character(len=*), dimension(:), intent(in) :: attribute_names, attribute_values
    integer, intent(in) :: number_of_contents

    character(len=STRING_LENGTH) :: time_manip
    integer :: field_index

    field_index=get_field_index_from_name(attribute_names, "time_manipulation")
    if (field_index .gt. 0) then
      time_manip=retrieve_string_value(attribute_values(field_index), STRING_DATA_TYPE)
      if (time_manip == "instantaneous") then
        building_config%file_writers(current_building_file_writer)%contents(number_of_contents)%time_manipulation_type=&
             INSTANTANEOUS_TYPE
      else if (time_manip == "averaged") then
        building_config%file_writers(current_building_file_writer)%contents(number_of_contents)%time_manipulation_type=&
             TIME_AVERAGED_TYPE
      else if (time_manip == "none") then
        building_config%file_writers(current_building_file_writer)%contents(number_of_contents)%time_manipulation_type=&
             NONE_TYPE
      else
        call log_log(LOG_ERROR, "Time manipulation '"//trim(time_manip)//"' option not recognised")
      end if
    else
      call log_log(LOG_ERROR, "Inclusion to file writer requires time manipulation")
    end if

    field_index=get_field_index_from_name(attribute_names, "output_frequency")
    if (field_index .gt. 0) then
      building_config%file_writers(current_building_file_writer)%&
           contents(number_of_contents)%output_time_frequency=&
           conv_to_real(retrieve_string_value(attribute_values(field_index), DOUBLE_DATA_TYPE))
    else if (building_config%file_writers(current_building_file_writer)%contents(number_of_contents)%time_manipulation_type &
         == NONE_TYPE) then
      building_config%file_writers(current_building_file_writer)%contents(number_of_contents)%output_time_frequency=-1.0
    else
      call log_log(LOG_ERROR, "Inclusion to file writer requires an output frequency")
    end if
  end subroutine add_include_group_or_field_to_file_writer

  subroutine define_file_writer(attribute_names, attribute_values)     
    character(len=*), dimension(:), intent(in) :: attribute_names, attribute_values

    integer :: field_index


    if (current_building_file_writer .gt. size(building_config%file_writers)) call extend_file_writer_array()

    field_index=get_field_index_from_name(attribute_names, "name")
    if (field_index .gt. 0) then
      building_config%file_writers(current_building_file_writer)%file_name=&
           retrieve_string_value(attribute_values(field_index), STRING_DATA_TYPE)
    else
      call log_log(LOG_ERROR, "File writer requires a file name")
    end if

    field_index=get_field_index_from_name(attribute_names, "write_time_frequency")
    if (field_index .gt. 0) then    
      building_config%file_writers(current_building_file_writer)%write_time_frequency=&
           conv_to_real(retrieve_string_value(attribute_values(field_index), DOUBLE_DATA_TYPE))
      building_config%file_writers(current_building_file_writer)%write_on_model_time=.true.      
    else
      field_index=get_field_index_from_name(attribute_names, "write_timestep_frequency")
      if (field_index .gt. 0) then
        building_config%file_writers(current_building_file_writer)%write_timestep_frequency=&
             conv_to_real(retrieve_string_value(attribute_values(field_index), INTEGER_DATA_TYPE))
        building_config%file_writers(current_building_file_writer)%write_on_model_time=.false. 
      else
        call log_log(LOG_ERROR, "File writer requires either a write time frequency or write timestep frequency")
      end if
    end if

    field_index=get_field_index_from_name(attribute_names, "title")
    if (field_index .gt. 0) then 
      building_config%file_writers(current_building_file_writer)%title=&
           retrieve_string_value(attribute_values(field_index), STRING_DATA_TYPE)
    else
      building_config%file_writers(current_building_file_writer)%title=DEFAULT_FILE_TITLE
    end if    

    field_index=get_field_index_from_name(attribute_names, "write_on_terminate")
    if (field_index .gt. 0) then 
      building_config%file_writers(current_building_file_writer)%write_on_terminate=&
           retrieve_string_value(attribute_values(field_index), STRING_DATA_TYPE) == "true"
      if (building_config%file_writers(current_building_file_writer)%write_on_model_time) &
        call log_log(LOG_ERROR, "Inconsitent settings.  write_on_terminate cannot be used with write_on_model_time")
    else
      building_config%file_writers(current_building_file_writer)%write_on_terminate=.false.
    end if

    field_index=get_field_index_from_name(attribute_names, "store_state")
    if (field_index .gt. 0) then 
      building_config%file_writers(current_building_file_writer)%include_in_io_state_write=&
           retrieve_string_value(attribute_values(field_index), STRING_DATA_TYPE) == "true"
    else
      building_config%file_writers(current_building_file_writer)%include_in_io_state_write=.true.
    end if

    field_index=get_field_index_from_name(attribute_names, "write_precision")
    if (field_index .gt. 0) then
      building_config%file_writers(current_building_file_writer)%write_precision=&
           merge(NF90_REAL, NF90_DOUBLE, &
                 conv_to_uppercase(retrieve_string_value(attribute_values(field_index), STRING_DATA_TYPE)) == "FLOAT")
      call log_master_newline()
      call log_master_log(LOG_INFO, "Data will be written in FLOAT precision for file: "//&
                   trim(building_config%file_writers(current_building_file_writer)%file_name))
    else
      building_config%file_writers(current_building_file_writer)%write_precision = NF90_DOUBLE
    end if

    
    building_config%file_writers(current_building_file_writer)%number_of_contents=0
    allocate(building_config%file_writers(current_building_file_writer)%contents(DATA_SIZE_STRIDE))    
  end subroutine define_file_writer  

  !> Defines a new data handling rule
  !! @param element_name The name of the XML element
  !! @param number_of_attributes Number of XML attributes associated with this element
  !! @param attribute_names Each attribute name (same location as attribute value)
  !! @param attribute_values Each attribute value (same location as attribute name)
  subroutine define_diagnostic(attribute_names, attribute_values)   
    character(len=*), dimension(:), intent(in) :: attribute_names, attribute_values
   
    integer :: field_index, type_field_index, data_field_index
    character(len=STRING_LENGTH) :: field_type_str, field_data_type_str, size_definitions

    if (current_building_diagnostic .gt. size(building_config%diagnostics)) call extend_diagnostics_array()

    field_index=get_field_index_from_name(attribute_names, "field")
    type_field_index=get_field_index_from_name(attribute_names, "type")
    data_field_index=get_field_index_from_name(attribute_names, "data_type")

    if (field_index == 0 .or. type_field_index == 0 .or. data_field_index == 0) then
      call log_log(LOG_ERROR, "Each diagnostic definition requires a field name, field type and data type")
    else
      building_config%diagnostics(current_building_diagnostic)%name=&
           retrieve_string_value(attribute_values(field_index), STRING_DATA_TYPE)
      building_config%diagnostics(current_building_diagnostic)%namespace=data_handling_namespace
      field_type_str=retrieve_string_value(attribute_values(type_field_index), STRING_DATA_TYPE)
      building_config%diagnostics(current_building_diagnostic)%field_type=get_field_type_from_attribute(field_type_str)
      if (building_config%diagnostics(current_building_diagnostic)%field_type==0) then
        call log_log(LOG_ERROR, "The field type of '"//trim(field_type_str)//"' is not recognised")
      end if
      field_data_type_str=retrieve_string_value(attribute_values(data_field_index), STRING_DATA_TYPE)
      building_config%diagnostics(current_building_diagnostic)%data_type=get_field_datatype_from_attribute(field_data_type_str)
      if (building_config%diagnostics(current_building_diagnostic)%data_type==0) then
        call log_log(LOG_ERROR, "The field data type of '"//trim(field_data_type_str)//"' is not recognised")
      end if
      field_index=get_field_index_from_name(attribute_names, "units")
      if (field_index .ne. 0) then
        building_config%diagnostics(current_building_diagnostic)%units=&
             retrieve_string_value(attribute_values(field_index), STRING_DATA_TYPE)
      else
        building_config%diagnostics(current_building_diagnostic)%units=""
      end if
      if (building_config%diagnostics(current_building_diagnostic)%field_type == ARRAY_FIELD_TYPE .or. &
           building_config%diagnostics(current_building_diagnostic)%field_type == MAP_FIELD_TYPE) then
        field_index=get_field_index_from_name(attribute_names, "size")
        if (field_index .ne. 0) then
          size_definitions=retrieve_string_value(attribute_values(field_index), STRING_DATA_TYPE)      
          building_config%diagnostics(current_building_diagnostic)%dimensions=process_sizing_definition(size_definitions, &
               building_config%diagnostics(current_building_diagnostic)%dim_size_defns)
        else
          call log_log(LOG_ERROR, "A diagnostic of field type array or map requires sizing a definition")
        end if
        field_index=get_field_index_from_name(attribute_names, "collective")
        if (field_index .ne. 0) then
          building_config%diagnostics(current_building_diagnostic)%collective=&
               retrieve_string_value(attribute_values(field_index), STRING_DATA_TYPE) == "true" 
        else
          building_config%diagnostics(current_building_diagnostic)%collective=.false.
        end if
      else
        building_config%diagnostics(current_building_diagnostic)%dimensions=0
        building_config%diagnostics(current_building_diagnostic)%collective=.false.
      end if
    end if
  end subroutine define_diagnostic

  integer function process_sizing_definition(size_definitions, individual_str_defn)
    character(len=*), intent(in) :: size_definitions
    character(len=STRING_LENGTH), intent(out) :: individual_str_defn(4)

    integer :: comma_index, sizing_index, cp

    cp=1
    sizing_index=1
    comma_index=index(size_definitions(cp:), ",")
    do while (comma_index .ne. 0)
      comma_index=comma_index+cp-1
      individual_str_defn(sizing_index)= trim(size_definitions(cp:comma_index-1))
      cp=comma_index+1
      sizing_index=sizing_index+1
      comma_index=index(size_definitions(cp:), ",")
      if (sizing_index .gt. 4) then
        call log_log(LOG_ERROR, "Can only have a maximum of four diagnostic field sizing dimensions")
      end if
    end do
    if (cp .le. len(size_definitions)) then
      individual_str_defn(sizing_index)=trim(size_definitions(cp:))
      sizing_index=sizing_index+1
    end if    
    process_sizing_definition=sizing_index-1
  end function process_sizing_definition

  !> Process XML into a field description by identifying the attributes of the field and storing these in a refined
  !! configuration format
  !! @param element_name The name of the XML element
  !! @param number_of_attributes Number of XML attributes associated with this element
  !! @param attribute_names Each attribute name (same location as attribute value)
  !! @param attribute_values Each attribute value (same location as attribute name)
  subroutine process_xml_into_field_description(attribute_names, attribute_values)    
    character(len=*), dimension(:), intent(in) :: attribute_names, attribute_values

    character(len=STRING_LENGTH) :: field_type_str, field_data_type_str, sizing_defn_str
    integer :: name_field_index, type_field_index, data_field_index, field_index, optional_field_index, idx
    

    name_field_index=get_field_index_from_name(attribute_names, "name")
    type_field_index=get_field_index_from_name(attribute_names, "type")
    data_field_index=get_field_index_from_name(attribute_names, "data_type")

    if (name_field_index == 0 .or. type_field_index == 0 .or. data_field_index == 0) then
      call log_log(LOG_ERROR, "Each data field definition requires a name, field type and data type")
    else
      if (current_building_field .eq. size(building_config%data_definitions(current_building_definition)%fields)) &
           call extend_field_array()

      building_config%data_definitions(current_building_definition)%fields(current_building_field)%name=&
           retrieve_string_value(attribute_values(name_field_index), STRING_DATA_TYPE)

      call c_add_string(data_field_names, &
           building_config%data_definitions(current_building_definition)%fields(current_building_field)%name)

      field_type_str=retrieve_string_value(attribute_values(type_field_index), STRING_DATA_TYPE)
      building_config%data_definitions(current_building_definition)%fields(current_building_field)%field_type=&
           get_field_type_from_attribute(field_type_str)

      if (building_config%data_definitions(current_building_definition)%fields(current_building_field)%field_type==0) then
        call log_log(LOG_ERROR, "The field type of '"//trim(field_type_str)//"' is not recognised")
      end if

      building_config%data_definitions(current_building_definition)%fields(current_building_field)%dimensions=0

      if (building_config%data_definitions(current_building_definition)%fields(current_building_field)%field_type == &
           ARRAY_FIELD_TYPE .or. &
           building_config%data_definitions(current_building_definition)%fields(current_building_field)%field_type == &
           MAP_FIELD_TYPE) then
        idx=get_field_index_from_name(attribute_names, "size")
        if (idx .ne. 0) then
          sizing_defn_str=retrieve_string_value(attribute_values(idx), STRING_DATA_TYPE)
          building_config%data_definitions(current_building_definition)%fields(current_building_field)%dimensions=&
               process_sizing_definition(sizing_defn_str, &
               building_config%data_definitions(current_building_definition)%fields(current_building_field)%dim_size_defns)
        end if

        field_index=get_field_index_from_name(attribute_names, "collective")
        if (field_index .ne. 0) then
          building_config%data_definitions(current_building_definition)%fields(current_building_field)%collective=&
               retrieve_string_value(attribute_values(field_index), STRING_DATA_TYPE) == "true" 
        else
          building_config%data_definitions(current_building_definition)%fields(current_building_field)%collective=.false.
        end if
      else
        building_config%data_definitions(current_building_definition)%fields(current_building_field)%collective=.false.
      end if

      field_index=get_field_index_from_name(attribute_names, "units")
      if (field_index .ne. 0) then
        building_config%data_definitions(current_building_definition)%fields(current_building_field)%units=&
             retrieve_string_value(attribute_values(field_index), STRING_DATA_TYPE)
      else
        building_config%data_definitions(current_building_definition)%fields(current_building_field)%units=""
      end if
            
      field_data_type_str=retrieve_string_value(attribute_values(data_field_index), STRING_DATA_TYPE)
      building_config%data_definitions(current_building_definition)%fields(current_building_field)%data_type=&
           get_field_datatype_from_attribute(field_data_type_str)

      if (building_config%data_definitions(current_building_definition)%fields(current_building_field)%data_type==0) then
        call log_log(LOG_ERROR, "The field data type of '"//trim(field_data_type_str)//"' is not recognised")
      end if      

      if (building_config%data_definitions(current_building_definition)%fields(current_building_field)%field_type == &
           MAP_FIELD_TYPE .and. building_config%data_definitions(current_building_definition)%fields(&
           current_building_field)%data_type/=STRING_DATA_TYPE) then
        call log_log(LOG_ERROR, "A map field type must have a data type of ""string""")
      end if

      optional_field_index=get_field_index_from_name(attribute_names, "optional")
      if (optional_field_index .ne. 0) then
        if (retrieve_string_value(attribute_values(optional_field_index), STRING_DATA_TYPE) == "true") then
          building_config%data_definitions(current_building_definition)%fields(current_building_field)%optional=.true.
        else
          building_config%data_definitions(current_building_definition)%fields(current_building_field)%optional=.false.
        end if
      else
        building_config%data_definitions(current_building_definition)%fields(current_building_field)%optional=.false.
      end if      
      current_building_field=current_building_field+1
    end if
  end subroutine process_xml_into_field_description

  integer function get_field_type_from_attribute(field_type_str)
    character(len=*), intent(in) :: field_type_str

    get_field_type_from_attribute=0
    if (field_type_str == "scalar") get_field_type_from_attribute=SCALAR_FIELD_TYPE
    if (field_type_str == "array") get_field_type_from_attribute=ARRAY_FIELD_TYPE
    if (field_type_str == "map") get_field_type_from_attribute=MAP_FIELD_TYPE
  end function get_field_type_from_attribute  

  integer function get_field_datatype_from_attribute(field_data_type_str)
    character(len=*), intent(in) :: field_data_type_str

    get_field_datatype_from_attribute=0
    if (field_data_type_str == "integer") get_field_datatype_from_attribute=INTEGER_DATA_TYPE
    if (field_data_type_str == "boolean") get_field_datatype_from_attribute=BOOLEAN_DATA_TYPE
    if (field_data_type_str == "string") get_field_datatype_from_attribute=STRING_DATA_TYPE
    if (field_data_type_str == "float") get_field_datatype_from_attribute=FLOAT_DATA_TYPE
    if (field_data_type_str == "double") get_field_datatype_from_attribute=DOUBLE_DATA_TYPE
  end function get_field_datatype_from_attribute  

  !> Replaces specific characters in a string and returns a new string with this replaced by nothing (i.e. removed)
  !! @param original_string The original string the process
  !! @param new_string New string which contains the string after replacement has been performed
  !! @param to_replace The substring to replace by nothing (i.e. remove) from original_string
  subroutine replace_characters_in_string(original_string, new_string, to_replace)
    character(len=*), intent(in) :: original_string, to_replace
    character(len=*), intent(out) :: new_string

    integer :: current_index, string_len, occurance

    current_index=1
    string_len=len(original_string)
    new_string=""
    do while (current_index .lt. string_len)
      occurance=index(original_string(current_index:), to_replace)
      if (occurance .eq. 0) then
        occurance=len(original_string)
      else
        occurance=occurance+current_index
      end if      
      new_string=trim(new_string)//trim(original_string(current_index:occurance-len(to_replace)-1))
      current_index=current_index+occurance+len(to_replace)-2
    end do
  end subroutine replace_characters_in_string

  character(len=STRING_LENGTH) function retrieve_string_value(original_string, field_value_type)
    character(len=*), intent(in) :: original_string
    integer, intent(in) :: field_value_type

    integer :: last_char
    character(len=STRING_LENGTH) :: lookup_key

    call replace_characters_in_string(original_string, retrieve_string_value, """")

    retrieve_string_value=adjustl(retrieve_string_value)
    if (retrieve_string_value(1:1)=="{") then
      last_char=len_trim(retrieve_string_value)
      if (retrieve_string_value(last_char:last_char)=="}") then
        lookup_key=retrieve_string_value(2:last_char-1)
        if (options_has_key(options_database, lookup_key)) then
          if (field_value_type==INTEGER_DATA_TYPE) then
            retrieve_string_value=conv_to_string(options_get_integer(options_database, lookup_key))
          else if (field_value_type==BOOLEAN_DATA_TYPE) then
            retrieve_string_value=conv_to_string(options_get_logical(options_database, lookup_key))
          else if (field_value_type==DOUBLE_DATA_TYPE) then
            retrieve_string_value=conv_to_string(options_get_real(options_database, lookup_key))
          else if (field_value_type==STRING_DATA_TYPE) then
            retrieve_string_value=options_get_string(options_database, lookup_key)
          end if          
        else
          call log_log(LOG_ERROR, "Can not find IO configuration key '"//trim(lookup_key)//"' in the options database")
        end if
      end if      
    end if    
  end function retrieve_string_value  

  !> Given the name of an attribute will return the index of this in the names collection or 0 if it is not found
  !! @param attribute_names Collection of attribute names
  !! @param search_name The name to search for
  integer function get_field_index_from_name(attribute_names, search_name)
    character(len=*), intent(in) :: search_name       
    character(len=*), dimension(:) :: attribute_names

    integer :: i, size_of_names

    size_of_names=size(attribute_names)
    do i=1,size_of_names
      if (attribute_names(i) == search_name) then
        get_field_index_from_name=i
        return
      end if      
    end do  
    get_field_index_from_name=0
  end function get_field_index_from_name  

  !> Extends the array of inter io communications from its current suze to current size+data_stride+current size deficit
  !! @param io_configuration The IO server configuration state
  !! @param inter_io_size The target number of elements in the array
  subroutine extend_inter_io_comm_array(io_configuration)
    type(io_configuration_type), intent(inout) :: io_configuration

    type(io_configuration_inter_communication_description), dimension(:), allocatable :: temp_descriptions

    allocate(temp_descriptions(lbound(io_configuration%inter_io_communications,1):&
         ubound(io_configuration%inter_io_communications,1)+DATA_SIZE_STRIDE))
    temp_descriptions(lbound(io_configuration%inter_io_communications,1):&
         ubound(io_configuration%inter_io_communications,1)) = io_configuration%inter_io_communications
    call move_alloc(from=temp_descriptions,to=io_configuration%inter_io_communications)
  end subroutine extend_inter_io_comm_array

  subroutine extend_file_writer_contents_array()
    type(io_configuration_file_writer_facet_type), dimension(:), allocatable :: temp_filewriter_contents

    allocate(temp_filewriter_contents(lbound(building_config%file_writers(current_building_file_writer)%contents,1): &
         ubound(building_config%file_writers(current_building_file_writer)%contents,1)+DATA_SIZE_STRIDE)) 
    temp_filewriter_contents(lbound(building_config%file_writers(current_building_file_writer)%contents,1):&
         ubound(building_config%file_writers(current_building_file_writer)%contents,1)) = &
         building_config%file_writers(current_building_file_writer)%contents
    call move_alloc(from=temp_filewriter_contents,to=building_config%file_writers(current_building_file_writer)%contents)
  end subroutine extend_file_writer_contents_array

  subroutine extend_file_writer_array()
    type(io_configuration_file_writer_type), dimension(:), allocatable :: temp_filewriter

    allocate(temp_filewriter(lbound(building_config%file_writers,1): ubound(building_config%file_writers,1)+DATA_SIZE_STRIDE)) 
    temp_filewriter(lbound(building_config%file_writers,1):ubound(building_config%file_writers,1)) = building_config%file_writers
    call move_alloc(from=temp_filewriter,to=building_config%file_writers)
  end subroutine extend_file_writer_array  

  !> Extends the rules array of a specific rule from the current size to the current size + data size stride
  subroutine extend_diagnostics_array()
    type(io_configuration_diagnostic_field_type), dimension(:), allocatable :: temp_diagnostics

    allocate(temp_diagnostics(lbound(building_config%diagnostics,1): ubound(building_config%diagnostics,1)+DATA_SIZE_STRIDE)) 
    temp_diagnostics(lbound(building_config%diagnostics,1):ubound(building_config%diagnostics,1)) = building_config%diagnostics
    call move_alloc(from=temp_diagnostics,to=building_config%diagnostics)
  end subroutine extend_diagnostics_array  

  subroutine extend_groups_array()
    type(io_configuration_group_type), dimension(:), allocatable :: temp_groups

    allocate(temp_groups(lbound(building_config%groups,1): ubound(building_config%groups,1)+DATA_SIZE_STRIDE)) 
    temp_groups(lbound(building_config%groups,1):ubound(building_config%groups,1)) = building_config%groups
    call move_alloc(from=temp_groups,to=building_config%groups)
  end subroutine extend_groups_array  

  !> Extends the fields array of the current data definition from the current size to the current size + data size stride
  subroutine extend_field_array()
    type(io_configuration_field_type), dimension(:), allocatable :: temp_fields

    allocate(temp_fields(lbound(building_config%data_definitions(current_building_definition)%fields,1):&
         ubound(building_config%data_definitions(current_building_definition)%fields,1)+DATA_SIZE_STRIDE)) 
    temp_fields(lbound(building_config%data_definitions(current_building_definition)%fields,1):&
         ubound(building_config%data_definitions(current_building_definition)%fields,1)) = &
         building_config%data_definitions(current_building_definition)%fields
    call move_alloc(from=temp_fields,to=building_config%data_definitions(current_building_definition)%fields)
  end subroutine extend_field_array

  !> Extends the data definitions array from the current size to the current size + data size stride
  subroutine extend_data_definition_array()
    type(io_configuration_data_definition_type), dimension(:), allocatable :: temp_data_definitions

    allocate(temp_data_definitions(lbound(building_config%data_definitions, 1):&
         ubound(building_config%data_definitions,1)+DATA_SIZE_STRIDE)) 
    temp_data_definitions(lbound(building_config%data_definitions, 1):&
         ubound(building_config%data_definitions, 1)) = building_config%data_definitions
    call move_alloc(from=temp_data_definitions,to=building_config%data_definitions)
  end subroutine extend_data_definition_array 

  !> Extends the data definitions array from the current size to the current size + data size stride
  !! @param io_configuration IO server configuration state
  subroutine extend_registered_moncs_array(io_configuration)
    type(io_configuration_type), intent(inout) :: io_configuration

    type(io_configuration_registered_monc_type), dimension(:), allocatable :: temp_registered_moncs

    allocate(temp_registered_moncs(lbound(io_configuration%registered_moncs, 1):&
         io_configuration%number_of_moncs+MONC_SIZE_STRIDE)) 
    temp_registered_moncs(lbound(io_configuration%registered_moncs, 1):&
         ubound(io_configuration%registered_moncs, 1)) = io_configuration%registered_moncs
    call move_alloc(from=temp_registered_moncs,to=io_configuration%registered_moncs)
  end subroutine extend_registered_moncs_array

  !> Retrieves a specific data definition from the configuration which matches a key
  !! @param io_configuration IO server configuration state
  !! @param key Data definition key that we are searching for
  !! @returns Index of the found data definition or zero if none is found
  integer function retrieve_data_definition(io_configuration, key)
    type(io_configuration_type), intent(inout) :: io_configuration
    character(len=*), intent(in) :: key

    integer :: i

    do i=1,io_configuration%number_of_data_definitions
      if (io_configuration%data_definitions(i)%name .eq. key) then
        retrieve_data_definition=i
        return
      end if      
    end do    
    retrieve_data_definition=0
  end function retrieve_data_definition

  !> Retrieves a specific MONC definition from the configuration which matches a source PID
  !! @param io_configuration IO server configuration state
  !! @param source MONC process ID to search for
  !! @param monc_defn The MONC definition, if it is found, is written into here
  !! @returns Whether a matching MONC definition is found or not
  logical function retrieve_monc_definition(io_configuration, source, monc_defn)
    type(io_configuration_type), intent(inout) :: io_configuration
    integer, intent(in) :: source
    type(io_configuration_registered_monc_type), intent(out) :: monc_defn

    class(*), pointer :: generic
    integer :: location

    generic=>c_get_generic(io_configuration%monc_to_index, conv_to_string(source))
    if (associated(generic)) then
      location=conv_to_integer(generic, .false.)
      monc_defn=io_configuration%registered_moncs(location)
      retrieve_monc_definition=.true.
    else
      retrieve_monc_definition=.false.
    end if
  end function retrieve_monc_definition

  !> Builds up the data definition description type from the structured definitions in the IO configuration
  !! @param io_configuration IO server configuration
  !! @returns Array of data definition descriptions ready to be sent to a MONC process when it registers
  function build_definition_description_type_from_configuration(io_configuration)
    type(io_configuration_type), intent(inout) :: io_configuration
    type(definition_description_type), dimension(:), allocatable :: build_definition_description_type_from_configuration

    integer :: i

    allocate(build_definition_description_type_from_configuration(io_configuration%number_of_data_definitions))

    do i=1, io_configuration%number_of_data_definitions
      build_definition_description_type_from_configuration(i)%definition_name=io_configuration%data_definitions(i)%name
      build_definition_description_type_from_configuration(i)%send_on_terminate=&
           io_configuration%data_definitions(i)%send_on_terminate
      build_definition_description_type_from_configuration(i)%number_fields=&
           io_configuration%data_definitions(i)%number_of_data_fields
      build_definition_description_type_from_configuration(i)%frequency=io_configuration%data_definitions(i)%frequency
    end do
  end function build_definition_description_type_from_configuration

  !> Builds up the field definition description type from the structured definitions in the IO configuration
  !! @param io_configuration IO server configuration
  !! @returns Array of field definition descriptions ready to be sent to a MONC process when it registers
  function build_field_description_type_from_configuration(io_configuration)
    type(io_configuration_type), intent(inout) :: io_configuration
    type(field_description_type), dimension(:), allocatable :: build_field_description_type_from_configuration

    integer :: i, j, field_index

    allocate(build_field_description_type_from_configuration(get_total_number_fields(io_configuration)))
    field_index=1
    do i=1, io_configuration%number_of_data_definitions
      do j=1, io_configuration%data_definitions(i)%number_of_data_fields
        build_field_description_type_from_configuration(field_index)%definition_name=io_configuration%data_definitions(i)%name
        build_field_description_type_from_configuration(field_index)%field_name=&
             io_configuration%data_definitions(i)%fields(j)%name
        build_field_description_type_from_configuration(field_index)%field_type=&
             io_configuration%data_definitions(i)%fields(j)%field_type
        build_field_description_type_from_configuration(field_index)%data_type=&
             io_configuration%data_definitions(i)%fields(j)%data_type
        build_field_description_type_from_configuration(field_index)%optional=&
             io_configuration%data_definitions(i)%fields(j)%optional
        field_index=field_index+1
      end do      
    end do    
  end function build_field_description_type_from_configuration  

  !> Retrieves the total number of fields held in all data definitions
  !! @param io_configuration The IO server configuration
  !! @returns Total number of fields
  integer function get_total_number_fields(io_configuration)
    type(io_configuration_type), intent(inout) :: io_configuration

    integer :: i

    get_total_number_fields=0
    do i=1, io_configuration%number_of_data_definitions
      get_total_number_fields=get_total_number_fields+io_configuration%data_definitions(i)%number_of_data_fields
    end do    
  end function get_total_number_fields
  
  !> Retrieves the number of field dimensions that a specific field has from a MONC process within a data definition. This
  !! is determined by the MONC process when it registers with the IO server and can be different from one to another
  !! @param io_configuration The IO server configuration
  !! @param field_name The field name to look up
  !! @param source The pid of the source MONC process
  !! @param data_id The data definition id
  !! @returns The number of dimensions for this field, within the data definition from the source MONC process
  integer function get_number_field_dimensions(io_configuration, field_name, source, data_id)
    type(io_configuration_type), intent(inout) :: io_configuration
    character(len=*), intent(in) :: field_name
    integer, intent(in) :: source, data_id

    integer :: monc_location

    monc_location=c_get_integer(io_configuration%monc_to_index, conv_to_string(source))
    get_number_field_dimensions=c_get_integer(io_configuration%registered_moncs(monc_location)%dimensions(data_id), field_name)
  end function get_number_field_dimensions

  !> Retrieves the data value (wrapper) by field name or null if no entry was found in the provided collection
  !! @param collection A hashmap to search for this data value in
  !! @param field_name The field name to search for
  !! @returns The corresponding data value or null if none is found
  function get_data_value_from_hashmap_by_field_name(collection, field_name)
    type(hashmap_type), intent(inout) :: collection
    character(len=*), intent(in) :: field_name
    type(data_values_type), pointer :: get_data_value_from_hashmap_by_field_name

    class(*), pointer :: generic

    generic=>c_get_generic(collection, field_name)
    if (associated(generic)) then
      select type(generic)
        type is (data_values_type)
          get_data_value_from_hashmap_by_field_name=>generic
      end select
    else
      get_data_value_from_hashmap_by_field_name=>null()
    end if
  end function get_data_value_from_hashmap_by_field_name  

  !> Retrieves the data value (wrapper) by field name or null if no entry was found in the provided collection
  !! @param collection A map to search for this data value in
  !! @param field_name The field name to search for
  !! @returns The corresponding data value or null if none is found
  function get_data_value_from_map_by_field_name(collection, field_name)
    type(map_type), intent(inout) :: collection
    character(len=*), intent(in) :: field_name
    type(data_values_type), pointer :: get_data_value_from_map_by_field_name

    class(*), pointer :: generic

    generic=>c_get_generic(collection, field_name)
    if (associated(generic)) then
      select type(generic)
      type is (data_values_type)
        get_data_value_from_map_by_field_name=>generic
      end select
    else
      get_data_value_from_map_by_field_name=>null()
    end if
  end function get_data_value_from_map_by_field_name

  !> Retrieves the data value (wrapper) by field name or null if no entry was found in the provided map entry
  !! @param map_entry A map entry to convert into data wrapper
  !! @returns The corresponding data value or null if none is found
  function get_data_value_from_map_entry(map_entry)
    type(mapentry_type), intent(in) :: map_entry
    type(data_values_type), pointer :: get_data_value_from_map_entry

    class(*), pointer :: generic

    generic=>c_get_generic(map_entry)
    if (associated(generic)) then
      select type(generic)
      type is (data_values_type)
        get_data_value_from_map_entry=>generic
      end select
    else
      get_data_value_from_map_entry=>null()
    end if
  end function get_data_value_from_map_entry

  !> A helper function to get the location of a MONC's configuration in the IO data structure
  !! @param source Source index of the MONC process
  !! @returns Index that that MONC corresponds to
  integer function get_monc_location(io_configuration, source)
    type(io_configuration_type), intent(inout) :: io_configuration
    integer, intent(in) :: source

    get_monc_location=c_get_integer(io_configuration%monc_to_index, conv_to_string(source))
  end function get_monc_location

  !> Retrieves the diagnostics field configuration corresponding to a specific field name and returns whether one was found or not
  !! @param io_configuration The current IO server configuration
  !! @param field_name The name of the diagnostics field we are searching for
  !! @param diagnostic_config The found diagnostics is written into here if located
  !! @returns Whether a corresponding diagnostics field was found or not
  logical function get_diagnostic_field_configuration(io_configuration, field_name, field_namespace, diagnostic_config)
    type(io_configuration_type), intent(inout) :: io_configuration
    character(len=*), intent(in) :: field_name, field_namespace
    type(io_configuration_diagnostic_field_type), intent(out) :: diagnostic_config

    integer :: i

    do i=1, size(io_configuration%diagnostics)
      if (io_configuration%diagnostics(i)%name == field_name .and. &
           io_configuration%diagnostics(i)%namespace == field_namespace) then
        diagnostic_config=io_configuration%diagnostics(i)
        get_diagnostic_field_configuration=.true.
        return
      end if
    end do
    get_diagnostic_field_configuration=.false.
  end function get_diagnostic_field_configuration

  !> Retrieves the prognostic field configuration corresponding to a specific field name and returns whether one was found or not
  !! @param io_configuration The current IO server configuration
  !! @param field_name The name of the prognostics field we are searching for
  !! @param prognostic_config The found prognostic is written into here if located
  !! @returns Whether a corresponding prognostic field was found or not
  logical function get_prognostic_field_configuration(io_configuration, field_name, field_namespace, &
       prognostic_config, prognostic_containing_data_defn)
    type(io_configuration_type), intent(inout) :: io_configuration
    character(len=*), intent(in) :: field_name, field_namespace
    type(io_configuration_field_type), intent(out) :: prognostic_config
    type(io_configuration_data_definition_type), intent(out), optional :: prognostic_containing_data_defn

    integer :: i, j
    do i=1, io_configuration%number_of_data_definitions
      do j=1, io_configuration%data_definitions(i)%number_of_data_fields
        if (io_configuration%data_definitions(i)%fields(j)%name == field_name .and. &
             io_configuration%data_definitions(i)%fields(j)%namespace == field_namespace) then
          prognostic_config=io_configuration%data_definitions(i)%fields(j)
          if (present(prognostic_containing_data_defn)) then
            prognostic_containing_data_defn=io_configuration%data_definitions(i)
          end if
          get_prognostic_field_configuration=.true.
          return
        end if
      end do
    end do
    get_prognostic_field_configuration=.false.
  end function get_prognostic_field_configuration  

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
end module configuration_parser_mod
