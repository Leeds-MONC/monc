!> Parses a configuration file and loads the contents into the options database which can 
!! then be interogated by components in the model
module configuration_file_parser_mod
  use datadefn_mod, only : l_config_double
  use collections_mod, only : hashmap_type
  use conversions_mod, only : conv_to_logical, conv_to_integer, conv_to_real, &
       conv_is_logical, conv_is_integer, conv_is_real, conv_single_real_to_double, &
       string_to_double
  use optionsdatabase_mod, only : options_add, options_get_string, options_has_key, &
       options_get_array_size, options_remove_key
  use logging_mod, only : LOG_ERROR, log_master_log
  implicit none

#ifndef TEST_MODE
  private
#endif

  ! Ids used in the opening and reading in of configuration files
  integer, parameter :: USER_FILE_ID=15, GLOBAL_FILE_ID=16

  public parse_configuration_file
contains

  !> Parses a specific configuration and adds the contents into the options database
  !! @param options_database The options database
  !! @param user_configuration_file Filename of the initial (user) configuration file to open
  !!        and process
  subroutine parse_configuration_file(options_database, user_configuration_file)
    type(hashmap_type), intent(inout) :: options_database
    character(*), intent(in) :: user_configuration_file

    call process_configuration_file(options_database, user_configuration_file, .true., &
         USER_FILE_ID)
  end subroutine parse_configuration_file

  !> Will actually open a specific file and read it in line by line, parsing this and storing 
  !! the configuration options
  !! @param options_database The options database
  !! @param filename Name of file to open
  !! @param is_user_file Whether or not this is the user configuration file or not
  !! @param file_id The ID to use in the open call
  recursive subroutine process_configuration_file(options_database, filename, is_user_file,&
       file_id)
    type(hashmap_type), intent(inout) :: options_database
    character(*), intent(in) :: filename
    integer, intent(in) :: file_id
    logical, intent(in) :: is_user_file

    integer :: file_status, comment_point
    logical :: continue_parsing, found_global
    character(len=10000) :: raw_line

    found_global=.false.
    continue_parsing=.true.
    open(unit=file_id, file=filename, status='old', access='sequential', form='formatted',&
         action='read', iostat=file_status)
    if (file_status .ne. 0) then
      call log_master_log(LOG_ERROR, "Configuration file "//trim(filename)//" does not exist")
    end if

    do while (continue_parsing)
      read(file_id, '(A)', iostat=file_status) raw_line
      if (is_iostat_end(file_status)) then
        continue_parsing=.false.
      else
        !> Remove in-line commments.  If there's a comma in the comment, an array is assumed
        !  and errors will occur.  See has_multiple_values.
        !> Makes a section of store_configuration redundant.
        comment_point=index(raw_line,'#')
        if (comment_point .eq. 0) comment_point=index(raw_line,'!')
        if (comment_point .ne. 0) raw_line = raw_line(:comment_point - 1)
        
        !> Left adjust
        raw_line=adjustl(raw_line)
        if (len_trim(raw_line) .gt. 0) then
          call process_configuration_line(options_database, raw_line, is_user_file, &
               found_global)
        end if
      end if
    end do
    close(file_id)
  end subroutine process_configuration_file

  !> Processes a line from the configuration file, breaks it up into its key and value and 
  !! depending upon the specifics of the line it will store it in the options database in a
  !! variety of ways
  !! @param options_database The options database
  !! @param raw_line The raw text line to process
  !! @param is_user_file Whether this line was read from the user file or not
  !! @param found_global Whether or not we have found the global configuration file
  recursive subroutine process_configuration_line(options_database, raw_line, is_user_file, &
       found_global)
    type(hashmap_type), intent(inout) :: options_database
    character(*), intent(in) :: raw_line
    logical, intent(in) :: is_user_file
    logical, intent(inout) :: found_global
    
    integer :: mode, start_split, end_split
    character(len=10000) :: config_key, config_value

    call get_mode_and_split_points_from_line(raw_line, mode, start_split, end_split)

    if (mode .ge. 1 .and. raw_line(1:1) .ne. '#' .and. raw_line(1:1) .ne. '!') then
      config_key=raw_line(1:start_split)
      config_value=adjustl(raw_line(end_split:))     
      !> If multiple values exist in comma-separated format (a common event) or
      !  The key already exists with multiple values, treat as an array.
      !  The second case is important if, for instance, the global_config defaults
      !  to an array with multiple values, but the option may be acceptably entered
      !  as a single value.  This ensures the full key will be removed.
      if (has_multiple_values(config_value) .or.    &
          options_has_key(options_database, trim(config_key)//"a_size")) then
         call process_configuration_array(options_database, config_key, config_value, mode)        
      else
         if (is_key_array_index_specifier(config_key) .or. mode .eq. 2) then
            call handle_array_element_set(options_database, config_key, config_value, mode)
         else
            call store_configuration(options_database, config_key, config_value)
            if (is_user_file .and. .not. found_global) &
                 found_global=parse_global_configuration_if_available(options_database)
         end if
      end if
   end if
 end subroutine process_configuration_line

  !> Parses the global configuration file if it is available and calls on to add all of this
  !! to the options database
  !! @param options_database The options database
  !! @returns Whether the global configuration file was processed
  logical function parse_global_configuration_if_available(options_database)
    type(hashmap_type), intent(inout) :: options_database

    if (options_has_key(options_database, "global_configuration")) then
      parse_global_configuration_if_available=.true.
      call process_configuration_file(options_database, &
           trim(options_get_string(options_database, "global_configuration")), .false., &
           GLOBAL_FILE_ID)
    else
      parse_global_configuration_if_available=.false.
    end if
  end function parse_global_configuration_if_available

  !> Determines if a specific string contains multiple values such as str1, str2, str3
  !! @param configuration_value The string value
  !! @returns Whether or not the string is a vector of values (and needs to be split up)
  logical function has_multiple_values(configuration_value)
    character(*), intent(in) :: configuration_value

    has_multiple_values=scan(configuration_value, ",") .ne. 0
  end function has_multiple_values  

  !> Processes a line to determine the mode (replace or additive) and where the split point
  !! is between the key and value
  !! @param raw_line The raw configuration line to process
  !! @mode The mode, 1=replace, 2=additive and 0=none (unparsable)
  !! @param start_split The last element before the split, i.e. the key is up to and including 
  !!this element
  !! @param end_split The element after the split, i.e. the value is this key to the end
  subroutine get_mode_and_split_points_from_line(raw_line, mode, start_split, end_split)
    character(*), intent(in) :: raw_line
    integer, intent(out) :: mode, start_split, end_split

    integer :: split_point

    split_point=index(raw_line, "+=")
    if (split_point .eq. 0) split_point=index(raw_line, "=+")
    if (split_point .ne. 0) then
      mode=2
      start_split=split_point-1
      end_split=split_point+2
    else
      split_point=index(raw_line, "=")
      if (split_point .ne. 0) then
        mode=1
        start_split=split_point-1
        end_split=split_point+1
      else
        mode=0
      end if
    end if
  end subroutine get_mode_and_split_points_from_line

  !> Handles setting a specific array element, when the key has something like k(n) - 
  !! n being the index to set
  !! @param options_database The options database
  !! @param config_key The configuration key including the index
  !! @param config_value The configuration value to set
  !! @param mode Whether this is replace or additive
  subroutine handle_array_element_set(options_database, config_key, config_value, mode)
    type(hashmap_type), intent(inout) :: options_database
    character(*), intent(in) :: config_key, config_value
    integer, intent(in) :: mode

    integer :: array_index, key_end_index
    
    if (mode .eq. 2) then
      array_index=options_get_array_size(options_database, config_key)+1
      key_end_index=len(config_key)
    else
      array_index=get_key_array_index(config_key)
      key_end_index=scan(config_key,"(")-1
    end if

    call store_configuration(options_database, config_key(:key_end_index), &
         trim(adjustl(config_value)), array_index)
  end subroutine handle_array_element_set

  !> Given a configuration key of the form k(n), this returns the n
  !! @param config_key The configuration key to parse and extract the index from
  !! @returns The index encoded in the key
  integer function get_key_array_index(config_key)
    character(*), intent(in) :: config_key

    integer :: open_brace_index, close_brace_index

    open_brace_index=scan(config_key,"(")
    close_brace_index=scan(config_key,")")

    if (close_brace_index - open_brace_index .lt. 2) then
      call log_master_log(LOG_ERROR, "Array element format for key "//&
           trim(config_key)//" but no element provided")
    end if

    get_key_array_index=conv_to_integer(config_key(open_brace_index+1:close_brace_index-1))    
  end function get_key_array_index  

  !> Determines whether a configuration key represents a specific array element, i.e. is of 
  !! the form k(n)
  !! @param config_key The configuration key to check
  !! @returns Whether the key represents a specific element or not
  logical function is_key_array_index_specifier(config_key)
    character(*), intent(in) :: config_key

    integer :: loc

    loc=scan(config_key,"(")
    if (loc .ne. 0) then
      loc=scan(config_key,")")
      if (loc .ne. 0) then
        is_key_array_index_specifier=.true.
        return
      end if
    end if
    is_key_array_index_specifier=.false.
  end function is_key_array_index_specifier  

  !> Will process a configuration array of values such as v1,v2,v3,v4
  !! @param options_database The options database
  !! @param config_key The configuration key
  !! @param config_value The values each separated by comma
  !! @mode The mode to use for storage, 1=insert, 2=additive
  subroutine process_configuration_array(options_database, config_key, config_value, mode)
    type(hashmap_type), intent(inout) :: options_database
    character(*), intent(in) :: config_key, config_value
    integer, intent(in) :: mode

    character(len=len(config_value)) :: raw_value
    character(len=len(config_key)) :: parsed_config_key

    integer :: comma_posn, index

    if (mode == 1) then
      call options_remove_key(options_database, config_key)
      if (is_key_array_index_specifier(config_key)) then
        index=get_key_array_index(config_key)
        parsed_config_key=config_key(:scan(config_key,"(")-1)
      else
        parsed_config_key=config_key
        index=1
      end if
    else if (mode==2) then
      index=options_get_array_size(options_database, config_key)+1
    end if
    
    raw_value=config_value
    comma_posn=scan(raw_value, ",")
    do while (comma_posn .gt. 0)
      call store_configuration(options_database, parsed_config_key, &
           trim(adjustl(raw_value(1:comma_posn-1))), index)
      raw_value=raw_value(comma_posn+1:)
      comma_posn=scan(raw_value, ",")
      index=index+1
    end do
    call store_configuration(options_database, parsed_config_key, &
         trim(adjustl(raw_value(1:))), index)
  end subroutine process_configuration_array  

  !> Stores a specific configuration by determining the type of a value and calling on to 
  !! the options database for storage
  !! @param options_database The options database
  !! @param config_key The configuration key
  !! @param config_value The configuration value
  !! @param array_index The index in the array to insert the element into
  subroutine store_configuration(options_database, config_key, config_value, array_index)
    type(hashmap_type), intent(inout) :: options_database
    character(*), intent(in) :: config_key, config_value
    integer, intent(in), optional :: array_index

    integer :: comment_location
    character(len=len(config_value)) :: parsed_value

    comment_location=scan(config_value,"#")
    if (comment_location == 0) comment_location=scan(config_value,"!")
    if (comment_location .gt. 0) then
      parsed_value=config_value(:comment_location-1)
    else
      parsed_value=config_value
    end if

    if (conv_is_logical(trim(parsed_value))) then
      if (present(array_index)) then
        call options_add(options_database, trim(config_key), &
             conv_to_logical(trim(parsed_value)), &
             array_index=array_index)
      else
        call options_add(options_database, trim(config_key), &
             conv_to_logical(trim(parsed_value)))
      end if
    else if (conv_is_integer(parsed_value)) then
      if (present(array_index)) then
        call options_add(options_database, trim(config_key), &
             conv_to_integer(trim(parsed_value)), &
             array_index=array_index)
      else
        call options_add(options_database, trim(config_key), &
             conv_to_integer(trim(parsed_value)))
      end if
    else if (conv_is_real(parsed_value)) then
      if (present(array_index)) then
        if (.not. l_config_double) then
          call options_add(options_database, trim(config_key), &
                           conv_single_real_to_double(conv_to_real(trim(parsed_value))), &
                           array_index=array_index)
        else
          call options_add(options_database, trim(config_key), &
               string_to_double(trim(parsed_value)), array_index=array_index)
        end if
      else
        if (.not. l_config_double) then
          call options_add(options_database, trim(config_key), &
                           conv_single_real_to_double(conv_to_real(trim(parsed_value))))
        else    
          call options_add(options_database, trim(config_key), string_to_double(trim(parsed_value)))
        end if
      end if
    else
      if (present(array_index)) then
        call options_add(options_database, trim(config_key), &
             trim(remove_string_quotation(parsed_value)), array_index=array_index)
      else
        call options_add(options_database, trim(config_key), &
             trim(remove_string_quotation(parsed_value)))
      end if
    end if
  end subroutine store_configuration

  !> Removes quotations from a string if these are included, regardless of before it will
  !!  return the contents i.e. "abc" = abc, hg"abc"re = abc
  !! @param string_value The string value which may or may not contain string quotation marks
  !! @returns Either the quotation marks stripped out or the string unchanged if there were 
  !! no quotes
  function remove_string_quotation(string_value)
    character(len=*), intent(in) :: string_value
    character(len=len(string_value)) :: remove_string_quotation

    integer :: quotation_index_start, quotation_index_end

    quotation_index_start=scan(string_value, """")
    if (quotation_index_start .gt. 0) then
      quotation_index_end=scan(string_value(quotation_index_start+1:), """")+&
           quotation_index_start
      if (quotation_index_end .gt. 0) then
        remove_string_quotation=string_value(quotation_index_start+1:quotation_index_end-1)
      else
        remove_string_quotation=string_value
      end if
    else
      remove_string_quotation=string_value
    end if
  end function remove_string_quotation  
end module configuration_file_parser_mod
