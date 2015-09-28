!> A SAX parser for XML files. This is used to parse the description of the data and rules.
!! Being a SAX parser it works in a callback fashion, as these often do, and will call the
!! start and end subroutines along with tag details. It still requires some additional
!! stability work to handle ill formatted XML and error checking.
module sax_xml_parser_mod
  use datadefn_mod, only : STRING_LENGTH
  use logging_mod, only : LOG_WARN, LOG_ERROR, log_log
  implicit none

#ifndef TEST_MODE
  private
#endif

  interface
     !> The start element callback interface (on opening of XML tag)
     subroutine start_element_callback_interface(element_name, number_of_attributes, attribute_names, attribute_values)
       character(len=*), intent(in) :: element_name       
       character(len=*), dimension(:), intent(in) :: attribute_names, attribute_values
       integer, intent(in) :: number_of_attributes
     end subroutine start_element_callback_interface

     !> The end element callback interface (on closing of XML tag, this is not called if an opening tag self closes)
     subroutine end_element_callback_interface(element_name)
       character(len=*), intent(in) :: element_name
     end subroutine end_element_callback_interface
  end interface
  
  public xml_parse
contains

  !> Parses some raw XML
  !! raw_contents The raw (unparsed) XML string
  !! start_element_callback Subroutine to call for each XML start element
  !! end_element_callback Subroutine to call for each XML end element
  subroutine xml_parse(raw_contents, start_element_callback, end_element_callback)
    character, dimension(:), intent(in) :: raw_contents
    procedure(start_element_callback_interface) :: start_element_callback
    procedure(end_element_callback_interface) :: end_element_callback

    character(len=size(raw_contents)) :: string_to_process

    integer :: current_index, start_index, end_index, i
    current_index=1

    ! Here we copy the raw contents array into a string so that string intrinsics can be used on it
    do i=1, size(raw_contents)
      string_to_process(i:i)=raw_contents(i)
    end do
    do while (current_index .lt. len(string_to_process))
      start_index=index(string_to_process(current_index:),"<")
      if (start_index .eq. 0) exit
      start_index=start_index+current_index-1
      end_index=index(string_to_process(start_index:),">")
      if (end_index .eq. 0) exit
      end_index=end_index+start_index-1
      call process_individual_tag(string_to_process, start_element_callback, end_element_callback, start_index, end_index)
      current_index=end_index
    end do    
  end subroutine xml_parse

  !> Processes an individual XML tag. This deduces whether it is a start or end tag, the name and any additional attributes
  !! @param raw_contents The raw XML string
  !! @param start_element_callback Subroutine to call on XML element opening tags
  !! @param end_element_callback Subroutine to call on XML element closing tags
  !! @param start_index The start index in the raw contents to go from (the start of this tag)
  !! @param end_index The end index in the raw contents to parse to (the end of this tag)
  subroutine process_individual_tag(raw_contents, start_element_callback, end_element_callback, start_index, end_index)
    character(len=*), intent(in) :: raw_contents
    procedure(start_element_callback_interface) :: start_element_callback
    procedure(end_element_callback_interface) :: end_element_callback
    integer, intent(in) :: start_index, end_index

    character(len=STRING_LENGTH) :: tag_name
    character(len=STRING_LENGTH), dimension(:), allocatable :: attribute_names, attribute_values
    logical :: start_tag
    integer :: name_start_index, name_end_index, number_attributes, attribute_index, i, attribute_start_posn, new_start_posn

    if (raw_contents(start_index+1:start_index+1) .eq. "!") return

    start_tag=.not. raw_contents(start_index+1:start_index+1) .eq. "/"
    name_start_index = start_index+1
    if (.not. start_tag) name_start_index=name_start_index+1
    name_end_index=index(raw_contents(name_start_index:end_index), " ")
    if (name_end_index .eq. 0) then
      name_end_index=end_index-1
    else
      name_end_index=name_end_index+name_start_index-1
    end if
      
    tag_name=raw_contents(name_start_index : name_end_index)    

    if (.not. start_tag) then
      call end_element_callback(tag_name)
    else
      number_attributes=occurances_of_substring(raw_contents(name_end_index+1:end_index), "=")
      attribute_index=1
      attribute_start_posn=1
      allocate(attribute_names(number_attributes), attribute_values(number_attributes))
      attribute_names(:)=""
      attribute_values(:)=""
      do i=1,number_attributes
        new_start_posn=get_attribute(raw_contents(name_end_index+1:end_index), attribute_start_posn, attribute_names, &
             attribute_values, i)
        attribute_start_posn=new_start_posn
      end do
      call start_element_callback(tag_name, number_attributes, attribute_names, attribute_values)
      deallocate(attribute_names, attribute_values)
      if (raw_contents(end_index-1:end_index) .eq. "/>") call end_element_callback(tag_name)
    end if
  end subroutine process_individual_tag

  !> Retrieves the "next" attribute from the XML tag and returns the position after this attribute to search from next. 
  !! This contains some additional complexities to deal with a number of different formatted tags
  !! @param contents The XML string
  !! @param start_index The index to start from when getting the next attribute
  !! @param attribute_names Collection of attribute names that the extraced name is appended to
  !! @param attribute_values Collection of attribute values that the extacted value is appended to
  !! @param attribute_index The index to write the name and value of this attribute into
  integer function get_attribute(contents, start_index, attribute_names, attribute_values, attribute_index)
    character(len=*), intent(in) :: contents
    integer, intent(in) :: start_index, attribute_index
    character(len=*), dimension(:), allocatable, intent(inout) :: attribute_names, attribute_values

    integer :: equals_posn, end_oftag_absolute, equals_absolute, open_quote, close_quote, begin_index, space_point

    close_quote=0
    equals_posn=index(contents(start_index:), "=")
    if (equals_posn .ne. 0) then
      ! Currently each attribute requires a specified value
      equals_absolute=equals_posn+start_index-1
      if (index(trim(adjustl(contents(start_index:equals_absolute))), " ") .ne. 0) then
        ! This warns of and eliminates any garbage before the attribute name, such as 't abc' sets 'abc' as the name
        call log_log(LOG_WARN, "Ignorning leading garbage in attribute name '"//contents(start_index:equals_absolute)//"'")
        begin_index=start_index+index(trim(adjustl(contents(start_index:equals_absolute))), " ")
      else 
        begin_index=start_index
      end if
      open_quote=index(contents(equals_absolute:), """")
      if (open_quote .ne. 0) close_quote=index(contents(equals_absolute+open_quote:), """")

      if (close_quote == 0) then
        ! No quote therefore must be quoteless and check for whitespace, also check for termination tag and use whichever is closed
        space_point=index(contents(equals_absolute+open_quote:), " ")        
        close_quote=index(contents(equals_absolute+open_quote:), "/>")-1
        if (close_quote .lt. 0) close_quote=index(contents(equals_absolute+open_quote:), ">")-1
        if (space_point .ne. 0 .and. space_point .lt. close_quote) close_quote=space_point
      else
        if (open_quote .gt. 2) then
          ! Deals with matching over to the next tag as no quotes around this value
          if (len(trim(adjustl(contents(equals_absolute+1:equals_absolute+open_quote-1)))) .gt. 0) then        
            close_quote=index(trim(adjustl(contents(equals_absolute+1:equals_absolute+open_quote-1))), " ")
            open_quote=0
          end if
        end if
      end if
      end_oftag_absolute=close_quote+equals_absolute+open_quote

      attribute_names(attribute_index)=trim(adjustl(contents(begin_index:equals_absolute-1)))
      attribute_values(attribute_index)=trim(adjustl(contents(equals_absolute+1:end_oftag_absolute-1)))
      if (len_trim(attribute_names(attribute_index)) == 0 .or. len_trim(attribute_values(attribute_index)) == 0) then
        call log_log(LOG_ERROR, "Empty IO server XML configuration name or value")
      end if
      get_attribute=end_oftag_absolute
    else
      get_attribute=0
    end if
  end function get_attribute  

  !> Returns the number of times a specific substring can be found in a string
  !! @param string The whole string to search
  !! @param substring The substring to search for and count the number occurances of
  integer function occurances_of_substring(string, substring)
    character(len=*), intent(in) :: string, substring

    integer :: current_index, found_index, sub_len

    sub_len=len(substring)
    occurances_of_substring=0
    current_index=1
    found_index=1

    do while (found_index .gt. 0)
      found_index = index(string(current_index:), substring)
      if (found_index .gt. 0) then
        occurances_of_substring=occurances_of_substring+1
        current_index=current_index+found_index+sub_len
      end if      
    end do    
  end function occurances_of_substring  
end module sax_xml_parser_mod
