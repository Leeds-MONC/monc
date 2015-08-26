!> Entry point for the IO server when run as a separate executable. This is currently used as a
!! driver for testing the aspects of our server
program iodriver
  use configuration_parser_mod, only : io_configuration_type, configuration_parse
  implicit none

  character(len=1000) :: xml, temp_line
  type(io_configuration_type) :: parsed_configuration
  integer :: k

  xml=""
  k=0
  open (unit = 2, file = "description.xml")
  do while (k == 0)
    read(2,"(A)",iostat=k) temp_line
    if (k == 0) then
      xml=trim(xml)//trim(temp_line)//new_line("A")
    end if    
  end do
  close(2)

  call configuration_parse(trim(xml), parsed_configuration)

  write(*,*) "Data fields", parsed_configuration%number_of_data_field, "Handling rules", parsed_configuration%number_of_rules
  do k=1,parsed_configuration%number_of_data_field
    write(*,*) "Field",k,parsed_configuration%fields(k)%name
  end do
  do k=1, parsed_configuration%number_of_rules
    write(*,*) "Rule ",k,parsed_configuration%rules(k)%name, "triggers",parsed_configuration%rules(k)%number_of_triggers, &
         "actions", parsed_configuration%rules(k)%number_of_actions
  end do
end program iodriver
