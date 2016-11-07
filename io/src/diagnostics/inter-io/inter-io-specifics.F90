!> Inter IO server communication specific functionality. This manages all of the communication that might happen between
!! different IO servers.
module inter_io_specifics_mod
  use datadefn_mod, only : DEFAULT_PRECISION, STRING_LENGTH
  use configuration_parser_mod, only : DATA_SIZE_STRIDE, handle_recv_data_from_io_server, io_configuration_type, &
       io_configuration_data_definition_type, io_configuration_inter_communication_description, extend_inter_io_comm_array       
  use io_server_client_mod, only : data_sizing_description_type, pack_scalar_field, pack_array_field
  use forthread_mod, only : forthread_rwlock_rdlock, forthread_rwlock_wrlock, forthread_rwlock_unlock, &
       forthread_rwlock_init, forthread_rwlock_destroy, forthread_mutex_init
  use threadpool_mod, only : check_thread_status
  implicit none

#ifndef TEST_MODE
  private
#endif

  abstract interface
     subroutine handle_completion(io_configuration, values, field_name, timestep)
       import DEFAULT_PRECISION, STRING_LENGTH, io_configuration_type
       type(io_configuration_type), intent(inout) :: io_configuration
       real(DEFAULT_PRECISION), dimension(:) :: values
       character(len=STRING_LENGTH) :: field_name
       integer :: timestep
     end subroutine handle_completion
  end interface

  integer :: starting_tag=20

  public handle_completion, register_inter_io_communication, package_inter_io_communication_message, &
       unpackage_inter_io_communication_message, find_inter_io_from_name
contains  

  !> Registers an inter IO communication operation
  !! @param io_configuration IO server configuration state
  !! @param message_tag The inter IO message tag
  !! @param handling_procedure Callback procedure back to the inter IO functionality when a message arrives
  !! @param name The name of the inter IO operation
  subroutine register_inter_io_communication(io_configuration, message_tag, handling_procedure, name)
    type(io_configuration_type), intent(inout) :: io_configuration
    integer, intent(in) :: message_tag
    procedure(handle_recv_data_from_io_server) :: handling_procedure
    character(len=*), intent(in) :: name

    integer :: new_index

    io_configuration%number_inter_io_communications=io_configuration%number_inter_io_communications+1

    if (io_configuration%number_inter_io_communications .gt. size(io_configuration%inter_io_communications)) then        
      call extend_inter_io_comm_array(io_configuration)
    end if

    new_index=io_configuration%number_inter_io_communications
    io_configuration%inter_io_communications(new_index)%message_tag=message_tag+starting_tag
    io_configuration%inter_io_communications(new_index)%handling_procedure=>handling_procedure
    io_configuration%inter_io_communications(new_index)%name=name
  end subroutine register_inter_io_communication

  !> Locates a the index of an inter IO entry from the operator name or returns 0 if none is found
  !! @param io_configuration The configuration of the IO server
  !! @param name Name of the inter IO communication operation
  !! @returns The index of the inter IO descriptor or 0 if none is found
  integer function find_inter_io_from_name(io_configuration, name)
    type(io_configuration_type), intent(inout) :: io_configuration
    character(len=*), intent(in) :: name

    integer :: i

    do i=1, io_configuration%number_inter_io_communications
      if (io_configuration%inter_io_communications(i)%name .eq. name) then
        find_inter_io_from_name=i
        return
      end if      
    end do
    find_inter_io_from_name=0
  end function find_inter_io_from_name  

  !> Packages up fields into an io binary message (allocated here) which is used for sending
  !! @param field_name The field name to package
  !! @param timestep The timestep to package
  !! @param field_values The field values to package
  !! @param other_int Optional other integer to package
  !! @returns The binary data ready for communication
  function package_inter_io_communication_message(field_name, timestep, field_values, other_int)
    character(len=STRING_LENGTH), intent(in) :: field_name
    integer, intent(in) :: timestep
    integer, intent(in), optional :: other_int
    real(kind=DEFAULT_PRECISION), dimension(:), intent(in) :: field_values
    character, dimension(:), allocatable :: package_inter_io_communication_message

    integer :: current_location
    
    allocate(package_inter_io_communication_message(STRING_LENGTH+8+(8*size(field_values))))
    
    current_location=pack_scalar_field(package_inter_io_communication_message, 1, string_value=field_name)
    current_location=pack_scalar_field(package_inter_io_communication_message, current_location, int_value=timestep)
    if (present(other_int)) then
      current_location=pack_scalar_field(package_inter_io_communication_message, current_location, int_value=other_int)
    end if
    current_location=pack_array_field(package_inter_io_communication_message, current_location, real_array_1d=field_values)
  end function package_inter_io_communication_message

  !> Unpackages some binary data into its individual fields. The field values are allocated here and the size is the remaining
  !! message size after everything else has been accounted for
  !! @param data_buffer The data buffer to unpackage
  !! @param field_name Written into from the buffer
  !! @param timestep Written into from the buffer
  !! @param field_values Remaining message size, allocated here
  !! @param other_int Optional other integer written in from the buffer
  subroutine unpackage_inter_io_communication_message(data_buffer, field_name, timestep, field_values, other_int)
    character, dimension(:), intent(in) :: data_buffer
    character(len=STRING_LENGTH), intent(out) :: field_name
    integer, intent(out) :: timestep
    integer, intent(out), optional :: other_int
    real(kind=DEFAULT_PRECISION), dimension(:), allocatable, intent(out) :: field_values
    
    integer :: values_entries

    values_entries=(size(data_buffer) - (STRING_LENGTH+8))/8

    allocate(field_values(values_entries))
    field_name=transfer(data_buffer(1:STRING_LENGTH), field_name)
    timestep=transfer(data_buffer(STRING_LENGTH+1:STRING_LENGTH+4), timestep)
    if (present(other_int)) other_int=transfer(data_buffer(STRING_LENGTH+5:STRING_LENGTH+8), other_int)
    field_values=transfer(data_buffer(STRING_LENGTH+9:), field_values)
  end subroutine unpackage_inter_io_communication_message  
end module inter_io_specifics_mod
