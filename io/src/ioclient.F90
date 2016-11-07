!> This defines some constants and procedures that are useful to the IO server and clients that call it. By using the contents
!! then the client can guarantee consistency against what the server expects
module io_server_client_mod
  use datadefn_mod, only : DEFAULT_PRECISION, SINGLE_PRECISION, DOUBLE_PRECISION, STRING_LENGTH
  use collections_mod, only : hashmap_type, iterator_type, mapentry_type, c_get_iterator, c_has_next, &
       c_next_mapentry, c_get_generic
  use conversions_mod, only : conv_to_string
  use mpi, only : MPI_CHARACTER, MPI_INT, MPI_LOGICAL, MPI_REAL, MPI_DOUBLE_PRECISION, MPI_ADDRESS_KIND
  implicit none

#ifndef TEST_MODE
  private
#endif

  !< Data structure used to hold a sizing description of a field
  type data_sizing_description_type
     character(len=STRING_LENGTH) :: field_name !< Name of the field that this describes
     integer :: dimensions, dim_sizes(4) !< The number of dimensions and size in each dimension
  end type data_sizing_description_type

  type field_description_type
     character(len=STRING_LENGTH) :: definition_name, field_name
     integer :: field_type, data_type
     logical :: optional
  end type field_description_type

  type definition_description_type
     character(len=STRING_LENGTH) :: definition_name
     logical :: send_on_terminate
     integer :: number_fields, frequency
  end type definition_description_type  

  ! Constants used in sending and receiving IO data
  integer, parameter :: COMMAND_TAG=9, DATA_TAG=10, REGISTER_COMMAND=1, DEREGISTER_COMMAND=3, DATA_COMMAND_START=4, &
       INTER_IO_COMMUNICATION=2

  !< Field type identifiers
  integer, parameter :: SCALAR_FIELD_TYPE = 1, ARRAY_FIELD_TYPE=2, MAP_FIELD_TYPE=3
  !< Field data type identifiers
  integer, parameter :: INTEGER_DATA_TYPE = 1, BOOLEAN_DATA_TYPE=2, STRING_DATA_TYPE=3, FLOAT_DATA_TYPE=4, &
       DOUBLE_DATA_TYPE=5

  character(len=STRING_LENGTH), parameter :: LOCAL_SIZES_KEY="local_sizes", LOCAL_START_POINTS_KEY="local_start_points", &
       LOCAL_END_POINTS_KEY="local_end_points", NUMBER_Q_INDICIES_KEY="num_q_indicies"

  public COMMAND_TAG, DATA_TAG, REGISTER_COMMAND, DEREGISTER_COMMAND, DATA_COMMAND_START, INTER_IO_COMMUNICATION, &
       SCALAR_FIELD_TYPE, ARRAY_FIELD_TYPE, MAP_FIELD_TYPE, INTEGER_DATA_TYPE, BOOLEAN_DATA_TYPE, &
       STRING_DATA_TYPE, FLOAT_DATA_TYPE, LOCAL_SIZES_KEY, LOCAL_START_POINTS_KEY, LOCAL_END_POINTS_KEY, NUMBER_Q_INDICIES_KEY, &
       DOUBLE_DATA_TYPE, data_sizing_description_type, populate_mpi_type_extents, append_mpi_datatype, &
       get_mpi_datatype_from_internal_representation, definition_description_type, field_description_type, &
       build_mpi_type_data_sizing_description, build_mpi_type_field_description, build_mpi_type_definition_description, &
       pack_scalar_field, pack_array_field, pack_map_field, get_data_description_from_name

contains

  !> Provides the type extents of the types that we are using in construction of the MPI data type. This is done
  !! once for fast look up in the actual construction phase.
  !! @returns An array of type extents keyed by the DATA_TYPE constants of the configuration parser module)
  function populate_mpi_type_extents()
    integer :: populate_mpi_type_extents(5)

    integer :: ierr
    integer(kind=MPI_ADDRESS_KIND) :: large_number_extents(5)
    
    call mpi_type_extent(MPI_INT, large_number_extents(INTEGER_DATA_TYPE), ierr)
    call mpi_type_extent(MPI_LOGICAL, large_number_extents(BOOLEAN_DATA_TYPE), ierr)
    call mpi_type_extent(MPI_CHARACTER, large_number_extents(STRING_DATA_TYPE), ierr)
    call mpi_type_extent(MPI_REAL, large_number_extents(FLOAT_DATA_TYPE), ierr)
    call mpi_type_extent(MPI_DOUBLE_PRECISION, large_number_extents(DOUBLE_DATA_TYPE), ierr)

    populate_mpi_type_extents=int(large_number_extents)
  end function populate_mpi_type_extents  

  !> Builds the MPI data type for sending data descriptions to registree MONCs
  integer function build_mpi_type_definition_description()
    integer :: new_type, ierr, block_counts(4), old_types(4), offsets(4)
    integer(MPI_ADDRESS_KIND) :: num_addr, base_addr

    type(definition_description_type) :: basic_type

    call mpi_get_address(basic_type, base_addr, ierr)
    old_types(1) = MPI_CHARACTER
    block_counts(1) = STRING_LENGTH
    offsets(1)=0

    call mpi_get_address(basic_type%send_on_terminate, num_addr, ierr)
    old_types(2) = MPI_LOGICAL
    block_counts(2) = 1
    offsets(2)=int(num_addr-base_addr)

    call mpi_get_address(basic_type%number_fields, num_addr, ierr)
    old_types(3) = MPI_INT
    block_counts(3) = 1
    offsets(3)=int(num_addr-base_addr)

    call mpi_get_address(basic_type%frequency, num_addr, ierr)
    old_types(4) = MPI_INT
    block_counts(4) = 1
    offsets(4)=int(num_addr-base_addr)

    call mpi_type_struct(4, block_counts, offsets, old_types, new_type, ierr) 
    call mpi_type_commit(new_type, ierr)
    build_mpi_type_definition_description=new_type
  end function build_mpi_type_definition_description

  !> Builds the MPI data type for sending field descriptions to registree MONCs
  integer function build_mpi_type_field_description()
    integer :: new_type, ierr, old_types(5), block_counts(5), offsets(5)
    integer(MPI_ADDRESS_KIND) :: num_addr, base_addr

    type(field_description_type) :: basic_type

    call mpi_get_address(basic_type, base_addr, ierr)    
    old_types(1) = MPI_CHARACTER
    block_counts(1) = STRING_LENGTH
    offsets(1)=0

    call mpi_get_address(basic_type%field_name, num_addr, ierr)    
    old_types(2) = MPI_CHARACTER
    block_counts(2) = STRING_LENGTH
    offsets(2)=int(num_addr-base_addr)

    call mpi_get_address(basic_type%field_type, num_addr, ierr)
    old_types(3) = MPI_INT
    block_counts(3) = 1
    offsets(3)=int(num_addr-base_addr)
    
    call mpi_get_address(basic_type%data_type, num_addr, ierr)
    old_types(4) = MPI_INT
    block_counts(4) = 1
    offsets(4)=int(num_addr-base_addr)

    call mpi_get_address(basic_type%optional, num_addr, ierr)
    old_types(5) = MPI_LOGICAL
    block_counts(5) = 1
    offsets(5)=int(num_addr-base_addr)

    call mpi_type_struct(5, block_counts, offsets, old_types, new_type, ierr) 
    call mpi_type_commit(new_type, ierr)
    build_mpi_type_field_description=new_type
  end function build_mpi_type_field_description

  !> Builds the MPI type used for sending to the IO server a description of the data, namely the size
  !! of the arrays on this process
  !! @return The handle of the MPI type
  integer function build_mpi_type_data_sizing_description()
    integer :: new_type, ierr, block_counts(3), old_types(3), offsets(3)
    integer(kind=MPI_ADDRESS_KIND) :: num_addr, base_addr

    type(data_sizing_description_type) :: basic_type

    call mpi_get_address(basic_type, base_addr, ierr)
    old_types(1) = MPI_CHARACTER
    block_counts(1) = STRING_LENGTH 
    offsets(1)=0

    call mpi_get_address(basic_type%dimensions, num_addr, ierr)
    old_types(2) = MPI_INT
    block_counts(2) = 1
    offsets(2)=int(num_addr-base_addr)

    call mpi_get_address(basic_type%dim_sizes, num_addr, ierr)    
    old_types(3) = MPI_INT
    block_counts(3) = 4
    offsets(3)=int(num_addr-base_addr)

    call mpi_type_struct(3, block_counts, offsets, old_types, new_type, ierr) 
    call mpi_type_commit(new_type, ierr)
    build_mpi_type_data_sizing_description=new_type
  end function build_mpi_type_data_sizing_description

  !> Appends the MPI datatype details to the block counts, old types and offsets arrays. This will lump together multiple
  !! concurrent fields with the same type
  !! @param field_start Number of field that we are starting at to build
  !! @param field_end Number of field that we are ending at in this build
  !! @param field_array_sizes The extra size that is taken up not just by the data but the fact that some members are arrays (array size)
  !! @param data_type IO server representation of the current field data type
  !! @param type_extents Array of MPI type extents
  !! @param prev_data_type IO server representation of the previous distinct field data type
  !! @param type_index The number in the MPI type structure that we are creating here
  !! @param old_types Array of MPI types that we use for each member in the structure
  !! @param offsets Starting offsets for each MPI type structure element
  !! @param block_counts Number of type blocks for each MPI type structure element
  subroutine append_mpi_datatype(field_start, field_end, field_array_sizes, data_type, &
       type_extents, prev_data_type, type_index, old_types, offsets, block_counts)
    integer, intent(in) :: field_start, field_end, field_array_sizes, data_type, type_index, prev_data_type, type_extents(5)
    integer, intent(inout) :: old_types(20), offsets(20), block_counts(20)

    block_counts(type_index)=(field_end-field_start) + 1 + field_array_sizes
    old_types(type_index)=get_mpi_datatype_from_internal_representation(data_type)
    if (type_index == 1) then
      offsets(1)=0
    else
      offsets(type_index)=offsets(type_index-1)+type_extents(prev_data_type) * block_counts(type_index-1)
    end if    
  end subroutine append_mpi_datatype

  !> Gets the MPI datatype from out internal representation of the field data type (as in the configuration
  !! parse module)
  !! @param Our internal type code representation of the data type
  !! @returns The corresponding MPI data type
  integer function get_mpi_datatype_from_internal_representation(type_code)
    integer, intent(in) :: type_code

    if (type_code==INTEGER_DATA_TYPE) then
      get_mpi_datatype_from_internal_representation=MPI_INT
    else if (type_code==BOOLEAN_DATA_TYPE) then
      get_mpi_datatype_from_internal_representation=MPI_LOGICAL
    else if (type_code==STRING_DATA_TYPE) then
      get_mpi_datatype_from_internal_representation=MPI_CHARACTER
    else if (type_code==FLOAT_DATA_TYPE) then
      get_mpi_datatype_from_internal_representation=MPI_REAL
    else if (type_code==DOUBLE_DATA_TYPE) then
      get_mpi_datatype_from_internal_representation=MPI_DOUBLE_PRECISION
    end if
  end function get_mpi_datatype_from_internal_representation

    !> Packs a map into the send buffer
  !! @param buffer The buffer to pack the field into
  !! @param start_offset The starting offset to write into the buffer
  !! @param map_to_pack The map to pack into the buffer
  !! @returns The next location in the buffer to write to (next start offset)
  integer function pack_map_field(buffer, start_offset, map_to_pack)
    character, dimension(:), intent(inout) :: buffer
    integer, intent(in) :: start_offset
    type(hashmap_type) :: map_to_pack

    integer :: i, target_end, current_offset
    character(len=STRING_LENGTH) :: temp_string
    character(len=STRING_LENGTH), pointer :: sized_raw_character
    class(*), pointer :: raw_data, raw_to_string
    type(iterator_type) :: map_iterator
    type(mapentry_type) :: specific_mapentry

    current_offset=start_offset
    map_iterator=c_get_iterator(map_to_pack)
    do while (c_has_next(map_iterator))
      specific_mapentry=c_next_mapentry(map_iterator)
      temp_string=specific_mapentry%key
      target_end=current_offset+STRING_LENGTH-1
      buffer(current_offset:target_end)=transfer(temp_string, buffer(current_offset:target_end))
      current_offset=target_end+1

      raw_data=>c_get_generic(specific_mapentry)
      raw_to_string=>raw_data
      select type (raw_data)
      type is(integer)
        temp_string=conv_to_string(raw_data)
      type is(real(4))
        temp_string=conv_to_string(raw_data)
      type is (real(8))
        temp_string=conv_to_string(raw_data)
      type is(logical)
        temp_string=conv_to_string(raw_data)
      type is(character(len=*))
        sized_raw_character=>conv_to_string(raw_to_string, .false., STRING_LENGTH)
        temp_string=sized_raw_character
      end select
      target_end=current_offset+STRING_LENGTH-1
      buffer(current_offset:target_end)=transfer(temp_string, buffer(current_offset:target_end))
      current_offset=target_end+1
    end do
    pack_map_field=current_offset
  end function pack_map_field

  !> Packs an array field into the sending buffer
  !! @param buffer The buffer to pack the field into
  !! @param start_offset The starting offset to write into the buffer
  !! @param int_value (Optional) integer array values to pack
  !! @param real_value (Optional) default precision real array values to pack
  !! @returns The next location in the buffer to write to (next start offset)
  integer function pack_array_field(buffer, start_offset, int_array, real_array_1d, real_array_2d, real_array_3d, real_array_4d)
    character, dimension(:), intent(inout) :: buffer
    integer, intent(in) :: start_offset
    integer, dimension(:), intent(in), optional :: int_array
    real(kind=DEFAULT_PRECISION), dimension(:), intent(in), optional :: real_array_1d
    real(kind=DEFAULT_PRECISION), dimension(:,:), intent(in), optional :: real_array_2d
    real(kind=DEFAULT_PRECISION), dimension(:,:,:), intent(in), optional :: real_array_3d
    real(kind=DEFAULT_PRECISION), dimension(:,:,:,:), intent(in), optional :: real_array_4d

    integer :: target_end

    if (present(int_array)) then
      target_end=start_offset+kind(int_array)*size(int_array)-1    
      buffer(start_offset:target_end) = transfer(int_array, buffer(start_offset:target_end))
    else if (present(real_array_1d)) then
      target_end=start_offset+kind(real_array_1d)*size(real_array_1d)-1    
      buffer(start_offset:target_end) = transfer(real_array_1d, buffer(start_offset:target_end))
    else if (present(real_array_2d)) then
      target_end=start_offset+kind(real_array_2d)*size(real_array_2d)-1    
      buffer(start_offset:target_end) = transfer(real_array_2d, buffer(start_offset:target_end))
    else if (present(real_array_3d)) then
      target_end=start_offset+kind(real_array_3d)*size(real_array_3d)-1    
      buffer(start_offset:target_end) = transfer(real_array_3d, buffer(start_offset:target_end))
    else if (present(real_array_4d)) then
      target_end=start_offset+kind(real_array_4d)*size(real_array_4d)-1    
      buffer(start_offset:target_end) = transfer(real_array_4d, buffer(start_offset:target_end))
    end if
    pack_array_field=target_end+1
  end function pack_array_field

  !> Packs the data of a scalar field into a buffer
  !! @param buffer The buffer to pack the field into
  !! @param start_offset The starting offset to write into the buffer
  !! @param int_value (Optional) integer scalar value to pack
  !! @param real_value (Optional) default precision real scalar value to pack
  !! @param single_real_value (Optional) single precision real scalar value to pack
  !! @param double_real_value (Optional) double precision real scalar value to pack
  !! @returns The next location in the buffer to write to (next start offset)
  integer function pack_scalar_field(buffer, start_offset, int_value, real_value, single_real_value, double_real_value, &
       string_value, logical_value)
    character, dimension(:), intent(inout) :: buffer
    integer, intent(in) :: start_offset
    integer, intent(in), optional :: int_value
    real(kind=DEFAULT_PRECISION), intent(in), optional :: real_value
    real(kind=SINGLE_PRECISION), intent(in), optional :: single_real_value
    real(kind=DOUBLE_PRECISION), intent(in), optional :: double_real_value
    character(len=*), intent(in), optional :: string_value
    logical, intent(in), optional :: logical_value

    integer :: target_end
    character(len=STRING_LENGTH) :: string_to_insert

    if (present(int_value)) then
      target_end=start_offset+kind(int_value)-1
      buffer(start_offset:target_end) = transfer(int_value, buffer(start_offset:target_end))
    else if (present(real_value)) then
      target_end=start_offset+kind(real_value)-1
      buffer(start_offset:target_end) = transfer(real_value, buffer(start_offset:target_end))
    else if (present(single_real_value)) then
      target_end=start_offset+kind(single_real_value)-1
      buffer(start_offset:target_end) = transfer(single_real_value, buffer(start_offset:target_end))
    else if (present(double_real_value)) then
      target_end=start_offset+kind(double_real_value)-1
      buffer(start_offset:target_end) = transfer(double_real_value, buffer(start_offset:target_end))
    else if (present(string_value)) then
      target_end=start_offset+STRING_LENGTH-1
      string_to_insert=string_value
      buffer(start_offset:target_end) = transfer(string_to_insert, buffer(start_offset:target_end))
    else if (present(logical_value)) then
      target_end=start_offset+kind(logical_value)-1
      buffer(start_offset:target_end) = transfer(logical_value, buffer(start_offset:target_end))
    else 
      target_end=start_offset-1
    end if
    pack_scalar_field=target_end+1
  end function pack_scalar_field

  !> Look up the data description that corresponds to a specific field keyed by its name
  !! @param descriptions The data descriptions
  !! @param name The field name to look up
  !! @param field_description The resulting field description if the field is found
  !! @returns Whether or not the field, based upon its name lookup, was found
  logical function get_data_description_from_name(descriptions, name, field_description)
    type(data_sizing_description_type), dimension(:), intent(in) :: descriptions
    type(data_sizing_description_type), intent(out), optional :: field_description
    character(len=*), intent(in) :: name

    integer :: i
    do i=1,size(descriptions)
      if (descriptions(i)%field_name == name) then
        if (present(field_description)) field_description=descriptions(i)
        get_data_description_from_name=.true.
        return
      end if
    end do
    get_data_description_from_name=.false.
  end function get_data_description_from_name
end module io_server_client_mod
