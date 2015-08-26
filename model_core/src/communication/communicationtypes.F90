!> Contains the types used for communication, holding the state of communications and 
!! supporting activities. These are held in a separate module to allow for dependencies from 
!! multiple areas of the code base which might otherwise result in circular dependencies if 
!! they lived in a specific module of functionality
module communication_types_mod
  use datadefn_mod, only : DEFAULT_PRECISION
  implicit none

#ifndef TEST_MODE
  private
#endif

  ! A wrapper type (as F doesn't allow an array of 3D arrays) to point to field data
  type field_data_wrapper_type
     real(kind=DEFAULT_PRECISION), dimension(:,:,:), pointer :: data
  end type field_data_wrapper_type      

  !> Describes the neighbours of a process in a specific dimension and contains the 
  !! communication buffers associated with these
  type neighbour_description_type
     integer :: pid, halo_pages=0, halo_corners=0, dimension, recv_size, send_size, &
          recv_corner_size, send_corner_size
     real(kind=DEFAULT_PRECISION), dimension(:,:,:), allocatable :: send_halo_buffer, &
          recv_halo_buffer, send_corner_buffer, recv_corner_buffer
  end type neighbour_description_type

  !> Maintains the state of a halo swap and contains buffers, neighbours etc
  type halo_communication_type
     integer :: number_distinct_neighbours, fields_per_cell, halo_depth, cell_match(3,4)
     type(neighbour_description_type), dimension(:), allocatable :: halo_swap_neighbours     
     integer, dimension(:), allocatable :: send_requests, recv_requests
     logical :: initialised = .false., swap_in_progress=.false., involve_corners=.false.
  end type halo_communication_type

  public halo_communication_type, neighbour_description_type, field_data_wrapper_type
end module communication_types_mod
