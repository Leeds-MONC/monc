!> Performs halo swapping. In the parallel case this is between neighbouring processes
!! and in the serial case it still needs to wrap the halos around for the boundary conditions. This module
!! determines the policy of halo swapping (i.e. the fields to communicate) and the halo communication
!! module is used to provide the actual mechanism.
module haloswapper_mod
  use monc_component_mod, only : component_descriptor_type
  use state_mod, only : model_state_type
  use prognostics_mod, only : prognostic_field_type
  use grids_mod, only : local_grid_type, X_INDEX, Y_INDEX, Z_INDEX
  use logging_mod, only : LOG_DEBUG, log_get_logging_level, log_log
  use conversions_mod, only : conv_to_string
  use datadefn_mod, only : DEFAULT_PRECISION, PRECISION_TYPE
  use communication_types_mod, only : halo_communication_type, neighbour_description_type, &
       field_data_wrapper_type
  use halo_communication_mod, only : copy_buffer_to_field, copy_field_to_buffer, &
       perform_local_data_copy_for_field,&
       init_halo_communication, finalise_halo_communication, blocking_halo_swap, &
       copy_corner_to_buffer, copy_buffer_to_corner
  use mpi, only : MPI_REQUEST_NULL, MPI_STATUSES_IGNORE
  use optionsdatabase_mod, only : options_get_integer
  implicit none

#ifndef TEST_MODE
  private
#endif  
  !< First call tracked for use in displaying debugging information only once
  logical, private :: first_call = .true. 
  type(halo_communication_type), save :: halo_swap_state

  public haloswapper_get_descriptor

contains

  !> Provides a description of this component for the core to register
  !! @returns The descriptor containing registration information for this component
  type(component_descriptor_type) function haloswapper_get_descriptor()
    haloswapper_get_descriptor%name="halo_swapper"
    haloswapper_get_descriptor%version=0.1
    haloswapper_get_descriptor%initialisation=>initialisation_callback
    haloswapper_get_descriptor%timestep=>timestep_callback
    haloswapper_get_descriptor%finalisation=>finalisation_callback
  end function haloswapper_get_descriptor  

  !> Initialisation callback hook which will set up the halo swapping state and cache some 
  !!precalculated data for fast(er) halo swapping at each timestep
  !! @param current_state The current model state
  subroutine initialisation_callback(current_state)
    type(model_state_type), target, intent(inout) :: current_state
    integer :: halo_depth

    ! get halo_depth and pass it to the halo_swapping routines
    halo_depth = options_get_integer(current_state%options_database, "halo_depth")
    call init_halo_communication(current_state, get_fields_per_halo_cell, halo_swap_state, &
         halo_depth, .true.)
  end subroutine initialisation_callback  

  !> Timestep callback hook which performs the halo swapping for each prognostic field
  !!
  !! In parallel this is performed with MPI communication calls and wrapping around. In serial
  !! still need to wrap data around
  !! @param current_state The current model state_mod
  subroutine timestep_callback(current_state)
    type(model_state_type), target, intent(inout) :: current_state

    call blocking_halo_swap(current_state, halo_swap_state, copy_fields_to_halo_buffer, &
         perform_local_data_copy_for_all_prognostics, copy_halo_buffer_to_field, &
         copy_corners_to_halo_buffer, copy_halo_buffer_to_corners)    
    call display_debugging_info_if_needed(halo_swap_state%number_distinct_neighbours, &
         halo_swap_state%fields_per_cell, current_state%parallel%my_rank)

  end subroutine timestep_callback
  
  !> The finalisation callback hook which will clean up and free the memory associated with the 
  !! halo swapping
  !! @param current_state The current model state
  subroutine finalisation_callback(current_state)
    type(model_state_type), target, intent(inout) :: current_state

    call finalise_halo_communication(halo_swap_state)
  end subroutine finalisation_callback 

  !> Copies the halo buffer to halo location in a corner for a halo cell/column and corner 
  !! location.
  !! The copies are performed for each prognostic field.
  !! @param current_state The current model state
  !! @param neighbour_description The halo swapping description of the neighbour we are 
  !!                              accessing the buffer of
  !! @param corner_loc The location of the corner
  !! @param x_target_index The target index for the x dimension we are receiving for
  !! @param y_target_index The target index for the y dimension we are receiving for
  !! @param neighbour_location The location in the local neighbour data stores of this neighbour
  !! @param current_page The current, next, halo swap page to read from (all previous have been
  !!                     read and copied already)
  !! @param source_data Optional source data which is read from into send buffers and written
  !!                    into by receieve buffers
  subroutine copy_halo_buffer_to_corners(current_state, neighbour_description, corner_loc, &
       x_target_index, y_target_index, neighbour_location, current_page, source_data)
    
    type(model_state_type), intent(inout) :: current_state
    integer, intent(in) :: corner_loc, x_target_index, y_target_index, neighbour_location
    integer, intent(inout) :: current_page(:)
    type(neighbour_description_type), intent(inout) :: neighbour_description
    type(field_data_wrapper_type), dimension(:), intent(in), optional :: source_data

    integer :: page_bookmark, i

    page_bookmark = current_page(neighbour_location)
#ifdef U_ACTIVE
    call copy_buffer_to_corner(current_state%local_grid, neighbour_description%recv_corner_buffer,&
         current_state%u%data, corner_loc, x_target_index, y_target_index, page_bookmark)
    page_bookmark=page_bookmark+1
    call copy_buffer_to_corner(current_state%local_grid, neighbour_description%recv_corner_buffer, &
         current_state%zu%data, corner_loc, x_target_index, y_target_index, page_bookmark)
    page_bookmark=page_bookmark+1
#endif
#ifdef V_ACTIVE
    call copy_buffer_to_corner(current_state%local_grid, neighbour_description%recv_corner_buffer, &
         current_state%v%data, corner_loc, x_target_index, y_target_index, page_bookmark)
    page_bookmark=page_bookmark+1
    call copy_buffer_to_corner(current_state%local_grid, neighbour_description%recv_corner_buffer, &
         current_state%zv%data, corner_loc, x_target_index, y_target_index, page_bookmark)
    page_bookmark=page_bookmark+1
#endif
#ifdef W_ACTIVE
    call copy_buffer_to_corner(current_state%local_grid, neighbour_description%recv_corner_buffer, &
         current_state%w%data, corner_loc, x_target_index, y_target_index, page_bookmark)
    page_bookmark=page_bookmark+1
    call copy_buffer_to_corner(current_state%local_grid, neighbour_description%recv_corner_buffer, &
         current_state%zw%data, corner_loc, x_target_index, y_target_index, page_bookmark)
    page_bookmark=page_bookmark+1
#endif
    if (current_state%th%active) then 
      call copy_buffer_to_corner(current_state%local_grid, neighbour_description%recv_corner_buffer, &
           current_state%th%data, corner_loc, x_target_index, y_target_index, page_bookmark)
      page_bookmark=page_bookmark+1
      call copy_buffer_to_corner(current_state%local_grid, neighbour_description%recv_corner_buffer, &
           current_state%zth%data, corner_loc, x_target_index, y_target_index, page_bookmark)
      page_bookmark=page_bookmark+1
    end if
    do i=1,current_state%number_q_fields
      if (current_state%q(i)%active) then
        call copy_buffer_to_corner(current_state%local_grid, neighbour_description%recv_corner_buffer, &
             current_state%q(i)%data, corner_loc, x_target_index, y_target_index, page_bookmark)
        page_bookmark=page_bookmark + 1
        call copy_buffer_to_corner(current_state%local_grid, neighbour_description%recv_corner_buffer, &
             current_state%zq(i)%data, corner_loc, x_target_index, y_target_index, page_bookmark)
        page_bookmark=page_bookmark + 1        
      end if
    end do
    current_page(neighbour_location)=page_bookmark
  end subroutine copy_halo_buffer_to_corners

  !> Copies the halo buffer to halo location in a field for a specific dimension and halo cell/column.
  !! The copies are performed for each prognostic field.
  !! @param current_state The current model state
  !! @param neighbour_description The halo swapping description of the neighbour we are accessing the buffer of
  !! @param dim The dimension we receive for
  !! @param target_index The target index for the dimension we are receiving for
  !! @param neighbour_location The location in the local neighbour data stores of this neighbour
  !! @param current_page The current, next, halo swap page to read from (all previous have been read and copied
  !!        already)
  !! @param source_data Optional source data which is read from into send buffers and written into by receieve 
  !!        buffers
  subroutine copy_halo_buffer_to_field(current_state, neighbour_description, dim, target_index, &
       neighbour_location, current_page, source_data)
    type(model_state_type), intent(inout) :: current_state
    integer, intent(in) :: dim, target_index, neighbour_location
    integer, intent(inout) :: current_page(:)
    type(neighbour_description_type), intent(inout) :: neighbour_description
    type(field_data_wrapper_type), dimension(:), intent(in), optional :: source_data

    integer :: page_bookmark, i

    page_bookmark = current_page(neighbour_location)
#ifdef U_ACTIVE
    call copy_buffer_to_field(current_state%local_grid, neighbour_description%recv_halo_buffer, &
         current_state%u%data,  dim, target_index, page_bookmark)
    page_bookmark=page_bookmark+1
    call copy_buffer_to_field(current_state%local_grid, neighbour_description%recv_halo_buffer, &
         current_state%zu%data, dim, target_index, page_bookmark)
    page_bookmark=page_bookmark+1
#endif
#ifdef V_ACTIVE
    call copy_buffer_to_field(current_state%local_grid, neighbour_description%recv_halo_buffer, &
         current_state%v%data, dim, target_index, page_bookmark)
    page_bookmark=page_bookmark+1
    call copy_buffer_to_field(current_state%local_grid, neighbour_description%recv_halo_buffer, &
         current_state%zv%data, dim, target_index, page_bookmark)
    page_bookmark=page_bookmark+1
#endif
#ifdef W_ACTIVE
    call copy_buffer_to_field(current_state%local_grid, neighbour_description%recv_halo_buffer, &
         current_state%w%data, dim, target_index, page_bookmark)
    page_bookmark=page_bookmark+1
    call copy_buffer_to_field(current_state%local_grid, neighbour_description%recv_halo_buffer, &
         current_state%zw%data, dim, target_index, page_bookmark)
    page_bookmark=page_bookmark+1
#endif
    if (current_state%th%active) then 
       call copy_buffer_to_field(current_state%local_grid, neighbour_description%recv_halo_buffer, &
            current_state%th%data, dim, target_index, page_bookmark)
       page_bookmark=page_bookmark+1
       call copy_buffer_to_field(current_state%local_grid, neighbour_description%recv_halo_buffer, &
            current_state%zth%data, dim, target_index, page_bookmark)
       page_bookmark=page_bookmark+1
    end if
    do i=1,current_state%number_q_fields
       if (current_state%q(i)%active) then
          call copy_buffer_to_field(current_state%local_grid, neighbour_description%recv_halo_buffer, &
               current_state%q(i)%data, dim, target_index, page_bookmark)
          page_bookmark=page_bookmark+1
          call copy_buffer_to_field(current_state%local_grid, neighbour_description%recv_halo_buffer, &
               current_state%zq(i)%data, dim, target_index, page_bookmark)
          page_bookmark=page_bookmark+1        
       end if
    end do
    current_page(neighbour_location)=page_bookmark
  end subroutine copy_halo_buffer_to_field

  !> Copies the prognostic field data to halo buffers for a specific process in a dimension and
  !! halo cell
  !! @param current_state The current model state
  !! @param neighbour_descriptions Description of the neighbour halo swapping status
  !! @param dim Dimension to copy from
  !! @param source_index The source index of the dimension we are reading from in the prognostic
  !!        field
  !! @param pid_location Location of the neighbouring process in the local stored data structures
  !! @param current_page The current (next) buffer page to copy into
  !! @param source_data Optional source data which is read from into send buffers and written 
  !!        into by receieve buffers
  subroutine copy_fields_to_halo_buffer(current_state, neighbour_description, dim, source_index,&
       pid_location, current_page, source_data)
    
    type(model_state_type), intent(inout) :: current_state
    integer, intent(in) :: dim, pid_location, source_index
    integer, intent(inout) :: current_page(:)
    type(neighbour_description_type), intent(inout) :: neighbour_description
    type(field_data_wrapper_type), dimension(:), intent(in), optional :: source_data

    integer :: page_bookmark, i

    page_bookmark = current_page(pid_location)

#ifdef U_ACTIVE
    call copy_field_to_buffer(current_state%local_grid, neighbour_description%send_halo_buffer,&
         current_state%u%data, dim, source_index, page_bookmark)
    page_bookmark = page_bookmark + 1
    call copy_field_to_buffer(current_state%local_grid, neighbour_description%send_halo_buffer, &
         current_state%zu%data, dim, source_index, page_bookmark)
    page_bookmark = page_bookmark + 1
#endif
#ifdef V_ACTIVE
    call copy_field_to_buffer(current_state%local_grid, neighbour_description%send_halo_buffer, &
         current_state%v%data, dim, source_index, page_bookmark)
    page_bookmark = page_bookmark + 1
    call copy_field_to_buffer(current_state%local_grid, neighbour_description%send_halo_buffer, &
         current_state%zv%data,  dim, source_index, page_bookmark)
    page_bookmark = page_bookmark + 1
#endif
#ifdef W_ACTIVE
    call copy_field_to_buffer(current_state%local_grid, neighbour_description%send_halo_buffer, &
         current_state%w%data, dim, source_index, page_bookmark)
    page_bookmark = page_bookmark+1
    call copy_field_to_buffer(current_state%local_grid, neighbour_description%send_halo_buffer, &
         current_state%zw%data,  dim, source_index, page_bookmark)
    page_bookmark = page_bookmark+1
#endif
    if (current_state%th%active) then 
       call copy_field_to_buffer(current_state%local_grid, &
            neighbour_description%send_halo_buffer, current_state%th%data,  dim, source_index,&
            page_bookmark)
      page_bookmark = page_bookmark+1
      call copy_field_to_buffer(current_state%local_grid, &
           neighbour_description%send_halo_buffer, current_state%zth%data, dim, source_index, &
           page_bookmark)
      page_bookmark = page_bookmark+1
   end if
   
    do i = 1,current_state%number_q_fields
       if (current_state%q(i)%active) then
          call copy_field_to_buffer(current_state%local_grid, &
               neighbour_description%send_halo_buffer, current_state%q(i)%data, dim, &
               source_index, page_bookmark)
          page_bookmark = page_bookmark + 1
        call copy_field_to_buffer(current_state%local_grid, &
             neighbour_description%send_halo_buffer, &
             current_state%zq(i)%data, dim, source_index, page_bookmark)
        page_bookmark = page_bookmark + 1        
      end if
    end do
    current_page(pid_location) = page_bookmark
  end subroutine copy_fields_to_halo_buffer

  !> Copies the prognostic corner field data to halo buffers for a specific process in a 
  !! dimension and halo cell
  !! @param current_state The current model state
  !! @param neighbour_descriptions Description of the neighbour halo swapping status
  !! @param corner_loc Location of the corner
  !! @param x_source_index The X source index of the dimension we are reading from in the 
  !! prognostic field
  !! @param y_source_index The Y source index of the dimension we are reading from in the 
  !! prognostic field
  !! @param pid_location Location of the neighbouring process in the local stored data structures
  !! @param current_page The current (next) buffer page to copy into
  !! @param source_data Optional source data which is read from into send buffers and written
  !! into by receieve buffers
  subroutine copy_corners_to_halo_buffer(current_state, neighbour_description, corner_loc, &
       x_source_index, y_source_index, pid_location, current_page, source_data)

    type(model_state_type), intent(inout) :: current_state
    integer, intent(in) :: corner_loc, pid_location, x_source_index, y_source_index
    integer, intent(inout) :: current_page(:)
    type(neighbour_description_type), intent(inout) :: neighbour_description
    type(field_data_wrapper_type), dimension(:), intent(in), optional :: source_data

    integer :: page_bookmark, i

    page_bookmark = current_page(pid_location)

#ifdef U_ACTIVE
    call copy_corner_to_buffer(current_state%local_grid, &
         neighbour_description%send_corner_buffer,&
         current_state%u%data, corner_loc, x_source_index, y_source_index, page_bookmark)
    page_bookmark=page_bookmark+1
    call copy_corner_to_buffer(current_state%local_grid, &
         neighbour_description%send_corner_buffer, &
         current_state%zu%data, corner_loc, x_source_index, y_source_index, page_bookmark)
    page_bookmark=page_bookmark+1
#endif
#ifdef V_ACTIVE
    call copy_corner_to_buffer(current_state%local_grid, &
         neighbour_description%send_corner_buffer, &
         current_state%v%data, corner_loc, x_source_index, y_source_index, page_bookmark)
    page_bookmark=page_bookmark+1
    call copy_corner_to_buffer(current_state%local_grid, &
         neighbour_description%send_corner_buffer, &
         current_state%zv%data, corner_loc, x_source_index, y_source_index, page_bookmark)
    page_bookmark=page_bookmark+1
#endif
#ifdef W_ACTIVE
    call copy_corner_to_buffer(current_state%local_grid, &
         neighbour_description%send_corner_buffer, &
         current_state%w%data, corner_loc, x_source_index, y_source_index, page_bookmark)
    page_bookmark=page_bookmark+1
    call copy_corner_to_buffer(current_state%local_grid, &
         neighbour_description%send_corner_buffer, &
         current_state%zw%data, corner_loc, x_source_index, y_source_index, page_bookmark)
    page_bookmark=page_bookmark+1
#endif
    if (current_state%th%active) then 
      call copy_corner_to_buffer(current_state%local_grid, &
           neighbour_description%send_corner_buffer,&
           current_state%th%data,  corner_loc, x_source_index, y_source_index, page_bookmark)
      page_bookmark=page_bookmark+1
      call copy_corner_to_buffer(current_state%local_grid, &
           neighbour_description%send_corner_buffer, &
           current_state%zth%data, corner_loc, x_source_index, y_source_index, page_bookmark)
      page_bookmark=page_bookmark+1
    end if
    do i=1,current_state%number_q_fields
      if (current_state%q(i)%active) then
        call copy_corner_to_buffer(current_state%local_grid, &
             neighbour_description%send_corner_buffer, &
             current_state%q(i)%data, corner_loc, x_source_index, y_source_index, page_bookmark)
        page_bookmark=page_bookmark+1
        call copy_corner_to_buffer(current_state%local_grid, &
             neighbour_description%send_corner_buffer, &
             current_state%zq(i)%data, corner_loc, x_source_index, y_source_index, page_bookmark)
        page_bookmark=page_bookmark+1        
      end if
    end do
    current_page(pid_location)=page_bookmark
  end subroutine copy_corners_to_halo_buffer

  !> Deduces the number of fields per halo cell. This depends upon what fields are active in the model
  !! @param current_state The current model state
  integer function get_fields_per_halo_cell(current_state)
    type(model_state_type), intent(inout) :: current_state

    integer :: i

    get_fields_per_halo_cell=0

#ifdef U_ACTIVE
    get_fields_per_halo_cell=get_fields_per_halo_cell+2
#endif
#ifdef V_ACTIVE
    get_fields_per_halo_cell=get_fields_per_halo_cell+2
#endif
#ifdef W_ACTIVE
    get_fields_per_halo_cell=get_fields_per_halo_cell+2
#endif
    if (current_state%th%active) then 
      get_fields_per_halo_cell=get_fields_per_halo_cell+2
    end if    
    do i=1,current_state%number_q_fields
      if (current_state%q(i)%active) then
        get_fields_per_halo_cell=get_fields_per_halo_cell+2
      end if
    end do
  end function get_fields_per_halo_cell    

  !> Displays some debugging information about who is sending what if that logging_mod level is selected
  !! @param neighbour_counts Number of neighbours that I have
  !! @param included_fields Number of prognostic fields that we include in this halo swap
  !! @param rank My current rank
  subroutine display_debugging_info_if_needed(neighbour_counts, included_fields, rank)
    integer, intent(in) :: neighbour_counts, rank, included_fields

    if (.not. first_call .or. log_get_logging_level() .lt. LOG_DEBUG) return
    call log_log(LOG_DEBUG, "Rank "//trim(conv_to_string(rank))//": "//trim(conv_to_string(neighbour_counts))//&
         " neighbours per timestep over "//trim(conv_to_string(included_fields))//" fields")
    first_call = .false.
  end subroutine display_debugging_info_if_needed

  !> Will do any local copying of data required for the boundary conditions. I.e. if all columns in a slice
  !! are on a process then it will copy in the y dimension
  !! @param current_state The current model state_mod
  !! @param copy_counts Number of local copies performed
  !! @param source_data Optional source data which is read from into send buffers and written into by receieve
  !!                     buffers
  subroutine perform_local_data_copy_for_all_prognostics(current_state, halo_depth, involve_corners, source_data)
    type(model_state_type), intent(inout) :: current_state
    integer, intent(in) :: halo_depth
    logical, intent(in) :: involve_corners
    type(field_data_wrapper_type), dimension(:), intent(in), optional :: source_data

    integer :: i

#ifdef U_ACTIVE
    call perform_local_data_copy_for_field(current_state%u%data, current_state%local_grid, &
         current_state%parallel%my_rank, halo_depth, involve_corners)
    call perform_local_data_copy_for_field(current_state%zu%data, current_state%local_grid, &
         current_state%parallel%my_rank, halo_depth, involve_corners)
#endif    
#ifdef V_ACTIVE
    call perform_local_data_copy_for_field(current_state%v%data, current_state%local_grid, &
         current_state%parallel%my_rank, halo_depth, involve_corners)
    call perform_local_data_copy_for_field(current_state%zv%data, current_state%local_grid, &
         current_state%parallel%my_rank, halo_depth, involve_corners)
#endif
#ifdef W_ACTIVE
    call perform_local_data_copy_for_field(current_state%w%data, current_state%local_grid, &
         current_state%parallel%my_rank, halo_depth, involve_corners)
    call perform_local_data_copy_for_field(current_state%zw%data, current_state%local_grid, &
         current_state%parallel%my_rank, halo_depth, involve_corners)
#endif
    if (current_state%th%active) then
      call perform_local_data_copy_for_field(current_state%th%data, current_state%local_grid, &
           current_state%parallel%my_rank, halo_depth, involve_corners)
      call perform_local_data_copy_for_field(current_state%zth%data, current_state%local_grid, &
           current_state%parallel%my_rank, halo_depth, involve_corners)
    end if
    do i=1,current_state%number_q_fields
      if (current_state%q(i)%active) then
        call perform_local_data_copy_for_field(current_state%q(i)%data, current_state%local_grid, &
             current_state%parallel%my_rank, halo_depth, involve_corners)
        call perform_local_data_copy_for_field(current_state%zq(i)%data, current_state%local_grid, &
             current_state%parallel%my_rank, halo_depth, involve_corners)
      end if
    end do
  end subroutine perform_local_data_copy_for_all_prognostics
end module haloswapper_mod
