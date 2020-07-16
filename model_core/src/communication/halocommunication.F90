!> Provides the mechanism for halo swapping. This module contains the functionality required to 
!! determine what messages get sent where, pack data up en mass, send and receive from 
!! neighbouring processes and unpack into the appropriate data locations. There is also some 
!! functionality to help local copying of data to corresponding locations. The idea is that the 
!! caller will determine the policy (i.e. exactly what fields are to be communicated) through 
!! procedure arguments and this mechanism can be used again and again. It implements bulk sending
!! of all field data in one large message to reduce communication overhead
module halo_communication_mod
  use datadefn_mod, only : DEFAULT_PRECISION, PRECISION_TYPE
  use communication_types_mod, only : halo_communication_type, neighbour_description_type, &
       field_data_wrapper_type
  use prognostics_mod, only : prognostic_field_type
  use state_mod, only : model_state_type
  use conversions_mod, only : conv_to_string
  use logging_mod, only : LOG_ERROR, LOG_WARN, log_log,LOG_DEBUG
  use grids_mod, only : local_grid_type, X_INDEX, Y_INDEX, Z_INDEX
  use mpi, only : MPI_REQUEST_NULL, MPI_STATUSES_IGNORE
  implicit none

#ifndef TEST_MODE
  private
#endif  

  !> Procedure interfaces used to determine the policy (i.e. the fields) of halo swapping and 
  !  data copying
  interface
     integer function get_fields_per_halo_cell_proc_interface(current_state)
       import model_state_type
       type(model_state_type), intent(inout) :: current_state
     end function get_fields_per_halo_cell_proc_interface

     subroutine copy_fields_to_halo_buffer_proc_interface(current_state, neighbour_description, &
          dim, source_index, &
          pid_location, current_page, source_data)
       import model_state_type, neighbour_description_type, field_data_wrapper_type
       type(model_state_type), intent(inout) :: current_state
       integer, intent(in) :: dim, pid_location, source_index
       integer, intent(inout) :: current_page(:)
       type(neighbour_description_type), intent(inout) :: neighbour_description
       type(field_data_wrapper_type), dimension(:), intent(in), optional :: source_data
     end subroutine copy_fields_to_halo_buffer_proc_interface

     subroutine copy_corners_to_halo_buffer_proc_interface(current_state, neighbour_description, &
          dim, x_source_index, &
          y_source_index, pid_location, current_page, source_data)
       import model_state_type, neighbour_description_type, field_data_wrapper_type
       type(model_state_type), intent(inout) :: current_state
       integer, intent(in) :: dim, pid_location, x_source_index, y_source_index
       integer, intent(inout) :: current_page(:)
       type(neighbour_description_type), intent(inout) :: neighbour_description
       type(field_data_wrapper_type), dimension(:), intent(in), optional :: source_data
     end subroutine copy_corners_to_halo_buffer_proc_interface

     subroutine copy_halo_buffer_to_field_proc_interface(current_state, neighbour_description, &
          dim, target_index, &
          neighbour_location, current_page, source_data)
       import model_state_type, neighbour_description_type, field_data_wrapper_type
       type(model_state_type), intent(inout) :: current_state
       integer, intent(in) :: dim, target_index, neighbour_location
       integer, intent(inout) :: current_page(:)
       type(neighbour_description_type), intent(inout) :: neighbour_description
       type(field_data_wrapper_type), dimension(:), intent(in), optional :: source_data
     end subroutine copy_halo_buffer_to_field_proc_interface

     subroutine copy_halo_buffer_to_corner_proc_interface(current_state, neighbour_description,&
          corner_loc, x_target_index, &
          y_target_index, neighbour_location, current_page, source_data)
       import model_state_type, neighbour_description_type, field_data_wrapper_type
       type(model_state_type), intent(inout) :: current_state
       integer, intent(in) :: corner_loc, x_target_index, y_target_index, neighbour_location
       integer, intent(inout) :: current_page(:)
       type(neighbour_description_type), intent(inout) :: neighbour_description
       type(field_data_wrapper_type), dimension(:), intent(in), optional :: source_data
     end subroutine copy_halo_buffer_to_corner_proc_interface

     subroutine perform_local_data_copy_proc_interface(current_state, halo_depth, &
          involve_corners, source_data)
       import model_state_type, field_data_wrapper_type
       type(model_state_type), intent(inout) :: current_state
       integer, intent(in) :: halo_depth
       logical, intent(in) :: involve_corners
       type(field_data_wrapper_type), dimension(:), intent(in), optional :: source_data
     end subroutine perform_local_data_copy_proc_interface     
  end interface
  
  public copy_buffer_to_field, copy_field_to_buffer, copy_corner_to_buffer, &
       copy_buffer_to_corner, perform_local_data_copy_for_field, init_halo_communication, &
       finalise_halo_communication, initiate_nonblocking_halo_swap, &
       complete_nonblocking_halo_swap, blocking_halo_swap, get_single_field_per_halo_cell       
contains

  !> Performs the entire halo swap operation, this is simply a wrapper around the nonblocking 
  !! initiate and complete procedures and saves the programmer from calling these directly if
  !! they do not wish to interleave any computation
  !! @param current_state The current model state
  !! @param halo_swap_state The halo swapping state
  !! @param copy_to_halo_buffer Procedure pointer which is called to copy the data into the
  !!                            halo send buffer
  !! @param perform_local_data_copy Procedure pointer which performs local data copying 
  !!                                (where the neighbour is the local process)
  !! @param copy_from_halo_buffer Procedure pointer which copies received data from the halo 
  !!                              buffer into the field
  !! @param copy_corners_to_halo_buffer Optional procedure pointer which copies data in corners
  !!                                    to the halo send buffer
  !! @param copy_from_halo_buffer_to_corner Optional procedure pointer which copies
  !!                                        received data into halo corners
  !! @param source_data Optional source data which is read from into send buffers and written
  !!                    into by receieve buffers
  subroutine blocking_halo_swap(current_state, halo_swap_state, copy_to_halo_buffer, &
       perform_local_data_copy, copy_from_halo_buffer, copy_corners_to_halo_buffer, &
       copy_from_halo_buffer_to_corner, source_data)
    
    type(model_state_type), intent(inout) :: current_state
    type(halo_communication_type), intent(inout) :: halo_swap_state
    procedure(copy_fields_to_halo_buffer_proc_interface) :: copy_to_halo_buffer
    procedure(copy_halo_buffer_to_field_proc_interface) :: copy_from_halo_buffer
    procedure(perform_local_data_copy_proc_interface) :: perform_local_data_copy
    procedure(copy_corners_to_halo_buffer_proc_interface), optional :: &
         copy_corners_to_halo_buffer
    procedure(copy_halo_buffer_to_corner_proc_interface), optional :: &
         copy_from_halo_buffer_to_corner
    type(field_data_wrapper_type), dimension(:), intent(in), optional :: source_data
    
    if ((present(copy_corners_to_halo_buffer) .and. .not. &
         present(copy_from_halo_buffer_to_corner)) .or. &
         (.not. present(copy_corners_to_halo_buffer) .and. &
         present(copy_from_halo_buffer_to_corner))) then
       call log_log(LOG_ERROR, "You must either provie no or both corner optional arguments to the halo swap function")
    end if
    
    if (present(source_data)) then
       if (present(copy_corners_to_halo_buffer)) then
          call initiate_nonblocking_halo_swap(current_state, halo_swap_state, &
               copy_to_halo_buffer, copy_corners_to_halo_buffer, source_data)
       else
          call initiate_nonblocking_halo_swap(current_state, halo_swap_state, &
               copy_to_halo_buffer, source_data=source_data)
       end if
      if (present(copy_from_halo_buffer_to_corner)) then
         call complete_nonblocking_halo_swap(current_state, halo_swap_state, &
              perform_local_data_copy, copy_from_halo_buffer, copy_from_halo_buffer_to_corner, &
              source_data)
      else
        call complete_nonblocking_halo_swap(current_state, halo_swap_state, &
             perform_local_data_copy, copy_from_halo_buffer, source_data=source_data)
      end if      
    else
      if (present(copy_corners_to_halo_buffer)) then
        call initiate_nonblocking_halo_swap(current_state, halo_swap_state, copy_to_halo_buffer,&
             copy_corners_to_halo_buffer)
      else
        call initiate_nonblocking_halo_swap(current_state, halo_swap_state, copy_to_halo_buffer)
      end if
      if (present(copy_from_halo_buffer_to_corner)) then
         call complete_nonblocking_halo_swap(current_state, halo_swap_state, &
              perform_local_data_copy, &
             copy_from_halo_buffer, copy_from_halo_buffer_to_corner)
      else
        call complete_nonblocking_halo_swap(current_state, halo_swap_state, &
             perform_local_data_copy, copy_from_halo_buffer)
      end if      
    end if    
  end subroutine blocking_halo_swap
  
  !> This completes a nonblocking halo swap and it is only during this call that the data fields
  !! themselves are modified. This will perform local data copying, wait for all outstanding 
  !! receives to complete, copy the received buffer data into the data and then wait for all 
  !! outstanding sends to complete
  !! @param current_state The current model state
  !! @param halo_swap_state The halo swapping state
  !! @param perform_local_data_copy Procedure pointer which performs local data copying (where 
  !!                                the neighbour is the local process)
  !! @param copy_from_halo_buffer Procedure pointer which copies received data from the halo 
  !!                              buffer into the field
  !! @param copy_from_halo_buffer_to_corner Optional procedure pointer which copies halo data 
  !!                                        into field corners
  !! @param source_data Optional source data which is read from into send buffers and written 
  !!                    into by receieve buffers
  subroutine complete_nonblocking_halo_swap(current_state, halo_swap_state, &
       perform_local_data_copy, copy_from_halo_buffer, copy_from_halo_buffer_to_corner, &
       source_data)
    type(model_state_type), intent(inout) :: current_state
    type(halo_communication_type), intent(inout) :: halo_swap_state
    procedure(copy_halo_buffer_to_field_proc_interface) :: copy_from_halo_buffer
    procedure(perform_local_data_copy_proc_interface) :: perform_local_data_copy
    procedure(copy_halo_buffer_to_corner_proc_interface), optional :: &
         copy_from_halo_buffer_to_corner
    type(field_data_wrapper_type), dimension(:), intent(in), optional :: source_data

    integer :: ierr

    if ((present(copy_from_halo_buffer_to_corner) .and. .not. halo_swap_state%involve_corners) &
         .or. (.not. present(copy_from_halo_buffer_to_corner) .and. &
         halo_swap_state%involve_corners)) then
      call log_log(LOG_WARN, "Inconsistent halo swap corner state and corner subroutine call arguments")
    end if    

    if (present(source_data)) then
      call perform_local_data_copy(current_state, halo_swap_state%halo_depth, &
           halo_swap_state%involve_corners, source_data)
    else
      call perform_local_data_copy(current_state, halo_swap_state%halo_depth, &
           halo_swap_state%involve_corners)
    end if
    if (halo_swap_state%number_distinct_neighbours .gt. 0) then
      call mpi_waitall(size(halo_swap_state%recv_requests), halo_swap_state%recv_requests, &
           MPI_STATUSES_IGNORE, ierr)
      if (present(source_data)) then
        if (halo_swap_state%involve_corners .and. present(copy_from_halo_buffer_to_corner)) then
          call copy_buffer_data_for_prognostics(current_state, halo_swap_state, &
               copy_from_halo_buffer, copy_from_halo_buffer_to_corner, source_data)
        else
          call copy_buffer_data_for_prognostics(current_state, halo_swap_state, &
               copy_from_halo_buffer, source_data=source_data)
        end if
      else
        if (halo_swap_state%involve_corners .and. present(copy_from_halo_buffer_to_corner)) then
          call copy_buffer_data_for_prognostics(current_state, halo_swap_state, &
               copy_from_halo_buffer, copy_from_halo_buffer_to_corner)
        else
          call copy_buffer_data_for_prognostics(current_state, halo_swap_state, &
               copy_from_halo_buffer)
        end if
      end if
      call mpi_waitall(size(halo_swap_state%send_requests), halo_swap_state%send_requests, &
           MPI_STATUSES_IGNORE, ierr)
    end if
    halo_swap_state%swap_in_progress=.false.
  end subroutine complete_nonblocking_halo_swap

  !> Initiates a non blocking halo swap, this registers the receive requests, copies data into
  !! the send buffer and then registers send requests for these. As far as this call is 
  !! concerned, the data is immutable - no modifications will take place to it until the
  !! matching completion call is made
  !! @param current_state The current model state
  !! @param halo_swap_state The halo swapping state
  !! @param copy_to_halo_buffer Procedure pointer which is called to copy the data into the halo 
  !!                            send buffer
  !! @param copy_corners_to_halo_buffer Optional procedure pointer which copies field corners 
  !!                                    into halo buffer
  !! @param source_data Optional source data which is read from into send buffers and written 
  !!                    into by receieve buffers
  subroutine initiate_nonblocking_halo_swap(current_state, halo_swap_state, copy_to_halo_buffer, &
       copy_corners_to_halo_buffer, source_data)
    
    type(model_state_type), intent(inout) :: current_state
    type(halo_communication_type), intent(inout) :: halo_swap_state
    procedure(copy_fields_to_halo_buffer_proc_interface) :: copy_to_halo_buffer
    procedure(copy_corners_to_halo_buffer_proc_interface), optional :: copy_corners_to_halo_buffer
    type(field_data_wrapper_type), dimension(:), intent(in), optional :: source_data
    
    halo_swap_state%swap_in_progress = .true.
    if (halo_swap_state%number_distinct_neighbours .gt. 0) then
       halo_swap_state%send_requests(:) = MPI_REQUEST_NULL
       halo_swap_state%recv_requests(:) = MPI_REQUEST_NULL
       
       if ((present(copy_corners_to_halo_buffer) .and. .not. halo_swap_state%involve_corners) &
            .or.  (.not. present(copy_corners_to_halo_buffer) .and. &
            halo_swap_state%involve_corners)) then
          call log_log(LOG_WARN, "Inconsistent halo swap corner state and corner subroutine call arguments")
       end if

       ! we call recv before send
       call recv_all_halos(current_state, halo_swap_state)
       
       if (present(source_data)) then
          if (halo_swap_state%involve_corners .and. present(copy_corners_to_halo_buffer)) then
             call send_all_halos(current_state, halo_swap_state, copy_to_halo_buffer, &
                  copy_corners_to_halo_buffer, source_data)
          else
             call send_all_halos(current_state, halo_swap_state, copy_to_halo_buffer, &
                  source_data=source_data)
          end if
       else
          if (halo_swap_state%involve_corners .and. present(copy_corners_to_halo_buffer)) then
             call send_all_halos(current_state, halo_swap_state, copy_to_halo_buffer, &
                  copy_corners_to_halo_buffer)
          else
             call send_all_halos(current_state, halo_swap_state, copy_to_halo_buffer)
          end if
       end if
    end if
  end subroutine initiate_nonblocking_halo_swap  

  !> Initialises a halo swapping state, by determining the neighbours, size of data in each swap
  !! and allocating the required memory - which are communication buffers and request handles.
  !! All this information is placed into the returned state
  !! @param current_state The current model state
  !! @param get_fields_per_halo_cell Procedure pointer to the procedure which determines how many
  !!                                 fields per halo cell there are
  !! @param involve_corners Whether this involves corners or not
  !! @returns A halo swapping state which can be used for halo swapping operations
  subroutine init_halo_communication(current_state, get_fields_per_halo_cell, halo_state, &
       halo_depth, involve_corners)
    type(model_state_type), intent(inout) :: current_state
    procedure(get_fields_per_halo_cell_proc_interface) :: get_fields_per_halo_cell
    logical, intent(in) :: involve_corners
    integer, intent(in) :: halo_depth
    type(halo_communication_type), intent(out) :: halo_state

    integer :: number_comm_requests

    halo_state%involve_corners = involve_corners
    halo_state%halo_depth = halo_depth
    halo_state%number_distinct_neighbours = get_number_of_processes_involved_in_communication(&
         current_state%local_grid, current_state%parallel%my_rank, involve_corners)
    if (halo_state%number_distinct_neighbours .gt. 0) then
      allocate(halo_state%halo_swap_neighbours(halo_state%number_distinct_neighbours))
      halo_state%halo_swap_neighbours = populate_halo_swap_neighbours(current_state%local_grid, &
           current_state%parallel%my_rank, halo_state%number_distinct_neighbours, involve_corners)
      call deduce_halo_pages_per_neighbour(current_state, halo_state%halo_swap_neighbours, &
           halo_state%number_distinct_neighbours, get_fields_per_halo_cell, &
           halo_state%fields_per_cell, halo_depth)
      if (involve_corners) call deduce_halo_corners_per_neighbour(current_state, &
           halo_state%halo_swap_neighbours, &
           halo_state%number_distinct_neighbours, halo_state%fields_per_cell)
      call allocate_halo_buffers_for_each_neighbour(current_state%local_grid, &
           halo_state%number_distinct_neighbours, &
           halo_state%halo_swap_neighbours)
      call determine_recv_and_send_sizes(current_state%local_grid, &
           halo_state%halo_swap_neighbours, &
           halo_state%number_distinct_neighbours, involve_corners)
      call generate_recv_field_buffer_matches(current_state, halo_state%halo_depth, &
           halo_state%cell_match)
      
      ! required for nonblocking MPI communcations
      number_comm_requests = get_number_communication_requests(halo_state%halo_swap_neighbours, &
           halo_state%number_distinct_neighbours)
      allocate(halo_state%send_requests(number_comm_requests), &
           halo_state%recv_requests(number_comm_requests))
      halo_state%send_requests(:) = MPI_REQUEST_NULL
      halo_state%recv_requests(:) = MPI_REQUEST_NULL
    end if
    halo_state%initialised = .true.
  end subroutine init_halo_communication

  !> Finalises the halo swap represented by the state by freeing up all the allocated memory
  !! @param halo_swap_state State of the specific halo swap
  subroutine finalise_halo_communication(halo_swap_state)
    type(halo_communication_type), intent(inout) :: halo_swap_state

    integer :: i

    ! TODO - issue cancel of outstanding requests here
    if (allocated(halo_swap_state%send_requests)) deallocate(halo_swap_state%send_requests)
    if (allocated(halo_swap_state%recv_requests)) deallocate(halo_swap_state%recv_requests)

    do i=1,halo_swap_state%number_distinct_neighbours
      if (allocated(halo_swap_state%halo_swap_neighbours(i)%send_halo_buffer)) &
           deallocate(halo_swap_state%halo_swap_neighbours(i)%send_halo_buffer)
      if (allocated(halo_swap_state%halo_swap_neighbours(i)%recv_halo_buffer)) &
           deallocate(halo_swap_state%halo_swap_neighbours(i)%recv_halo_buffer)
    end do
    if (allocated(halo_swap_state%halo_swap_neighbours)) deallocate(halo_swap_state%halo_swap_neighbours)
    halo_swap_state%initialised=.false.
  end subroutine finalise_halo_communication

  !> Copies the received buffer for a specific field to the corresponding corner of that field
  !! @param local_grid Description of the local grid
  !! @param halo_buffer Raw halo buffer data
  !! @param field_data Raw prognostic field data
  !! @param corner_loc Location of the corner
  !! @param x_target_index The X target index that we copy into
  !! @param y_target_index The Y target index that we copy into
  !! @param halo_page The halo page to read from
  subroutine copy_buffer_to_corner(local_grid, halo_buffer, field_data, corner_loc, &
       x_target_index, y_target_index, halo_page)
    type(local_grid_type), intent(inout) :: local_grid
    integer, intent(in) :: corner_loc, x_target_index, y_target_index, halo_page
    real(kind=DEFAULT_PRECISION), dimension(:,:,:), intent(inout) :: halo_buffer
    real(kind=DEFAULT_PRECISION), dimension(:,:,:), intent(inout) :: field_data

    field_data(local_grid%local_domain_start_index(Z_INDEX):&
         local_grid%local_domain_end_index(Z_INDEX),&
         y_target_index, x_target_index)=halo_buffer(:,4,halo_page)
    field_data(local_grid%local_domain_start_index(Z_INDEX):&
         local_grid%local_domain_end_index(Z_INDEX),&
         merge(y_target_index-1, y_target_index+1, corner_loc .lt. 3), x_target_index)=&
         halo_buffer(:,3,halo_page)
    field_data(local_grid%local_domain_start_index(Z_INDEX):&
         local_grid%local_domain_end_index(Z_INDEX),&
         y_target_index, merge(x_target_index-1, x_target_index+1, corner_loc == 1 .or.&
         corner_loc == 3))= halo_buffer(:,2,halo_page)
    field_data(local_grid%local_domain_start_index(Z_INDEX):&
         local_grid%local_domain_end_index(Z_INDEX),&
         merge(y_target_index-1, y_target_index+1, corner_loc .lt. 3), &
         merge(x_target_index-1, x_target_index+1, corner_loc == 1 .or.&
         corner_loc == 3))= halo_buffer(:,1,halo_page)    
  end subroutine copy_buffer_to_corner

  !> Copies the received buffer for a specific field to the corresponding halo data of that
  !! prognostic field
  !! @param local_grid Description of the local grid
  !! @param halo_buffer Raw halo buffer data
  !! @param field_data Raw prognostic field data
  !! @param dim The dimension to copy into
  !! @param target_index The target index in the dimension that we copy into
  !! @param halo_page The halo page to read from
  subroutine copy_buffer_to_field(local_grid, halo_buffer, field_data, dim, target_index, &
       halo_page)
    type(local_grid_type), intent(inout) :: local_grid
    integer, intent(in) :: dim, target_index, halo_page
    real(kind=DEFAULT_PRECISION), dimension(:,:,:), intent(inout) :: halo_buffer
    real(kind=DEFAULT_PRECISION), dimension(:,:,:), intent(inout) :: field_data

    ! If the neighbours are the same then reverse our placement of the data due to wrapping
    ! around and order of
    ! messages being sent. This is not an issue if the neighbours are different    
    if (dim == X_INDEX) then
      field_data(local_grid%local_domain_start_index(Z_INDEX):&
           local_grid%local_domain_end_index(Z_INDEX),&
           local_grid%local_domain_start_index(Y_INDEX):&
           local_grid%local_domain_end_index(Y_INDEX), &
           target_index) = halo_buffer(:,:,halo_page)
    else
      field_data(local_grid%local_domain_start_index(Z_INDEX):&
           local_grid%local_domain_end_index(Z_INDEX), target_index, &
           local_grid%local_domain_start_index(X_INDEX):&
           local_grid%local_domain_end_index(X_INDEX)) = &
           halo_buffer(:,:,halo_page)
    end if
  end subroutine copy_buffer_to_field  

  !> Copies prognostic field data to send buffer for specific field, dimension, halo cell
  !! @param local_grid Description of the local grid
  !! @param halo_buffer Raw halo_buffer data that is written to
  !! @param field_data Raw prognostic field data that is read from
  !! @param dim Dimension that we are reading from
  !! @param source_index The index in the read dimension
  !! @param halo_page The halo buffer page to write to
  subroutine copy_field_to_buffer(local_grid, halo_buffer, field_data, dim, source_index, &
       halo_page)
    type(local_grid_type), intent(inout) :: local_grid
    integer, intent(in) :: dim, source_index, halo_page
    real(kind=DEFAULT_PRECISION), dimension(:,:,:), intent(inout) :: halo_buffer
    real(kind=DEFAULT_PRECISION), dimension(:,:,:), intent(inout) :: field_data
    
    if (dim == X_INDEX) then
       halo_buffer(:,:,halo_page) = field_data(local_grid%local_domain_start_index(Z_INDEX):&
            local_grid%local_domain_end_index(Z_INDEX), &
            local_grid%local_domain_start_index(Y_INDEX):&
            local_grid%local_domain_end_index(Y_INDEX), source_index)
    else
       halo_buffer(:,:,halo_page)=field_data(local_grid%local_domain_start_index(Z_INDEX):&
            local_grid%local_domain_end_index(Z_INDEX), source_index, &
            local_grid%local_domain_start_index(X_INDEX):&
            local_grid%local_domain_end_index(X_INDEX))
    end if
  end subroutine copy_field_to_buffer

  !> Copies prognostic field corner data to send buffer for specific field
  !! @param local_grid Description of the local grid
  !! @param halo_buffer Raw halo_buffer data that is written to
  !! @param field_data Raw prognostic field data that is read from
  !! @param corner_loc Location of the corner
  !! @param x_source_index The X index in the read dimension
  !! @param y_source_index The Y index in the read dimension
  !! @param halo_page The halo buffer page to write to
  subroutine copy_corner_to_buffer(local_grid, halo_buffer, field_data, corner_loc, &
       x_source_index, y_source_index, halo_page)
    type(local_grid_type), intent(inout) :: local_grid
    integer, intent(in) :: corner_loc, x_source_index, y_source_index, halo_page
    real(kind=DEFAULT_PRECISION), dimension(:,:,:), intent(inout) :: halo_buffer
    real(kind=DEFAULT_PRECISION), dimension(:,:,:), intent(inout) :: field_data

    !! TODO: hardcoded size of the corners
    halo_buffer(:,1,halo_page) = field_data(:,y_source_index,x_source_index)
    halo_buffer(:,2,halo_page) = field_data(:, merge(y_source_index-1, y_source_index+1,   &
         corner_loc .lt. 3), x_source_index)
    halo_buffer(:,3,halo_page) = field_data(:, y_source_index,                             &
         merge(x_source_index-1,x_source_index+1, corner_loc == 1 .or. corner_loc == 3))
    halo_buffer(:,4,halo_page) = field_data(:, merge(y_source_index-1, y_source_index+1,   &
         corner_loc .lt. 3), merge(x_source_index-1,x_source_index+1, corner_loc == 1 .or. &
         corner_loc == 3))
  end subroutine copy_corner_to_buffer

  !> Will perform a a local copy for the halo data of a field
  !! @param field_data The field data to perform the copy for
  !! @param local_grid Local grid information
  !! @param my_rank My process id
  subroutine perform_local_data_copy_for_field(field_data, local_grid, my_rank, halo_depth, &
       involve_corners)
    real(kind=DEFAULT_PRECISION), dimension(:,:,:), intent(inout) :: field_data
    type(local_grid_type), intent(inout) :: local_grid
    integer, intent(in) :: my_rank, halo_depth
    logical, intent(in) :: involve_corners

    call perform_local_data_copy_for_dimension(Y_INDEX, my_rank, halo_depth, local_grid, &
         field_data)
    call perform_local_data_copy_for_dimension(X_INDEX, my_rank, halo_depth, local_grid, &
         field_data)
    if (involve_corners) call perform_local_data_copy_for_corners(my_rank, local_grid, field_data)
  end subroutine perform_local_data_copy_for_field

  !--------------------------------------------------------------------------
  ! Private procedures acting as helpers
  !--------------------------------------------------------------------------

  !> Determines the overall number of communication requests, which is made up of normal halo 
  !! swaps and potentially corner
  !! swaps too if that is enabled
  !! @param halo_swap_neighbours Neighbouring halo swap state
  !! @param number_distinct_neighbours The number of distinct neighours that I will swap with
  !! @returns The number of communication requests
  integer function get_number_communication_requests(halo_swap_neighbours, &
       number_distinct_neighbours)
    integer, intent(in) :: number_distinct_neighbours
    type(neighbour_description_type), dimension(:), allocatable :: halo_swap_neighbours

    integer :: i

    get_number_communication_requests=0
    do i=1, number_distinct_neighbours
      if (halo_swap_neighbours(i)%recv_size .gt. 0) &
           get_number_communication_requests=get_number_communication_requests+1
      if (halo_swap_neighbours(i)%recv_corner_size .gt. 0)&
           get_number_communication_requests=get_number_communication_requests+1
    end do
  end function get_number_communication_requests  

  !> Determines the amount (in elements) of data that each neighbour will be sent and I will 
  !! receive from in a halo swap
  !! @param local_grid The local grid
  !! @param halo_swap_neighbours Structure describing state of halo swap neighbours
  !! @param number_distinct_neighbours The number of distinct neighbours I swap with
  !! @param involve_corners Whether or not to involve corners in a halo swap
  subroutine determine_recv_and_send_sizes(local_grid, halo_swap_neighbours, &
       number_distinct_neighbours, involve_corners)
    type(local_grid_type), intent(inout) :: local_grid
    integer, intent(in) :: number_distinct_neighbours
    logical, intent(in) :: involve_corners
    type(neighbour_description_type), dimension(:), allocatable :: halo_swap_neighbours
    
    integer :: i, normal_size, corner_size

    do i=1, number_distinct_neighbours
      if (halo_swap_neighbours(i)%halo_pages .gt. 0) then
        if (halo_swap_neighbours(i)%dimension == 0) then
          call log_log(LOG_ERROR, "Halo swapping with neighbour needed but dimension is 0 which suggests corner only")
        end if        
        normal_size=local_grid%size(Z_INDEX) * merge(local_grid%size(Y_INDEX), local_grid%size(X_INDEX), &
             halo_swap_neighbours(i)%dimension==X_INDEX)*halo_swap_neighbours(i)%halo_pages
      else
        normal_size=0
      end if
      halo_swap_neighbours(i)%recv_size=normal_size
      halo_swap_neighbours(i)%send_size=normal_size
      if (involve_corners .and. halo_swap_neighbours(i)%halo_corners .gt. 0) then
        ! For the moment assume both halos are the same neighbour - hence the 4, otherwise should call determine_halo_corner_element_sizes
        corner_size=local_grid%size(Z_INDEX)*4*halo_swap_neighbours(i)%halo_corners
      else
        corner_size=0
      end if
        halo_swap_neighbours(i)%recv_corner_size=corner_size
        halo_swap_neighbours(i)%send_corner_size=corner_size
    end do
  end subroutine determine_recv_and_send_sizes

  !> Determine the halo corner size in elements
  !! @param local_grid The local grid
  !! @param pid The process id that this calculation is based upon
  !! @reterns The number of elements involved in a specific corner halo swap
  integer function determine_halo_corner_size(local_grid)
    type(local_grid_type), intent(inout) :: local_grid

    determine_halo_corner_size = local_grid%halo_size(X_INDEX)*local_grid%halo_size(Y_INDEX)*&
         local_grid%size(Z_INDEX)    
  end function determine_halo_corner_size
 
  !> For a specific process id this determines the number of halo swap corner elements to involve
  !! in a communication
  !! @param local_grid The local grid
  !! @param pid The process id that this calculation is based upon
  !! @reterns The number of elements involved in a specific corner halo swap
  integer function determine_halo_corner_element_sizes(local_grid, pid)
    type(local_grid_type), intent(inout) :: local_grid
    integer, intent(in) :: pid

    integer :: i,j
    determine_halo_corner_element_sizes=0

    do i=1,size(local_grid%corner_neighbours, 2)
      do j=1,size(local_grid%corner_neighbours, 1)
        if (local_grid%corner_neighbours(j,i) == pid) then
          determine_halo_corner_element_sizes=&
               determine_halo_corner_element_sizes+local_grid%size(Z_INDEX)
          ! If second halo then there are 3 corner elements, therefore add extra two to this
          if (i==2) determine_halo_corner_element_sizes=determine_halo_corner_element_sizes+&
               local_grid%size(Z_INDEX)*2
        end if
      end do      
    end do    
  end function determine_halo_corner_element_sizes  

  !> Deduces the number of distinct neighbours that will be involved in a halo swap. This
  !! information is used to then allocate the appropriate amount of memory to store the
  !! neighbour halo swapping data structure
  !! @param local_grid Description of the local grid
  !! @param my_rank My global PID
  !! @param include_corners Whether to include corners or not
  integer function get_number_of_processes_involved_in_communication(local_grid, my_rank, &
       include_corners)
    type(local_grid_type), intent(inout) :: local_grid
    integer, intent(in) :: my_rank
    logical, intent(in) :: include_corners

    integer :: i, j, temp_neighbour_pids(merge(16, 8, include_corners))

    temp_neighbour_pids(:)=-1
    get_number_of_processes_involved_in_communication=0
    do i=2,3
      do j=1,4
        if (local_grid%neighbours(i,j) .ne. my_rank .and. .not. &
             has_pid_already_been_seen(temp_neighbour_pids, &
             local_grid%neighbours(i,j))) then
          get_number_of_processes_involved_in_communication = &
               get_number_of_processes_involved_in_communication+1
          temp_neighbour_pids(get_number_of_processes_involved_in_communication) = &
               local_grid%neighbours(i,j)
        end if        
      end do      
    end do

    if (include_corners) then
      do i=1,size(local_grid%corner_neighbours, 2)
        do j=1,size(local_grid%corner_neighbours, 1)
          if (local_grid%corner_neighbours(j,i) .ne. my_rank .and. .not. &
               has_pid_already_been_seen(temp_neighbour_pids, &
               local_grid%corner_neighbours(j,i))) then
            get_number_of_processes_involved_in_communication = &
                 get_number_of_processes_involved_in_communication+1
            temp_neighbour_pids(get_number_of_processes_involved_in_communication) = &
                 local_grid%corner_neighbours(j,i)
          end if          
        end do
      end do      
    end if
  end function get_number_of_processes_involved_in_communication

  !> Will populate the halo swap neighbour data strutures with appropriate neighbour
  !! pid and dimension numbers
  !! @param local_grid Description of the local grid
  !! @param my_rank My global PID
  !! @param number_distinct_neighbours The number of distinct neighbours that I have
  !! @param include_corners Whether to include corners or not
  function populate_halo_swap_neighbours(local_grid, my_rank, number_distinct_neighbours, &
       involve_corners)
    type(local_grid_type), intent(inout) :: local_grid
    integer, intent(in) :: my_rank, number_distinct_neighbours
    logical, intent(in) :: involve_corners

    type(neighbour_description_type), dimension(number_distinct_neighbours) :: &
         populate_halo_swap_neighbours
    integer :: i, j, current_pid_location, temp_neighbour_pids(merge(16, 8, involve_corners))

    current_pid_location=0
    temp_neighbour_pids(:)=-1
    do i=2,3
      do j=1,4
        if (local_grid%neighbours(i,j) .ne. my_rank .and. .not. &
             has_pid_already_been_seen(temp_neighbour_pids, &
             local_grid%neighbours(i,j))) then
          current_pid_location=current_pid_location+1
          populate_halo_swap_neighbours(current_pid_location)%pid=local_grid%neighbours(i,j)
          temp_neighbour_pids(current_pid_location)=local_grid%neighbours(i,j)
          populate_halo_swap_neighbours(current_pid_location)%dimension=i
        end if
      end do
    end do

    if (involve_corners) then
      do i=1,size(local_grid%corner_neighbours, 2)
        do j=1,size(local_grid%corner_neighbours, 1)
          if (local_grid%corner_neighbours(j,i) .ne. my_rank .and. .not. &
               has_pid_already_been_seen(temp_neighbour_pids, &
               local_grid%corner_neighbours(j,i))) then
            current_pid_location=current_pid_location+1
            populate_halo_swap_neighbours(current_pid_location)%pid = &
                 local_grid%corner_neighbours(j,i)
            temp_neighbour_pids(current_pid_location)=local_grid%corner_neighbours(j,i)
            populate_halo_swap_neighbours(current_pid_location)%dimension=0
          end if
        end do
      end do
    end if
  end function populate_halo_swap_neighbours

  !> Deduces the number of halo pages per neighbour halo swap and places this information in the appropriate data
  !! structures. We call a "page" of data the contiguous data of a field that we are going to send, such as 
  !! halo 1 of w, halo 1 of zw and halo 2 of w (assuming these go to the same neighbour as 1)
  !! @param current_state The current model state
  !! @param halo_swap_neighbours My neighbouring PIDs
  !! @param number_distinct_neighbours The number of distinct neighbours that I have
  !! @param get_fields_per_halo_cell Procedure pointer to get the number of fields per halo cell
  !! @param fields_per_cell The number of fields per cell is written into here
  subroutine deduce_halo_pages_per_neighbour(current_state, halo_swap_neighbours, &
       number_distinct_neighbours, get_fields_per_halo_cell, fields_per_cell, halo_depth)
    type(model_state_type), intent(inout) :: current_state
    type(neighbour_description_type), dimension(:), allocatable :: halo_swap_neighbours
    integer, intent(in) :: number_distinct_neighbours, halo_depth
    procedure(get_fields_per_halo_cell_proc_interface) :: get_fields_per_halo_cell
    integer, intent(out) :: fields_per_cell

    integer :: i, j, pid_location, halo_start, halo_end

    fields_per_cell = get_fields_per_halo_cell(current_state)
    halo_start = merge(2, 1, halo_depth==1)
    halo_end = merge(3, 4, halo_depth==1)

    ! i moves in x and y. z is 1
    do i = 2, 3
      do j = halo_start, halo_end
        if (current_state%parallel%my_rank .ne. current_state%local_grid%neighbours(i,j)) then
          pid_location = get_pid_neighbour_location(halo_swap_neighbours, &
               current_state%local_grid%neighbours(i,j), number_distinct_neighbours)
          halo_swap_neighbours(pid_location)%halo_pages = &
               halo_swap_neighbours(pid_location)%halo_pages + fields_per_cell
        end if
      end do
    end do
  end subroutine deduce_halo_pages_per_neighbour

  !> Determines the number of halo corners to swap between specific neighours, this is similar
  !! to deducing the number of halo
  !! pages per neighbour, only it is for corners
  !! @param current_state The current model state
  !! @param halo_swap_neighbours Description of neighbouring halo swap state
  !! @param number_distinct_neighbours The number of distinct neighbours that I swap with
  !! @param fields_per_cell The number of fields per cell is written into here
  subroutine deduce_halo_corners_per_neighbour(current_state, halo_swap_neighbours, &
       number_distinct_neighbours, fields_per_cell)
    type(model_state_type), intent(inout) :: current_state
    type(neighbour_description_type), dimension(:), allocatable :: halo_swap_neighbours
    integer, intent(in) :: number_distinct_neighbours, fields_per_cell

    integer :: i, j, pid_location
    ! i moves in x ,j in y
    do i=1,size(current_state%local_grid%corner_neighbours, 2)
      do j=1,size(current_state%local_grid%corner_neighbours, 1)
        if (current_state%parallel%my_rank .ne. &
             current_state%local_grid%corner_neighbours(j,i)) then
           pid_location = get_pid_neighbour_location(halo_swap_neighbours, &
                current_state%local_grid%corner_neighbours(j,i), number_distinct_neighbours)
           halo_swap_neighbours(pid_location)%halo_corners = &
                halo_swap_neighbours(pid_location)%halo_corners + fields_per_cell
        end if        
      end do      
    end do    
  end subroutine deduce_halo_corners_per_neighbour  

  !> Allocates the locally stored halo buffers (send and receive) for each neighbouring process
  !! @param local_grid Description of the local grid
  !! @param halo_swap_neighbours My neighbouring PIDs
  !! @param number_distinct_neighbours The number of distinct neighbours that I have
  subroutine allocate_halo_buffers_for_each_neighbour(local_grid, number_distinct_neighbours, &
       halo_swap_neighbours)
    type(local_grid_type), intent(inout) :: local_grid
    integer, intent(in) :: number_distinct_neighbours
    type(neighbour_description_type), dimension(:), allocatable, intent(inout) :: &
         halo_swap_neighbours

    integer :: i
    
    !! AH - test code to see if this corrects the allocation error on restart with gcc 7
     do i=1,number_distinct_neighbours
       if (allocated(halo_swap_neighbours(i)%send_halo_buffer)) &
            deallocate(halo_swap_neighbours(i)%send_halo_buffer)
       if (allocated(halo_swap_neighbours(i)%recv_halo_buffer)) &
            deallocate(halo_swap_neighbours(i)%recv_halo_buffer)
       if (allocated(halo_swap_neighbours(i)%send_corner_buffer)) &
           deallocate(halo_swap_neighbours(i)%send_corner_buffer) 
       if (allocated(halo_swap_neighbours(i)%recv_corner_buffer)) &
           deallocate(halo_swap_neighbours(i)%recv_corner_buffer)   
     end do
    !!

    do i=1,number_distinct_neighbours
      if (halo_swap_neighbours(i)%halo_pages .gt. 0) then
         !depending on the direction of the swapping, the send and recv buffer size would change
        allocate(halo_swap_neighbours(i)%send_halo_buffer(local_grid%size(Z_INDEX), &
             merge(local_grid%size(Y_INDEX), local_grid%size(X_INDEX), &
             halo_swap_neighbours(i)%dimension==X_INDEX), halo_swap_neighbours(i)%halo_pages))
        allocate(halo_swap_neighbours(i)%recv_halo_buffer(local_grid%size(Z_INDEX), &
             merge(local_grid%size(Y_INDEX), local_grid%size(X_INDEX), &
             halo_swap_neighbours(i)%dimension==X_INDEX), halo_swap_neighbours(i)%halo_pages))
      end if
      if (halo_swap_neighbours(i)%halo_corners .gt. 0) then
         ! is 4 because of the 4 cells to swap since the halo_depth is 2(2 in x and 2 in y)??
        allocate(halo_swap_neighbours(i)%send_corner_buffer(local_grid%size(Z_INDEX), 4, &
             halo_swap_neighbours(i)%halo_corners))
        allocate(halo_swap_neighbours(i)%recv_corner_buffer(local_grid%size(Z_INDEX), 4, &
             halo_swap_neighbours(i)%halo_corners))
      end if
    end do    
  end subroutine allocate_halo_buffers_for_each_neighbour

  !> Precalculates the received buffer to field halo cell matches for each dimension and called 
  !! from the initialisation stage
  !! @param current_state The current model state
  !! @param halo_depth The halo depth
  !! @param cell_match The matching cells are written into here
  subroutine generate_recv_field_buffer_matches(current_state, halo_depth, cell_match)
    type(model_state_type), intent(inout) :: current_state
    integer, intent(in) :: halo_depth
    integer, intent(out) :: cell_match(:,:)

    logical, dimension(3) :: same_neighbours
    integer :: i,j

    same_neighbours = retrieve_same_neighbour_information(current_state%local_grid)
    
    do i = 2,3
       if (halo_depth == 1) then
          cell_match(i,1)=1
          cell_match(i,2)=merge(2, 3, .not. same_neighbours(i))
          cell_match(i,3)=merge(3, 2, .not. same_neighbours(i))
          cell_match(i,4)=4
       else
          do j = 1, halo_depth
             cell_match(i,j) = merge(j, j+halo_depth, .not. same_neighbours(i))
             cell_match(i,j+halo_depth) = merge(j+halo_depth, j, .not. same_neighbours(i))
          end do
       !          cell_match(i,j) = merge(2, 4, .not. same_neighbours(i))
       !          cell_match(i,j) = merge(3, 1, .not. same_neighbours(i))
!          cell_match(i,j) = merge(4, 2, .not. same_neighbours(i))
       end if
    end do
  end subroutine generate_recv_field_buffer_matches

  !> Registers receive requests for all prognostic fields from the appropriate neighbouring 
  !! processes (that we have already deduced in the initialisation stage.)
  !! @param current_state The current model state
  !! @param halo_swap_state The halo swap state
  subroutine recv_all_halos(current_state, halo_swap_state)
    type(model_state_type), intent(inout) :: current_state
    type(halo_communication_type), intent(inout) :: halo_swap_state
    
    integer :: i, request_counter, ierr

    request_counter = 1 

    do i = 1, halo_swap_state%number_distinct_neighbours
       if (halo_swap_state%halo_swap_neighbours(i)%recv_size .gt. 0) then
          call mpi_irecv(halo_swap_state%halo_swap_neighbours(i)%recv_halo_buffer, &
               halo_swap_state%halo_swap_neighbours(i)%recv_size, PRECISION_TYPE, &
               halo_swap_state%halo_swap_neighbours(i)%pid, 0, &
               current_state%parallel%neighbour_comm, &
               halo_swap_state%recv_requests(request_counter), ierr)
          request_counter = request_counter + 1
       end if
       if (halo_swap_state%halo_swap_neighbours(i)%recv_corner_size .gt. 0) then
          call mpi_irecv(halo_swap_state%halo_swap_neighbours(i)%recv_corner_buffer, &
               halo_swap_state%halo_swap_neighbours(i)%recv_corner_size, PRECISION_TYPE, &
               halo_swap_state%halo_swap_neighbours(i)%pid, 0, &
               current_state%parallel%neighbour_comm, &
               halo_swap_state%recv_requests(request_counter), ierr)
          request_counter = request_counter + 1
       end if
    end do
  end subroutine recv_all_halos

  !> Copies all applicable bits of the prognostics into a send buffer for each neighbour and 
  !! then issues asynchronous sends of this data
  !! @param current_state The current model state
  !! @param halo_swap_state The halo swap state
  !! @param copy_fields_to_halo_buffer Procedure pointer to copy the data from the fields into 
  !!                                   the send buffer
  !! @param copy_corner_fields_to_halo_buffer Optional procedure pointer to copy corner fields 
  !!                                          into halo send buffer
  !! @param source_data Optional field data that can be copied from (passed to the user procedure)
  subroutine send_all_halos(current_state, halo_swap_state, copy_fields_to_halo_buffer, &
       copy_corner_fields_to_halo_buffer, source_data)

    type(model_state_type), intent(inout) :: current_state
    type(halo_communication_type), intent(inout) :: halo_swap_state
    procedure(copy_fields_to_halo_buffer_proc_interface) :: copy_fields_to_halo_buffer
    procedure(copy_corners_to_halo_buffer_proc_interface), optional :: &
         copy_corner_fields_to_halo_buffer
    type(field_data_wrapper_type), dimension(:), intent(in), optional :: source_data

    integer :: i, j, ierr,hstart, hend, pid_location, source_index, request_number, &
         x_source_index, y_source_index, current_page(halo_swap_state%number_distinct_neighbours),&
         halo_depth
    ! halo_size(Y_INDEX) = halo_size(X_INDEX), just pick one
    halo_depth = current_state%local_grid%halo_size(Y_INDEX)
    current_page(:) = 1
    request_number = 1
    ! TODO: hardcoded to halodepth 1 or 2
    hstart = merge(2, 1, halo_swap_state%halo_depth==1)
    hend = merge(3, halo_depth*2, halo_swap_state%halo_depth==1)
    
    do i=2, 3
       do j=hstart, hend
          if (current_state%parallel%my_rank .ne. current_state%local_grid%neighbours(i,j)) then
             
             if (j==1) then
                source_index = current_state%local_grid%local_domain_start_index(i)
             else if (j==2) then
                source_index = current_state%local_grid%local_domain_start_index(i) + &
                     merge(1, 0, halo_swap_state%halo_depth .ne. 1)
             else if (j==3) then
                source_index = current_state%local_grid%local_domain_end_index(i) - &
                     merge(1, 0, halo_swap_state%halo_depth .ne. 1)
             else if (j==4) then
                source_index = current_state%local_grid%local_domain_end_index(i)
             end if
             
             pid_location = get_pid_neighbour_location(halo_swap_state%halo_swap_neighbours, &
                  current_state%local_grid%neighbours(i,j), &
                  halo_swap_state%number_distinct_neighbours)

             if (present(source_data)) then
                call copy_fields_to_halo_buffer(current_state, &
                     halo_swap_state%halo_swap_neighbours(pid_location), i, source_index, &
                     pid_location, current_page, source_data)
             else
                call copy_fields_to_halo_buffer(current_state, &
                     halo_swap_state%halo_swap_neighbours(pid_location), i, source_index, &
                     pid_location, current_page)
             end if
             ! call log_log(LOG_DEBUG, "PID ="//trim(conv_to_String(&
             !current_state%parallel%my_rank))//" source_index = "//&
             !      trim(conv_to_string(source_index))//" PID location ="//trim(&
             !conv_to_string(pid_location))//&
             !      " i= "//trim(conv_to_string(i))//" j="//trim(conv_to_string(j)))
          end if
       end do
    end do
    
    if (present(copy_corner_fields_to_halo_buffer)) then
       current_page(:)=1
       do j = 1, size(current_state%local_grid%corner_neighbours, 1)
          if (current_state%parallel%my_rank .ne. &
               current_state%local_grid%corner_neighbours(j,1)) then
             x_source_index = merge(current_state%local_grid%local_domain_start_index(X_INDEX)+1,&
                  current_state%local_grid%local_domain_end_index(X_INDEX)-1, j==1 .or. j==3)
             y_source_index  =merge(current_state%local_grid%local_domain_start_index(Y_INDEX)+1,&
                  current_state%local_grid%local_domain_end_index(Y_INDEX)-1, j==1 .or. j==2)
             pid_location = get_pid_neighbour_location(halo_swap_state%halo_swap_neighbours, &
                  current_state%local_grid%corner_neighbours(j,1), &
                  halo_swap_state%number_distinct_neighbours)
             if (present(source_data)) then
                call copy_corner_fields_to_halo_buffer(current_state, &
                     halo_swap_state%halo_swap_neighbours(pid_location), j, &
                     x_source_index, y_source_index, pid_location, current_page, source_data)
             else
                call copy_corner_fields_to_halo_buffer(current_state, &
                     halo_swap_state%halo_swap_neighbours(pid_location), j, &
                     x_source_index, y_source_index, pid_location, current_page)
             end if
          end if
       end do
    end if
    
    do i=1,halo_swap_state%number_distinct_neighbours
       if (halo_swap_state%halo_swap_neighbours(i)%send_size .gt. 0) then
          call mpi_isend(halo_swap_state%halo_swap_neighbours(i)%send_halo_buffer, &
               halo_swap_state%halo_swap_neighbours(i)%send_size, PRECISION_TYPE, &
               halo_swap_state%halo_swap_neighbours(i)%pid, 0, &
               current_state%parallel%neighbour_comm, &
               halo_swap_state%send_requests(request_number), ierr)
          request_number = request_number+1
       end if
       if (halo_swap_state%halo_swap_neighbours(i)%send_corner_size .gt. 0) then
          call mpi_isend(halo_swap_state%halo_swap_neighbours(i)%send_corner_buffer, &
               halo_swap_state%halo_swap_neighbours(i)%send_corner_size, PRECISION_TYPE, &
               halo_swap_state%halo_swap_neighbours(i)%pid, 0, &
               current_state%parallel%neighbour_comm, &
               halo_swap_state%send_requests(request_number), ierr)
          request_number = request_number+1
       end if
    end do
  end subroutine send_all_halos
  
  !> Copies the received data (held in buffers) from neighbours into the correct halo location
  !! in the prognostic fields
  !! @param current_state The current model state
  !! @param halo_swap_state The halo swapping state
  !! @param copy_halo_buffer_to_field Procedure pointer which copies the halo buffer into the data field
  !! @param copy_halo_buffer_to_corner Optional procedure pointer to copy halo buffer into field corners
  !! @param source_data Optional source data which is read from into send buffers and written into by receieve buffers
  subroutine copy_buffer_data_for_prognostics(current_state, halo_swap_state, copy_halo_buffer_to_field, &
       copy_halo_buffer_to_corner, source_data)
    type(model_state_type), intent(inout) :: current_state
    type(halo_communication_type), intent(inout) :: halo_swap_state
    procedure(copy_halo_buffer_to_field_proc_interface) :: copy_halo_buffer_to_field
    procedure(copy_halo_buffer_to_corner_proc_interface), optional :: copy_halo_buffer_to_corner
    type(field_data_wrapper_type), dimension(:), intent(in), optional :: source_data

    integer :: i, j, hstart, hend, pid_location, target_index, x_target_index, &
         y_target_index, current_page(halo_swap_state%number_distinct_neighbours)

    hstart=merge(2, 1, halo_swap_state%halo_depth==1)
    hend=merge(3, 4, halo_swap_state%halo_depth==1)

    current_page(:)=1
    do i=2, 3
      do j=hstart, hend
        if (current_state%parallel%my_rank .ne. current_state%local_grid%neighbours(i,j)) then
          if (j==halo_swap_state%cell_match(i, 1)) then
            target_index=1
          else if (j==halo_swap_state%cell_match(i, 2)) then
            target_index=2
          else if (j==halo_swap_state%cell_match(i, 3)) then
            target_index=current_state%local_grid%local_domain_end_index(i)+1
          else if (j==halo_swap_state%cell_match(i, 4)) then
            target_index=current_state%local_grid%local_domain_end_index(i)+2
          end if
          pid_location=get_pid_neighbour_location(halo_swap_state%halo_swap_neighbours, &
               current_state%local_grid%neighbours(i,j),&
               halo_swap_state%number_distinct_neighbours)
          if (present(source_data)) then
            call copy_halo_buffer_to_field(current_state, &
                 halo_swap_state%halo_swap_neighbours(pid_location), i, target_index,&
                 pid_location, current_page, source_data)
          else
            call copy_halo_buffer_to_field(current_state, &
                 halo_swap_state%halo_swap_neighbours(pid_location), i, target_index,&
                 pid_location, current_page)
          end if
        end if
      end do
    end do
    if (present(copy_halo_buffer_to_corner)) then
      current_page(:)=1
      do j=size(current_state%local_grid%corner_neighbours, 1),1,-1
        if (current_state%parallel%my_rank .ne. &
             current_state%local_grid%corner_neighbours(j,1)) then
          x_target_index=merge(2, current_state%local_grid%local_domain_end_index(X_INDEX)+1,&
               j==1 .or. j==3)
          y_target_index=merge(2, current_state%local_grid%local_domain_end_index(Y_INDEX)+1, &
               j==1 .or. j==2)
          pid_location=get_pid_neighbour_location(halo_swap_state%halo_swap_neighbours, &
               current_state%local_grid%corner_neighbours(j,1), &
               halo_swap_state%number_distinct_neighbours)
          if (present(source_data)) then
            call copy_halo_buffer_to_corner(current_state, &
                 halo_swap_state%halo_swap_neighbours(pid_location), j, &
                 x_target_index, y_target_index, pid_location, current_page, source_data)
          else
            call copy_halo_buffer_to_corner(current_state, &
                 halo_swap_state%halo_swap_neighbours(pid_location), j, &
                 x_target_index, y_target_index, pid_location, current_page)
          end if
        end if
      end do
    end if
  end subroutine copy_buffer_data_for_prognostics

  !> Performs a local data copy for corners when the neighbour is local (me)
  !! @param my_rank My global rank
  !! @param local_grid The local grid
  !! @param field_data The field data top copy into
  subroutine perform_local_data_copy_for_corners(my_rank, local_grid, field_data)
    type(local_grid_type), intent(inout) :: local_grid
    real(kind=DEFAULT_PRECISION), dimension(:,:,:), intent(inout) :: field_data
    integer, intent(in) :: my_rank

    integer :: i, y_source_index, x_source_index, y_target_index, x_target_index

    do i=1,size(local_grid%corner_neighbours, 1)
      if (my_rank .eq. local_grid%corner_neighbours(i,1)) then
        if (i==1) then
          y_source_index=local_grid%local_domain_end_index(Y_INDEX)-1
          x_source_index=local_grid%local_domain_end_index(X_INDEX)-1
          y_target_index=1
          x_target_index=1
        else if (i==2) then
          y_source_index=local_grid%local_domain_end_index(Y_INDEX)-1
          x_source_index=local_grid%local_domain_start_index(X_INDEX)
          y_target_index=1
          x_target_index=local_grid%local_domain_end_index(X_INDEX)+1
        else if (i==3) then
          y_source_index=local_grid%local_domain_start_index(Y_INDEX)
          x_source_index=local_grid%local_domain_end_index(X_INDEX)-1
          y_target_index=local_grid%local_domain_end_index(Y_INDEX)+1
          x_target_index=1
        else if (i==4) then
          y_source_index=local_grid%local_domain_start_index(Y_INDEX)
          x_source_index=local_grid%local_domain_start_index(X_INDEX)
          y_target_index=local_grid%local_domain_end_index(Y_INDEX)+1
          x_target_index=local_grid%local_domain_end_index(X_INDEX)+1
        end if

        field_data(local_grid%local_domain_start_index(Z_INDEX):&
             local_grid%local_domain_end_index(Z_INDEX),&
             y_target_index, x_target_index)=&
             field_data(local_grid%local_domain_start_index(Z_INDEX):&
             local_grid%local_domain_end_index(Z_INDEX),y_source_index, x_source_index)
        
        field_data(local_grid%local_domain_start_index(Z_INDEX):&
             local_grid%local_domain_end_index(Z_INDEX),&
             y_target_index+1, x_target_index)=&
             field_data(local_grid%local_domain_start_index(Z_INDEX):&
             local_grid%local_domain_end_index(Z_INDEX),y_source_index+1, x_source_index)

        field_data(local_grid%local_domain_start_index(Z_INDEX):&
             local_grid%local_domain_end_index(Z_INDEX),&
             y_target_index, x_target_index+1)=&
             field_data(local_grid%local_domain_start_index(Z_INDEX):&
             local_grid%local_domain_end_index(Z_INDEX),y_source_index, x_source_index+1)
        
        field_data(local_grid%local_domain_start_index(Z_INDEX):&
             local_grid%local_domain_end_index(Z_INDEX),&
             y_target_index+1, x_target_index+1)=&
             field_data(local_grid%local_domain_start_index(Z_INDEX):&
             local_grid%local_domain_end_index(Z_INDEX),y_source_index+1, x_source_index+1)       
      end if      
    end do    
  end subroutine perform_local_data_copy_for_corners  

  !> Performs a local data copy for a specific dimension of a prognostic field
  !! @param dim The dimension
  !! @param my_rank My process Id
  !! @param local_grid The local grid description
  !! @param field_data The field data that we are going to copy from and to
  subroutine perform_local_data_copy_for_dimension(dim, my_rank, halo_depth, local_grid, &
       field_data)
    type(local_grid_type), intent(inout) :: local_grid
    real(kind=DEFAULT_PRECISION), dimension(:,:,:), intent(inout) :: field_data
    integer, intent(in) :: dim, my_rank, halo_depth

    integer i, target_index, source_index, hstart, hend

    hstart=merge(2,1, halo_depth==1)
    hend=merge(3,4, halo_depth==1)

    do i=hstart, hend
      if (local_grid%neighbours(dim,i) .eq. my_rank) then 
        if (i==1) then
          target_index=1
          source_index=local_grid%local_domain_end_index(dim)-1
        else if (i==2) then 
          target_index=2
          source_index=local_grid%local_domain_end_index(dim)
        else if (i==3) then 
          target_index=local_grid%local_domain_end_index(dim)+1
          source_index=local_grid%local_domain_start_index(dim)
        else if (i==4) then
          target_index=local_grid%local_domain_end_index(dim)+2
          source_index=local_grid%local_domain_start_index(dim)+1
        end if
        if (dim == X_INDEX) then
          field_data(local_grid%local_domain_start_index(Z_INDEX):&
               local_grid%local_domain_end_index(Z_INDEX),&
               local_grid%local_domain_start_index(Y_INDEX):&
               local_grid%local_domain_end_index(Y_INDEX), target_index) = &
               field_data(local_grid%local_domain_start_index(Z_INDEX):&
               local_grid%local_domain_end_index(Z_INDEX),&
               local_grid%local_domain_start_index(Y_INDEX):&
               local_grid%local_domain_end_index(Y_INDEX), source_index)
        else
          field_data(local_grid%local_domain_start_index(Z_INDEX):&
               local_grid%local_domain_end_index(Z_INDEX),&
               target_index, local_grid%local_domain_start_index(X_INDEX):&
               local_grid%local_domain_end_index(X_INDEX)) = &
               field_data(local_grid%local_domain_start_index(Z_INDEX):&
               local_grid%local_domain_end_index(Z_INDEX),&
               source_index, local_grid%local_domain_start_index(X_INDEX):&
               local_grid%local_domain_end_index(X_INDEX))
        end if
      end if
    end do
  end subroutine perform_local_data_copy_for_dimension

  !> Retrieves whether we have the same neighbours for L and R halo swaps in each dimension
  !!
  !! This is because if we do then the data storage needs to be reversed - i.e. the first two 
  !! messages sent to this process are Left and the second two are Right, but due to wrapping
  !! around, the L messages need to be put on the right halo and R messages on the left halo
  !! @param neighbours The neighbouring processes in each dimension
  function retrieve_same_neighbour_information(local_grid)
    type(local_grid_type), intent(inout) :: local_grid
!    integer, dimension(:,:), intent(in) :: neighbours
    logical, dimension(3) :: retrieve_same_neighbour_information

    integer :: i, nd1, nd2
    
    retrieve_same_neighbour_information=(/.true., .true., .true./)

    ! halo_size in X and Y are the same, therefore it does not matter which one we take
    ! we multiply by 2 since there are 2 sides Up&Down or Left&Right
    do i = 1,local_grid%halo_size(Y_INDEX)*2
       if (i==1) then
          nd1=local_grid%neighbours(Y_INDEX,i)
          nd2=local_grid%neighbours(X_INDEX,i)
       else
          if (local_grid%neighbours(Y_INDEX,i) .ne. nd1) &
               retrieve_same_neighbour_information(Y_INDEX) = .false.
          if (local_grid%neighbours(X_INDEX,i) .ne. nd2) &
               retrieve_same_neighbour_information(X_INDEX) = .false.
       end if
    end do
  end function retrieve_same_neighbour_information

  !> Returns whether or not a specific process id has already been "seen" by searching a list
  !! of already seen process ids. 
  !! @param temp_neighbour_pids The process Ids already seen
  !! @param pid The PID to search for
  logical function has_pid_already_been_seen(temp_neighbour_pids, pid)
    integer, intent(in) :: pid, temp_neighbour_pids(8)

    integer :: i

    has_pid_already_been_seen=.true.
    do i=1,8
      if (temp_neighbour_pids(i) == pid) return
      if (temp_neighbour_pids(i) == -1) then
        has_pid_already_been_seen=.false.
        return
      end if      
    end do
    has_pid_already_been_seen=.false.
  end function has_pid_already_been_seen

  !> Given the process id of a neighbour this determines the location in the data structure
  !! of corresponding data for that. Note that this is an O(n) operation so ideally call once
  !! and reuse the results. If no process id is found then returns -1
  !! @param halo_swap_neighbours PIDs of the neighbouring processes to me
  !! @param pid The process id to find
  !! @param number_distinct_neighbours The number of distinct neighbours I have
  integer function get_pid_neighbour_location(halo_swap_neighbours, pid, &
       number_distinct_neighbours)
    type(neighbour_description_type), dimension(:), allocatable :: halo_swap_neighbours
    integer, intent(in) :: pid, number_distinct_neighbours

    integer :: i

    do i=1, number_distinct_neighbours
      if (halo_swap_neighbours(i)%pid == pid) then
        get_pid_neighbour_location = i
        return
      end if      
    end do
    ! Not found
    get_pid_neighbour_location=-1
  end function get_pid_neighbour_location

  !> A very common function, which returns a single field per halo cell which is used to halo
  !! swap just one field
  !! @param current_state The current model state
  integer function get_single_field_per_halo_cell(current_state)
    type(model_state_type), intent(inout) :: current_state    

    get_single_field_per_halo_cell=1
  end function get_single_field_per_halo_cell  
end module halo_communication_mod
