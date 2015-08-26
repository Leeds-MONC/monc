!> Pressure solver which uses a tridiagonal algorithm operating on the pressure terms in Fourier space. It uses the pencil FFT
!! module for 3D FFTs in pencil decomposition. These use FFTW for the actual FFT kernel.
module fftsolver_mod
  use datadefn_mod, only : DEFAULT_PRECISION
  use grids_mod, only : X_INDEX, Y_INDEX, Z_INDEX
  use state_mod, only : model_state_type
  use monc_component_mod, only : component_descriptor_type
  use pencil_fft_mod, only : initialise_pencil_fft, finalise_pencil_fft, perform_forward_3dfft, perform_backwards_3dfft
  use communication_types_mod, only : halo_communication_type, neighbour_description_type, field_data_wrapper_type
  use halo_communication_mod, only : copy_buffer_to_field, copy_field_to_buffer, perform_local_data_copy_for_field, &
       init_halo_communication, finalise_halo_communication, blocking_halo_swap, get_single_field_per_halo_cell
  use registry_mod, only : is_component_enabled
  use logging_mod, only : LOG_ERROR, log_master_log
  use mpi, only : MPI_REQUEST_NULL, MPI_STATUSES_IGNORE
  implicit none

#ifndef TEST_MODE
  private
#endif

  real(kind=DEFAULT_PRECISION), dimension(:), allocatable :: cos_x, cos_y
  real(kind=DEFAULT_PRECISION) :: PI
  real(kind=DEFAULT_PRECISION), dimension(:,:,:), allocatable :: pt
  integer :: fourier_space_sizes(3)
  type(halo_communication_type), save :: halo_swap_state

  public fftsolver_get_descriptor
contains

  !> Descriptor of this component for registration
  !! @returns The fft solver component descriptor
  type(component_descriptor_type) function fftsolver_get_descriptor()
    fftsolver_get_descriptor%name="fftsolver"
    fftsolver_get_descriptor%version=0.1
    fftsolver_get_descriptor%initialisation=>initialisation_callback
    fftsolver_get_descriptor%timestep=>timestep_callback
    fftsolver_get_descriptor%finalisation=>finalisation_callback
  end function fftsolver_get_descriptor

  !> This initialisation callback sets up the pencil fft module, allocates data for the fourier space pressure term (might be
  !! different size than p) and populates cos x and y
  subroutine initialisation_callback(current_state)
    type(model_state_type), target, intent(inout) :: current_state

    integer :: i, my_y_start, my_x_start

    if (.not. is_component_enabled(current_state%options_database, "diverr")) then
      call log_master_log(LOG_ERROR, "The FFT solver component requires the diverr component to be enabled")
    end if 

    PI=4.0_DEFAULT_PRECISION*atan(1.0_DEFAULT_PRECISION)

    fourier_space_sizes=initialise_pencil_fft(current_state, my_y_start, my_x_start)

    call init_halo_communication(current_state, get_single_field_per_halo_cell, halo_swap_state, 1, .false.)

    allocate(pt(fourier_space_sizes(Z_INDEX), fourier_space_sizes(Y_INDEX), fourier_space_sizes(X_INDEX)))

#ifdef U_ACTIVE
    allocate(cos_x(fourier_space_sizes(X_INDEX)))
    do i=1, fourier_space_sizes(X_INDEX)
      cos_x(i)=(2.0_DEFAULT_PRECISION/(current_state%global_grid%configuration%horizontal%dx*&
           current_state%global_grid%configuration%horizontal%dx))&
           *(cos(2.0_DEFAULT_PRECISION*PI*real(((i-1)+(my_x_start-1))/2, kind=DEFAULT_PRECISION)/&
           real(current_state%global_grid%size(X_INDEX), kind=DEFAULT_PRECISION))-1.0_DEFAULT_PRECISION)
    end do    
#endif

#ifdef V_ACTIVE
    allocate(cos_y(fourier_space_sizes(Y_INDEX)))
    do i=1, fourier_space_sizes(Y_INDEX)
      cos_y(i)=(2.0_DEFAULT_PRECISION/(current_state%global_grid%configuration%horizontal%dy*&
           current_state%global_grid%configuration%horizontal%dy))&
           *(cos(2.0_DEFAULT_PRECISION*PI*real(((i-1)+(my_y_start-1))/2, kind=DEFAULT_PRECISION)/&
           real(current_state%global_grid%size(Y_INDEX), kind=DEFAULT_PRECISION))-1.0_DEFAULT_PRECISION)
    end do    
#endif
    current_state%psrce_x_hs_send_request=MPI_REQUEST_NULL
    current_state%psrce_y_hs_send_request=MPI_REQUEST_NULL
    current_state%psrce_x_hs_recv_request=MPI_REQUEST_NULL
    current_state%psrce_y_hs_recv_request=MPI_REQUEST_NULL 
  end subroutine initialisation_callback  

  !> Timestep call back, which will transform to Fourier space, do a tridiagonal solve and then back into time space
  !! @param current_state The current model state
  subroutine timestep_callback(current_state)
    type(model_state_type), target, intent(inout) :: current_state

    integer :: start_loc(3), end_loc(3), i

    do i=1,3
      start_loc(i)=current_state%local_grid%local_domain_start_index(i)
      end_loc(i)=current_state%local_grid%local_domain_end_index(i)
    end do

    call complete_psrce_calculation(current_state, current_state%local_grid%halo_size(Y_INDEX), &
         current_state%local_grid%halo_size(X_INDEX))

#ifdef FFT_TEST_MODE
    current_state%p%data=real(current_state%parallel%my_rank, kind=8) + 1.d0
#endif

    call perform_forward_3dfft(current_state, current_state%p%data(start_loc(Z_INDEX):end_loc(Z_INDEX), &
         start_loc(Y_INDEX):end_loc(Y_INDEX), start_loc(X_INDEX):end_loc(X_INDEX)), pt)

#ifndef FFT_TEST_MODE
    call tridiagonal_solver(current_state, pt, fourier_space_sizes)
#endif
    
    ! Here it is complex space, distributed as needed. The other option is real space but with +1 for the last process in Y
    ! as require that complex number - not sure the best approach. We can extract the values here from complex anyway so
    ! worth trying that first. Operating on the size of pt rather than grid sizes. Downside of real space xform is then we have +1
    ! in the forwards transformation    
    call perform_backwards_3dfft(current_state, pt, current_state%p%data(start_loc(Z_INDEX):end_loc(Z_INDEX), &
         start_loc(Y_INDEX):end_loc(Y_INDEX), start_loc(X_INDEX):end_loc(X_INDEX)))   

#ifdef FFT_TEST_MODE
    if (current_state%parallel%my_rank == 1) then
      write(*,*) current_state%parallel%my_rank, current_state%p%data(:,3,3)
    end if
#endif

    call blocking_halo_swap(current_state, halo_swap_state, copy_p_to_halo_buffer, &
         perform_local_data_copy_for_p, copy_halo_buffer_to_p)
    
  end subroutine timestep_callback  

  !> Called at MONC finalisation, will call to the pencil fft module to clean itself up and free the pressure term
  !! @param current_state The current model state
  subroutine finalisation_callback(current_state)
    type(model_state_type), target, intent(inout) :: current_state

    call finalise_pencil_fft(current_state%parallel%monc_communicator)
    deallocate(pt)
    if (allocated(cos_x)) deallocate(cos_x)
    if (allocated(cos_y)) deallocate(cos_y)
  end subroutine finalisation_callback

  !> Completes the psrce calculation by waiting on all outstanding psrce communications to complete and then combine the
  !! received values with the P field for U and V
  !! @param current_state The current model state
  !! @param y_halo_size The halo size in the Y dimension
  !! @param x_halo_size The halo size in the X dimension
  subroutine complete_psrce_calculation(current_state, y_halo_size, x_halo_size)
    type(model_state_type), target, intent(inout) :: current_state
    integer, intent(in) :: y_halo_size, x_halo_size

    integer :: ierr, combined_handles(2), i, j, k

    combined_handles(1)=current_state%psrce_x_hs_recv_request
    combined_handles(2)=current_state%psrce_y_hs_recv_request
    call mpi_waitall(2, combined_handles, MPI_STATUSES_IGNORE, ierr)

    do j=current_state%local_grid%local_domain_start_index(Y_INDEX), current_state%local_grid%local_domain_end_index(Y_INDEX)
      do k=2,current_state%local_grid%size(Z_INDEX)
#ifdef U_ACTIVE
        current_state%p%data(k,j,x_halo_size+1)=current_state%p%data(k,j,x_halo_size+1)-&
               current_state%psrce_recv_buffer_x(k-1,j-x_halo_size)
#endif
#ifdef V_ACTIVE
        if (j .gt. y_halo_size+1) current_state%p%data(k, j, x_halo_size+1)=current_state%p%data(k, j, x_halo_size+1)-&
             current_state%global_grid%configuration%horizontal%cy * current_state%sv%data(k, j-1, x_halo_size+1)
#endif
      end do
    end do
      
#ifdef V_ACTIVE
    do i=current_state%local_grid%local_domain_start_index(X_INDEX), current_state%local_grid%local_domain_end_index(X_INDEX)
      do k=2,current_state%local_grid%size(Z_INDEX)
        current_state%p%data(k,y_halo_size+1,i)=current_state%p%data(k,y_halo_size+1,i)-&
             current_state%psrce_recv_buffer_y(k-1,i-y_halo_size)
      end do
    end do
#endif

    combined_handles(1)=current_state%psrce_x_hs_send_request
    combined_handles(2)=current_state%psrce_y_hs_send_request
    call mpi_waitall(2, combined_handles, MPI_STATUSES_IGNORE, ierr)
  end subroutine complete_psrce_calculation 
  
  !> The tridiagonal solver which runs in Fourier space on the pressure terms. Note that because we are going
  !! from n reals to n/2+1 complex numbers (and into their reals for this solver) then there might be more data 
  !! in the pressure term than in p. Therefore use the space sizes which are provided as an argument to this
  !! procedure rather than the current state local grid dimensions
  !! @param current_state The current model state
  !! @param pressure_term The pressure term in Fourier space
  !! @param fourier_space_sizes The size of the pressure term in each dimension
  subroutine tridiagonal_solver(current_state, pressure_term, fourier_space_sizes)
    type(model_state_type), target, intent(inout) :: current_state
    real(kind=DEFAULT_PRECISION), dimension(:,:,:), intent(inout) :: pressure_term
    integer, intent(in) :: fourier_space_sizes(3)

    integer :: i, j, k, j_start
    real(kind=DEFAULT_PRECISION) :: b(fourier_space_sizes(Z_INDEX)), b1(fourier_space_sizes(Z_INDEX)),   &
        s(fourier_space_sizes(Z_INDEX)), s1(fourier_space_sizes(Z_INDEX)), cij
    
    do i=1,fourier_space_sizes(X_INDEX)
      j_start=merge(3, 1, current_state%local_grid%start(Y_INDEX) ==1 .and. current_state%local_grid%start(X_INDEX) == 1 &
           .and. i .lt. 3)
      do j=j_start,fourier_space_sizes(Y_INDEX)

        cij=cos_x(i)+cos_y(j)        
        b(2)=cij-current_state%global_grid%configuration%vertical%cza(2)
        b1(2)=1.0_DEFAULT_PRECISION/b(2)
        s1(2)=pressure_term(2,j,i)

        do k=3,fourier_space_sizes(Z_INDEX)
           b(k)=cij+current_state%global_grid%configuration%vertical%czg(k)
           b1(k)=1.0_DEFAULT_PRECISION/(b(k)-current_state%global_grid%configuration%vertical%czh(k)*b1(k-1))
           s1(k)=pressure_term(k,j,i)-current_state%global_grid%configuration%vertical%czb(k)*s1(k-1)*b1(k-1)
        end do

        pressure_term(fourier_space_sizes(Z_INDEX),j,i)=s1(fourier_space_sizes(Z_INDEX))* b1(fourier_space_sizes(Z_INDEX))

        do k=fourier_space_sizes(Z_INDEX)-1,2,-1
          pressure_term(k,j,i)=(s1(k)-current_state%global_grid%configuration%vertical%cza(k)*&
               pressure_term(k+1,j,i))*b1(k)
        end do                
      end do
    end do    

    ! Handle the zero wavenumber, which will only be on processes where the X and Y  dimension starts at 1
    if (current_state%local_grid%start(Y_INDEX) == 1 .and. current_state%local_grid%start(X_INDEX) == 1) then
      do i=1,2
        do j=1,2
          s(2)=pressure_term(2,j,i)
          pressure_term(1,j,i)=0.0_DEFAULT_PRECISION
          pressure_term(2,j,i)=0.0_DEFAULT_PRECISION
          do k=3, fourier_space_sizes(Z_INDEX)
            s(k)=pressure_term(k,j,i)
            pressure_term(k,j,i)=(s(k-1)-current_state%global_grid%configuration%vertical%czg(k-1)*pressure_term(k-1,j,i)-&
                 current_state%global_grid%configuration%vertical%czb(k-1)*pressure_term(k-2,j,i))/&
                 current_state%global_grid%configuration%vertical%cza(k-1)
          end do          
        end do        
      end do           
    end if
  end subroutine tridiagonal_solver

  !> Copies the p field data to halo buffers for a specific process in a dimension and halo cell
  !! @param current_state The current model state
  !! @param neighbour_descriptions Description of the neighbour halo swapping status
  !! @param dim Dimension to copy from
  !! @param source_index The source index of the dimension we are reading from in the prognostic field
  !! @param pid_location Location of the neighbouring process in the local stored data structures
  !! @param current_page The current (next) buffer page to copy into
  !! @param source_data Optional source data which is read from
  subroutine copy_p_to_halo_buffer(current_state, neighbour_description, dim, source_index, &
       pid_location, current_page, source_data)
    type(model_state_type), intent(inout) :: current_state
    integer, intent(in) :: dim, pid_location, source_index
    integer, intent(inout) :: current_page(:)
    type(neighbour_description_type), intent(inout) :: neighbour_description
    type(field_data_wrapper_type), dimension(:), intent(in), optional :: source_data

    call copy_field_to_buffer(current_state%local_grid, neighbour_description%send_halo_buffer, current_state%p%data, &
         dim, source_index, current_page(pid_location))

    current_page(pid_location)=current_page(pid_location)+1
  end subroutine copy_p_to_halo_buffer

  !! @param dim The dimension we receive for
  !! @param target_index The target index for the dimension we are receiving for
  !! @param neighbour_location The location in the local neighbour data stores of this neighbour
  !! @param current_page The current, next, halo swap page to read from (all previous have been read and copied already)
  !! @param source_data Optional source data which is written into
  subroutine copy_halo_buffer_to_p(current_state, neighbour_description, dim, target_index, &
       neighbour_location, current_page, source_data)
    type(model_state_type), intent(inout) :: current_state
    integer, intent(in) :: dim, target_index, neighbour_location
    integer, intent(inout) :: current_page(:)
    type(neighbour_description_type), intent(inout) :: neighbour_description
    type(field_data_wrapper_type), dimension(:), intent(in), optional :: source_data

    call copy_buffer_to_field(current_state%local_grid, neighbour_description%recv_halo_buffer, current_state%p%data, &
         dim, target_index, current_page(neighbour_location))

    current_page(neighbour_location)=current_page(neighbour_location)+1
  end subroutine copy_halo_buffer_to_p

  !> Does local data copying for P variable halo swap
  !! @param current_state The current model state_mod
  !! @param source_data Optional source data which is written into
  subroutine perform_local_data_copy_for_p(current_state, halo_depth, involve_corners, source_data)
    type(model_state_type), intent(inout) :: current_state
    integer, intent(in) :: halo_depth
    logical, intent(in) :: involve_corners
    type(field_data_wrapper_type), dimension(:), intent(in), optional :: source_data

    call perform_local_data_copy_for_field(current_state%p%data, current_state%local_grid, &
         current_state%parallel%my_rank, halo_depth, involve_corners)
  end subroutine perform_local_data_copy_for_p
end module fftsolver_mod
