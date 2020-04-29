!>  This sets the lower boundary conditions for theta and the q variables
module lowerbc_mod
  use monc_component_mod, only : component_descriptor_type
  use state_mod, only : FORWARD_STEPPING, PRESCRIBED_SURFACE_FLUXES, PRESCRIBED_SURFACE_VALUES, &
   model_state_type
  use grids_mod, only : Z_INDEX, Y_INDEX, X_INDEX, vertical_grid_configuration_type
  use datadefn_mod, only : DEFAULT_PRECISION, PRECISION_TYPE
  use prognostics_mod, only : prognostic_field_type
  use science_constants_mod, only : von_karman_constant, smallp, alphah, betah, betam, pi, &
       z0, z0th, convective_limit, gammah, gammam
  use logging_mod, only : LOG_ERROR, LOG_WARN, log_log
  use registry_mod, only : is_component_enabled
  use logging_mod, only : LOG_ERROR, log_master_log
  use q_indices_mod, only: get_q_index, standard_q_names
  use mpi, only: MPI_REQUEST_NULL, MPI_STATUSES_IGNORE
  implicit none

#ifndef TEST_MODE
  private
#endif

  integer, parameter :: CONVERGENCE_SUCCESS=1, CONVERGENCE_RICHARDSON_TOO_LARGE=2, CONVERGENCE_FAILURE=3
   
  real(kind=DEFAULT_PRECISION), parameter :: smth = 0.05_DEFAULT_PRECISION,& ! Smoothing between iterations
        tolm=1.0E-4_DEFAULT_PRECISION,  tolt=1.0E-4_DEFAULT_PRECISION ! Convergence tollerance for u and t star

  real(kind=DEFAULT_PRECISION) :: tstrcona, rhmbc, ddbc, ddbc_x4, eecon, r2ddbc, rcmbc, tstrconb, &
       x4con, xx0con, y2con, yy0con, viscous_courant_coefficient

  real(kind=DEFAULT_PRECISION), dimension(:,:,:), allocatable :: x_wrapping_send_buffer, y_wrapping_send_buffer, &
       x_wrapping_recv_buffer, y_wrapping_recv_buffer

  integer :: iqv  ! index for vapour
  integer :: wrapping_comm_requests(4), y_wrapping_target_id, x_wrapping_target_id

  public lowerbc_get_descriptor
contains

  !> Descriptor of this component for registration
  !! @returns The diverr component descriptor
  type(component_descriptor_type) function lowerbc_get_descriptor()
    lowerbc_get_descriptor%name="lower_bc"
    lowerbc_get_descriptor%version=0.1
    lowerbc_get_descriptor%initialisation=>initialisation_callback
    lowerbc_get_descriptor%timestep=>timestep_callback
    lowerbc_get_descriptor%finalisation=>finalisation_callback
  end function lowerbc_get_descriptor
  
  subroutine initialisation_callback(current_state)
    type(model_state_type), target, intent(inout) :: current_state

    real(kind=DEFAULT_PRECISION) :: bhbc
    integer :: num_wrapped_fields

    ! Adhill - this check is only required so that the vis_ and diff_coefficients
    !          are allocated in their respective components
    if (.not. is_component_enabled(current_state%options_database, "diffusion")) then
      call log_master_log(LOG_ERROR, "Lowerbc requires the diffusion component to be enabled")
    end if    
    if (.not. is_component_enabled(current_state%options_database, "viscosity")) then
      call log_master_log(LOG_ERROR, "Lowerbc requires the viscosity component to be enabled")
    end if

    call allocate_applicable_fields(current_state)

    wrapping_comm_requests=MPI_REQUEST_NULL

    num_wrapped_fields=0
    if (current_state%th%active) num_wrapped_fields=1
    num_wrapped_fields=num_wrapped_fields+current_state%number_q_fields

    if (num_wrapped_fields .gt. 0) then
      if (current_state%parallel%my_coords(Y_INDEX) == 0  .or. &
           current_state%parallel%my_coords(Y_INDEX) == current_state%parallel%dim_sizes(Y_INDEX)-1) then        
        if (current_state%parallel%my_coords(Y_INDEX) == 0) then
          y_wrapping_target_id=current_state%local_grid%neighbours(Y_INDEX, 1)
        else
          y_wrapping_target_id=current_state%local_grid%neighbours(Y_INDEX, 3)
        end if
        if (current_state%parallel%my_rank .ne. y_wrapping_target_id) then
          allocate(y_wrapping_send_buffer(current_state%local_grid%size(X_INDEX)+4, 2, num_wrapped_fields), &
               y_wrapping_recv_buffer(current_state%local_grid%size(X_INDEX)+4, 2, num_wrapped_fields))
        end if
      end if

      if (current_state%parallel%my_coords(X_INDEX) == 0  .or. &
           current_state%parallel%my_coords(X_INDEX) == current_state%parallel%dim_sizes(X_INDEX)-1) then        
        if (current_state%parallel%my_coords(X_INDEX) == 0) then
          x_wrapping_target_id=current_state%local_grid%neighbours(X_INDEX, 1)
        else
          x_wrapping_target_id=current_state%local_grid%neighbours(X_INDEX, 3)
        end if
        if (current_state%parallel%my_rank .ne. x_wrapping_target_id) then
          allocate(x_wrapping_send_buffer(current_state%local_grid%size(Y_INDEX)+4, 2, num_wrapped_fields), &
               x_wrapping_recv_buffer(current_state%local_grid%size(Y_INDEX)+4, 2, num_wrapped_fields))
        end if
      end if
    end if

    viscous_courant_coefficient=1.0_DEFAULT_PRECISION/current_state%global_grid%configuration%vertical%dzn(2)**2+&
         1.0_DEFAULT_PRECISION/(current_state%global_grid%configuration%horizontal%dx*&
         current_state%global_grid%configuration%horizontal%dx)+&
         1.0_DEFAULT_PRECISION/(current_state%global_grid%configuration%horizontal%dy*&
         current_state%global_grid%configuration%horizontal%dy)
    
    if ( current_state%use_surface_boundary_conditions .and.  &
         current_state%type_of_surface_boundary_conditions == PRESCRIBED_SURFACE_VALUES) then
      ! variables below are only required when PRESCRIBED_SURFACE_VALUES are used. 
       tstrcona=von_karman_constant/alphah*current_state%global_grid%configuration%vertical%zlogth
       bhbc=alphah*current_state%global_grid%configuration%vertical%zlogth
       rhmbc=betah*(current_state%global_grid%configuration%vertical%zn(2)+z0-z0th)/&
            (betam*current_state%global_grid%configuration%vertical%zn(2))
       ddbc=current_state%global_grid%configuration%vertical%zlogm*(bhbc-&
            rhmbc*current_state%global_grid%configuration%vertical%zlogm)
       ddbc_x4=4.*ddbc
       r2ddbc=0.5_DEFAULT_PRECISION/ddbc
       eecon=2.0_DEFAULT_PRECISION*rhmbc*current_state%global_grid%configuration%vertical%zlogm-bhbc
       rcmbc=1.0_DEFAULT_PRECISION/current_state%cmbc
       tstrconb=von_karman_constant/alphah
       x4con=gammam*(current_state%global_grid%configuration%vertical%zn(2)+z0)
       xx0con=gammam*z0
       y2con=gammah*(current_state%global_grid%configuration%vertical%zn(2)+z0)
       yy0con=gammah*z0th
    else
       tstrcona=0.0
       bhbc=0.0
       rhmbc=0.0
       ddbc=0.0
       ddbc_x4=0.0
       r2ddbc=0.0
       eecon=0.0
       rcmbc=0.0
       tstrconb=0.0
       x4con=0.0
       xx0con=0.0
       y2con=0.0
       yy0con=0.0
    endif 
    
    ! Determine vapour index
    if (.not. current_state%passive_q) then 
       iqv = get_q_index(standard_q_names%VAPOUR, 'lowerbc')
    endif

  end subroutine initialisation_callback

  subroutine finalisation_callback(current_state)
    type(model_state_type), target, intent(inout) :: current_state

    if (allocated(x_wrapping_send_buffer)) deallocate(x_wrapping_send_buffer)
    if (allocated(y_wrapping_send_buffer)) deallocate(y_wrapping_send_buffer)
    if (allocated(x_wrapping_recv_buffer)) deallocate(x_wrapping_recv_buffer)
    if (allocated(y_wrapping_recv_buffer)) deallocate(y_wrapping_recv_buffer)
  end subroutine finalisation_callback

  subroutine allocate_applicable_fields(current_state)
    type(model_state_type), target, intent(inout) :: current_state

    integer :: z_size, y_size, x_size, i

    z_size=current_state%local_grid%size(Z_INDEX) + current_state%local_grid%halo_size(Z_INDEX) * 2
    y_size=current_state%local_grid%size(Y_INDEX) + current_state%local_grid%halo_size(Y_INDEX) * 2
    x_size=current_state%local_grid%size(X_INDEX) + current_state%local_grid%halo_size(X_INDEX) * 2

    allocate(current_state%dis%data(z_size, y_size, x_size), &
         current_state%dis_th%data(z_size, y_size, x_size), current_state%disq(current_state%number_q_fields))

    do i=1,current_state%number_q_fields
      allocate(current_state%disq(i)%data(z_size, y_size, x_size))
    end do   
  end subroutine allocate_applicable_fields

  subroutine timestep_callback(current_state)
    type(model_state_type), target, intent(inout) :: current_state
    
    integer :: current_y_index, current_x_index

    current_y_index=current_state%column_local_y
    current_x_index=current_state%column_local_x

    if (current_state%field_stepping == FORWARD_STEPPING) then
      call compute_lower_boundary_conditions(current_state, current_y_index, current_x_index, &
           current_state%u, current_state%v, current_state%th, current_state%th, current_state%q, current_state%q)
    else
      if (current_state%scalar_stepping == FORWARD_STEPPING) then
        call compute_lower_boundary_conditions(current_state, current_y_index, current_x_index, &
           current_state%zu, current_state%zv, current_state%th, current_state%zth, current_state%q, current_state%zq)
      else
        call compute_lower_boundary_conditions(current_state, current_y_index, current_x_index, &
           current_state%zu, current_state%zv, current_state%zth, current_state%zth, current_state%zq, current_state%zq)
      end if
    end if    
  end subroutine timestep_callback  

  subroutine compute_lower_boundary_conditions(current_state, current_y_index, current_x_index, zu, zv, zth, th, zq, q)
    type(model_state_type), target, intent(inout) :: current_state
    type(prognostic_field_type), intent(inout) :: zu, zv, th, zth
    type(prognostic_field_type), dimension(:), intent(inout) :: q, zq
    integer, intent(in) :: current_y_index, current_x_index

    integer :: n
    real(kind=DEFAULT_PRECISION) :: horizontal_velocity_at_k2

    if (.not. current_state%use_viscosity_and_diffusion .or. .not. current_state%use_surface_boundary_conditions) then
       call simple_boundary_values(current_state, current_y_index, current_x_index, th, q)
    else
       ! one level in for the halo-column
       if (.not. (current_state%column_local_y .lt. current_state%local_grid%local_domain_start_index(Y_INDEX)-1 .or.&
         current_state%column_local_x .lt. current_state%local_grid%local_domain_start_index(X_INDEX)-1 .or.&
         current_state%column_local_y .gt. current_state%local_grid%local_domain_end_index(Y_INDEX)+1 .or.&
         current_state%column_local_x .gt. current_state%local_grid%local_domain_end_index(X_INDEX)+1)) then

          !if (.not. current_state%halo_column) then
          ! Include one halo to ensure that the halo is set in tvdadvection. This is done using the 
          ! logic from the timestep callback in tvdadvection in the timestep callback above
          horizontal_velocity_at_k2=0.0_DEFAULT_PRECISION
#ifdef U_ACTIVE
         horizontal_velocity_at_k2=(0.5_DEFAULT_PRECISION*(current_state%zu%data(2,current_y_index,current_x_index)+&
              zu%data(2,current_y_index,current_x_index-1))+ current_state%ugal)**2
#endif
#ifdef V_ACTIVE
         horizontal_velocity_at_k2=horizontal_velocity_at_k2+(0.5_DEFAULT_PRECISION*(zv%data(&
              2,current_y_index,current_x_index)+zv%data(2,current_y_index-1,current_x_index))+current_state%vgal)**2
#endif
         horizontal_velocity_at_k2=sqrt(horizontal_velocity_at_k2)+smallp      

         if (current_state%type_of_surface_boundary_conditions == PRESCRIBED_SURFACE_FLUXES) then
            call compute_using_fixed_surface_fluxes(current_state, current_y_index, current_x_index, &
                 horizontal_velocity_at_k2, th, q)
         else
            call compute_using_fixed_surface_temperature(current_state, current_y_index, current_x_index, &
               horizontal_velocity_at_k2, zth, th, zq, q)
         end if

         current_state%dis%data(1, current_y_index, current_x_index)=0.0_DEFAULT_PRECISION
         current_state%dis_th%data(1, current_y_index, current_x_index)=0.0_DEFAULT_PRECISION

         if (current_state%backscatter) then
            do n=1, current_state%number_q_fields
               current_state%disq(n)%data(1,current_y_index,current_x_index)=0.0_DEFAULT_PRECISION
            end do
        end if
        
        !-----------------------
        ! _return viscous number
        !-----------------------

        current_state%cvis=max(current_state%cvis,max(current_state%vis_coefficient%data(1, current_y_index, current_x_index),&
             current_state%diff_coefficient%data(1, current_y_index, current_x_index))*viscous_courant_coefficient)
        !            CVIS will be multiplied by DTM_X4 in TESTCFL
     else if (current_x_index == 1 .and. current_y_index == 1) then
        call register_async_wrapping_recv_requests(current_state)
     else if (current_x_index == current_state%local_grid%local_domain_end_index(X_INDEX)+&
             current_state%local_grid%halo_size(X_INDEX) .and. current_y_index == &
             current_state%local_grid%local_domain_end_index(Y_INDEX)+current_state%local_grid%halo_size(Y_INDEX)) then
        call complete_async_wrapping(current_state, zth, zq)
     end if
  end if
  end subroutine compute_lower_boundary_conditions

  !> Registers asynchronous wrapping recv requests as needed
  !! @param current_state The current model state
  subroutine register_async_wrapping_recv_requests(current_state)
    type(model_state_type), target, intent(inout) :: current_state

    integer :: ierr

    if (allocated(y_wrapping_recv_buffer)) then
      call mpi_irecv(y_wrapping_recv_buffer, size(y_wrapping_recv_buffer), PRECISION_TYPE, &
           y_wrapping_target_id, 0, current_state%parallel%neighbour_comm, wrapping_comm_requests(1), ierr)
    end if
    if (allocated(x_wrapping_recv_buffer)) then
      call mpi_irecv(x_wrapping_recv_buffer, size(x_wrapping_recv_buffer), PRECISION_TYPE, &
           x_wrapping_target_id, 0, current_state%parallel%neighbour_comm, wrapping_comm_requests(3), ierr)
    end if
  end subroutine register_async_wrapping_recv_requests
  
  !> Completes the asynchronous wrapping if required for periodic boundary conditions
  !! @param current_state The current model state
  !! @param zth Temperature field
  !! @param zq Q fields
  subroutine complete_async_wrapping(current_state, zth, zq)
    type(model_state_type), target, intent(inout) :: current_state
    type(prognostic_field_type), intent(inout) :: zth
    type(prognostic_field_type), dimension(:), intent(inout) :: zq

    integer :: ierr, n

    if (allocated(x_wrapping_send_buffer) .or. allocated(y_wrapping_send_buffer)) then
      if (allocated(y_wrapping_send_buffer)) then
        if (current_state%parallel%my_coords(Y_INDEX) == 0) then
          call package_y_wrapping_send_buffer(current_state, zth, zq, current_state%local_grid%local_domain_start_index(Y_INDEX),&
               current_state%local_grid%local_domain_start_index(Y_INDEX)+1)        
        else
          call package_y_wrapping_send_buffer(current_state, zth, zq, current_state%local_grid%local_domain_end_index(Y_INDEX)-1,&
               current_state%local_grid%local_domain_end_index(Y_INDEX))        
        end if        
        call mpi_isend(y_wrapping_send_buffer, size(y_wrapping_send_buffer), PRECISION_TYPE, &
             y_wrapping_target_id, 0, current_state%parallel%neighbour_comm, &
             wrapping_comm_requests(2), ierr)       
      end if
      if (allocated(x_wrapping_send_buffer)) then
        if (current_state%parallel%my_coords(X_INDEX) == 0) then
          call package_x_wrapping_send_buffer(current_state, zth, zq, current_state%local_grid%local_domain_start_index(X_INDEX),&
               current_state%local_grid%local_domain_start_index(X_INDEX)+1)        
        else
          call package_x_wrapping_send_buffer(current_state, zth, zq, current_state%local_grid%local_domain_end_index(X_INDEX)-1,&
               current_state%local_grid%local_domain_end_index(X_INDEX))        
        end if        
        call mpi_isend(x_wrapping_send_buffer, size(x_wrapping_send_buffer), PRECISION_TYPE, &
             x_wrapping_target_id, 0, current_state%parallel%neighbour_comm, &
             wrapping_comm_requests(4), ierr)        
      end if

      ! If send buffer is allocated then recv buffer is allocated, therefore only test the send buffer here and assume recv
      call mpi_waitall(4, wrapping_comm_requests, MPI_STATUSES_IGNORE, ierr)
      wrapping_comm_requests=MPI_REQUEST_NULL
      if (allocated(y_wrapping_recv_buffer)) then
        if (current_state%parallel%my_coords(Y_INDEX) == 0) then
          call unpackage_y_wrapping_recv_buffer(current_state, zth, zq, 1, 2)
        else
          call unpackage_y_wrapping_recv_buffer(current_state, zth, zq, &
               current_state%local_grid%local_domain_end_index(Y_INDEX)+1, &
               current_state%local_grid%local_domain_end_index(Y_INDEX)+2)          
        end if
      end if
      if (allocated(x_wrapping_recv_buffer)) then
        if (current_state%parallel%my_coords(X_INDEX) == 0) then
          call unpackage_x_wrapping_recv_buffer(current_state, zth, zq, 1, 2)
        else
          call unpackage_x_wrapping_recv_buffer(current_state, zth, zq, &
               current_state%local_grid%local_domain_end_index(X_INDEX)+1, &
               current_state%local_grid%local_domain_end_index(X_INDEX)+2)          
        end if
      end if
    end if 
    if (current_state%parallel%my_rank == y_wrapping_target_id) then
      if (current_state%th%active) then
        zth%data(1,1,:)=zth%data(1, current_state%local_grid%local_domain_end_index(Y_INDEX)-1, :)
        zth%data(1,2,:)=zth%data(1, current_state%local_grid%local_domain_end_index(Y_INDEX), :)
        zth%data(1,current_state%local_grid%local_domain_end_index(Y_INDEX)+1,:)=&
             zth%data(1, current_state%local_grid%local_domain_start_index(Y_INDEX), :)
        zth%data(1,current_state%local_grid%local_domain_end_index(Y_INDEX)+2,:)=&
             zth%data(1, current_state%local_grid%local_domain_start_index(Y_INDEX)+1, :)
      end if
      if (current_state%number_q_fields .gt. 0) then
        do n=1, current_state%number_q_fields
          zq(n)%data(1,1,:)=zq(n)%data(1, current_state%local_grid%local_domain_end_index(Y_INDEX)-1, :)
          zq(n)%data(1,2,:)=zq(n)%data(1, current_state%local_grid%local_domain_end_index(Y_INDEX), :)
          zq(n)%data(1,current_state%local_grid%local_domain_end_index(Y_INDEX)+1,:)=&
             zq(n)%data(1, current_state%local_grid%local_domain_start_index(Y_INDEX), :)
          zq(n)%data(1,current_state%local_grid%local_domain_end_index(Y_INDEX)+2,:)=&
             zq(n)%data(1, current_state%local_grid%local_domain_start_index(Y_INDEX)+1, :)
        end do
      end if
    end if

    if (current_state%parallel%my_rank == x_wrapping_target_id) then
      if (current_state%th%active) then
        zth%data(1,:,1)=zth%data(1,:,current_state%local_grid%local_domain_end_index(X_INDEX)-1)
        zth%data(1,:,2)=zth%data(1,:,current_state%local_grid%local_domain_end_index(X_INDEX))
        zth%data(1,:,current_state%local_grid%local_domain_end_index(X_INDEX)+1)=&
             zth%data(1,:,current_state%local_grid%local_domain_start_index(X_INDEX))
        zth%data(1,:,current_state%local_grid%local_domain_end_index(X_INDEX)+2)=&
             zth%data(1,:,current_state%local_grid%local_domain_start_index(X_INDEX)+1)
      end if
      if (current_state%number_q_fields .gt. 0) then
        do n=1, current_state%number_q_fields
          zq(n)%data(1,:,1)=zq(n)%data(1,:,current_state%local_grid%local_domain_end_index(X_INDEX)-1)
          zq(n)%data(1,:,2)=zq(n)%data(1,:,current_state%local_grid%local_domain_end_index(X_INDEX))
          zq(n)%data(1,:,current_state%local_grid%local_domain_end_index(X_INDEX)+1)=&
               zq(n)%data(1,:,current_state%local_grid%local_domain_start_index(X_INDEX))
          zq(n)%data(1,:,current_state%local_grid%local_domain_end_index(X_INDEX)+2)=&
               zq(n)%data(1,:,current_state%local_grid%local_domain_start_index(X_INDEX)+1)
        end do
      end if
    end if
  end subroutine complete_async_wrapping

  !> Packages theta and Q fields (if enabled) into the send buffer for Y
  !! @param current_state The current model state
  !! @param zth Temperature field
  !! @param zq Q fields
  !! @param first_y_index The first Y index to read from the data field
  !! @param second_y_index The second Y index to read from the data field
  subroutine package_y_wrapping_send_buffer(current_state, zth, zq, first_y_index, second_y_index)
    type(model_state_type), target, intent(inout) :: current_state
    type(prognostic_field_type), intent(inout) :: zth
    type(prognostic_field_type), dimension(:), intent(inout) :: zq
    integer, intent(in) :: first_y_index, second_y_index

    integer :: index_start, n

    index_start=0
    if (current_state%th%active) then
      y_wrapping_send_buffer(:,1,1)=zth%data(1,first_y_index,:)
      y_wrapping_send_buffer(:,2,1)=zth%data(1,second_y_index,:)
      index_start=index_start+1
    end if
    if (current_state%number_q_fields .gt. 0) then
      do n=1, current_state%number_q_fields
        y_wrapping_send_buffer(:,1,index_start+n)=zq(n)%data(1,first_y_index,:)
        y_wrapping_send_buffer(:,2,index_start+n)=zq(n)%data(1,second_y_index,:)
      end do
    end if
  end subroutine package_y_wrapping_send_buffer

  !> Packages theta and Q fields (if enabled) into the send buffer for X
  !! @param current_state The current model state
  !! @param zth Temperature field
  !! @param zq Q fields
  !! @param first_x_index The first X index to read from the data field
  !! @param second_x_index The second X index to read from the data field
  subroutine package_x_wrapping_send_buffer(current_state, zth, zq, first_x_index, second_x_index)
    type(model_state_type), target, intent(inout) :: current_state
    type(prognostic_field_type), intent(inout) :: zth
    type(prognostic_field_type), dimension(:), intent(inout) :: zq
    integer, intent(in) :: first_x_index, second_x_index

    integer :: index_start, n
    
    index_start=0
    if (current_state%th%active) then
      x_wrapping_send_buffer(:,1,1)=zth%data(1,:,first_x_index)
      x_wrapping_send_buffer(:,2,1)=zth%data(1,:,second_x_index)
      index_start=index_start+1
    end if
    if (current_state%number_q_fields .gt. 0) then
      do n=1, current_state%number_q_fields
        x_wrapping_send_buffer(:,1,index_start+n)= zq(n)%data(1,:,first_x_index)
        x_wrapping_send_buffer(:,2,index_start+n)= zq(n)%data(1,:,second_x_index)
      end do
    end if
  end subroutine package_x_wrapping_send_buffer

  !> Unpackages theta and Q fields from the receive buffer into the fields themselves (if enabled) for Y
  !! @param current_state The current model state
  !! @param zth Temperature field
  !! @param zq Q fields
  !! @param first_y_index The first Y index to read from the data field
  !! @param second_y_index The second Y index to read from the data field
  subroutine unpackage_y_wrapping_recv_buffer(current_state, zth, zq, first_y_index, second_y_index)
    type(model_state_type), target, intent(inout) :: current_state
    type(prognostic_field_type), intent(inout) :: zth
    type(prognostic_field_type), dimension(:), intent(inout) :: zq
    integer, intent(in) :: first_y_index, second_y_index

    integer :: index_start, n

    index_start=0
    if (current_state%th%active) then
      zth%data(1,first_y_index,:)=y_wrapping_recv_buffer(:,1,1)
      zth%data(1,second_y_index,:)=y_wrapping_recv_buffer(:,2,1)
      index_start=index_start+1
    end if
    if (current_state%number_q_fields .gt. 0) then
      do n=1, current_state%number_q_fields
        zq(n)%data(1,first_y_index,:)=y_wrapping_recv_buffer(:,1,index_start+n)
        zq(n)%data(1,second_y_index,:)=y_wrapping_recv_buffer(:,2,index_start+n)
      end do
    end if
  end subroutine unpackage_y_wrapping_recv_buffer

  !> Unpackages theta and Q fields from the receive buffer into the fields themselves (if enabled) for X
  !! @param current_state The current model state
  !! @param zth Temperature field
  !! @param zq Q fields
  !! @param first_x_index The first X index to read from the data field
  !! @param second_x_index The second X index to read from the data field
  subroutine unpackage_x_wrapping_recv_buffer(current_state, zth, zq, first_x_index, second_x_index)
    type(model_state_type), target, intent(inout) :: current_state
    type(prognostic_field_type), intent(inout) :: zth
    type(prognostic_field_type), dimension(:), intent(inout) :: zq
    integer, intent(in) :: first_x_index, second_x_index

    integer :: index_start, n

    index_start=0
    if (current_state%th%active) then
      zth%data(1,:,first_x_index)=x_wrapping_recv_buffer(:,1,1)
      zth%data(1,:,second_x_index)=x_wrapping_recv_buffer(:,2,1)
      index_start=index_start+1
    end if
    if (current_state%number_q_fields .gt. 0) then
      do n=1, current_state%number_q_fields
        zq(n)%data(1,:,first_x_index)=x_wrapping_recv_buffer(:,1,index_start+n)
        zq(n)%data(1,:,second_x_index)=x_wrapping_recv_buffer(:,2,index_start+n)
      end do
    end if
  end subroutine unpackage_x_wrapping_recv_buffer

  subroutine handle_convective_fluxes(current_state, current_y_index, current_x_index, horizontal_velocity_at_k2, th, q)
    type(model_state_type), target, intent(inout) :: current_state
    type(prognostic_field_type), intent(inout) :: th
    type(prognostic_field_type), dimension(:), intent(inout) :: q
    integer, intent(in) :: current_y_index, current_x_index
    real(kind=DEFAULT_PRECISION), intent(in) :: horizontal_velocity_at_k2

    integer :: n
    real(kind=DEFAULT_PRECISION) :: ustr

    ustr=look(current_state, horizontal_velocity_at_k2)

    current_state%vis_coefficient%data(1, current_y_index, current_x_index)=&
         current_state%global_grid%configuration%vertical%czn*ustr**2/ horizontal_velocity_at_k2
    current_state%diff_coefficient%data(1, current_y_index, current_x_index)=&
         (von_karman_constant*current_state%global_grid%configuration%vertical%czn*ustr/alphah)/&
         (current_state%global_grid%configuration%vertical%zlogth- 2.*log(&
         (1.+sqrt(1.+gammah*von_karman_constant*current_state%fbuoy*(current_state%global_grid%configuration%vertical%zn(2)+z0)&
         /ustr**3))/ (1.+sqrt(1.+gammah*von_karman_constant*current_state%fbuoy*z0th/ustr**3))))
    if (current_state%th%active) th%data(1, current_y_index, current_x_index)= &
         (current_state%surface_temperature_flux*current_state%global_grid%configuration%vertical%dzn(2)/&
         current_state%diff_coefficient%data(1, current_y_index, current_x_index))+th%data(2, current_y_index, current_x_index)-& 
         current_state%global_grid%configuration%vertical%thref(1)+& 
         current_state%global_grid%configuration%vertical%thref(2)

    ! Surface Flux of vapour
    if (current_state%number_q_fields .gt. 0) then
       n=iqv
       q(n)%data(1, current_y_index, current_x_index)=q(n)%data(2, current_y_index, current_x_index)+&
            current_state%surface_vapour_flux*current_state%global_grid%configuration%vertical%dzn(2)/&
            current_state%diff_coefficient%data(1, current_y_index, current_x_index)
    endif
  end subroutine handle_convective_fluxes

  real(kind=DEFAULT_PRECISION) function look(current_state, vel)
    type(model_state_type), target, intent(inout) :: current_state
    real(kind=DEFAULT_PRECISION), intent(in) :: vel     ! Horizontal speed at lowest model level

    real(kind=DEFAULT_PRECISION) :: lookup_real_posn
    integer :: lookup_int_posn

    lookup_real_posn=1.0_DEFAULT_PRECISION+real(current_state%lookup_table_entries-1)*&
         log(vel/current_state%velmin)*current_state%aloginv
    lookup_int_posn=int(lookup_real_posn)                                                         

    if (lookup_int_posn .ge. 1) then                                                       
      if (lookup_int_posn .lt. current_state%lookup_table_entries) then      ! Linear interpolation
        look=current_state%lookup_table_ustr(lookup_int_posn)+ (lookup_real_posn-real(lookup_int_posn))*&
             (current_state%lookup_table_ustr(lookup_int_posn+1)-current_state%lookup_table_ustr(lookup_int_posn))
      else     ! Near neutral
        look=vel*current_state%cneut
      end if
    else       ! Nearly free convection                                     
      look=vel**(-convective_limit)*current_state%cfc
    end if
  end function look

! check allocations and initialisations

  subroutine handle_neutral_fluxes(current_state, current_y_index, current_x_index, horizontal_velocity_at_k2, th, q)
    type(model_state_type), target, intent(inout) :: current_state
    type(prognostic_field_type), intent(inout) :: th
    type(prognostic_field_type), dimension(:), intent(inout) :: q
    integer, intent(in) :: current_y_index, current_x_index
    real(kind=DEFAULT_PRECISION) :: horizontal_velocity_at_k2

    integer :: n
    real(kind=DEFAULT_PRECISION) :: ustr
    
    ustr=horizontal_velocity_at_k2*current_state%global_grid%configuration%vertical%vk_on_zlogm
    current_state%vis_coefficient%data(1, current_y_index, current_x_index)=current_state%global_grid%configuration%vertical%czn*&
         ustr**2/horizontal_velocity_at_k2
    current_state%diff_coefficient%data(1, current_y_index, current_x_index)=&
         current_state%vis_coefficient%data(1, current_y_index, current_x_index)*&
         current_state%global_grid%configuration%vertical%zlogm/(alphah*current_state%global_grid%configuration%vertical%zlogth)
    if (current_state%th%active) th%data(1, current_y_index, current_x_index)=  &
         (current_state%surface_temperature_flux*current_state%global_grid%configuration%vertical%dzn(2)/&
         current_state%diff_coefficient%data(1, current_y_index, current_x_index))+th%data(2, current_y_index, current_x_index)-& 
         current_state%global_grid%configuration%vertical%thref(1)+& 
         current_state%global_grid%configuration%vertical%thref(2)

    ! Flux of vapour
    if (current_state%number_q_fields .gt. 0) then
       n=iqv
       q(n)%data(1, current_y_index, current_x_index)=q(n)%data(2, current_y_index, current_x_index)+&
            current_state%surface_vapour_flux*current_state%global_grid%configuration%vertical%dzn(2)/&
            current_state%diff_coefficient%data(1, current_y_index, current_x_index)
    endif
  end subroutine handle_neutral_fluxes

  subroutine handle_stable_fluxes(current_state, current_y_index, current_x_index, horizontal_velocity_at_k2, th, q)
    type(model_state_type), target, intent(inout) :: current_state
    type(prognostic_field_type), intent(inout) :: th
    type(prognostic_field_type), dimension(:), intent(inout) :: q
    integer, intent(in) :: current_y_index, current_x_index
    real(kind=DEFAULT_PRECISION) :: horizontal_velocity_at_k2

    real(kind=DEFAULT_PRECISION) :: ustr
    integer :: n

    ustr=horizontal_velocity_at_k2*current_state%global_grid%configuration%vertical%vk_on_zlogm
    if((current_state%fbuoy - 1.E-9_DEFAULT_PRECISION) .lt. -4.0_DEFAULT_PRECISION*&
         von_karman_constant**2*horizontal_velocity_at_k2**3/ (27.0_DEFAULT_PRECISION*betam*&
         current_state%global_grid%configuration%vertical%zn(2)*current_state%global_grid%configuration%vertical%zlogm**2)) then
      ! Too stable for turbulence
      current_state%vis_coefficient%data(1, current_y_index, current_x_index)=0.0_DEFAULT_PRECISION
      current_state%diff_coefficient%data(1, current_y_index, current_x_index)=0.0_DEFAULT_PRECISION
      ! Level 1 values of l_th and l_q set to be harmless for advection scheme
      if (current_state%th%active) th%data(1, current_y_index, current_x_index)=th%data(2, current_y_index, current_x_index)
      if (current_state%number_q_fields .gt. 0) then
        do n=1, current_state%number_q_fields
          q(n)%data(1, current_y_index, current_x_index)=q(n)%data(2, current_y_index, current_x_index)
        end do
      end if
    else
      ustr=ustr/3.0_DEFAULT_PRECISION*(1.0_DEFAULT_PRECISION-2.0_DEFAULT_PRECISION*&
           cos((acos(-27.0_DEFAULT_PRECISION*betam*von_karman_constant*current_state%global_grid%configuration%vertical%zn(2)*&
           current_state%fbuoy/(current_state%global_grid%configuration%vertical%zlogm*&
           2.0_DEFAULT_PRECISION*ustr**3)-1.0_DEFAULT_PRECISION)+ 2.0_DEFAULT_PRECISION*pi)/3.0_DEFAULT_PRECISION))
      current_state%vis_coefficient%data(1, current_y_index, current_x_index)=&
           current_state%global_grid%configuration%vertical%czn*ustr**2/horizontal_velocity_at_k2
      current_state%diff_coefficient%data(1, current_y_index, current_x_index)=&
           current_state%vis_coefficient%data(1, current_y_index, current_x_index)*&
           (current_state%global_grid%configuration%vertical%zlogm-betam*current_state%global_grid%configuration%vertical%zn(2)*&
           von_karman_constant*current_state%fbuoy/ustr**3)/(alphah*current_state%global_grid%configuration%vertical%zlogth-betah*&
           von_karman_constant*current_state%fbuoy* (current_state%global_grid%configuration%vertical%zn(2)+ z0-z0th)/ustr**3)
      if (current_state%th%active) th%data(1, current_y_index, current_x_index)= &
           (current_state%surface_temperature_flux*current_state%global_grid%configuration%vertical%dzn(2)/&
           current_state%diff_coefficient%data(1, current_y_index, current_x_index))+th%data(2, current_y_index, current_x_index)-& 
           current_state%global_grid%configuration%vertical%thref(1)+& 
           current_state%global_grid%configuration%vertical%thref(2)


      !Flux of vapour
      if (current_state%number_q_fields .gt. 0) then
         n=iqv
         q(n)%data(1, current_y_index, current_x_index)=q(n)%data(2, current_y_index, current_x_index)+&
              current_state%surface_vapour_flux*current_state%global_grid%configuration%vertical%dzn(2)/&
              current_state%diff_coefficient%data(1, current_y_index, current_x_index)
      endif
    end if
  end subroutine handle_stable_fluxes

  ! set surface_boundary_flux in init == FBUOY
  subroutine compute_using_fixed_surface_fluxes(current_state, current_y_index, current_x_index, horizontal_velocity_at_k2, th, q)
    type(model_state_type), target, intent(inout) :: current_state
    type(prognostic_field_type), intent(inout) :: th
    type(prognostic_field_type), dimension(:), intent(inout) :: q
    integer, intent(in) :: current_y_index, current_x_index
    real(kind=DEFAULT_PRECISION) :: horizontal_velocity_at_k2

    if (current_state%fbuoy .gt. 0.0_DEFAULT_PRECISION) then
      call handle_convective_fluxes(current_state, current_y_index, current_x_index, horizontal_velocity_at_k2, th, q)
    else if (current_state%fbuoy .eq. 0.0_DEFAULT_PRECISION) then
      call handle_neutral_fluxes(current_state, current_y_index, current_x_index, horizontal_velocity_at_k2, th, q)
    else
      call handle_stable_fluxes(current_state, current_y_index, current_x_index, horizontal_velocity_at_k2, th, q)
    end if
  end subroutine compute_using_fixed_surface_fluxes
  

  subroutine compute_using_fixed_surface_temperature(current_state, current_y_index, current_x_index, horizontal_velocity_at_k2, &
       zth, th, zq, q)
    type(model_state_type), target, intent(inout) :: current_state
    type(prognostic_field_type), intent(inout) :: th, zth
    type(prognostic_field_type), dimension(:), intent(inout) :: q, zq
    real(kind=DEFAULT_PRECISION), intent(in) :: horizontal_velocity_at_k2
    integer, intent(in) :: current_y_index, current_x_index

    real(kind=DEFAULT_PRECISION) :: dthv_surf, ustr, thvstr
    integer :: convergence_status, n

    if (current_state%passive_q) then ! i.e. q is not active 
      ! Assuming no liquid water at level 2
      dthv_surf = zth%data(2, current_y_index, current_x_index) + &
         current_state%global_grid%configuration%vertical%thref(2) - current_state%theta_virtual_surf
    else   
      dthv_surf=zth%data(2, current_y_index, current_x_index) + current_state%global_grid%configuration%vertical%thref(2)&
         *(1.0_DEFAULT_PRECISION + current_state%cq(current_state%water_vapour_mixing_ratio_index)*&
         zq(current_state%water_vapour_mixing_ratio_index)%data(2,current_y_index,current_x_index)) - &
         current_state%theta_virtual_surf                                                                        
    end if
    convergence_status = mostbc(current_state, horizontal_velocity_at_k2, dthv_surf,&
         current_state%global_grid%configuration%vertical%zn(2), ustr, thvstr)
                                                    
    current_state%vis_coefficient%data(1, current_y_index, current_x_index)=&
         current_state%global_grid%configuration%vertical%czn*ustr**2/horizontal_velocity_at_k2
    current_state%diff_coefficient%data(1, current_y_index, current_x_index)=&
         current_state%global_grid%configuration%vertical%czn*ustr*thvstr/(dthv_surf+smallp)
    zth%data(1, current_y_index, current_x_index)=2.0*current_state%theta_surf-zth%data(2, current_y_index, current_x_index)-&
         (current_state%global_grid%configuration%vertical%thref(2)+current_state%global_grid%configuration%vertical%thref(1))
    th%data(1, current_y_index, current_x_index)=2.0*current_state%theta_surf-th%data(2, current_y_index, current_x_index)-&
         (current_state%global_grid%configuration%vertical%thref(2)+current_state%global_grid%configuration%vertical%thref(1))
    
    if (current_state%number_q_fields .gt. 0) then
       n=iqv
       zq(n)%data(1, current_y_index, current_x_index)=zq(n)%data(2, current_y_index, current_x_index)
       q(n)%data(1, current_y_index, current_x_index)=q(n)%data(2, current_y_index, current_x_index)
       if (.not. current_state%passive_q) then
          zq(current_state%water_vapour_mixing_ratio_index)%data(1,current_y_index,current_x_index)=&
               2.0_DEFAULT_PRECISION*current_state%surface_vapour_mixing_ratio-&
               zq(current_state%water_vapour_mixing_ratio_index)%data(2,current_y_index,current_x_index)
          q(current_state%water_vapour_mixing_ratio_index)%data(1,current_y_index,current_x_index)=&
               2.0_DEFAULT_PRECISION*current_state%surface_vapour_mixing_ratio-&
               q(current_state%water_vapour_mixing_ratio_index)%data(2,current_y_index,current_x_index)
       endif
    end if    
  end subroutine compute_using_fixed_surface_temperature  

  subroutine simple_boundary_values(current_state, current_y_index, current_x_index, th, q)
    type(model_state_type), target, intent(inout) :: current_state
    type(prognostic_field_type), intent(inout) :: th
    type(prognostic_field_type), dimension(:), intent(inout) :: q
    integer, intent(in) :: current_y_index, current_x_index

    integer :: n

    current_state%vis_coefficient%data(1, current_y_index, current_x_index)=0.0_DEFAULT_PRECISION
    current_state%diff_coefficient%data(1, current_y_index, current_x_index)=0.0_DEFAULT_PRECISION
    if(current_state%backscatter) current_state%dis%data(1, current_y_index, current_x_index)=0.0_DEFAULT_PRECISION
    if(current_state%backscatter .and. current_state%th%active) &
         current_state%dis_th%data(1, current_y_index, current_x_index)=0.0_DEFAULT_PRECISION
    if (current_state%th%active) then
      th%data(1, current_y_index, current_x_index) = th%data(2, current_y_index, current_x_index)
    end if
    if (current_state%number_q_fields .gt. 0) then
      do n=1, current_state%number_q_fields
        q(n)%data(1, current_y_index, current_x_index) = q(n)%data(2, current_y_index, current_x_index)
        if (current_state%backscatter) current_state%disq(n)%data(1, current_y_index, current_x_index) = 0.0_DEFAULT_PRECISION
      end do
    end if
  end subroutine simple_boundary_values

  !> Solves the Monin-Obukhov equations in the case of specified surface values of temperature and mixing ratio,combined
  !! into a specified value of virtual temperature. It is a modified version of the subroutine described in Bull and
  !! Derbyshire (TDN197) based on the assumption that the similarity functions and roughness lengths for temperature and 
  !! mixing ratio are the same. In that case, all the original theory can be used if we replace temperature by virtual 
  !! temperature.
  !! The form of the non-dimensionalised wind shear used is 1.0 + BETAM z/L  for the stable case and
  !! (1.0 - GAMMAM z/L)**1/4 for the unstable case. The form of the non-dimensionalised temperature gradient used is
  !! 1.0 + BETAH z/L  for the stable case and (1.0 - GAMMAH z/L)**1/2 for the unstable case
  !! @param current_state The current model state
  !! @param delu The wind speed at the lowest grid point
  !! @param delt The virtual potential temperature difference between the lowest grid point and the surface
  !! @param z1 The height of the lowest grid point ABOVE the roughness length Z0
  !! @param ustrdg The diagnosed value of friction velocity
  !! @param tstrdg The diagnosed value of surface virtual temperature scale
  !! @returns The convergence criteria - success, richardson number too large or failure
  integer function mostbc(current_state, delu, delt, z1, ustrdg, tstrdg)
    type(model_state_type), target, intent(inout) :: current_state
    real(kind=DEFAULT_PRECISION), intent(in)  :: delu, delt, z1
    real(kind=DEFAULT_PRECISION), intent(out) :: ustrdg, tstrdg

    if (delt .lt. 0.0_DEFAULT_PRECISION) then
      if (delu .le. smallp) then 
        ustrdg=0.0_DEFAULT_PRECISION
        tstrdg=tstrcona*delt
      else
         ! The unstable case
         mostbc=solve_monin_obukhov_unstable_case(delu, delt, current_state%ellmocon, &
           ustrdg, tstrdg, current_state%global_grid%configuration%vertical)
      end if
    else if (delt .gt. 0.0_DEFAULT_PRECISION) then
      ! The stable case
      mostbc=solve_monin_obukhov_stable_case(delu, delt, current_state%global_grid%configuration%vertical%zlogm, &
           current_state%cmbc, ustrdg, tstrdg)
    else
      ! Trivial neutral case
        ustrdg=current_state%global_grid%configuration%vertical%vk_on_zlogm*delu
        tstrdg=0.0_DEFAULT_PRECISION 
        mostbc=CONVERGENCE_SUCCESS
    end if

    if (mostbc .ne. CONVERGENCE_SUCCESS) then
      if (mostbc .eq. CONVERGENCE_RICHARDSON_TOO_LARGE) then
        !call log_log(LOG_WARN, "Richardson number greater than critical value")
      else if(mostbc .eq. CONVERGENCE_FAILURE) then
        call log_log(LOG_ERROR, "Convergence failure after 200 iterations")
      end if
    end if   
  end function mostbc

  integer function solve_monin_obukhov_unstable_case(delu, delt, ellmocon, ustrdg, tstrdg, vertical_grid)
    real(kind=DEFAULT_PRECISION), intent(in) :: delu, delt, ellmocon
    real(kind=DEFAULT_PRECISION), intent(out) :: ustrdg, tstrdg
    type(vertical_grid_configuration_type), intent(inout) :: vertical_grid

    integer :: i
    real(kind=DEFAULT_PRECISION) :: ellmo, psim, psih, x4, xx, xx0, y2, yy, yy0, err_ustr, err_tstr, &
         ustrl, tstrl, & ! U and T star at start of iteration
         ustrpl, tstrpl ! U and T star at end of iteration
    
    ! First set initial values
    ustrl=vertical_grid%vk_on_zlogm*delu
    tstrl=tstrcona*delt

    ! Now start iteration
    do i=1, 200
      ellmo=ustrl*ustrl*ellmocon/tstrl

      ! Test for possible square root of negative quantity
      x4=1.0_DEFAULT_PRECISION-x4con/ellmo      
      if (x4 .lt. 0.0_DEFAULT_PRECISION) call log_log(LOG_ERROR, "Negative square root in x4")        

      xx=sqrt(sqrt(x4))
      xx0=sqrt(sqrt(1.0_DEFAULT_PRECISION-xx0con / ellmo))        
      psim=2.*( log((xx+1.0_DEFAULT_PRECISION)/(xx0+1.0_DEFAULT_PRECISION))-atan(xx)+atan(xx0) )+&
           log((xx*xx+1.0_DEFAULT_PRECISION)/(xx0*xx0+1.0_DEFAULT_PRECISION))
      ustrpl=von_karman_constant*delu/(vertical_grid%zlogm-psim)

      ! Test for possible square root of negative quantity
      y2=1.-y2con/ellmo
      if (y2 .lt. 0.0_DEFAULT_PRECISION) call log_log(LOG_ERROR, "Negative square root in y2")
      yy=sqrt(y2)
      yy0=sqrt(1.0_DEFAULT_PRECISION-yy0con/ellmo)
      psih=2.*log((1.0_DEFAULT_PRECISION+yy)/(1.0_DEFAULT_PRECISION+yy0))
      tstrpl=tstrconb*delt/(vertical_grid%zlogth-psih)
      err_ustr=abs((ustrpl-ustrl)/ ustrl)
      err_tstr=abs((tstrpl-tstrl)/ tstrl)                                       
      if ((err_tstr .lt. tolt) .and. (err_ustr .lt. tolm))  then                                                 
        ustrdg=ustrpl
        tstrdg=tstrpl
        solve_monin_obukhov_unstable_case=CONVERGENCE_SUCCESS
        return
      else                                                                 
        ustrl=(1.0_DEFAULT_PRECISION-smth)*ustrpl+smth*ustrl
        tstrl=(1.0_DEFAULT_PRECISION-smth)*tstrpl+smth*tstrl
      end if
    end do
    solve_monin_obukhov_unstable_case=CONVERGENCE_FAILURE
  end function solve_monin_obukhov_unstable_case

  integer function solve_monin_obukhov_stable_case(delu, delt, zlogm, cmbc, ustrdg, tstrdg)
    real(kind=DEFAULT_PRECISION), intent(in) :: delu, delt, zlogm, cmbc
    real(kind=DEFAULT_PRECISION), intent(out) :: ustrdg, tstrdg

    real(kind=DEFAULT_PRECISION) :: am, ah, ee, ff, det

    am=von_karman_constant*delu
    ah=von_karman_constant*delt
    ee=am*eecon
    ff=ah*cmbc-rhmbc*am*am  !
    det=ee*ee-DDBC_X4*ff
    solve_monin_obukhov_stable_case=CONVERGENCE_SUCCESS
    ! Test for laminar flow
    if (ff .gt. 0.0_DEFAULT_PRECISION) then
      if ((ee .lt. 0.0_DEFAULT_PRECISION).and.(det .gt. 0.0_DEFAULT_PRECISION)) then
        ustrdg=(-ee+sqrt(det))*r2ddbc
        tstrdg=ustrdg*(am-zlogm*ustrdg)*rcmbc
      else
        solve_monin_obukhov_stable_case=CONVERGENCE_RICHARDSON_TOO_LARGE        
        ustrdg=0.0_DEFAULT_PRECISION
        tstrdg=0.0_DEFAULT_PRECISION
      end if      
    else if (ddbc .eq. 0.0_DEFAULT_PRECISION) then
      ! Degenerate case
      ustrdg=-ff/ee
      tstrdg=delt*ustrdg/delu
    else
      ! Solve quadratic for USTRDG
      ustrdg=(-ee+sqrt(det))*r2ddbc
      tstrdg=ustrdg*(am-zlogm*ustrdg)*rcmbc
    end if
  end function solve_monin_obukhov_stable_case
end module lowerbc_mod
