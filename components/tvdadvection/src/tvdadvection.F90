!> Implements TVD advection for prognostic fields
module tvdadvection_mod
  use monc_component_mod, only : COMPONENT_ARRAY_FIELD_TYPE, COMPONENT_DOUBLE_DATA_TYPE, component_descriptor_type, &
       component_field_value_type, component_field_information_type
  use stencil_mod, only : grid_stencil_type, interpolate_to_dual, create_stencil, free_stencil
  use state_mod, only : model_state_type, parallel_state_type, FORWARD_STEPPING
  use grids_mod, only : local_grid_type, global_grid_type, X_INDEX, Y_INDEX, Z_INDEX
  use prognostics_mod, only : prognostic_field_type, prognostic_field_ptr_type
  use ultimateflux_mod, only : ultflx
  use datadefn_mod, only : DEFAULT_PRECISION, PRECISION_TYPE
  use optionsdatabase_mod, only : options_get_string
  use collections_mod, only : map_type
  use mpi, only : MPI_REQUEST_NULL, MPI_STATUS_IGNORE
  implicit none

#ifndef TEST_MODE
  private
#endif

  type(grid_stencil_type), save :: star_stencil
  integer, save :: u_index=0, v_index=0, w_index=0
  logical :: advect_flow, advect_th, advect_q
  real(kind=DEFAULT_PRECISION), dimension(:), allocatable :: flux_x, flux_y, flux_z, u_advection, v_advection, &
       w_advection, th_advection
  real(kind=DEFAULT_PRECISION), dimension(:,:), allocatable :: q_advection
  type(prognostic_field_type), dimension(:), allocatable :: interpolated_fields

  public tvdadvection_get_descriptor

contains

  !> Provides a description of this component for the core to register
  !! @returns The descriptor containing registration information for this component
  type(component_descriptor_type) function tvdadvection_get_descriptor()
    tvdadvection_get_descriptor%name="tvd_advection"
    tvdadvection_get_descriptor%version=0.1
    tvdadvection_get_descriptor%initialisation=>initialisation_callback
    tvdadvection_get_descriptor%timestep=>timestep_callback
    tvdadvection_get_descriptor%finalisation=>finalisation_callback

    tvdadvection_get_descriptor%field_value_retrieval=>field_value_retrieval_callback
    tvdadvection_get_descriptor%field_information_retrieval=>field_information_retrieval_callback
    allocate(tvdadvection_get_descriptor%published_fields(5))
    tvdadvection_get_descriptor%published_fields(1)="u_advection"
    tvdadvection_get_descriptor%published_fields(2)="v_advection"
    tvdadvection_get_descriptor%published_fields(3)="w_advection"
    tvdadvection_get_descriptor%published_fields(4)="th_advection"
    tvdadvection_get_descriptor%published_fields(5)="q_advection"
  end function tvdadvection_get_descriptor

  !> Field information retrieval callback, this returns information for a specific components published field
  !! @param current_state Current model state
  !! @param name The name of the field to retrieve information for
  !! @param field_information Populated with information about the field
  subroutine field_information_retrieval_callback(current_state, name, field_information)
    type(model_state_type), target, intent(inout) :: current_state
    character(len=*), intent(in) :: name
    type(component_field_information_type), intent(out) :: field_information

    ! Field description is the same regardless of the specific field being retrieved
    field_information%field_type=COMPONENT_ARRAY_FIELD_TYPE
    field_information%data_type=COMPONENT_DOUBLE_DATA_TYPE
    if (name .eq. "q_advection") then
      field_information%number_dimensions=2
    else
      field_information%number_dimensions=1
    end if
    field_information%dimension_sizes(1)=current_state%local_grid%size(Z_INDEX)
    if (name .eq. "q_advection") field_information%dimension_sizes(2)=current_state%number_q_fields
    field_information%enabled=.true.
  end subroutine field_information_retrieval_callback

  !> Field value retrieval callback, this returns the value of a specific published field
  !! @param current_state Current model state
  !! @param name The name of the field to retrieve the value for
  !! @param field_value Populated with the value of the field
  subroutine field_value_retrieval_callback(current_state, name, field_value)
    type(model_state_type), target, intent(inout) :: current_state
    character(len=*), intent(in) :: name
    type(component_field_value_type), intent(out) :: field_value
    
    if (name .eq. "u_advection") then
      allocate(field_value%real_1d_array(size(u_advection)), source=u_advection)
    else if (name .eq. "v_advection") then
      allocate(field_value%real_1d_array(size(v_advection)), source=v_advection)
    else if (name .eq. "w_advection") then
      allocate(field_value%real_1d_array(size(w_advection)), source=w_advection)
    else if (name .eq. "th_advection") then
      allocate(field_value%real_1d_array(size(th_advection)), source=th_advection)
    else if (name .eq. "q_advection") then
      allocate(field_value%real_2d_array(size(q_advection, 1), size(q_advection, 2)), source=q_advection)
    end if
  end subroutine field_value_retrieval_callback

  !> Sets up the stencil_mod (used in interpolation) and allocates data for the flux fields
  !! @param current_state The current model state_mod
  subroutine initialisation_callback(current_state)
    type(model_state_type), target, intent(inout) :: current_state

    type(prognostic_field_ptr_type), dimension(3) :: fields
    integer, dimension(3, 2) :: sizes
    integer :: num_fields
    logical :: xdim, ydim

    xdim=.false.
    ydim=.false.
    num_fields=0
#ifdef U_ACTIVE    
    xdim=.true.
    num_fields = num_fields + 1
    fields(num_fields)%ptr => current_state%u
    sizes(num_fields,:) = (/ 2, 2 /) ! need um2 therefore -2 (applies to all interpolations)
    u_index = num_fields    
#endif

#ifdef V_ACTIVE
    ydim=.true.
    num_fields = num_fields + 1     
    fields(num_fields)%ptr => current_state%v
    sizes(num_fields,:) = (/ 1, 1 /)
    v_index=num_fields
#endif

#ifdef W_ACTIVE    
    num_fields = num_fields + 1
    fields(num_fields)%ptr => current_state%w
    sizes(num_fields,:) = (/ 1, 1 /)
    w_index=num_fields
#endif
    ! Allocate from 0, as any inactive dimensions will issue 0 to the ultimate flux which ignores the field
    allocate(interpolated_fields(0:num_fields))
#ifdef U_ACTIVE
    allocate(interpolated_fields(u_index)%data(current_state%global_grid%size(Z_INDEX), -1:3, -1:3))
    interpolated_fields(u_index)%active=.true.
#endif
#ifdef V_ACTIVE
    allocate(interpolated_fields(v_index)%data(current_state%global_grid%size(Z_INDEX), 0:2, 0:2))
    interpolated_fields(v_index)%active=.true.
#endif
#ifdef W_ACTIVE
    allocate(interpolated_fields(w_index)%data(current_state%global_grid%size(Z_INDEX), 0:2, 0:2))
    interpolated_fields(w_index)%active=.true.
#endif

    star_stencil = create_stencil(current_state%local_grid, fields, num_fields, 3, sizes, xdim, ydim)
    allocate(flux_y(current_state%global_grid%size(Z_INDEX)))
    allocate(flux_z(current_state%global_grid%size(Z_INDEX)))
    allocate(flux_x(current_state%global_grid%size(Z_INDEX)))
    allocate(u_advection(current_state%global_grid%size(Z_INDEX)), v_advection(current_state%global_grid%size(Z_INDEX)), &
         w_advection(current_state%global_grid%size(Z_INDEX)), th_advection(current_state%global_grid%size(Z_INDEX)), &
         q_advection(current_state%global_grid%size(Z_INDEX), current_state%number_q_fields))
    advect_flow=determine_if_advection_here(options_get_string(current_state%options_database, "advection_flow_fields"))    
    advect_th=determine_if_advection_here(options_get_string(current_state%options_database, "advection_theta_field"))
    advect_q=determine_if_advection_here(options_get_string(current_state%options_database, "advection_q_fields"))   
  end subroutine initialisation_callback

  !> Frees up the memory associated with the advection
  !! @param current_state The current model state
  subroutine finalisation_callback(current_state)
    type(model_state_type), target, intent(inout) :: current_state

    call free_stencil(star_stencil)
    if (allocated(flux_x)) deallocate(flux_x)
    if (allocated(flux_y)) deallocate(flux_y)
    if (allocated(flux_z)) deallocate(flux_z)
    if (allocated(interpolated_fields)) deallocate(interpolated_fields)
    if (allocated(u_advection)) deallocate(u_advection)
    if (allocated(v_advection)) deallocate(v_advection)
    if (allocated(w_advection)) deallocate(w_advection)
    if (allocated(th_advection)) deallocate(th_advection)
    if (allocated(q_advection)) deallocate(q_advection)
  end subroutine finalisation_callback  

  !> Timestep callback hook which performs the TVD advection for each prognostic field
  !! @param current_state The current model state_mod
  subroutine timestep_callback(current_state)
    type(model_state_type), target, intent(inout) :: current_state

    if (current_state%halo_column) then
      if (.not. ((current_state%column_local_y == current_state%local_grid%halo_size(Y_INDEX) .and. &
           current_state%column_local_x .le. current_state%local_grid%local_domain_end_index(X_INDEX) .and. &
           current_state%column_local_x .ge. current_state%local_grid%local_domain_start_index(X_INDEX)-1) .or. &
           (current_state%column_local_x == current_state%local_grid%halo_size(X_INDEX) .and. &
           current_state%column_local_y .ge. current_state%local_grid%local_domain_start_index(Y_INDEX) &
           .and. current_state%column_local_y .le. current_state%local_grid%local_domain_end_index(Y_INDEX)) )) return
    end if

    if (advect_flow) call advect_flow_fields(current_state)
    if (advect_th) call advect_theta(current_state)
    if (advect_q) call advect_q_fields(current_state)
  end subroutine timestep_callback

  !> Will advect the flow fields
  !! @param current_state The current model state_mod
  subroutine advect_flow_fields(current_state)
    type(model_state_type), intent(inout) :: current_state

    real(kind=DEFAULT_PRECISION) :: dtm

    dtm = current_state%dtm*2.0_DEFAULT_PRECISION
    if (current_state%momentum_stepping == FORWARD_STEPPING) dtm=current_state%dtm

#ifdef U_ACTIVE
    call advect_u(current_state%column_local_y, current_state%column_local_x, dtm, current_state%u, &
         current_state%v, current_state%w, current_state%zu, current_state%su, current_state%global_grid, &
         current_state%local_grid, current_state%parallel, current_state%halo_column, current_state%field_stepping)
#endif

#ifdef V_ACTIVE
    call advect_v(current_state%column_local_y, current_state%column_local_x, dtm, current_state%u, &
         current_state%v, current_state%w, current_state%zv, current_state%sv, current_state%global_grid, &
         current_state%local_grid, current_state%parallel, current_state%halo_column, current_state%field_stepping)
#endif

#ifdef W_ACTIVE
    call advect_w(current_state%column_local_y, current_state%column_local_x, dtm, current_state%u, &
         current_state%v, current_state%w, current_state%zw, current_state%sw, current_state%global_grid,&
         current_state%local_grid, current_state%parallel, current_state%halo_column, current_state%field_stepping)
#endif
  end subroutine advect_flow_fields

  !> Advects the Q fields
  !! @param current_state The current model state_mod
  subroutine advect_q_fields(current_state)
    type(model_state_type), intent(inout) :: current_state

    integer :: i
    real(kind=DEFAULT_PRECISION) :: dtm

    dtm = current_state%dtm*2.0_DEFAULT_PRECISION
    if (current_state%scalar_stepping == FORWARD_STEPPING) dtm=current_state%dtm

    do i=1,current_state%number_q_fields
      if (current_state%q(i)%active) then           
        call advect_scalar_field(current_state%column_local_y, current_state%column_local_x, dtm, current_state%u, &
             current_state%v, current_state%w, current_state%zq(i), current_state%q(i), current_state%sq(i), &
             current_state%global_grid, current_state%local_grid, current_state%parallel, &
             current_state%halo_column, current_state%field_stepping)
        q_advection(:,i)=current_state%sq(i)%data(:, current_state%column_local_y, current_state%column_local_x)
      end if
    end do
  end subroutine advect_q_fields

  !> Advects the theta field if it is active
  !! @param current_state The current model state_mod
  subroutine advect_theta(current_state)
    type(model_state_type), intent(inout) :: current_state

    real(kind=DEFAULT_PRECISION) :: dtm

    dtm = current_state%dtm*2.0_DEFAULT_PRECISION
    if (current_state%scalar_stepping == FORWARD_STEPPING) dtm=current_state%dtm

    if (current_state%th%active) then
      call advect_scalar_field(current_state%column_local_y, current_state%column_local_x, dtm, current_state%u,&
           current_state%v, current_state%w, current_state%zth, current_state%th, current_state%sth, current_state%global_grid,&
           current_state%local_grid, current_state%parallel, current_state%halo_column, current_state%field_stepping)
      th_advection=current_state%sth%data(:, current_state%column_local_y, current_state%column_local_x)
    end if
  end subroutine advect_theta

  !> Advects a single scalar field
  subroutine advect_scalar_field(y_local_index, x_local_index, dt, u, v, w, z_q_field, q_field, source_field, &
       global_grid, local_grid, parallel, halo_column, field_stepping)

    integer, intent(in) ::y_local_index, x_local_index, field_stepping
    real(kind=DEFAULT_PRECISION), intent(in) ::dt   ! timestep (s)
    logical, intent(in) :: halo_column
    type(prognostic_field_type), intent(inout) :: u, w, v, z_q_field, q_field, source_field
    type(global_grid_type), intent(inout) :: global_grid
    type(local_grid_type), intent(inout) :: local_grid
    type(parallel_state_type), intent(inout) :: parallel

    if (.not. allocated(q_field%flux_previous_x)) allocate(q_field%flux_previous_x(local_grid%size(Z_INDEX), &
         local_grid%size(Y_INDEX)+4))
    if (.not. allocated(q_field%flux_previous_y)) allocate(q_field%flux_previous_y(local_grid%size(Z_INDEX)))

    call register_y_flux_wrap_recv_if_required(y_local_index, q_field, parallel, local_grid)

    if (field_stepping == FORWARD_STEPPING) then
      call ultflx(y_local_index, x_local_index, u, v, w, y_local_index, x_local_index, q_field, local_grid, &
           global_grid%configuration, parallel, 0, dt, &
           flux_y, flux_z, flux_x, q_field%flux_previous_x(:,y_local_index), global_grid%configuration%vertical%rdz,&
           global_grid%configuration%vertical%rdzn, global_grid%configuration%vertical%dzn, 2, local_grid%size(Z_INDEX))
    else
      call ultflx(y_local_index, x_local_index, u, v, w, y_local_index, x_local_index, z_q_field, local_grid, &
           global_grid%configuration, parallel, 0, dt, &
           flux_y, flux_z, flux_x, q_field%flux_previous_x(:,y_local_index), global_grid%configuration%vertical%rdz,&
           global_grid%configuration%vertical%rdzn, global_grid%configuration%vertical%dzn, 2, local_grid%size(Z_INDEX))
    end if

    call complete_y_flux_wrap_recv_if_required(y_local_index, q_field, parallel, local_grid)

    if (.not. halo_column) then
      call differentiate_face_values(y_local_index, x_local_index, u, v, w, y_local_index, x_local_index, source_field, &
           local_grid, global_grid, q_field%flux_previous_y, q_field%flux_previous_x(:,y_local_index), &
           global_grid%configuration%vertical%tzc1, global_grid%configuration%vertical%tzc2, .true.)
    end if

    if (y_local_index .lt. local_grid%local_domain_end_index(Y_INDEX)) q_field%flux_previous_y(:) = flux_y(:)
    call register_y_flux_wrap_send_if_required(y_local_index, q_field, parallel, local_grid)
    call complete_y_flux_wrap_send_if_required(y_local_index, q_field, parallel, local_grid)
  end subroutine advect_scalar_field

  !> Advects the U flow field
  subroutine advect_u(y_local_index, x_local_index, dt, u, v, w, zf, su, global_grid, local_grid, parallel, &
       halo_column, field_stepping)

    integer, intent(in) :: y_local_index, x_local_index, field_stepping
    real(kind=DEFAULT_PRECISION), intent(in) ::dt   ! timestep (s)
    logical, intent(in) :: halo_column
    type(prognostic_field_type), intent(inout) :: u, w, v, zf, su
    type(global_grid_type), intent(inout) :: global_grid
    type(local_grid_type), intent(inout) :: local_grid
    type(parallel_state_type), intent(inout) :: parallel

    if (.not. allocated(u%flux_previous_x)) allocate(u%flux_previous_x(local_grid%size(Z_INDEX), local_grid%size(Y_INDEX)+4))
    if (.not. allocated(u%flux_previous_y)) allocate(u%flux_previous_y(local_grid%size(Z_INDEX)))

    call register_y_flux_wrap_recv_if_required(y_local_index, u, parallel, local_grid)

    call interpolate_to_dual(local_grid, u, star_stencil, x_local_index, y_local_index, interpolated_fields, u_index)

    if (field_stepping == FORWARD_STEPPING) then
      call ultflx(1, 1, interpolated_fields(u_index), interpolated_fields(v_index), &
           interpolated_fields(w_index), y_local_index, x_local_index, u, local_grid, global_grid%configuration, parallel, 0, &
           dt, flux_y, flux_z, flux_x, &
           u%flux_previous_x(:,y_local_index), global_grid%configuration%vertical%rdz, &
           global_grid%configuration%vertical%rdzn, global_grid%configuration%vertical%dzn, 2, local_grid%size(Z_INDEX))
    else
      ! Centred stepping
      call ultflx(1, 1, interpolated_fields(u_index), interpolated_fields(v_index), &
           interpolated_fields(w_index), y_local_index, x_local_index, zf, local_grid, global_grid%configuration, parallel, 0, &
           dt, flux_y, flux_z, flux_x, &
           u%flux_previous_x(:,y_local_index), global_grid%configuration%vertical%rdz, &
           global_grid%configuration%vertical%rdzn, global_grid%configuration%vertical%dzn, 2, local_grid%size(Z_INDEX))
    end if

    call complete_y_flux_wrap_recv_if_required(y_local_index, u, parallel, local_grid)

    if (.not. halo_column) then
      call differentiate_face_values(1, 1, interpolated_fields(u_index), interpolated_fields(v_index), &
           interpolated_fields(w_index), y_local_index, x_local_index, su, local_grid, global_grid, &
           u%flux_previous_y, u%flux_previous_x(:,y_local_index), &
           global_grid%configuration%vertical%tzc1, global_grid%configuration%vertical%tzc2, .true.)
      u_advection=su%data(:, y_local_index, x_local_index)
    end if

    if (y_local_index .lt. local_grid%local_domain_end_index(Y_INDEX)) u%flux_previous_y(:) = flux_y(:)
    call register_y_flux_wrap_send_if_required(y_local_index, u, parallel, local_grid)
    call complete_y_flux_wrap_send_if_required(y_local_index, u, parallel, local_grid)
  end subroutine advect_u

  !> Advects the V flow field
  subroutine advect_v(y_local_index, x_local_index, dt, u, v, w, zf, sv, global_grid, local_grid, parallel, &
       halo_column, field_stepping)

    integer, intent(in) ::y_local_index, x_local_index, field_stepping
    real(kind=DEFAULT_PRECISION), intent(in) ::dt   ! timestep (s)
    logical, intent(in) :: halo_column
    type(prognostic_field_type), intent(inout) :: u, w, v, zf, sv
    type(global_grid_type), intent(inout) :: global_grid
    type(local_grid_type), intent(inout) :: local_grid
    type(parallel_state_type), intent(inout) :: parallel

    if (.not. allocated(v%flux_previous_x)) allocate(v%flux_previous_x(local_grid%size(Z_INDEX), local_grid%size(Y_INDEX)+4))
    if (.not. allocated(v%flux_previous_y)) allocate(v%flux_previous_y(local_grid%size(Z_INDEX)))

    call register_y_flux_wrap_recv_if_required(y_local_index, v, parallel, local_grid)

    call interpolate_to_dual(local_grid, v, star_stencil, x_local_index, y_local_index, interpolated_fields, v_index)

    if (field_stepping == FORWARD_STEPPING) then
      call ultflx(1, 1, interpolated_fields(u_index), interpolated_fields(v_index), &
           interpolated_fields(w_index), y_local_index, x_local_index, v, local_grid, global_grid%configuration, parallel, 0, &
           dt, flux_y, flux_z, flux_x, &
           v%flux_previous_x(:,y_local_index), global_grid%configuration%vertical%rdz, &
           global_grid%configuration%vertical%rdzn, global_grid%configuration%vertical%dzn, 2, local_grid%size(Z_INDEX))
    else
      ! Centred stepping
      call ultflx(1, 1, interpolated_fields(u_index), interpolated_fields(v_index), &
           interpolated_fields(w_index), y_local_index, x_local_index, zf, local_grid, global_grid%configuration, parallel, 0, &
           dt, flux_y, flux_z, flux_x, &
           v%flux_previous_x(:,y_local_index), global_grid%configuration%vertical%rdz, &
           global_grid%configuration%vertical%rdzn, global_grid%configuration%vertical%dzn, 2, local_grid%size(Z_INDEX))
    end if

    call complete_y_flux_wrap_recv_if_required(y_local_index, v, parallel, local_grid)

    if (.not. halo_column) then
      call differentiate_face_values(1, 1, interpolated_fields(u_index), interpolated_fields(v_index), &
           interpolated_fields(w_index), y_local_index, x_local_index, sv, local_grid, global_grid, &
           v%flux_previous_y, v%flux_previous_x(:,y_local_index), &           
           global_grid%configuration%vertical%tzc1, global_grid%configuration%vertical%tzc2, .true.)
      v_advection=sv%data(:, y_local_index, x_local_index)
    end if

    if (y_local_index .lt. local_grid%local_domain_end_index(Y_INDEX)) v%flux_previous_y(:) = flux_y(:)
    call register_y_flux_wrap_send_if_required(y_local_index, v, parallel, local_grid)
    call complete_y_flux_wrap_send_if_required(y_local_index, v, parallel, local_grid)
  end subroutine advect_v

  !> Advects the W flow field
  subroutine advect_w(y_local_index, x_local_index, dt, u, v, w, zf, sw, global_grid, local_grid, parallel, &
       halo_column, field_stepping)

    integer, intent(in) ::y_local_index, x_local_index, field_stepping
    real(kind=DEFAULT_PRECISION), intent(in) ::dt   ! timestep (s)
    logical, intent(in) :: halo_column
    type(prognostic_field_type), intent(inout) :: u, w, v, zf, sw
    type(global_grid_type), intent(inout) :: global_grid
    type(local_grid_type), intent(inout) :: local_grid
    type(parallel_state_type), intent(inout) :: parallel

    if (.not. allocated(w%flux_previous_x)) allocate(w%flux_previous_x(local_grid%size(Z_INDEX), local_grid%size(Y_INDEX)+4))
    if (.not. allocated(w%flux_previous_y)) allocate(w%flux_previous_y(local_grid%size(Z_INDEX)))

    call register_y_flux_wrap_recv_if_required(y_local_index, w, parallel, local_grid)

    call interpolate_to_dual(local_grid, w, star_stencil, x_local_index, y_local_index, interpolated_fields, w_index)

    if (field_stepping == FORWARD_STEPPING) then
      call ultflx(1, 1, interpolated_fields(u_index), interpolated_fields(v_index), &
           interpolated_fields(w_index), y_local_index, x_local_index, w, local_grid, global_grid%configuration, parallel, 1, &
           dt, flux_y, flux_z, flux_x,&
           w%flux_previous_x(:,y_local_index), global_grid%configuration%vertical%rdzn, &
           global_grid%configuration%vertical%rdz, global_grid%configuration%vertical%dz, 1, local_grid%size(Z_INDEX)-1)
    else
      ! Centred stepping
      call ultflx(1, 1, interpolated_fields(u_index), interpolated_fields(v_index), &
           interpolated_fields(w_index), y_local_index, x_local_index, zf, local_grid, global_grid%configuration, parallel, 1, &
           dt, flux_y, flux_z, flux_x,&
           w%flux_previous_x(:,y_local_index), global_grid%configuration%vertical%rdzn, &
           global_grid%configuration%vertical%rdz, global_grid%configuration%vertical%dz, 1, local_grid%size(Z_INDEX)-1)
    end if

    call complete_y_flux_wrap_recv_if_required(y_local_index, w, parallel, local_grid)

    if (.not. halo_column) then
      call differentiate_face_values(1, 1, interpolated_fields(u_index), interpolated_fields(v_index),&
           interpolated_fields(w_index), y_local_index, x_local_index, sw, local_grid, global_grid, &
           w%flux_previous_y, w%flux_previous_x(:,y_local_index),&
           global_grid%configuration%vertical%tzd1, global_grid%configuration%vertical%tzd2, .false.)
      w_advection=sw%data(:, y_local_index, x_local_index)
    end if
    if (y_local_index .lt. local_grid%local_domain_end_index(Y_INDEX)) w%flux_previous_y(:) = flux_y(:)
    call register_y_flux_wrap_send_if_required(y_local_index, w, parallel, local_grid)
    call complete_y_flux_wrap_send_if_required(y_local_index, w, parallel, local_grid)
  end subroutine advect_w

  !> Differentiates face values to update the source field
  subroutine differentiate_face_values(y_flow_index, x_flow_index, u, v, w, y_source_index, x_source_index, source_field, &
       local_grid, global_grid, flux_y_previous, flux_x_previous, tzc1, tzc2, differentiate_top)

    integer, intent(in) :: y_flow_index, x_flow_index, y_source_index, x_source_index
    logical, intent(in) :: differentiate_top
    real(kind=DEFAULT_PRECISION), intent(in), dimension(*) :: tzc1, tzc2
    type(prognostic_field_type), intent(inout) :: u, w, v
    type(prognostic_field_type), intent(inout) :: source_field
    type(global_grid_type), intent(inout) :: global_grid
    type(local_grid_type), intent(inout) :: local_grid
    real(kind=DEFAULT_PRECISION), intent(in), dimension(:) :: flux_y_previous, flux_x_previous

    integer :: k

    do k=2,local_grid%size(Z_INDEX)-1
      source_field%data(k, y_source_index, x_source_index)=0.0_DEFAULT_PRECISION
#ifdef V_ACTIVE
      source_field%data(k, y_source_index, x_source_index)=source_field%data(k, y_source_index, x_source_index)+&
           (v%data(k, y_flow_index-1, x_flow_index)* flux_y_previous(k) - v%data(k, y_flow_index, x_flow_index)*flux_y(k))*&
           global_grid%configuration%horizontal%cy
#endif
#ifdef W_ACTIVE
      source_field%data(k, y_source_index, x_source_index)=source_field%data(k, y_source_index, x_source_index)+&
           4.0_DEFAULT_PRECISION*(w%data(k-1, y_flow_index, x_flow_index)* flux_z(k)*tzc1(k) - &
           w%data(k, y_flow_index, x_flow_index)*flux_z(k+1)*tzc2(k))
#endif
#ifdef U_ACTIVE
      source_field%data(k, y_source_index, x_source_index)=source_field%data(k, y_source_index, x_source_index)+&
           (u%data(k, y_flow_index, x_flow_index-1)* flux_x(k) - u%data(k, y_flow_index, x_flow_index)*flux_x_previous(k))*&
           global_grid%configuration%horizontal%cx
#endif
    end do
    if (differentiate_top) then
      k=local_grid%size(Z_INDEX)
      source_field%data(k, y_source_index, x_source_index)=0.0_DEFAULT_PRECISION
#ifdef V_ACTIVE
      source_field%data(k, y_source_index, x_source_index)=source_field%data(k, y_source_index, x_source_index)+&
           (v%data(k, y_flow_index-1, x_flow_index)* flux_y_previous(k) - v%data(k, y_flow_index, x_flow_index)*flux_y(k))*&
           global_grid%configuration%horizontal%cy
#endif
#ifdef W_ACTIVE
      source_field%data(k, y_source_index, x_source_index)=source_field%data(k, y_source_index, x_source_index)+&
           4.0_DEFAULT_PRECISION*tzc1(k)* w%data(k-1, y_flow_index, x_flow_index)*flux_z(k)
#endif
#ifdef U_ACTIVE
      source_field%data(k, y_source_index, x_source_index)=source_field%data(k, y_source_index, x_source_index)+&
           (u%data(k, y_flow_index, x_flow_index-1)* flux_x(k) -u%data(k, y_flow_index, x_flow_index)*flux_x_previous(k))*&
           global_grid%configuration%horizontal%cx
#endif
    end if
  end subroutine differentiate_face_values

  !> Completes the Y flux MPI asynchronous send if required
  !! @param y_local_index The local index in Y
  !! @param field The prognostic field
  !! @param parallel Parallel system description
  !! @param local_grid The local grid description
  subroutine complete_y_flux_wrap_send_if_required(y_local_index, field, parallel, local_grid)
    integer, intent(in) :: y_local_index
    type(prognostic_field_type), intent(inout) :: field
    type(parallel_state_type), intent(inout) :: parallel
    type(local_grid_type), intent(inout) :: local_grid

    integer :: ierr

    if (y_local_index == local_grid%local_domain_end_index(Y_INDEX) .and. parallel%my_coords(Y_INDEX) == 0 .and. &
         field%async_flux_handle .ne. MPI_REQUEST_NULL) then
      call mpi_wait(field%async_flux_handle, MPI_STATUS_IGNORE, ierr)
    end if
  end subroutine complete_y_flux_wrap_send_if_required

  !> Registers an asynchronous send for the Y flux if required.
  !!
  !! This is done after the second y is computed and we have until the entire Y dimension is completed
  !! until the communication must be complete. If the wrap around process is the same (one process in Y dimension)
  !! then just issues a local copy to the buffer.
  !! @param y_local_index The local index in Y
  !! @param field The prognostic field
  !! @param parallel Parallel system description
  !! @param local_grid The local grid description
  subroutine register_y_flux_wrap_send_if_required(y_local_index, field, parallel, local_grid)
    integer, intent(in) :: y_local_index
    type(prognostic_field_type), intent(inout) :: field
    type(parallel_state_type), intent(inout) :: parallel
    type(local_grid_type), intent(inout) :: local_grid

    integer :: ierr

    if (y_local_index == local_grid%local_domain_start_index(Y_INDEX)-1 .and. parallel%my_coords(Y_INDEX) == 0) then
      if (.not. allocated(field%flux_y_buffer)) allocate(field%flux_y_buffer(local_grid%size(Z_INDEX)))
      field%flux_y_buffer(:) = flux_y(:)
      if (parallel%my_rank .ne. local_grid%neighbours(Y_INDEX,1)) then      
        call mpi_isend(field%flux_y_buffer, local_grid%size(Z_INDEX), PRECISION_TYPE, local_grid%neighbours(Y_INDEX,1), 0, &
             parallel%neighbour_comm, field%async_flux_handle, ierr)
      end if
    end if
  end subroutine register_y_flux_wrap_send_if_required

  !> Completes the Y flux MPI asynchronous recieve if required. If the wrap around process is the same (one process
  !! in the y dimension) then just issues a local copy
  !! @param y_local_index The local index in Y
  !! @param field The prognostic field
  !! @param parallel Parallel system description
  !! @param local_grid The local grid description
  subroutine complete_y_flux_wrap_recv_if_required(y_local_index, field, parallel, local_grid)
    integer, intent(in) :: y_local_index
    type(prognostic_field_type), intent(inout) :: field
    type(parallel_state_type), intent(inout) :: parallel
    type(local_grid_type), intent(inout) :: local_grid

    integer :: ierr

    if (y_local_index == local_grid%local_domain_end_index(Y_INDEX) .and. parallel%my_coords(Y_INDEX) == &
         parallel%dim_sizes(Y_INDEX)-1) then
      if (field%async_flux_handle .ne. MPI_REQUEST_NULL) then
        call mpi_wait(field%async_flux_handle, MPI_STATUS_IGNORE, ierr)
      end if
      flux_y(:) = field%flux_y_buffer(:)
    end if
  end subroutine complete_y_flux_wrap_recv_if_required

  !> Registers an MPI asynchronous receive for the flux if required.
  !!
  !! This is registered at the start and we have until the last column in Y until it must be completed. No 
  !! communication is registered if this is a local operation
  !! @param y_local_index The local index in Y
  !! @param field The prognostic field
  !! @param parallel Parallel system description
  !! @param local_grid The local grid description
  subroutine register_y_flux_wrap_recv_if_required(y_local_index, field, parallel, local_grid)
    integer, intent(in) :: y_local_index
    type(prognostic_field_type), intent(inout) :: field
    type(parallel_state_type), intent(inout) :: parallel
    type(local_grid_type), intent(inout) :: local_grid

    integer :: ierr

    if (y_local_index == local_grid%local_domain_start_index(Y_INDEX) .and. parallel%my_coords(Y_INDEX) == &
         parallel%dim_sizes(Y_INDEX)-1) then
      if (parallel%my_rank .ne. local_grid%neighbours(Y_INDEX,3)) then
        if (.not. allocated(field%flux_y_buffer)) allocate(field%flux_y_buffer(local_grid%size(Z_INDEX)))
        call mpi_irecv(field%flux_y_buffer, local_grid%size(Z_INDEX), PRECISION_TYPE, local_grid%neighbours(Y_INDEX,3), 0, &
             parallel%neighbour_comm, field%async_flux_handle, ierr)
      end if
    end if
  end subroutine register_y_flux_wrap_recv_if_required

  !> Parses a field string (read in from the configuration file) and determines whether this algorithm should be used
  !! for advecting that field
  !! @param field The string configuration of field advection
  !! @returns Whether or not the field is advected here
  logical function determine_if_advection_here(field)
    character(len=*), intent(in) :: field

    if (len_trim(field) .ne. 0) then
      if (trim(field) .eq. "tvd" .or. trim(field) .eq. "any") then
        determine_if_advection_here=.true.
      else
        determine_if_advection_here=.false.
      end if
    else
      determine_if_advection_here=.true.
    end if
  end function determine_if_advection_here
end module tvdadvection_mod
