!> Computes the viscosity dynamics for the U,V,W source terms
module viscosity_mod
  use datadefn_mod, only : DEFAULT_PRECISION
  use monc_component_mod, only : COMPONENT_ARRAY_FIELD_TYPE, &
      COMPONENT_DOUBLE_DATA_TYPE, component_descriptor_type, &
      component_field_value_type, component_field_information_type
  use optionsdatabase_mod, only : options_get_integer
  use state_mod, only : FORWARD_STEPPING, model_state_type
  use prognostics_mod, only : prognostic_field_type
  use grids_mod, only : Z_INDEX, X_INDEX, Y_INDEX
  use communication_types_mod, only : halo_communication_type, neighbour_description_type, field_data_wrapper_type
  use halo_communication_mod, only : copy_buffer_to_field, perform_local_data_copy_for_field, complete_nonblocking_halo_swap, &
       copy_buffer_to_corner
  use registry_mod, only : is_component_enabled
  use logging_mod, only : LOG_ERROR, LOG_WARN, log_master_log

  implicit none

#ifndef TEST_MODE
  private
#endif

  real(kind=DEFAULT_PRECISION), dimension(:), allocatable :: u_viscosity, v_viscosity, w_viscosity

  ! Local tendency diagnostic variables for this component
  ! 3D tendency fields and logicals for their use
  real(kind=DEFAULT_PRECISION), dimension(:,:,:), allocatable ::     &
       tend_3d_u, tend_3d_v, tend_3d_w
  logical :: l_tend_3d_u, l_tend_3d_v, l_tend_3d_w
  ! Local mean tendency profile fields and logicals for their use
  real(kind=DEFAULT_PRECISION), dimension(:), allocatable ::         &
       tend_pr_tot_u, tend_pr_tot_v, tend_pr_tot_w
  logical :: l_tend_pr_tot_u, l_tend_pr_tot_v, l_tend_pr_tot_w

  integer :: diagnostic_generation_frequency

  public viscosity_get_descriptor

contains

  !> Provides the descriptor back to the caller and is used in component registration
  !! @returns The termination check component descriptor
  type(component_descriptor_type) function viscosity_get_descriptor()
    viscosity_get_descriptor%name="viscosity"
    viscosity_get_descriptor%version=0.1
    viscosity_get_descriptor%initialisation=>initialisation_callback
    viscosity_get_descriptor%timestep=>timestep_callback
    viscosity_get_descriptor%finalisation=>finalisation_callback

    viscosity_get_descriptor%field_value_retrieval=>field_value_retrieval_callback
    viscosity_get_descriptor%field_information_retrieval=>field_information_retrieval_callback
    allocate(viscosity_get_descriptor%published_fields(3+3+3))

    viscosity_get_descriptor%published_fields(1)="u_viscosity"
    viscosity_get_descriptor%published_fields(2)="v_viscosity"
    viscosity_get_descriptor%published_fields(3)="w_viscosity"

    viscosity_get_descriptor%published_fields(3+1)= "tend_u_viscosity_3d_local"
    viscosity_get_descriptor%published_fields(3+2)= "tend_v_viscosity_3d_local"
    viscosity_get_descriptor%published_fields(3+3)= "tend_w_viscosity_3d_local"

    viscosity_get_descriptor%published_fields(3+3+1)= "tend_u_viscosity_profile_total_local"
    viscosity_get_descriptor%published_fields(3+3+2)= "tend_v_viscosity_profile_total_local"
    viscosity_get_descriptor%published_fields(3+3+3)= "tend_w_viscosity_profile_total_local"

  end function viscosity_get_descriptor

  !> Field information retrieval callback, this returns information for a specific components published field
  !! @param current_state Current model state
  !! @param name The name of the field to retrieve information for
  !! @param field_information Populated with information about the field
  subroutine field_information_retrieval_callback(current_state, name, field_information)
    type(model_state_type), target, intent(inout) :: current_state
    character(len=*), intent(in) :: name
    type(component_field_information_type), intent(out) :: field_information
    integer :: strcomp

    ! Field description is the same regardless of the specific field being retrieved
    field_information%field_type=COMPONENT_ARRAY_FIELD_TYPE
    field_information%data_type=COMPONENT_DOUBLE_DATA_TYPE
    field_information%number_dimensions=1
    field_information%dimension_sizes(1)=current_state%local_grid%size(Z_INDEX)
    field_information%enabled=.true.

    ! Field information for 3d
    strcomp=INDEX(name, "viscosity_3d_local")
    if (strcomp .ne. 0) then
      field_information%field_type=COMPONENT_ARRAY_FIELD_TYPE
      field_information%number_dimensions=3
      field_information%dimension_sizes(1)=current_state%local_grid%size(Z_INDEX)
      field_information%dimension_sizes(2)=current_state%local_grid%size(Y_INDEX)
      field_information%dimension_sizes(3)=current_state%local_grid%size(X_INDEX)
      field_information%data_type=COMPONENT_DOUBLE_DATA_TYPE

      if      (name .eq. "tend_u_viscosity_3d_local") then
        field_information%enabled=l_tend_3d_u
      else if (name .eq. "tend_v_viscosity_3d_local") then
        field_information%enabled=l_tend_3d_v
      else if (name .eq. "tend_w_viscosity_3d_local") then
        field_information%enabled=l_tend_3d_w
      else
        field_information%enabled=.true.
      end if

    end if !end 3d check

    ! Field information for profiles
    strcomp=INDEX(name, "viscosity_profile_total_local")
    if (strcomp .ne. 0) then
      field_information%field_type=COMPONENT_ARRAY_FIELD_TYPE
      field_information%number_dimensions=1
      field_information%dimension_sizes(1)=current_state%local_grid%size(Z_INDEX)
      field_information%data_type=COMPONENT_DOUBLE_DATA_TYPE

      if      (name .eq. "tend_u_viscosity_profile_total_local") then
        field_information%enabled=l_tend_pr_tot_u
      else if (name .eq. "tend_v_viscosity_profile_total_local") then
        field_information%enabled=l_tend_pr_tot_v
      else if (name .eq. "tend_w_viscosity_profile_total_local") then
        field_information%enabled=l_tend_pr_tot_w
      else
        field_information%enabled=.true.
      end if

    end if !end profile check

  end subroutine field_information_retrieval_callback


  !> Field value retrieval callback, this returns the value of a specific published field
  !! @param current_state Current model state
  !! @param name The name of the field to retrieve the value for
  !! @param field_value Populated with the value of the field
  subroutine field_value_retrieval_callback(current_state, name, field_value)
    type(model_state_type), target, intent(inout) :: current_state
    character(len=*), intent(in) :: name
    type(component_field_value_type), intent(out) :: field_value
    
    if (name .eq. "u_viscosity") then
      allocate(field_value%real_1d_array(size(u_viscosity)), source=u_viscosity)
    else if (name .eq. "v_viscosity") then
      allocate(field_value%real_1d_array(size(v_viscosity)), source=v_viscosity)
    else if (name .eq. "w_viscosity") then
      allocate(field_value%real_1d_array(size(w_viscosity)), source=w_viscosity)

    ! 3d Tendency Fields
    else if (name .eq. "tend_u_viscosity_3d_local" .and. allocated(tend_3d_u)) then
      call set_published_field_value(field_value, real_3d_field=tend_3d_u)
    else if (name .eq. "tend_v_viscosity_3d_local" .and. allocated(tend_3d_v)) then
      call set_published_field_value(field_value, real_3d_field=tend_3d_v)
    else if (name .eq. "tend_w_viscosity_3d_local" .and. allocated(tend_3d_w)) then
      call set_published_field_value(field_value, real_3d_field=tend_3d_w)

    ! Profile Tendency Fields
    else if (name .eq. "tend_u_viscosity_profile_total_local" .and. allocated(tend_pr_tot_u)) then
      call set_published_field_value(field_value, real_1d_field=tend_pr_tot_u)
    else if (name .eq. "tend_v_viscosity_profile_total_local" .and. allocated(tend_pr_tot_v)) then
      call set_published_field_value(field_value, real_1d_field=tend_pr_tot_v)
    else if (name .eq. "tend_w_viscosity_profile_total_local" .and. allocated(tend_pr_tot_w)) then
      call set_published_field_value(field_value, real_1d_field=tend_pr_tot_w)
    end if


  end subroutine field_value_retrieval_callback

  !> Sets up the stencil_mod (used in interpolation) and allocates data for the flux fields
  !! @param current_state The current model state_mod
  subroutine initialisation_callback(current_state)
    type(model_state_type), target, intent(inout) :: current_state

    integer :: z_size, y_size, x_size

    z_size=current_state%local_grid%size(Z_INDEX) + current_state%local_grid%halo_size(Z_INDEX) * 2
    y_size=current_state%local_grid%size(Y_INDEX) + current_state%local_grid%halo_size(Y_INDEX) * 2
    x_size=current_state%local_grid%size(X_INDEX) + current_state%local_grid%halo_size(X_INDEX) * 2
    allocate(current_state%vis_coefficient%data(z_size, y_size, x_size))

    z_size=current_state%global_grid%size(Z_INDEX)
    allocate(u_viscosity(z_size), v_viscosity(z_size), w_viscosity(z_size))

    ! Tendency Logicals
    l_tend_pr_tot_u   = current_state%u%active 
    l_tend_pr_tot_v   = current_state%v%active
    l_tend_pr_tot_w   = current_state%w%active

    l_tend_3d_u   = current_state%u%active .or. l_tend_pr_tot_u
    l_tend_3d_v   = current_state%v%active .or. l_tend_pr_tot_v
    l_tend_3d_w   = current_state%w%active .or. l_tend_pr_tot_w

    ! Allocate 3d tendency fields upon availability
    if (l_tend_3d_u) then
      allocate( tend_3d_u(current_state%local_grid%size(Z_INDEX),  &
                          current_state%local_grid%size(Y_INDEX),  &
                          current_state%local_grid%size(X_INDEX)   )    )
    endif
    if (l_tend_3d_v) then
      allocate( tend_3d_v(current_state%local_grid%size(Z_INDEX),  &
                          current_state%local_grid%size(Y_INDEX),  &
                          current_state%local_grid%size(X_INDEX)   )    )
    endif
    if (l_tend_3d_w) then
      allocate( tend_3d_w(current_state%local_grid%size(Z_INDEX),  &
                          current_state%local_grid%size(Y_INDEX),  &
                          current_state%local_grid%size(X_INDEX)   )    )
    endif

    ! Allocate profile tendency fields upon availability
    if (l_tend_pr_tot_u) then
      allocate( tend_pr_tot_u(current_state%local_grid%size(Z_INDEX)) )
    endif
    if (l_tend_pr_tot_v) then
      allocate( tend_pr_tot_v(current_state%local_grid%size(Z_INDEX)) )
    endif
    if (l_tend_pr_tot_w) then
      allocate( tend_pr_tot_w(current_state%local_grid%size(Z_INDEX)) )
    endif

    ! Save the sampling_frequency to force diagnostic calculation on select time steps
    diagnostic_generation_frequency=options_get_integer(current_state%options_database, "sampling_frequency")

  end subroutine initialisation_callback

  subroutine finalisation_callback(current_state)
    type(model_state_type), target, intent(inout) :: current_state

    if (allocated(u_viscosity)) deallocate(u_viscosity)
    if (allocated(v_viscosity)) deallocate(v_viscosity)
    if (allocated(w_viscosity)) deallocate(w_viscosity)

    if (allocated(tend_3d_u)) deallocate(tend_3d_u)
    if (allocated(tend_3d_v)) deallocate(tend_3d_v)
    if (allocated(tend_3d_w)) deallocate(tend_3d_w)

    if (allocated(tend_pr_tot_u)) deallocate(tend_pr_tot_u)
    if (allocated(tend_pr_tot_v)) deallocate(tend_pr_tot_v)
    if (allocated(tend_pr_tot_w)) deallocate(tend_pr_tot_w)

  end subroutine finalisation_callback  

  !> At each timestep will compute the viscosity U,V,W source terms
  !! @param current_state The current model state
  subroutine timestep_callback(current_state)
    type(model_state_type), target, intent(inout) :: current_state

    integer :: local_y, locaL_x, k, target_x_index, target_y_index
    real(kind=DEFAULT_PRECISION), dimension(current_state%local_grid%size(Z_INDEX)) :: tau12, tau12_ym1, tau12m1, &
         tau11, tau22, tau22_yp1, tau33, tau23_ym1, tau11p1, tau13, tau13m1, tau23

    local_y=current_state%column_local_y
    local_x=current_state%column_local_x
    target_y_index=local_y-current_state%local_grid%halo_size(Y_INDEX)
    target_x_index=local_x-current_state%local_grid%halo_size(X_INDEX)

    ! Zero profile tendency totals on first instance in the sum
    if (current_state%first_timestep_column) then
      if (l_tend_pr_tot_u) then
        tend_pr_tot_u(:)= 0.0_DEFAULT_PRECISION
      endif
      if (l_tend_pr_tot_v) then
        tend_pr_tot_v(:)= 0.0_DEFAULT_PRECISION
      endif
      if (l_tend_pr_tot_w) then
        tend_pr_tot_w(:)= 0.0_DEFAULT_PRECISION
      endif
    endif  ! zero totals

    if (.not. current_state%use_viscosity_and_diffusion .or. current_state%halo_column) return

    if (current_state%viscosity_halo_swap_state%swap_in_progress) then
      ! If there is a viscosity halo swap in progress then complete it
      call complete_nonblocking_halo_swap(current_state, current_state%viscosity_halo_swap_state, &
           perform_local_data_copy_for_vis, copy_halo_buffer_to_vis, copy_halo_buffer_to_vis_corners)
    end if

    if (mod(current_state%timestep, diagnostic_generation_frequency) == 0) then
      call save_precomponent_tendencies(current_state, local_x, local_y, target_x_index, target_y_index)
    end if

    if (current_state%field_stepping == FORWARD_STEPPING) then
      call calculate_tau(current_state, local_y, local_x, current_state%u, current_state%v, current_state%w, tau12, tau12_ym1, &
           tau12m1, tau11, tau22, tau22_yp1, tau33, tau11p1, tau13, tau13m1, tau23, tau23_ym1)
    else
      call calculate_tau(current_state, local_y, local_x, current_state%zu, current_state%zv, current_state%zw, tau12, tau12_ym1, &
           tau12m1, tau11, tau22, tau22_yp1, tau33, tau11p1, tau13, tau13m1, tau23, tau23_ym1)
    end if
    call calculate_viscous_sources(current_state, local_y, local_x, tau12, tau12_ym1, tau12m1, tau11, tau22, tau22_yp1, tau33, &
         tau11p1, tau13, tau13m1, tau23, tau23_ym1)

    if (mod(current_state%timestep, diagnostic_generation_frequency) == 0) then
      call compute_component_tendencies(current_state, local_x, local_y, target_x_index, target_y_index)
    end if

  end subroutine timestep_callback

  !> Calculates the viscous sources based upon TAU for the U,V and W fields
  !! @param current_state The current model state
  !! @param local_y Local Y index
  !! @param local_x Local X index
  subroutine calculate_viscous_sources(current_state, local_y, local_x, tau12, tau12_ym1, tau12m1, &
       tau11, tau22, tau22_yp1, tau33, tau11p1, tau13, tau13m1, tau23, tau23_ym1)
    type(model_state_type), target, intent(inout) :: current_state
    integer, intent(in) :: local_y, local_x
    real(kind=DEFAULT_PRECISION), dimension(:), intent(in) :: tau12, tau12_ym1, tau12m1, tau11, tau22, tau22_yp1, tau33, &
         tau11p1, tau13, tau13m1, tau23, tau23_ym1

    integer :: k

    do k=2, current_state%local_grid%size(Z_INDEX)
#ifdef U_ACTIVE
      u_viscosity(k)=((tau11p1(k)-tau11(k))*current_state%global_grid%configuration%horizontal%cx+(tau12(k)-tau12_ym1(k))*&
           current_state%global_grid%configuration%horizontal%cy+(tau13(k)-tau13(k-1))*&
           current_state%global_grid%configuration%vertical%rdz(k))/current_state%global_grid%configuration%vertical%rhon(k)
      current_state%su%data(k, local_y, local_x)=current_state%su%data(k, local_y, local_x)+u_viscosity(k)
#endif
#ifdef V_ACTIVE
      v_viscosity(k)=((tau12(k)-tau12m1(k))*current_state%global_grid%configuration%horizontal%cx+(tau22_yp1(k)-tau22(k))*&
           current_state%global_grid%configuration%horizontal%cy+(tau23(k)-tau23(k-1))*&
           current_state%global_grid%configuration%vertical%rdz(k))/current_state%global_grid%configuration%vertical%rhon(k)
      current_state%sv%data(k, local_y, local_x)=current_state%sv%data(k, local_y, local_x)+v_viscosity(k)
#endif
    end do
#ifdef W_ACTIVE
    do k=2, current_state%local_grid%size(Z_INDEX)-1
      w_viscosity(k)=((tau13(k)-tau13m1(k))*current_state%global_grid%configuration%horizontal%cx+(tau23(k)-tau23_ym1(k))*&
           current_state%global_grid%configuration%horizontal%cy+(tau33(k+1)-tau33(k))*&
           current_state%global_grid%configuration%vertical%rdzn(k+1))/current_state%global_grid%configuration%vertical%rho(k)
      current_state%sw%data(k, local_y, local_x)=current_state%sw%data(k, local_y, local_x)+w_viscosity(k)
    end do
#endif
  end subroutine calculate_viscous_sources  

  !> Calculated TAU which is used in computing the viscous source terms
  !! @param current_state The current model state
  !! @param local_y Local Y index
  !! @param local_x Local X index
  subroutine calculate_tau(current_state, local_y, local_x, zu, zv, zw, tau12, tau12_ym1, tau12m1, tau11, &
       tau22, tau22_yp1, tau33, tau11p1, tau13, tau13m1, tau23, tau23_ym1)
    type(model_state_type), target, intent(inout) :: current_state
    type(prognostic_field_type), intent(inout) :: zu, zv, zw
    integer, intent(in) :: local_y, local_x
    real(kind=DEFAULT_PRECISION), dimension(:), intent(out) :: tau12, tau12m1, tau11, tau22, tau22_yp1, tau33, &
         tau11p1, tau13, tau13m1, tau23, tau23_ym1, tau12_ym1

    integer :: k
    real(kind=DEFAULT_PRECISION) :: vistmp, vis12, vis12a, visonp2, visonp2a, vis13, vis23
    
    ! Do p levels and w-levels
    do k=1, current_state%local_grid%size(Z_INDEX)     
      if (k .gt. 1) then
        vistmp=current_state%vis_coefficient%data(k, local_y, local_x)+current_state%vis_coefficient%data(k, local_y+1, local_x)+&
             current_state%vis_coefficient%data(k-1, local_y, local_x)+current_state%vis_coefficient%data(k-1, local_y+1, local_x)
        vis12=0.125_DEFAULT_PRECISION*(vistmp + current_state%vis_coefficient%data(k, local_y, local_x+1)+&
             current_state%vis_coefficient%data(k, local_y+1, local_x+1)+&
             current_state%vis_coefficient%data(k-1, local_y, local_x+1)+&
             current_state%vis_coefficient%data(k-1, local_y+1, local_x+1))
        tau12(k)=0.0_DEFAULT_PRECISION
#ifdef U_ACTIVE
        tau12(k)=(zu%data(k, local_y+1, local_x)-zu%data(k, local_y, local_x))*&
             current_state%global_grid%configuration%horizontal%cy
#endif
#ifdef V_ACTIVE
        tau12(k)=tau12(k)+(zv%data(k, local_y, local_x+1)-zv%data(k, local_y, local_x))*&
             current_state%global_grid%configuration%horizontal%cx
#endif
        tau12(k)=tau12(k)*current_state%global_grid%configuration%vertical%rhon(k)*vis12
        
        vis12a=0.125_DEFAULT_PRECISION*(vistmp +current_state%vis_coefficient%data(k, local_y, local_x-1)+&
             current_state%vis_coefficient%data(k, local_y+1, local_x-1)+&
             current_state%vis_coefficient%data(k-1, local_y, local_x-1)+&
             current_state%vis_coefficient%data(k-1, local_y+1, local_x-1))
        tau12m1(k)=0.0_DEFAULT_PRECISION
#ifdef U_ACTIVE
        tau12m1(k)=(zu%data(k, local_y+1, local_x-1)-zu%data(k, local_y, local_x-1))*&
             current_state%global_grid%configuration%horizontal%cy
#endif
#ifdef V_ACTIVE
        tau12m1(k)=tau12m1(k)+(zv%data(k, local_y, local_x)-zv%data(k, local_y, local_x-1))*&
             current_state%global_grid%configuration%horizontal%cx
#endif
        tau12m1(k)=tau12m1(k)*current_state%global_grid%configuration%vertical%rhon(k)*vis12a

        vistmp=current_state%vis_coefficient%data(k, local_y-1, local_x)+current_state%vis_coefficient%data(k, local_y, local_x)+&
             current_state%vis_coefficient%data(k-1, local_y-1, local_x)+current_state%vis_coefficient%data(k-1, local_y, local_x)
        vis12=0.125_DEFAULT_PRECISION*(vistmp + current_state%vis_coefficient%data(k, local_y-1, local_x+1)+&
             current_state%vis_coefficient%data(k, local_y, local_x+1)+&
             current_state%vis_coefficient%data(k-1, local_y-1, local_x+1)+&
             current_state%vis_coefficient%data(k-1, local_y, local_x+1))
        tau12_ym1(k)=0.0_DEFAULT_PRECISION
#ifdef U_ACTIVE
        tau12_ym1(k)=(zu%data(k, local_y, local_x)-zu%data(k, local_y-1, local_x))*&
             current_state%global_grid%configuration%horizontal%cy
#endif
#ifdef V_ACTIVE
        tau12_ym1(k)=tau12_ym1(k)+(zv%data(k, local_y-1, local_x+1)-zv%data(k, local_y-1, local_x))*&
             current_state%global_grid%configuration%horizontal%cx
#endif
        tau12_ym1(k)=tau12_ym1(k)*current_state%global_grid%configuration%vertical%rhon(k)*vis12

        visonp2=current_state%global_grid%configuration%vertical%rhon(k)*&
             (current_state%vis_coefficient%data(k, local_y, local_x)+current_state%vis_coefficient%data(k-1, local_y, local_x))
#ifdef U_ACTIVE
        tau11(k)=visonp2*(zu%data(k, local_y, local_x)-zu%data(k, local_y, local_x-1))*&
             current_state%global_grid%configuration%horizontal%cx
#else
        tau11(k)=0.0_DEFAULT_PRECISION
#endif
#ifdef V_ACTIVE
        tau22(k)=visonp2*(zv%data(k, local_y, local_x)-zv%data(k, local_y-1, local_x))*&
             current_state%global_grid%configuration%horizontal%cy
#else
        tau22(k)=0.0_DEFAULT_PRECISION
#endif
#ifdef W_ACTIVE
        tau33(k)=visonp2*(zw%data(k, local_y, local_x)-zw%data(k-1, local_y, local_x))*&
             current_state%global_grid%configuration%vertical%rdz(k)
#else
        tau33(k)=0.0_DEFAULT_PRECISION
#endif

#ifdef V_ACTIVE
        visonp2=current_state%global_grid%configuration%vertical%rhon(k)*&
             (current_state%vis_coefficient%data(k, local_y+1, local_x)+&
             current_state%vis_coefficient%data(k-1, local_y+1, local_x))
        tau22_yp1(k)=visonp2*(zv%data(k, local_y+1, local_x)-zv%data(k, local_y, local_x))*&
             current_state%global_grid%configuration%horizontal%cy
#endif

#ifdef U_ACTIVE
        visonp2a=current_state%global_grid%configuration%vertical%rhon(k)*&
             (current_state%vis_coefficient%data(k, local_y, local_x+1)+&
             current_state%vis_coefficient%data(k-1, local_y, local_x+1))
        tau11p1(k)=visonp2a*(current_state%zu%data(k, local_y, local_x+1)-&
             zu%data(k, local_y, local_x))*current_state%global_grid%configuration%horizontal%cx
#endif
      else
        tau12(k)=0.0_DEFAULT_PRECISION
        tau12_ym1(k)=0.0_DEFAULT_PRECISION
        tau12m1(k)=0.0_DEFAULT_PRECISION
        tau11(k)=0.0_DEFAULT_PRECISION
        tau22(k)=0.0_DEFAULT_PRECISION
        tau22_yp1(k)=0.0_DEFAULT_PRECISION
        tau33(k)=0.0_DEFAULT_PRECISION
        tau11p1(k)=0.0_DEFAULT_PRECISION
      end if
      if (k .lt. current_state%local_grid%size(Z_INDEX)) then
        vis13=0.5_DEFAULT_PRECISION*(current_state%vis_coefficient%data(k, local_y, local_x)+&
             current_state%vis_coefficient%data(k, local_y, local_x+1))
        tau13(k)=0.0_DEFAULT_PRECISION
#ifdef U_ACTIVE
        tau13(k)=(zu%data(k+1, local_y, local_x)-zu%data(k, local_y, local_x))*&
             current_state%global_grid%configuration%vertical%rdzn(k+1)
#endif
#ifdef W_ACTIVE
        tau13(k)=tau13(k)+(zw%data(k, local_y, local_x+1)-zw%data(k, local_y, local_x))*&
             current_state%global_grid%configuration%horizontal%cx
#endif
        tau13(k)=tau13(k)*current_state%global_grid%configuration%vertical%rho(k)*vis13

        vis13=0.5_DEFAULT_PRECISION*(current_state%vis_coefficient%data(k, local_y, local_x-1)+&
             current_state%vis_coefficient%data(k, local_y, local_x))
        tau13m1(k)=0.0_DEFAULT_PRECISION
#ifdef U_ACTIVE
        tau13m1(k)=(zu%data(k+1, local_y, local_x-1)-zu%data(k, local_y, local_x-1))*&
             current_state%global_grid%configuration%vertical%rdzn(k+1)
#endif
#ifdef W_ACTIVE
        tau13m1(k)=tau13m1(k)+(zw%data(k, local_y, local_x)-zw%data(k, local_y, local_x-1))*&
             current_state%global_grid%configuration%horizontal%cx
#endif
        tau13m1(k)=tau13m1(k)*current_state%global_grid%configuration%vertical%rho(k)*vis13

        vis23=0.5_DEFAULT_PRECISION*(current_state%vis_coefficient%data(k, local_y, local_x)+&
             current_state%vis_coefficient%data(k, local_y+1, local_x))
        tau23(k)=0.0_DEFAULT_PRECISION
#ifdef W_ACTIVE
        tau23(k)=(zw%data(k, local_y+1, local_x)-zw%data(k, local_y, local_x))*&
             current_state%global_grid%configuration%horizontal%cy
#endif
#ifdef V_ACTIVE
        tau23(k)=tau23(k)+(zv%data(k+1, local_y, local_x)-zv%data(k, local_y, local_x))*&
             current_state%global_grid%configuration%vertical%rdzn(k+1)
#endif
        tau23(k)=tau23(k)*current_state%global_grid%configuration%vertical%rho(k)*vis23

        vis23=0.5_DEFAULT_PRECISION*(current_state%vis_coefficient%data(k, local_y-1, local_x)+&
             current_state%vis_coefficient%data(k, local_y, local_x))
        tau23_ym1(k)=0.0_DEFAULT_PRECISION
#ifdef W_ACTIVE
        tau23_ym1(k)=(zw%data(k, local_y, local_x)-zw%data(k, local_y-1, local_x))*&
             current_state%global_grid%configuration%horizontal%cy
#endif
#ifdef V_ACTIVE
        tau23_ym1(k)=tau23_ym1(k)+(zv%data(k+1, local_y-1, local_x)-zv%data(k, local_y-1, local_x))*&
             current_state%global_grid%configuration%vertical%rdzn(k+1)
#endif
        tau23_ym1(k)=tau23_ym1(k)*current_state%global_grid%configuration%vertical%rho(k)*vis23
      else
        tau13(k)=0.0_DEFAULT_PRECISION
        tau13m1(k)=0.0_DEFAULT_PRECISION
        tau23(k)=0.0_DEFAULT_PRECISION
        tau23_ym1(k)=0.0_DEFAULT_PRECISION
      end if
    end do
  end subroutine calculate_tau  

  !> Does local data copying for viscosity coefficient variable halo swap
  !! @param current_state The current model state_mod
  !! @param source_data Optional source data which is written into
  subroutine perform_local_data_copy_for_vis(current_state, halo_depth, involve_corners, source_data)
    type(model_state_type), intent(inout) :: current_state
    integer, intent(in) :: halo_depth
    logical, intent(in) :: involve_corners
    type(field_data_wrapper_type), dimension(:), intent(in), optional :: source_data

    call perform_local_data_copy_for_field(current_state%vis_coefficient%data, current_state%local_grid, &
         current_state%parallel%my_rank, halo_depth, involve_corners)
  end subroutine perform_local_data_copy_for_vis

  !> Copies the halo buffer to halo location for the viscosity coefficient field
  !! @param current_state The current model state
  !! @param neighbour_description The halo swapping description of the neighbour we are accessing the buffer of
  !! @param dim The dimension we receive for
  !! @param target_index The target index for the dimension we are receiving for
  !! @param neighbour_location The location in the local neighbour data stores of this neighbour
  !! @param current_page The current, next, halo swap page to read from (all previous have been read and copied already)
  !! @param source_data Optional source data which is written into
  subroutine copy_halo_buffer_to_vis(current_state, neighbour_description, dim, target_index, &
       neighbour_location, current_page, source_data)
    type(model_state_type), intent(inout) :: current_state
    integer, intent(in) :: dim, target_index, neighbour_location
    integer, intent(inout) :: current_page(:)
    type(neighbour_description_type), intent(inout) :: neighbour_description
    type(field_data_wrapper_type), dimension(:), intent(in), optional :: source_data

    call copy_buffer_to_field(current_state%local_grid, neighbour_description%recv_halo_buffer, &
         current_state%vis_coefficient%data, dim, target_index, current_page(neighbour_location))

    current_page(neighbour_location)=current_page(neighbour_location)+1
  end subroutine copy_halo_buffer_to_vis

  !> Copies the corner halo buffer to the viscosity coefficient field corners
  !! @param current_state The current model state
  !! @param neighbour_description The halo swapping description of the neighbour we are accessing the buffer of
  !! @param corner_loc The corner location
  !! @param x_target_index The X target index for the dimension we are receiving for
  !! @param y_target_index The Y target index for the dimension we are receiving for
  !! @param neighbour_location The location in the local neighbour data stores of this neighbour
  !! @param current_page The current, next, halo swap page to read from (all previous have been read and copied already)
  !! @param source_data Optional source data which is written into
  subroutine copy_halo_buffer_to_vis_corners(current_state, neighbour_description, corner_loc, x_target_index, &
       y_target_index, neighbour_location, current_page, source_data)
    type(model_state_type), intent(inout) :: current_state
    integer, intent(in) :: corner_loc, x_target_index, y_target_index, neighbour_location
    integer, intent(inout) :: current_page(:)
    type(neighbour_description_type), intent(inout) :: neighbour_description
    type(field_data_wrapper_type), dimension(:), intent(in), optional :: source_data

    call copy_buffer_to_corner(current_state%local_grid, neighbour_description%recv_corner_buffer, &
         current_state%vis_coefficient%data, corner_loc, x_target_index, y_target_index, current_page(neighbour_location))

    current_page(neighbour_location)=current_page(neighbour_location)+1
  end subroutine copy_halo_buffer_to_vis_corners  

   !> Save the 3d tendencies coming into the component.
  !! @param current_state Current model state
  !! @param cxn The current slice, x, index
  !! @param cyn The current column, y, index.
  !! @param txn target_x_index
  !! @param tyn target_y_index
  subroutine save_precomponent_tendencies(current_state, cxn, cyn, txn, tyn)
    type(model_state_type), target, intent(in) :: current_state
    integer, intent(in) ::  cxn, cyn, txn, tyn

    ! Save 3d tendency fields upon request (of 3d or profiles) and availability
    if (l_tend_3d_u) then
      tend_3d_u(:,tyn,txn)=current_state%su%data(:,cyn,cxn)
    endif
    if (l_tend_3d_v) then
      tend_3d_v(:,tyn,txn)=current_state%sv%data(:,cyn,cxn)
    endif
    if (l_tend_3d_w) then
      tend_3d_w(:,tyn,txn)=current_state%sw%data(:,cyn,cxn)
    endif

  end subroutine save_precomponent_tendencies


   !> Computation of component tendencies
  !! @param current_state Current model state
  !! @param cxn The current slice, x, index
  !! @param cyn The current column, y, index.
  !! @param txn target_x_index
  !! @param tyn target_y_index
  subroutine compute_component_tendencies(current_state, cxn, cyn, txn, tyn)
    type(model_state_type), target, intent(inout) :: current_state
    integer, intent(in) ::  cxn, cyn, txn, tyn

    ! Calculate change in tendency due to component
    if (l_tend_3d_u) then
      tend_3d_u(:,tyn,txn)=current_state%su%data(:,cyn,cxn) - tend_3d_u(:,tyn,txn)
    endif
    if (l_tend_3d_v) then
      tend_3d_v(:,tyn,txn)=current_state%sv%data(:,cyn,cxn) - tend_3d_v(:,tyn,txn)
    endif
    if (l_tend_3d_w) then
      tend_3d_w(:,tyn,txn)=current_state%sw%data(:,cyn,cxn) - tend_3d_w(:,tyn,txn)
    endif

   ! Add local tendency fields to the profile total
    if (l_tend_pr_tot_u) then
      tend_pr_tot_u(:)=tend_pr_tot_u(:) + tend_3d_u(:,tyn,txn)
    endif
    if (l_tend_pr_tot_v) then
      tend_pr_tot_v(:)=tend_pr_tot_v(:) + tend_3d_v(:,tyn,txn)
    endif
    if (l_tend_pr_tot_w) then
      tend_pr_tot_w(:)=tend_pr_tot_w(:) + tend_3d_w(:,tyn,txn)
    endif

  end subroutine compute_component_tendencies


  !> Sets the published field value from the temporary diagnostic values held by this component.
  !! @param field_value Populated with the value of the field
  !! @param real_1d_field Optional one dimensional real of values to publish
  !! @param real_2d_field Optional two dimensional real of values to publish
  subroutine set_published_field_value(field_value, real_1d_field, real_2d_field, real_3d_field)
    type(component_field_value_type), intent(inout) :: field_value
    real(kind=DEFAULT_PRECISION), dimension(:), optional :: real_1d_field
    real(kind=DEFAULT_PRECISION), dimension(:,:), optional :: real_2d_field
    real(kind=DEFAULT_PRECISION), dimension(:,:,:), optional :: real_3d_field

    if (present(real_1d_field)) then
      allocate(field_value%real_1d_array(size(real_1d_field)), source=real_1d_field)
    else if (present(real_2d_field)) then
      allocate(field_value%real_2d_array(size(real_2d_field, 1), size(real_2d_field, 2)), source=real_2d_field)
    else if (present(real_3d_field)) then
      allocate(field_value%real_3d_array(size(real_3d_field, 1), size(real_3d_field, 2), size(real_3d_field, 3)), &
               source=real_3d_field)
    end if
  end subroutine set_published_field_value

end module viscosity_mod
