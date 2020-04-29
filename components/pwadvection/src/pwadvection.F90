!> Piacsek-Williams advection scheme
module pwadvection_mod
  use datadefn_mod, only : DEFAULT_PRECISION
  use monc_component_mod, only : COMPONENT_ARRAY_FIELD_TYPE, &
      COMPONENT_DOUBLE_DATA_TYPE, component_descriptor_type, &
      component_field_value_type, component_field_information_type
  use optionsdatabase_mod, only : options_get_string, options_get_integer
  use collections_mod, only : map_type
  use state_mod, only : model_state_type
  use grids_mod, only : Z_INDEX, Y_INDEX, X_INDEX
  use q_indices_mod, only: get_q_index, standard_q_names

implicit none

#ifndef TEST_MODE
  private
#endif

  logical :: advect_flow, advect_th, advect_q

  logical :: l_toplevel=.true.

  ! Local tendency diagnostic variables for this component
  ! 3D tendency fields and logicals for their use
  real(kind=DEFAULT_PRECISION), dimension(:,:,:), allocatable ::     &
       tend_3d_u, tend_3d_v, tend_3d_w, tend_3d_th,tend_3d_qv,       &
       tend_3d_ql,tend_3d_qi,tend_3d_qr,tend_3d_qs,tend_3d_qg,       &
       tend_3d_tabs
  logical :: l_tend_3d_u, l_tend_3d_v, l_tend_3d_w, l_tend_3d_th,l_tend_3d_qv,       &
             l_tend_3d_ql,l_tend_3d_qi,l_tend_3d_qr,l_tend_3d_qs,l_tend_3d_qg,       &
             l_tend_3d_tabs
  ! Local mean tendency profile fields and logicals for their use
  real(kind=DEFAULT_PRECISION), dimension(:), allocatable ::                             &
       tend_pr_tot_u, tend_pr_tot_v, tend_pr_tot_w, tend_pr_tot_th,tend_pr_tot_qv,       &
       tend_pr_tot_ql,tend_pr_tot_qi,tend_pr_tot_qr,tend_pr_tot_qs,tend_pr_tot_qg,       &
       tend_pr_tot_tabs
  logical :: l_tend_pr_tot_u, l_tend_pr_tot_v, l_tend_pr_tot_w, l_tend_pr_tot_th,l_tend_pr_tot_qv,       &
             l_tend_pr_tot_ql,l_tend_pr_tot_qi,l_tend_pr_tot_qr,l_tend_pr_tot_qs,l_tend_pr_tot_qg,       &
             l_tend_pr_tot_tabs
  ! q indices
  integer :: iqv=0, iql=0, iqr=0, iqi=0, iqs=0, iqg=0
  integer :: diagnostic_generation_frequency


 public pwadvection_get_descriptor
contains

  !> Provides a description of this component for the core to register
  !! @returns The descriptor containing registration information for this component
  type(component_descriptor_type) function pwadvection_get_descriptor()
    pwadvection_get_descriptor%name="pw_advection"
    pwadvection_get_descriptor%version=0.1
    pwadvection_get_descriptor%initialisation=>initialisation_callback
    pwadvection_get_descriptor%timestep=>timestep_callback

    pwadvection_get_descriptor%field_value_retrieval=>field_value_retrieval_callback
    pwadvection_get_descriptor%field_information_retrieval=>field_information_retrieval_callback
    allocate(pwadvection_get_descriptor%published_fields(11+11))

    pwadvection_get_descriptor%published_fields(1)= "tend_u_pwadvection_3d_local"
    pwadvection_get_descriptor%published_fields(2)= "tend_v_pwadvection_3d_local"
    pwadvection_get_descriptor%published_fields(3)= "tend_w_pwadvection_3d_local"
    pwadvection_get_descriptor%published_fields(4)= "tend_th_pwadvection_3d_local"
    pwadvection_get_descriptor%published_fields(5)= "tend_qv_pwadvection_3d_local"
    pwadvection_get_descriptor%published_fields(6)= "tend_ql_pwadvection_3d_local"
    pwadvection_get_descriptor%published_fields(7)= "tend_qi_pwadvection_3d_local"
    pwadvection_get_descriptor%published_fields(8)= "tend_qr_pwadvection_3d_local"
    pwadvection_get_descriptor%published_fields(9)= "tend_qs_pwadvection_3d_local"
    pwadvection_get_descriptor%published_fields(10)="tend_qg_pwadvection_3d_local"
    pwadvection_get_descriptor%published_fields(11)="tend_tabs_pwadvection_3d_local"

    pwadvection_get_descriptor%published_fields(11+1)= "tend_u_pwadvection_profile_total_local"
    pwadvection_get_descriptor%published_fields(11+2)= "tend_v_pwadvection_profile_total_local"
    pwadvection_get_descriptor%published_fields(11+3)= "tend_w_pwadvection_profile_total_local"
    pwadvection_get_descriptor%published_fields(11+4)= "tend_th_pwadvection_profile_total_local"
    pwadvection_get_descriptor%published_fields(11+5)= "tend_qv_pwadvection_profile_total_local"
    pwadvection_get_descriptor%published_fields(11+6)= "tend_ql_pwadvection_profile_total_local"
    pwadvection_get_descriptor%published_fields(11+7)= "tend_qi_pwadvection_profile_total_local"
    pwadvection_get_descriptor%published_fields(11+8)= "tend_qr_pwadvection_profile_total_local"
    pwadvection_get_descriptor%published_fields(11+9)= "tend_qs_pwadvection_profile_total_local"
    pwadvection_get_descriptor%published_fields(11+10)="tend_qg_pwadvection_profile_total_local"
    pwadvection_get_descriptor%published_fields(11+11)="tend_tabs_pwadvection_profile_total_local"

  end function pwadvection_get_descriptor

  !> Initialisation callback, will set up the configuration of this advection scheme
  !! @param current_state The current model state
  subroutine initialisation_callback(current_state)
    type(model_state_type), target, intent(inout) :: current_state
    logical :: l_qdiag

    advect_flow=determine_if_advection_here(options_get_string(current_state%options_database, "advection_flow_fields"))
    advect_th=determine_if_advection_here(options_get_string(current_state%options_database, "advection_theta_field"))
    advect_q=determine_if_advection_here(options_get_string(current_state%options_database, "advection_q_fields"))

    ! Set tendency diagnostic logicals based on availability 
    ! Need to use 3d tendencies to compute the profiles, so they will be allocated 
    !      in the case where profiles are available
    l_qdiag =  (.not. current_state%passive_q .and. current_state%number_q_fields .gt. 0) .and. advect_q

    l_tend_pr_tot_u   = current_state%u%active .and. advect_flow 
    l_tend_pr_tot_v   = current_state%v%active .and. advect_flow
    l_tend_pr_tot_w   = current_state%w%active .and. advect_flow
    l_tend_pr_tot_th  = current_state%th%active .and. advect_th
    l_tend_pr_tot_qv  = l_qdiag .and. current_state%number_q_fields .ge. 1
    l_tend_pr_tot_ql  = l_qdiag .and. current_state%number_q_fields .ge. 2
    l_tend_pr_tot_qi  = l_qdiag .and. current_state%number_q_fields .ge. 11
    l_tend_pr_tot_qr  = l_qdiag .and. current_state%number_q_fields .ge. 11
    l_tend_pr_tot_qs  = l_qdiag .and. current_state%number_q_fields .ge. 11
    l_tend_pr_tot_qg  = l_qdiag .and. current_state%number_q_fields .ge. 11
    l_tend_pr_tot_tabs  = l_tend_pr_tot_th

    l_tend_3d_u   = (current_state%u%active .and. advect_flow) .or. l_tend_pr_tot_u 
    l_tend_3d_v   = (current_state%v%active .and. advect_flow) .or. l_tend_pr_tot_v
    l_tend_3d_w   = (current_state%w%active .and. advect_flow) .or. l_tend_pr_tot_w
    l_tend_3d_th  = (current_state%th%active .and. advect_th)  .or. l_tend_pr_tot_th
    l_tend_3d_qv  = (l_qdiag .and. current_state%number_q_fields .ge. 1) .or. l_tend_pr_tot_qv
    l_tend_3d_ql  = (l_qdiag .and. current_state%number_q_fields .ge. 2) .or. l_tend_pr_tot_ql
    l_tend_3d_qi  = (l_qdiag .and. current_state%number_q_fields .ge. 11) .or. l_tend_pr_tot_qi
    l_tend_3d_qr  = (l_qdiag .and. current_state%number_q_fields .ge. 11) .or. l_tend_pr_tot_qr
    l_tend_3d_qs  = (l_qdiag .and. current_state%number_q_fields .ge. 11) .or. l_tend_pr_tot_qs
    l_tend_3d_qg  = (l_qdiag .and. current_state%number_q_fields .ge. 11) .or. l_tend_pr_tot_qg
    l_tend_3d_tabs = l_tend_3d_th

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
    if (l_tend_3d_th) then
      allocate( tend_3d_th(current_state%local_grid%size(Z_INDEX),  &
                           current_state%local_grid%size(Y_INDEX),  &
                           current_state%local_grid%size(X_INDEX)   )    )
    endif
    if (l_tend_3d_qv) then
      iqv=get_q_index(standard_q_names%VAPOUR, 'pw_advection')
      allocate( tend_3d_qv(current_state%local_grid%size(Z_INDEX),  &
                           current_state%local_grid%size(Y_INDEX),  &
                           current_state%local_grid%size(X_INDEX)   )    )
    endif
    if (l_tend_3d_ql) then
      iql=get_q_index(standard_q_names%CLOUD_LIQUID_MASS, 'pw_advection')
      allocate( tend_3d_ql(current_state%local_grid%size(Z_INDEX),  &
                           current_state%local_grid%size(Y_INDEX),  &
                           current_state%local_grid%size(X_INDEX)   )    )
    endif
    if (l_tend_3d_qi) then
      iqi=get_q_index(standard_q_names%ICE_MASS, 'pw_advection')
      allocate( tend_3d_qi(current_state%local_grid%size(Z_INDEX),  &
                           current_state%local_grid%size(Y_INDEX),  &
                           current_state%local_grid%size(X_INDEX)   )    )
    endif
    if (l_tend_3d_qr) then
      iqr=get_q_index(standard_q_names%RAIN_MASS, 'pw_advection')
      allocate( tend_3d_qr(current_state%local_grid%size(Z_INDEX),  &
                           current_state%local_grid%size(Y_INDEX),  &
                           current_state%local_grid%size(X_INDEX)   )    )
    endif
    if (l_tend_3d_qs) then
      iqs=get_q_index(standard_q_names%SNOW_MASS, 'pw_advection')
      allocate( tend_3d_qs(current_state%local_grid%size(Z_INDEX),  &
                           current_state%local_grid%size(Y_INDEX),  &
                           current_state%local_grid%size(X_INDEX)   )    )
    endif
    if (l_tend_3d_qg) then
      iqg=get_q_index(standard_q_names%GRAUPEL_MASS, 'pw_advection')
      allocate( tend_3d_qg(current_state%local_grid%size(Z_INDEX),  &
                           current_state%local_grid%size(Y_INDEX),  &
                           current_state%local_grid%size(X_INDEX)   )    )
    endif
    if (l_tend_3d_tabs) then
      allocate( tend_3d_tabs(current_state%local_grid%size(Z_INDEX),  &
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
    if (l_tend_pr_tot_th) then
      allocate( tend_pr_tot_th(current_state%local_grid%size(Z_INDEX)) )
    endif
    if (l_tend_pr_tot_qv) then
      allocate( tend_pr_tot_qv(current_state%local_grid%size(Z_INDEX)) )
    endif
    if (l_tend_pr_tot_ql) then
      allocate( tend_pr_tot_ql(current_state%local_grid%size(Z_INDEX)) )
    endif
    if (l_tend_pr_tot_qi) then
      allocate( tend_pr_tot_qi(current_state%local_grid%size(Z_INDEX)) )
    endif
    if (l_tend_pr_tot_qr) then
      allocate( tend_pr_tot_qr(current_state%local_grid%size(Z_INDEX)) )
    endif
    if (l_tend_pr_tot_qs) then
      allocate( tend_pr_tot_qs(current_state%local_grid%size(Z_INDEX)) )
    endif
    if (l_tend_pr_tot_qg) then
      allocate( tend_pr_tot_qg(current_state%local_grid%size(Z_INDEX)) )
    endif
    if (l_tend_pr_tot_tabs) then
      allocate( tend_pr_tot_tabs(current_state%local_grid%size(Z_INDEX)) )
    endif

    ! Save the sampling_frequency to force diagnostic calculation on select time steps
    diagnostic_generation_frequency=options_get_integer(current_state%options_database, "sampling_frequency")

  end subroutine initialisation_callback


  subroutine finalisation_callback(current_state)
    type(model_state_type), target, intent(inout) :: current_state

    if (allocated(tend_3d_u)) deallocate(tend_3d_u)
    if (allocated(tend_3d_v)) deallocate(tend_3d_v)
    if (allocated(tend_3d_w)) deallocate(tend_3d_w)
    if (allocated(tend_3d_th)) deallocate(tend_3d_th)
    if (allocated(tend_3d_qv)) deallocate(tend_3d_qv)
    if (allocated(tend_3d_ql)) deallocate(tend_3d_ql)
    if (allocated(tend_3d_qi)) deallocate(tend_3d_qi)
    if (allocated(tend_3d_qr)) deallocate(tend_3d_qr)
    if (allocated(tend_3d_qs)) deallocate(tend_3d_qs)
    if (allocated(tend_3d_qg)) deallocate(tend_3d_qg)
    if (allocated(tend_3d_tabs)) deallocate(tend_3d_tabs)

    if (allocated(tend_pr_tot_u)) deallocate(tend_pr_tot_u)
    if (allocated(tend_pr_tot_v)) deallocate(tend_pr_tot_v)
    if (allocated(tend_pr_tot_w)) deallocate(tend_pr_tot_w)
    if (allocated(tend_pr_tot_th)) deallocate(tend_pr_tot_th)
    if (allocated(tend_pr_tot_qv)) deallocate(tend_pr_tot_qv)
    if (allocated(tend_pr_tot_ql)) deallocate(tend_pr_tot_ql)
    if (allocated(tend_pr_tot_qi)) deallocate(tend_pr_tot_qi)
    if (allocated(tend_pr_tot_qr)) deallocate(tend_pr_tot_qr)
    if (allocated(tend_pr_tot_qs)) deallocate(tend_pr_tot_qs)
    if (allocated(tend_pr_tot_qg)) deallocate(tend_pr_tot_qg)
    if (allocated(tend_pr_tot_tabs)) deallocate(tend_pr_tot_tabs)

  end subroutine finalisation_callback


  !> Called per column of data, this will perform Piacsek-Williams advection on the applicable fields for non halo data
  !! @param current_state The current model state
  !! @param target_(x/y)_index This is the index with the halos subtracted. This is needed so that diagnostic does 
  !!                           not include halos and to prevent array out-of-bounds
  subroutine timestep_callback(current_state)
    type(model_state_type), target, intent(inout) :: current_state

    integer :: current_x_index, current_y_index, target_x_index, target_y_index

    current_x_index=current_state%column_local_x
    current_y_index=current_state%column_local_y
    target_y_index=current_y_index-current_state%local_grid%halo_size(Y_INDEX)
    target_x_index=current_x_index-current_state%local_grid%halo_size(X_INDEX)


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
      if (l_tend_pr_tot_th) then
        tend_pr_tot_th(:)=0.0_DEFAULT_PRECISION
      endif
      if (l_tend_pr_tot_qv) then
        tend_pr_tot_qv(:)=0.0_DEFAULT_PRECISION
      endif
      if (l_tend_pr_tot_ql) then
        tend_pr_tot_ql(:)=0.0_DEFAULT_PRECISION
      endif
      if (l_tend_pr_tot_qi) then
        tend_pr_tot_qi(:)=0.0_DEFAULT_PRECISION
      endif
      if (l_tend_pr_tot_qr) then
        tend_pr_tot_qr(:)=0.0_DEFAULT_PRECISION
      endif
      if (l_tend_pr_tot_qs) then
        tend_pr_tot_qs(:)=0.0_DEFAULT_PRECISION
      endif
      if (l_tend_pr_tot_qg) then
        tend_pr_tot_qg(:)=0.0_DEFAULT_PRECISION
      endif
      if (l_tend_pr_tot_tabs) then
        tend_pr_tot_tabs(:)=0.0_DEFAULT_PRECISION
      endif
    endif  ! zero totals


    if (current_state%halo_column) return
    
    if (mod(current_state%timestep, diagnostic_generation_frequency) == 0) then
      call save_precomponent_tendencies(current_state, current_x_index, current_y_index, target_x_index, target_y_index)
    end if

    if (advect_flow) call advect_flow_fields(current_state, current_x_index, current_y_index)
    if (advect_th) call advect_th_field(current_state, current_x_index, current_y_index)
    if (advect_q) call advect_q_field(current_state, current_x_index, current_y_index)

    if (mod(current_state%timestep, diagnostic_generation_frequency) == 0) then
      call compute_component_tendencies(current_state, current_x_index, current_y_index, target_x_index, target_y_index)
    end if

  end subroutine timestep_callback

  !> Advects the q fields in a column
  !! @param current_state The current model state
  !! @param current_x_index The current slice, x, index
  !! @param current_y_index The current column, y, index
  subroutine advect_q_field(current_state, current_x_index, current_y_index)
    type(model_state_type), target, intent(inout) :: current_state
    integer, intent(in) ::  current_x_index, current_y_index

    integer :: k, n

    do n=1,current_state%number_q_fields
      do k=2,current_state%local_grid%size(Z_INDEX)-1
#ifdef U_ACTIVE
        current_state%sq(n)%data(k, current_y_index, current_x_index)=&
             current_state%sq(n)%data(k, current_y_index, current_x_index)+current_state%global_grid%configuration%horizontal%cx*&
             0.5_DEFAULT_PRECISION*(current_state%u%data(k, current_y_index, current_x_index-1)*&
             current_state%q(n)%data(k, current_y_index, current_x_index-1)-&
             current_state%u%data(k, current_y_index, current_x_index)*&
             current_state%q(n)%data(k, current_y_index, current_x_index+1))
#endif
#ifdef V_ACTIVE
        current_state%sq(n)%data(k, current_y_index, current_x_index)=&
             current_state%sq(n)%data(k, current_y_index, current_x_index)+&
             current_state%global_grid%configuration%horizontal%cy*0.5_DEFAULT_PRECISION*&
             (current_state%v%data(k, current_y_index-1, current_x_index)*&
             current_state%q(n)%data(k, current_y_index-1, current_x_index)-&
             current_state%v%data(k, current_y_index, current_x_index)*&
             current_state%q(n)%data(k, current_y_index+1, current_x_index))
#endif
#ifdef W_ACTIVE
        current_state%sq(n)%data(k, current_y_index, current_x_index)=&
             current_state%sq(n)%data(k, current_y_index, current_x_index)+&
             2.0_DEFAULT_PRECISION*(current_state%global_grid%configuration%vertical%tzc1(k)*&
             current_state%w%data(k-1, current_y_index, current_x_index)*&
             current_state%q(n)%data(k-1, current_y_index, current_x_index)-&
             current_state%global_grid%configuration%vertical%tzc2(k)*&
             current_state%w%data(k, current_y_index, current_x_index)*&
             current_state%q(n)%data(k+1, current_y_index, current_x_index))
#endif
      end do

      if (l_toplevel)then
      k=current_state%local_grid%size(Z_INDEX)
#ifdef U_ACTIVE
      current_state%sq(n)%data(k, current_y_index, current_x_index)=&
           current_state%sq(n)%data(k, current_y_index, current_x_index)+current_state%global_grid%configuration%horizontal%cx*&
           0.5_DEFAULT_PRECISION*(current_state%u%data(k, current_y_index, current_x_index-1)*&
           current_state%q(n)%data(k, current_y_index, current_x_index-1)-&
           current_state%u%data(k, current_y_index, current_x_index)*&
           current_state%q(n)%data(k, current_y_index, current_x_index+1))
#endif
#ifdef V_ACTIVE
      current_state%sq(n)%data(k, current_y_index, current_x_index)=&
           current_state%sq(n)%data(k, current_y_index, current_x_index)+current_state%global_grid%configuration%horizontal%cy*&
           0.5_DEFAULT_PRECISION*(current_state%v%data(k, current_y_index-1, current_x_index)*&
           current_state%q(n)%data(k, current_y_index-1, current_x_index)-&
           current_state%v%data(k, current_y_index, current_x_index)*&
           current_state%q(n)%data(k, current_y_index+1, current_x_index))
#endif
#ifdef W_ACTIVE
      current_state%sq(n)%data(k, current_y_index, current_x_index)=&
           current_state%sq(n)%data(k, current_y_index, current_x_index)+&
           current_state%global_grid%configuration%vertical%tzc1(k)*2.0_DEFAULT_PRECISION*&
           current_state%w%data(k-1, current_y_index, current_x_index)*&
           current_state%q(n)%data(k-1, current_y_index, current_x_index)
#endif
    end if
    end do
  end subroutine advect_q_field

  !> Advects the theta field in a column
  !! @param current_state The current model state
  !! @param current_x_index The current slice, x, index
  !! @param current_y_index The current column, y, index
  subroutine advect_th_field(current_state, current_x_index, current_y_index)
    type(model_state_type), target, intent(inout) :: current_state
    integer, intent(in) ::  current_x_index, current_y_index

    integer :: k

    if (current_state%th%active) then

      do k=2,current_state%local_grid%size(Z_INDEX)-1
#ifdef U_ACTIVE
        current_state%sth%data(k, current_y_index, current_x_index)= &   !current_state%sth%data(k, current_y_index, current_x_index)+&
             current_state%global_grid%configuration%horizontal%cx*&
             0.5_DEFAULT_PRECISION*(current_state%u%data(k, current_y_index, current_x_index-1)*&
             current_state%th%data(k, current_y_index, current_x_index-1)-&
             current_state%u%data(k, current_y_index, current_x_index)*&
             current_state%th%data(k, current_y_index, current_x_index+1))
#endif
#ifdef V_ACTIVE
        current_state%sth%data(k, current_y_index, current_x_index)=current_state%sth%data(k, current_y_index, current_x_index)+&
             current_state%global_grid%configuration%horizontal%cy*0.5_DEFAULT_PRECISION*&
             (current_state%v%data(k, current_y_index-1, current_x_index)*&
             current_state%th%data(k, current_y_index-1, current_x_index)-&
             current_state%v%data(k, current_y_index, current_x_index)*&
             current_state%th%data(k, current_y_index+1, current_x_index))
#endif
#ifdef W_ACTIVE
        current_state%sth%data(k, current_y_index, current_x_index)=current_state%sth%data(k, current_y_index, current_x_index)+&
             2.0_DEFAULT_PRECISION*(current_state%global_grid%configuration%vertical%tzc1(k)*&
             current_state%w%data(k-1, current_y_index, current_x_index)*&
             current_state%th%data(k-1, current_y_index, current_x_index)-&
             current_state%global_grid%configuration%vertical%tzc2(k)*&
             current_state%w%data(k, current_y_index, current_x_index)*&
             current_state%th%data(k+1, current_y_index, current_x_index))
#endif
      end do

      if (l_toplevel)then

        k=current_state%local_grid%size(Z_INDEX)
#ifdef U_ACTIVE
        current_state%sth%data(k, current_y_index, current_x_index)= &  !current_state%sth%data(k, current_y_index, current_x_index)+&
             current_state%global_grid%configuration%horizontal%cx*&
             0.5_DEFAULT_PRECISION*(current_state%u%data(k, current_y_index, current_x_index-1)*&
             current_state%th%data(k, current_y_index, current_x_index-1)-&
             current_state%u%data(k, current_y_index, current_x_index)*&
             current_state%th%data(k, current_y_index, current_x_index+1))
#endif
#ifdef V_ACTIVE
        current_state%sth%data(k, current_y_index, current_x_index)=current_state%sth%data(k, current_y_index, current_x_index)+&
             current_state%global_grid%configuration%horizontal%cy*0.5_DEFAULT_PRECISION*&
             (current_state%v%data(k, current_y_index-1, current_x_index)*&
             current_state%th%data(k, current_y_index-1, current_x_index)-&
             current_state%v%data(k, current_y_index, current_x_index)*&
             current_state%th%data(k, current_y_index+1, current_x_index))
#endif
#ifdef W_ACTIVE
        current_state%sth%data(k, current_y_index, current_x_index)=current_state%sth%data(k, current_y_index, current_x_index)+&
             current_state%global_grid%configuration%vertical%tzc1(k)*2.0_DEFAULT_PRECISION*&
             current_state%w%data(k-1, current_y_index, current_x_index)*current_state%th%data(k-1, current_y_index, &
             current_x_index)
#endif
      end if
   end if
  end subroutine advect_th_field

  !> Advects the flow fields depending upon which fields are active in the model in a column
  !! @param current_state The current model state
  !! @param current_x_index The current slice, x, index
  !! @param current_y_index The current column, y, index
  subroutine advect_flow_fields(current_state, current_x_index, current_y_index)
    type(model_state_type), target, intent(inout) :: current_state
    integer, intent(in) ::  current_x_index, current_y_index

    integer :: k

    do k=2,current_state%local_grid%size(Z_INDEX)-1
#ifdef U_ACTIVE
      current_state%su%data(k, current_y_index, current_x_index)=&
           current_state%global_grid%configuration%horizontal%tcx*(current_state%u%data(k, current_y_index, current_x_index-1)*&
           (current_state%u%data(k, current_y_index, current_x_index)+&
           current_state%u%data(k, current_y_index, current_x_index-1))-&
           current_state%u%data(k, current_y_index, current_x_index+1)*&
           (current_state%u%data(k, current_y_index, current_x_index)+&
           current_state%u%data(k, current_y_index, current_x_index+1)))
#ifdef V_ACTIVE
      current_state%su%data(k, current_y_index, current_x_index)=current_state%su%data(k, current_y_index, current_x_index)+&
           current_state%global_grid%configuration%horizontal%tcy*(current_state%u%data(k, current_y_index-1, current_x_index)*&
           (current_state%v%data(k, current_y_index-1, current_x_index)+&
           current_state%v%data(k, current_y_index-1, current_x_index+1))-&
           current_state%u%data(k, current_y_index+1, current_x_index)*&
           (current_state%v%data(k, current_y_index, current_x_index)+&
           current_state%v%data(k, current_y_index, current_x_index+1)))
#endif
#ifdef W_ACTIVE
      current_state%su%data(k, current_y_index, current_x_index)=current_state%su%data(k, current_y_index, current_x_index)+&
           (current_state%global_grid%configuration%vertical%tzc1(k)*current_state%u%data(k-1, current_y_index, current_x_index)*&
           (current_state%w%data(k-1, current_y_index, current_x_index)+&
           current_state%w%data(k-1, current_y_index, current_x_index+1))-&
           current_state%global_grid%configuration%vertical%tzc2(k)*current_state%u%data(k+1, current_y_index, current_x_index)*&
           (current_state%w%data(k, current_y_index, current_x_index)+&
           current_state%w%data(k, current_y_index, current_x_index+1)))
#endif
#endif

#ifdef V_ACTIVE
      current_state%sv%data(k, current_y_index, current_x_index)=current_state%global_grid%configuration%horizontal%tcy*(&
           current_state%v%data(k, current_y_index-1, current_x_index)*&
           (current_state%v%data(k, current_y_index, current_x_index)+&
           current_state%v%data(k, current_y_index-1, current_x_index))-&
           current_state%v%data(k, current_y_index+1, current_x_index)*&
           (current_state%v%data(k, current_y_index, current_x_index)+&
           current_state%v%data(k, current_y_index+1, current_x_index)))
#ifdef U_ACTIVE
      current_state%sv%data(k, current_y_index, current_x_index)=current_state%sv%data(k, current_y_index, current_x_index)+&
           current_state%global_grid%configuration%horizontal%tcx*(current_state%v%data(k, current_y_index, current_x_index-1)*&
           (current_state%u%data(k, current_y_index, current_x_index-1)+&
           current_state%u%data(k, current_y_index+1, current_x_index-1))-&
           current_state%v%data(k, current_y_index, current_x_index+1)*&
           (current_state%u%data(k, current_y_index, current_x_index)+&
           current_state%u%data(k, current_y_index+1, current_x_index)))
#endif
#ifdef W_ACTIVE
      current_state%sv%data(k, current_y_index, current_x_index)=current_state%sv%data(k, current_y_index, current_x_index)+&
           (current_state%global_grid%configuration%vertical%tzc1(k)*current_state%v%data(k-1, current_y_index, current_x_index)*&
           (current_state%w%data(k-1, current_y_index, current_x_index)+&
           current_state%w%data(k-1, current_y_index+1, current_x_index))-&
           current_state%global_grid%configuration%vertical%tzc2(k)*current_state%v%data(k+1, current_y_index, current_x_index)*&
           (current_state%w%data(k, current_y_index, current_x_index)+&
           current_state%w%data(k, current_y_index+1, current_x_index)))
#endif
#endif

#ifdef W_ACTIVE
      current_state%sw%data(k, current_y_index, current_x_index)=(current_state%global_grid%configuration%vertical%tzd1(k)*&
           current_state%w%data(k-1, current_y_index, current_x_index)*&
           (current_state%w%data(k, current_y_index, current_x_index)+&
           current_state%w%data(k-1, current_y_index, current_x_index))-&
           current_state%global_grid%configuration%vertical%tzd2(k)*current_state%w%data(k+1, current_y_index, current_x_index)*&
           (current_state%w%data(k, current_y_index, current_x_index)+&
           current_state%w%data(k+1, current_y_index, current_x_index)))
#ifdef U_ACTIVE
      current_state%sw%data(k, current_y_index, current_x_index)=current_state%sw%data(k, current_y_index, current_x_index)+&
           current_state%global_grid%configuration%horizontal%tcx*(current_state%w%data(k, current_y_index, current_x_index-1)*&
           (current_state%u%data(k, current_y_index, current_x_index-1)+&
           current_state%u%data(k+1, current_y_index, current_x_index-1))-&
           current_state%w%data(k, current_y_index, current_x_index+1)*&
           (current_state%u%data(k, current_y_index, current_x_index)+&
           current_state%u%data(k+1, current_y_index, current_x_index)))
#endif
#ifdef V_ACTIVE
      current_state%sw%data(k, current_y_index, current_x_index)=current_state%sw%data(k, current_y_index, current_x_index)+&
           current_state%global_grid%configuration%horizontal%tcy*(current_state%w%data(k, current_y_index-1, current_x_index)*&
           (current_state%v%data(k, current_y_index-1, current_x_index)+&
           current_state%v%data(k+1, current_y_index-1, current_x_index))-&
           current_state%w%data(k, current_y_index+1, current_x_index)*&
           (current_state%v%data(k, current_y_index, current_x_index)+&
           current_state%v%data(k+1, current_y_index, current_x_index)))
#endif
#endif
    end do

    if (l_toplevel)then
    k=current_state%local_grid%size(Z_INDEX)
#ifdef U_ACTIVE
    current_state%su%data(k, current_y_index, current_x_index)=current_state%global_grid%configuration%horizontal%tcx*&
         (current_state%u%data(k, current_y_index, current_x_index-1)*&
         (current_state%u%data(k, current_y_index, current_x_index)+&
         current_state%u%data(k, current_y_index, current_x_index-1))-&
         current_state%u%data(k, current_y_index, current_x_index+1)*&
         (current_state%u%data(k, current_y_index, current_x_index)+&
         current_state%u%data(k, current_y_index, current_x_index+1)))
#ifdef V_ACTIVE
    current_state%su%data(k, current_y_index, current_x_index)=current_state%su%data(k, current_y_index, current_x_index)+&
         current_state%global_grid%configuration%horizontal%tcy*(current_state%u%data(k, current_y_index-1, current_x_index)*&
         (current_state%v%data(k, current_y_index-1, current_x_index)+&
         current_state%v%data(k, current_y_index-1, current_x_index+1))-&
         current_state%u%data(k, current_y_index+1, current_x_index)*&
         (current_state%v%data(k, current_y_index, current_x_index)+&
         current_state%v%data(k, current_y_index, current_x_index+1)))
#endif
#ifdef W_ACTIVE
    current_state%su%data(k, current_y_index, current_x_index)=current_state%su%data(k, current_y_index, current_x_index)+&
         current_state%global_grid%configuration%vertical%tzc1(k)*current_state%u%data(k-1, current_y_index, current_x_index)*&
         (current_state%w%data(k-1, current_y_index, current_x_index)+&
         current_state%w%data(k-1, current_y_index, current_x_index+1))
#endif
#endif

#ifdef V_ACTIVE
    current_state%sv%data(k, current_y_index, current_x_index)=current_state%global_grid%configuration%horizontal%tcy*&
         (current_state%v%data(k, current_y_index-1, current_x_index)*&
         (current_state%v%data(k, current_y_index, current_x_index)+&
         current_state%v%data(k, current_y_index-1, current_x_index))-&
         current_state%v%data(k, current_y_index+1, current_x_index)*&
         (current_state%v%data(k, current_y_index, current_x_index)+&
         current_state%v%data(k, current_y_index+1, current_x_index)))
#ifdef U_ACTIVE
    current_state%sv%data(k, current_y_index, current_x_index)=current_state%sv%data(k, current_y_index, current_x_index)+&
         current_state%global_grid%configuration%horizontal%tcx*(current_state%v%data(k, current_y_index, current_x_index-1)*&
         (current_state%u%data(k, current_y_index, current_x_index-1)+&
         current_state%u%data(k, current_y_index+1, current_x_index-1))-&
         current_state%v%data(k, current_y_index, current_x_index+1)*&
         (current_state%u%data(k, current_y_index, current_x_index)+&
         current_state%u%data(k, current_y_index+1, current_x_index)))
#endif
#ifdef W_ACTIVE
    current_state%sv%data(k, current_y_index, current_x_index)=current_state%sv%data(k, current_y_index, current_x_index)+&
         current_state%global_grid%configuration%vertical%tzc1(k)*current_state%v%data(k-1, current_y_index, current_x_index)*&
         (current_state%w%data(k-1, current_y_index, current_x_index)+&
         current_state%w%data(k-1, current_y_index+1, current_x_index))
#endif
#endif
 end if
  end subroutine advect_flow_fields
  

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
    if (l_tend_3d_th) then
      tend_3d_th(:,tyn,txn)=current_state%sth%data(:,cyn,cxn)
    endif
    if (l_tend_3d_qv) then
      tend_3d_qv(:,tyn,txn)=current_state%sq(iqv)%data(:,cyn,cxn)
    endif
    if (l_tend_3d_ql) then
      tend_3d_ql(:,tyn,txn)=current_state%sq(iql)%data(:,cyn,cxn)
    endif
    if (l_tend_3d_qi) then
      tend_3d_qi(:,tyn,txn)=current_state%sq(iqi)%data(:,cyn,cxn)
    endif
    if (l_tend_3d_qr) then
      tend_3d_qr(:,tyn,txn)=current_state%sq(iqr)%data(:,cyn,cxn)
    endif
    if (l_tend_3d_qs) then
      tend_3d_qs(:,tyn,txn)=current_state%sq(iqs)%data(:,cyn,cxn)
    endif
    if (l_tend_3d_qg) then
      tend_3d_qg(:,tyn,txn)=current_state%sq(iqg)%data(:,cyn,cxn)
    endif
    if (l_tend_3d_tabs) then
      tend_3d_tabs(:,tyn,txn)=current_state%sth%data(:,cyn,cxn) * current_state%global_grid%configuration%vertical%rprefrcp(:)
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
      tend_3d_u(:,tyn,txn)=current_state%su%data(:,cyn,cxn)       - tend_3d_u(:,tyn,txn)
    endif
    if (l_tend_3d_v) then
      tend_3d_v(:,tyn,txn)=current_state%sv%data(:,cyn,cxn)       - tend_3d_v(:,tyn,txn)
    endif
    if (l_tend_3d_w) then
      tend_3d_w(:,tyn,txn)=current_state%sw%data(:,cyn,cxn)       - tend_3d_w(:,tyn,txn)
    endif
    if (l_tend_3d_th) then
      tend_3d_th(:,tyn,txn)=current_state%sth%data(:,cyn,cxn)     - tend_3d_th(:,tyn,txn)
    endif
    if (l_tend_3d_qv) then
      tend_3d_qv(:,tyn,txn)=current_state%sq(iqv)%data(:,cyn,cxn) - tend_3d_qv(:,tyn,txn)
    endif
    if (l_tend_3d_ql) then
      tend_3d_ql(:,tyn,txn)=current_state%sq(iql)%data(:,cyn,cxn) - tend_3d_ql(:,tyn,txn)
    endif
    if (l_tend_3d_qi) then
      tend_3d_qi(:,tyn,txn)=current_state%sq(iqi)%data(:,cyn,cxn) - tend_3d_qi(:,tyn,txn)
    endif
    if (l_tend_3d_qr) then
      tend_3d_qr(:,tyn,txn)=current_state%sq(iqr)%data(:,cyn,cxn) - tend_3d_qr(:,tyn,txn)
    endif
    if (l_tend_3d_qs) then
      tend_3d_qs(:,tyn,txn)=current_state%sq(iqs)%data(:,cyn,cxn) - tend_3d_qs(:,tyn,txn)
    endif
    if (l_tend_3d_qg) then
      tend_3d_qg(:,tyn,txn)=current_state%sq(iqg)%data(:,cyn,cxn) - tend_3d_qg(:,tyn,txn)
    endif
    if (l_tend_3d_tabs) then
      tend_3d_tabs(:,tyn,txn)=   &
       current_state%sth%data(:,cyn,cxn) * current_state%global_grid%configuration%vertical%rprefrcp(:)   &
        - tend_3d_tabs(:,tyn,txn)
    endif

   ! Add local tendency fields to the profile total 
    if (l_tend_pr_tot_u) then
      tend_pr_tot_u(:)=tend_pr_tot_u(:)   + tend_3d_u(:,tyn,txn)
    endif
    if (l_tend_pr_tot_v) then
      tend_pr_tot_v(:)=tend_pr_tot_v(:)   + tend_3d_v(:,tyn,txn)
    endif
    if (l_tend_pr_tot_w) then
      tend_pr_tot_w(:)=tend_pr_tot_w(:)   + tend_3d_w(:,tyn,txn)
    endif
    if (l_tend_pr_tot_th) then
      tend_pr_tot_th(:)=tend_pr_tot_th(:) + tend_3d_th(:,tyn,txn)
    endif
    if (l_tend_pr_tot_qv) then
      tend_pr_tot_qv(:)=tend_pr_tot_qv(:) + tend_3d_qv(:,tyn,txn)
    endif
    if (l_tend_pr_tot_ql) then
      tend_pr_tot_ql(:)=tend_pr_tot_ql(:) + tend_3d_ql(:,tyn,txn)
    endif
    if (l_tend_pr_tot_qi) then
      tend_pr_tot_qi(:)=tend_pr_tot_qi(:) + tend_3d_qi(:,tyn,txn)
    endif
    if (l_tend_pr_tot_qr) then
      tend_pr_tot_qr(:)=tend_pr_tot_qr(:) + tend_3d_qr(:,tyn,txn)
    endif
    if (l_tend_pr_tot_qs) then
      tend_pr_tot_qs(:)=tend_pr_tot_qs(:) + tend_3d_qs(:,tyn,txn)
    endif
    if (l_tend_pr_tot_qg) then
      tend_pr_tot_qg(:)=tend_pr_tot_qg(:) + tend_3d_qg(:,tyn,txn)
    endif
    if (l_tend_pr_tot_tabs) then
      tend_pr_tot_tabs(:)=tend_pr_tot_tabs(:) + tend_3d_tabs(:,tyn,txn)
    endif

  end subroutine compute_component_tendencies
  

   !> Field information retrieval callback, this returns information for a specific component's published field
  !! @param current_state Current model state
  !! @param name The name of the field to retrieve information for
  !! @param field_information Populated with information about the field
  !! @param strcomp Starting index within 1st argument string that matches substring (2nd argument); 0 if not a match.
  subroutine field_information_retrieval_callback(current_state, name, field_information)
    type(model_state_type), target, intent(inout) :: current_state
    character(len=*), intent(in) :: name
    type(component_field_information_type), intent(out) :: field_information
    integer :: strcomp

    ! Field information for 3d
    strcomp=INDEX(name, "pwadvection_3d_local")
    if (strcomp .ne. 0) then
      field_information%field_type=COMPONENT_ARRAY_FIELD_TYPE
      field_information%number_dimensions=3
      field_information%dimension_sizes(1)=current_state%local_grid%size(Z_INDEX)
      field_information%dimension_sizes(2)=current_state%local_grid%size(Y_INDEX)
      field_information%dimension_sizes(3)=current_state%local_grid%size(X_INDEX)
      field_information%data_type=COMPONENT_DOUBLE_DATA_TYPE
  
      if      (name .eq. "tend_u_pwadvection_3d_local") then
        field_information%enabled=l_tend_3d_u
      else if (name .eq. "tend_v_pwadvection_3d_local") then
        field_information%enabled=l_tend_3d_v
      else if (name .eq. "tend_w_pwadvection_3d_local") then
        field_information%enabled=l_tend_3d_w
      else if (name .eq. "tend_th_pwadvection_3d_local") then
        field_information%enabled=l_tend_3d_th
      else if (name .eq. "tend_qv_pwadvection_3d_local") then
        field_information%enabled=l_tend_3d_qv
      else if (name .eq. "tend_ql_pwadvection_3d_local") then
        field_information%enabled=l_tend_3d_ql
      else if (name .eq. "tend_qi_pwadvection_3d_local") then
        field_information%enabled=l_tend_3d_qi
      else if (name .eq. "tend_qr_pwadvection_3d_local") then
        field_information%enabled=l_tend_3d_qr
      else if (name .eq. "tend_qs_pwadvection_3d_local") then
        field_information%enabled=l_tend_3d_qs
      else if (name .eq. "tend_qg_pwadvection_3d_local") then
        field_information%enabled=l_tend_3d_qg
      else if (name .eq. "tend_tabs_pwadvection_3d_local") then
        field_information%enabled=l_tend_3d_tabs
      else 
        field_information%enabled=.true.
      end if
  
    end if !end 3d check

    ! Field information for profiles
    strcomp=INDEX(name, "pwadvection_profile_total_local")
    if (strcomp .ne. 0) then
      field_information%field_type=COMPONENT_ARRAY_FIELD_TYPE
      field_information%number_dimensions=1
      field_information%dimension_sizes(1)=current_state%local_grid%size(Z_INDEX)
      field_information%data_type=COMPONENT_DOUBLE_DATA_TYPE
  
      if      (name .eq. "tend_u_pwadvection_profile_total_local") then
        field_information%enabled=l_tend_pr_tot_u
      else if (name .eq. "tend_v_pwadvection_profile_total_local") then
        field_information%enabled=l_tend_pr_tot_v
      else if (name .eq. "tend_w_pwadvection_profile_total_local") then
        field_information%enabled=l_tend_pr_tot_w
      else if (name .eq. "tend_th_pwadvection_profile_total_local") then
        field_information%enabled=l_tend_pr_tot_th
      else if (name .eq. "tend_qv_pwadvection_profile_total_local") then
        field_information%enabled=l_tend_pr_tot_qv
      else if (name .eq. "tend_ql_pwadvection_profile_total_local") then
        field_information%enabled=l_tend_pr_tot_ql
      else if (name .eq. "tend_qi_pwadvection_profile_total_local") then
        field_information%enabled=l_tend_pr_tot_qi
      else if (name .eq. "tend_qr_pwadvection_profile_total_local") then
        field_information%enabled=l_tend_pr_tot_qr
      else if (name .eq. "tend_qs_pwadvection_profile_total_local") then
        field_information%enabled=l_tend_pr_tot_qs
      else if (name .eq. "tend_qg_pwadvection_profile_total_local") then
        field_information%enabled=l_tend_pr_tot_qg
      else if (name .eq. "tend_tabs_pwadvection_profile_total_local") then
        field_information%enabled=l_tend_pr_tot_tabs
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
    
    ! 3d Tendency Fields
    if (name .eq. "tend_u_pwadvection_3d_local" .and. allocated(tend_3d_u)) then
      call set_published_field_value(field_value, real_3d_field=tend_3d_u)
    else if (name .eq. "tend_v_pwadvection_3d_local" .and. allocated(tend_3d_v)) then
      call set_published_field_value(field_value, real_3d_field=tend_3d_v)
    else if (name .eq. "tend_w_pwadvection_3d_local" .and. allocated(tend_3d_w)) then
      call set_published_field_value(field_value, real_3d_field=tend_3d_w)
    else if (name .eq. "tend_th_pwadvection_3d_local" .and. allocated(tend_3d_th)) then
      call set_published_field_value(field_value, real_3d_field=tend_3d_th)
    else if (name .eq. "tend_qv_pwadvection_3d_local" .and. allocated(tend_3d_qv)) then
      call set_published_field_value(field_value, real_3d_field=tend_3d_qv)
    else if (name .eq. "tend_ql_pwadvection_3d_local" .and. allocated(tend_3d_ql)) then
      call set_published_field_value(field_value, real_3d_field=tend_3d_ql)
    else if (name .eq. "tend_qi_pwadvection_3d_local" .and. allocated(tend_3d_qi)) then
      call set_published_field_value(field_value, real_3d_field=tend_3d_qi)
    else if (name .eq. "tend_qr_pwadvection_3d_local" .and. allocated(tend_3d_qr)) then
      call set_published_field_value(field_value, real_3d_field=tend_3d_qr)
    else if (name .eq. "tend_qs_pwadvection_3d_local" .and. allocated(tend_3d_qs)) then
      call set_published_field_value(field_value, real_3d_field=tend_3d_qs)
    else if (name .eq. "tend_qg_pwadvection_3d_local" .and. allocated(tend_3d_qg)) then
      call set_published_field_value(field_value, real_3d_field=tend_3d_qg)
    else if (name .eq. "tend_tabs_pwadvection_3d_local" .and. allocated(tend_3d_tabs)) then
      call set_published_field_value(field_value, real_3d_field=tend_3d_tabs)

    ! Profile Tendency Fields	
    else if (name .eq. "tend_u_pwadvection_profile_total_local" .and. allocated(tend_pr_tot_u)) then
      call set_published_field_value(field_value, real_1d_field=tend_pr_tot_u)
    else if (name .eq. "tend_v_pwadvection_profile_total_local" .and. allocated(tend_pr_tot_v)) then
      call set_published_field_value(field_value, real_1d_field=tend_pr_tot_v)
    else if (name .eq. "tend_w_pwadvection_profile_total_local" .and. allocated(tend_pr_tot_w)) then
      call set_published_field_value(field_value, real_1d_field=tend_pr_tot_w)
    else if (name .eq. "tend_th_pwadvection_profile_total_local" .and. allocated(tend_pr_tot_th)) then
      call set_published_field_value(field_value, real_1d_field=tend_pr_tot_th)
    else if (name .eq. "tend_qv_pwadvection_profile_total_local" .and. allocated(tend_pr_tot_qv)) then
      call set_published_field_value(field_value, real_1d_field=tend_pr_tot_qv)
    else if (name .eq. "tend_ql_pwadvection_profile_total_local" .and. allocated(tend_pr_tot_ql)) then
      call set_published_field_value(field_value, real_1d_field=tend_pr_tot_ql)
    else if (name .eq. "tend_qi_pwadvection_profile_total_local" .and. allocated(tend_pr_tot_qi)) then
      call set_published_field_value(field_value, real_1d_field=tend_pr_tot_qi)
    else if (name .eq. "tend_qr_pwadvection_profile_total_local" .and. allocated(tend_pr_tot_qr)) then
      call set_published_field_value(field_value, real_1d_field=tend_pr_tot_qr)
    else if (name .eq. "tend_qs_pwadvection_profile_total_local" .and. allocated(tend_pr_tot_qs)) then
      call set_published_field_value(field_value, real_1d_field=tend_pr_tot_qs)
    else if (name .eq. "tend_qg_pwadvection_profile_total_local" .and. allocated(tend_pr_tot_qg)) then
      call set_published_field_value(field_value, real_1d_field=tend_pr_tot_qg)
    else if (name .eq. "tend_tabs_pwadvection_profile_total_local" .and. allocated(tend_pr_tot_tabs)) then
      call set_published_field_value(field_value, real_1d_field=tend_pr_tot_tabs)
    end if

  end subroutine field_value_retrieval_callback  

  
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

  
  !> Parses a field string (read in from the configuration file) and determines whether this algorithm should be used
  !! for advecting that field
  !! @param field The string configuration of field advection
  !! @returns Whether or not the field is advected here
  logical function determine_if_advection_here(field)
    character(len=*), intent(in) :: field

    if (len_trim(field) .ne. 0) then
      if (trim(field) .eq. "pw" .or. trim(field) .eq. "any") then
        determine_if_advection_here=.true.
      else
        determine_if_advection_here=.false.
      end if
    else
      determine_if_advection_here=.true.
    end if
  end function determine_if_advection_here

end module pwadvection_mod

