!> Implements TVD advection for prognostic fields
module tvdadvection_mod
  use monc_component_mod, only : COMPONENT_ARRAY_FIELD_TYPE, &
      COMPONENT_DOUBLE_DATA_TYPE, component_descriptor_type, &
      component_field_value_type, component_field_information_type
  use stencil_mod, only : grid_stencil_type, interpolate_to_dual, create_stencil, free_stencil
  use state_mod, only : model_state_type, parallel_state_type, FORWARD_STEPPING
  use grids_mod, only : local_grid_type, global_grid_type, X_INDEX, Y_INDEX, Z_INDEX
  use prognostics_mod, only : prognostic_field_type, prognostic_field_ptr_type
  use ultimateflux_mod, only : ultflx
  use datadefn_mod, only : DEFAULT_PRECISION, PRECISION_TYPE
  use optionsdatabase_mod, only : options_get_string, options_get_integer
  use collections_mod, only : map_type
  use mpi, only : MPI_REQUEST_NULL, MPI_STATUS_IGNORE
  use q_indices_mod, only: get_q_index, standard_q_names
  ! Some tvd diagnostic terms 
  use def_tvd_diagnostic_terms, only: tvd_dgs_terms 
  use registry_mod, only : is_component_enabled

  implicit none

#ifndef TEST_MODE
  private
#endif

  type(grid_stencil_type), save :: star_stencil
  integer, save :: u_index=0, v_index=0, w_index=0
  logical :: advect_flow, advect_th, advect_q, advect_tracer
  real(kind=DEFAULT_PRECISION), dimension(:), allocatable :: flux_x, flux_y, flux_z, u_advection, v_advection, &
       w_advection, th_advection
  real(kind=DEFAULT_PRECISION), dimension(:,:), allocatable :: q_advection
  real(kind=DEFAULT_PRECISION), dimension(:,:), allocatable :: tracer_advection
  
  type(prognostic_field_type), dimension(:), allocatable :: interpolated_fields

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
       tend_pr_tot_u, tend_pr_tot_v, tend_pr_tot_w,tend_pr_tot_th,tend_pr_tot_qv,        &
       tend_pr_tot_ql,tend_pr_tot_qi,tend_pr_tot_qr,tend_pr_tot_qs,tend_pr_tot_qg,       &
       tend_pr_tot_tabs
  logical :: l_tend_pr_tot_u, l_tend_pr_tot_v, l_tend_pr_tot_w,l_tend_pr_tot_th,l_tend_pr_tot_qv,        &
             l_tend_pr_tot_ql,l_tend_pr_tot_qi,l_tend_pr_tot_qr,l_tend_pr_tot_qs,l_tend_pr_tot_qg,       &
             l_tend_pr_tot_tabs
  ! q indices
  integer :: iqv=0, iql=0, iqr=0, iqi=0, iqs=0, iqg=0

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
    allocate(tvdadvection_get_descriptor%published_fields(5+11+11+1))
    tvdadvection_get_descriptor%published_fields(1)="u_advection"
    tvdadvection_get_descriptor%published_fields(2)="v_advection"
    tvdadvection_get_descriptor%published_fields(3)="w_advection"
    tvdadvection_get_descriptor%published_fields(4)="th_advection"
    tvdadvection_get_descriptor%published_fields(5)="q_advection"
    

    tvdadvection_get_descriptor%published_fields(5+1)="tend_u_tvdadvection_3d_local"
    tvdadvection_get_descriptor%published_fields(5+2)="tend_v_tvdadvection_3d_local"
    tvdadvection_get_descriptor%published_fields(5+3)="tend_w_tvdadvection_3d_local"
    tvdadvection_get_descriptor%published_fields(5+4)="tend_th_tvdadvection_3d_local"
    tvdadvection_get_descriptor%published_fields(5+5)="tend_qv_tvdadvection_3d_local"
    tvdadvection_get_descriptor%published_fields(5+6)="tend_ql_tvdadvection_3d_local"
    tvdadvection_get_descriptor%published_fields(5+7)="tend_qi_tvdadvection_3d_local"
    tvdadvection_get_descriptor%published_fields(5+8)="tend_qr_tvdadvection_3d_local"
    tvdadvection_get_descriptor%published_fields(5+9)="tend_qs_tvdadvection_3d_local"
    tvdadvection_get_descriptor%published_fields(5+10)="tend_qg_tvdadvection_3d_local"
    tvdadvection_get_descriptor%published_fields(5+11)="tend_tabs_tvdadvection_3d_local"

    tvdadvection_get_descriptor%published_fields(5+11+1)="tend_u_tvdadvection_profile_total_local"
    tvdadvection_get_descriptor%published_fields(5+11+2)="tend_v_tvdadvection_profile_total_local"
    tvdadvection_get_descriptor%published_fields(5+11+3)="tend_w_tvdadvection_profile_total_local"
    tvdadvection_get_descriptor%published_fields(5+11+4)="tend_th_tvdadvection_profile_total_local"
    tvdadvection_get_descriptor%published_fields(5+11+5)="tend_qv_tvdadvection_profile_total_local"
    tvdadvection_get_descriptor%published_fields(5+11+6)="tend_ql_tvdadvection_profile_total_local"
    tvdadvection_get_descriptor%published_fields(5+11+7)="tend_qi_tvdadvection_profile_total_local"
    tvdadvection_get_descriptor%published_fields(5+11+8)="tend_qr_tvdadvection_profile_total_local"
    tvdadvection_get_descriptor%published_fields(5+11+9)="tend_qs_tvdadvection_profile_total_local"
    tvdadvection_get_descriptor%published_fields(5+11+10)="tend_qg_tvdadvection_profile_total_local"
    tvdadvection_get_descriptor%published_fields(5+11+11)="tend_tabs_tvdadvection_profile_total_local"
    
    tvdadvection_get_descriptor%published_fields(5+11+11+1)="tracer_advection"


  end function tvdadvection_get_descriptor

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
    if (name .eq. "q_advection") then
      field_information%number_dimensions=2
    else if (name .eq. "tracer_advection") then
      field_information%number_dimensions=2
    else
      field_information%number_dimensions=1
    end if
    field_information%dimension_sizes(1)=current_state%local_grid%size(Z_INDEX)
    if (name .eq. "q_advection") field_information%dimension_sizes(2)=current_state%number_q_fields
    if (name .eq. "tracer_advection") field_information%dimension_sizes(2)=current_state%n_tracers
    field_information%enabled=.true.

    ! Field information for 3d
    strcomp=INDEX(name, "tvdadvection_3d_local")
    if (strcomp .ne. 0) then
      field_information%field_type=COMPONENT_ARRAY_FIELD_TYPE
      field_information%number_dimensions=3
      field_information%dimension_sizes(1)=current_state%local_grid%size(Z_INDEX)
      field_information%dimension_sizes(2)=current_state%local_grid%size(Y_INDEX)
      field_information%dimension_sizes(3)=current_state%local_grid%size(X_INDEX)
      field_information%data_type=COMPONENT_DOUBLE_DATA_TYPE

      if      (name .eq. "tend_u_tvdadvection_3d_local") then
        field_information%enabled=l_tend_3d_u
      else if (name .eq. "tend_v_tvdadvection_3d_local") then
        field_information%enabled=l_tend_3d_v
      else if (name .eq. "tend_w_tvdadvection_3d_local") then
        field_information%enabled=l_tend_3d_w
      else if (name .eq. "tend_th_tvdadvection_3d_local") then
        field_information%enabled=l_tend_3d_th
      else if (name .eq. "tend_qv_tvdadvection_3d_local") then
        field_information%enabled=l_tend_3d_qv
      else if (name .eq. "tend_ql_tvdadvection_3d_local") then
        field_information%enabled=l_tend_3d_ql
      else if (name .eq. "tend_qi_tvdadvection_3d_local") then
        field_information%enabled=l_tend_3d_qi
      else if (name .eq. "tend_qr_tvdadvection_3d_local") then
        field_information%enabled=l_tend_3d_qr
      else if (name .eq. "tend_qs_tvdadvection_3d_local") then
        field_information%enabled=l_tend_3d_qs
      else if (name .eq. "tend_qg_tvdadvection_3d_local") then
        field_information%enabled=l_tend_3d_qg
      else if (name .eq. "tend_tabs_tvdadvection_3d_local") then
        field_information%enabled=l_tend_3d_tabs
      else
        field_information%enabled=.true.
      end if

    end if !end 3d check

    ! Field information for profiles
    strcomp=INDEX(name, "tvdadvection_profile_total_local")
    if (strcomp .ne. 0) then
      field_information%field_type=COMPONENT_ARRAY_FIELD_TYPE
      field_information%number_dimensions=1
      field_information%dimension_sizes(1)=current_state%local_grid%size(Z_INDEX)
      field_information%data_type=COMPONENT_DOUBLE_DATA_TYPE

      if      (name .eq. "tend_u_tvdadvection_profile_total_local") then
        field_information%enabled=l_tend_pr_tot_u
      else if (name .eq. "tend_v_tvdadvection_profile_total_local") then
        field_information%enabled=l_tend_pr_tot_v
      else if (name .eq. "tend_w_tvdadvection_profile_total_local") then
        field_information%enabled=l_tend_pr_tot_w
      else if (name .eq. "tend_th_tvdadvection_profile_total_local") then
        field_information%enabled=l_tend_pr_tot_th
      else if (name .eq. "tend_qv_tvdadvection_profile_total_local") then
        field_information%enabled=l_tend_pr_tot_qv
      else if (name .eq. "tend_ql_tvdadvection_profile_total_local") then
        field_information%enabled=l_tend_pr_tot_ql
      else if (name .eq. "tend_qi_tvdadvection_profile_total_local") then
        field_information%enabled=l_tend_pr_tot_qi
      else if (name .eq. "tend_qr_tvdadvection_profile_total_local") then
        field_information%enabled=l_tend_pr_tot_qr
      else if (name .eq. "tend_qs_tvdadvection_profile_total_local") then
        field_information%enabled=l_tend_pr_tot_qs
      else if (name .eq. "tend_qg_tvdadvection_profile_total_local") then
        field_information%enabled=l_tend_pr_tot_qg
      else if (name .eq. "tend_tabs_tvdadvection_profile_total_local") then
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
    else if (name .eq. "tracer_advection") then
      allocate(field_value%real_2d_array(size(tracer_advection, 1), size(tracer_advection, 2)), source=tracer_advection)
    end if

    ! 3d Tendency Fields
    if      (name .eq. "tend_u_tvdadvection_3d_local" .and. allocated(tend_3d_u)) then
      call set_published_field_value(field_value, real_3d_field=tend_3d_u)
    else if (name .eq. "tend_v_tvdadvection_3d_local" .and. allocated(tend_3d_v)) then
      call set_published_field_value(field_value, real_3d_field=tend_3d_v)
    else if (name .eq. "tend_w_tvdadvection_3d_local" .and. allocated(tend_3d_w)) then
      call set_published_field_value(field_value, real_3d_field=tend_3d_w)
    else if (name .eq. "tend_th_tvdadvection_3d_local" .and. allocated(tend_3d_th)) then
      call set_published_field_value(field_value, real_3d_field=tend_3d_th)
    else if (name .eq. "tend_qv_tvdadvection_3d_local" .and. allocated(tend_3d_qv)) then
      call set_published_field_value(field_value, real_3d_field=tend_3d_qv)
    else if (name .eq. "tend_ql_tvdadvection_3d_local" .and. allocated(tend_3d_ql)) then
      call set_published_field_value(field_value, real_3d_field=tend_3d_ql)
    else if (name .eq. "tend_qi_tvdadvection_3d_local" .and. allocated(tend_3d_qi)) then
      call set_published_field_value(field_value, real_3d_field=tend_3d_qi)
    else if (name .eq. "tend_qr_tvdadvection_3d_local" .and. allocated(tend_3d_qr)) then
      call set_published_field_value(field_value, real_3d_field=tend_3d_qr)
    else if (name .eq. "tend_qs_tvdadvection_3d_local" .and. allocated(tend_3d_qs)) then
      call set_published_field_value(field_value, real_3d_field=tend_3d_qs)
    else if (name .eq. "tend_qg_tvdadvection_3d_local" .and. allocated(tend_3d_qg)) then
      call set_published_field_value(field_value, real_3d_field=tend_3d_qg)
    else if (name .eq. "tend_tabs_tvdadvection_3d_local" .and. allocated(tend_3d_tabs)) then
      call set_published_field_value(field_value, real_3d_field=tend_3d_tabs)

    ! Profile Tendency Fields
    else if (name .eq. "tend_u_tvdadvection_profile_total_local" .and. allocated(tend_pr_tot_u)) then
      call set_published_field_value(field_value, real_1d_field=tend_pr_tot_u)
    else if (name .eq. "tend_v_tvdadvection_profile_total_local" .and. allocated(tend_pr_tot_v)) then
      call set_published_field_value(field_value, real_1d_field=tend_pr_tot_v)
    else if (name .eq. "tend_w_tvdadvection_profile_total_local" .and. allocated(tend_pr_tot_w)) then
      call set_published_field_value(field_value, real_1d_field=tend_pr_tot_w)
    else if (name .eq. "tend_th_tvdadvection_profile_total_local" .and. allocated(tend_pr_tot_th)) then
      call set_published_field_value(field_value, real_1d_field=tend_pr_tot_th)
    else if (name .eq. "tend_qv_tvdadvection_profile_total_local" .and. allocated(tend_pr_tot_qv)) then
      call set_published_field_value(field_value, real_1d_field=tend_pr_tot_qv)
    else if (name .eq. "tend_ql_tvdadvection_profile_total_local" .and. allocated(tend_pr_tot_ql)) then
      call set_published_field_value(field_value, real_1d_field=tend_pr_tot_ql)
    else if (name .eq. "tend_qi_tvdadvection_profile_total_local" .and. allocated(tend_pr_tot_qi)) then
      call set_published_field_value(field_value, real_1d_field=tend_pr_tot_qi)
    else if (name .eq. "tend_qr_tvdadvection_profile_total_local" .and. allocated(tend_pr_tot_qr)) then
      call set_published_field_value(field_value, real_1d_field=tend_pr_tot_qr)
    else if (name .eq. "tend_qs_tvdadvection_profile_total_local" .and. allocated(tend_pr_tot_qs)) then
      call set_published_field_value(field_value, real_1d_field=tend_pr_tot_qs)
    else if (name .eq. "tend_qg_tvdadvection_profile_total_local" .and. allocated(tend_pr_tot_qg)) then
      call set_published_field_value(field_value, real_1d_field=tend_pr_tot_qg)
    else if (name .eq. "tend_tabs_tvdadvection_profile_total_local" .and. allocated(tend_pr_tot_tabs)) then
      call set_published_field_value(field_value, real_1d_field=tend_pr_tot_tabs)
    end if

  end subroutine field_value_retrieval_callback

  !> Sets up the stencil_mod (used in interpolation) and allocates data for the flux fields
  !! @param current_state The current model state_mod
  subroutine initialisation_callback(current_state)
    type(model_state_type), target, intent(inout) :: current_state
    logical :: l_qdiag

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
         q_advection(current_state%global_grid%size(Z_INDEX), current_state%number_q_fields), &
         tracer_advection(current_state%global_grid%size(Z_INDEX), current_state%n_tracers))

    advect_flow=determine_if_advection_here(options_get_string(current_state%options_database, "advection_flow_fields"))    
    advect_th=determine_if_advection_here(options_get_string(current_state%options_database, "advection_theta_field"))
    advect_q=determine_if_advection_here(options_get_string(current_state%options_database, "advection_q_fields"))
    advect_tracer=determine_if_advection_here(options_get_string(current_state%options_database, "advection_theta_field"))

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
      iqv=get_q_index(standard_q_names%VAPOUR, 'tvd_advection')
      allocate( tend_3d_qv(current_state%local_grid%size(Z_INDEX),  &
                           current_state%local_grid%size(Y_INDEX),  &
                           current_state%local_grid%size(X_INDEX)   )    )
    endif
    if (l_tend_3d_ql) then
      iql=get_q_index(standard_q_names%CLOUD_LIQUID_MASS, 'tvd_advection')
      allocate( tend_3d_ql(current_state%local_grid%size(Z_INDEX),  &
                           current_state%local_grid%size(Y_INDEX),  &
                           current_state%local_grid%size(X_INDEX)   )    )
    endif
    if (l_tend_3d_qi) then
      iqi=get_q_index(standard_q_names%ICE_MASS, 'tvd_advection')
      allocate( tend_3d_qi(current_state%local_grid%size(Z_INDEX),  &
                           current_state%local_grid%size(Y_INDEX),  &
                           current_state%local_grid%size(X_INDEX)   )    )
    endif
    if (l_tend_3d_qr) then
      iqr=get_q_index(standard_q_names%RAIN_MASS, 'tvd_advection')
      allocate( tend_3d_qr(current_state%local_grid%size(Z_INDEX),  &
                           current_state%local_grid%size(Y_INDEX),  &
                           current_state%local_grid%size(X_INDEX)   )    )
    endif
    if (l_tend_3d_qs) then
      iqs=get_q_index(standard_q_names%SNOW_MASS, 'tvd_advection')
      allocate( tend_3d_qs(current_state%local_grid%size(Z_INDEX),  &
                           current_state%local_grid%size(Y_INDEX),  &
                           current_state%local_grid%size(X_INDEX)   )    )
    endif
    if (l_tend_3d_qg) then
      iqg=get_q_index(standard_q_names%GRAUPEL_MASS, 'tvd_advection')
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
    if (allocated(tracer_advection)) deallocate(tracer_advection)

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

  !> Timestep callback hook which performs the TVD advection for each prognostic field
  !! @param current_state The current model state_mod
  subroutine timestep_callback(current_state)
    type(model_state_type), target, intent(inout) :: current_state
    integer :: current_x_index, current_y_index, target_x_index, target_y_index
    logical :: calculate_diagnostics

    calculate_diagnostics = current_state%diagnostic_sample_timestep &
                            .and. .not. current_state%halo_column

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


    if (current_state%halo_column) then
      if (.not. ((current_state%column_local_y == current_state%local_grid%halo_size(Y_INDEX) .and. &
           current_state%column_local_x .le. current_state%local_grid%local_domain_end_index(X_INDEX) .and. &
           current_state%column_local_x .ge. current_state%local_grid%local_domain_start_index(X_INDEX)-1) .or. &
           (current_state%column_local_x == current_state%local_grid%halo_size(X_INDEX) .and. &
           current_state%column_local_y .ge. current_state%local_grid%local_domain_start_index(Y_INDEX) &
           .and. current_state%column_local_y .le. current_state%local_grid%local_domain_end_index(Y_INDEX)) )) return
    end if

    if (calculate_diagnostics) &
        call save_precomponent_tendencies(current_state, current_x_index, current_y_index, target_x_index, target_y_index)

    if (advect_flow) call advect_flow_fields(current_state)
    if (advect_th) call advect_theta(current_state)
    if (advect_q) call advect_q_fields(current_state)
    if (advect_tracer) call advect_tracer_fields(current_state)

    if (calculate_diagnostics) &
        call compute_component_tendencies(current_state, current_x_index, current_y_index, target_x_index, target_y_index)

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
    if (is_component_enabled(current_state%options_database, "profile_diagnostics")) then
       ! NOTE: flux_z is declared at the top of module and then passed into ultflx, through argument 
       !       list in advect_scalar_field.
       tvd_dgs_terms%adv_u_dgs(:, current_state%column_local_y, current_state%column_local_x) =  &
            flux_z(:)
    endif
#endif

#ifdef V_ACTIVE
    call advect_v(current_state%column_local_y, current_state%column_local_x, dtm, current_state%u, &
         current_state%v, current_state%w, current_state%zv, current_state%sv, current_state%global_grid, &
         current_state%local_grid, current_state%parallel, current_state%halo_column, current_state%field_stepping)
    if (is_component_enabled(current_state%options_database, "profile_diagnostics")) then
       ! NOTE: flux_z is declared at the top of module and then passed into ultflx, through argument 
       !       list in advect_scalar_field.
       tvd_dgs_terms%adv_v_dgs(:, current_state%column_local_y, current_state%column_local_x) =  &
            flux_z(:)
    endif
#endif

#ifdef W_ACTIVE
    call advect_w(current_state%column_local_y, current_state%column_local_x, dtm, current_state%u, &
         current_state%v, current_state%w, current_state%zw, current_state%sw, current_state%global_grid,&
         current_state%local_grid, current_state%parallel, current_state%halo_column, current_state%field_stepping)
    if (is_component_enabled(current_state%options_database, "profile_diagnostics")) then
       ! NOTE: flux_z is declared at the top of module and then passed into ultflx, through argument 
       !       list in advect_scalar_field.
       tvd_dgs_terms%adv_w_dgs(:, current_state%column_local_y, current_state%column_local_x) =  &
            flux_z(:)
    endif
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
        if (is_component_enabled(current_state%options_database, "profile_diagnostics")) then
           ! NOTE: flux_z is declared at the top of module and then passed into ultflx, through argument 
           !       list in advect_scalar_field.
           tvd_dgs_terms%adv_q_dgs(:, current_state%column_local_y, current_state%column_local_x, i) =  &
                flux_z(:)
        endif
           
      end if
    end do
  end subroutine advect_q_fields

  !> Advects the tracer fields
  !! @param current_state The current model state_mod
  subroutine advect_tracer_fields(current_state)
    type(model_state_type), intent(inout) :: current_state

    integer :: i
    real(kind=DEFAULT_PRECISION) :: dtm

    dtm = current_state%dtm*2.0_DEFAULT_PRECISION
    if (current_state%scalar_stepping == FORWARD_STEPPING) dtm=current_state%dtm

    do i=1,current_state%n_tracers
      call advect_scalar_field(current_state%column_local_y, current_state%column_local_x, dtm, current_state%u, &
           current_state%v, current_state%w, current_state%ztracer(i), current_state%tracer(i), current_state%stracer(i), &
           current_state%global_grid, current_state%local_grid, current_state%parallel, &
           current_state%halo_column, current_state%field_stepping)
      tracer_advection(:,i)=current_state%stracer(i)%data(:, current_state%column_local_y, current_state%column_local_x)          
    end do
  end subroutine advect_tracer_fields

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
      if (is_component_enabled(current_state%options_database, "profile_diagnostics")) then
           ! NOTE: flux_z is declared at the top of module and then passed into ultflx, through argument 
           !       list in advect_scalar_field.
         tvd_dgs_terms%adv_th_dgs(:, current_state%column_local_y, current_state%column_local_x) =  &
              flux_z(:)
        endif
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


  !> Sets the published field value from the temporary diagnostic values held by
  !this component.
  !! @param field_value Populated with the value of the field
  !! @param real_1d_field Optional one dimensional real of values to publish
  !! @param real_2d_field Optional two dimensional real of values to publish
  subroutine set_published_field_value(field_value, real_1d_field, real_2d_field, real_3d_field)
    type(component_field_value_type), intent(inout) :: field_value
    real(kind=DEFAULT_PRECISION), dimension(:), optional :: real_1d_field
    real(kind=DEFAULT_PRECISION), dimension(:,:), optional :: real_2d_field
    real(kind=DEFAULT_PRECISION), dimension(:,:,:), optional :: real_3d_field


    if (present(real_1d_field)) then
      allocate(field_value%real_1d_array(size(real_1d_field)),source=real_1d_field)
    else if (present(real_2d_field)) then
      allocate(field_value%real_2d_array(size(real_2d_field, 1),size(real_2d_field, 2)), source=real_2d_field)
    else if (present(real_3d_field)) then
      allocate(field_value%real_3d_array(size(real_3d_field, 1),size(real_3d_field, 2), size(real_3d_field, 3)), &
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
