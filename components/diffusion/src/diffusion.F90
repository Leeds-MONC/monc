!> Diffusion on the TH and Q fields
module diffusion_mod
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
  use q_indices_mod, only: get_q_index, standard_q_names

  implicit none

#ifndef TEST_MODE
  private
#endif

  real(kind=DEFAULT_PRECISION), dimension(:), allocatable :: th_diffusion
  real(kind=DEFAULT_PRECISION), dimension(:,:), allocatable :: q_diffusion
  real(kind=DEFAULT_PRECISION), dimension(:,:), allocatable :: tracer_diffusion

  ! Local tendency diagnostic variables for this component
  ! 3D tendency fields and logicals for their use
  real(kind=DEFAULT_PRECISION), dimension(:,:,:), allocatable ::     &
       tend_3d_th,tend_3d_qv,tend_3d_ql,tend_3d_qi,tend_3d_qr,tend_3d_qs,tend_3d_qg,tend_3d_tabs
  logical :: l_tend_3d_th,l_tend_3d_qv,l_tend_3d_ql,l_tend_3d_qi,l_tend_3d_qr,l_tend_3d_qs,l_tend_3d_qg,l_tend_3d_tabs
  ! Local mean tendency profile fields and logicals for their use
  real(kind=DEFAULT_PRECISION), dimension(:), allocatable ::         &
       tend_pr_tot_th,tend_pr_tot_qv,tend_pr_tot_ql,tend_pr_tot_qi,tend_pr_tot_qr,tend_pr_tot_qs,tend_pr_tot_qg,  &
       tend_pr_tot_tabs
  logical :: l_tend_pr_tot_th,l_tend_pr_tot_qv,l_tend_pr_tot_ql,l_tend_pr_tot_qi,l_tend_pr_tot_qr,  &
             l_tend_pr_tot_qs,l_tend_pr_tot_qg,l_tend_pr_tot_tabs
  ! q indices
  integer :: iqv=0, iql=0, iqr=0, iqi=0, iqs=0, iqg=0

  public diffusion_get_descriptor

contains

  !> Provides the descriptor back to the caller and is used in component registration
  !! @returns The termination check component descriptor
  type(component_descriptor_type) function diffusion_get_descriptor()
    diffusion_get_descriptor%name="diffusion"
    diffusion_get_descriptor%version=0.1
    diffusion_get_descriptor%initialisation=>initialisation_callback
    diffusion_get_descriptor%timestep=>timestep_callback
    diffusion_get_descriptor%finalisation=>finalisation_callback

    diffusion_get_descriptor%field_value_retrieval=>field_value_retrieval_callback
    diffusion_get_descriptor%field_information_retrieval=>field_information_retrieval_callback
    allocate(diffusion_get_descriptor%published_fields(2+8+8+1))

    diffusion_get_descriptor%published_fields(1)="th_diffusion"
    diffusion_get_descriptor%published_fields(2)="q_diffusion"

    diffusion_get_descriptor%published_fields(2+1)="tend_th_diffusion_3d_local"
    diffusion_get_descriptor%published_fields(2+2)="tend_qv_diffusion_3d_local"
    diffusion_get_descriptor%published_fields(2+3)="tend_ql_diffusion_3d_local"
    diffusion_get_descriptor%published_fields(2+4)="tend_qi_diffusion_3d_local"
    diffusion_get_descriptor%published_fields(2+5)="tend_qr_diffusion_3d_local"
    diffusion_get_descriptor%published_fields(2+6)="tend_qs_diffusion_3d_local"
    diffusion_get_descriptor%published_fields(2+7)="tend_qg_diffusion_3d_local"
    diffusion_get_descriptor%published_fields(2+8)="tend_tabs_diffusion_3d_local"

    diffusion_get_descriptor%published_fields(2+8+1)="tend_th_diffusion_profile_total_local"
    diffusion_get_descriptor%published_fields(2+8+2)="tend_qv_diffusion_profile_total_local"
    diffusion_get_descriptor%published_fields(2+8+3)="tend_ql_diffusion_profile_total_local"
    diffusion_get_descriptor%published_fields(2+8+4)="tend_qi_diffusion_profile_total_local"
    diffusion_get_descriptor%published_fields(2+8+5)="tend_qr_diffusion_profile_total_local"
    diffusion_get_descriptor%published_fields(2+8+6)="tend_qs_diffusion_profile_total_local"
    diffusion_get_descriptor%published_fields(2+8+7)="tend_qg_diffusion_profile_total_local"
    diffusion_get_descriptor%published_fields(2+8+8)="tend_tabs_diffusion_profile_total_local"

    diffusion_get_descriptor%published_fields(2+8+8+1)="tracer_diffusion"

  end function diffusion_get_descriptor


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
    if (name .eq. "q_diffusion") then
      field_information%number_dimensions=2
    else if (name .eq. "tracer_diffusion") then
      field_information%number_dimensions=2
    else
      field_information%number_dimensions=1
    end if
    field_information%dimension_sizes(1)=current_state%local_grid%size(Z_INDEX)
    if (name .eq. "q_diffusion") field_information%dimension_sizes(2)=current_state%number_q_fields
    if (name .eq. "tracer_diffusion") field_information%dimension_sizes(2)=current_state%n_tracers
    field_information%enabled=.true.

    ! Field information for 3d
    strcomp=INDEX(name, "diffusion_3d_local")
    if (strcomp .ne. 0) then
      field_information%field_type=COMPONENT_ARRAY_FIELD_TYPE
      field_information%number_dimensions=3
      field_information%dimension_sizes(1)=current_state%local_grid%size(Z_INDEX)
      field_information%dimension_sizes(2)=current_state%local_grid%size(Y_INDEX)
      field_information%dimension_sizes(3)=current_state%local_grid%size(X_INDEX)
      field_information%data_type=COMPONENT_DOUBLE_DATA_TYPE

      if      (name .eq. "tend_th_diffusion_3d_local") then
        field_information%enabled=l_tend_3d_th
      else if (name .eq. "tend_qv_diffusion_3d_local") then
        field_information%enabled=l_tend_3d_qv
      else if (name .eq. "tend_ql_diffusion_3d_local") then
        field_information%enabled=l_tend_3d_ql
      else if (name .eq. "tend_qi_diffusion_3d_local") then
        field_information%enabled=l_tend_3d_qi
      else if (name .eq. "tend_qr_diffusion_3d_local") then
        field_information%enabled=l_tend_3d_qr
      else if (name .eq. "tend_qs_diffusion_3d_local") then
        field_information%enabled=l_tend_3d_qs
      else if (name .eq. "tend_qg_diffusion_3d_local") then
        field_information%enabled=l_tend_3d_qg
      else if (name .eq. "tend_tabs_diffusion_3d_local") then
        field_information%enabled=l_tend_3d_tabs
      else
        field_information%enabled=.true.
      end if

    end if !end 3d check

    ! Field information for profiles
    strcomp=INDEX(name, "diffusion_profile_total_local")
    if (strcomp .ne. 0) then
      field_information%field_type=COMPONENT_ARRAY_FIELD_TYPE
      field_information%number_dimensions=1
      field_information%dimension_sizes(1)=current_state%local_grid%size(Z_INDEX)
      field_information%data_type=COMPONENT_DOUBLE_DATA_TYPE

      if      (name .eq. "tend_th_diffusion_profile_total_local") then
        field_information%enabled=l_tend_pr_tot_th
      else if (name .eq. "tend_qv_diffusion_profile_total_local") then
        field_information%enabled=l_tend_pr_tot_qv
      else if (name .eq. "tend_ql_diffusion_profile_total_local") then
        field_information%enabled=l_tend_pr_tot_ql
      else if (name .eq. "tend_qi_diffusion_profile_total_local") then
        field_information%enabled=l_tend_pr_tot_qi
      else if (name .eq. "tend_qr_diffusion_profile_total_local") then
        field_information%enabled=l_tend_pr_tot_qr
      else if (name .eq. "tend_qs_diffusion_profile_total_local") then
        field_information%enabled=l_tend_pr_tot_qs
      else if (name .eq. "tend_qg_diffusion_profile_total_local") then
        field_information%enabled=l_tend_pr_tot_qg
      else if (name .eq. "tend_tabs_diffusion_profile_total_local") then
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
    
    if (name .eq. "th_diffusion") then
      call set_published_field_value(field_value, real_1d_field=th_diffusion)
    else if (name .eq. "q_diffusion") then
      call set_published_field_value(field_value, real_2d_field=q_diffusion)
    else if (name .eq. "tracer_diffusion") then
      call set_published_field_value(field_value, real_2d_field=tracer_diffusion)

    ! 3d Tendency Fields
    else if (name .eq. "tend_th_diffusion_3d_local" .and. allocated(tend_3d_th)) then
      call set_published_field_value(field_value, real_3d_field=tend_3d_th)
    else if (name .eq. "tend_qv_diffusion_3d_local" .and. allocated(tend_3d_qv)) then
      call set_published_field_value(field_value, real_3d_field=tend_3d_qv)
    else if (name .eq. "tend_ql_diffusion_3d_local" .and. allocated(tend_3d_ql)) then
      call set_published_field_value(field_value, real_3d_field=tend_3d_ql)
    else if (name .eq. "tend_qi_diffusion_3d_local" .and. allocated(tend_3d_qi)) then
      call set_published_field_value(field_value, real_3d_field=tend_3d_qi)
    else if (name .eq. "tend_qr_diffusion_3d_local" .and. allocated(tend_3d_qr)) then
      call set_published_field_value(field_value, real_3d_field=tend_3d_qr)
    else if (name .eq. "tend_qs_diffusion_3d_local" .and. allocated(tend_3d_qs)) then
      call set_published_field_value(field_value, real_3d_field=tend_3d_qs)
    else if (name .eq. "tend_qg_diffusion_3d_local" .and. allocated(tend_3d_qg)) then
      call set_published_field_value(field_value, real_3d_field=tend_3d_qg)
    else if (name .eq. "tend_tabs_diffusion_3d_local" .and. allocated(tend_3d_tabs)) then
      call set_published_field_value(field_value, real_3d_field=tend_3d_tabs)

    ! Profile Tendency Fields
    else if (name .eq. "tend_th_diffusion_profile_total_local" .and. allocated(tend_pr_tot_th)) then
      call set_published_field_value(field_value, real_1d_field=tend_pr_tot_th)
    else if (name .eq. "tend_qv_diffusion_profile_total_local" .and. allocated(tend_pr_tot_qv)) then
      call set_published_field_value(field_value, real_1d_field=tend_pr_tot_qv)
    else if (name .eq. "tend_ql_diffusion_profile_total_local" .and. allocated(tend_pr_tot_ql)) then
      call set_published_field_value(field_value, real_1d_field=tend_pr_tot_ql)
    else if (name .eq. "tend_qi_diffusion_profile_total_local" .and. allocated(tend_pr_tot_qi)) then
      call set_published_field_value(field_value, real_1d_field=tend_pr_tot_qi)
    else if (name .eq. "tend_qr_diffusion_profile_total_local" .and. allocated(tend_pr_tot_qr)) then
      call set_published_field_value(field_value, real_1d_field=tend_pr_tot_qr)
    else if (name .eq. "tend_qs_diffusion_profile_total_local" .and. allocated(tend_pr_tot_qs)) then
      call set_published_field_value(field_value, real_1d_field=tend_pr_tot_qs)
    else if (name .eq. "tend_qg_diffusion_profile_total_local" .and. allocated(tend_pr_tot_qg)) then
      call set_published_field_value(field_value, real_1d_field=tend_pr_tot_qg)
    else if (name .eq. "tend_tabs_diffusion_profile_total_local" .and. allocated(tend_pr_tot_tabs)) then
      call set_published_field_value(field_value, real_1d_field=tend_pr_tot_tabs)
    end if

  end subroutine field_value_retrieval_callback

  !> Sets up the stencil_mod (used in interpolation) and allocates data for the flux fields
  !! @param current_state The current model state_mod
  subroutine initialisation_callback(current_state)
    type(model_state_type), target, intent(inout) :: current_state

    integer :: z_size, y_size, x_size
    logical :: l_qdiag

    z_size=current_state%local_grid%size(Z_INDEX) + current_state%local_grid%halo_size(Z_INDEX) * 2
    y_size=current_state%local_grid%size(Y_INDEX) + current_state%local_grid%halo_size(Y_INDEX) * 2
    x_size=current_state%local_grid%size(X_INDEX) + current_state%local_grid%halo_size(X_INDEX) * 2
    allocate(current_state%diff_coefficient%data(z_size, y_size, x_size))

    z_size=current_state%global_grid%size(Z_INDEX)
    allocate(th_diffusion(z_size))
    allocate(q_diffusion(z_size, current_state%number_q_fields))
    if (current_state%n_tracers > 0) allocate(tracer_diffusion(z_size, current_state%n_tracers))

    ! Set tendency diagnostic logicals based on availability
    ! Need to use 3d tendencies to compute the profiles, so they will be allocated
    !      in the case where profiles are available
    l_qdiag =  (.not. current_state%passive_q .and. current_state%number_q_fields .gt. 0)

    l_tend_pr_tot_th  = current_state%th%active 
    l_tend_pr_tot_qv  = l_qdiag .and. current_state%number_q_fields .ge. 1
    l_tend_pr_tot_ql  = l_qdiag .and. current_state%number_q_fields .ge. 2
    l_tend_pr_tot_qi  = l_qdiag .and. current_state%number_q_fields .ge. 11
    l_tend_pr_tot_qr  = l_qdiag .and. current_state%number_q_fields .ge. 11
    l_tend_pr_tot_qs  = l_qdiag .and. current_state%number_q_fields .ge. 11
    l_tend_pr_tot_qg  = l_qdiag .and. current_state%number_q_fields .ge. 11
    l_tend_pr_tot_tabs  = l_tend_pr_tot_th

    l_tend_3d_th  = current_state%th%active .or. l_tend_pr_tot_th
    l_tend_3d_qv  = (l_qdiag .and. current_state%number_q_fields .ge. 1) .or. l_tend_pr_tot_qv
    l_tend_3d_ql  = (l_qdiag .and. current_state%number_q_fields .ge. 2) .or. l_tend_pr_tot_ql
    l_tend_3d_qi  = (l_qdiag .and. current_state%number_q_fields .ge. 11) .or. l_tend_pr_tot_qi
    l_tend_3d_qr  = (l_qdiag .and. current_state%number_q_fields .ge. 11) .or. l_tend_pr_tot_qr
    l_tend_3d_qs  = (l_qdiag .and. current_state%number_q_fields .ge. 11) .or. l_tend_pr_tot_qs
    l_tend_3d_qg  = (l_qdiag .and. current_state%number_q_fields .ge. 11) .or. l_tend_pr_tot_qg
    l_tend_3d_tabs = l_tend_3d_th

    ! Allocate 3d tendency fields upon availability
    if (l_tend_3d_th) then
      allocate( tend_3d_th(current_state%local_grid%size(Z_INDEX),  &
                           current_state%local_grid%size(Y_INDEX),  &
                           current_state%local_grid%size(X_INDEX)   )    )
    endif
    if (l_tend_3d_qv) then
      iqv=get_q_index(standard_q_names%VAPOUR, 'diffusion')
      allocate( tend_3d_qv(current_state%local_grid%size(Z_INDEX),  &
                           current_state%local_grid%size(Y_INDEX),  &
                           current_state%local_grid%size(X_INDEX)   )    )
    endif
    if (l_tend_3d_ql) then
      iql=get_q_index(standard_q_names%CLOUD_LIQUID_MASS, 'diffusion')
      allocate( tend_3d_ql(current_state%local_grid%size(Z_INDEX),  &
                           current_state%local_grid%size(Y_INDEX),  &
                           current_state%local_grid%size(X_INDEX)   )    )
    endif
    if (l_tend_3d_qi) then
      iqi=get_q_index(standard_q_names%ICE_MASS, 'diffusion')
      allocate( tend_3d_qi(current_state%local_grid%size(Z_INDEX),  &
                           current_state%local_grid%size(Y_INDEX),  &
                           current_state%local_grid%size(X_INDEX)   )    )
    endif
    if (l_tend_3d_qr) then
      iqr=get_q_index(standard_q_names%RAIN_MASS, 'diffusion')
      allocate( tend_3d_qr(current_state%local_grid%size(Z_INDEX),  &
                           current_state%local_grid%size(Y_INDEX),  &
                           current_state%local_grid%size(X_INDEX)   )    )
    endif
    if (l_tend_3d_qs) then
      iqs=get_q_index(standard_q_names%SNOW_MASS, 'diffusion')
      allocate( tend_3d_qs(current_state%local_grid%size(Z_INDEX),  &
                           current_state%local_grid%size(Y_INDEX),  &
                           current_state%local_grid%size(X_INDEX)   )    )
    endif
    if (l_tend_3d_qg) then
      iqg=get_q_index(standard_q_names%GRAUPEL_MASS, 'diffusion')
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


  subroutine finalisation_callback(current_state)
    type(model_state_type), target, intent(inout) :: current_state

    if (allocated(th_diffusion)) deallocate(th_diffusion)
    if (allocated(q_diffusion)) deallocate(q_diffusion)
    if (allocated(tracer_diffusion)) deallocate(tracer_diffusion)

    if (allocated(tend_3d_th)) deallocate(tend_3d_th)
    if (allocated(tend_3d_qv)) deallocate(tend_3d_qv)
    if (allocated(tend_3d_ql)) deallocate(tend_3d_ql)
    if (allocated(tend_3d_qi)) deallocate(tend_3d_qi)
    if (allocated(tend_3d_qr)) deallocate(tend_3d_qr)
    if (allocated(tend_3d_qs)) deallocate(tend_3d_qs)
    if (allocated(tend_3d_qg)) deallocate(tend_3d_qg)
    if (allocated(tend_3d_tabs)) deallocate(tend_3d_tabs)

    if (allocated(tend_pr_tot_th)) deallocate(tend_pr_tot_th)
    if (allocated(tend_pr_tot_qv)) deallocate(tend_pr_tot_qv)
    if (allocated(tend_pr_tot_ql)) deallocate(tend_pr_tot_ql)
    if (allocated(tend_pr_tot_qi)) deallocate(tend_pr_tot_qi)
    if (allocated(tend_pr_tot_qr)) deallocate(tend_pr_tot_qr)
    if (allocated(tend_pr_tot_qs)) deallocate(tend_pr_tot_qs)
    if (allocated(tend_pr_tot_qg)) deallocate(tend_pr_tot_qg)
    if (allocated(tend_pr_tot_tabs)) deallocate(tend_pr_tot_tabs)

  end subroutine finalisation_callback

  !> At each timestep will compute the diffusion source terms for TH and Q fields per column if these fields are active
  !! @param current_state The current model state
  subroutine timestep_callback(current_state)
    type(model_state_type), target, intent(inout) :: current_state

    integer :: local_y, local_x, target_x_index, target_y_index

    logical :: calculate_diagnostics

    calculate_diagnostics = current_state%diagnostic_sample_timestep

    local_y=current_state%column_local_y
    local_x=current_state%column_local_x
    target_y_index=local_y-current_state%local_grid%halo_size(Y_INDEX)
    target_x_index=local_x-current_state%local_grid%halo_size(X_INDEX)

    ! Zero profile tendency totals on first instance in the sum
    if (current_state%first_timestep_column) then
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

    if (.not. current_state%use_viscosity_and_diffusion .or. current_state%halo_column) return

    if (current_state%diffusion_halo_swap_state%swap_in_progress) then
      ! If there is a diffusion halo swap in progress then complete it
      call complete_nonblocking_halo_swap(current_state, current_state%diffusion_halo_swap_state, &
           perform_local_data_copy_for_diff, copy_halo_buffer_to_diff, copy_halo_buffer_to_diff_corners)
    end if

    if (calculate_diagnostics) call save_precomponent_tendencies(current_state, local_x, local_y, target_x_index, target_y_index)

    if (current_state%th%active) call perform_th_diffusion(current_state, local_y, local_x)
    if (current_state%number_q_fields .gt. 0) call perform_q_diffusion(current_state, local_y, local_x)
    if (current_state%n_radioactive_tracers .gt. 0) call perform_tracer_diffusion(current_state, local_y, local_x)

    if (calculate_diagnostics) call compute_component_tendencies(current_state, local_x, local_y, target_x_index, target_y_index)

  end subroutine timestep_callback

  !> Computes the diffusion source terms for each tracer field
  !! @param current_state The current model state
  !! @param local_y Local Y index
  !! @param local_x Local X index
  subroutine perform_tracer_diffusion(current_state, local_y, local_x)
    type(model_state_type), target, intent(inout) :: current_state
    integer, intent(in) :: local_y, local_x

    integer :: n

    do n=current_state%radioactive_tracer_index, current_state%radioactive_tracer_index + current_state%n_radioactive_tracers - 1 
      call general_diffusion(current_state, current_state%ztracer(n), current_state%stracer(n), local_y, local_x, &
        tracer_diffusion(:,n))
    end do
  end subroutine perform_tracer_diffusion  


  !> Computes the diffusion source terms for each Q field
  !! @param current_state The current model state
  !! @param local_y Local Y index
  !! @param local_x Local X index
  subroutine perform_q_diffusion(current_state, local_y, local_x)
    type(model_state_type), target, intent(inout) :: current_state
    integer, intent(in) :: local_y, local_x

    integer :: n

    do n=1, current_state%number_q_fields
      call general_diffusion(current_state, current_state%zq(n), current_state%sq(n), local_y, local_x, q_diffusion(:,n))
    end do
  end subroutine perform_q_diffusion  

  !> Computes the diffusion source terms for the theta field
  !! @param current_state The current model state
  !! @param local_y Local Y index
  !! @param local_x Local X index
  subroutine perform_th_diffusion(current_state, local_y, local_x)
    type(model_state_type), target, intent(inout) :: current_state
    integer, intent(in) :: local_y, local_x

    integer :: k
    real(kind=DEFAULT_PRECISION) :: adjustment

    call general_diffusion(current_state, current_state%zth, current_state%sth, local_y, local_x, th_diffusion)

    if (current_state%use_anelastic_equations) then
      ! This code only needs to be executed if anelastic, otherwise THREF is constant and the increment calculated here is zero
      do k=2, current_state%local_grid%size(Z_INDEX)
        adjustment=(current_state%global_grid%configuration%vertical%cza(k)*&
             current_state%global_grid%configuration%vertical%dthref(k)*&
             current_state%diff_coefficient%data(k, local_y, local_x) - current_state%global_grid%configuration%vertical%czb(k)*&
             current_state%global_grid%configuration%vertical%dthref(k-1)*&
             current_state%diff_coefficient%data(k-1, local_y, local_x))
        current_state%sth%data(k, local_y, local_x)=current_state%sth%data(k, local_y, local_x)+adjustment
        th_diffusion(k)=th_diffusion(k)+adjustment
      end do
    end if
  end subroutine perform_th_diffusion

  !> General diffusion computation for any field which is provided as arguments. Works in a column
  !! @param current_state The current model state
  !! @param field The field to take values from, typically zth or zq(n)
  !! @param source_field The source target field to update, typically sth or sq(n)
  !! @param local_y Local Y index
  !! @param local_x Local X index
  subroutine general_diffusion(current_state, field, source_field, local_y, local_x, diagnostics)
    type(model_state_type), target, intent(inout) :: current_state
    type(prognostic_field_type), intent(inout) :: field, source_field
    integer, intent(in) :: local_y, local_x
    real(kind=DEFAULT_PRECISION), dimension(:), intent(inout), optional :: diagnostics

    real(kind=DEFAULT_PRECISION) :: term
    integer :: k

    do k=2, current_state%local_grid%size(Z_INDEX)
      term=current_state%global_grid%configuration%horizontal%cx2*0.25_DEFAULT_PRECISION*&
           (((current_state%diff_coefficient%data(k, local_y, local_x)+&
           current_state%diff_coefficient%data(k, local_y, local_x-1))+&
           (current_state%diff_coefficient%data(k-1, local_y, local_x)+&
           current_state%diff_coefficient%data(k-1, local_y, local_x-1)))&
           *(field%data(k, local_y, local_x-1)-field%data(k, local_y, local_x)) -&
           ((current_state%diff_coefficient%data(k, local_y, local_x+1)+&
           current_state%diff_coefficient%data(k, local_y, local_x))+&
           (current_state%diff_coefficient%data(k-1, local_y, local_x+1)+&
           current_state%diff_coefficient%data(k-1, local_y, local_x)))&
           *(field%data(k, local_y, local_x)-field%data(k, local_y, local_x+1)) )&
           +current_state%global_grid%configuration%horizontal%cy2*0.25_DEFAULT_PRECISION*(&
           ((current_state%diff_coefficient%data(k, local_y, local_x)+&
           current_state%diff_coefficient%data(k, local_y-1, local_x))+&
           (current_state%diff_coefficient%data(k-1, local_y, local_x)+&
           current_state%diff_coefficient%data(k-1, local_y-1, local_x)))&
           *(field%data(k, local_y-1, local_x)-field%data(k, local_y, local_x)) -&
           ((current_state%diff_coefficient%data(k, local_y+1, local_x)+&
           current_state%diff_coefficient%data(k, local_y, local_x))+&
           (current_state%diff_coefficient%data(k-1, local_y+1, local_x)+&
           current_state%diff_coefficient%data(k-1, local_y, local_x)))&
           *(field%data(k, local_y, local_x)-field%data(k, local_y+1, local_x)) )&
           +( current_state%global_grid%configuration%vertical%czb(k)*&
           current_state%diff_coefficient%data(k-1, local_y, local_x)*&
           (field%data(k-1, local_y, local_x)-field%data(k, local_y, local_x)))

      if (k .lt. current_state%local_grid%size(Z_INDEX)) then
        term=term - current_state%global_grid%configuration%vertical%cza(k)*&
             current_state%diff_coefficient%data(k, local_y, local_x)*&
             (field%data(k, local_y, local_x)-field%data(k+1, local_y, local_x))
      end if
      source_field%data(k, local_y, local_x)=source_field%data(k, local_y, local_x)+term
      if (present(diagnostics)) diagnostics(k)=term
    end do
  end subroutine general_diffusion

    !> Does local data copying for diffusion coefficient variable halo swap
  !! @param current_state The current model state_mod
  !! @param source_data Optional source data which is written into
  subroutine perform_local_data_copy_for_diff(current_state, halo_depth, involve_corners, source_data)
    type(model_state_type), intent(inout) :: current_state
    integer, intent(in) :: halo_depth
    logical, intent(in) :: involve_corners
    type(field_data_wrapper_type), dimension(:), intent(in), optional :: source_data

    call perform_local_data_copy_for_field(current_state%diff_coefficient%data, current_state%local_grid, &
         current_state%parallel%my_rank, halo_depth, involve_corners)
  end subroutine perform_local_data_copy_for_diff

  !> Copies the halo buffer to halo location for the diffusion coefficient field
  !! @param current_state The current model state
  !! @param neighbour_description The halo swapping description of the neighbour we are accessing the buffer of
  !! @param dim The dimension we receive for
  !! @param target_index The target index for the dimension we are receiving for
  !! @param neighbour_location The location in the local neighbour data stores of this neighbour
  !! @param current_page The current, next, halo swap page to read from (all previous have been read and copied already)
  !! @param source_data Optional source data which is written into
  subroutine copy_halo_buffer_to_diff(current_state, neighbour_description, dim, target_index, &
       neighbour_location, current_page, source_data)
    type(model_state_type), intent(inout) :: current_state
    integer, intent(in) :: dim, target_index, neighbour_location
    integer, intent(inout) :: current_page(:)
    type(neighbour_description_type), intent(inout) :: neighbour_description
    type(field_data_wrapper_type), dimension(:), intent(in), optional :: source_data

    call copy_buffer_to_field(current_state%local_grid, neighbour_description%recv_halo_buffer, &
         current_state%diff_coefficient%data, dim, target_index, current_page(neighbour_location))

    current_page(neighbour_location)=current_page(neighbour_location)+1
  end subroutine copy_halo_buffer_to_diff

  !> Copies the corner halo buffer to the diffusion coefficient field corners
  !! @param current_state The current model state
  !! @param neighbour_description The halo swapping description of the neighbour we are accessing the buffer of
  !! @param corner_loc The corner location
  !! @param x_target_index The X target index for the dimension we are receiving for
  !! @param y_target_index The Y target index for the dimension we are receiving for
  !! @param neighbour_location The location in the local neighbour data stores of this neighbour
  !! @param current_page The current, next, halo swap page to read from (all previous have been read and copied already)
  !! @param source_data Optional source data which is written into
  subroutine copy_halo_buffer_to_diff_corners(current_state, neighbour_description, corner_loc, x_target_index, &
       y_target_index, neighbour_location, current_page, source_data)
    type(model_state_type), intent(inout) :: current_state
    integer, intent(in) :: corner_loc, x_target_index, y_target_index, neighbour_location
    integer, intent(inout) :: current_page(:)
    type(neighbour_description_type), intent(inout) :: neighbour_description
    type(field_data_wrapper_type), dimension(:), intent(in), optional :: source_data

    call copy_buffer_to_corner(current_state%local_grid, neighbour_description%recv_corner_buffer, &
         current_state%diff_coefficient%data, corner_loc, x_target_index, y_target_index, current_page(neighbour_location))

    current_page(neighbour_location)=current_page(neighbour_location)+1
  end subroutine copy_halo_buffer_to_diff_corners  

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


end module diffusion_mod
