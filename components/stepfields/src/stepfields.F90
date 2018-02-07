!> Does the field stepping
!! Stepping is called at the end of processing a column and steps the x-2 column
module stepfields_mod
  use monc_component_mod, only : COMPONENT_ARRAY_FIELD_TYPE, &
      COMPONENT_DOUBLE_DATA_TYPE, component_descriptor_type, &
      component_field_value_type, component_field_information_type
  use optionsdatabase_mod, only : options_get_integer
  use state_mod, only : FORWARD_STEPPING, CENTRED_STEPPING, model_state_type
  use prognostics_mod, only : prognostic_field_type
  use grids_mod, only : local_grid_type, X_INDEX, Y_INDEX, Z_INDEX
  use datadefn_mod, only : DEFAULT_PRECISION
  use registry_mod, only : is_component_enabled
  use science_constants_mod, only : rlargep
  use q_indices_mod, only: get_q_index, standard_q_names

  implicit none

#ifndef TEST_MODE
  private
#endif

  logical :: determine_flow_minmax=.false., cfl_is_enabled

 
  real(kind=DEFAULT_PRECISION), allocatable :: resetq_min(:)
  logical :: l_nonconservative_positive_q=.true.

  ! Local tendency diagnostic variables for this component
  ! 3D tendency fields and logicals for their use
  real(kind=DEFAULT_PRECISION), dimension(:,:,:), allocatable ::     &
       tend_3d_th,tend_3d_qv,       &
       tend_3d_ql,tend_3d_qi,tend_3d_qr,tend_3d_qs,tend_3d_qg,       &
       tend_3d_tabs
  logical :: l_tend_3d_th,l_tend_3d_qv,       &
             l_tend_3d_ql,l_tend_3d_qi,l_tend_3d_qr,l_tend_3d_qs,l_tend_3d_qg,       &
             l_tend_3d_tabs
  ! Local mean tendency profile fields and logicals for their use
  real(kind=DEFAULT_PRECISION), dimension(:), allocatable ::         &
       tend_pr_tot_th,tend_pr_tot_qv,       &
       tend_pr_tot_ql,tend_pr_tot_qi,tend_pr_tot_qr,tend_pr_tot_qs,tend_pr_tot_qg,       &
       tend_pr_tot_tabs
  logical :: l_tend_pr_tot_th,l_tend_pr_tot_qv,       &
             l_tend_pr_tot_ql,l_tend_pr_tot_qi,l_tend_pr_tot_qr,l_tend_pr_tot_qs,l_tend_pr_tot_qg,       &
             l_tend_pr_tot_tabs
  ! q indices
  integer :: iqv=0, iql=0, iqr=0, iqi=0, iqs=0, iqg=0
  integer :: diagnostic_generation_frequency


  public stepfields_get_descriptor

contains

  !> Provides the descriptor back to the caller and is used in component registration
  !! @returns The KidReader component descriptor
  type(component_descriptor_type) function stepfields_get_descriptor()
    stepfields_get_descriptor%name="stepfields"
    stepfields_get_descriptor%version=0.1
    stepfields_get_descriptor%initialisation=>initialisation_callback
    stepfields_get_descriptor%timestep=>timestep_callback
    stepfields_get_descriptor%finalisation=>finalisation_callback

    stepfields_get_descriptor%field_value_retrieval=>field_value_retrieval_callback
    stepfields_get_descriptor%field_information_retrieval=>field_information_retrieval_callback
    allocate(stepfields_get_descriptor%published_fields(8+8))

    stepfields_get_descriptor%published_fields(1)="tend_th_total_3d_local"
    stepfields_get_descriptor%published_fields(2)="tend_qv_total_3d_local"
    stepfields_get_descriptor%published_fields(3)="tend_ql_total_3d_local"
    stepfields_get_descriptor%published_fields(4)="tend_qi_total_3d_local"
    stepfields_get_descriptor%published_fields(5)="tend_qr_total_3d_local"
    stepfields_get_descriptor%published_fields(6)="tend_qs_total_3d_local"
    stepfields_get_descriptor%published_fields(7)="tend_qg_total_3d_local"
    stepfields_get_descriptor%published_fields(8)="tend_tabs_total_3d_local"

    stepfields_get_descriptor%published_fields(8+1)="tend_th_total_profile_total_local"
    stepfields_get_descriptor%published_fields(8+2)="tend_qv_total_profile_total_local"
    stepfields_get_descriptor%published_fields(8+3)="tend_ql_total_profile_total_local"
    stepfields_get_descriptor%published_fields(8+4)="tend_qi_total_profile_total_local"
    stepfields_get_descriptor%published_fields(8+5)="tend_qr_total_profile_total_local"
    stepfields_get_descriptor%published_fields(8+6)="tend_qs_total_profile_total_local"
    stepfields_get_descriptor%published_fields(8+7)="tend_qg_total_profile_total_local"
    stepfields_get_descriptor%published_fields(8+8)="tend_tabs_total_profile_total_local"

  end function stepfields_get_descriptor

  !> Initialisation callback
  !! @param current_state The current model state
  subroutine initialisation_callback(current_state)
    type(model_state_type), target, intent(inout) :: current_state
    logical :: l_qdiag

    allocate(resetq_min(current_state%number_q_fields))
    cfl_is_enabled=is_component_enabled(current_state%options_database, "cfltest") 
    if (cfl_is_enabled) call reset_local_minmax_values(current_state)

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
      iqv=get_q_index(standard_q_names%VAPOUR, 'stepfields')
      allocate( tend_3d_qv(current_state%local_grid%size(Z_INDEX),  &
                           current_state%local_grid%size(Y_INDEX),  &
                           current_state%local_grid%size(X_INDEX)   )    )
    endif
    if (l_tend_3d_ql) then
      iql=get_q_index(standard_q_names%CLOUD_LIQUID_MASS, 'stepfields')
      allocate( tend_3d_ql(current_state%local_grid%size(Z_INDEX),  &
                           current_state%local_grid%size(Y_INDEX),  &
                           current_state%local_grid%size(X_INDEX)   )    )
    endif
    if (l_tend_3d_qi) then
      iqi=get_q_index(standard_q_names%ICE_MASS, 'stepfields')
      allocate( tend_3d_qi(current_state%local_grid%size(Z_INDEX),  &
                           current_state%local_grid%size(Y_INDEX),  &
                           current_state%local_grid%size(X_INDEX)   )    )
    endif
    if (l_tend_3d_qr) then
      iqr=get_q_index(standard_q_names%RAIN_MASS, 'stepfields')
      allocate( tend_3d_qr(current_state%local_grid%size(Z_INDEX),  &
                           current_state%local_grid%size(Y_INDEX),  &
                           current_state%local_grid%size(X_INDEX)   )    )
    endif
    if (l_tend_3d_qs) then
      iqs=get_q_index(standard_q_names%SNOW_MASS, 'stepfields')
      allocate( tend_3d_qs(current_state%local_grid%size(Z_INDEX),  &
                           current_state%local_grid%size(Y_INDEX),  &
                           current_state%local_grid%size(X_INDEX)   )    )
    endif
    if (l_tend_3d_qg) then
      iqg=get_q_index(standard_q_names%GRAUPEL_MASS, 'stepfields')
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

    ! Save the sampling_frequency to force diagnostic calculation on select time steps
    diagnostic_generation_frequency=options_get_integer(current_state%options_database, "sampling_frequency")

  end subroutine initialisation_callback


  !> Finalisation callback
  !! @param current_state The current model state
  subroutine finalisation_callback(current_state)
    type(model_state_type), target, intent(inout) :: current_state

    deallocate(resetq_min)

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

  !> Called at each timestep and will perform swapping and smoothing as required
  !! @param current_state The current model state_mod
  subroutine timestep_callback(current_state)
    type(model_state_type), target, intent(inout) :: current_state

    integer :: iq
    integer :: current_x_index, current_y_index, target_x_index, target_y_index

    current_x_index=current_state%column_local_x
    current_y_index=current_state%column_local_y
    target_y_index=current_y_index-current_state%local_grid%halo_size(Y_INDEX)
    target_x_index=current_x_index-current_state%local_grid%halo_size(X_INDEX)


    if (cfl_is_enabled .and. current_state%first_timestep_column) then
      if (mod(current_state%timestep, current_state%cfl_frequency) == 1 .or. &
         current_state%timestep-current_state%start_timestep .le. current_state%cfl_frequency) then
        determine_flow_minmax=.true.
        call reset_local_minmax_values(current_state)
      else
        determine_flow_minmax=.false.
      end if
    end if

    if (.not. current_state%halo_column) then
      if (determine_flow_minmax .and. cfl_is_enabled) &
         call determine_local_flow_minmax(current_state, current_state%column_local_y,  current_state%column_local_x)
      call step_all_fields(current_state)
    end if

    ! Remove negative rounding errors
    if (l_nonconservative_positive_q)then
      do iq=1,current_state%number_q_fields
        if (current_state%first_timestep_column)then
          resetq_min(iq)=minval(current_state%zq(iq)%data(:,current_state%column_local_y,  current_state%column_local_x))
        else
          resetq_min(iq)=min(resetq_min(iq),&
             minval(current_state%zq(iq)%data(:,current_state%column_local_y,  current_state%column_local_x)))
        end if
        call remove_negative_rounding_errors_for_single_field(current_state%column_local_x, current_state%column_local_y, &
             current_state%column_local_x-2, current_state%column_local_y-1, current_state%zq(iq), current_state%local_grid)
      end do
    end if


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

    if (mod(current_state%timestep, diagnostic_generation_frequency) == 0 .and. .not. current_state%halo_column) then
      call compute_component_tendencies(current_state, current_x_index, current_y_index, target_x_index, target_y_index)
    end if

  end subroutine timestep_callback

  !> Steps all fields
  !! @param current_state The current model state_mod
  subroutine step_all_fields(current_state)
    type(model_state_type), target, intent(inout) :: current_state

    integer :: x_prev, y_prev, i, k
    real(kind=DEFAULT_PRECISION) :: c1, c2

    x_prev = current_state%column_local_x-2
    y_prev = current_state%column_local_y-1

    c1 = 1.0_DEFAULT_PRECISION - 2.0_DEFAULT_PRECISION*current_state%tsmth                                               
    c2 = current_state%tsmth

#ifdef U_ACTIVE   
    call step_single_field(current_state%column_local_x,  current_state%column_local_y, &
         x_prev, y_prev, current_state%u, current_state%zu, current_state%su, current_state%local_grid, .true., &
         current_state%field_stepping, current_state%dtm, current_state%ugal, c1, c2, .false., current_state%savu)
#endif
#ifdef V_ACTIVE
    call step_single_field(current_state%column_local_x,  current_state%column_local_y, &
         x_prev, y_prev, current_state%v, current_state%zv, current_state%sv, current_state%local_grid, .true., &
         current_state%field_stepping, current_state%dtm, current_state%vgal, c1, c2, .false., current_state%savv)
#endif
#ifdef W_ACTIVE
    call step_single_field(current_state%column_local_x,  current_state%column_local_y, &
         x_prev, y_prev, current_state%w, current_state%zw, current_state%sw, current_state%local_grid, .false., &
         current_state%field_stepping, current_state%dtm, real(0., kind=DEFAULT_PRECISION), c1, c2, .false., current_state%savw)
#endif
    if (current_state%th%active) then       
       call step_single_field(current_state%column_local_x,  current_state%column_local_y, &
         x_prev, y_prev, current_state%th, current_state%zth, current_state%sth, current_state%local_grid, .false., &
         current_state%field_stepping, current_state%dtm, real(0., kind=DEFAULT_PRECISION), c1, c2, &
         current_state%field_stepping == CENTRED_STEPPING)
    endif
    do i=1,current_state%number_q_fields
      if (current_state%q(i)%active) then         
        call step_single_field(current_state%column_local_x,  current_state%column_local_y, x_prev, y_prev, &
             current_state%q(i), current_state%zq(i), current_state%sq(i), current_state%local_grid, .false., &
             current_state%field_stepping, current_state%dtm, real(0., kind=DEFAULT_PRECISION), c1, c2, &
             current_state%field_stepping == CENTRED_STEPPING)
      end if
    end do
  end subroutine step_all_fields 

  !> Determines the minimum and maximum values for the local flow field. These are before the stepping, and are all reduced
  !! later on in the cfl test
  !! @param current_state The current model state
  !! @param local_y The local y index
  !! @param local_x The local x index
  subroutine determine_local_flow_minmax(current_state, local_y, local_x)
    type(model_state_type), target, intent(inout) :: current_state
    integer, intent(in) :: local_x, local_y
    
    integer :: k
    
    do k=2, current_state%local_grid%local_domain_end_index(Z_INDEX)
#ifdef U_ACTIVE
      current_state%local_zumax = max(current_state%local_zumax, current_state%zu%data(k,local_y,local_x))
      current_state%local_zumin = min(current_state%local_zumin, current_state%zu%data(k,local_y,local_x))
#endif
#ifdef V_ACTIVE
      current_state%local_zvmax = max(current_state%local_zvmax, current_state%zv%data(k,local_y,local_x))
      current_state%local_zvmin = min(current_state%local_zvmin, current_state%zv%data(k,local_y,local_x))
#endif
#ifdef W_ACTIVE
      if (k .lt. current_state%local_grid%local_domain_end_index(Z_INDEX)) then
        current_state%abswmax(k) = max(current_state%abswmax(k), abs(current_state%zw%data(k,local_y,local_x)))
      end if
#endif
    end do    
  end subroutine determine_local_flow_minmax

  !> Resets the local min and max values for the flow fields
  !! @param current_state The current model state
  subroutine reset_local_minmax_values(current_state)
    type(model_state_type), intent(inout), target :: current_state

    ! Reset the local values for the next timestep
    current_state%local_zumin=rlargep
    current_state%local_zumax=-rlargep
    current_state%local_zvmin=rlargep
    current_state%local_zvmax=-rlargep
    current_state%abswmax=-rlargep
  end subroutine reset_local_minmax_values  

  !> Removes the negative rounding errors from a specific single field. This works two columns behind and then catches up
  !! on the last column
  !! @param x_local_index The current local x index
  !! @param y_local_index The current local y index
  !! @param x_prev The previous x index to step
  !! @param y_prev The previous y index to step
  !! @param field The prognostic field
  !! @param local_grid Description of the local grid
  subroutine remove_negative_rounding_errors_for_single_field(x_local_index, y_local_index, x_prev, y_prev, field, local_grid)
    integer, intent(in) :: x_local_index, y_local_index, x_prev, y_prev
    type(local_grid_type), intent(inout) :: local_grid
    type(prognostic_field_type), intent(inout) :: field

    if (x_prev .ge. local_grid%local_domain_start_index(X_INDEX)) then
      call remove_negative_rounding_errors_in_slice(y_local_index, x_prev, y_prev, field, local_grid)
    end if
    if (x_local_index == local_grid%local_domain_end_index(X_INDEX)) then
      if (x_local_index .gt. 1) then
        call remove_negative_rounding_errors_in_slice(y_local_index, x_local_index-1, y_prev, field, local_grid)
      end if
      call remove_negative_rounding_errors_in_slice(y_local_index, x_local_index, y_prev, field, local_grid)
    end if    
  end subroutine remove_negative_rounding_errors_for_single_field

  !> Removes the negative rounding errors from a slice of a single field. This works two columns behind and then catches up
  !! on the last column
  !! @param y_local_index The current local y index
  !! @param x_prev The previous x index to step
  !! @param y_prev The previous y index to step
  !! @param field The prognostic field
  !! @param local_grid Description of the local grid
  subroutine remove_negative_rounding_errors_in_slice(y_local_index, x_prev, y_prev, field, local_grid)
    integer, intent(in) :: y_local_index, x_prev, y_prev
    type(local_grid_type), intent(inout) :: local_grid
    type(prognostic_field_type), intent(inout) :: field

    if (y_prev .ge. local_grid%local_domain_start_index(Y_INDEX)) then
      where (field%data(:, y_prev,  x_prev) < 0.0_DEFAULT_PRECISION)
        field%data(:, y_prev,  x_prev)=0.0_DEFAULT_PRECISION
      end where
    end if
    if (y_local_index == local_grid%local_domain_end_index(Y_INDEX)) then
      where (field%data(:, y_local_index,  x_prev) < 0.0_DEFAULT_PRECISION)
        field%data(:, y_local_index,  x_prev)=0.0_DEFAULT_PRECISION
      end where
    end if    
  end subroutine remove_negative_rounding_errors_in_slice  

  !> Steps a single specific field. This will step on the yth column of the x-2 slice and x-1 and x if this is the last slice
  !! @param x_local_index The current local x index
  !! @param y_local_index The current local y index
  !! @param x_prev The previous x index to step
  !! @param y_prev The previous y index to step
  !! @param field The prognostic field
  !! @param zfield Z prognostic field
  !! @param sfield Source terms for the prognostic field
  !! @param local_grid Description of the local grid
  !! @param flow_field Whether or not this is a flow field
  !! @param direction The stepping direction (centred or forward)
  !! @param dtm The delta time per timestep
  !! @param gal Galilean transformation
  !! @param c1 Constant to use in smoothing
  !! @param c2 Constant to use in smoothing
  !! @param do_timesmoothing Whether timesmoothing using Robert filter should be done on the field
  !! @param sav Optional sav field
  subroutine step_single_field(x_local_index, y_local_index, x_prev, y_prev, field, zfield, sfield, local_grid,&
       flow_field, direction, dtm, gal, c1, c2, do_timesmoothing, sav)
    integer, intent(in) :: x_local_index, y_local_index, x_prev, y_prev, direction
    real(kind=DEFAULT_PRECISION), intent(in) :: dtm, gal
    logical, intent(in) :: flow_field, do_timesmoothing
    type(local_grid_type), intent(inout) :: local_grid
    type(prognostic_field_type), intent(inout) :: field, zfield, sfield
    real(kind=DEFAULT_PRECISION), intent(in) :: c1, c2
    type(prognostic_field_type), optional, intent(inout) :: sav

    if (x_prev .ge. local_grid%local_domain_start_index(X_INDEX)) then
      if (present(sav)) then
        call step_column_in_slice(y_local_index, x_prev, y_prev, field, zfield, sfield, local_grid, &
             flow_field, direction, dtm, gal, c1, c2, do_timesmoothing, sav)
      else               
        call step_column_in_slice(y_local_index, x_prev, y_prev, field, zfield, sfield, local_grid, &
             flow_field, direction, dtm, gal, c1, c2, do_timesmoothing)
      end if
    end if

    if (x_local_index == local_grid%local_domain_end_index(X_INDEX)) then
      ! If this is the last slice then process x-1 (if applicable) and x too
      if (x_local_index .gt. 1) then
        if (present(sav)) then
          call step_column_in_slice(y_local_index, x_local_index-1, y_prev, field, zfield, sfield, local_grid, &
               flow_field, direction, dtm, gal, c1, c2, do_timesmoothing, sav)
        else
          call step_column_in_slice(y_local_index, x_local_index-1, y_prev, field, zfield, sfield, local_grid, &
               flow_field, direction, dtm, gal, c1, c2, do_timesmoothing)
        end if
      end if
      if (present(sav)) then
        call step_column_in_slice(y_local_index, x_local_index, y_prev, field, zfield, sfield, local_grid, &
             flow_field, direction, dtm, gal, c1, c2, do_timesmoothing, sav)     
      else
        call step_column_in_slice(y_local_index, x_local_index, y_prev, field, zfield, sfield, local_grid, &
             flow_field, direction, dtm, gal, c1, c2, do_timesmoothing)     
      end if
    end if
  end subroutine step_single_field

  !> Performs initial timesmoothing for a theta or Q field using Robert filter. This is finished off in swapsmooth
  !! @param field The field to smooth
  !! @param zfield The zfield to use in smoothing
  !! @param local_grid Description of the local grid
  !! @param x_index The X index to work on
  !! @param y_index The Y index to work on
  !! @param c1 Constant to use in smoothing
  !! @param c2 Constant to use in smoothing
  subroutine perform_timesmooth_for_field(field, zfield, local_grid, x_index, y_index, c1, c2)
    type(prognostic_field_type), intent(inout) :: field, zfield
    type(local_grid_type), intent(inout) :: local_grid
    integer, intent(in) :: x_index, y_index
    real(kind=DEFAULT_PRECISION), intent(in) :: c1, c2

    integer :: k

    do k=1,local_grid%size(Z_INDEX)
      field%data(k, y_index, x_index)=c1*field%data(k, y_index, x_index)+c2*zfield%data(k, y_index, x_index)
    end do    
  end subroutine perform_timesmooth_for_field 

  !> Will step a column in a specific slice. If y_prev is large enough then will step the y-1 column and if this
  !! is the last column of the slice then will also step the current column
  !! @param x_local_index The current local x index
  !! @param y_local_index The current local y index
  !! @param x_prev The previous x index to step
  !! @param y_prev The previous y index to step
  !! @param field The prognostic field
  !! @param zfield Z prognostic field
  !! @param sfield Source terms for the prognostic field
  !! @param local_grid Description of the local grid
  !! @param flow_field Whether or not this is a flow field
  !! @param direction The stepping direction (centred or forward)
  !! @param dtm The delta time per timestep
  !! @param gal Galilean transformation
  !! @param c1 Constant to use in smoothing
  !! @param c2 Constant to use in smoothing
  !! @param do_timesmoothing Whether timesmoothing using Robert filter should be done on the field
  !! @param sav Optional sav field
  subroutine step_column_in_slice(y_local_index, x_prev, y_prev, field, zfield, sfield, local_grid,&
       flow_field, direction, dtm, gal,  c1, c2, do_timesmoothing, sav)
    integer, intent(in) :: y_local_index, x_prev, y_prev, direction
    real(kind=DEFAULT_PRECISION), intent(in) :: dtm, gal
    logical, intent(in) :: flow_field, do_timesmoothing
    type(local_grid_type), intent(inout) :: local_grid
    type(prognostic_field_type), intent(inout) :: field, zfield, sfield
    real(kind=DEFAULT_PRECISION), intent(in) :: c1, c2
    type(prognostic_field_type), optional, intent(inout) :: sav

    if (y_prev .ge. local_grid%local_domain_start_index(Y_INDEX)) then
      if (do_timesmoothing) then
        call perform_timesmooth_for_field(field, zfield, local_grid, x_prev, y_prev, c1, c2)
      end if      
      if (present(sav)) then
        call step_field(x_prev, y_prev, field, zfield, sfield, local_grid, flow_field, direction, dtm, gal, sav)
      else
        call step_field(x_prev, y_prev, field, zfield, sfield, local_grid, flow_field, direction, dtm, gal)
      end if
    end if

    if (y_local_index == local_grid%local_domain_end_index(Y_INDEX)) then
      if (do_timesmoothing) then
        call perform_timesmooth_for_field(field, zfield, local_grid, x_prev, y_local_index, c1, c2)
      end if 
      if (present(sav)) then
        call step_field(x_prev, y_local_index, field, zfield, sfield, local_grid, flow_field, direction, dtm, gal, sav)
      else
        call step_field(x_prev, y_local_index, field, zfield, sfield, local_grid, flow_field, direction, dtm, gal)
      end if
    end if
  end subroutine step_column_in_slice

  !> Will do the actual field stepping
  !! @param flow_field Whether or not we are stepping a flow field
  !! @param direction 1=forward, 0=centred
  !! @param x_index The local X slice index
  !! @param y_index The local Y column index
  !! @param kkp Points in the vertical column
  !! @param dtm The model timestep
  !! @param field The prognostic field
  !! @param zfield The prognostic z field (which goes to timestep t+1)
  !! @param xfield The tendency of the field 
  !! @param gal The galilean transformation
  subroutine step_field(x_local_index, y_local_index, field, zfield, sfield, local_grid, flow_field, direction, dtm, gal, sav)
    integer, intent(in) :: x_local_index, y_local_index, direction
    real(kind=DEFAULT_PRECISION), intent(in) :: dtm, gal
    logical, intent(in) :: flow_field
    type(local_grid_type), intent(inout) :: local_grid
    type(prognostic_field_type), intent(inout) :: field, zfield, sfield
    type(prognostic_field_type), optional, intent(inout) :: sav

    integer :: k
    real(kind=DEFAULT_PRECISION) :: actual_gal, dtm_x2

    dtm_x2 = 2.0_DEFAULT_PRECISION * dtm

    actual_gal = merge(gal, real(0.0_DEFAULT_PRECISION, kind=DEFAULT_PRECISION), flow_field)

    sfield%data(1,y_local_index, x_local_index)=0.0_DEFAULT_PRECISION

    do k=1,local_grid%size(Z_INDEX)
      ! Save the Z field which is used in the Robert filter
      if (present(sav) .and. direction .eq. CENTRED_STEPPING) &
           sav%data(k,y_local_index, x_local_index) = zfield%data(k, y_local_index, x_local_index) + actual_gal
      if (flow_field) field%data(k, y_local_index, x_local_index) = actual_gal + field%data(k, y_local_index, x_local_index)
      if (direction == FORWARD_STEPPING) then
        zfield%data(k, y_local_index, x_local_index) = field%data(k, y_local_index, x_local_index) + dtm * &
             sfield%data(k, y_local_index, x_local_index)
      else
        zfield%data(k, y_local_index, x_local_index) = actual_gal+zfield%data(k, y_local_index, x_local_index)+dtm_x2*&
             sfield%data(k, y_local_index, x_local_index)
      end if
    end do
  end subroutine step_field

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
      tend_3d_tabs(:,tyn,txn)=   &
          current_state%sth%data(:,cyn,cxn) * current_state%global_grid%configuration%vertical%rprefrcp(:)
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
    strcomp=INDEX(name, "_total_3d_local")
    if (strcomp .ne. 0) then
      field_information%field_type=COMPONENT_ARRAY_FIELD_TYPE
      field_information%number_dimensions=3
      field_information%dimension_sizes(1)=current_state%local_grid%size(Z_INDEX)
      field_information%dimension_sizes(2)=current_state%local_grid%size(Y_INDEX)
      field_information%dimension_sizes(3)=current_state%local_grid%size(X_INDEX)
      field_information%data_type=COMPONENT_DOUBLE_DATA_TYPE

      if      (name .eq. "tend_th_total_3d_local") then
        field_information%enabled=l_tend_3d_th
      else if (name .eq. "tend_qv_total_3d_local") then
        field_information%enabled=l_tend_3d_qv
      else if (name .eq. "tend_ql_total_3d_local") then
        field_information%enabled=l_tend_3d_ql
      else if (name .eq. "tend_qi_total_3d_local") then
        field_information%enabled=l_tend_3d_qi
      else if (name .eq. "tend_qr_total_3d_local") then
        field_information%enabled=l_tend_3d_qr
      else if (name .eq. "tend_qs_total_3d_local") then
        field_information%enabled=l_tend_3d_qs
      else if (name .eq. "tend_qg_total_3d_local") then
        field_information%enabled=l_tend_3d_qg
      else if (name .eq. "tend_tabs_total_3d_local") then
        field_information%enabled=l_tend_3d_tabs
      else
        field_information%enabled=.true.
      end if

    end if !end 3d check

    ! Field information for profiles
   strcomp=INDEX(name, "_total_profile_total_local")
    if (strcomp .ne. 0) then
      field_information%field_type=COMPONENT_ARRAY_FIELD_TYPE
      field_information%number_dimensions=1
      field_information%dimension_sizes(1)=current_state%local_grid%size(Z_INDEX)
      field_information%data_type=COMPONENT_DOUBLE_DATA_TYPE

      if      (name .eq. "tend_th_total_profile_total_local") then
        field_information%enabled=l_tend_pr_tot_th
      else if (name .eq. "tend_qv_total_profile_total_local") then
        field_information%enabled=l_tend_pr_tot_qv
      else if (name .eq. "tend_ql_total_profile_total_local") then
        field_information%enabled=l_tend_pr_tot_ql
      else if (name .eq. "tend_qi_total_profile_total_local") then
        field_information%enabled=l_tend_pr_tot_qi
      else if (name .eq. "tend_qr_total_profile_total_local") then
        field_information%enabled=l_tend_pr_tot_qr
      else if (name .eq. "tend_qs_total_profile_total_local") then
        field_information%enabled=l_tend_pr_tot_qs
      else if (name .eq. "tend_qg_total_profile_total_local") then
        field_information%enabled=l_tend_pr_tot_qg
      else if (name .eq. "tend_tabs_total_profile_total_local") then
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
    if      (name .eq. "tend_th_total_3d_local" .and. allocated(tend_3d_th)) then
      call set_published_field_value(field_value, real_3d_field=tend_3d_th)
    else if (name .eq. "tend_qv_total_3d_local" .and. allocated(tend_3d_qv)) then
      call set_published_field_value(field_value, real_3d_field=tend_3d_qv)
    else if (name .eq. "tend_ql_total_3d_local" .and. allocated(tend_3d_ql)) then
      call set_published_field_value(field_value, real_3d_field=tend_3d_ql)
    else if (name .eq. "tend_qi_total_3d_local" .and. allocated(tend_3d_qi)) then
      call set_published_field_value(field_value, real_3d_field=tend_3d_qi)
    else if (name .eq. "tend_qr_total_3d_local" .and. allocated(tend_3d_qr)) then
      call set_published_field_value(field_value, real_3d_field=tend_3d_qr)
    else if (name .eq. "tend_qs_total_3d_local" .and. allocated(tend_3d_qs)) then
      call set_published_field_value(field_value, real_3d_field=tend_3d_qs)
    else if (name .eq. "tend_qg_total_3d_local" .and. allocated(tend_3d_qg)) then
      call set_published_field_value(field_value, real_3d_field=tend_3d_qg)
    else if (name .eq. "tend_tabs_total_3d_local" .and. allocated(tend_3d_tabs)) then
      call set_published_field_value(field_value, real_3d_field=tend_3d_tabs)

    ! Profile Tendency Fields
    else if (name .eq. "tend_th_total_profile_total_local" .and. allocated(tend_pr_tot_th)) then
      call set_published_field_value(field_value, real_1d_field=tend_pr_tot_th)
    else if (name .eq. "tend_qv_total_profile_total_local" .and. allocated(tend_pr_tot_qv)) then
      call set_published_field_value(field_value, real_1d_field=tend_pr_tot_qv)
    else if (name .eq. "tend_ql_total_profile_total_local" .and. allocated(tend_pr_tot_ql)) then
      call set_published_field_value(field_value, real_1d_field=tend_pr_tot_ql)
    else if (name .eq. "tend_qi_total_profile_total_local" .and. allocated(tend_pr_tot_qi)) then
      call set_published_field_value(field_value, real_1d_field=tend_pr_tot_qi)
    else if (name .eq. "tend_qr_total_profile_total_local" .and. allocated(tend_pr_tot_qr)) then
      call set_published_field_value(field_value, real_1d_field=tend_pr_tot_qr)
    else if (name .eq. "tend_qs_total_profile_total_local" .and. allocated(tend_pr_tot_qs)) then
      call set_published_field_value(field_value, real_1d_field=tend_pr_tot_qs)
    else if (name .eq. "tend_qg_total_profile_total_local" .and. allocated(tend_pr_tot_qg)) then
      call set_published_field_value(field_value, real_1d_field=tend_pr_tot_qg)
    else if (name .eq. "tend_tabs_total_profile_total_local" .and. allocated(tend_pr_tot_tabs)) then
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

end module stepfields_mod
