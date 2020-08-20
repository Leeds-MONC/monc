!> Stepping of the pressure field. Completes the time-stepping of the velocity fields
!! by adding the pressure term (dp/dx_i). In addition, ensures that l_zu and l_zv satisfy the
!! Galilean-transformed boundary condition. This does not do the flow field _p terms which are only
!! needed for diagnostics, nore does it do field halo swapping which is again only needed for diagnostics
module pstep_mod
  use monc_component_mod, only : COMPONENT_ARRAY_FIELD_TYPE, &
      COMPONENT_DOUBLE_DATA_TYPE, component_descriptor_type, &
      component_field_value_type, component_field_information_type
  use optionsdatabase_mod, only : options_get_integer
  use state_mod, only : model_state_type, CENTRED_STEPPING
  use grids_mod, only : Z_INDEX, X_INDEX, Y_INDEX
  use datadefn_mod, only : DEFAULT_PRECISION
  use registry_mod, only : is_component_enabled
  use logging_mod, only : LOG_ERROR, log_master_log

  implicit none

#ifndef TEST_MODE
  private
#endif

  ! Local tendency diagnostic variables for this component
  !
  ! This one is a little different.  tendp_ terms are after the tendency due to the pressure term,
  ! and the tend_ terms will collect the total tendency for the time step.
  !
  ! 3D tendency fields and logicals for their use
  real(kind=DEFAULT_PRECISION), dimension(:,:,:), allocatable ::     &
               tendp_3d_u,   tendp_3d_v,   tendp_3d_w,   tend_3d_u,   tend_3d_v,   tend_3d_w
  logical :: l_tendp_3d_u, l_tendp_3d_v, l_tendp_3d_w, l_tend_3d_u, l_tend_3d_v, l_tend_3d_w
  ! Local mean tendency profile fields and logicals for their use
  real(kind=DEFAULT_PRECISION), dimension(:), allocatable ::         &
               tendp_pr_tot_u,   tendp_pr_tot_v,   tendp_pr_tot_w,   tend_pr_tot_u,   tend_pr_tot_v,   tend_pr_tot_w
  logical :: l_tendp_pr_tot_u, l_tendp_pr_tot_v, l_tendp_pr_tot_w, l_tend_pr_tot_u, l_tend_pr_tot_v, l_tend_pr_tot_w

  public pstep_get_descriptor

contains

  !> Descriptor of this component for registration
  !! @returns The pstep component descriptor
  type(component_descriptor_type) function pstep_get_descriptor()
    pstep_get_descriptor%name="pstep"
    pstep_get_descriptor%version=0.1
    pstep_get_descriptor%initialisation=>initialisation_callback
    pstep_get_descriptor%timestep=>timestep_callback
    pstep_get_descriptor%finalisation=>finalisation_callback

    pstep_get_descriptor%field_value_retrieval=>field_value_retrieval_callback
    pstep_get_descriptor%field_information_retrieval=>field_information_retrieval_callback
    allocate(pstep_get_descriptor%published_fields(6+6))

    pstep_get_descriptor%published_fields(1)= "tend_u_pstep_3d_local"
    pstep_get_descriptor%published_fields(2)= "tend_v_pstep_3d_local"
    pstep_get_descriptor%published_fields(3)= "tend_w_pstep_3d_local"
    pstep_get_descriptor%published_fields(4)= "tend_u_total_3d_local"
    pstep_get_descriptor%published_fields(5)= "tend_v_total_3d_local"
    pstep_get_descriptor%published_fields(6)= "tend_w_total_3d_local"

    pstep_get_descriptor%published_fields(6+1)= "tend_u_pstep_profile_total_local"
    pstep_get_descriptor%published_fields(6+2)= "tend_v_pstep_profile_total_local"
    pstep_get_descriptor%published_fields(6+3)= "tend_w_pstep_profile_total_local"
    pstep_get_descriptor%published_fields(6+4)= "tend_u_total_profile_total_local"
    pstep_get_descriptor%published_fields(6+5)= "tend_v_total_profile_total_local"
    pstep_get_descriptor%published_fields(6+6)= "tend_w_total_profile_total_local"

  end function pstep_get_descriptor

  !> Field information retrieval callback, this returns information for a specific components published field
  !! @param current_state Current model state
  !! @param name The name of the field to retrieve information for
  !! @param field_information Populated with information about the field
  subroutine field_information_retrieval_callback(current_state, name, field_information)
    type(model_state_type), target, intent(inout) :: current_state
    character(len=*), intent(in) :: name
    type(component_field_information_type), intent(out) :: field_information
    integer :: strcomp

    ! Field information for 3d
    strcomp=INDEX(name, "_pstep_3d_local")
    if (strcomp .ne. 0) then
      field_information%field_type=COMPONENT_ARRAY_FIELD_TYPE
      field_information%number_dimensions=3
      field_information%dimension_sizes(1)=current_state%local_grid%size(Z_INDEX)
      field_information%dimension_sizes(2)=current_state%local_grid%size(Y_INDEX)
      field_information%dimension_sizes(3)=current_state%local_grid%size(X_INDEX)
      field_information%data_type=COMPONENT_DOUBLE_DATA_TYPE

      if      (name .eq. "tend_u_pstep_3d_local") then
        field_information%enabled=l_tend_3d_u
      else if (name .eq. "tend_v_pstep_3d_local") then
        field_information%enabled=l_tend_3d_v
      else if (name .eq. "tend_w_pstep_3d_local") then
        field_information%enabled=l_tend_3d_w
      else
        field_information%enabled=.true.
      end if

    end if !end 3d check

    ! Field information for profiles
    strcomp=INDEX(name, "_pstep_profile_total_local")
    if (strcomp .ne. 0) then
      field_information%field_type=COMPONENT_ARRAY_FIELD_TYPE
      field_information%number_dimensions=1
      field_information%dimension_sizes(1)=current_state%local_grid%size(Z_INDEX)
      field_information%data_type=COMPONENT_DOUBLE_DATA_TYPE

      if      (name .eq. "tend_u_pstep_profile_total_local") then
        field_information%enabled=l_tend_pr_tot_u
      else if (name .eq. "tend_v_pstep_profile_total_local") then
        field_information%enabled=l_tend_pr_tot_v
      else if (name .eq. "tend_w_pstep_profile_total_local") then
        field_information%enabled=l_tend_pr_tot_w
      else
        field_information%enabled=.true.
      end if

    end if !end profile check

    ! Field information for 3d totals
    strcomp=INDEX(name, "_total_3d_local")
    if (strcomp .ne. 0) then
      field_information%field_type=COMPONENT_ARRAY_FIELD_TYPE
      field_information%number_dimensions=3
      field_information%dimension_sizes(1)=current_state%local_grid%size(Z_INDEX)
      field_information%dimension_sizes(2)=current_state%local_grid%size(Y_INDEX)
      field_information%dimension_sizes(3)=current_state%local_grid%size(X_INDEX)
      field_information%data_type=COMPONENT_DOUBLE_DATA_TYPE

      if      (name .eq. "tend_u_total_3d_local") then
        field_information%enabled=l_tend_3d_u
      else if (name .eq. "tend_v_total_3d_local") then
        field_information%enabled=l_tend_3d_v
      else if (name .eq. "tend_w_total_3d_local") then
        field_information%enabled=l_tend_3d_w
      else
        field_information%enabled=.true.
      end if

    end if !end 3d check

    ! Field information for total profiles
    strcomp=INDEX(name, "_total_profile_total_local")
    if (strcomp .ne. 0) then
      field_information%field_type=COMPONENT_ARRAY_FIELD_TYPE
      field_information%number_dimensions=1
      field_information%dimension_sizes(1)=current_state%local_grid%size(Z_INDEX)
      field_information%data_type=COMPONENT_DOUBLE_DATA_TYPE

      if      (name .eq. "tend_u_total_profile_total_local") then
        field_information%enabled=l_tend_pr_tot_u
      else if (name .eq. "tend_v_total_profile_total_local") then
        field_information%enabled=l_tend_pr_tot_v
      else if (name .eq. "tend_w_total_profile_total_local") then
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
    ! 3d Tendency Fields

    if      (name .eq. "tend_u_pstep_3d_local" .and. allocated(tendp_3d_u)) then
      call set_published_field_value(field_value, real_3d_field=tendp_3d_u)
    else if (name .eq. "tend_v_pstep_3d_local" .and. allocated(tendp_3d_v)) then
      call set_published_field_value(field_value, real_3d_field=tendp_3d_v)
    else if (name .eq. "tend_w_pstep_3d_local" .and. allocated(tendp_3d_w)) then
      call set_published_field_value(field_value, real_3d_field=tendp_3d_w)

    ! Profile Tendency Fields
    else if (name .eq. "tend_u_pstep_profile_total_local" .and. allocated(tendp_pr_tot_u)) then
      call set_published_field_value(field_value, real_1d_field=tendp_pr_tot_u)
    else if (name .eq. "tend_v_pstep_profile_total_local" .and. allocated(tendp_pr_tot_v)) then
      call set_published_field_value(field_value, real_1d_field=tendp_pr_tot_v)
    else if (name .eq. "tend_w_pstep_profile_total_local" .and. allocated(tendp_pr_tot_w)) then
      call set_published_field_value(field_value, real_1d_field=tendp_pr_tot_w)

    ! 3d Tendency Fields
    else if (name .eq. "tend_u_total_3d_local" .and. allocated(tend_3d_u)) then
      call set_published_field_value(field_value, real_3d_field=tend_3d_u)
    else if (name .eq. "tend_v_total_3d_local" .and. allocated(tend_3d_v)) then
      call set_published_field_value(field_value, real_3d_field=tend_3d_v)
    else if (name .eq. "tend_w_total_3d_local" .and. allocated(tend_3d_w)) then
      call set_published_field_value(field_value, real_3d_field=tend_3d_w)

    ! Profile Tendency Fields
    else if (name .eq. "tend_u_total_profile_total_local" .and. allocated(tend_pr_tot_u)) then
      call set_published_field_value(field_value, real_1d_field=tend_pr_tot_u)
    else if (name .eq. "tend_v_total_profile_total_local" .and. allocated(tend_pr_tot_v)) then
      call set_published_field_value(field_value, real_1d_field=tend_pr_tot_v)
    else if (name .eq. "tend_w_total_profile_total_local" .and. allocated(tend_pr_tot_w)) then
      call set_published_field_value(field_value, real_1d_field=tend_pr_tot_w)
    end if

  end subroutine field_value_retrieval_callback

  !> Initialisation callback hook which will check the diverr component is enabled (as this allocates p)
  !! @param current_state The current model state
  subroutine initialisation_callback(current_state)
    type(model_state_type), target, intent(inout) :: current_state

    if (.not. is_component_enabled(current_state%options_database, "diverr")) then
      call log_master_log(LOG_ERROR, "The pstep component requires the diverr component to be enabled")
    end if

    ! Tendency Logicals
    l_tendp_pr_tot_u   = current_state%u%active
    l_tendp_pr_tot_v   = current_state%v%active
    l_tendp_pr_tot_w   = current_state%w%active
    l_tend_pr_tot_u   = current_state%u%active
    l_tend_pr_tot_v   = current_state%v%active
    l_tend_pr_tot_w   = current_state%w%active

    l_tendp_3d_u   = current_state%u%active .or. l_tendp_pr_tot_u
    l_tendp_3d_v   = current_state%v%active .or. l_tendp_pr_tot_v
    l_tendp_3d_w   = current_state%w%active .or. l_tendp_pr_tot_w
    l_tend_3d_u   = current_state%u%active .or. l_tend_pr_tot_u
    l_tend_3d_v   = current_state%v%active .or. l_tend_pr_tot_v
    l_tend_3d_w   = current_state%w%active .or. l_tend_pr_tot_w

    ! Allocate 3d tendency fields upon availability
    if (l_tendp_3d_u) then
      allocate( tendp_3d_u(current_state%local_grid%size(Z_INDEX),  &
                           current_state%local_grid%size(Y_INDEX),  &
                           current_state%local_grid%size(X_INDEX)   )    )
    endif
    if (l_tendp_3d_v) then
      allocate( tendp_3d_v(current_state%local_grid%size(Z_INDEX),  &
                           current_state%local_grid%size(Y_INDEX),  &
                           current_state%local_grid%size(X_INDEX)   )    )
    endif
    if (l_tendp_3d_w) then
      allocate( tendp_3d_w(current_state%local_grid%size(Z_INDEX),  &
                           current_state%local_grid%size(Y_INDEX),  &
                           current_state%local_grid%size(X_INDEX)   )    )
    endif
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
    if (l_tendp_pr_tot_u) then
      allocate( tendp_pr_tot_u(current_state%local_grid%size(Z_INDEX)) )
    endif
    if (l_tendp_pr_tot_v) then
      allocate( tendp_pr_tot_v(current_state%local_grid%size(Z_INDEX)) )
    endif
    if (l_tendp_pr_tot_w) then
      allocate( tendp_pr_tot_w(current_state%local_grid%size(Z_INDEX)) )
    endif
    if (l_tend_pr_tot_u) then
      allocate( tend_pr_tot_u(current_state%local_grid%size(Z_INDEX)) )
    endif
    if (l_tend_pr_tot_v) then
      allocate( tend_pr_tot_v(current_state%local_grid%size(Z_INDEX)) )
    endif
    if (l_tend_pr_tot_w) then
      allocate( tend_pr_tot_w(current_state%local_grid%size(Z_INDEX)) )
    endif

  end subroutine initialisation_callback  


  subroutine finalisation_callback(current_state)
    type(model_state_type), target, intent(inout) :: current_state

    if (allocated(tendp_3d_u)) deallocate(tendp_3d_u)
    if (allocated(tendp_3d_v)) deallocate(tendp_3d_v)
    if (allocated(tendp_3d_w)) deallocate(tendp_3d_w)
    if (allocated(tend_3d_u)) deallocate(tend_3d_u)
    if (allocated(tend_3d_v)) deallocate(tend_3d_v)
    if (allocated(tend_3d_w)) deallocate(tend_3d_w)

    if (allocated(tendp_pr_tot_u)) deallocate(tendp_pr_tot_u)
    if (allocated(tendp_pr_tot_v)) deallocate(tendp_pr_tot_v)
    if (allocated(tendp_pr_tot_w)) deallocate(tendp_pr_tot_w)
    if (allocated(tend_pr_tot_u)) deallocate(tend_pr_tot_u)
    if (allocated(tend_pr_tot_v)) deallocate(tend_pr_tot_v)
    if (allocated(tend_pr_tot_w)) deallocate(tend_pr_tot_w)

  end subroutine finalisation_callback

  !> Called each timestep, this will step the pressure field for the non halo columns
  !! @param current_state The current model state
  subroutine timestep_callback(current_state)
    type(model_state_type), target, intent(inout) :: current_state
    integer :: local_y, local_x, target_x_index, target_y_index
    logical :: calculate_diagnostics

    calculate_diagnostics = current_state%diagnostic_sample_timestep &
                            .and. .not. current_state%halo_column

    local_y=current_state%column_local_y
    local_x=current_state%column_local_x
    target_y_index=local_y-current_state%local_grid%halo_size(Y_INDEX)
    target_x_index=local_x-current_state%local_grid%halo_size(X_INDEX)

    ! Zero profile tendency totals on first instance in the sum
    if (current_state%first_timestep_column) then
      if (l_tendp_pr_tot_u) then
        tendp_pr_tot_u(:)= 0.0_DEFAULT_PRECISION
      endif
      if (l_tendp_pr_tot_v) then
        tendp_pr_tot_v(:)= 0.0_DEFAULT_PRECISION
      endif
      if (l_tendp_pr_tot_w) then
        tendp_pr_tot_w(:)= 0.0_DEFAULT_PRECISION
      endif
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
    
    if (current_state%galilean_transformation) call perform_galilean_transformation(current_state, &
         current_state%column_local_y, current_state%column_local_x)

    if (calculate_diagnostics) call save_precomponent_tendencies(current_state, local_x, local_y, target_x_index, target_y_index)

    if (.not. current_state%halo_column) call step_pressure_field(current_state)

    if (calculate_diagnostics) call compute_component_tendencies(current_state, local_x, local_y, target_x_index, target_y_index)

  end subroutine timestep_callback


  !> Does the actual stepping of the pressure field
  !! @param current_state The current model state
  subroutine step_pressure_field(current_state)
    type(model_state_type), target, intent(inout) :: current_state

    integer :: k, x_index, y_index
    real(kind=DEFAULT_PRECISION) :: dtmtmp

    x_index=current_state%column_local_x
    y_index=current_state%column_local_y

    dtmtmp=merge(current_state%dtm, 0.5_DEFAULT_PRECISION*current_state%dtm, current_state%field_stepping == CENTRED_STEPPING)
    do k=2,current_state%local_grid%size(Z_INDEX)

#ifdef U_ACTIVE
      current_state%zu%data(k, y_index, x_index)= current_state%zu%data(k, y_index, x_index)+ 2.0_DEFAULT_PRECISION*&
           current_state%global_grid%configuration%horizontal%cx*dtmtmp*(current_state%p%data(k, y_index, x_index)-&
           current_state%p%data(k, y_index, x_index+1))      
#endif
#ifdef V_ACTIVE
      current_state%zv%data(k, y_index, x_index)=&
           current_state%zv%data(k, y_index, x_index)+2.0_DEFAULT_PRECISION*&
           current_state%global_grid%configuration%horizontal%cy*dtmtmp*&
           (current_state%p%data(k, y_index, x_index) - current_state%p%data(k, y_index+1, x_index))
#endif
#ifdef W_ACTIVE
      if (k .lt. current_state%local_grid%size(Z_INDEX)) then
        current_state%zw%data(k, y_index, x_index)=current_state%zw%data(k, y_index, x_index)+2.0_DEFAULT_PRECISION*&
             current_state%global_grid%configuration%vertical%rdzn(k+1)*dtmtmp*(current_state%p%data(k, y_index, x_index)-&
             current_state%p%data(k+1, y_index, x_index))
      end if
#endif
    end do
    if (current_state%use_viscosity_and_diffusion .and. current_state%use_surface_boundary_conditions) then
#ifdef U_ACTIVE
      current_state%zu%data(1, y_index, x_index)=-current_state%zu%data(2, y_index, x_index)-&
           2.0_DEFAULT_PRECISION*current_state%ugal
#endif
#ifdef V_ACTIVE
      current_state%zv%data(1, y_index, x_index)=-current_state%zv%data(2, y_index, x_index)-&
           2.0_DEFAULT_PRECISION*current_state%vgal
#endif
    else
#ifdef U_ACTIVE
      current_state%zu%data(1, y_index, x_index)=current_state%zu%data(2, y_index, x_index)
#endif
#ifdef V_ACTIVE
      current_state%zv%data(1, y_index, x_index)=current_state%zv%data(2, y_index, x_index)
#endif
    end if
  end subroutine step_pressure_field  

  !> Performs Galilean transformation of flow current and z fields.
  !! @param current_state The current model state
  !! @param y_index Local y index which we are working with
  !! @param x_index Local x index which we are working with
  subroutine perform_galilean_transformation(current_state, y_index, x_index)
    type(model_state_type), target, intent(inout) :: current_state
    integer, intent(in) :: y_index, x_index

    integer :: k

    do k=1,current_state%local_grid%size(Z_INDEX)
#ifdef U_ACTIVE
      current_state%u%data(k, y_index, x_index)= current_state%u%data(k, y_index, x_index)-current_state%ugal
      current_state%zu%data(k, y_index, x_index)= current_state%zu%data(k, y_index, x_index)-current_state%ugal
#endif
#ifdef V_ACTIVE
      current_state%v%data(k, y_index, x_index)= current_state%v%data(k, y_index, x_index)-current_state%vgal
      current_state%zv%data(k, y_index, x_index)= current_state%zv%data(k, y_index, x_index)-current_state%vgal
#endif
    end do
  end subroutine perform_galilean_transformation


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
    if (l_tendp_3d_u) then
      tendp_3d_u(:,tyn,txn)=current_state%zu%data(:,cyn,cxn)
    endif
    if (l_tendp_3d_v) then
      tendp_3d_v(:,tyn,txn)=current_state%zv%data(:,cyn,cxn)
    endif
    if (l_tendp_3d_w) then
      tendp_3d_w(:,tyn,txn)=current_state%zw%data(:,cyn,cxn)
    endif
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
    real(kind=DEFAULT_PRECISION) :: dtmtmp

    dtmtmp=2.0_DEFAULT_PRECISION* &
         merge(current_state%dtm, 0.5_DEFAULT_PRECISION*current_state%dtm, current_state%field_stepping == CENTRED_STEPPING)

    ! Calculate change in tendency due to component
    if (l_tendp_3d_u) then
      tendp_3d_u(:,tyn,txn)=(current_state%zu%data(:,cyn,cxn) - tendp_3d_u(:,tyn,txn))/dtmtmp
    endif
    if (l_tendp_3d_v) then
      tendp_3d_v(:,tyn,txn)=(current_state%zv%data(:,cyn,cxn) - tendp_3d_v(:,tyn,txn))/dtmtmp
    endif
    if (l_tendp_3d_w) then
      tendp_3d_w(:,tyn,txn)=(current_state%zw%data(:,cyn,cxn) - tendp_3d_w(:,tyn,txn))/dtmtmp
    endif
    if (l_tend_3d_u) then
      tend_3d_u(:,tyn,txn)=tend_3d_u(:,tyn,txn) + tendp_3d_u(:,tyn,txn)
    endif
    if (l_tend_3d_v) then
      tend_3d_v(:,tyn,txn)=tend_3d_v(:,tyn,txn) + tendp_3d_v(:,tyn,txn)
    endif
    if (l_tend_3d_w) then
      tend_3d_w(:,tyn,txn)=tend_3d_w(:,tyn,txn) + tendp_3d_w(:,tyn,txn)
    endif

   ! Add local tendency fields to the profile total
    if (l_tendp_pr_tot_u) then
      tendp_pr_tot_u(:)=tendp_pr_tot_u(:) + tendp_3d_u(:,tyn,txn)
    endif
    if (l_tendp_pr_tot_v) then
      tendp_pr_tot_v(:)=tendp_pr_tot_v(:) + tendp_3d_v(:,tyn,txn)
    endif
    if (l_tendp_pr_tot_w) then
      tendp_pr_tot_w(:)=tendp_pr_tot_w(:) + tendp_3d_w(:,tyn,txn)
    endif
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

end module pstep_mod
