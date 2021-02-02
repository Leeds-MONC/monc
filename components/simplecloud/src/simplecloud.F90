!> A very simple saturation adjustment scheme without any microphysics
module simplecloud_mod
  use datadefn_mod, only : DEFAULT_PRECISION
  use monc_component_mod, only : COMPONENT_ARRAY_FIELD_TYPE, &
      COMPONENT_DOUBLE_DATA_TYPE, component_descriptor_type, &
      component_field_value_type, component_field_information_type
  use optionsdatabase_mod, only : options_get_real, options_get_integer
  use state_mod, only : FORWARD_STEPPING, model_state_type
  use grids_mod, only : Z_INDEX, X_INDEX, Y_INDEX
  use science_constants_mod, only : r_over_cp, rlvap_over_cp
  use saturation_mod, only: qsaturation, dqwsatdt
  use q_indices_mod, only: get_q_index, standard_q_names
  use registry_mod, only : is_component_enabled
  use logging_mod, only : LOG_ERROR, log_master_log

implicit none

#ifndef TEST_MODE
  private
#endif
    
  ! Indices for vapour and cloud
  integer :: iqv, iql

  integer :: k_cloudmax ! max k index for height
  real(kind=DEFAULT_PRECISION) :: max_height_cloud

  ! Local tendency diagnostic variables for this component
  ! 3D tendency fields and logicals for their use
  real(kind=DEFAULT_PRECISION), dimension(:,:,:), allocatable ::     &
                   tend_3d_th,      tend_3d_qv,      tend_3d_ql,      tend_3d_tabs
  logical ::     l_tend_3d_th,    l_tend_3d_qv,    l_tend_3d_ql,    l_tend_3d_tabs
  ! Local mean tendency profile fields and logicals for their use
  real(kind=DEFAULT_PRECISION), dimension(:), allocatable ::         &
               tend_pr_tot_th,  tend_pr_tot_qv,  tend_pr_tot_ql,  tend_pr_tot_tabs
  logical :: l_tend_pr_tot_th,l_tend_pr_tot_qv,l_tend_pr_tot_ql,l_tend_pr_tot_tabs

  public simplecloud_get_descriptor

contains

  !> Provides the descriptor back to the caller and is used in component registration
  !! @returns The termination check component descriptor
  type(component_descriptor_type) function simplecloud_get_descriptor()
    simplecloud_get_descriptor%name="simplecloud"
    simplecloud_get_descriptor%version=0.1
    simplecloud_get_descriptor%initialisation=>initialisation_callback
    simplecloud_get_descriptor%timestep=>timestep_callback
    simplecloud_get_descriptor%finalisation=>finalisation_callback

    simplecloud_get_descriptor%field_value_retrieval=>field_value_retrieval_callback
    simplecloud_get_descriptor%field_information_retrieval=>field_information_retrieval_callback
    allocate(simplecloud_get_descriptor%published_fields(4+4))

    simplecloud_get_descriptor%published_fields(1)="tend_th_simplecloud_3d_local"
    simplecloud_get_descriptor%published_fields(2)="tend_qv_simplecloud_3d_local"
    simplecloud_get_descriptor%published_fields(3)="tend_ql_simplecloud_3d_local"
    simplecloud_get_descriptor%published_fields(4)="tend_tabs_simplecloud_3d_local"

    simplecloud_get_descriptor%published_fields(4+1)="tend_th_simplecloud_profile_total_local"
    simplecloud_get_descriptor%published_fields(4+2)="tend_qv_simplecloud_profile_total_local"
    simplecloud_get_descriptor%published_fields(4+3)="tend_ql_simplecloud_profile_total_local"
    simplecloud_get_descriptor%published_fields(4+4)="tend_tabs_simplecloud_profile_total_local"

  end function simplecloud_get_descriptor


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
    strcomp=INDEX(name, "simplecloud_3d_local")
    if (strcomp .ne. 0) then
      field_information%field_type=COMPONENT_ARRAY_FIELD_TYPE
      field_information%number_dimensions=3
      field_information%dimension_sizes(1)=current_state%local_grid%size(Z_INDEX)
      field_information%dimension_sizes(2)=current_state%local_grid%size(Y_INDEX)
      field_information%dimension_sizes(3)=current_state%local_grid%size(X_INDEX)
      field_information%data_type=COMPONENT_DOUBLE_DATA_TYPE

      if      (name .eq. "tend_th_simplecloud_3d_local") then
        field_information%enabled=l_tend_3d_th
      else if (name .eq. "tend_qv_simplecloud_3d_local") then
        field_information%enabled=l_tend_3d_qv
      else if (name .eq. "tend_ql_simplecloud_3d_local") then
        field_information%enabled=l_tend_3d_ql
      else if (name .eq. "tend_tabs_simplecloud_3d_local") then
        field_information%enabled=l_tend_3d_tabs
      else
        field_information%enabled=.true.
      end if

    end if !end 3d check

    ! Field information for profiles
    strcomp=INDEX(name, "simplecloud_profile_total_local")
    if (strcomp .ne. 0) then
      field_information%field_type=COMPONENT_ARRAY_FIELD_TYPE
      field_information%number_dimensions=1
      field_information%dimension_sizes(1)=current_state%local_grid%size(Z_INDEX)
      field_information%data_type=COMPONENT_DOUBLE_DATA_TYPE

      if      (name .eq. "tend_th_simplecloud_profile_total_local") then
        field_information%enabled=l_tend_pr_tot_th
      else if (name .eq. "tend_qv_simplecloud_profile_total_local") then
        field_information%enabled=l_tend_pr_tot_qv
      else if (name .eq. "tend_ql_simplecloud_profile_total_local") then
        field_information%enabled=l_tend_pr_tot_ql
      else if (name .eq. "tend_tabs_simplecloud_profile_total_local") then
        field_information%enabled=l_tend_pr_tot_tabs
      else
        field_information%enabled=.true.
      end if

    end if !end profile check

  end subroutine field_information_retrieval_callback


  !> The initialisation callback sets up the moisture fields
  !! @param current_state The current model state
  !> Field value retrieval callback, this returns the value of a specific published field
  !! @param current_state Current model state
  !! @param name The name of the field to retrieve the value for
  !! @param field_value Populated with the value of the field
  subroutine field_value_retrieval_callback(current_state, name, field_value)
    type(model_state_type), target, intent(inout) :: current_state
    character(len=*), intent(in) :: name
    type(component_field_value_type), intent(out) :: field_value

    ! 3d Tendency Fields
    if      (name .eq. "tend_th_simplecloud_3d_local" .and. allocated(tend_3d_th)) then
      call set_published_field_value(field_value, real_3d_field=tend_3d_th)
    else if (name .eq. "tend_qv_simplecloud_3d_local" .and. allocated(tend_3d_qv)) then
      call set_published_field_value(field_value, real_3d_field=tend_3d_qv)
    else if (name .eq. "tend_ql_simplecloud_3d_local" .and. allocated(tend_3d_ql)) then
      call set_published_field_value(field_value, real_3d_field=tend_3d_ql)
    else if (name .eq. "tend_tabs_simplecloud_3d_local" .and. allocated(tend_3d_tabs)) then
      call set_published_field_value(field_value, real_3d_field=tend_3d_tabs)

    ! Profile Tendency Fields
    else if (name .eq. "tend_th_simplecloud_profile_total_local" .and. allocated(tend_pr_tot_th)) then
      call set_published_field_value(field_value, real_1d_field=tend_pr_tot_th)
    else if (name .eq. "tend_qv_simplecloud_profile_total_local" .and. allocated(tend_pr_tot_qv)) then
      call set_published_field_value(field_value, real_1d_field=tend_pr_tot_qv)
    else if (name .eq. "tend_ql_simplecloud_profile_total_local" .and. allocated(tend_pr_tot_ql)) then
      call set_published_field_value(field_value, real_1d_field=tend_pr_tot_ql)
    else if (name .eq. "tend_tabs_simplecloud_profile_total_local" .and. allocated(tend_pr_tot_tabs)) then
      call set_published_field_value(field_value, real_1d_field=tend_pr_tot_tabs)
    end if

  end subroutine field_value_retrieval_callback


  subroutine initialisation_callback(current_state)
    type(model_state_type), target, intent(inout) :: current_state

    integer :: k ! look counter
    logical :: l_qdiag

    if (is_component_enabled(current_state%options_database, "casim")) then
      call log_master_log(LOG_ERROR, "Casim and Simplecloud are enabled, this does not work yet. Please disable one")
    end if 

    iqv=get_q_index(standard_q_names%VAPOUR, 'simplecloud')
    iql=get_q_index(standard_q_names%CLOUD_LIQUID_MASS, 'simplecloud')

    ! set buoyancy coefficient (value for vapour should be set
    ! elsewhere for a moist model
    if (.not. allocated(current_state%cq))then
      allocate(current_state%cq(current_state%number_q_fields))
      current_state%cq=0.0_DEFAULT_PRECISION
    end if
    current_state%cq(iql) = -1.0 

    max_height_cloud=options_get_real(current_state%options_database, "max_height_cloud")
    do k=2, current_state%local_grid%size(Z_INDEX)-1
      if (current_state%global_grid%configuration%vertical%zn(k) > max_height_cloud) exit
    end do
    k_cloudmax=k-1

    ! Set tendency diagnostic logicals based on availability
    ! Need to use 3d tendencies to compute the profiles, so they will be allocated
    !      in the case where profiles are available
    l_qdiag =  (.not. current_state%passive_q .and. current_state%number_q_fields .gt. 0)

    l_tend_pr_tot_th  = current_state%th%active
    l_tend_pr_tot_qv  = l_qdiag .and. current_state%number_q_fields .ge. 1
    l_tend_pr_tot_ql  = l_qdiag .and. current_state%number_q_fields .ge. 2
    l_tend_pr_tot_tabs  = l_tend_pr_tot_th

    l_tend_3d_th  = current_state%th%active .or. l_tend_pr_tot_th
    l_tend_3d_qv  = (l_qdiag .and. current_state%number_q_fields .ge. 1) .or. l_tend_pr_tot_qv
    l_tend_3d_ql  = (l_qdiag .and. current_state%number_q_fields .ge. 2) .or. l_tend_pr_tot_ql
    l_tend_3d_tabs = l_tend_3d_th

    ! Allocate 3d tendency fields upon availability
    if (l_tend_3d_th) then
      allocate( tend_3d_th(current_state%local_grid%size(Z_INDEX),  &
                           current_state%local_grid%size(Y_INDEX),  &
                           current_state%local_grid%size(X_INDEX)   )    )
    endif
    if (l_tend_3d_qv) then
      allocate( tend_3d_qv(current_state%local_grid%size(Z_INDEX),  &
                           current_state%local_grid%size(Y_INDEX),  &
                           current_state%local_grid%size(X_INDEX)   )    )
    endif
    if (l_tend_3d_ql) then
      allocate( tend_3d_ql(current_state%local_grid%size(Z_INDEX),  &
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
    if (l_tend_pr_tot_tabs) then
      allocate( tend_pr_tot_tabs(current_state%local_grid%size(Z_INDEX)) )
    endif

  end subroutine initialisation_callback  


  subroutine finalisation_callback(current_state)
    type(model_state_type), target, intent(inout) :: current_state

    if (allocated(tend_3d_th)) deallocate(tend_3d_th)
    if (allocated(tend_3d_qv)) deallocate(tend_3d_qv)
    if (allocated(tend_3d_ql)) deallocate(tend_3d_ql)
    if (allocated(tend_3d_tabs)) deallocate(tend_3d_tabs)

    if (allocated(tend_pr_tot_th)) deallocate(tend_pr_tot_th)
    if (allocated(tend_pr_tot_qv)) deallocate(tend_pr_tot_qv)
    if (allocated(tend_pr_tot_ql)) deallocate(tend_pr_tot_ql)
    if (allocated(tend_pr_tot_tabs)) deallocate(tend_pr_tot_tabs)

  end subroutine finalisation_callback


  !> Called for each column per timestep this will apply a forcing term 
  !> to the aerosol fields
  !! @param current_state The current model state
  subroutine timestep_callback(current_state)
    type(model_state_type), target, intent(inout) :: current_state

    real(DEFAULT_PRECISION) :: TdegK   ! Temperature in Kelvin
    real(DEFAULT_PRECISION) :: Pmb     ! Pressure in mb
    real(DEFAULT_PRECISION) :: exner   ! Exner pressure
    real(DEFAULT_PRECISION) :: one_over_exner ! Reciprocal of Exner pressure
    real(DEFAULT_PRECISION) :: qv,qc   ! Shorthand for vapour and cloud mass mixing ratio
    real(DEFAULT_PRECISION) :: qs      ! Saturation mixing ratio
    real(DEFAULT_PRECISION) :: dqsdT   ! Rate of change of qs with temperature
    real(DEFAULT_PRECISION) :: qsatfac ! Multiplicative factor
    real(DEFAULT_PRECISION) :: dmass   ! Mass transfer mixing ratio
    
    integer :: k          ! Loop counter
    integer :: icol, jcol ! Shorthand column indices

    real(DEFAULT_PRECISION) :: dtm  ! Local timestep variable
    integer :: target_x_index, target_y_index
    logical :: calculate_diagnostics

    calculate_diagnostics = current_state%diagnostic_sample_timestep

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
      if (l_tend_pr_tot_tabs) then
        tend_pr_tot_tabs(:)=0.0_DEFAULT_PRECISION
      endif
    endif  ! zero totals


    if (current_state%halo_column) return


    dtm = current_state%dtm*2.0
    if (current_state%field_stepping == FORWARD_STEPPING) dtm=current_state%dtm! Should this be revised to scalar_stepping

    icol=current_state%column_local_x
    jcol=current_state%column_local_y
    target_y_index=jcol-current_state%local_grid%halo_size(Y_INDEX)
    target_x_index=icol-current_state%local_grid%halo_size(X_INDEX)

    if (calculate_diagnostics) call save_precomponent_tendencies(current_state, icol, jcol, target_x_index, target_y_index)

    do k=2,k_cloudmax

      exner = current_state%global_grid%configuration%vertical%rprefrcp(k)
      one_over_exner = current_state%global_grid%configuration%vertical%prefrcp(k)
      Pmb   = (current_state%global_grid%configuration%vertical%prefn(k)/100.) 

      if (current_state%field_stepping == FORWARD_STEPPING) then  ! Should this be revised to scalar_stepping
        qv    = current_state%q(iqv)%data(k, jcol, icol) + current_state%sq(iqv)%data(k, jcol, icol)*dtm
        qc    = current_state%q(iql)%data(k, jcol, icol) + current_state%sq(iql)%data(k, jcol, icol)*dtm
        TdegK = (current_state%th%data(k, jcol, icol) + current_state%sth%data(k, jcol, icol)*dtm &
           + current_state%global_grid%configuration%vertical%thref(k))*exner
      else
        qv    = current_state%zq(iqv)%data(k, jcol, icol) + current_state%sq(iqv)%data(k, jcol, icol)*dtm
        qc    = current_state%zq(iql)%data(k, jcol, icol) + current_state%sq(iql)%data(k, jcol, icol)*dtm
        TdegK = (current_state%zth%data(k, jcol, icol) + current_state%sth%data(k, jcol, icol)*dtm &
           + current_state%global_grid%configuration%vertical%thref(k))*exner
      end if
      ! Calculate the cloud/vapour increments

      qs = qsaturation(TdegK, Pmb)

      if (qv > qs .or. qc >0.0)then
        dqsdT = dqwsatdt(qs, TdegK)
        
        qsatfac = 1.0/(1.0 + rlvap_over_cp*dqsdT)
      
        dmass = MAX (-qc,(qv-qs)*qsatfac)/dtm
      
        current_state%sq(iqv)%data(k, jcol, icol) = current_state%sq(iqv)%data(k, jcol, icol) - dmass
        current_state%sq(iql)%data(k, jcol, icol) = current_state%sq(iql)%data(k, jcol, icol) + dmass

        current_state%sth%data(k, jcol, icol) = current_state%sth%data(k, jcol, icol) &
           + rlvap_over_cp*dmass*one_over_exner

      end if

    end do

    ! If there's any cloud above then evaporate it
    do k=k_cloudmax+1, current_state%local_grid%size(Z_INDEX)
      if (current_state%scalar_stepping == FORWARD_STEPPING) then
        qv    = current_state%q(iqv)%data(k, jcol, icol) + current_state%sq(iqv)%data(k, jcol, icol)*dtm
        qc    = current_state%q(iql)%data(k, jcol, icol) + current_state%sq(iql)%data(k, jcol, icol)*dtm
        TdegK = (current_state%th%data(k, jcol, icol) + current_state%sth%data(k, jcol, icol)*dtm &
           + current_state%global_grid%configuration%vertical%thref(k))*exner
      else
        qv    = current_state%zq(iqv)%data(k, jcol, icol) + current_state%sq(iqv)%data(k, jcol, icol)*dtm
        qc    = current_state%zq(iql)%data(k, jcol, icol) + current_state%sq(iql)%data(k, jcol, icol)*dtm
        TdegK = (current_state%zth%data(k, jcol, icol) + current_state%sth%data(k, jcol, icol)*dtm &
           + current_state%global_grid%configuration%vertical%thref(k))*exner
      end if
      if (qc >0.0)then
        dmass = -qc/dtm
      
        current_state%sq(iqv)%data(k, jcol, icol) = current_state%sq(iqv)%data(k, jcol, icol) - dmass
        current_state%sq(iql)%data(k, jcol, icol) = current_state%sq(iql)%data(k, jcol, icol) + dmass

        current_state%sth%data(k, jcol, icol) = current_state%sth%data(k, jcol, icol) &
           + rlvap_over_cp*dmass*one_over_exner

      end if
    end do

    if (calculate_diagnostics) call compute_component_tendencies(current_state, icol, jcol, target_x_index, target_y_index)

  end subroutine timestep_callback


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

end module simplecloud_mod
