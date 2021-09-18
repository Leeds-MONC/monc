!> Simple exponential scheme to calculate the longwave radiation associated
!> with cloud. The scheme is based on the methods used in GASS intercomparison
!> cases, e.g. DYCOMS, ISDAC.
!>
!> This scheme depends on exp_lw, lwtop_in, lwbase_in, which are
!> set in the config file.
module lwrad_exponential_mod
  use datadefn_mod, only : DEFAULT_PRECISION
  use monc_component_mod, only : COMPONENT_ARRAY_FIELD_TYPE, &
      COMPONENT_DOUBLE_DATA_TYPE, component_descriptor_type, &
      component_field_value_type, component_field_information_type
  use optionsdatabase_mod, only : options_get_real
  use state_mod, only : FORWARD_STEPPING, model_state_type
  use grids_mod, only : Z_INDEX, Y_INDEX, X_INDEX
  use science_constants_mod, only : cp
  use q_indices_mod, only: get_q_index, standard_q_names
  use registry_mod, only : is_component_enabled
  use logging_mod, only : LOG_ERROR, log_master_log
  use optionsdatabase_mod, only : options_get_integer

  implicit none

  ! Indices for cloud
  integer :: iql
  ! Index for the top of the domain
  integer :: k_top, x_local, y_local

  ! local column longwave flux variables, calculated using prescribed values from
  ! the global/user_config file
  real(DEFAULT_PRECISION), allocatable ::  lwrad_flux_top(:), lwrad_flux_base(:), lwrad_flux(:)
  ! local column liquid water content
  real(DEFAULT_PRECISION), allocatable ::  qc_col(:)
  ! local density factor and raidiation factor used in conversions
  real(DEFAULT_PRECISION), allocatable ::  density_factor(:), radiation_factor(:)
  ! declare radiative heating rate for longwave. This is declared as a 3D array so
  ! that the heating rates can be stored for diagnostics or use when the radiation
  ! timestep is longer than the model timestep
  real(DEFAULT_PRECISION), allocatable :: sth_lw(:,:,:)

  ! declare radiative flux variables that are read for global or user_config
  real(DEFAULT_PRECISION) :: longwave_exp_decay, cltop_longwave_flux, clbase_longwave_flux

  ! Local tendency diagnostic variables for this component
  ! 3D tendency fields and logicals for their use
  real(kind=DEFAULT_PRECISION), dimension(:,:,:), allocatable :: tend_3d_th, tend_3d_tabs
  logical :: l_tend_3d_th, l_tend_3d_tabs
  ! Local mean tendency profile fields and logicals for their use
  real(kind=DEFAULT_PRECISION), dimension(:), allocatable ::     tend_pr_tot_th, tend_pr_tot_tabs
  logical :: l_tend_pr_tot_th, l_tend_pr_tot_tabs

  public lwrad_exponential_get_descriptor
contains

  !> Provides the descriptor back to the caller and is used in component registration
  !! @returns The termination check component descriptor
  type(component_descriptor_type) function lwrad_exponential_get_descriptor()
    lwrad_exponential_get_descriptor%name="lwrad_exponential"
    lwrad_exponential_get_descriptor%version=0.1
    lwrad_exponential_get_descriptor%initialisation=>initialisation_callback
    lwrad_exponential_get_descriptor%timestep=>timestep_callback
    lwrad_exponential_get_descriptor%finalisation=>finalisation_callback

    lwrad_exponential_get_descriptor%field_value_retrieval=>field_value_retrieval_callback
    lwrad_exponential_get_descriptor%field_information_retrieval=>field_information_retrieval_callback
    allocate(lwrad_exponential_get_descriptor%published_fields(2+2))

    lwrad_exponential_get_descriptor%published_fields(1)= "tend_th_lwrad_exponential_3d_local"
    lwrad_exponential_get_descriptor%published_fields(2)= "tend_tabs_lwrad_exponential_3d_local"

    lwrad_exponential_get_descriptor%published_fields(2+1)= "tend_th_lwrad_exponential_profile_total_local"
    lwrad_exponential_get_descriptor%published_fields(2+2)= "tend_tabs_lwrad_exponential_profile_total_local"
  end function lwrad_exponential_get_descriptor

  !> The initialisation callback sets up the prescribed longwave fluxes and the
  !> exponential decay factor
  !! @param current_state The current model state
  subroutine initialisation_callback(current_state)
    type(model_state_type), target, intent(inout) :: current_state

    integer :: kkp, k ! look counter

    
    k_top = current_state%local_grid%size(Z_INDEX) + current_state%local_grid%halo_size(Z_INDEX) * 2
    x_local = current_state%local_grid%size(X_INDEX) + current_state%local_grid%halo_size(X_INDEX) * 2
    y_local = current_state%local_grid%size(Y_INDEX) + current_state%local_grid%halo_size(X_INDEX) * 2

    if (is_component_enabled(current_state%options_database, "socrates_couple")) then
       call log_master_log &
            (LOG_ERROR, "Socrates and lwrad_exponential both enabled, switch off one in config - STOP")
    endif
    
    allocate(lwrad_flux_top(k_top))
    allocate(lwrad_flux_base(k_top))
    allocate(lwrad_flux(k_top))
    allocate(qc_col(k_top))
    allocate(density_factor(k_top))
    allocate(radiation_factor(k_top))
    ! NOTE: this may be the wrong declaration as sth_lw may need to declared on the whole domain
    allocate(sth_lw(k_top,y_local,x_local))

    iql=get_q_index(standard_q_names%CLOUD_LIQUID_MASS, 'lwrad_exponential')

    longwave_exp_decay=options_get_real(current_state%options_database, "longwave_exp_decay")
    cltop_longwave_flux=options_get_real(current_state%options_database, "cltop_longwave_flux")
    clbase_longwave_flux=options_get_real(current_state%options_database, "clbase_longwave_flux")

    density_factor(:) = current_state%global_grid%configuration%vertical%rhon(:)*                  &
         current_state%global_grid%configuration%vertical%dz(:)

    radiation_factor(2:k_top) = 1.0/(density_factor(2:k_top)*cp)
    radiation_factor(1) = radiation_factor(2)

    ! Set tendency diagnostic logicals based on availability
    ! Need to use 3d tendencies to compute the profiles, so they will be allocated
    !      in the case where profiles are available
    l_tend_pr_tot_th  = current_state%th%active 
    l_tend_pr_tot_tabs  = l_tend_pr_tot_th

    l_tend_3d_th  = current_state%th%active .or. l_tend_pr_tot_th
    l_tend_3d_tabs = l_tend_3d_th

    ! Allocate 3d tendency fields upon availability
    if (l_tend_3d_th) then
      allocate( tend_3d_th(current_state%local_grid%size(Z_INDEX),  &
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
    if (l_tend_pr_tot_tabs) then
      allocate( tend_pr_tot_tabs(current_state%local_grid%size(Z_INDEX)) )
    endif

  end subroutine initialisation_callback

  !> Called for each column per timestep this will apply a forcing term
  !> to the aerosol fields
  !! @param current_state The current model state
  subroutine timestep_callback(current_state)
    type(model_state_type), target, intent(inout) :: current_state

    integer :: k          ! Loop counter
    integer :: icol, jcol ! Shorthand column indices
    integer :: target_x_index, target_y_index ! Column indices with halos subtracted
    real(DEFAULT_PRECISION) :: dtm  ! Local timestep variable

    if (current_state%halo_column) return

    dtm = current_state%dtm*2.0
    if (current_state%field_stepping == FORWARD_STEPPING) dtm=current_state%dtm! Should this be revised to scalar_stepping

    icol=current_state%column_local_x
    jcol=current_state%column_local_y
    target_y_index = jcol - current_state%local_grid%halo_size(Y_INDEX)
    target_x_index = icol - current_state%local_grid%halo_size(X_INDEX)

    ! Zero profile tendency totals on first instance in the sum
    if (current_state%first_nonhalo_timestep_column) then
      if (l_tend_pr_tot_th) then
        tend_pr_tot_th(:) = 0.0_DEFAULT_PRECISION
      endif
      if (l_tend_pr_tot_tabs) then
        tend_pr_tot_tabs(:) = 0.0_DEFAULT_PRECISION
      endif
    endif  ! zero totals

    ! set the column liquid water content
    if (current_state%field_stepping == FORWARD_STEPPING) then  ! Should this be revised to scalar_stepping
       qc_col(:) = current_state%q(iql)%data(:, jcol, icol) + current_state%sq(iql)%data(:, jcol, icol)*dtm
    else
       qc_col(:)= current_state%zq(iql)%data(:, jcol, icol) + current_state%sq(iql)%data(:, jcol, icol)*dtm
    end if

    ! initialise the flux top and base to 0.0
    lwrad_flux_top(:) = 0.0
    lwrad_flux_base(:) = 0.0

    ! First workout cloud top cooling
    lwrad_flux_top(k_top)=cltop_longwave_flux !units Wm-2

    do k = k_top-1, 1, -1
       if (qc_col(k+1) > 1.e-10) then
          lwrad_flux_top(k) = lwrad_flux_top(k+1)*                                    &
               exp(-qc_col(k+1)*density_factor(k+1)*longwave_exp_decay)
       else
          lwrad_flux_top(k) = lwrad_flux_top(k+1)
       endif
    enddo

    ! Second workout the cloud base warming
    do k = 2, k_top
       if (qc_col(k) > 1.e-10) then
          lwrad_flux_base(k) = lwrad_flux_base(k-1)*                                  &
               exp(-qc_col(k)*density_factor(k)*longwave_exp_decay)
       else
          lwrad_flux_base(k) = lwrad_flux_base(k-1)
       endif
    enddo

    ! Next, sum the fluxes
    do k = 1, k_top
       lwrad_flux(k) = lwrad_flux_base(k) + lwrad_flux_top(k)
    enddo

    ! workout radiative heating rate, and stored on the processors local grid - is this correct??
    ! or should the declaration be on the global grid?
    do k = 2, k_top
       sth_lw(k, jcol, icol) = -(lwrad_flux(k) - lwrad_flux(k-1))*radiation_factor(k)
    enddo

    ! update the current_state sth
    current_state%sth%data(:, jcol, icol) = current_state%sth%data(:, jcol, icol) + sth_lw(:, jcol, icol)

    ! Zero profile tendency totals on first instance in the sum
    if (current_state%first_timestep_column) then
      tend_pr_tot_th(:) = 0.0_DEFAULT_PRECISION
      tend_pr_tot_tabs(:) = 0.0_DEFAULT_PRECISION
    endif  ! zero totals

    if (current_state%diagnostic_sample_timestep .and. .not. current_state%halo_column) then
      call compute_component_tendencies(current_state, icol, jcol, target_x_index, target_y_index)
    end if

  end subroutine timestep_callback

  subroutine finalisation_callback(current_state)

    type(model_state_type), target, intent(inout) :: current_state
    deallocate(lwrad_flux_top,lwrad_flux_base, lwrad_flux, density_factor, radiation_factor, sth_lw)

    if (allocated(tend_3d_th)) deallocate(tend_3d_th)
    if (allocated(tend_3d_tabs)) deallocate(tend_3d_tabs)

    if (allocated(tend_pr_tot_th)) deallocate(tend_pr_tot_th)
    if (allocated(tend_pr_tot_tabs)) deallocate(tend_pr_tot_tabs)

  end subroutine finalisation_callback

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
      tend_3d_th(:,tyn,txn) = sth_lw(:,cyn,cxn)
    endif
    if (l_tend_3d_tabs) then
      tend_3d_tabs(:,tyn,txn) = sth_lw(:,cyn,cxn) * current_state%global_grid%configuration%vertical%rprefrcp(:)
    endif

   ! Add local tendency fields to the profile total
    if (l_tend_pr_tot_th) then
      tend_pr_tot_th(:)=tend_pr_tot_th(:) + tend_3d_th(:,tyn,txn)
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
    strcomp=INDEX(name, "lwrad_exponential_3d_local")
    if (strcomp .ne. 0) then
      field_information%field_type=COMPONENT_ARRAY_FIELD_TYPE
      field_information%number_dimensions=3
      field_information%dimension_sizes(1)=current_state%local_grid%size(Z_INDEX)
      field_information%dimension_sizes(2)=current_state%local_grid%size(Y_INDEX)
      field_information%dimension_sizes(3)=current_state%local_grid%size(X_INDEX)
      field_information%data_type=COMPONENT_DOUBLE_DATA_TYPE

      if (name .eq. "tend_th_lwrad_exponential_3d_local") then
        field_information%enabled=l_tend_3d_th
      else if (name .eq. "tend_tabs_lwrad_exponential_3d_local") then
        field_information%enabled=l_tend_3d_tabs
      else
        field_information%enabled=.true.
      end if

    end if !end 3d check

    ! Field information for profiles
    strcomp=INDEX(name, "lwrad_exponential_profile_total_local")
    if (strcomp .ne. 0) then
      field_information%field_type=COMPONENT_ARRAY_FIELD_TYPE
      field_information%number_dimensions=1
      field_information%dimension_sizes(1)=current_state%local_grid%size(Z_INDEX)
      field_information%data_type=COMPONENT_DOUBLE_DATA_TYPE

      if (name .eq. "tend_th_lwrad_exponential_profile_total_local") then
        field_information%enabled=l_tend_pr_tot_th
      else if (name .eq. "tend_tabs_lwrad_exponential_profile_total_local") then
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
    if (name .eq. "tend_th_lwrad_exponential_3d_local" .and. allocated(tend_3d_th)) then
      call set_published_field_value(field_value, real_3d_field=tend_3d_th)
    else if (name .eq. "tend_tabs_lwrad_exponential_3d_local" .and. allocated(tend_3d_tabs)) then
      call set_published_field_value(field_value, real_3d_field=tend_3d_tabs)

    ! Profile Tendency Fields
    else if (name .eq. "tend_th_lwrad_exponential_profile_total_local" .and. allocated(tend_pr_tot_th)) then
      call set_published_field_value(field_value, real_1d_field=tend_pr_tot_th)
    else if (name .eq. "tend_tabs_lwrad_exponential_profile_total_local" .and. allocated(tend_pr_tot_tabs)) then
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


end module lwrad_exponential_mod
