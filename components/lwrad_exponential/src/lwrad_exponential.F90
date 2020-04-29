!> Simple exponential scheme to calculate the longwave radiation associated
!> with cloud. The scheme is based on the methods used in GASS intercomparison
!> cases, e.g. DYCOMS, ISDAC.
!>
!> This scheme depends on exp_lw, lwtop_in, lwbase_in, which are
!> set in the config file.
module lwrad_exponential_mod
  use datadefn_mod, only : DEFAULT_PRECISION
  use monc_component_mod, only : component_descriptor_type
  use optionsdatabase_mod, only : options_get_real
  use state_mod, only : FORWARD_STEPPING, model_state_type
  use grids_mod, only : Z_INDEX, Y_INDEX, X_INDEX
  use science_constants_mod, only : cp
  use q_indices_mod, only: get_q_index, standard_q_names
  use registry_mod, only : is_component_enabled
  use logging_mod, only : LOG_ERROR, log_master_log

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

  end subroutine initialisation_callback

  !> Called for each column per timestep this will apply a forcing term
  !> to the aerosol fields
  !! @param current_state The current model state
  subroutine timestep_callback(current_state)
    type(model_state_type), target, intent(inout) :: current_state

    integer :: k          ! Loop counter
    integer :: icol, jcol ! Shorthand column indices

    real(DEFAULT_PRECISION) :: dtm  ! Local timestep variable

    if (current_state%halo_column) return

    dtm = current_state%dtm*2.0
    if (current_state%field_stepping == FORWARD_STEPPING) dtm=current_state%dtm! Should this be revised to scalar_stepping

    icol=current_state%column_local_x
    jcol=current_state%column_local_y

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

  end subroutine timestep_callback

  subroutine finalisation_callback(current_state)

    type(model_state_type), target, intent(inout) :: current_state
    deallocate(lwrad_flux_top,lwrad_flux_base, lwrad_flux, density_factor, radiation_factor, sth_lw)

  end subroutine finalisation_callback

end module lwrad_exponential_mod
