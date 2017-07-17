!> A very simple saturation adjustment scheme without any microphysics
module simplecloud_mod
  use datadefn_mod, only : DEFAULT_PRECISION
  use monc_component_mod, only : component_descriptor_type
  use optionsdatabase_mod, only : options_get_real
  use state_mod, only : FORWARD_STEPPING, model_state_type
  use grids_mod, only : Z_INDEX
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

    public simplecloud_get_descriptor
contains

  !> Provides the descriptor back to the caller and is used in component registration
  !! @returns The termination check component descriptor
  type(component_descriptor_type) function simplecloud_get_descriptor()
    simplecloud_get_descriptor%name="simplecloud"
    simplecloud_get_descriptor%version=0.1
    simplecloud_get_descriptor%initialisation=>initialisation_callback
    simplecloud_get_descriptor%timestep=>timestep_callback
  end function simplecloud_get_descriptor

  !> The initialisation callback sets up the moisture fields
  !! @param current_state The current model state
  subroutine initialisation_callback(current_state)
    type(model_state_type), target, intent(inout) :: current_state

    integer :: k ! look counter

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

  end subroutine initialisation_callback  

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

    if (current_state%halo_column) return

    dtm = current_state%dtm*2.0
    if (current_state%field_stepping == FORWARD_STEPPING) dtm=current_state%dtm! Should this be revised to scalar_stepping

    icol=current_state%column_local_x
    jcol=current_state%column_local_y

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

  end subroutine timestep_callback

end module simplecloud_mod
