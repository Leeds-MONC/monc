!> Conditionally averaged diagnostics, Part 2 of 2.
!  Using the totalled values from Part 1, averages the results over the horizontal domain,
!  creating profiles of the requested diagnostics under the requested conditions and .not. conditions.
!  The "area" diagnostic is a special case, recording the fraction of the domain meeting the associated
!  condition and .not. condition.
module conditional_diagnostics_whole_mod
  use monc_component_mod, only : component_descriptor_type
  use state_mod, only : model_state_type

  use conditional_diagnostics_column_mod, only : CondDiags_tot, ncond, ndiag, gpts_total, requested_area
  use grids_mod, only : Z_INDEX
  use datadefn_mod, only : PRECISION_TYPE, DEFAULT_PRECISION
  use mpi, only : MPI_SUM, MPI_IN_PLACE, MPI_INT, MPI_REAL, MPI_DOUBLE, MPI_Comm
  use missing_data_mod, only: rmdi
  use optionsdatabase_mod, only : options_get_integer

  implicit none

#ifndef TEST_MODE
  private
#endif

  public conditional_diagnostics_whole_get_descriptor

contains

  !> Provides registry information for the component
  !! @returns The component descriptor that describes this component
  type(component_descriptor_type) function conditional_diagnostics_whole_get_descriptor()
    conditional_diagnostics_whole_get_descriptor%name="conditional_diagnostics_whole"
    conditional_diagnostics_whole_get_descriptor%version=0.1
    conditional_diagnostics_whole_get_descriptor%initialisation=>initialisation_callback
    conditional_diagnostics_whole_get_descriptor%timestep=>timestep_callback
    conditional_diagnostics_whole_get_descriptor%finalisation=>finalisation_callback
  end function conditional_diagnostics_whole_get_descriptor

  !> Initialisation hook: currently doesn't need to do anything
  !! @param current_state The current model state
  subroutine initialisation_callback(current_state)
    type(model_state_type), target, intent(inout) :: current_state

  end subroutine initialisation_callback


  !> The timestep hook will perform averaging of the conditional diagnostics
  !! @param current_state The current model state
  subroutine timestep_callback(current_state)
    type(model_state_type), target, intent(inout) :: current_state
    integer :: k,cnc,dnc,ierr
    real(kind=DEFAULT_PRECISION) :: temp
    logical :: calculate_diagnostics

    calculate_diagnostics = current_state%diagnostic_sample_timestep

    !> Decide if conditions are appropriate to proceed with calculations
    if (.not. calculate_diagnostics) return

    !> Sum conditional diagnostics total array (horizontally), placing the result on process 0
    !! Reduction call on process 0 requires special MPI_IN_PLACE handling
    if (current_state%parallel%my_rank == 0) then
      call mpi_reduce(MPI_IN_PLACE , CondDiags_tot, ncond*2*ndiag*current_state%local_grid%size(Z_INDEX), &
                      PRECISION_TYPE, MPI_SUM, 0, current_state%parallel%monc_communicator, ierr)
    else
      call mpi_reduce(CondDiags_tot, CondDiags_tot, ncond*2*ndiag*current_state%local_grid%size(Z_INDEX), &
                      PRECISION_TYPE, MPI_SUM, 0, current_state%parallel%monc_communicator, ierr)
    end if

    !> Average summed diagnostics over the domain by dividing the total diagnostic for each condition 
    !! by the total number of points for the associated conditions.
    !! This is NOT done for the area diagnostic, identified by the requested_area position in the array.
    !! Note: CondDiags_tot(k, ncond, ndiag)
    do dnc = 1, ndiag
      if (dnc /= requested_area) then
        do cnc = 1, ncond
          do k = 2, current_state%local_grid%size(Z_INDEX) - 1
            temp = CondDiags_tot(k,cnc,requested_area)
            if (temp .gt. 0) then
              CondDiags_tot(k,cnc,dnc) = CondDiags_tot(k,cnc,dnc) / temp
            else
              CondDiags_tot(k,cnc,dnc) = rmdi
            end if
          end do ! k
        end do ! cnc over ncond*2
      end if ! check requested_area
    end do ! dnc over ndiag

    !> Convert the area total number of points for each condition to fraction of the horizontal domain.
    CondDiags_tot(:,:,requested_area) = CondDiags_tot(:,:,requested_area) / gpts_total

    !> Apply missing data mask to top/bottom
    CondDiags_tot(1,:,:) = rmdi
    CondDiags_tot(current_state%local_grid%size(Z_INDEX),:,:) = rmdi

    !> Since the xml handling of CondDiags_tot will perform a sum over processes, divide by the number
    !! of proceses.
    CondDiags_tot = CondDiags_tot / current_state%parallel%processes

    !> Broadcast the process-fractional solution to all processes.
    call mpi_bcast(CondDiags_tot,  ncond*2*ndiag*current_state%local_grid%size(Z_INDEX), &
                   PRECISION_TYPE, 0, current_state%parallel%monc_communicator, ierr)

  end subroutine timestep_callback

  
  !> Called on termination: currently doesn't need to do anything
  !! @param current_state The current model state
  subroutine finalisation_callback(current_state)
    type(model_state_type), target, intent(inout) :: current_state

  end subroutine finalisation_callback

end module conditional_diagnostics_whole_mod
