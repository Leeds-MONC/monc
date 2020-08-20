!> This contains the CFL test. It will perform the local advective CFL and Galilean transfromation calculations,
!! compute the global values of these, then check the cfl criterion and determine an absolute new dtm. Depending upon the
!! maximum and increment values this is then smoothed into a new dtm value. The value of dtm is physically set to this new
!! value at the start of the next timestep. 
module cfltest_mod
  use datadefn_mod, only : DEFAULT_PRECISION, PRECISION_TYPE
  use monc_component_mod, only : component_descriptor_type
  use state_mod, only : model_state_type, parallel_state_type
  use collections_mod, only : map_type
  use logging_mod, only : LOG_WARN, LOG_DEBUG, LOG_ERROR, LOG_INFO, &
                          log_log, log_get_logging_level, log_master_newline, &
                          log_master_log, log_is_master
  use conversions_mod, only : conv_to_string
  use optionsdatabase_mod, only : options_get_integer, options_get_real, options_get_logical
  use grids_mod, only : Z_INDEX, Y_INDEX, X_INDEX
  use mpi, only : MPI_MAX, MPI_MIN  
  implicit none

#ifndef TEST_MODE
  private
#endif

  !! Configuration options - all are optional and have default values
  real(kind=DEFAULT_PRECISION) :: tollerance, cvismax, cvelmax, dtmmax, dtmmin, rincmax
  logical l_monitor_cfl, l_constant_dtm

  public cfltest_get_descriptor
contains

  !> Provides the descriptor back to the caller and is used in component registration
  !! @returns The termination check component descriptor
  type(component_descriptor_type) function cfltest_get_descriptor()
    cfltest_get_descriptor%name="cfltest"
    cfltest_get_descriptor%version=0.1
    cfltest_get_descriptor%initialisation=>initialisation_callback
    cfltest_get_descriptor%timestep=>timestep_callback
  end function cfltest_get_descriptor

  !> Called at initialisation, will read in configuration and use either configured or default values. 
  !! @param current_state The current model state
  subroutine initialisation_callback(current_state)
    type(model_state_type), intent(inout), target :: current_state

    current_state%cfl_frequency=options_get_integer(current_state%options_database, "cfl_frequency")
    tollerance=options_get_real(current_state%options_database, "cfl_tollerance")
    cvismax=options_get_real(current_state%options_database, "cfl_cvismax")
    cvelmax=options_get_real(current_state%options_database, "cfl_cvelmax")
    dtmmax=options_get_real(current_state%options_database, "cfl_dtmmax")
    dtmmin=options_get_real(current_state%options_database, "cfl_dtmmin")
    rincmax=options_get_real(current_state%options_database, "cfl_rincmax")

    l_monitor_cfl = options_get_logical(current_state%options_database,"cfl_monitor")
    l_constant_dtm = options_get_logical(current_state%options_database,"l_constant_dtm")

    allocate(current_state%abswmax(current_state%local_grid%local_domain_end_index(Z_INDEX)))
  end subroutine initialisation_callback

  !> Called at each timestep, this will only do the CFL computation every nncfl timesteps (or every timestep up to nncfl) but
  !! will ratchet up to the absolute (target) dtm as needed.
  !! @param current_state The current model state
  subroutine timestep_callback(current_state)
    type(model_state_type), intent(inout), target :: current_state

    real(kind=DEFAULT_PRECISION) :: cfl_number

    if (current_state%normal_step) then

      if ((mod(current_state%timestep, current_state%cfl_frequency) == 1 .or. &
           current_state%timestep-current_state%start_timestep .le. current_state%cfl_frequency) &
          .or. current_state%timestep .ge. (current_state%last_cfl_timestep + current_state%cfl_frequency)) then

        current_state%last_cfl_timestep = current_state%timestep
        current_state%cvel=0.0_DEFAULT_PRECISION
        current_state%cvel_x=0.0_DEFAULT_PRECISION
        current_state%cvel_y=0.0_DEFAULT_PRECISION
        current_state%cvel_z=0.0_DEFAULT_PRECISION
  
        call perform_cfl_and_galilean_transformation_calculation(current_state)
  
        current_state%cvel=(current_state%cvel_x*current_state%global_grid%configuration%horizontal%cx+current_state%cvel_y*&
             current_state%global_grid%configuration%horizontal%cy+current_state%cvel_z)*current_state%dtm
        current_state%cvis=current_state%cvis*(current_state%dtm * 4)
  
        cfl_number=current_state%cvis/cvismax+current_state%cvel/cvelmax
  
        current_state%absolute_new_dtm=current_state%dtm
        current_state%update_dtm=.false.
        if (cfl_number .gt. 0.0_DEFAULT_PRECISION) then
          if (cfl_number .lt. (1.0_DEFAULT_PRECISION-tollerance) .or. cfl_number .gt. (1.0_DEFAULT_PRECISION+tollerance)) then
            current_state%absolute_new_dtm=current_state%dtm/cfl_number
          end if
        end if
      end if ! evaluate the cfl 

      call update_dtm_based_on_absolute(current_state, cfl_number)

      if (current_state%time_basis .or. current_state%force_output_on_interval) &
         call evaluate_time_basis(current_state) 

    else ! do not evaluate if taking reduced NON-normal_step

      current_state%update_dtm=.false.

    end if ! end check normal_step

    current_state%cvis=0.0_DEFAULT_PRECISION

  end subroutine timestep_callback

  !> Updates the (new) dtm value, which is actioned after time step completion, based upon the absolute value. This is incremented
  !! towards the absolute if that is too large a step, and even if the absolute value has not been updated in this timestep,
  !! this ratcheting will still occur if needed.
  !! @param current_state The current model state
  subroutine update_dtm_based_on_absolute(current_state, cfl_number)
    type(model_state_type), intent(inout), target :: current_state
    real(kind=DEFAULT_PRECISION), intent(in)      :: cfl_number

    if (current_state%dtm .ne. current_state%absolute_new_dtm .and. &
         (current_state%dtm .ne. dtmmax .or. current_state%absolute_new_dtm .lt. dtmmax)) then

      current_state%update_dtm=.true.

      current_state%dtm_new=min(current_state%dtm*(1.0_DEFAULT_PRECISION+rincmax), current_state%absolute_new_dtm, dtmmax)

      if (l_monitor_cfl) then
        call log_master_log(LOG_INFO, " --- CFL Monitoring Information --- ")
        if (l_constant_dtm) &
          call log_master_log(LOG_INFO, " *** l_constant_dtm=.true. - NOT CHANGING ***")
        call log_master_log(LOG_INFO, "dtm changed from "//trim(conv_to_string(current_state%dtm, 5))&
                //" to "//trim(conv_to_string(current_state%dtm_new, 5)))
        if (cfl_number .gt. 0.0) then 
          call log_master_log(LOG_INFO, "cfl_number :  "//trim(conv_to_string(cfl_number))//"  (change divisor)")
          call log_master_log(LOG_INFO, "cvis       :  "//trim(conv_to_string(current_state%cvis)) )
          call log_master_log(LOG_INFO, "cvel       :  "//trim(conv_to_string(current_state%cvel)))
!//", "//trim(conv_to_string(current_state%cvel_x))//", "//trim(conv_to_string(current_state%cvel_y))//", "//trim(conv_to_string(current_state%cvel_z)) )
        else
          call log_master_log(LOG_INFO, "dtm change due to ratcheting only. Target dtm unchanged.")
        end if
        call log_master_log(LOG_INFO, "target dtm :  "//trim(conv_to_string(current_state%absolute_new_dtm)) )
        call log_master_newline()
      end if ! l_monitor_cfl


      !! --- Diagnostic Writing -----------------
      call log_master_log(LOG_DEBUG, "dtm changed from "//trim(conv_to_string(current_state%dtm, 5))//" to "//&
             trim(conv_to_string(current_state%dtm_new, 5)))

      if (current_state%dtm_new .lt. dtmmin .and. .not. l_constant_dtm) then
        call log_master_log(LOG_ERROR, "Timestep too small, dtmnew="//&
                                       trim(conv_to_string(current_state%dtm_new, 5))//&
                                       ", dtmmin="//trim(conv_to_string(dtmmin, 5)))
      end if

      if (current_state%dtm_new .lt. current_state%dtm .and. l_constant_dtm) then
        call log_master_log(LOG_WARN, "The CFL check would like to reduce the model timestep "//&
                                "to a value below the specified constant dtm: dtmnew="//&
                                trim(conv_to_string(current_state%dtm_new, 5))//&
                                ", constant dtm="//trim(conv_to_string(current_state%dtm, 5)))
      end if
    end if

    if (l_constant_dtm) then
      current_state%dtm = options_get_real(current_state%options_database, "dtm")
      current_state%dtm_new = current_state%dtm
      current_state%update_dtm=.false.
    end if
  end subroutine update_dtm_based_on_absolute  

  !> Performs the CFL and Galilean transformation calculations. First locally and then will determine the global value
  !! of each calculation. If U, V or W are not active then these are set to 0 and the calculation use these values
  !! @param current_state The current model state
  subroutine perform_cfl_and_galilean_transformation_calculation(current_state)
    type(model_state_type), intent(inout), target :: current_state

    integer :: k
    real(kind=DEFAULT_PRECISION) :: global_zumin, global_zumax, global_zvmin, &
         global_zvmax, global_cvel_z, global_cvis

#ifdef U_ACTIVE
    current_state%local_zumin=current_state%local_zumin+current_state%ugal               ! _undo Gal-trfm                        
    current_state%local_zumax=current_state%local_zumax+current_state%ugal
#else
    current_state%local_zumin=0.0_DEFAULT_PRECISION
    current_state%local_zumax=0.0_DEFAULT_PRECISION
#endif
#ifdef V_ACTIVE
    current_state%local_zvmin=current_state%local_zvmin+current_state%vgal               ! _undo Gal-trfm                        
    current_state%local_zvmax=current_state%local_zvmax+current_state%vgal
#else
    current_state%local_zvmin=0.0_DEFAULT_PRECISION
    current_state%local_zvmax=0.0_DEFAULT_PRECISION
#endif
#ifdef W_ACTIVE
    current_state%local_cvel_z=current_state%cvel_z
    do k=2,current_state%local_grid%local_domain_end_index(Z_INDEX)-1
      !  CVELZ will be multiplied by DTM in TESTCFL
      current_state%local_cvel_z=max(current_state%local_cvel_z, &
           current_state%abswmax(k)*current_state%global_grid%configuration%vertical%rdzn(k+1))
    end do
#else
    current_state%local_cvel_z=0.0_DEFAULT_PRECISION
#endif
    call get_global_values(current_state%local_zumin, current_state%local_zumax, current_state%local_zvmin, &
         current_state%local_zvmax, current_state%local_cvel_z, current_state%cvis, &
         global_zumin, global_zumax, global_zvmin, global_zvmax, global_cvel_z, global_cvis, current_state%parallel)

    if (current_state%galilean_transformation) then
      if (.not.current_state%fix_ugal)current_state%ugal=0.5_DEFAULT_PRECISION*(global_zumin+global_zumax)
      if (.not.current_state%fix_vgal)current_state%vgal=0.5_DEFAULT_PRECISION*(global_zvmin+global_zvmax)
    else
      current_state%ugal=0.0_DEFAULT_PRECISION
      current_state%vgal=0.0_DEFAULT_PRECISION
    end if
    current_state%cvel_z=global_cvel_z
    current_state%cvel_x=max(abs(global_zumax-current_state%ugal), abs(global_zumin-current_state%ugal))
    current_state%cvel_y=max(abs(global_zvmax-current_state%vgal), abs(global_zvmin-current_state%vgal))
    current_state%cvis=global_cvis
  end subroutine perform_cfl_and_galilean_transformation_calculation 

  !> Gets the global reduction values based upon the local contributions of CFL and Galilean transformations provided
  !! @param local_zumin Local contribution to ZU minimum
  !! @param local_zumax Local contribution to ZU maximum
  !! @param local_zvmin Local contribution to ZV minimum
  !! @param local_zvmax Local contribution to ZV maximum
  !! @param local_cvel_z Local contribution to component of Courant number in z
  !! @param local_cvis Local contribution to viscous Courant number
  !! @param global_zumin Globally reduced value of ZU minimum
  !! @param global_zumax Globally reduced value of ZU maximum
  !! @param global_zvmin Globally reduced value of ZV minimum
  !! @param global_zvmax Globally reduced value of ZV maximum
  !! @param global_cvel_z Globally reduced value of contribution to component of Courant number in z
  !! @param global_cvis Globally reduced value of contribution to viscous Courant number
  !! @param parallel_state The parallel state
  subroutine get_global_values(local_zumin, local_zumax, local_zvmin, local_zvmax, local_cvel_z, local_cvis, &
       global_zumin, global_zumax, global_zvmin, global_zvmax, global_cvel_z, global_cvis, parallel_state)
    type(parallel_state_type), intent(inout) :: parallel_state
    real(kind=DEFAULT_PRECISION), intent(in) :: local_zumin, local_zumax, local_zvmin, local_zvmax, local_cvel_z, local_cvis
    real(kind=DEFAULT_PRECISION), intent(out) :: global_zumin, global_zumax, global_zvmin, global_zvmax, global_cvel_z, global_cvis

    integer :: ierr

    call mpi_allreduce(local_zumax, global_zumax, 1, PRECISION_TYPE, MPI_MAX, parallel_state%monc_communicator, ierr)
    call mpi_allreduce(local_zvmax, global_zvmax, 1, PRECISION_TYPE, MPI_MAX, parallel_state%monc_communicator, ierr)
    call mpi_allreduce(local_cvel_z, global_cvel_z, 1, PRECISION_TYPE, MPI_MAX, parallel_state%monc_communicator, ierr)
    call mpi_allreduce(local_cvis, global_cvis, 1, PRECISION_TYPE, MPI_MAX, parallel_state%monc_communicator, ierr)
    call mpi_allreduce(local_zumin, global_zumin, 1, PRECISION_TYPE, MPI_MIN, parallel_state%monc_communicator, ierr)
    call mpi_allreduce(local_zvmin, global_zvmin, 1, PRECISION_TYPE, MPI_MIN, parallel_state%monc_communicator, ierr)
  end subroutine get_global_values


  !> Handle timestep adjustment for time_basis=.true.
  !! @param current_state The current model state
  subroutine evaluate_time_basis(current_state)
    type(model_state_type), intent(inout), target :: current_state

    real(kind=DEFAULT_PRECISION) :: projected_time, dtm_trial
    integer :: sample_nts, next_sample_time, ts_to_next_cfl, next_step, interval

    next_sample_time = minval(current_state%sampling(:)%next_time)

    ! Handle timestep adjustment for time_basis=.true.
    ! Reduces timestep when approaching the next sample time.
    if (current_state%time_basis) then
      ts_to_next_cfl = current_state%cfl_frequency &
                       - mod(current_state%timestep, current_state%cfl_frequency)
      projected_time = current_state%time + current_state%dtm &
                       + (current_state%dtm_new * ts_to_next_cfl) + dtmmin
      if ( next_sample_time .gt. 0 .and. projected_time .ge. next_sample_time ) then
        sample_nts = max(1, ceiling( (next_sample_time                             &
                                      - (current_state%time + current_state%dtm) ) &
                                    / current_state%dtm_new ) )
        current_state%dtm_new = (next_sample_time - (current_state%time + current_state%dtm) ) &
                                / sample_nts
        if (l_constant_dtm) then
          current_state%dtm = options_get_real(current_state%options_database, "dtm")
          current_state%dtm_new = current_state%dtm
          sample_nts = nint((next_sample_time - (current_state%time + current_state%dtm) ) &
                           / current_state%dtm_new )
        end if

        ! Record the next sampling step for intervals matching the next time
        where(next_sample_time .eq. current_state%sampling(:)%next_time)        &
          current_state%sampling(:)%next_step = current_state%timestep + sample_nts

        current_state%update_dtm = .true.
        current_state%normal_step = .false.

      end if ! next_sample_time to be exceeded

    ! Handle force_output_on_interval by finding if time step can be reduced to hit 
    ! output interval exactly on a sampling interval.
    else if (current_state%force_output_on_interval) then

      next_step = maxval(current_state%sampling(:)%next_step, &
                  (next_sample_time .eq. current_state%sampling(:)%next_time))
      interval = maxval(current_state%sampling(:)%interval, &
                  (next_step .eq. current_state%sampling(:)%next_step) .and. &
                  (next_sample_time .eq. current_state%sampling(:)%next_time))
      sample_nts = next_step - current_state%timestep + interval

      if (sample_nts .eq. 0) sample_nts = interval
      projected_time = (current_state%time + current_state%dtm) + sample_nts*current_state%dtm_new
      dtm_trial = (next_sample_time - (current_state%time + current_state%dtm) ) &
                                / sample_nts

      if (dtm_trial .gt. current_state%dtm) then
        return
      else ! dtm should be reduced to hit output time
        current_state%dtm_new = dtm_trial
        current_state%update_dtm = .true.
        current_state%normal_step = .false.
      end if

    end if ! time_basis or force_output_on_interval

  end subroutine evaluate_time_basis


end module cfltest_mod
