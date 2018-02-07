!> Calculates the Smagorinsky eddy viscosity and diffusivity at l_w-points
module smagorinsky_mod
  use datadefn_mod, only : DEFAULT_PRECISION
  use monc_component_mod, only : component_descriptor_type
  use state_mod, only : FORWARD_STEPPING, model_state_type
  use grids_mod, only : Z_INDEX
  use optionsdatabase_mod, only : options_get_real
  use prognostics_mod, only : prognostic_field_type
  use science_constants_mod, only : smallp, rlvap_over_cp, ratio_mol_wts, G, cp
  use communication_types_mod, only : halo_communication_type, neighbour_description_type, field_data_wrapper_type
  use halo_communication_mod, only : copy_field_to_buffer, init_halo_communication, finalise_halo_communication, &
       initiate_nonblocking_halo_swap, get_single_field_per_halo_cell, copy_corner_to_buffer
  use registry_mod, only : is_component_enabled
  use logging_mod, only : LOG_ERROR, LOG_INFO, log_master_log
  use q_indices_mod, only: get_q_index, standard_q_names
  implicit none

#ifndef TEST_MODE
  private
#endif

  integer, parameter :: RICHARDSON_NUMBER_CALCULATION=2
  real(kind=DEFAULT_PRECISION) :: eps, repsh, thcona, thconb, thconap1, suba, subb, subc, subg, subh, subr, pr_n, ric, ricinv

  public smagorinsky_get_descriptor, calculate_half_squared_strain_rate, &
         calculate_richardson_number, calculate_thermal_dissipation_rate

  integer :: iqv, iql ! index for water vapour and liquid

contains

  !> Provides the descriptor back to the caller and is used in component registration
  !! @returns The termination check component descriptor
  type(component_descriptor_type) function smagorinsky_get_descriptor()
    smagorinsky_get_descriptor%name="smagorinsky"
    smagorinsky_get_descriptor%version=0.1
    smagorinsky_get_descriptor%initialisation=>initialisation_callback
    smagorinsky_get_descriptor%timestep=>timestep_callback
    smagorinsky_get_descriptor%finalisation=>finalisation_callback
  end function smagorinsky_get_descriptor

  !> Initialisation call back which will read in the coriolis configuration and set up the geostrophic winds
  !! @param current_state The current model state
  subroutine initialisation_callback(current_state)
    type(model_state_type), target, intent(inout) :: current_state


    if (.not. is_component_enabled(current_state%options_database, "diffusion")) then
      call log_master_log(LOG_ERROR, "Smagorinsky requires the diffusion component to be enabled")
    end if    
    if (.not. is_component_enabled(current_state%options_database, "viscosity")) then
      call log_master_log(LOG_ERROR, "Smagorinsky requires the viscosity component to be enabled")
    end if
    if (.not. current_state%use_viscosity_and_diffusion) then 
       call log_master_log(LOG_ERROR, "Smagorinsky requires use_viscosity_and_diffusion=.true. or monc will fail")
    endif

    if (.not. is_component_enabled(current_state%options_database, "lower_bc")) then
       call log_master_log(LOG_INFO, "LOWERBC is disabled, zero diff and vis_coeff on level 1")
       current_state%vis_coefficient%data(1,:,:)=0.0_DEFAULT_PRECISION
       current_state%diff_coefficient%data(1,:,:)=0.0_DEFAULT_PRECISION
    endif
       
    eps=0.01_DEFAULT_PRECISION
    repsh=0.5_DEFAULT_PRECISION/eps
    thcona=ratio_mol_wts*current_state%thref0
    thconb=thcona-current_state%thref0
    thconap1=thcona

    subb=options_get_real(current_state%options_database, "smag-subb")
    subc=options_get_real(current_state%options_database, "smag-subc")

    subg=1.2_DEFAULT_PRECISION
    subh=0.0_DEFAULT_PRECISION
    subr=4.0_DEFAULT_PRECISION
    pr_n=0.7_DEFAULT_PRECISION
    suba=1.0_DEFAULT_PRECISION/pr_n
    ric=0.25_DEFAULT_PRECISION
    ricinv=1.0_DEFAULT_PRECISION/ric
    if (current_state%use_viscosity_and_diffusion) then
      call init_halo_communication(current_state, get_single_field_per_halo_cell, &
           current_state%viscosity_halo_swap_state, 2, .true.)
      call init_halo_communication(current_state, get_single_field_per_halo_cell, &
           current_state%diffusion_halo_swap_state, 2, .true.)
    end if    

    if (.not. current_state%passive_q) then 
       iqv = get_q_index(standard_q_names%VAPOUR, 'smagorinsky')
       current_state%water_vapour_mixing_ratio_index=iqv
       
       iql = get_q_index(standard_q_names%CLOUD_LIQUID_MASS, 'smagorinsky')
       current_state%liquid_water_mixing_ratio_index=iql
    endif
       
  end subroutine initialisation_callback
  
  !> For each none halo cell this will calculate the subgrid terms for viscosity and diffusion
  !! @param current_state The current model state
  subroutine timestep_callback(current_state)
    type(model_state_type), target, intent(inout) :: current_state

    real(kind=DEFAULT_PRECISION), dimension(current_state%local_grid%size(Z_INDEX)) :: richardson_number, ssq

    if (.not. current_state%halo_column) then
      if (.not. current_state%use_viscosity_and_diffusion) then
        current_state%vis_coefficient%data(2:,:,:)=0.0_DEFAULT_PRECISION
        current_state%diff_coefficient%data(2:,:,:)=0.0_DEFAULT_PRECISION
        current_state%cvis=0.0_DEFAULT_PRECISION
      else
        if (current_state%field_stepping == FORWARD_STEPPING) then
          ssq=calculate_half_squared_strain_rate(current_state, current_state%u, current_state%v, current_state%w)
          richardson_number=calculate_richardson_number(current_state, ssq, current_state%th, current_state%q)
        else
          ssq=calculate_half_squared_strain_rate(current_state, current_state%zu, current_state%zv, current_state%zw)        
          richardson_number=calculate_richardson_number(current_state, ssq, current_state%zth, current_state%zq)
        end if
        call setfri(current_state, richardson_number, ssq)
        if (is_component_enabled(current_state%options_database, "cfltest")) then
           call update_viscous_number(current_state)
        endif
      end if
    end if
    if (current_state%last_timestep_column) then
      ! need to ensure not already in progress
      call initiate_nonblocking_halo_swap(current_state, current_state%diffusion_halo_swap_state, &
           copy_diff_to_halo_buffer, copy_diff_corners_to_halo_buffer)
      call initiate_nonblocking_halo_swap(current_state, current_state%viscosity_halo_swap_state, &
           copy_vis_to_halo_buffer, copy_vis_corners_to_halo_buffer)
    end if
  end subroutine timestep_callback

  !> Called when the model is finishing up, will finalise the halo communications represented by the state
  !! @param current_state The current model state
  subroutine finalisation_callback(current_state)
    type(model_state_type), target, intent(inout) :: current_state

    call finalise_halo_communication(current_state%viscosity_halo_swap_state)
    call finalise_halo_communication(current_state%diffusion_halo_swap_state)
  end subroutine finalisation_callback  

  !> Update viscous number based upon the viscosity and diffusivity coefficients
  !! @param current_state The current model state
  subroutine update_viscous_number(current_state)
    type(model_state_type), target, intent(inout) :: current_state

    integer :: k

    if (mod(current_state%timestep, current_state%cfl_frequency) == 1 .or. &
         current_state%timestep-current_state%start_timestep .le. current_state%cfl_frequency) then
      do k=2, current_state%local_grid%size(Z_INDEX)-1
        current_state%cvis=max(current_state%cvis, max(current_state%vis_coefficient%data(k, current_state%column_local_y, &
             current_state%column_local_x),current_state%diff_coefficient%data(k, current_state%column_local_y, &
             current_state%column_local_x))*(current_state%global_grid%configuration%vertical%rdzn(k+1)**2+&
             current_state%global_grid%configuration%horizontal%cx2+current_state%global_grid%configuration%horizontal%cy2))
      end do
    end if
  end subroutine update_viscous_number  

  !> Calculates the eddy viscosity (VIS) and diffusivity (DIFF) depending on the Richardson Number (RI) and half
  !! squared strain rate
  !! @param current_state The current model state
  !! @param richardson_number Richardson number
  !! @param ssq The half squared strain rate
  subroutine setfri(current_state, richardson_number, ssq)
    type(model_state_type), target, intent(inout) :: current_state
    real(kind=DEFAULT_PRECISION), dimension(:), intent(in) :: richardson_number, ssq

    integer :: k
    real(kind=DEFAULT_PRECISION) :: rifac, sctmp

    do k=2,current_state%local_grid%size(Z_INDEX)-1
      if (richardson_number(k) .le. 0.0_DEFAULT_PRECISION) then
        current_state%vis_coefficient%data(k, current_state%column_local_y, current_state%column_local_x)=&
             sqrt(1.0_DEFAULT_PRECISION-subc*richardson_number(k))
        current_state%diff_coefficient%data(k, current_state%column_local_y, current_state%column_local_x)=&
             suba*sqrt(1.-subb*richardson_number(k))        
      else if ((richardson_number(k) .gt. 0.0_DEFAULT_PRECISION) .and. (richardson_number(k) .lt. ric)) then
        rifac=(1.-richardson_number(k)*ricinv)**4
        current_state%vis_coefficient%data(k, current_state%column_local_y, current_state%column_local_x)=&
             rifac*(1.-subh*richardson_number(k))
        current_state%diff_coefficient%data(k, current_state%column_local_y, current_state%column_local_x)=&
             rifac*suba*(1.-subg*richardson_number(k))       
      else
        current_state%vis_coefficient%data(k, current_state%column_local_y, current_state%column_local_x)=0.0_DEFAULT_PRECISION
        current_state%diff_coefficient%data(k, current_state%column_local_y, current_state%column_local_x)=0.0_DEFAULT_PRECISION
      end if
      sctmp=current_state%global_grid%configuration%vertical%rneutml_sq(k)*sqrt(ssq(k))
      current_state%vis_coefficient%data(k, current_state%column_local_y, current_state%column_local_x)=&
           current_state%vis_coefficient%data(k, current_state%column_local_y, current_state%column_local_x)*sctmp
      current_state%diff_coefficient%data(k, current_state%column_local_y, current_state%column_local_x)=&
           current_state%diff_coefficient%data(k, current_state%column_local_y, current_state%column_local_x)*sctmp     
    end do

    current_state%vis_coefficient%data(current_state%local_grid%size(Z_INDEX), &
         current_state%column_local_y, current_state%column_local_x)= 0.0_DEFAULT_PRECISION
    current_state%diff_coefficient%data(current_state%local_grid%size(Z_INDEX), &
         current_state%column_local_y, current_state%column_local_x)= 0.0_DEFAULT_PRECISION
  end subroutine setfri  

  !> Calculates the richardson number depending upon the setup of the model and the method selected
  !! @param current_state The current model state
  !! @param ssq The half squared strain rate
  !! @returns The Richardson number for a specific column
  function calculate_richardson_number(current_state, ssq, th, q)
    type(model_state_type), target, intent(inout) :: current_state
    real(kind=DEFAULT_PRECISION), dimension(:), intent(in) :: ssq
    type(prognostic_field_type), intent(inout) :: th
    type(prognostic_field_type), dimension(:), intent(inout) :: q
    real(kind=DEFAULT_PRECISION), dimension(current_state%local_grid%size(Z_INDEX)) :: calculate_richardson_number

    integer :: k

    if (.not. current_state%passive_th) then
      if (.not. current_state%passive_q) then
        if (RICHARDSON_NUMBER_CALCULATION .eq. 2) then
          calculate_richardson_number=moist_ri_2(current_state, ssq, th, q)
        else
          calculate_richardson_number=moist_ri_1(current_state, ssq, th, q)
        end if
      else
        do k=2, current_state%local_grid%size(Z_INDEX)-1
          calculate_richardson_number(k) = (current_state%global_grid%configuration%vertical%dthref(k) + &
               current_state%th%data(k+1, current_state%column_local_y, current_state%column_local_x) -&
               current_state%th%data(k, current_state%column_local_y, current_state%column_local_x))* &
               current_state%global_grid%configuration%vertical%rdzn(k+1)*&
               current_state%global_grid%configuration%vertical%buoy_co(k) / ssq(k) 
        end do
      end if
    else
      calculate_richardson_number=0.0_DEFAULT_PRECISION        
    end if
  end function calculate_richardson_number

  !> Calculates another "moist version" of the Richardson number based on "change in subgrid buoyancy flux 
  !! when a fraction EPS is exchanged"
  !! @param current_state The current model state
  !! @param ssq The half squared strain rate
  !! @returns The Richardson number calculated here for a specific column
  function moist_ri_2(current_state, ssq, th, q)
    type(model_state_type), target, intent(inout) :: current_state
    real(kind=DEFAULT_PRECISION), dimension(:), intent(in) :: ssq
    type(prognostic_field_type), intent(inout) :: th
    type(prognostic_field_type), dimension(:), intent(inout) :: q
    real(kind=DEFAULT_PRECISION), dimension(current_state%local_grid%size(Z_INDEX)) :: moist_ri_2

    integer :: k
    real(kind=DEFAULT_PRECISION) :: tmpinv, scltmp, qt_k, qt_kp1, thlpr_k, thlpr_kp1, qlli, qlui, &
         thvli, thvui, thlpr_un, thlpr_ln, qtun, qtln, thvln, thvun, qlln, qlun

    do K=2,current_state%local_grid%size(Z_INDEX)-1
      tmpinv = 1.0_DEFAULT_PRECISION/ssq(k)
      scltmp=current_state%global_grid%configuration%vertical%rdzn(k+1)*&
           current_state%global_grid%configuration%vertical%buoy_co(k)
      if (current_state%use_anelastic_equations) then
        thcona=ratio_mol_wts*current_state%global_grid%configuration%vertical%thref(k)
        thconb=0.5_DEFAULT_PRECISION*(ratio_mol_wts-1.0_DEFAULT_PRECISION)*&
             (current_state%global_grid%configuration%vertical%thref(k)+&
             current_state%global_grid%configuration%vertical%thref(k+1))
        thconap1=ratio_mol_wts*current_state%global_grid%configuration%vertical%thref(k+1)
      end if    
      
      qt_k = q(current_state%water_vapour_mixing_ratio_index)%data(&
           k, current_state%column_local_y, current_state%column_local_x)+&
           q(current_state%liquid_water_mixing_ratio_index)%data(k, current_state%column_local_y, current_state%column_local_x)
      qt_kp1 = q(current_state%water_vapour_mixing_ratio_index)%data(&
           k+1, current_state%column_local_y, current_state%column_local_x)+&
           q(current_state%liquid_water_mixing_ratio_index)%data(k+1, current_state%column_local_y, current_state%column_local_x)
      thlpr_k = th%data(k, current_state%column_local_y, current_state%column_local_x) -&
           rlvap_over_cp*q(current_state%liquid_water_mixing_ratio_index)%data(&
           k, current_state%column_local_y, current_state%column_local_x)*&
           current_state%global_grid%configuration%vertical%prefrcp(k)
      ! calculate theta_l at constant (level K) pressure
      thlpr_kp1 = th%data(k+1, current_state%column_local_y, current_state%column_local_x) - &
           rlvap_over_cp*q(current_state%liquid_water_mixing_ratio_index)%data(&
           k+1, current_state%column_local_y, current_state%column_local_x)*&
           current_state%global_grid%configuration%vertical%prefrcp(k)

      !.......CALCULATE QL AND BUOYANCY OF LEVELS K AND K+1 (LOWER AND UPPER)
      !.......QL CALC AT CONSTANT PRESSURE
      qlli = max(0.0_DEFAULT_PRECISION, (qt_k-(current_state%global_grid%configuration%vertical%qsat(k) +&
           current_state%global_grid%configuration%vertical%dqsatdt(k)*&
           (current_state%global_grid%configuration%vertical%rprefrcp(k)*thlpr_k-&
           current_state%global_grid%configuration%vertical%tstarpr(k))))*&
           current_state%global_grid%configuration%vertical%qsatfac(k))
      qlui = max(0.0_DEFAULT_PRECISION, (qt_kp1-(current_state%global_grid%configuration%vertical%qsat(k) +&
           current_state%global_grid%configuration%vertical%dqsatdt(k)*&
           (current_state%global_grid%configuration%vertical%rprefrcp(k)*(thlpr_kp1+&
           current_state%global_grid%configuration%vertical%thref(k+1))-&
           current_state%global_grid%configuration%vertical%tstarpr(k)-&
           current_state%global_grid%configuration%vertical%tref(k))))*&
           current_state%global_grid%configuration%vertical%qsatfac(k))
      thvli = thlpr_k+current_state%global_grid%configuration%vertical%thref(k) +&
           ( rlvap_over_cp*current_state%global_grid%configuration%vertical%prefrcp(k)-thcona ) * qlli                         
      thvui = thlpr_kp1+current_state%global_grid%configuration%vertical%thref(k+1) +&
           ( rlvap_over_cp*current_state%global_grid%configuration%vertical%prefrcp(k)-thconap1 ) * qlui                       

      !.......MIX CONSERVATIVE QUANTITIES AT CONSTANT PRESSURE
      !.......NOTE THAT FOR COMPUTATIONAL CONVENIENCE THE RELEVANT LEVEL'S
      !.......THREF HAS BEEN OMITTED FROM ALL THV'S, THLPR_LN AND THLPR_UN AS
      !.......ONLY THE DIFFERENCE BETWEEN THV'S AT THE SAME LEVEL IS REQUIRED.
      thlpr_un = thlpr_kp1 - eps*(thlpr_kp1 - thlpr_k + current_state%global_grid%configuration%vertical%dthref(k))
      thlpr_ln = thlpr_k   + eps*(thlpr_kp1 - thlpr_k + current_state%global_grid%configuration%vertical%dthref(k))
      qtun     = qt_kp1    - eps*(qt_kp1 - qt_k) 
      qtln     = qt_k      + eps*(qt_kp1 - qt_k)
      !.......CALCULATE LIQUID WATER AND BUOYANCY OF MIXED AIR
      qlln = max(0.0_DEFAULT_PRECISION,( qtln - ( current_state%global_grid%configuration%vertical%qsat(k) +&
           current_state%global_grid%configuration%vertical%dqsatdt(k)*&
           (current_state%global_grid%configuration%vertical%rprefrcp(k)*&
           thlpr_ln-current_state%global_grid%configuration%vertical%tstarpr(k)))&
           )*current_state%global_grid%configuration%vertical%qsatfac(k))
      qlun = max(0.0_DEFAULT_PRECISION,( qtun - ( current_state%global_grid%configuration%vertical%qsat(k) + &
           current_state%global_grid%configuration%vertical%dqsatdt(k)*&
           (current_state%global_grid%configuration%vertical%rprefrcp(k)*&
           (thlpr_un+current_state%global_grid%configuration%vertical%thref(k+1))-&
           current_state%global_grid%configuration%vertical%tstarpr(k)-current_state%global_grid%configuration%vertical%tref(k)))&
           )*current_state%global_grid%configuration%vertical%qsatfac(k))
      thvln = thlpr_ln+current_state%global_grid%configuration%vertical%thref(k) +&
           ( rlvap_over_cp*current_state%global_grid%configuration%vertical%prefrcp(k)-thcona) * qlln
      thvun = thlpr_un+current_state%global_grid%configuration%vertical%thref(k+1) +&
           ( rlvap_over_cp*current_state%global_grid%configuration%vertical%prefrcp(k)-thconap1) * qlun

      !.....CHANGE IN BUOYANCY DUE TO FRAC MIXING IS USED TO DETERMINE l_ri
      moist_ri_2(k) = scltmp*tmpinv*(thconb*(qt_kp1-qt_k) +((thvln-thvli)-(thvun-thvui))*repsh)
    end do
  end function moist_ri_2

  !> Calculates numerator of Richarson number as the difference in buoyancy between a parcel at K+1 and a parcel lifted from
  !! K to K+1. Condensation is calculated by assuming that TL=T-(Lv/CP)*QL+gz/CP and QT = QV+QL are conserved on lifting.
  !! Microphysics conversions are assumed to happen slowly compared to turbulent processes so are, therefore, neglected
  !! in the calculation.
  !! @param current_state The current model state
  !! @param ssq The half squared strain rate
  !! @returns The Richardson number calculated here for a specific column
  function moist_ri_1(current_state, ssq, th, q)
    type(model_state_type), target, intent(inout) :: current_state
    real(kind=DEFAULT_PRECISION), dimension(:), intent(in) :: ssq
    type(prognostic_field_type), intent(inout) :: th
    type(prognostic_field_type), dimension(:), intent(inout) :: q
    real(kind=DEFAULT_PRECISION), dimension(current_state%local_grid%size(Z_INDEX)) :: moist_ri_1

    integer :: k
    real(kind=DEFAULT_PRECISION) :: liquid_water_at_kp1, temperature_at_kp1, dtlref, tmpinv, &
         liquid_water_static_energy_temperature, total_water_substance, total_water_substance_p1, scltmp, qloading

    ! Set up the liquid water static energy & total water
    do k=2, current_state%local_grid%size(Z_INDEX)
      liquid_water_static_energy_temperature = &
           th%data(k, current_state%column_local_y, current_state%column_local_x)*&
           current_state%global_grid%configuration%vertical%rprefrcp(k) - rlvap_over_cp*&
           q(current_state%liquid_water_mixing_ratio_index)%data(k, current_state%column_local_y, current_state%column_local_x)
      if (k==2) then
        total_water_substance = q(current_state%water_vapour_mixing_ratio_index)%data(&
             k, current_state%column_local_y, current_state%column_local_x) + &
             q(current_state%liquid_water_mixing_ratio_index)%data(k, current_state%column_local_y, current_state%column_local_x)
      else
        total_water_substance=total_water_substance_p1
      end if
      total_water_substance_p1 = q(current_state%water_vapour_mixing_ratio_index)%data(&
           k+1, current_state%column_local_y, current_state%column_local_x) + &
           q(current_state%liquid_water_mixing_ratio_index)%data(k+1, current_state%column_local_y, current_state%column_local_x)
      if (k .le. current_state%local_grid%size(Z_INDEX)-1) then
        tmpinv = 1.0_DEFAULT_PRECISION/ssq(k)
        dtlref = current_state%global_grid%configuration%vertical%tref(k+1)-&
             current_state%global_grid%configuration%vertical%tref(k)+(G/cp)*&
             current_state%global_grid%configuration%vertical%dzn(k+1)
        temperature_at_kp1=current_state%global_grid%configuration%vertical%tstarpr(k+1)+&
             current_state%global_grid%configuration%vertical%tref(k+1)+&
             (liquid_water_static_energy_temperature+rlvap_over_cp*(total_water_substance-&
             current_state%global_grid%configuration%vertical%qsat(k+1)) - &
             dtlref - current_state%global_grid%configuration%vertical%tstarpr(k+1) )*&
             current_state%global_grid%configuration%vertical%qsatfac(k+1)
        liquid_water_at_kp1=total_water_substance -( current_state%global_grid%configuration%vertical%qsat(k+1)+&
             current_state%global_grid%configuration%vertical%dqsatdt(k+1)*(temperature_at_kp1- &
             current_state%global_grid%configuration%vertical%tstarpr(k+1)-&
             current_state%global_grid%configuration%vertical%tref(k+1)) )
        if (liquid_water_at_kp1 .le. 0.0_DEFAULT_PRECISION) then
          liquid_water_at_kp1=0.0_DEFAULT_PRECISION
          temperature_at_kp1=current_state%global_grid%configuration%vertical%tref(k+1)+&
               (liquid_water_static_energy_temperature-dtlref)
        end if
        scltmp=current_state%global_grid%configuration%vertical%rdzn(k+1)*&
             current_state%global_grid%configuration%vertical%buoy_co(k)
        ! Calculate the Richardson number taking into account the water
        ! loading terms from the other hydrometeors in the buoyancy
        ! calculation
        qloading = current_state%cq(current_state%water_vapour_mixing_ratio_index)*&
             (total_water_substance_p1-total_water_substance)-ratio_mol_wts*&
             (q(current_state%liquid_water_mixing_ratio_index)%data(&
             k+1, current_state%column_local_y, current_state%column_local_x)-liquid_water_at_kp1)
        !if(IRAINP.gt.0) qloading = qloading + current_state%cq(IQR)*(l_q(J,K+1,IQR)-l_q(J,K,IQR))
        !if(ISNOWP.gt.0) qloading = qloading + current_state%cq(IQS)*(l_q(J,K+1,IQS)-l_q(J,K,IQS))
        !if(ICECLP.gt.0) qloading = qloading + current_state%cq(IQI)*(l_q(J,K+1,IQI)-l_q(J,K,IQI))
        !if(IGRAUP.gt.0) qloading = qloading + currnet_state%cq(IQG)*(l_q(J,K+1,IQG)-l_q(J,K,IQG))

        moist_ri_1(k) = ((th%data(k+1, current_state%column_local_y, current_state%column_local_x)+&
             current_state%global_grid%configuration%vertical%thref(k+1) - temperature_at_kp1*&
             current_state%global_grid%configuration%vertical%prefrcp(k+1))&
             +current_state%global_grid%configuration%vertical%thref(k+1)*qloading)*scltmp*tmpinv
      end if
    end do
  end function moist_ri_1  
  
  !> Calculates the half squared strain rate on l_w-points which is used in determining the Richardson number. 
  !! CX=1./DX, CY=1./DY, RDZ(K)=1./DZ(K), RDZN(K) =1./DZN(K)
  !! _SSQ= 0.5*^^DU_I/DX_J+DU_J/DX_I^^**2
  !! _SSQIJ= (DU_I/DX_J+DU_J/DX_I)**2
  !! _Hence SSQ= SUM(I,J) {0.5*(SSQIJ)}
  !! @param current_state The current model state
  !! @returns The hard squared strain rate for a specific column
  function calculate_half_squared_strain_rate(current_state, u, v, w)
    type(model_state_type), target, intent(inout) :: current_state
    type(prognostic_field_type), intent(inout) :: u, v, w
    real(kind=DEFAULT_PRECISION) :: calculate_half_squared_strain_rate(current_state%local_grid%size(Z_INDEX))

    integer :: k
    real(kind=DEFAULT_PRECISION) :: ssq11, ssq22, ssq33, ssq13, ssq23, ssq12

    do k=2,current_state%local_grid%size(Z_INDEX)-1   
#ifdef U_ACTIVE    
      ssq11=current_state%global_grid%configuration%horizontal%cx2*(&
           (u%data(k+1, current_state%column_local_y, current_state%column_local_x)-&
           u%data(k+1, current_state%column_local_y, current_state%column_local_x-1))**2+&
           (u%data(k, current_state%column_local_y, current_state%column_local_x)-&
           u%data(k, current_state%column_local_y, current_state%column_local_x-1))**2)
#else
      ssq11=0.0_DEFAULT_PRECISION
#endif
#ifdef V_ACTIVE
      ssq22=current_state%global_grid%configuration%horizontal%cy2*(&
           (v%data(k+1, current_state%column_local_y, current_state%column_local_x)-&
           v%data(k+1, current_state%column_local_y-1, current_state%column_local_x))**2+&
           (v%data(k, current_state%column_local_y, current_state%column_local_x)-&
           v%data(k, current_state%column_local_y-1, current_state%column_local_x))**2)
#else
      ssq22=0.0_DEFAULT_PRECISION
#endif
#ifdef W_ACTIVE
      ssq33=((w%data(k, current_state%column_local_y, current_state%column_local_x)-&
           w%data(k-1, current_state%column_local_y, current_state%column_local_x))*&
           current_state%global_grid%configuration%vertical%rdz(k))**2 +&
           ((w%data(k+1, current_state%column_local_y, current_state%column_local_x)-&
           w%data(k, current_state%column_local_y, current_state%column_local_x))*&
           current_state%global_grid%configuration%vertical%rdz(k+1))**2
#else
      ssq33=0.0_DEFAULT_PRECISION
#endif
#if defined(U_ACTIVE) && defined(W_ACTIVE)
      ! Average over 2 points U and W
      ssq13=(((u%data(k+1, current_state%column_local_y, current_state%column_local_x)-&
           u%data(k, current_state%column_local_y, current_state%column_local_x))*&
           current_state%global_grid%configuration%vertical%rdzn(k+1)+&
           (w%data(k, current_state%column_local_y, current_state%column_local_x+1)-&
           w%data(k, current_state%column_local_y, current_state%column_local_x))*&
           current_state%global_grid%configuration%horizontal%cx)**2+&
           ((u%data(k+1, current_state%column_local_y, current_state%column_local_x-1)-&
           u%data(k, current_state%column_local_y, current_state%column_local_x-1))*&
           current_state%global_grid%configuration%vertical%rdzn(k+1)+&
           (w%data(k, current_state%column_local_y, current_state%column_local_x)-&
           w%data(k, current_state%column_local_y, current_state%column_local_x-1))*&
           current_state%global_grid%configuration%horizontal%cx)**2)*0.5_DEFAULT_PRECISION
#elif defined(U_ACTIVE) && !defined(W_ACTIVE)
      ssq13=(((u%data(k+1, current_state%column_local_y, current_state%column_local_x)-&
           u%data(k, current_state%column_local_y, current_state%column_local_x))*&
           current_state%global_grid%configuration%vertical%rdzn(k+1))**2+&
           ((u%data(k+1, current_state%column_local_y, current_state%column_local_x-1)-&
           u%data(k, current_state%column_local_y, current_state%column_local_x-1))*&
           current_state%global_grid%configuration%vertical%rdzn(k+1))**2)*0.5_DEFAULT_PRECISION
#elif !defined(U_ACTIVE) && defined(W_ACTIVE)
      ssq13=(((w%data(k, current_state%column_local_y, current_state%column_local_x+1)-&
           w%data(k, current_state%column_local_y, current_state%column_local_x))*&
           current_state%global_grid%configuration%horizontal%cx)**2+&
           ((w%data(k, current_state%column_local_y, current_state%column_local_x)-&
           w%data(k, current_state%column_local_y, current_state%column_local_x-1))*&
           current_state%global_grid%configuration%horizontal%cx)**2)*0.5_DEFAULT_PRECISION
#else
      ssq13=0.0_DEFAULT_PRECISION
#endif
#if defined(W_ACTIVE) && defined(V_ACTIVE)
      ! Average over 2 points W and V
      ssq23=(((w%data(k, current_state%column_local_y, current_state%column_local_x)-&
           w%data(k, current_state%column_local_y-1, current_state%column_local_x))*&
           current_state%global_grid%configuration%horizontal%cy+&
           (v%data(k+1, current_state%column_local_y-1, current_state%column_local_x)-&
           v%data(k, current_state%column_local_y-1, current_state%column_local_x))*&
           current_state%global_grid%configuration%vertical%rdzn(k+1))**2+&
           ((w%data(k, current_state%column_local_y+1, current_state%column_local_x)-&
           w%data(k, current_state%column_local_y, current_state%column_local_x))*&
           current_state%global_grid%configuration%horizontal%cy+&
           (v%data(k+1, current_state%column_local_y, current_state%column_local_x)-&
           v%data(k, current_state%column_local_y, current_state%column_local_x))*&
           current_state%global_grid%configuration%vertical%rdzn(k+1))**2)*0.5_DEFAULT_PRECISION
#elif defined(W_ACTIVE) && !defined(V_ACTIVE)
      ssq23=(((w%data(k, current_state%column_local_y, current_state%column_local_x)-&
           w%data(k, current_state%column_local_y-1, current_state%column_local_x))*&
           current_state%global_grid%configuration%horizontal%cy)**2+&
           ((w%data(k, current_state%column_local_y+1, current_state%column_local_x)-&
           w%data(k, current_state%column_local_y, current_state%column_local_x))*&
           current_state%global_grid%configuration%horizontal%cy)**2)*0.5_DEFAULT_PRECISION
#elif !defined(W_ACTIVE) && defined(V_ACTIVE)
      ssq23=(((v%data(k+1, current_state%column_local_y-1, current_state%column_local_x)-&
           v%data(k, current_state%column_local_y-1, current_state%column_local_x))*&
           current_state%global_grid%configuration%vertical%rdzn(k+1))**2+&
           ((v%data(k+1, current_state%column_local_y, current_state%column_local_x)-&
           v%data(k, current_state%column_local_y, current_state%column_local_x))*&
           current_state%global_grid%configuration%vertical%rdzn(k+1))**2)*0.5_DEFAULT_PRECISION
#else
      ssq23=0.0_DEFAULT_PRECISION
#endif

#if defined(U_ACTIVE) && defined(V_ACTIVE)
      ! Average over 8 points from U and V
      ssq12=(((((u%data(k, current_state%column_local_y, current_state%column_local_x-1)-&
           u%data(k, current_state%column_local_y-1, current_state%column_local_x-1))*&
           current_state%global_grid%configuration%horizontal%cy+&
           (v%data(k, current_state%column_local_y-1, current_state%column_local_x)-&
           v%data(k, current_state%column_local_y-1, current_state%column_local_x-1))*&
           current_state%global_grid%configuration%horizontal%cx)**2 +&
           ((u%data(k, current_state%column_local_y+1, current_state%column_local_x-1)-&
           u%data(k, current_state%column_local_y, current_state%column_local_x-1))*&
           current_state%global_grid%configuration%horizontal%cy+&
           (v%data(k, current_state%column_local_y, current_state%column_local_x)-&
           v%data(k, current_state%column_local_y, current_state%column_local_x-1))*&
           current_state%global_grid%configuration%horizontal%cx)**2) +(&
           ((u%data(k, current_state%column_local_y, current_state%column_local_x)-&
           u%data(k, current_state%column_local_y-1, current_state%column_local_x))*&
           current_state%global_grid%configuration%horizontal%cy+&
           (v%data(k, current_state%column_local_y-1, current_state%column_local_x+1)-&
           v%data(k, current_state%column_local_y-1, current_state%column_local_x))*&
           current_state%global_grid%configuration%horizontal%cx)**2 +&
           ((u%data(k, current_state%column_local_y+1, current_state%column_local_x)-&
           u%data(k, current_state%column_local_y, current_state%column_local_x))*&
           current_state%global_grid%configuration%horizontal%cy+&
           (v%data(k, current_state%column_local_y, current_state%column_local_x+1)-&
           v%data(k, current_state%column_local_y, current_state%column_local_x))*&
           current_state%global_grid%configuration%horizontal%cx)**2))+((&
           ((u%data(k+1, current_state%column_local_y, current_state%column_local_x-1)-&
           u%data(k+1, current_state%column_local_y-1, current_state%column_local_x-1))*&
           current_state%global_grid%configuration%horizontal%cy+&
           (v%data(k+1, current_state%column_local_y-1, current_state%column_local_x)-&
           v%data(k+1, current_state%column_local_y-1, current_state%column_local_x-1))*&
           current_state%global_grid%configuration%horizontal%cx)**2+&
           ((u%data(k+1, current_state%column_local_y+1, current_state%column_local_x-1)-&
           u%data(k+1, current_state%column_local_y, current_state%column_local_x-1))*&
           current_state%global_grid%configuration%horizontal%cy+&
           (v%data(k+1, current_state%column_local_y, current_state%column_local_x)-&
           v%data(k+1, current_state%column_local_y, current_state%column_local_x-1))*&
           current_state%global_grid%configuration%horizontal%cx)**2)+(&
           ((u%data(k+1, current_state%column_local_y, current_state%column_local_x)-&
           u%data(k+1, current_state%column_local_y-1, current_state%column_local_x))*&
           current_state%global_grid%configuration%horizontal%cy+&
           (v%data(k+1, current_state%column_local_y-1, current_state%column_local_x+1)-&
           v%data(k+1, current_state%column_local_y-1, current_state%column_local_x))*&
           current_state%global_grid%configuration%horizontal%cx)**2+&
           ((u%data(k+1, current_state%column_local_y+1, current_state%column_local_x)-&
           u%data(k+1, current_state%column_local_y, current_state%column_local_x))*&
           current_state%global_grid%configuration%horizontal%cy+&
           (v%data(k+1, current_state%column_local_y, current_state%column_local_x+1)-&
           v%data(k+1, current_state%column_local_y, current_state%column_local_x))*&
           current_state%global_grid%configuration%horizontal%cx)**2)))*0.125_DEFAULT_PRECISION

#elif defined(U_ACTIVE) && !defined(V_ACTIVE)

      ssq12=(((((u%data(k, current_state%column_local_y, current_state%column_local_x-1)-&
           u%data(k, current_state%column_local_y-1, current_state%column_local_x-1))*&
           current_state%global_grid%configuration%horizontal%cy)**2 +&
           ((u%data(k, current_state%column_local_y+1, current_state%column_local_x-1)-&
           u%data(k, current_state%column_local_y, current_state%column_local_x-1))*&
           current_state%global_grid%configuration%horizontal%cy)**2) +(&
           ((u%data(k, current_state%column_local_y, current_state%column_local_x)-&
           u%data(k, current_state%column_local_y-1, current_state%column_local_x))*&
           current_state%global_grid%configuration%horizontal%cy)**2 +&
           ((u%data(k, current_state%column_local_y+1, current_state%column_local_x)-&
           u%data(k, current_state%column_local_y, current_state%column_local_x))*&
           current_state%global_grid%configuration%horizontal%cy)**2))+((&
           ((u%data(k+1, current_state%column_local_y, current_state%column_local_x-1)-&
           u%data(k+1, current_state%column_local_y-1, current_state%column_local_x-1))*&
           current_state%global_grid%configuration%horizontal%cy)**2+&
           ((u%data(k+1, current_state%column_local_y+1, current_state%column_local_x-1)-&
           u%data(k+1, current_state%column_local_y, current_state%column_local_x-1))*&
           current_state%global_grid%configuration%horizontal%cy)**2)+(&
           ((u%data(k+1, current_state%column_local_y, current_state%column_local_x)-&
           u%data(k+1, current_state%column_local_y-1, current_state%column_local_x))*&
           current_state%global_grid%configuration%horizontal%cy)**2+&
           ((u%data(k+1, current_state%column_local_y+1, current_state%column_local_x)-&
           u%data(k+1, current_state%column_local_y, current_state%column_local_x))*&
           current_state%global_grid%configuration%horizontal%cy)**2)))*0.125_DEFAULT_PRECISION

#elif !defined(U_ACTIVE) && defined(V_ACTIVE)

      ssq12=(((((v%data(k, current_state%column_local_y-1, current_state%column_local_x)-&
           v%data(k, current_state%column_local_y-1, current_state%column_local_x-1))*&
           current_state%global_grid%configuration%horizontal%cx)**2 +&
           ((v%data(k, current_state%column_local_y, current_state%column_local_x)-&
           v%data(k, current_state%column_local_y, current_state%column_local_x-1))*&
           current_state%global_grid%configuration%horizontal%cx)**2) +&
           (((v%data(k, current_state%column_local_y-1, current_state%column_local_x+1)-&
           v%data(k, current_state%column_local_y-1, current_state%column_local_x))*&
           current_state%global_grid%configuration%horizontal%cx)**2 +&
           ((v%data(k, current_state%column_local_y, current_state%column_local_x+1)-&
           v%data(k, current_state%column_local_y, current_state%column_local_x))*&
           current_state%global_grid%configuration%horizontal%cx)**2))+((&
           ((v%data(k+1, current_state%column_local_y-1, current_state%column_local_x)-&
           v%data(k+1, current_state%column_local_y-1, current_state%column_local_x-1))*&
           current_state%global_grid%configuration%horizontal%cx)**2+&
           ((v%data(k+1, current_state%column_local_y, current_state%column_local_x)-&
           v%data(k+1, current_state%column_local_y, current_state%column_local_x-1))*&
           current_state%global_grid%configuration%horizontal%cx)**2)+(&
           ((v%data(k+1, current_state%column_local_y-1, current_state%column_local_x+1)-&
           v%data(k+1, current_state%column_local_y-1, current_state%column_local_x))*&
           current_state%global_grid%configuration%horizontal%cx)**2+&
           ((v%data(k+1, current_state%column_local_y, current_state%column_local_x+1)-&
           v%data(k+1, current_state%column_local_y, current_state%column_local_x))*&
           current_state%global_grid%configuration%horizontal%cx)**2)))*0.125_DEFAULT_PRECISION
#else
      ssq12=0.0_DEFAULT_PRECISION
#endif
      calculate_half_squared_strain_rate(k)=ssq11+ssq22+ssq33+ssq13+ssq23+ssq12+smallp
    end do
  end function calculate_half_squared_strain_rate

!===========================================
  !! @param current_state The current model state
  !! @returns The thermal dissiptation rate for a specific column: backscatter is NOT included -- exercise for partners
  function calculate_thermal_dissipation_rate(current_state, th)
    type(model_state_type), target, intent(inout) :: current_state
    type(prognostic_field_type), intent(inout) :: th
    real(kind=DEFAULT_PRECISION) :: calculate_thermal_dissipation_rate(current_state%local_grid%size(Z_INDEX))

    integer :: k

    do k=2,current_state%local_grid%size(Z_INDEX)-1
      calculate_thermal_dissipation_rate(k) = &
        current_state%diff_coefficient%data(k, current_state%column_local_y, current_state%column_local_x) * ( &
          (current_state%global_grid%configuration%vertical%rdzn(k+1) * &
            (current_state%th%data(k+1, current_state%column_local_y, current_state%column_local_x) - &
             current_state%th%data(k, current_state%column_local_y, current_state%column_local_x) - &
             current_state%global_grid%configuration%vertical%dthref(k) ) ) ** 2 + &
          0.25 * current_state%global_grid%configuration%horizontal%cx2 * &
            ( (current_state%th%data(k, current_state%column_local_y, current_state%column_local_x+1) - &
               current_state%th%data(k, current_state%column_local_y, current_state%column_local_x) ) ** 2 + &
              (current_state%th%data(k, current_state%column_local_y, current_state%column_local_x) - &
               current_state%th%data(k, current_state%column_local_y, current_state%column_local_x-1) ) ** 2 + &
              (current_state%th%data(k+1, current_state%column_local_y, current_state%column_local_x+1) - &
               current_state%th%data(k+1, current_state%column_local_y, current_state%column_local_x) ) **2 + &
              (current_state%th%data(k+1, current_state%column_local_y, current_state%column_local_x) - &
               current_state%th%data(k+1, current_state%column_local_y, current_state%column_local_x-1) ) ** 2 ) + &
          0.25 * current_state%global_grid%configuration%horizontal%cy2 * &
            ( (current_state%th%data(k, current_state%column_local_y+1, current_state%column_local_x) - &
               current_state%th%data(k, current_state%column_local_y, current_state%column_local_x) ) ** 2 + &
              (current_state%th%data(k, current_state%column_local_y, current_state%column_local_x) - &
               current_state%th%data(k, current_state%column_local_y-1, current_state%column_local_x) ) ** 2 + &
              (current_state%th%data(k+1, current_state%column_local_y+1, current_state%column_local_x) - &
               current_state%th%data(k+1, current_state%column_local_y, current_state%column_local_x) ) **2 + &
              (current_state%th%data(k+1, current_state%column_local_y, current_state%column_local_x) - &
               current_state%th%data(k+1, current_state%column_local_y-1, current_state%column_local_x) ) ** 2 ) )
    end do

  end function calculate_thermal_dissipation_rate

!===========================================

  !> Copies the viscosity field data to halo buffers for a specific process in a dimension and halo cell
  !! @param current_state The current model state
  !! @param neighbour_descriptions Description of the neighbour halo swapping status
  !! @param dim Dimension to copy from
  !! @param source_index The source index of the dimension we are reading from in the prognostic field
  !! @param pid_location Location of the neighbouring process in the local stored data structures
  !! @param current_page The current (next) buffer page to copy into
  !! @param source_data Optional source data which is read from
  subroutine copy_vis_to_halo_buffer(current_state, neighbour_description, dim, source_index, &
       pid_location, current_page, source_data)
    type(model_state_type), intent(inout) :: current_state
    integer, intent(in) :: dim, pid_location, source_index
    integer, intent(inout) :: current_page(:)
    type(neighbour_description_type), intent(inout) :: neighbour_description
    type(field_data_wrapper_type), dimension(:), intent(in), optional :: source_data

    call copy_field_to_buffer(current_state%local_grid, neighbour_description%send_halo_buffer, &
         current_state%vis_coefficient%data, dim, source_index, current_page(pid_location))

    current_page(pid_location)=current_page(pid_location)+1
  end subroutine copy_vis_to_halo_buffer

  !> Copies the viscosity field corner data to halo buffers for a specific process
  !! @param current_state The current model state
  !! @param neighbour_descriptions Description of the neighbour halo swapping status
  !! @param corner_loc The corner location
  !! @param x_source_index The X source index of the dimension we are reading from in the prognostic field
  !! @param y_source_index The Y source index of the dimension we are reading from in the prognostic field
  !! @param pid_location Location of the neighbouring process in the local stored data structures
  !! @param current_page The current (next) buffer page to copy into
  !! @param source_data Optional source data which is read from
  subroutine copy_vis_corners_to_halo_buffer(current_state, neighbour_description, corner_loc, x_source_index, y_source_index, &
       pid_location, current_page, source_data)
    type(model_state_type), intent(inout) :: current_state
    integer, intent(in) :: corner_loc, pid_location, x_source_index, y_source_index
    integer, intent(inout) :: current_page(:)
    type(neighbour_description_type), intent(inout) :: neighbour_description
    type(field_data_wrapper_type), dimension(:), intent(in), optional :: source_data

    call copy_corner_to_buffer(current_state%local_grid, neighbour_description%send_corner_buffer, &
         current_state%vis_coefficient%data, corner_loc, x_source_index, y_source_index, current_page(pid_location))

    current_page(pid_location)=current_page(pid_location)+1
  end subroutine copy_vis_corners_to_halo_buffer 

  !> Copies the diffusion field data to halo buffers for a specific process in a dimension and halo cell
  !! @param current_state The current model state
  !! @param neighbour_descriptions Description of the neighbour halo swapping status
  !! @param dim Dimension to copy from
  !! @param source_index The source index of the dimension we are reading from in the prognostic field
  !! @param pid_location Location of the neighbouring process in the local stored data structures
  !! @param current_page The current (next) buffer page to copy into
  !! @param source_data Optional source data which is read from
  subroutine copy_diff_to_halo_buffer(current_state, neighbour_description, dim, source_index, &
       pid_location, current_page, source_data)
    type(model_state_type), intent(inout) :: current_state
    integer, intent(in) :: dim, pid_location, source_index
    integer, intent(inout) :: current_page(:)
    type(neighbour_description_type), intent(inout) :: neighbour_description
    type(field_data_wrapper_type), dimension(:), intent(in), optional :: source_data

    call copy_field_to_buffer(current_state%local_grid, neighbour_description%send_halo_buffer, &
         current_state%diff_coefficient%data, dim, source_index, current_page(pid_location))

    current_page(pid_location)=current_page(pid_location)+1
  end subroutine copy_diff_to_halo_buffer

  !> Copies the diffusion field corner data to halo buffers for a specific process
  !! @param current_state The current model state
  !! @param neighbour_descriptions Description of the neighbour halo swapping status
  !! @param corner_loc The corner location
  !! @param x_source_index The X source index of the dimension we are reading from in the prognostic field
  !! @param y_source_index The Y source index of the dimension we are reading from in the prognostic field
  !! @param pid_location Location of the neighbouring process in the local stored data structures
  !! @param current_page The current (next) buffer page to copy into
  !! @param source_data Optional source data which is read from
  subroutine copy_diff_corners_to_halo_buffer(current_state, neighbour_description, corner_loc, x_source_index, y_source_index, &
       pid_location, current_page, source_data)
    type(model_state_type), intent(inout) :: current_state
    integer, intent(in) :: corner_loc, pid_location, x_source_index, y_source_index
    integer, intent(inout) :: current_page(:)
    type(neighbour_description_type), intent(inout) :: neighbour_description
    type(field_data_wrapper_type), dimension(:), intent(in), optional :: source_data

    call copy_corner_to_buffer(current_state%local_grid, neighbour_description%send_corner_buffer, &
         current_state%diff_coefficient%data, corner_loc, x_source_index, y_source_index, current_page(pid_location))

    current_page(pid_location)=current_page(pid_location)+1
  end subroutine copy_diff_corners_to_halo_buffer  
end module smagorinsky_mod
