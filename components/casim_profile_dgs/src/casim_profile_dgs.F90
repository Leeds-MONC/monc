module casim_profile_dgs_mod
  use monc_component_mod, only :  COMPONENT_DOUBLE_DATA_TYPE, COMPONENT_ARRAY_FIELD_TYPE, &
       component_descriptor_type, component_field_value_type, component_field_information_type
  use grids_mod, only : Z_INDEX, Y_INDEX, X_INDEX
  use registry_mod, only : is_component_enabled
  use logging_mod, only : LOG_ERROR, log_master_log
  use state_mod, only : model_state_type
  use datadefn_mod, only : DEFAULT_PRECISION

  ! needs casim modules
  use mphys_switches, only: l_warm
  ! need casdiags for the diag switches
  use generic_diagnostic_variables, ONLY: casdiags
  ! and casim component structure, which contains the process rate data
  use casim_monc_dgs_space, only: casim_monc_dgs

  implicit none

#ifndef TEST_MODE
  private
#endif

  integer :: total_points, iqv, iql
  real(kind=DEFAULT_PRECISION), dimension(:), allocatable ::     &
       ! local process rate totals
       phomc_tot, pinuc_tot, pidep_tot, psdep_tot, piacw_tot, psacw_tot, psacr_tot, pisub_tot,   &
       pssub_tot, pimlt_tot, psmlt_tot, psaut_tot, psaci_tot, praut_tot, pracw_tot, prevp_tot,   &
       pgacw_tot, pgacs_tot, pgmlt_tot, pgsub_tot, psedi_tot, pseds_tot, psedr_tot, psedg_tot,   &
       psedl_tot, pcond_tot, &
       ! local total tendencies
       dth_mphys_tot, dth_cond_evap_tot, dqv_mphys_tot, dqv_cond_evap_tot, &
       dqc_mphys_tot, dqr_mphys_tot, dqi_mphys_tot, dqs_mphys_tot, dqg_mphys_tot

  public casim_profile_dgs_get_descriptor

contains

  !> Provides the component descriptor for the core to register
  !! @returns The descriptor describing this component
  type(component_descriptor_type) function casim_profile_dgs_get_descriptor()
    casim_profile_dgs_get_descriptor%name="casim_profile_dgs"
    casim_profile_dgs_get_descriptor%version=0.1

    casim_profile_dgs_get_descriptor%initialisation=>initialisation_callback
    casim_profile_dgs_get_descriptor%timestep=>timestep_callback

    casim_profile_dgs_get_descriptor%field_value_retrieval=>field_value_retrieval_callback
    casim_profile_dgs_get_descriptor%field_information_retrieval=>field_information_retrieval_callback

    allocate(casim_profile_dgs_get_descriptor%published_fields(35))

    casim_profile_dgs_get_descriptor%published_fields(1)="phomc_total"
    casim_profile_dgs_get_descriptor%published_fields(2)="pinuc_total"
    casim_profile_dgs_get_descriptor%published_fields(3)="pidep_total"
    casim_profile_dgs_get_descriptor%published_fields(4)="psdep_total"
    casim_profile_dgs_get_descriptor%published_fields(5)="piacw_total"
    casim_profile_dgs_get_descriptor%published_fields(6)="psacw_total"
    casim_profile_dgs_get_descriptor%published_fields(7)="psacr_total"
    casim_profile_dgs_get_descriptor%published_fields(8)="pisub_total"
    casim_profile_dgs_get_descriptor%published_fields(9)="pssub_total"
    casim_profile_dgs_get_descriptor%published_fields(10)="pimlt_total"
    casim_profile_dgs_get_descriptor%published_fields(11)="psmlt_total"
    casim_profile_dgs_get_descriptor%published_fields(12)="psaut_total"
    casim_profile_dgs_get_descriptor%published_fields(13)="psaci_total"
    casim_profile_dgs_get_descriptor%published_fields(14)="praut_total"
    casim_profile_dgs_get_descriptor%published_fields(15)="pracw_total"
    casim_profile_dgs_get_descriptor%published_fields(16)="prevp_total"
    casim_profile_dgs_get_descriptor%published_fields(17)="pgacw_total"
    casim_profile_dgs_get_descriptor%published_fields(18)="pgacs_total"
    casim_profile_dgs_get_descriptor%published_fields(19)="pgmlt_total"
    casim_profile_dgs_get_descriptor%published_fields(20)="pgsub_total"
    casim_profile_dgs_get_descriptor%published_fields(21)="psedi_total"
    casim_profile_dgs_get_descriptor%published_fields(22)="pseds_total"
    casim_profile_dgs_get_descriptor%published_fields(23)="psedr_total"
    casim_profile_dgs_get_descriptor%published_fields(24)="psedg_total"
    casim_profile_dgs_get_descriptor%published_fields(25)="psedl_total"
    casim_profile_dgs_get_descriptor%published_fields(26)="pcond_total"
    casim_profile_dgs_get_descriptor%published_fields(27)="dth_mphys_total"
    casim_profile_dgs_get_descriptor%published_fields(28)="dth_cond_evap_total"
    casim_profile_dgs_get_descriptor%published_fields(29)="dqv_mphys_total"
    casim_profile_dgs_get_descriptor%published_fields(30)="dqv_cond_evap_total"
    casim_profile_dgs_get_descriptor%published_fields(31)="dqc_mphys_total"
    casim_profile_dgs_get_descriptor%published_fields(32)="dqr_mphys_total"
    casim_profile_dgs_get_descriptor%published_fields(33)="dqi_mphys_total"
    casim_profile_dgs_get_descriptor%published_fields(34)="dqs_mphys_total"
    casim_profile_dgs_get_descriptor%published_fields(35)="dqg_mphys_total"
    
  end function casim_profile_dgs_get_descriptor

  subroutine initialisation_callback(current_state)
    type(model_state_type), target, intent(inout) :: current_state

    integer :: k

    if (.not. is_component_enabled(current_state%options_database, "casim")) then
      call log_master_log(LOG_ERROR, "Casim profile diags requested but casim not enabled - check config, STOP")
    end if 

    ! allocate local arrays for the horizontal wind averages
    allocate(psedl_tot(current_state%local_grid%size(Z_INDEX)), &
         pcond_tot(current_state%local_grid%size(Z_INDEX)), &
         praut_tot(current_state%local_grid%size(Z_INDEX)), &
         pracw_tot(current_state%local_grid%size(Z_INDEX)), & 
         prevp_tot(current_state%local_grid%size(Z_INDEX)), &
         psedr_tot(current_state%local_grid%size(Z_INDEX)), &
         dth_mphys_tot(current_state%local_grid%size(Z_INDEX)),     &
         dth_cond_evap_tot(current_state%local_grid%size(Z_INDEX)), &
         dqv_mphys_tot(current_state%local_grid%size(Z_INDEX)),     &
         dqv_cond_evap_tot(current_state%local_grid%size(Z_INDEX)), &
         dqc_mphys_tot(current_state%local_grid%size(Z_INDEX)),     &
         dqr_mphys_tot(current_state%local_grid%size(Z_INDEX)))
    if (.not. l_warm) then 
       allocate(phomc_tot(current_state%local_grid%size(Z_INDEX)), &
            pinuc_tot(current_state%local_grid%size(Z_INDEX)), &
            pidep_tot(current_state%local_grid%size(Z_INDEX)), &
            piacw_tot(current_state%local_grid%size(Z_INDEX)), &
            pisub_tot(current_state%local_grid%size(Z_INDEX)), &
            pimlt_tot(current_state%local_grid%size(Z_INDEX)), &
            psedi_tot(current_state%local_grid%size(Z_INDEX)), &
            psacw_tot(current_state%local_grid%size(Z_INDEX)), &
            psacr_tot(current_state%local_grid%size(Z_INDEX)), &
            pssub_tot(current_state%local_grid%size(Z_INDEX)), &
            psmlt_tot(current_state%local_grid%size(Z_INDEX)), &
            psaut_tot(current_state%local_grid%size(Z_INDEX)), & 
            psaci_tot(current_state%local_grid%size(Z_INDEX)), &
            psdep_tot(current_state%local_grid%size(Z_INDEX)), &
            pseds_tot(current_state%local_grid%size(Z_INDEX)), & 
            pgacw_tot(current_state%local_grid%size(Z_INDEX)), & 
            pgacs_tot(current_state%local_grid%size(Z_INDEX)), & 
            pgmlt_tot(current_state%local_grid%size(Z_INDEX)), & 
            pgsub_tot(current_state%local_grid%size(Z_INDEX)), & 
            psedg_tot(current_state%local_grid%size(Z_INDEX)), &
            dqi_mphys_tot(current_state%local_grid%size(Z_INDEX)), &
            dqs_mphys_tot(current_state%local_grid%size(Z_INDEX)), &
            dqg_mphys_tot(current_state%local_grid%size(Z_INDEX))) 
    endif
       
  end subroutine initialisation_callback  

  subroutine timestep_callback(current_state)
    type(model_state_type), target, intent(inout) :: current_state

    integer :: k, i
    integer :: icol, jcol, target_x_index, target_y_index
    
    icol=current_state%column_local_x
    jcol=current_state%column_local_y
    target_y_index=jcol-current_state%local_grid%halo_size(Y_INDEX)
    target_x_index=icol-current_state%local_grid%halo_size(X_INDEX)

    if (current_state%first_timestep_column) then
       psedl_tot(:)= 0.0_DEFAULT_PRECISION
       pcond_tot(:)= 0.0_DEFAULT_PRECISION
       praut_tot(:)= 0.0_DEFAULT_PRECISION 
       pracw_tot(:)= 0.0_DEFAULT_PRECISION  
       prevp_tot(:)= 0.0_DEFAULT_PRECISION 
       psedr_tot(:)= 0.0_DEFAULT_PRECISION
       dth_mphys_tot(:)= 0.0_DEFAULT_PRECISION
       dth_cond_evap_tot(:)= 0.0_DEFAULT_PRECISION
       dqv_mphys_tot(:)= 0.0_DEFAULT_PRECISION
       dqv_cond_evap_tot(:)= 0.0_DEFAULT_PRECISION
       dqc_mphys_tot(:)= 0.0_DEFAULT_PRECISION
       dqr_mphys_tot(:)= 0.0_DEFAULT_PRECISION
       if (.not. l_warm) then 
          phomc_tot(:)= 0.0_DEFAULT_PRECISION 
          pinuc_tot(:)= 0.0_DEFAULT_PRECISION 
          pidep_tot(:)= 0.0_DEFAULT_PRECISION 
          piacw_tot(:)= 0.0_DEFAULT_PRECISION 
          pisub_tot(:)= 0.0_DEFAULT_PRECISION 
          pimlt_tot(:)= 0.0_DEFAULT_PRECISION 
          psedi_tot(:)= 0.0_DEFAULT_PRECISION
          psacw_tot(:)= 0.0_DEFAULT_PRECISION 
          psacr_tot(:)= 0.0_DEFAULT_PRECISION 
          pssub_tot(:)= 0.0_DEFAULT_PRECISION 
          psmlt_tot(:)= 0.0_DEFAULT_PRECISION 
          psaut_tot(:)= 0.0_DEFAULT_PRECISION  
          psaci_tot(:)= 0.0_DEFAULT_PRECISION 
          psdep_tot(:)= 0.0_DEFAULT_PRECISION 
          pseds_tot(:)= 0.0_DEFAULT_PRECISION 
          pgacw_tot(:)= 0.0_DEFAULT_PRECISION  
          pgacs_tot(:)= 0.0_DEFAULT_PRECISION  
          pgmlt_tot(:)= 0.0_DEFAULT_PRECISION  
          pgsub_tot(:)= 0.0_DEFAULT_PRECISION
          psedg_tot(:)= 0.0_DEFAULT_PRECISION
          dqi_mphys_tot(:)= 0.0_DEFAULT_PRECISION
          dqs_mphys_tot(:)= 0.0_DEFAULT_PRECISION
          dqg_mphys_tot(:)= 0.0_DEFAULT_PRECISION 
       endif
    endif

    if (.not. current_state%halo_column) then
       if ( casdiags % l_psedl ) &
            psedl_tot(:)= psedl_tot(:) + casim_monc_dgs % psedl(:,target_y_index,target_x_index)
       if ( casdiags % l_pcond ) & 
            pcond_tot(:)= pcond_tot(:) + casim_monc_dgs % pcond(:,target_y_index,target_x_index)
       if (  casdiags % l_praut ) & 
            praut_tot(:)= praut_tot(:) + casim_monc_dgs % praut(:,target_y_index,target_x_index)
       if ( casdiags % l_pracw ) & 
            pracw_tot(:)= pracw_tot(:) + casim_monc_dgs % pracw(:,target_y_index,target_x_index)
       if ( casdiags % l_prevp ) & 
            prevp_tot(:)= prevp_tot(:) + casim_monc_dgs % prevp(:,target_y_index,target_x_index)
       if ( casdiags % l_psedr ) & 
            psedr_tot(:)= psedr_tot(:) + casim_monc_dgs % psedr(:,target_y_index,target_x_index)
       if ( casdiags % l_dth ) then
          dth_mphys_tot(:)= dth_mphys_tot(:) + &
               casim_monc_dgs %  dth_total(:,target_y_index,target_x_index)
          
          dth_cond_evap_tot(:)= dth_cond_evap_tot(:) + &
               casim_monc_dgs % dth_cond_evap(:,target_y_index,target_x_index)
       endif

       if ( casdiags % l_dqv ) then
          dqv_mphys_tot(:)= dqv_mphys_tot(:) + &
               casim_monc_dgs % dqv_total(:,target_y_index,target_x_index)
          dqv_cond_evap_tot(:)= dqv_cond_evap_tot(:) + &
               casim_monc_dgs % dqv_cond_evap(:,target_y_index,target_x_index)
       endif

       if ( casdiags % l_dqc ) &
            dqc_mphys_tot(:)= dqc_mphys_tot(:) + &
            casim_monc_dgs % dqc(:,target_y_index,target_x_index)
       if ( casdiags % l_dqr ) & 
            dqr_mphys_tot(:)= dqr_mphys_tot(:) + &
            casim_monc_dgs % dqr(:,target_y_index,target_x_index)
       
       if (.not. l_warm) then
          if ( casdiags % l_phomc ) & 
               phomc_tot(:)= phomc_tot(:) + &
               casim_monc_dgs % phomc(:,target_y_index,target_x_index)
          if ( casdiags % l_pinuc ) & 
               pinuc_tot(:)= pinuc_tot(:) + &
               casim_monc_dgs % pinuc(:,target_y_index,target_x_index)
          if ( casdiags % l_pidep ) & 
               pidep_tot(:)= pidep_tot(:) + &
               casim_monc_dgs % pidep(:,target_y_index,target_x_index)
          if ( casdiags % l_piacw ) & 
               piacw_tot(:)= piacw_tot(:) + &
               casim_monc_dgs % piacw(:,target_y_index,target_x_index)
          if ( casdiags % l_pisub ) &
               pisub_tot(:)= pisub_tot(:) + &
               casim_monc_dgs % pisub(:,target_y_index,target_x_index)
          if ( casdiags % l_pimlt ) &
               pimlt_tot(:)= pimlt_tot(:) + &
               casim_monc_dgs % pimlt(:,target_y_index,target_x_index)
          if ( casdiags % l_psedi ) & 
               psedi_tot(:)= psedi_tot(:) + &
               casim_monc_dgs % psedi(:,target_y_index,target_x_index)
          if ( casdiags % l_psacw ) & 
               psacw_tot(:)= psacw_tot(:) + &
               casim_monc_dgs % psacw(:,target_y_index,target_x_index)
          if ( casdiags % l_psacr ) & 
               psacr_tot(:)= psacr_tot(:) + &
               casim_monc_dgs % psacr(:,target_y_index,target_x_index)
          if ( casdiags % l_pssub ) & 
               pssub_tot(:)= pssub_tot(:) + &
               casim_monc_dgs % pssub(:,target_y_index,target_x_index)
          if ( casdiags % l_psmlt ) & 
               psmlt_tot(:)= psmlt_tot(:) + &
               casim_monc_dgs % psmlt(:,target_y_index,target_x_index)
          if ( casdiags % l_psaut ) &
               psaut_tot(:)= psaut_tot(:) + &
               casim_monc_dgs % psaut(:,target_y_index,target_x_index)
          if ( casdiags % l_psaci ) &
               psaci_tot(:)= psaci_tot(:) + &
               casim_monc_dgs % psaci(:,target_y_index,target_x_index)
          if ( casdiags % l_psdep ) & 
               psdep_tot(:)= psdep_tot(:) + &
               casim_monc_dgs % psdep(:,target_y_index,target_x_index)
          if ( casdiags % l_pseds ) & 
               pseds_tot(:)= pseds_tot(:) + &
               casim_monc_dgs % pseds(:,target_y_index,target_x_index)
          if ( casdiags % l_pgacw ) & 
               pgacw_tot(:)= pgacw_tot(:) + &
               casim_monc_dgs % pgacw(:,target_y_index,target_x_index)
          if ( casdiags % l_pgacs ) &
               pgacs_tot(:)= pgacs_tot(:) + &
               casim_monc_dgs % pgacs(:,target_y_index,target_x_index)
          if ( casdiags % l_pgmlt ) &
               pgmlt_tot(:)= pgmlt_tot(:) + &
               casim_monc_dgs % pgmlt(:,target_y_index,target_x_index)
          if ( casdiags % l_pgsub ) &
               pgsub_tot(:)= pgsub_tot(:) + &
               casim_monc_dgs % pgsub(:,target_y_index,target_x_index)
          if ( casdiags % l_psedg ) & 
               psedg_tot(:)= psedg_tot(:) + &
               casim_monc_dgs % psedg(:,target_y_index,target_x_index)
          if ( casdiags % l_dqi ) & 
               dqi_mphys_tot(:)= dqi_mphys_tot(:) + &
               casim_monc_dgs % dqi(:,target_y_index,target_x_index)
          if ( casdiags % l_dqs ) & 
               dqs_mphys_tot(:)= dqs_mphys_tot(:) + &
               casim_monc_dgs % dqs(:,target_y_index,target_x_index)
          if ( casdiags % l_dqg ) &
               dqg_mphys_tot(:)= dqg_mphys_tot(:) + &
               casim_monc_dgs % dqg(:,target_y_index,target_x_index)
       endif
       
    endif
  end subroutine timestep_callback  

  !> Field information retrieval callback, this returns information for a specific components published field
  !! @param current_state Current model state
  !! @param name The name of the field to retrieve information for
  !! @param field_information Populated with information about the field
  subroutine field_information_retrieval_callback(current_state, name, field_information)
    type(model_state_type), target, intent(inout) :: current_state
    character(len=*), intent(in) :: name
    type(component_field_information_type), intent(out) :: field_information

    field_information%field_type=COMPONENT_ARRAY_FIELD_TYPE
    field_information%number_dimensions=1
    field_information%dimension_sizes(1)=current_state%local_grid%size(Z_INDEX)
    field_information%data_type=COMPONENT_DOUBLE_DATA_TYPE
    if (l_warm) then
       if (name .eq. "pcond_total" .or. name .eq. "praut_total" &
            .or. name .eq. "pracw_total" .or. name .eq. "prevp_total" &
            .or. name .eq. "psedl_total" .or. name .eq. "psedr_total" &
            .or. name .eq. "dth_mphys_total" .or. name .eq. "dth_cond_evap_total" &
            .or. name .eq. "dqv_mphys_total" .or. name .eq. "dqv_cond_evap_total" &
            .or. name .eq. "dqc_mphys_total" .or. name .eq. "dqr_mphys_total") then
          field_information%enabled=.true.
       endif
    else
       field_information%enabled=.true.
    endif
           
  end subroutine field_information_retrieval_callback

  !> Field value retrieval callback, this returns the value of a specific published field
  !! @param current_state Current model state
  !! @param name The name of the field to retrieve the value for
  !! @param field_value Populated with the value of the field
  subroutine field_value_retrieval_callback(current_state, name, field_value)
    type(model_state_type), target, intent(inout) :: current_state
    character(len=*), intent(in) :: name
    type(component_field_value_type), intent(out) :: field_value
    
    integer :: k
    
    if (name .eq. "phomc_total") then
       allocate(field_value%real_1d_array(current_state%local_grid%size(Z_INDEX)))
       field_value%real_1d_array(:)=phomc_tot(:)
    else if (name .eq. "pinuc_total") then
       allocate(field_value%real_1d_array(current_state%local_grid%size(Z_INDEX)))
       field_value%real_1d_array(:)=pinuc_tot(:)
    else if (name .eq. "pidep_total") then
       allocate(field_value%real_1d_array(current_state%local_grid%size(Z_INDEX)))
       field_value%real_1d_array(:)=pidep_tot(:)
    else if (name .eq. "psdep_total") then
       allocate(field_value%real_1d_array(current_state%local_grid%size(Z_INDEX)))
       field_value%real_1d_array(:)=psdep_tot(:)
    else if (name .eq. "piacw_total") then
       allocate(field_value%real_1d_array(current_state%local_grid%size(Z_INDEX)))
       field_value%real_1d_array(:)=piacw_tot(:)
    else if (name .eq. "psacw_total") then
       allocate(field_value%real_1d_array(current_state%local_grid%size(Z_INDEX)))
       field_value%real_1d_array(:)=psacw_tot(:)
    else if (name .eq. "psacr_total") then
       allocate(field_value%real_1d_array(current_state%local_grid%size(Z_INDEX)))
       field_value%real_1d_array(:)=psacr_tot(:)
    else if (name .eq. "pisub_total") then
       allocate(field_value%real_1d_array(current_state%local_grid%size(Z_INDEX)))
       field_value%real_1d_array(:)=pisub_tot(:)
    else if (name .eq. "pssub_total") then
       allocate(field_value%real_1d_array(current_state%local_grid%size(Z_INDEX)))
       field_value%real_1d_array(:)=pssub_tot(:)
    else if (name .eq. "pimlt_total") then
       allocate(field_value%real_1d_array(current_state%local_grid%size(Z_INDEX)))
       field_value%real_1d_array(:)=pimlt_tot(:)
    else if (name .eq. "psmlt_total") then
       allocate(field_value%real_1d_array(current_state%local_grid%size(Z_INDEX)))
       field_value%real_1d_array(:)=psmlt_tot(:)
    else if (name .eq. "psaut_total") then
       allocate(field_value%real_1d_array(current_state%local_grid%size(Z_INDEX)))
       field_value%real_1d_array(:)=psaut_tot(:)
    else if (name .eq. "psaci_total") then
       allocate(field_value%real_1d_array(current_state%local_grid%size(Z_INDEX)))
       field_value%real_1d_array(:)=psaci_tot(:)
    else if (name .eq. "praut_total") then
       allocate(field_value%real_1d_array(current_state%local_grid%size(Z_INDEX)))
       field_value%real_1d_array(:)=praut_tot(:)
    else if (name .eq. "pracw_total") then
       allocate(field_value%real_1d_array(current_state%local_grid%size(Z_INDEX)))
       field_value%real_1d_array(:)=pracw_tot(:)
    else if (name .eq. "prevp_total") then
       allocate(field_value%real_1d_array(current_state%local_grid%size(Z_INDEX)))
       field_value%real_1d_array(:)=prevp_tot(:)
    else if (name .eq. "pgacw_total") then
       allocate(field_value%real_1d_array(current_state%local_grid%size(Z_INDEX)))
       field_value%real_1d_array(:)=pgacw_tot(:)
    else if (name .eq. "pgacs_total") then
       allocate(field_value%real_1d_array(current_state%local_grid%size(Z_INDEX)))
       field_value%real_1d_array(:)=pgacs_tot(:)
    else if (name .eq. "pgmlt_total") then
       allocate(field_value%real_1d_array(current_state%local_grid%size(Z_INDEX)))
       field_value%real_1d_array(:)=pgmlt_tot(:)
    else if (name .eq. "pgsub_total") then
       allocate(field_value%real_1d_array(current_state%local_grid%size(Z_INDEX)))
       field_value%real_1d_array(:)=pgsub_tot(:)
    else if (name .eq. "psedi_total") then
       allocate(field_value%real_1d_array(current_state%local_grid%size(Z_INDEX)))
       field_value%real_1d_array(:)=psedi_tot(:)
    else if (name .eq. "pseds_total") then
       allocate(field_value%real_1d_array(current_state%local_grid%size(Z_INDEX)))
       field_value%real_1d_array(:)=pseds_tot(:)
    else if (name .eq. "psedr_total") then
       allocate(field_value%real_1d_array(current_state%local_grid%size(Z_INDEX)))
       field_value%real_1d_array(:)=psedr_tot(:)
    else if (name .eq. "psedg_total") then
       allocate(field_value%real_1d_array(current_state%local_grid%size(Z_INDEX)))
       field_value%real_1d_array(:)=psedg_tot(:)
    else if (name .eq. "psedl_total") then
       allocate(field_value%real_1d_array(current_state%local_grid%size(Z_INDEX)))
       field_value%real_1d_array(:)=psedl_tot(:)
    else if (name .eq. "pcond_total") then
       allocate(field_value%real_1d_array(current_state%local_grid%size(Z_INDEX)))
       field_value%real_1d_array(:)=pcond_tot(:)
    else if (name .eq. "dth_mphys_total") then
       allocate(field_value%real_1d_array(current_state%local_grid%size(Z_INDEX)))
       field_value%real_1d_array(:)=dth_mphys_tot(:)
    else if (name .eq. "dth_cond_evap_total") then
       allocate(field_value%real_1d_array(current_state%local_grid%size(Z_INDEX)))
       field_value%real_1d_array(:)=dth_cond_evap_tot(:)
    else if (name .eq. "dqv_mphys_total") then
       allocate(field_value%real_1d_array(current_state%local_grid%size(Z_INDEX)))
       field_value%real_1d_array(:)=dqv_mphys_tot(:)
    else if (name .eq. "dqv_cond_evap_total") then
       allocate(field_value%real_1d_array(current_state%local_grid%size(Z_INDEX)))
       field_value%real_1d_array(:)=dqv_cond_evap_tot(:)
    else if (name .eq. "dqc_mphys_total") then
       allocate(field_value%real_1d_array(current_state%local_grid%size(Z_INDEX)))
       field_value%real_1d_array(:)=dqc_mphys_tot(:)
    else if (name .eq. "dqr_mphys_total") then
       allocate(field_value%real_1d_array(current_state%local_grid%size(Z_INDEX)))
       field_value%real_1d_array(:)=dqr_mphys_tot(:)
    else if (name .eq. "dqi_mphys_total") then
       allocate(field_value%real_1d_array(current_state%local_grid%size(Z_INDEX)))
       field_value%real_1d_array(:)=dqi_mphys_tot(:)
    else if (name .eq. "dqs_mphys_total") then
       allocate(field_value%real_1d_array(current_state%local_grid%size(Z_INDEX)))
       field_value%real_1d_array(:)=dqs_mphys_tot(:)
    else if (name .eq. "dqg_mphys_total") then
       allocate(field_value%real_1d_array(current_state%local_grid%size(Z_INDEX)))
       field_value%real_1d_array(:)=dqg_mphys_tot(:) 
    endif
  end subroutine field_value_retrieval_callback
end module casim_profile_dgs_mod
