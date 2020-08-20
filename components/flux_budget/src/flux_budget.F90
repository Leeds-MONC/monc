!> Flux budget component which produces diagnostic data for the flux aspects of the model
module flux_budget_mod
  use datadefn_mod, only : DEFAULT_PRECISION
  use collections_mod, only : hashmap_type, mapentry_type, iterator_type, c_contains, c_put_logical, c_size, c_get_logical, &
       c_get_iterator, c_has_next, c_next_mapentry
  use monc_component_mod, only : COMPONENT_ARRAY_FIELD_TYPE, COMPONENT_DOUBLE_DATA_TYPE, component_descriptor_type, &
       component_field_value_type, component_field_information_type
  use optionsdatabase_mod, only : options_get_real, options_get_integer
  use registry_mod, only : get_component_field_value, get_component_field_information, is_component_field_available
  use grids_mod, only : local_grid_type, global_grid_type, X_INDEX, Y_INDEX, Z_INDEX
  use state_mod, only : model_state_type
  use science_constants_mod, only: G, ratio_mol_wts
  implicit none

#ifndef TEST_MODE
  private
#endif

  real(kind=DEFAULT_PRECISION), dimension(:), allocatable :: th_flux_values, th_gradient, th_diff, th_buoyancy, th_tendency, &
       uw_advection, vw_advection, uw_viscosity, vw_viscosity, uw_buoyancy, vw_buoyancy, uw_tendency, vw_tendency, uw_w, &
       vw_w, tu_su, uu_advection, uu_viscosity, wu_u, tv_sv, vv_advection, vv_viscosity, wv_v, tw_sw, ww_advection, &
       ww_viscosity, ww_buoyancy, u_thetal, us_thetal, u_thetal_advection, u_thetal_viscosity_diffusion, wu_thetal, &
       v_thetal, vs_thetal, v_thetal_advection, v_thetal_viscosity_diffusion, wv_thetal, w_thetal, ws_thetal, &
       w_thetal_advection, w_thetal_viscosity_diffusion, w_thetal_buoyancy, ww_thetal, thetal_thetal, sthetal_thetal, &
       thetal_thetal_advection, thetal_thetal_diffusion, wthetal_thetal, u_mse, us_mse, u_mse_advection, &
       u_mse_viscosity_diffusion, wu_mse, v_mse, vs_mse, v_mse_advection, v_mse_viscosity_diffusion, wv_mse, w_mse, ws_mse, &
       w_mse_advection, w_mse_viscosity_diffusion, w_mse_buoyancy, ww_mse, mse_mse, smse_mse, mse_mse_advection, &
       mse_mse_diffusion, wmse_mse, us_qt, u_qt_advection, u_qt_viscosity_diffusion, wu_qt, vs_qt, v_qt_advection, &
       v_qt_viscosity_diffusion, wv_qt, w_qt, ws_qt, w_qt_advection, w_qt_viscosity_diffusion, w_qt_buoyancy, &
       ww_qt, qt_qt, sqt_qt, qt_qt_advection, qt_qt_diffusion, wqt_qt, sres, wke, buoy, wp, tend
  real(kind=DEFAULT_PRECISION) :: mflux, wmfcrit
  real(kind=DEFAULT_PRECISION), dimension(:,:), allocatable :: q_flux_values, q_gradient, q_diff, q_buoyancy, q_tendency
  type(hashmap_type) :: heat_flux_fields, q_flux_fields, uw_vw_fields, prognostic_budget_fields, thetal_fields, mse_fields, &
       qt_fields, scalar_fields, tke_fields

  logical :: some_theta_flux_diagnostics_enabled, some_q_flux_diagnostics_enabled, some_uw_vw_diagnostics_enabled, &
       some_prognostic_budget_diagnostics_enabled, some_thetal_diagnostics_enabled, some_mse_diagnostics_enabled, &
       some_qt_diagnostics_enabled, some_tke_diagnostics_enabled
       
  public flux_budget_get_descriptor
contains

  !> Provides the descriptor back to the caller and is used in component registration
  !! @returns The flux budget component descriptor
  type(component_descriptor_type) function flux_budget_get_descriptor()
    type(iterator_type) :: iterator
    type(mapentry_type) :: mapentry
    integer :: current_index, total_number_published_fields

    flux_budget_get_descriptor%name="flux_budget"
    flux_budget_get_descriptor%version=0.1
    flux_budget_get_descriptor%initialisation=>initialisation_callback
    flux_budget_get_descriptor%timestep=>timestep_callback
    flux_budget_get_descriptor%finalisation=>finalisation_callback

    flux_budget_get_descriptor%field_value_retrieval=>field_value_retrieval_callback
    flux_budget_get_descriptor%field_information_retrieval=>field_information_retrieval_callback
    call populate_field_names()
    total_number_published_fields=c_size(heat_flux_fields) + c_size(q_flux_fields) + c_size(uw_vw_fields) +       &
                                  c_size(prognostic_budget_fields) + c_size(thetal_fields) + c_size(mse_fields) + &
                                  c_size(qt_fields) + c_size(scalar_fields) + c_size(tke_fields)
    allocate(flux_budget_get_descriptor%published_fields(total_number_published_fields))

    current_index=1
    iterator=c_get_iterator(heat_flux_fields)
    do while (c_has_next(iterator))
      mapentry=c_next_mapentry(iterator)
      flux_budget_get_descriptor%published_fields(current_index)=mapentry%key
      current_index=current_index+1
    end do    
    iterator=c_get_iterator(q_flux_fields)
    do while (c_has_next(iterator))
      mapentry=c_next_mapentry(iterator)
      flux_budget_get_descriptor%published_fields(current_index)=mapentry%key
      current_index=current_index+1
    end do  
    iterator=c_get_iterator(uw_vw_fields)
    do while (c_has_next(iterator))
      mapentry=c_next_mapentry(iterator)
      flux_budget_get_descriptor%published_fields(current_index)=mapentry%key
      current_index=current_index+1
    end do 
    iterator=c_get_iterator(tke_fields)
    do while (c_has_next(iterator))
      mapentry=c_next_mapentry(iterator)
      flux_budget_get_descriptor%published_fields(current_index)=mapentry%key
      current_index=current_index+1
    end do 
    iterator=c_get_iterator(prognostic_budget_fields)
    do while (c_has_next(iterator))
      mapentry=c_next_mapentry(iterator)
      flux_budget_get_descriptor%published_fields(current_index)=mapentry%key
      current_index=current_index+1
    end do  
    iterator=c_get_iterator(thetal_fields)
    do while (c_has_next(iterator))
      mapentry=c_next_mapentry(iterator)
      flux_budget_get_descriptor%published_fields(current_index)=mapentry%key
      current_index=current_index+1
    end do  
    iterator=c_get_iterator(mse_fields)
    do while (c_has_next(iterator))
      mapentry=c_next_mapentry(iterator)
      flux_budget_get_descriptor%published_fields(current_index)=mapentry%key
      current_index=current_index+1
    end do  
    iterator=c_get_iterator(qt_fields)
    do while (c_has_next(iterator))
      mapentry=c_next_mapentry(iterator)
      flux_budget_get_descriptor%published_fields(current_index)=mapentry%key
      current_index=current_index+1
    end do 
    iterator=c_get_iterator(scalar_fields)
    do while (c_has_next(iterator))
      mapentry=c_next_mapentry(iterator)
      flux_budget_get_descriptor%published_fields(current_index)=mapentry%key
      current_index=current_index+1
    end do    
  end function flux_budget_get_descriptor 

  !> Initialisation call back
  !! @param current_state Current model state
  subroutine initialisation_callback(current_state)
    type(model_state_type), target, intent(inout) :: current_state

    call initialise_theta_flux_diagnostics(current_state)
    call initialise_q_flux_diagnostics(current_state)
    call initialise_uw_vw_diagnostics(current_state)
    call initialise_prognostic_budget_diagnostics(current_state)
    call initialise_thetal_diagnostics(current_state)
    call initialise_mse_diagnostics(current_state)
    call initialise_qt_diagnostics(current_state)
    call initialise_scalar_diagnostics(current_state)
    call initialise_tke_diagnostics(current_state)
  end subroutine initialisation_callback

  !> Timestep call back, this will deduce the diagnostics for the current (non halo) column
  !! @param current_state Current model state
  subroutine timestep_callback(current_state)
    type(model_state_type), target, intent(inout) :: current_state
    logical :: calculate_diagnostics

    calculate_diagnostics = current_state%diagnostic_sample_timestep

    if (calculate_diagnostics) then
      if (current_state%first_timestep_column) then
        if (some_theta_flux_diagnostics_enabled) call clear_theta_fluxes()
        if (some_q_flux_diagnostics_enabled) call clear_q_fluxes()
        if (some_uw_vw_diagnostics_enabled) call clear_uw_vw()
        if (some_prognostic_budget_diagnostics_enabled) call clear_prognostic_budgets()
        if (some_thetal_diagnostics_enabled) call clear_thetal()
        if (some_mse_diagnostics_enabled) call clear_mse()
        if (some_qt_diagnostics_enabled) call clear_qt()
        if (some_tke_diagnostics_enabled) call clear_tke()
        call clear_scalars()
      end if
      if (.not. current_state%halo_column) then
        if (some_theta_flux_diagnostics_enabled) call compute_theta_flux_for_column(current_state)
        if (some_q_flux_diagnostics_enabled) call compute_q_flux_for_column(current_state)
        if (some_uw_vw_diagnostics_enabled) call compute_uw_vw_for_column(current_state)
        if (some_prognostic_budget_diagnostics_enabled) call compute_prognostic_budgets_for_column(current_state)
        if (some_thetal_diagnostics_enabled) call compute_thetal_for_column(current_state)
        if (some_mse_diagnostics_enabled) call compute_mse_for_column(current_state)
        if (some_qt_diagnostics_enabled) call compute_qt_for_column(current_state)
        if (some_tke_diagnostics_enabled) call compute_tke_for_column(current_state)
        call compute_scalars_for_column(current_state)
      end if
    end if
  end subroutine timestep_callback

  !> Finalisation call back
  !! @param current_state Current model state
  subroutine finalisation_callback(current_state)
    type(model_state_type), target, intent(inout) :: current_state

    if (allocated(th_flux_values)) deallocate(th_flux_values)
    if (allocated(th_gradient)) deallocate(th_gradient)
    if (allocated(th_diff)) deallocate(th_diff)
    if (allocated(th_buoyancy)) deallocate(th_buoyancy)
    if (allocated(th_tendency)) deallocate(th_tendency)
    if (allocated(q_flux_values)) deallocate(q_flux_values)
    if (allocated(q_gradient)) deallocate(q_gradient)
    if (allocated(q_diff)) deallocate(q_diff)
    if (allocated(q_buoyancy)) deallocate(q_buoyancy)
    if (allocated(q_tendency)) deallocate(q_tendency)

    if (allocated(uw_advection)) deallocate(uw_advection)
    if (allocated(vw_advection)) deallocate(vw_advection)
    if (allocated(uw_viscosity)) deallocate(uw_viscosity)
    if (allocated(vw_viscosity)) deallocate(vw_viscosity)
    if (allocated(uw_buoyancy)) deallocate(uw_buoyancy)
    if (allocated(vw_buoyancy)) deallocate(vw_buoyancy)
    if (allocated(uw_tendency)) deallocate(uw_tendency)
    if (allocated(vw_tendency)) deallocate(vw_tendency)
    if (allocated(uw_w)) deallocate(uw_w)
    if (allocated(vw_w)) deallocate(vw_w)

    if (allocated(sres)) deallocate(sres)
    if (allocated(buoy)) deallocate(buoy)
    if (allocated(wke)) deallocate(wke)
    if (allocated(wp)) deallocate(wp)
    if (allocated(tend)) deallocate(tend)

    if (allocated(tu_su)) deallocate(tu_su)
    if (allocated(uu_advection)) deallocate(uu_advection)
    if (allocated(uu_viscosity)) deallocate(uu_viscosity)
    if (allocated(wu_u)) deallocate(wu_u)
    if (allocated(tv_sv)) deallocate(tv_sv)
    if (allocated(vv_advection)) deallocate(vv_advection)
    if (allocated(vv_viscosity)) deallocate(vv_viscosity)
    if (allocated(wv_v)) deallocate(wv_v)
    if (allocated(tw_sw)) deallocate(tw_sw)
    if (allocated(ww_advection)) deallocate(ww_advection)
    if (allocated(ww_viscosity)) deallocate(ww_viscosity)
    if (allocated(ww_buoyancy)) deallocate(ww_buoyancy)

    if (allocated(u_thetal)) deallocate(u_thetal)
    if (allocated(us_thetal)) deallocate(us_thetal)
    if (allocated(u_thetal_advection)) deallocate(u_thetal_advection)
    if (allocated(u_thetal_viscosity_diffusion)) deallocate(u_thetal_viscosity_diffusion)
    if (allocated(wu_thetal)) deallocate(wu_thetal)
    if (allocated(v_thetal)) deallocate(v_thetal)
    if (allocated(vs_thetal)) deallocate(vs_thetal)
    if (allocated(v_thetal_advection)) deallocate(v_thetal_advection)
    if (allocated(v_thetal_viscosity_diffusion)) deallocate(v_thetal_viscosity_diffusion)
    if (allocated(wv_thetal)) deallocate(wv_thetal)
    if (allocated(w_thetal)) deallocate(w_thetal)
    if (allocated(ws_thetal)) deallocate(ws_thetal)
    if (allocated(w_thetal_advection)) deallocate(w_thetal_advection)
    if (allocated(w_thetal_viscosity_diffusion)) deallocate(w_thetal_viscosity_diffusion)
    if (allocated(w_thetal_buoyancy)) deallocate(w_thetal_buoyancy)
    if (allocated(ww_thetal)) deallocate(ww_thetal)
    if (allocated(thetal_thetal)) deallocate(thetal_thetal)
    if (allocated(sthetal_thetal)) deallocate(sthetal_thetal)
    if (allocated(thetal_thetal_advection)) deallocate(thetal_thetal_advection)
    if (allocated(thetal_thetal_diffusion)) deallocate(thetal_thetal_diffusion)
    if (allocated(wthetal_thetal)) deallocate(wthetal_thetal)

    if (allocated(u_mse)) deallocate(u_mse)
    if (allocated(us_mse)) deallocate(us_mse)
    if (allocated(u_mse_advection)) deallocate(u_mse_advection)
    if (allocated(u_mse_viscosity_diffusion)) deallocate(u_mse_viscosity_diffusion)
    if (allocated(wu_mse)) deallocate(wu_mse)
    if (allocated(v_mse)) deallocate(v_mse)
    if (allocated(vs_mse)) deallocate(vs_mse)
    if (allocated(v_mse_advection)) deallocate(v_mse_advection)
    if (allocated(v_mse_viscosity_diffusion)) deallocate(v_mse_viscosity_diffusion)
    if (allocated(wv_mse)) deallocate(wv_mse)
    if (allocated(w_mse)) deallocate(w_mse)
    if (allocated(ws_mse)) deallocate(ws_mse)
    if (allocated(w_mse_advection)) deallocate(w_mse_advection)
    if (allocated(w_mse_viscosity_diffusion)) deallocate(w_mse_viscosity_diffusion)
    if (allocated(w_mse_buoyancy)) deallocate(w_mse_buoyancy)
    if (allocated(ww_mse)) deallocate(ww_mse)
    if (allocated(mse_mse)) deallocate(mse_mse)
    if (allocated(smse_mse)) deallocate(smse_mse)
    if (allocated(mse_mse_advection)) deallocate(mse_mse_advection)
    if (allocated(mse_mse_diffusion)) deallocate(mse_mse_diffusion)
    if (allocated(wmse_mse)) deallocate(wmse_mse)

    if (allocated(us_qt)) deallocate(us_qt)
    if (allocated(u_qt_advection)) deallocate(u_qt_advection)
    if (allocated(u_qt_viscosity_diffusion)) deallocate(u_qt_viscosity_diffusion)
    if (allocated(wu_qt)) deallocate(wu_qt)
    if (allocated(vs_qt)) deallocate(vs_qt)
    if (allocated(v_qt_advection)) deallocate(v_qt_advection)
    if (allocated(v_qt_viscosity_diffusion)) deallocate(v_qt_viscosity_diffusion)
    if (allocated(wv_qt)) deallocate(wv_qt)
    if (allocated(w_qt)) deallocate(w_qt)
    if (allocated(ws_qt)) deallocate(ws_qt)
    if (allocated(w_qt_advection)) deallocate(w_qt_advection)
    if (allocated(w_qt_viscosity_diffusion)) deallocate(w_qt_viscosity_diffusion)
    if (allocated(w_qt_buoyancy)) deallocate(w_qt_buoyancy)
    if (allocated(ww_qt)) deallocate(ww_qt)
    if (allocated(qt_qt)) deallocate(qt_qt)
    if (allocated(sqt_qt)) deallocate(sqt_qt)
    if (allocated(qt_qt_advection)) deallocate(qt_qt_advection)
    if (allocated(qt_qt_diffusion)) deallocate(qt_qt_diffusion)
    if (allocated(wqt_qt)) deallocate(wqt_qt)
  end subroutine finalisation_callback

    !> Populates the published field names in the appropriate map
  subroutine populate_field_names()
    call set_published_field_enabled_state(heat_flux_fields, "heat_flux_transport_local", .false.)
    call set_published_field_enabled_state(heat_flux_fields, "heat_flux_gradient_local", .false.)
    call set_published_field_enabled_state(heat_flux_fields, "heat_flux_dissipation_local", .false.)
    call set_published_field_enabled_state(heat_flux_fields, "heat_flux_buoyancy_local", .false.)
    call set_published_field_enabled_state(heat_flux_fields, "heat_flux_tendency_local", .false.)

    call set_published_field_enabled_state(q_flux_fields, "q_flux_transport_local", .false.)
    call set_published_field_enabled_state(q_flux_fields, "q_flux_gradient_local", .false.)
    call set_published_field_enabled_state(q_flux_fields, "q_flux_dissipation_local", .false.)
    call set_published_field_enabled_state(q_flux_fields, "q_flux_buoyancy_local", .false.)
    call set_published_field_enabled_state(q_flux_fields, "q_flux_tendency_local", .false.)

    call set_published_field_enabled_state(tke_fields, "resolved_shear_production_local", .false.)
    call set_published_field_enabled_state(tke_fields, "resolved_buoyant_production_local", .false.)
    call set_published_field_enabled_state(tke_fields, "resolved_turbulent_transport_local", .false.)
    call set_published_field_enabled_state(tke_fields, "resolved_pressure_transport_local", .false.)
    call set_published_field_enabled_state(tke_fields, "tke_tendency_local", .false.)

    call set_published_field_enabled_state(uw_vw_fields, "uw_advection_local", .false.)
    call set_published_field_enabled_state(uw_vw_fields, "vw_advection_local", .false.)
    call set_published_field_enabled_state(uw_vw_fields, "uw_viscosity_local", .false.)
    call set_published_field_enabled_state(uw_vw_fields, "vw_viscosity_local", .false.)
    call set_published_field_enabled_state(uw_vw_fields, "uw_buoyancy_local", .false.)
    call set_published_field_enabled_state(uw_vw_fields, "vw_buoyancy_local", .false.)
    call set_published_field_enabled_state(uw_vw_fields, "uw_tendency_local", .false.)
    call set_published_field_enabled_state(uw_vw_fields, "vw_tendency_local", .false.)
    call set_published_field_enabled_state(uw_vw_fields, "uw_w_local", .false.)
    call set_published_field_enabled_state(uw_vw_fields, "vw_w_local", .false.)

    call set_published_field_enabled_state(prognostic_budget_fields, "tu_su_local", .false.)
    call set_published_field_enabled_state(prognostic_budget_fields, "uu_advection_local", .false.)
    call set_published_field_enabled_state(prognostic_budget_fields, "uu_viscosity_local", .false.)
    call set_published_field_enabled_state(prognostic_budget_fields, "wu_u_local", .false.)
    call set_published_field_enabled_state(prognostic_budget_fields, "tv_sv_local", .false.)
    call set_published_field_enabled_state(prognostic_budget_fields, "vv_advection_local", .false.)
    call set_published_field_enabled_state(prognostic_budget_fields, "vv_viscosity_local", .false.)
    call set_published_field_enabled_state(prognostic_budget_fields, "wv_v_local", .false.)
    call set_published_field_enabled_state(prognostic_budget_fields, "tw_sw_local", .false.)
    call set_published_field_enabled_state(prognostic_budget_fields, "ww_advection_local", .false.)
    call set_published_field_enabled_state(prognostic_budget_fields, "ww_viscosity_local", .false.)
    call set_published_field_enabled_state(prognostic_budget_fields, "ww_buoyancy_local", .false.)

    call set_published_field_enabled_state(thetal_fields, "u_thetal_local", .false.)
    call set_published_field_enabled_state(thetal_fields, "us_thetal_local", .false.)
    call set_published_field_enabled_state(thetal_fields, "u_thetal_advection_local", .false.)
    call set_published_field_enabled_state(thetal_fields, "u_thetal_viscosity_diffusion_local", .false.)
    call set_published_field_enabled_state(thetal_fields, "wu_thetal_local", .false.)
    call set_published_field_enabled_state(thetal_fields, "v_thetal_local", .false.)
    call set_published_field_enabled_state(thetal_fields, "vs_thetal_local", .false.)
    call set_published_field_enabled_state(thetal_fields, "v_thetal_advection_local", .false.)
    call set_published_field_enabled_state(thetal_fields, "v_thetal_viscosity_diffusion_local", .false.)
    call set_published_field_enabled_state(thetal_fields, "wv_thetal_local", .false.)
    call set_published_field_enabled_state(thetal_fields, "w_thetal_local", .false.)
    call set_published_field_enabled_state(thetal_fields, "ws_thetal_local", .false.)
    call set_published_field_enabled_state(thetal_fields, "w_thetal_advection_local", .false.)
    call set_published_field_enabled_state(thetal_fields, "w_thetal_viscosity_diffusion_local", .false.)
    call set_published_field_enabled_state(thetal_fields, "ww_thetal_buoyancy_local", .false.)
    call set_published_field_enabled_state(thetal_fields, "ww_thetal_local", .false.)
    call set_published_field_enabled_state(thetal_fields, "thetal_thetal_local", .false.)
    call set_published_field_enabled_state(thetal_fields, "sthetal_thetal_local", .false.)
    call set_published_field_enabled_state(thetal_fields, "thetal_thetal_advection_local", .false.)
    call set_published_field_enabled_state(thetal_fields, "thetal_thetal_diffusion_local", .false.)
    call set_published_field_enabled_state(thetal_fields, "wthetal_thetal_local", .false.)

    call set_published_field_enabled_state(mse_fields, "u_mse_local", .false.)
    call set_published_field_enabled_state(mse_fields, "us_mse_local", .false.)
    call set_published_field_enabled_state(mse_fields, "u_mse_advection_local", .false.)
    call set_published_field_enabled_state(mse_fields, "u_mse_viscosity_diffusion_local", .false.)
    call set_published_field_enabled_state(mse_fields, "wu_mse_local", .false.)
    call set_published_field_enabled_state(mse_fields, "v_mse_local", .false.)
    call set_published_field_enabled_state(mse_fields, "vs_mse_local", .false.)
    call set_published_field_enabled_state(mse_fields, "v_mse_advection_local", .false.)
    call set_published_field_enabled_state(mse_fields, "v_mse_viscosity_diffusion_local", .false.)
    call set_published_field_enabled_state(mse_fields, "wv_mse_local", .false.)
    call set_published_field_enabled_state(mse_fields, "w_mse_local", .false.)
    call set_published_field_enabled_state(mse_fields, "ws_mse_local", .false.)
    call set_published_field_enabled_state(mse_fields, "w_mse_advection_local", .false.)
    call set_published_field_enabled_state(mse_fields, "w_mse_viscosity_diffusion_local", .false.)
    call set_published_field_enabled_state(mse_fields, "ww_mse_buoyancy_local", .false.)
    call set_published_field_enabled_state(mse_fields, "ww_mse_local", .false.)
    call set_published_field_enabled_state(mse_fields, "mse_mse_local", .false.)
    call set_published_field_enabled_state(mse_fields, "smse_mse_local", .false.)
    call set_published_field_enabled_state(mse_fields, "mse_mse_advection_local", .false.)
    call set_published_field_enabled_state(mse_fields, "mse_mse_diffusion_local", .false.)
    call set_published_field_enabled_state(mse_fields, "wmse_mse_local", .false.)

    call set_published_field_enabled_state(qt_fields, "us_qt_local", .false.)
    call set_published_field_enabled_state(qt_fields, "u_qt_advection_local", .false.)
    call set_published_field_enabled_state(qt_fields, "u_qt_viscosity_diffusion_local", .false.)
    call set_published_field_enabled_state(qt_fields, "wu_qt_local", .false.)
    call set_published_field_enabled_state(qt_fields, "vs_qt_local", .false.)
    call set_published_field_enabled_state(qt_fields, "v_qt_advection_local", .false.)
    call set_published_field_enabled_state(qt_fields, "v_qt_viscosity_diffusion_local", .false.)
    call set_published_field_enabled_state(qt_fields, "wv_qt_local", .false.)
    call set_published_field_enabled_state(qt_fields, "w_qt_local", .false.)
    call set_published_field_enabled_state(qt_fields, "ws_qt_local", .false.)
    call set_published_field_enabled_state(qt_fields, "w_qt_advection_local", .false.)
    call set_published_field_enabled_state(qt_fields, "w_qt_viscosity_diffusion_local", .false.)
    call set_published_field_enabled_state(qt_fields, "ww_qt_buoyancy_local", .false.)
    call set_published_field_enabled_state(qt_fields, "ww_qt_local", .false.)
    call set_published_field_enabled_state(qt_fields, "qt_qt_local", .false.)
    call set_published_field_enabled_state(qt_fields, "sqt_qt_local", .false.)
    call set_published_field_enabled_state(qt_fields, "qt_qt_advection_local", .false.)
    call set_published_field_enabled_state(qt_fields, "qt_qt_diffusion_local", .false.)
    call set_published_field_enabled_state(qt_fields, "wqt_qt_local", .false.)

    call set_published_field_enabled_state(scalar_fields, "mflux_local", .true.)
  end subroutine populate_field_names 

  !> Field information retrieval callback, this returns information for a specific components published field
  !! @param current_state Current model state
  !! @param name The name of the field to retrieve information for
  !! @param field_information Populated with information about the field
  subroutine field_information_retrieval_callback(current_state, name, field_information)
    type(model_state_type), target, intent(inout) :: current_state
    character(len=*), intent(in) :: name
    type(component_field_information_type), intent(out) :: field_information

    ! Field description is the same regardless of the specific field being retrieved
    field_information%field_type=COMPONENT_ARRAY_FIELD_TYPE
    field_information%data_type=COMPONENT_DOUBLE_DATA_TYPE

    if (is_field_heat_flux(name) .or. is_field_uw_vw(name) .or. is_field_prognostic_budget(name) &
        .or. is_field_thetal(name) .or. is_field_mse(name) .or. is_field_qt(name) .or. is_field_scalar(name) &
        .or. is_field_tke_flux(name)) then
      field_information%number_dimensions=1
      field_information%dimension_sizes(1)=current_state%local_grid%size(Z_INDEX)      
      if (is_field_heat_flux(name)) then
        field_information%enabled=get_published_field_enabled_state(heat_flux_fields, name)
      else if (is_field_uw_vw(name)) then
        field_information%enabled=get_published_field_enabled_state(uw_vw_fields, name)
      else if (is_field_prognostic_budget(name)) then
        field_information%enabled=get_published_field_enabled_state(prognostic_budget_fields, name)
      else if (is_field_thetal(name)) then
        field_information%enabled=get_published_field_enabled_state(thetal_fields, name)
      else if (is_field_tke_flux(name)) then
        field_information%enabled=get_published_field_enabled_state(tke_fields, name)
      else if (is_field_mse(name)) then
        field_information%enabled=get_published_field_enabled_state(mse_fields, name)
      else if (is_field_qt(name)) then
        field_information%enabled=get_published_field_enabled_state(qt_fields, name)
      else if (is_field_scalar(name)) then
        field_information%enabled=get_published_field_enabled_state(scalar_fields, name)
      end if
    else if (is_field_q_flux(name)) then
      field_information%number_dimensions=2
      field_information%dimension_sizes(1)=current_state%local_grid%size(Z_INDEX)
      field_information%dimension_sizes(2)=current_state%number_q_fields
      field_information%enabled=get_published_field_enabled_state(q_flux_fields, name)     
    end if
  end subroutine field_information_retrieval_callback

  !> Initialises the scalar diagnostics
  !! @param current_state The current model state
  subroutine initialise_scalar_diagnostics(current_state)
    type(model_state_type), target, intent(inout) :: current_state

    wmfcrit=options_get_real(current_state%options_database, "wmfcrit")
  end subroutine initialise_scalar_diagnostics  

  !> Initialises the qt diagnostics. For now we are assuming qt is the same as theta, which needs updating with moisture
  !! information
  !! @param current_state The current model state
  subroutine initialise_qt_diagnostics(current_state)
    type(model_state_type), target, intent(inout) :: current_state

    integer :: column_size
    logical :: us_qt_enabled, u_qt_advection_enabled, u_qt_viscosity_diffusion_enabled, &
         wu_qt_enabled, vs_qt_enabled, v_qt_advection_enabled, &
         v_qt_viscosity_diffusion_enabled, wv_qt_enabled, w_qt_enabled, ws_qt_enabled, &
         w_qt_advection_enabled, w_qt_viscosity_diffusion_enabled, w_qt_buoyancy_enabled, ww_qt_enabled, &
         qt_qt_enabled, sqt_qt_enabled, qt_qt_advection_enabled, &
         qt_qt_diffusion_enabled, wqt_qt_enabled

    column_size=current_state%local_grid%size(Z_INDEX)

    us_qt_enabled=current_state%u%active .and. current_state%th%active
    u_qt_advection_enabled=is_component_field_available("u_advection") .and. &
         is_component_field_available("th_advection") .and. us_qt_enabled
    u_qt_viscosity_diffusion_enabled=is_component_field_available("u_viscosity") .and. &
         is_component_field_available("th_diffusion") .and. us_qt_enabled
    wu_qt_enabled=current_state%w%active .and. us_qt_enabled

    some_qt_diagnostics_enabled=us_qt_enabled

    if (us_qt_enabled) then
      call set_published_field_enabled_state(qt_fields, "us_qt_local", .true.) 
      allocate(us_qt(column_size))
    end if
    if (u_qt_advection_enabled) then
      call set_published_field_enabled_state(qt_fields, "u_qt_advection_local", .true.) 
      allocate(u_qt_advection(column_size))
    end if
    if (u_qt_viscosity_diffusion_enabled) then
      call set_published_field_enabled_state(qt_fields, "u_qt_viscosity_diffusion_local", .true.) 
      allocate(u_qt_viscosity_diffusion(column_size))
    end if
    if (wu_qt_enabled) then
      call set_published_field_enabled_state(qt_fields, "wu_qt_local", .true.) 
      allocate(wu_qt(column_size))
    end if

    vs_qt_enabled=current_state%v%active .and. current_state%th%active
    v_qt_advection_enabled=is_component_field_available("v_advection") .and. &
         is_component_field_available("th_advection") .and. vs_qt_enabled
    v_qt_viscosity_diffusion_enabled=is_component_field_available("v_viscosity") .and. &
         is_component_field_available("th_diffusion") .and. vs_qt_enabled
    wv_qt_enabled=current_state%w%active .and. vs_qt_enabled

    some_qt_diagnostics_enabled=some_qt_diagnostics_enabled .or. vs_qt_enabled

    if (vs_qt_enabled) then
      call set_published_field_enabled_state(qt_fields, "vs_qt_local", .true.) 
      allocate(vs_qt(column_size))
    end if
    if (v_qt_advection_enabled) then
      call set_published_field_enabled_state(qt_fields, "v_qt_advection_local", .true.) 
      allocate(v_qt_advection(column_size))
    end if
    if (v_qt_viscosity_diffusion_enabled) then
      call set_published_field_enabled_state(qt_fields, "v_qt_viscosity_diffusion_local", .true.) 
      allocate(v_qt_viscosity_diffusion(column_size))
    end if
    if (wv_qt_enabled) then
      call set_published_field_enabled_state(qt_fields, "wv_qt_local", .true.) 
      allocate(wv_qt(column_size))
    end if

    w_qt_enabled=current_state%w%active .and. current_state%th%active
    ws_qt_enabled=w_qt_enabled
    w_qt_advection_enabled=is_component_field_available("w_advection") .and. &
         is_component_field_available("th_advection") .and. w_qt_enabled
    w_qt_viscosity_diffusion_enabled=is_component_field_available("w_viscosity") .and. &
         is_component_field_available("th_diffusion") .and. w_qt_enabled
    w_qt_buoyancy_enabled=current_state%th%active .and. is_component_field_available("w_buoyancy")
    ww_qt_enabled=w_qt_enabled

    some_qt_diagnostics_enabled=some_qt_diagnostics_enabled .or. w_qt_enabled

    if (w_qt_enabled) then
      call set_published_field_enabled_state(qt_fields, "w_qt_local", .true.) 
      allocate(w_qt(column_size))
    end if
    if (ws_qt_enabled) then
      call set_published_field_enabled_state(qt_fields, "ws_qt_local", .true.) 
      allocate(ws_qt(column_size))
    end if
    if (w_qt_advection_enabled) then
      call set_published_field_enabled_state(qt_fields, "w_qt_advection_local", .true.) 
      allocate(w_qt_advection(column_size))
    end if
    if (w_qt_viscosity_diffusion_enabled) then
      call set_published_field_enabled_state(qt_fields, "w_qt_viscosity_diffusion_local", .true.) 
      allocate(w_qt_viscosity_diffusion(column_size))
    end if
    if (w_qt_buoyancy_enabled) then
      call set_published_field_enabled_state(qt_fields, "w_qt_buoyancy_local", .true.) 
      allocate(w_qt_buoyancy(column_size))
    end if
    if (ww_qt_enabled) then
      call set_published_field_enabled_state(qt_fields, "ww_qt_local", .true.) 
      allocate(ww_qt(column_size))
    end if

    qt_qt_enabled=current_state%th%active
    sqt_qt_enabled=qt_qt_enabled
    qt_qt_advection_enabled=is_component_field_available("qt_advection") .and. qt_qt_enabled
    qt_qt_diffusion_enabled=is_component_field_available("qt_diffusion") .and. qt_qt_enabled
    wqt_qt_enabled=current_state%w%active .and. qt_qt_enabled

    some_qt_diagnostics_enabled=some_qt_diagnostics_enabled .or. qt_qt_enabled

    if (qt_qt_enabled) then
      call set_published_field_enabled_state(qt_fields, "qt_qt_local", .true.) 
      allocate(qt_qt(column_size))
    end if
    if (sqt_qt_enabled) then
      call set_published_field_enabled_state(qt_fields, "sqt_qt_local", .true.) 
      allocate(sqt_qt(column_size))
    end if
    if (qt_qt_advection_enabled) then
      call set_published_field_enabled_state(qt_fields, "qt_qt_advection_local", .true.) 
      allocate(qt_qt_advection(column_size))
    end if
    if (qt_qt_diffusion_enabled) then
      call set_published_field_enabled_state(qt_fields, "qt_qt_diffusion_local", .true.) 
      allocate(qt_qt_diffusion(column_size))
    end if
    if (wqt_qt_enabled) then
      call set_published_field_enabled_state(qt_fields, "wqt_qt_local", .true.) 
      allocate(wqt_qt(column_size))
    end if
  end subroutine initialise_qt_diagnostics

  !> Initialises the mse diagnostics. For now we are assuming mse is the same as theta, which needs updating with moisture
  !! information
  !! @param current_state The current model state
  subroutine initialise_mse_diagnostics(current_state)
    type(model_state_type), target, intent(inout) :: current_state

    integer :: column_size
    logical :: u_mse_enabled, us_mse_enabled, u_mse_advection_enabled, u_mse_viscosity_diffusion_enabled, &
         wu_mse_enabled, v_mse_enabled, vs_mse_enabled, v_mse_advection_enabled, &
         v_mse_viscosity_diffusion_enabled, wv_mse_enabled, w_mse_enabled, ws_mse_enabled, &
         w_mse_advection_enabled, w_mse_viscosity_diffusion_enabled, w_mse_buoyancy_enabled, ww_mse_enabled, &
         mse_mse_enabled, smse_mse_enabled, mse_mse_advection_enabled, &
         mse_mse_diffusion_enabled, wmse_mse_enabled

    column_size=current_state%local_grid%size(Z_INDEX)

    u_mse_enabled=current_state%u%active .and. current_state%th%active
    us_mse_enabled=u_mse_enabled
    u_mse_advection_enabled=is_component_field_available("u_advection") .and. &
         is_component_field_available("th_advection") .and. u_mse_enabled
    u_mse_viscosity_diffusion_enabled=is_component_field_available("u_viscosity") .and. &
         is_component_field_available("th_diffusion") .and. u_mse_enabled
    wu_mse_enabled=current_state%w%active .and. u_mse_enabled

    some_mse_diagnostics_enabled=u_mse_enabled

    if (u_mse_enabled) then
      call set_published_field_enabled_state(mse_fields, "u_mse_local", .true.) 
      allocate(u_mse(column_size))
    end if
    if (us_mse_enabled) then
      call set_published_field_enabled_state(mse_fields, "us_mse_local", .true.) 
      allocate(us_mse(column_size))
    end if
    if (u_mse_advection_enabled) then
      call set_published_field_enabled_state(mse_fields, "u_mse_advection_local", .true.) 
      allocate(u_mse_advection(column_size))
    end if
    if (u_mse_viscosity_diffusion_enabled) then
      call set_published_field_enabled_state(mse_fields, "u_mse_viscosity_diffusion_local", .true.) 
      allocate(u_mse_viscosity_diffusion(column_size))
    end if
    if (wu_mse_enabled) then
      call set_published_field_enabled_state(mse_fields, "wu_mse_local", .true.) 
      allocate(wu_mse(column_size))
    end if

    v_mse_enabled=current_state%v%active .and. current_state%th%active
    vs_mse_enabled=v_mse_enabled
    v_mse_advection_enabled=is_component_field_available("v_advection") .and. &
         is_component_field_available("th_advection") .and. v_mse_enabled
    v_mse_viscosity_diffusion_enabled=is_component_field_available("v_viscosity") .and. &
         is_component_field_available("th_diffusion") .and. v_mse_enabled
    wv_mse_enabled=current_state%w%active .and. v_mse_enabled

    some_mse_diagnostics_enabled=some_mse_diagnostics_enabled .or. v_mse_enabled

    if (v_mse_enabled) then
      call set_published_field_enabled_state(mse_fields, "v_mse_local", .true.) 
      allocate(v_mse(column_size))
    end if
    if (vs_mse_enabled) then
      call set_published_field_enabled_state(mse_fields, "vs_mse_local", .true.) 
      allocate(vs_mse(column_size))
    end if
    if (v_mse_advection_enabled) then
      call set_published_field_enabled_state(mse_fields, "v_mse_advection_local", .true.) 
      allocate(v_mse_advection(column_size))
    end if
    if (v_mse_viscosity_diffusion_enabled) then
      call set_published_field_enabled_state(mse_fields, "v_mse_viscosity_diffusion_local", .true.) 
      allocate(v_mse_viscosity_diffusion(column_size))
    end if
    if (wv_mse_enabled) then
      call set_published_field_enabled_state(mse_fields, "wv_mse_local", .true.) 
      allocate(wv_mse(column_size))
    end if

    w_mse_enabled=current_state%w%active .and. current_state%th%active
    ws_mse_enabled=w_mse_enabled
    w_mse_advection_enabled=is_component_field_available("w_advection") .and. &
         is_component_field_available("th_advection") .and. w_mse_enabled
    w_mse_viscosity_diffusion_enabled=is_component_field_available("w_viscosity") .and. &
         is_component_field_available("th_diffusion") .and. w_mse_enabled
    w_mse_buoyancy_enabled=current_state%th%active .and. is_component_field_available("w_buoyancy")
    ww_mse_enabled=w_mse_enabled

    some_mse_diagnostics_enabled=some_mse_diagnostics_enabled .or. w_mse_enabled

    if (w_mse_enabled) then
      call set_published_field_enabled_state(mse_fields, "w_mse_local", .true.) 
      allocate(w_mse(column_size))
    end if
    if (ws_mse_enabled) then
      call set_published_field_enabled_state(mse_fields, "ws_mse_local", .true.) 
      allocate(ws_mse(column_size))
    end if
    if (w_mse_advection_enabled) then
      call set_published_field_enabled_state(mse_fields, "w_mse_advection_local", .true.) 
      allocate(w_mse_advection(column_size))
    end if
    if (w_mse_viscosity_diffusion_enabled) then
      call set_published_field_enabled_state(mse_fields, "w_mse_viscosity_diffusion_local", .true.) 
      allocate(w_mse_viscosity_diffusion(column_size))
    end if
    if (w_mse_buoyancy_enabled) then
      call set_published_field_enabled_state(mse_fields, "w_mse_buoyancy_local", .true.) 
      allocate(w_mse_buoyancy(column_size))
    end if
    if (ww_mse_enabled) then
      call set_published_field_enabled_state(mse_fields, "ww_mse_local", .true.) 
      allocate(ww_mse(column_size))
    end if

    mse_mse_enabled=current_state%th%active
    smse_mse_enabled=mse_mse_enabled
    mse_mse_advection_enabled=is_component_field_available("mse_advection") .and. mse_mse_enabled
    mse_mse_diffusion_enabled=is_component_field_available("mse_diffusion") .and. mse_mse_enabled
    wmse_mse_enabled=current_state%w%active .and. mse_mse_enabled

    some_mse_diagnostics_enabled=some_mse_diagnostics_enabled .or. mse_mse_enabled

    if (mse_mse_enabled) then
      call set_published_field_enabled_state(mse_fields, "mse_mse_local", .true.) 
      allocate(mse_mse(column_size))
    end if
    if (smse_mse_enabled) then
      call set_published_field_enabled_state(mse_fields, "smse_mse_local", .true.) 
      allocate(smse_mse(column_size))
    end if
    if (mse_mse_advection_enabled) then
      call set_published_field_enabled_state(mse_fields, "mse_mse_advection_local", .true.) 
      allocate(mse_mse_advection(column_size))
    end if
    if (mse_mse_diffusion_enabled) then
      call set_published_field_enabled_state(mse_fields, "mse_mse_diffusion_local", .true.) 
      allocate(mse_mse_diffusion(column_size))
    end if
    if (wmse_mse_enabled) then
      call set_published_field_enabled_state(mse_fields, "wmse_mse_local", .true.) 
      allocate(wmse_mse(column_size))
    end if
  end subroutine initialise_mse_diagnostics

  !> Initialises the thetal diagnostics
  !! @param current_state The current model state
  subroutine initialise_thetal_diagnostics(current_state)
    type(model_state_type), target, intent(inout) :: current_state

    integer :: column_size
    logical :: u_thetal_enabled, us_thetal_enabled, u_thetal_advection_enabled, u_thetal_viscosity_diffusion_enabled, &
         wu_thetal_enabled, v_thetal_enabled, vs_thetal_enabled, v_thetal_advection_enabled, &
         v_thetal_viscosity_diffusion_enabled, wv_thetal_enabled, w_thetal_enabled, ws_thetal_enabled, &
         w_thetal_advection_enabled, w_thetal_viscosity_diffusion_enabled, w_thetal_buoyancy_enabled, ww_thetal_enabled, &
         thetal_thetal_enabled, sthetal_thetal_enabled, thetal_thetal_advection_enabled, &
         thetal_thetal_diffusion_enabled, wthetal_thetal_enabled

    column_size=current_state%local_grid%size(Z_INDEX)

    u_thetal_enabled=current_state%u%active .and. current_state%th%active
    us_thetal_enabled=u_thetal_enabled
    u_thetal_advection_enabled=is_component_field_available("u_advection") .and. &
         is_component_field_available("th_advection") .and. u_thetal_enabled
    u_thetal_viscosity_diffusion_enabled=is_component_field_available("u_viscosity") .and. &
         is_component_field_available("th_diffusion") .and. u_thetal_enabled
    wu_thetal_enabled=current_state%w%active .and. u_thetal_enabled

    some_thetal_diagnostics_enabled=u_thetal_enabled

    if (u_thetal_enabled) then
      call set_published_field_enabled_state(thetal_fields, "u_thetal_local", .true.) 
      allocate(u_thetal(column_size))
    end if
    if (us_thetal_enabled) then
      call set_published_field_enabled_state(thetal_fields, "us_thetal_local", .true.) 
      allocate(us_thetal(column_size))
    end if
    if (u_thetal_advection_enabled) then
      call set_published_field_enabled_state(thetal_fields, "u_thetal_advection_local", .true.) 
      allocate(u_thetal_advection(column_size))
    end if
    if (u_thetal_viscosity_diffusion_enabled) then
      call set_published_field_enabled_state(thetal_fields, "u_thetal_viscosity_diffusion_local", .true.) 
      allocate(u_thetal_viscosity_diffusion(column_size))
    end if
    if (wu_thetal_enabled) then
      call set_published_field_enabled_state(thetal_fields, "wu_thetal_local", .true.) 
      allocate(wu_thetal(column_size))
    end if

    v_thetal_enabled=current_state%v%active .and. current_state%th%active
    vs_thetal_enabled=v_thetal_enabled
    v_thetal_advection_enabled=is_component_field_available("v_advection") .and. &
         is_component_field_available("th_advection") .and. v_thetal_enabled
    v_thetal_viscosity_diffusion_enabled=is_component_field_available("v_viscosity") .and. &
         is_component_field_available("th_diffusion") .and. v_thetal_enabled
    wv_thetal_enabled=current_state%w%active .and. v_thetal_enabled

    some_thetal_diagnostics_enabled=some_thetal_diagnostics_enabled .or. v_thetal_enabled

    if (v_thetal_enabled) then
      call set_published_field_enabled_state(thetal_fields, "v_thetal_local", .true.) 
      allocate(v_thetal(column_size))
    end if
    if (vs_thetal_enabled) then
      call set_published_field_enabled_state(thetal_fields, "vs_thetal_local", .true.) 
      allocate(vs_thetal(column_size))
    end if
    if (v_thetal_advection_enabled) then
      call set_published_field_enabled_state(thetal_fields, "v_thetal_advection_local", .true.) 
      allocate(v_thetal_advection(column_size))
    end if
    if (v_thetal_viscosity_diffusion_enabled) then
      call set_published_field_enabled_state(thetal_fields, "v_thetal_viscosity_diffusion_local", .true.) 
      allocate(v_thetal_viscosity_diffusion(column_size))
    end if
    if (wv_thetal_enabled) then
      call set_published_field_enabled_state(thetal_fields, "wv_thetal_local", .true.) 
      allocate(wv_thetal(column_size))
    end if

    w_thetal_enabled=current_state%w%active .and. current_state%th%active
    ws_thetal_enabled=w_thetal_enabled
    w_thetal_advection_enabled=is_component_field_available("w_advection") .and. &
         is_component_field_available("th_advection") .and. w_thetal_enabled
    w_thetal_viscosity_diffusion_enabled=is_component_field_available("w_viscosity") .and. &
         is_component_field_available("th_diffusion") .and. w_thetal_enabled
    w_thetal_buoyancy_enabled=current_state%th%active .and. is_component_field_available("w_buoyancy")
    ww_thetal_enabled=w_thetal_enabled

    some_thetal_diagnostics_enabled=some_thetal_diagnostics_enabled .or. w_thetal_enabled

    if (w_thetal_enabled) then
      call set_published_field_enabled_state(thetal_fields, "w_thetal_local", .true.) 
      allocate(w_thetal(column_size))
    end if
    if (ws_thetal_enabled) then
      call set_published_field_enabled_state(thetal_fields, "ws_thetal_local", .true.) 
      allocate(ws_thetal(column_size))
    end if
    if (w_thetal_advection_enabled) then
      call set_published_field_enabled_state(thetal_fields, "w_thetal_advection_local", .true.) 
      allocate(w_thetal_advection(column_size))
    end if
    if (w_thetal_viscosity_diffusion_enabled) then
      call set_published_field_enabled_state(thetal_fields, "w_thetal_viscosity_diffusion_local", .true.) 
      allocate(w_thetal_viscosity_diffusion(column_size))
    end if
    if (w_thetal_buoyancy_enabled) then
      call set_published_field_enabled_state(thetal_fields, "w_thetal_buoyancy_local", .true.) 
      allocate(w_thetal_buoyancy(column_size))
    end if
    if (ww_thetal_enabled) then
      call set_published_field_enabled_state(thetal_fields, "ww_thetal_local", .true.) 
      allocate(ww_thetal(column_size))
    end if

    thetal_thetal_enabled=current_state%th%active
    sthetal_thetal_enabled=thetal_thetal_enabled
    thetal_thetal_advection_enabled=is_component_field_available("th_advection") .and. thetal_thetal_enabled
    thetal_thetal_diffusion_enabled=is_component_field_available("th_diffusion") .and. thetal_thetal_enabled
    wthetal_thetal_enabled=current_state%w%active .and. thetal_thetal_enabled

    some_thetal_diagnostics_enabled=some_thetal_diagnostics_enabled .or. thetal_thetal_enabled

    if (thetal_thetal_enabled) then
      call set_published_field_enabled_state(thetal_fields, "thetal_thetal_local", .true.) 
      allocate(thetal_thetal(column_size))
    end if
    if (sthetal_thetal_enabled) then
      call set_published_field_enabled_state(thetal_fields, "sthetal_thetal_local", .true.) 
      allocate(sthetal_thetal(column_size))
    end if
    if (thetal_thetal_advection_enabled) then
      call set_published_field_enabled_state(thetal_fields, "thetal_thetal_advection_local", .true.) 
      allocate(thetal_thetal_advection(column_size))
    end if
    if (thetal_thetal_diffusion_enabled) then
      call set_published_field_enabled_state(thetal_fields, "thetal_thetal_diffusion_local", .true.) 
      allocate(thetal_thetal_diffusion(column_size))
    end if
    if (wthetal_thetal_enabled) then
      call set_published_field_enabled_state(thetal_fields, "wthetal_thetal_local", .true.) 
      allocate(wthetal_thetal(column_size))
    end if
  end subroutine initialise_thetal_diagnostics  

  !> Initialises the prognostic (uu, vv, ww) budget diagnostics
  !! @param current_state The current model state
  subroutine initialise_prognostic_budget_diagnostics(current_state)
    type(model_state_type), target, intent(inout) :: current_state

    integer :: column_size
    logical :: tu_su_enabled, uu_advection_enabled, uu_viscosity_enabled, wu_u_enabled, tv_sv_enabled, vv_advection_enabled, &
         vv_viscosity_enabled, wv_v_enabled, tw_sw_enabled, ww_advection_enabled, ww_viscosity_enabled, ww_buoyancy_enabled

    tu_su_enabled=current_state%u%active
    uu_advection_enabled=is_component_field_available("u_advection") .and. current_state%u%active
    uu_viscosity_enabled=is_component_field_available("u_viscosity") .and. current_state%u%active
    wu_u_enabled=current_state%u%active .and. current_state%w%active
    tv_sv_enabled=current_state%v%active
    vv_advection_enabled=is_component_field_available("v_advection") .and. current_state%v%active
    vv_viscosity_enabled=is_component_field_available("v_viscosity") .and. current_state%v%active
    wv_v_enabled=current_state%v%active .and. current_state%w%active
    tw_sw_enabled=current_state%w%active
    ww_advection_enabled=is_component_field_available("w_advection") .and. current_state%w%active
    ww_viscosity_enabled=is_component_field_available("w_viscosity") .and. current_state%w%active
    ww_buoyancy_enabled=is_component_field_available("w_buoyancy") .and. current_state%w%active

    some_prognostic_budget_diagnostics_enabled=tu_su_enabled .or. tv_sv_enabled .or. tw_sw_enabled

    column_size=current_state%local_grid%size(Z_INDEX)
    
    if (tu_su_enabled) then
      call set_published_field_enabled_state(prognostic_budget_fields, "tu_su_local", .true.) 
      allocate(tu_su(column_size))
    end if
    if (uu_advection_enabled) then
      call set_published_field_enabled_state(prognostic_budget_fields, "uu_advection_local", .true.) 
      allocate(uu_advection(column_size))
    end if
    if (uu_viscosity_enabled) then
      call set_published_field_enabled_state(prognostic_budget_fields, "uu_viscosity_local", .true.) 
      allocate(uu_viscosity(column_size))
    end if
    if (wu_u_enabled) then
      call set_published_field_enabled_state(prognostic_budget_fields, "wu_u_local", .true.) 
      allocate(wu_u(column_size))
    end if
    if (tv_sv_enabled) then
      call set_published_field_enabled_state(prognostic_budget_fields, "tv_sv_local", .true.) 
      allocate(tv_sv(column_size))
    end if
    if (vv_advection_enabled) then
      call set_published_field_enabled_state(prognostic_budget_fields, "vv_advection_local", .true.) 
      allocate(vv_advection(column_size))
    end if
    if (vv_viscosity_enabled) then
      call set_published_field_enabled_state(prognostic_budget_fields, "vv_viscosity_local", .true.) 
      allocate(vv_viscosity(column_size))
    end if
    if (wv_v_enabled) then
      call set_published_field_enabled_state(prognostic_budget_fields, "wv_v_local", .true.) 
      allocate(wv_v(column_size))
    end if
    if (tw_sw_enabled) then
      call set_published_field_enabled_state(prognostic_budget_fields, "tw_sw_local", .true.) 
      allocate(tw_sw(column_size))
    end if
    if (ww_advection_enabled) then
      call set_published_field_enabled_state(prognostic_budget_fields, "ww_advection_local", .true.) 
      allocate(ww_advection(column_size))
    end if
    if (ww_viscosity_enabled) then
      call set_published_field_enabled_state(prognostic_budget_fields, "ww_viscosity_local", .true.) 
      allocate(ww_viscosity(column_size))
    end if
    if (ww_buoyancy_enabled) then
      call set_published_field_enabled_state(prognostic_budget_fields, "ww_buoyancy_local", .true.) 
      allocate(ww_buoyancy(column_size))
    end if    
  end subroutine initialise_prognostic_budget_diagnostics  

  !> Initialises the UW and VW diagnostics
  !! @param current_state Current model state
  subroutine initialise_uw_vw_diagnostics(current_state)
    type(model_state_type), target, intent(inout) :: current_state

    integer :: column_size
    logical :: uw_advection_term_enabled, vw_advection_term_enabled, uw_viscosity_term_enabled, &
       vw_viscosity_term_enabled, uw_buoyancy_term_enabled, vw_buoyancy_term_enabled, uw_tendency_term_enabled, &
       vw_tendency_term_enabled, uw_w_term_enabled, vw_w_term_enabled    

    uw_advection_term_enabled=is_component_field_available("w_advection") .and. is_component_field_available("u_advection") &
         .and. current_state%w%active
    vw_advection_term_enabled=is_component_field_available("w_advection") .and. is_component_field_available("v_advection") &
         .and. current_state%w%active
    uw_viscosity_term_enabled=is_component_field_available("w_viscosity") .and. is_component_field_available("u_viscosity") &
         .and. current_state%w%active
    vw_viscosity_term_enabled=is_component_field_available("w_viscosity") .and. is_component_field_available("v_viscosity") &
         .and. current_state%w%active
    uw_buoyancy_term_enabled=is_component_field_available("w_buoyancy")
    vw_buoyancy_term_enabled=is_component_field_available("w_buoyancy")
    uw_tendency_term_enabled=current_state%w%active .and. current_state%u%active
    vw_tendency_term_enabled=current_state%w%active .and. current_state%v%active
    uw_w_term_enabled=current_state%w%active .and. current_state%u%active
    vw_w_term_enabled=current_state%w%active .and. current_state%v%active

    some_uw_vw_diagnostics_enabled=uw_buoyancy_term_enabled .or. vw_buoyancy_term_enabled .or. uw_w_term_enabled &
      .or. vw_w_term_enabled .or. uw_advection_term_enabled .or. vw_advection_term_enabled .or. uw_viscosity_term_enabled .or. &
      vw_viscosity_term_enabled

    column_size=current_state%local_grid%size(Z_INDEX)

    if (uw_advection_term_enabled) then
      call set_published_field_enabled_state(uw_vw_fields, "uw_advection_local", .true.) 
      allocate(uw_advection(column_size))
    end if
    if (vw_advection_term_enabled) then
      call set_published_field_enabled_state(uw_vw_fields, "vw_advection_local", .true.)
      allocate(vw_advection(column_size))
    end if
    if (uw_viscosity_term_enabled) then
      call set_published_field_enabled_state(uw_vw_fields, "uw_viscosity_local", .true.)
      allocate(uw_viscosity(column_size))
    end if
    if (vw_viscosity_term_enabled) then
      call set_published_field_enabled_state(uw_vw_fields, "vw_viscosity_local", .true.)
      allocate(vw_viscosity(column_size))
    end if
    if (uw_buoyancy_term_enabled) then
      call set_published_field_enabled_state(uw_vw_fields, "uw_buoyancy_local", .true.)
      allocate(uw_buoyancy(column_size))
    end if
    if (vw_buoyancy_term_enabled) then
      call set_published_field_enabled_state(uw_vw_fields, "vw_buoyancy_local", .true.)
      allocate(vw_buoyancy(column_size))
    end if
    if (uw_tendency_term_enabled) then
      call set_published_field_enabled_state(uw_vw_fields, "uw_tendency_local", .true.)
      allocate(uw_tendency(column_size))
    end if
    if (vw_tendency_term_enabled) then
      call set_published_field_enabled_state(uw_vw_fields, "vw_tendency_local", .true.)
      allocate(vw_tendency(column_size))
    end if
    if (uw_w_term_enabled) then
      call set_published_field_enabled_state(uw_vw_fields, "uw_w_local", .true.)
      allocate(uw_w(column_size))
    end if
    if (vw_w_term_enabled) then
      call set_published_field_enabled_state(uw_vw_fields, "vw_w_local", .true.)
      allocate(vw_w(column_size))
    end if
  end subroutine initialise_uw_vw_diagnostics 

  !> Initialises the TKE diagnostics
  !! @param current_state The current model state
  subroutine initialise_TKE_diagnostics(current_state)
    type(model_state_type), target, intent(inout) :: current_state
    integer :: column_size
    logical :: sres_enabled, wp_enabled, wke_enabled, buoy_enabled, tend_enabled 

    column_size=current_state%local_grid%size(Z_INDEX)

    sres_enabled=current_state%u%active .and. current_state%v%active
    wke_enabled=current_state%w%active
    buoy_enabled=current_state%w%active
    wp_enabled=current_state%w%active
    tend_enabled=current_state%w%active

    some_tke_diagnostics_enabled=wp_enabled .or. sres_enabled .or. wke_enabled .or. buoy_enabled .or. tend_enabled

    if (wp_enabled) then
      call set_published_field_enabled_state(tke_fields, "resolved_pressure_transport_local", .true.) 
      allocate(wp(column_size))
    end if
    if (tend_enabled) then
      call set_published_field_enabled_state(tke_fields, "tke_tendency_local", .true.) 
      allocate(tend(column_size))
    end if
    if (sres_enabled) then
      call set_published_field_enabled_state(tke_fields, "resolved_shear_production_local", .true.) 
      allocate(sres(column_size))
    end if
    if (wke_enabled) then
      call set_published_field_enabled_state(tke_fields, "resolved_turbulent_transport_local", .true.) 
      allocate(wke(column_size))
    end if
    if (buoy_enabled) then
      call set_published_field_enabled_state(tke_fields, "resolved_buoyant_production_local", .true.) 
      allocate(buoy(column_size))
    end if
    
  end subroutine initialise_TKE_diagnostics   

  !> Initialises the Q field flux diagnostic areas and enabled flags depending upon the configuration of the model
  !! @param current_state The current model state
  subroutine initialise_q_flux_diagnostics(current_state)
    type(model_state_type), target, intent(inout) :: current_state

    logical :: q_flux_term_enabled, q_tendency_term_enabled, q_gradient_term_enabled, q_diff_enabled, &
       q_buoyancy_enabled

    q_flux_term_enabled=current_state%number_q_fields .gt. 0 .and. current_state%w%active
    q_tendency_term_enabled=current_state%number_q_fields .gt. 0 .and. current_state%w%active
    q_gradient_term_enabled=is_component_field_available("w_advection") .and. is_component_field_available("q_advection") &
         .and. current_state%w%active .and. current_state%number_q_fields .gt. 0
    q_diff_enabled=is_component_field_available("q_diffusion") .and. is_component_field_available("w_viscosity") &
         .and. current_state%w%active .and. current_state%number_q_fields .gt. 0
    q_buoyancy_enabled=is_component_field_available("w_buoyancy") .and. current_state%number_q_fields .gt. 0

    some_q_flux_diagnostics_enabled=q_flux_term_enabled .or. q_buoyancy_enabled
                    
    if (q_flux_term_enabled) then
      call set_published_field_enabled_state(q_flux_fields, "q_flux_transport_local", .true.)
      allocate(q_flux_values(current_state%local_grid%size(Z_INDEX), current_state%number_q_fields))
    end if
    if (q_tendency_term_enabled) then
      call set_published_field_enabled_state(q_flux_fields, "q_flux_tendency_local", .true.)
      allocate(q_tendency(current_state%local_grid%size(Z_INDEX), current_state%number_q_fields))
    end if
    if (q_gradient_term_enabled) then
      call set_published_field_enabled_state(q_flux_fields, "q_flux_gradient_local", .true.)
      allocate(q_gradient(current_state%local_grid%size(Z_INDEX), current_state%number_q_fields))
    end if
    if (q_diff_enabled) then
      call set_published_field_enabled_state(q_flux_fields, "q_flux_dissipation_local", .true.)
      allocate(q_diff(current_state%local_grid%size(Z_INDEX), current_state%number_q_fields))
    end if
    if (q_buoyancy_enabled) then
      call set_published_field_enabled_state(q_flux_fields, "q_flux_buoyancy_local", .true.)
      allocate(q_buoyancy(current_state%local_grid%size(Z_INDEX), current_state%number_q_fields))
    end if
  end subroutine initialise_q_flux_diagnostics

  !> Initialises the heat flux diagnostic areas and enabled flags depending upon the configuration of the model
  !! @param current_state The current model state
  subroutine initialise_theta_flux_diagnostics(current_state)
    type(model_state_type), target, intent(inout) :: current_state

    logical :: th_flux_term_enabled, th_tendency_term_enabled, th_diff_enabled, th_gradient_term_enabled, th_buoyancy_enabled

    th_flux_term_enabled=current_state%th%active .and. current_state%w%active
    th_tendency_term_enabled=current_state%th%active .and. current_state%w%active
    th_gradient_term_enabled=is_component_field_available("w_advection") .and. is_component_field_available("th_advection") &
         .and. current_state%w%active .and. current_state%th%active
    th_diff_enabled=is_component_field_available("th_diffusion") .and. is_component_field_available("w_viscosity") &
         .and. current_state%w%active .and. current_state%th%active
    th_buoyancy_enabled=is_component_field_available("w_buoyancy") .and. current_state%th%active

    some_theta_flux_diagnostics_enabled=th_flux_term_enabled .or. th_buoyancy_enabled

    if (th_flux_term_enabled) then
      call set_published_field_enabled_state(heat_flux_fields, "heat_flux_transport_local", .true.)
      allocate(th_flux_values(current_state%local_grid%size(Z_INDEX)))
    end if
    if (th_tendency_term_enabled) then
      call set_published_field_enabled_state(heat_flux_fields, "heat_flux_tendency_local", .true.)
      allocate(th_tendency(current_state%local_grid%size(Z_INDEX)))      
    end if
    if (th_diff_enabled) then
      call set_published_field_enabled_state(heat_flux_fields, "heat_flux_dissipation_local", .true.)
      allocate(th_diff(current_state%local_grid%size(Z_INDEX)))
    end if
    if (th_gradient_term_enabled) then
      call set_published_field_enabled_state(heat_flux_fields, "heat_flux_gradient_local", .true.)
      allocate(th_gradient(current_state%local_grid%size(Z_INDEX)))
    end if
    if (th_buoyancy_enabled) then
      call set_published_field_enabled_state(heat_flux_fields, "heat_flux_buoyancy_local", .true.)
      allocate(th_buoyancy(current_state%local_grid%size(Z_INDEX)))
    end if
  end subroutine initialise_theta_flux_diagnostics

  !> Clears the qt diagnostics
  subroutine clear_qt()
    if (allocated(us_qt)) us_qt=0.0_DEFAULT_PRECISION
    if (allocated(u_qt_advection)) u_qt_advection=0.0_DEFAULT_PRECISION
    if (allocated(u_qt_viscosity_diffusion)) u_qt_viscosity_diffusion=0.0_DEFAULT_PRECISION
    if (allocated(wu_qt)) wu_qt=0.0_DEFAULT_PRECISION
    if (allocated(vs_qt)) vs_qt=0.0_DEFAULT_PRECISION
    if (allocated(v_qt_advection)) v_qt_advection=0.0_DEFAULT_PRECISION
    if (allocated(v_qt_viscosity_diffusion)) v_qt_viscosity_diffusion=0.0_DEFAULT_PRECISION
    if (allocated(wv_qt)) wv_qt=0.0_DEFAULT_PRECISION
    if (allocated(w_qt)) w_qt=0.0_DEFAULT_PRECISION
    if (allocated(ws_qt)) ws_qt=0.0_DEFAULT_PRECISION
    if (allocated(w_qt_advection)) w_qt_advection=0.0_DEFAULT_PRECISION
    if (allocated(w_qt_viscosity_diffusion)) w_qt_viscosity_diffusion=0.0_DEFAULT_PRECISION
    if (allocated(w_qt_buoyancy)) w_qt_buoyancy=0.0_DEFAULT_PRECISION
    if (allocated(ww_qt)) ww_qt=0.0_DEFAULT_PRECISION
    if (allocated(qt_qt)) qt_qt=0.0_DEFAULT_PRECISION
    if (allocated(sqt_qt)) sqt_qt=0.0_DEFAULT_PRECISION
    if (allocated(qt_qt_advection)) qt_qt_advection=0.0_DEFAULT_PRECISION
    if (allocated(qt_qt_diffusion)) qt_qt_diffusion=0.0_DEFAULT_PRECISION
    if (allocated(wqt_qt)) wqt_qt=0.0_DEFAULT_PRECISION
  end subroutine clear_qt    

  !> Computes the qt diagnostics for a specific column. For now we are assuming qt is the same as theta, 
  !! which needs updating with moisture information as per resdgs
  !! @param current_state The current model state
  subroutine compute_qt_for_column(current_state)
    type(model_state_type), target, intent(inout) :: current_state

    real(kind=DEFAULT_PRECISION), dimension(current_state%local_grid%size(Z_INDEX)) :: upr, vpr, uprm1, vprm1, qtpr, qtprp1
    type(component_field_value_type) :: u_advection, u_viscosity, th_advection, th_diffusion, v_advection, v_viscosity, &
         w_advection, w_viscosity, w_buoyancy
    integer :: k

    if (is_component_field_available("u_advection")) u_advection=get_component_field_value(current_state, "u_advection")
    if (is_component_field_available("u_viscosity")) u_viscosity=get_component_field_value(current_state, "u_viscosity")
    if (is_component_field_available("th_advection")) th_advection=get_component_field_value(current_state, "th_advection")
    if (is_component_field_available("th_diffusion")) th_diffusion=get_component_field_value(current_state, "th_diffusion")
    if (is_component_field_available("v_advection")) v_advection=get_component_field_value(current_state, "v_advection")
    if (is_component_field_available("v_viscosity")) v_viscosity=get_component_field_value(current_state, "v_viscosity")
    if (is_component_field_available("w_advection")) w_advection=get_component_field_value(current_state, "w_advection")
    if (is_component_field_available("w_viscosity")) w_viscosity=get_component_field_value(current_state, "w_viscosity")
    if (is_component_field_available("w_buoyancy")) w_buoyancy=get_component_field_value(current_state, "w_buoyancy")

    do k=1, current_state%local_grid%size(Z_INDEX)
      upr(k)=current_state%u%data(k,current_state%column_local_y,current_state%column_local_x)
      uprm1(k)=current_state%u%data(k,current_state%column_local_y,current_state%column_local_x-1)
      if (allocated(current_state%global_grid%configuration%vertical%olubar)) then
        upr(k)=upr(k)-(current_state%global_grid%configuration%vertical%olubar(k)-current_state%ugal)
        uprm1(k)=uprm1(k)-(current_state%global_grid%configuration%vertical%olubar(k)-current_state%ugal)
      end if
      vpr(k)=current_state%v%data(k,current_state%column_local_y,current_state%column_local_x)
      vprm1(k)=current_state%v%data(k,current_state%column_local_y-1,current_state%column_local_x)
      if (allocated(current_state%global_grid%configuration%vertical%olvbar)) then
        vpr(k)=vpr(k)-(current_state%global_grid%configuration%vertical%olvbar(k)-current_state%vgal)
        vprm1(k)=vprm1(k)-(current_state%global_grid%configuration%vertical%olvbar(k)-current_state%vgal)
      end if

      qtpr(k)=current_state%th%data(k,current_state%column_local_y,current_state%column_local_x)
      qtprp1(k)=current_state%th%data(k,current_state%column_local_y,current_state%column_local_x+1)
      if (allocated(current_state%global_grid%configuration%vertical%olthbar)) then
        qtpr(k)=qtpr(k)-current_state%global_grid%configuration%vertical%olthbar(k)
        qtprp1(k)=qtprp1(k)-current_state%global_grid%configuration%vertical%olthbar(k)
      end if
    end do
    do k=2, current_state%local_grid%size(Z_INDEX)-1
      if (allocated(us_qt)) us_qt(k)=us_qt(k)+0.5*(upr(k)+uprm1(k))*&
           current_state%sth%data(k,current_state%column_local_y,current_state%column_local_x)+0.5*(qtpr(k)+qtprp1(k))*&
           current_state%su%data(k,current_state%column_local_y,current_state%column_local_x)
      if (allocated(u_qt_advection)) u_qt_advection(k)=u_qt_advection(k)+0.5*(upr(k)+uprm1(k))*&
           th_advection%real_1d_array(k)+0.5*(qtpr(k)+qtprp1(k))*u_advection%real_1d_array(k)
      if (allocated(u_qt_viscosity_diffusion)) u_qt_viscosity_diffusion(k)=u_qt_viscosity_diffusion(k)+0.5*&
           (qtpr(k)+qtprp1(k))*u_viscosity%real_1d_array(k)+0.5*(upr(k)+uprm1(k))*th_diffusion%real_1d_array(k)
      if (allocated(wu_qt)) wu_qt(k)=wu_qt(k)+&
           current_state%w%data(k,current_state%column_local_y,current_state%column_local_x)*0.25*(upr(k+1)+upr(k)+&
           uprm1(k+1)+uprm1(k))*0.5*(qtpr(k+1)+qtpr(k))

      if (allocated(vs_qt)) vs_qt(k)=vs_qt(k)+0.5*(vpr(k)+vprm1(k))*&
           current_state%sth%data(k,current_state%column_local_y,current_state%column_local_x)+0.5*(qtpr(k)+qtprp1(k))*&
           current_state%sv%data(k,current_state%column_local_y,current_state%column_local_x)
      if (allocated(v_qt_advection)) v_qt_advection(k)=v_qt_advection(k)+0.5*(vpr(k)+vprm1(k))*&
           th_advection%real_1d_array(k)+0.5*(qtpr(k)+qtprp1(k))*v_advection%real_1d_array(k)
      if (allocated(v_qt_viscosity_diffusion)) v_qt_viscosity_diffusion(k)=v_qt_viscosity_diffusion(k)+0.5*&
           (qtpr(k)+qtprp1(k))*v_viscosity%real_1d_array(k)+0.5*(vpr(k)+vprm1(k))*th_diffusion%real_1d_array(k)
      if (allocated(wv_qt)) wv_qt(k)=wv_qt(k)+&
           current_state%w%data(k,current_state%column_local_y,current_state%column_local_x)*0.25*(vpr(k+1)+vpr(k)+&
           vprm1(k+1)+vprm1(k))*0.5*(qtpr(k+1)+qtpr(k))

      if (allocated(w_qt)) w_qt(k)=w_qt(k)+&
           current_state%w%data(k,current_state%column_local_y,current_state%column_local_x)*0.5*(qtpr(k)+qtpr(k+1))
      if (allocated(ws_qt)) ws_qt(k)=ws_qt(k)+&
           current_state%w%data(k,current_state%column_local_y,current_state%column_local_x)*0.5*&
           (current_state%sth%data(k,current_state%column_local_y,current_state%column_local_x)+&
           current_state%sth%data(k+1,current_state%column_local_y,current_state%column_local_x))+0.5*(qtpr(k)+qtpr(k+1))*&
           current_state%sw%data(k,current_state%column_local_y,current_state%column_local_x)
      if (allocated(w_qt_advection)) w_qt_advection(k)=w_qt_advection(k)+&
           current_state%w%data(k,current_state%column_local_y,current_state%column_local_x)*0.5*(th_advection%real_1d_array(k)+&
           th_advection%real_1d_array(k+1))+0.5*(qtpr(k)+qtpr(k+1))*w_advection%real_1d_array(k)
      if (allocated(w_qt_viscosity_diffusion)) w_qt_viscosity_diffusion(k)=w_qt_viscosity_diffusion(k)+0.5*&
           (qtpr(k)+qtpr(k+1))*w_viscosity%real_1d_array(k)+&
           current_state%w%data(k,current_state%column_local_y,current_state%column_local_x)*0.5*&
           (th_diffusion%real_1d_array(k)+th_diffusion%real_1d_array(k+1))
      if (allocated(w_qt_buoyancy)) w_qt_buoyancy(k)=w_qt_buoyancy(k)+0.5*(qtpr(k+1)+qtpr(k))*&
           w_buoyancy%real_1d_array(k)
      if (allocated(wv_qt)) ww_qt(k)=ww_qt(k)+0.5*&
           (current_state%w%data(k,current_state%column_local_y,current_state%column_local_x)+&
           current_state%w%data(k-1,current_state%column_local_y,current_state%column_local_x))*0.5*(&
           current_state%w%data(k,current_state%column_local_y,current_state%column_local_x)+&
           current_state%w%data(k-1,current_state%column_local_y,current_state%column_local_x))*qtpr(k)

      if (allocated(qt_qt)) qt_qt(k)=qt_qt(k)+qtpr(k)*qtpr(k)
      if (allocated(sqt_qt)) sqt_qt(k)=sqt_qt(k)+2.0*qtpr(k)*&
           current_state%sth%data(k,current_state%column_local_y,current_state%column_local_x)
      if (allocated(qt_qt_advection)) qt_qt_advection(k)=qt_qt_advection(k)+2.0*qtpr(k)*&
           th_advection%real_1d_array(k)
      if (allocated(qt_qt_diffusion)) qt_qt_diffusion(k)=qt_qt_diffusion(k)+2.0*qtpr(k)*&
           th_diffusion%real_1d_array(k)
      if (allocated(wqt_qt)) wqt_qt(k)=wqt_qt(k)+&
           current_state%w%data(k,current_state%column_local_y,current_state%column_local_x)*0.5*(qtpr(k+1)+qtpr(k))*0.5*(&
           qtpr(k+1)+qtpr(k)) 
    end do

    if (allocated(u_advection%real_1d_array)) deallocate(u_advection%real_1d_array)
    if (allocated(u_viscosity%real_1d_array)) deallocate(u_viscosity%real_1d_array)
    if (allocated(th_advection%real_1d_array)) deallocate(th_advection%real_1d_array)
    if (allocated(th_diffusion%real_1d_array)) deallocate(th_diffusion%real_1d_array)
    if (allocated(v_advection%real_1d_array)) deallocate(v_advection%real_1d_array)
    if (allocated(v_viscosity%real_1d_array)) deallocate(v_viscosity%real_1d_array)
    if (allocated(w_advection%real_1d_array)) deallocate(w_advection%real_1d_array)
    if (allocated(w_viscosity%real_1d_array)) deallocate(w_viscosity%real_1d_array)
    if (allocated(w_buoyancy%real_1d_array)) deallocate(w_buoyancy%real_1d_array)
  end subroutine compute_qt_for_column

  !> Clears the scalar diagnostics
  subroutine clear_scalars()
    mflux=0.0_DEFAULT_PRECISION
  end subroutine clear_scalars
  
  !> Computes the scalar diagnostics for a specific column. 
  !! @param current_state The current model state
  subroutine compute_scalars_for_column(current_state)
    type(model_state_type), target, intent(inout) :: current_state

    integer :: k

    do k=2, current_state%local_grid%size(Z_INDEX)-1
      if (current_state%w%data(k, current_state%column_local_y,current_state%column_local_x) .gt. wmfcrit) then
        mflux=mflux+current_state%global_grid%configuration%vertical%rho(k)*&
             current_state%global_grid%configuration%vertical%dzn(k)*&
             current_state%w%data(k, current_state%column_local_y,current_state%column_local_x)
      end if
    end do
  end subroutine compute_scalars_for_column

  !> Clears the mse diagnostics
  subroutine clear_mse()
    if (allocated(u_mse)) u_mse=0.0_DEFAULT_PRECISION
    if (allocated(us_mse)) us_mse=0.0_DEFAULT_PRECISION
    if (allocated(u_mse_advection)) u_mse_advection=0.0_DEFAULT_PRECISION
    if (allocated(u_mse_viscosity_diffusion)) u_mse_viscosity_diffusion=0.0_DEFAULT_PRECISION
    if (allocated(wu_mse)) wu_mse=0.0_DEFAULT_PRECISION
    if (allocated(v_mse)) v_mse=0.0_DEFAULT_PRECISION
    if (allocated(vs_mse)) vs_mse=0.0_DEFAULT_PRECISION
    if (allocated(v_mse_advection)) v_mse_advection=0.0_DEFAULT_PRECISION
    if (allocated(v_mse_viscosity_diffusion)) v_mse_viscosity_diffusion=0.0_DEFAULT_PRECISION
    if (allocated(wv_mse)) wv_mse=0.0_DEFAULT_PRECISION
    if (allocated(w_mse)) w_mse=0.0_DEFAULT_PRECISION
    if (allocated(ws_mse)) ws_mse=0.0_DEFAULT_PRECISION
    if (allocated(w_mse_advection)) w_mse_advection=0.0_DEFAULT_PRECISION
    if (allocated(w_mse_viscosity_diffusion)) w_mse_viscosity_diffusion=0.0_DEFAULT_PRECISION
    if (allocated(w_mse_buoyancy)) w_mse_buoyancy=0.0_DEFAULT_PRECISION
    if (allocated(ww_mse)) ww_mse=0.0_DEFAULT_PRECISION
    if (allocated(mse_mse)) mse_mse=0.0_DEFAULT_PRECISION
    if (allocated(smse_mse)) smse_mse=0.0_DEFAULT_PRECISION
    if (allocated(mse_mse_advection)) mse_mse_advection=0.0_DEFAULT_PRECISION
    if (allocated(mse_mse_diffusion)) mse_mse_diffusion=0.0_DEFAULT_PRECISION
    if (allocated(wmse_mse)) wmse_mse=0.0_DEFAULT_PRECISION
  end subroutine clear_mse
  
  !> Computes the mse diagnostics for a specific column. For now we are assuming mse is the same as theta, 
  !! which needs updating with moisture information
  !! @param current_state The current model state
  subroutine compute_mse_for_column(current_state)
    type(model_state_type), target, intent(inout) :: current_state

    real(kind=DEFAULT_PRECISION), dimension(current_state%local_grid%size(Z_INDEX)) :: upr, vpr, uprm1, vprm1, msepr, mseprp1
    type(component_field_value_type) :: u_advection, u_viscosity, th_advection, th_diffusion, v_advection, v_viscosity, &
         w_advection, w_viscosity, w_buoyancy
    integer :: k

    if (is_component_field_available("u_advection")) u_advection=get_component_field_value(current_state, "u_advection")
    if (is_component_field_available("u_viscosity")) u_viscosity=get_component_field_value(current_state, "u_viscosity")
    if (is_component_field_available("th_advection")) th_advection=get_component_field_value(current_state, "th_advection")
    if (is_component_field_available("th_diffusion")) th_diffusion=get_component_field_value(current_state, "th_diffusion")
    if (is_component_field_available("v_advection")) v_advection=get_component_field_value(current_state, "v_advection")
    if (is_component_field_available("v_viscosity")) v_viscosity=get_component_field_value(current_state, "v_viscosity")
    if (is_component_field_available("w_advection")) w_advection=get_component_field_value(current_state, "w_advection")
    if (is_component_field_available("w_viscosity")) w_viscosity=get_component_field_value(current_state, "w_viscosity")
    if (is_component_field_available("w_buoyancy")) w_buoyancy=get_component_field_value(current_state, "w_buoyancy")

    do k=1, current_state%local_grid%size(Z_INDEX)
      upr(k)=current_state%u%data(k,current_state%column_local_y,current_state%column_local_x)
      uprm1(k)=current_state%u%data(k,current_state%column_local_y,current_state%column_local_x-1)
      if (allocated(current_state%global_grid%configuration%vertical%olubar)) then
        upr(k)=upr(k)-(current_state%global_grid%configuration%vertical%olubar(k)-current_state%ugal)
        uprm1(k)=uprm1(k)-(current_state%global_grid%configuration%vertical%olubar(k)-current_state%ugal)
      end if
      vpr(k)=current_state%v%data(k,current_state%column_local_y,current_state%column_local_x)
      vprm1(k)=current_state%v%data(k,current_state%column_local_y-1,current_state%column_local_x)
      if (allocated(current_state%global_grid%configuration%vertical%olvbar)) then
        vpr(k)=vpr(k)-(current_state%global_grid%configuration%vertical%olvbar(k)-current_state%vgal)
        vprm1(k)=vprm1(k)-(current_state%global_grid%configuration%vertical%olvbar(k)-current_state%vgal)
      end if
      msepr(k)=current_state%th%data(k,current_state%column_local_y,current_state%column_local_x)
      mseprp1(k)=current_state%th%data(k,current_state%column_local_y,current_state%column_local_x+1)
      if (allocated(current_state%global_grid%configuration%vertical%olthbar)) then
        msepr(k)=msepr(k)-current_state%global_grid%configuration%vertical%olthbar(k)
        mseprp1(k)=mseprp1(k)-current_state%global_grid%configuration%vertical%olthbar(k)
      end if
    end do
    do k=2, current_state%local_grid%size(Z_INDEX)-1
      if (allocated(u_mse)) u_mse(k)=u_mse(k)+0.5*(upr(k)+uprm1(k))*msepr(k)
      if (allocated(us_mse)) us_mse(k)=us_mse(k)+0.5*(upr(k)+uprm1(k))*&
           current_state%sth%data(k,current_state%column_local_y,current_state%column_local_x)+0.5*(msepr(k)+mseprp1(k))*&
           current_state%su%data(k,current_state%column_local_y,current_state%column_local_x)
      if (allocated(u_mse_advection)) u_mse_advection(k)=u_mse_advection(k)+0.5*(upr(k)+uprm1(k))*&
           th_advection%real_1d_array(k)+0.5*(msepr(k)+mseprp1(k))*u_advection%real_1d_array(k)
      if (allocated(u_mse_viscosity_diffusion)) u_mse_viscosity_diffusion(k)=u_mse_viscosity_diffusion(k)+0.5*&
           (msepr(k)+mseprp1(k))*u_viscosity%real_1d_array(k)+0.5*(upr(k)+uprm1(k))*th_diffusion%real_1d_array(k)
      if (allocated(wu_mse)) wu_mse(k)=wu_mse(k)+&
           current_state%w%data(k,current_state%column_local_y,current_state%column_local_x)*0.25*(upr(k+1)+upr(k)+&
           uprm1(k+1)+uprm1(k))*0.5*(msepr(k+1)+msepr(k))

      if (allocated(v_mse)) v_mse(k)=v_mse(k)+0.5*(vpr(k)+vprm1(k))*msepr(k)
      if (allocated(vs_mse)) vs_mse(k)=vs_mse(k)+0.5*(vpr(k)+vprm1(k))*&
           current_state%sth%data(k,current_state%column_local_y,current_state%column_local_x)+0.5*(msepr(k)+mseprp1(k))*&
           current_state%sv%data(k,current_state%column_local_y,current_state%column_local_x)
      if (allocated(v_mse_advection)) v_mse_advection(k)=v_mse_advection(k)+0.5*(vpr(k)+vprm1(k))*&
           th_advection%real_1d_array(k)+0.5*(msepr(k)+mseprp1(k))*v_advection%real_1d_array(k)
      if (allocated(v_mse_viscosity_diffusion)) v_mse_viscosity_diffusion(k)=v_mse_viscosity_diffusion(k)+0.5*&
           (msepr(k)+mseprp1(k))*v_viscosity%real_1d_array(k)+0.5*(vpr(k)+vprm1(k))*th_diffusion%real_1d_array(k)
      if (allocated(wv_mse)) wv_mse(k)=wv_mse(k)+&
           current_state%w%data(k,current_state%column_local_y,current_state%column_local_x)*0.25*(vpr(k+1)+vpr(k)+&
           vprm1(k+1)+vprm1(k))*0.5*(msepr(k+1)+msepr(k))

      if (allocated(w_mse)) w_mse(k)=w_mse(k)+&
           current_state%w%data(k,current_state%column_local_y,current_state%column_local_x)*0.5*(msepr(k)+msepr(k+1))
      if (allocated(ws_mse)) ws_mse(k)=ws_mse(k)+&
           current_state%w%data(k,current_state%column_local_y,current_state%column_local_x)*0.5*&
           (current_state%sth%data(k,current_state%column_local_y,current_state%column_local_x)+&
           current_state%sth%data(k+1,current_state%column_local_y,current_state%column_local_x))+0.5*(msepr(k)+msepr(k+1))*&
           current_state%sw%data(k,current_state%column_local_y,current_state%column_local_x)
      if (allocated(w_mse_advection)) w_mse_advection(k)=w_mse_advection(k)+&
           current_state%w%data(k,current_state%column_local_y,current_state%column_local_x)*0.5*(th_advection%real_1d_array(k)+&
           th_advection%real_1d_array(k+1))+0.5*(msepr(k)+msepr(k+1))*w_advection%real_1d_array(k)
      if (allocated(w_mse_viscosity_diffusion)) w_mse_viscosity_diffusion(k)=w_mse_viscosity_diffusion(k)+0.5*&
           (msepr(k)+msepr(k+1))*w_viscosity%real_1d_array(k)+&
           current_state%w%data(k,current_state%column_local_y,current_state%column_local_x)*0.5*&
           (th_diffusion%real_1d_array(k)+th_diffusion%real_1d_array(k+1))
      if (allocated(w_mse_buoyancy)) w_mse_buoyancy(k)=w_mse_buoyancy(k)+0.5*(msepr(k+1)+msepr(k))*&
           w_buoyancy%real_1d_array(k)
      if (allocated(wv_mse)) ww_mse(k)=ww_mse(k)+0.5*&
           (current_state%w%data(k,current_state%column_local_y,current_state%column_local_x)+&
           current_state%w%data(k-1,current_state%column_local_y,current_state%column_local_x))*0.5*(&
           current_state%w%data(k,current_state%column_local_y,current_state%column_local_x)+&
           current_state%w%data(k-1,current_state%column_local_y,current_state%column_local_x))*msepr(k)

      if (allocated(mse_mse)) mse_mse(k)=mse_mse(k)+msepr(k)*msepr(k)
      if (allocated(smse_mse)) smse_mse(k)=smse_mse(k)+2.0*msepr(k)*&
           current_state%sth%data(k,current_state%column_local_y,current_state%column_local_x)
      if (allocated(mse_mse_advection)) mse_mse_advection(k)=mse_mse_advection(k)+2.0*msepr(k)*&
           th_advection%real_1d_array(k)
      if (allocated(mse_mse_diffusion)) mse_mse_diffusion(k)=mse_mse_diffusion(k)+2.0*msepr(k)*&
           th_diffusion%real_1d_array(k)
      if (allocated(wmse_mse)) wmse_mse(k)=wmse_mse(k)+&
           current_state%w%data(k,current_state%column_local_y,current_state%column_local_x)*0.5*(msepr(k+1)+msepr(k))*0.5*(&
           msepr(k+1)+msepr(k)) 
    end do

    if (allocated(u_advection%real_1d_array)) deallocate(u_advection%real_1d_array)
    if (allocated(u_viscosity%real_1d_array)) deallocate(u_viscosity%real_1d_array)
    if (allocated(th_advection%real_1d_array)) deallocate(th_advection%real_1d_array)
    if (allocated(th_diffusion%real_1d_array)) deallocate(th_diffusion%real_1d_array)
    if (allocated(v_advection%real_1d_array)) deallocate(v_advection%real_1d_array)
    if (allocated(v_viscosity%real_1d_array)) deallocate(v_viscosity%real_1d_array)
    if (allocated(w_advection%real_1d_array)) deallocate(w_advection%real_1d_array)
    if (allocated(w_viscosity%real_1d_array)) deallocate(w_viscosity%real_1d_array)
    if (allocated(w_buoyancy%real_1d_array)) deallocate(w_buoyancy%real_1d_array)
  end subroutine compute_mse_for_column

  !> Clears the thetal diagnostics
  subroutine clear_thetal()
    if (allocated(u_thetal)) u_thetal=0.0_DEFAULT_PRECISION
    if (allocated(us_thetal)) us_thetal=0.0_DEFAULT_PRECISION
    if (allocated(u_thetal_advection)) u_thetal_advection=0.0_DEFAULT_PRECISION
    if (allocated(u_thetal_viscosity_diffusion)) u_thetal_viscosity_diffusion=0.0_DEFAULT_PRECISION
    if (allocated(wu_thetal)) wu_thetal=0.0_DEFAULT_PRECISION
    if (allocated(v_thetal)) v_thetal=0.0_DEFAULT_PRECISION
    if (allocated(vs_thetal)) vs_thetal=0.0_DEFAULT_PRECISION
    if (allocated(v_thetal_advection)) v_thetal_advection=0.0_DEFAULT_PRECISION
    if (allocated(v_thetal_viscosity_diffusion)) v_thetal_viscosity_diffusion=0.0_DEFAULT_PRECISION
    if (allocated(wv_thetal)) wv_thetal=0.0_DEFAULT_PRECISION
    if (allocated(w_thetal)) w_thetal=0.0_DEFAULT_PRECISION
    if (allocated(ws_thetal)) ws_thetal=0.0_DEFAULT_PRECISION
    if (allocated(w_thetal_advection)) w_thetal_advection=0.0_DEFAULT_PRECISION
    if (allocated(w_thetal_viscosity_diffusion)) w_thetal_viscosity_diffusion=0.0_DEFAULT_PRECISION
    if (allocated(w_thetal_buoyancy)) w_thetal_buoyancy=0.0_DEFAULT_PRECISION
    if (allocated(ww_thetal)) ww_thetal=0.0_DEFAULT_PRECISION
    if (allocated(thetal_thetal)) thetal_thetal=0.0_DEFAULT_PRECISION
    if (allocated(sthetal_thetal)) sthetal_thetal=0.0_DEFAULT_PRECISION
    if (allocated(thetal_thetal_advection)) thetal_thetal_advection=0.0_DEFAULT_PRECISION
    if (allocated(thetal_thetal_diffusion)) thetal_thetal_diffusion=0.0_DEFAULT_PRECISION
    if (allocated(wthetal_thetal)) wthetal_thetal=0.0_DEFAULT_PRECISION
  end subroutine clear_thetal

  !> Computes the thetal diagnostics for a specific column
  !! @param current_state The current model state
  subroutine compute_thetal_for_column(current_state)
    type(model_state_type), target, intent(inout) :: current_state

    real(kind=DEFAULT_PRECISION), dimension(current_state%local_grid%size(Z_INDEX)) :: upr, vpr, uprm1, vprm1, thlpr, thlprp1
    type(component_field_value_type) :: u_advection, u_viscosity, th_advection, th_diffusion, v_advection, v_viscosity, &
         w_advection, w_viscosity, w_buoyancy
    integer :: k
    
    if (is_component_field_available("u_advection")) u_advection=get_component_field_value(current_state, "u_advection")
    if (is_component_field_available("u_viscosity")) u_viscosity=get_component_field_value(current_state, "u_viscosity")
    if (is_component_field_available("th_advection")) th_advection=get_component_field_value(current_state, "th_advection")
    if (is_component_field_available("th_diffusion")) th_diffusion=get_component_field_value(current_state, "th_diffusion")
    if (is_component_field_available("v_advection")) v_advection=get_component_field_value(current_state, "v_advection")
    if (is_component_field_available("v_viscosity")) v_viscosity=get_component_field_value(current_state, "v_viscosity")
    if (is_component_field_available("w_advection")) w_advection=get_component_field_value(current_state, "w_advection")
    if (is_component_field_available("w_viscosity")) w_viscosity=get_component_field_value(current_state, "w_viscosity")
    if (is_component_field_available("w_buoyancy")) w_buoyancy=get_component_field_value(current_state, "w_buoyancy")

    do k=1, current_state%local_grid%size(Z_INDEX)
      upr(k)=current_state%u%data(k,current_state%column_local_y,current_state%column_local_x)
      uprm1(k)=current_state%u%data(k,current_state%column_local_y,current_state%column_local_x-1)
      if (allocated(current_state%global_grid%configuration%vertical%olubar)) then
        upr(k)=upr(k)-(current_state%global_grid%configuration%vertical%olubar(k)-current_state%ugal)
        uprm1(k)=uprm1(k)-(current_state%global_grid%configuration%vertical%olubar(k)-current_state%ugal)
      end if
      vpr(k)=current_state%v%data(k,current_state%column_local_y,current_state%column_local_x)
      vprm1(k)=current_state%v%data(k,current_state%column_local_y-1,current_state%column_local_x)
      if (allocated(current_state%global_grid%configuration%vertical%olvbar)) then
        vpr(k)=vpr(k)-(current_state%global_grid%configuration%vertical%olvbar(k)-current_state%vgal)
        vprm1(k)=vprm1(k)-(current_state%global_grid%configuration%vertical%olvbar(k)-current_state%vgal)
      end if
      thlpr(k)=current_state%th%data(k,current_state%column_local_y,current_state%column_local_x)
      thlprp1(k)=current_state%th%data(k,current_state%column_local_y,current_state%column_local_x+1)
      if (allocated(current_state%global_grid%configuration%vertical%olthbar)) then
        thlpr(k)=thlpr(k)-current_state%global_grid%configuration%vertical%olthbar(k)
        thlprp1(k)=thlprp1(k)-current_state%global_grid%configuration%vertical%olthbar(k)
      end if      
    end do
    do k=2, current_state%local_grid%size(Z_INDEX)-1
      if (allocated(u_thetal)) u_thetal(k)=u_thetal(k)+0.5*(upr(k)+uprm1(k))*thlpr(k)
      if (allocated(us_thetal)) us_thetal(k)=us_thetal(k)+0.5*(upr(k)+uprm1(k))*&
           current_state%sth%data(k,current_state%column_local_y,current_state%column_local_x)+0.5*(thlpr(k)+thlprp1(k))*&
           current_state%su%data(k,current_state%column_local_y,current_state%column_local_x)
      if (allocated(u_thetal_advection)) u_thetal_advection(k)=u_thetal_advection(k)+0.5*(upr(k)+uprm1(k))*&
           th_advection%real_1d_array(k)+0.5*(thlpr(k)+thlprp1(k))*u_advection%real_1d_array(k)
      if (allocated(u_thetal_viscosity_diffusion)) u_thetal_viscosity_diffusion(k)=u_thetal_viscosity_diffusion(k)+0.5*&
           (thlpr(k)+thlprp1(k))*u_viscosity%real_1d_array(k)+0.5*(upr(k)+uprm1(k))*th_diffusion%real_1d_array(k)
      if (allocated(wu_thetal)) wu_thetal(k)=wu_thetal(k)+&
           current_state%w%data(k,current_state%column_local_y,current_state%column_local_x)*0.25*(upr(k+1)+upr(k)+&
           uprm1(k+1)+uprm1(k))*0.5*(thlpr(k+1)+thlpr(k))

      if (allocated(v_thetal)) v_thetal(k)=v_thetal(k)+0.5*(vpr(k)+vprm1(k))*thlpr(k)
      if (allocated(vs_thetal)) vs_thetal(k)=vs_thetal(k)+0.5*(vpr(k)+vprm1(k))*&
           current_state%sth%data(k,current_state%column_local_y,current_state%column_local_x)+0.5*(thlpr(k)+thlprp1(k))*&
           current_state%sv%data(k,current_state%column_local_y,current_state%column_local_x)
      if (allocated(v_thetal_advection)) v_thetal_advection(k)=v_thetal_advection(k)+0.5*(vpr(k)+vprm1(k))*&
           th_advection%real_1d_array(k)+0.5*(thlpr(k)+thlprp1(k))*v_advection%real_1d_array(k)
      if (allocated(v_thetal_viscosity_diffusion)) v_thetal_viscosity_diffusion(k)=v_thetal_viscosity_diffusion(k)+0.5*&
           (thlpr(k)+thlprp1(k))*v_viscosity%real_1d_array(k)+0.5*(vpr(k)+vprm1(k))*th_diffusion%real_1d_array(k)
      if (allocated(wv_thetal)) wv_thetal(k)=wv_thetal(k)+&
           current_state%w%data(k,current_state%column_local_y,current_state%column_local_x)*0.25*(vpr(k+1)+vpr(k)+&
           vprm1(k+1)+vprm1(k))*0.5*(thlpr(k+1)+thlpr(k))

      if (allocated(w_thetal)) w_thetal(k)=w_thetal(k)+&
           current_state%w%data(k,current_state%column_local_y,current_state%column_local_x)*0.5*(thlpr(k)+thlpr(k+1))
      if (allocated(ws_thetal)) ws_thetal(k)=ws_thetal(k)+&
           current_state%w%data(k,current_state%column_local_y,current_state%column_local_x)*0.5*&
           (current_state%sth%data(k,current_state%column_local_y,current_state%column_local_x)+&
           current_state%sth%data(k+1,current_state%column_local_y,current_state%column_local_x))+0.5*(thlpr(k)+thlpr(k+1))*&
           current_state%sw%data(k,current_state%column_local_y,current_state%column_local_x)
      if (allocated(w_thetal_advection)) w_thetal_advection(k)=w_thetal_advection(k)+&
           current_state%w%data(k,current_state%column_local_y,current_state%column_local_x)*0.5*(th_advection%real_1d_array(k)+&
           th_advection%real_1d_array(k+1))+0.5*(thlpr(k)+thlpr(k+1))*w_advection%real_1d_array(k)
      if (allocated(w_thetal_viscosity_diffusion)) w_thetal_viscosity_diffusion(k)=w_thetal_viscosity_diffusion(k)+0.5*&
           (thlpr(k)+thlpr(k+1))*w_viscosity%real_1d_array(k)+&
           current_state%w%data(k,current_state%column_local_y,current_state%column_local_x)*0.5*&
           (th_diffusion%real_1d_array(k)+th_diffusion%real_1d_array(k+1))
      if (allocated(w_thetal_buoyancy)) w_thetal_buoyancy(k)=w_thetal_buoyancy(k)+0.5*(thlpr(k+1)+thlpr(k))*&
           w_buoyancy%real_1d_array(k)
      if (allocated(wv_thetal)) ww_thetal(k)=ww_thetal(k)+0.5*&
           (current_state%w%data(k,current_state%column_local_y,current_state%column_local_x)+&
           current_state%w%data(k-1,current_state%column_local_y,current_state%column_local_x))*0.5*(&
           current_state%w%data(k,current_state%column_local_y,current_state%column_local_x)+&
           current_state%w%data(k-1,current_state%column_local_y,current_state%column_local_x))*thlpr(k)

      if (allocated(thetal_thetal)) thetal_thetal(k)=thetal_thetal(k)+thlpr(k)*thlpr(k)
      if (allocated(sthetal_thetal)) sthetal_thetal(k)=sthetal_thetal(k)+2.0*thlpr(k)*&
           current_state%sth%data(k,current_state%column_local_y,current_state%column_local_x)
      if (allocated(thetal_thetal_advection)) thetal_thetal_advection(k)=thetal_thetal_advection(k)+2.0*thlpr(k)*&
           th_advection%real_1d_array(k)
      if (allocated(thetal_thetal_diffusion)) thetal_thetal_diffusion(k)=thetal_thetal_diffusion(k)+2.0*thlpr(k)*&
           th_diffusion%real_1d_array(k)
      if (allocated(wthetal_thetal)) wthetal_thetal(k)=wthetal_thetal(k)+&
           current_state%w%data(k,current_state%column_local_y,current_state%column_local_x)*0.5*(thlpr(k+1)+thlpr(k))*0.5*(&
           thlpr(k+1)+thlpr(k)) 
    end do

    if (allocated(u_advection%real_1d_array)) deallocate(u_advection%real_1d_array)
    if (allocated(u_viscosity%real_1d_array)) deallocate(u_viscosity%real_1d_array)
    if (allocated(th_advection%real_1d_array)) deallocate(th_advection%real_1d_array)
    if (allocated(th_diffusion%real_1d_array)) deallocate(th_diffusion%real_1d_array)
    if (allocated(v_advection%real_1d_array)) deallocate(v_advection%real_1d_array)
    if (allocated(v_viscosity%real_1d_array)) deallocate(v_viscosity%real_1d_array)
    if (allocated(w_advection%real_1d_array)) deallocate(w_advection%real_1d_array)
    if (allocated(w_viscosity%real_1d_array)) deallocate(w_viscosity%real_1d_array)
    if (allocated(w_buoyancy%real_1d_array)) deallocate(w_buoyancy%real_1d_array)
  end subroutine compute_thetal_for_column  

  !> Clears the prognostic (uu, vv, ww) budgets
  subroutine clear_prognostic_budgets()
    if (allocated(tu_su)) tu_su=0.0_DEFAULT_PRECISION
    if (allocated(uu_advection)) uu_advection=0.0_DEFAULT_PRECISION
    if (allocated(uu_viscosity)) uu_viscosity=0.0_DEFAULT_PRECISION
    if (allocated(wu_u)) wu_u=0.0_DEFAULT_PRECISION
    if (allocated(tv_sv)) tv_sv=0.0_DEFAULT_PRECISION
    if (allocated(vv_advection)) vv_advection=0.0_DEFAULT_PRECISION
    if (allocated(vv_viscosity)) vv_viscosity=0.0_DEFAULT_PRECISION
    if (allocated(wv_v)) wv_v=0.0_DEFAULT_PRECISION
    if (allocated(tw_sw)) tw_sw=0.0_DEFAULT_PRECISION
    if (allocated(ww_advection)) ww_advection=0.0_DEFAULT_PRECISION
    if (allocated(ww_viscosity)) ww_viscosity=0.0_DEFAULT_PRECISION
    if (allocated(ww_buoyancy)) ww_buoyancy=0.0_DEFAULT_PRECISION
  end subroutine clear_prognostic_budgets

  !> Computes the prognostic (uu, vv, ww) budgets for a specific column
  !! @param current_state The current model state
  subroutine compute_prognostic_budgets_for_column(current_state)
    type(model_state_type), target, intent(inout) :: current_state

    real(kind=DEFAULT_PRECISION), dimension(current_state%local_grid%size(Z_INDEX)) :: upr, vpr, uprm1, vprm1
    type(component_field_value_type) :: u_advection, u_viscosity, v_advection, v_viscosity, w_advection, w_viscosity, w_buoyancy
    integer :: k

    if (is_component_field_available("u_advection")) u_advection=get_component_field_value(current_state, "u_advection")
    if (is_component_field_available("u_viscosity")) u_viscosity=get_component_field_value(current_state, "u_viscosity")
    if (is_component_field_available("v_advection")) v_advection=get_component_field_value(current_state, "v_advection")
    if (is_component_field_available("v_viscosity")) v_viscosity=get_component_field_value(current_state, "v_viscosity")
    if (is_component_field_available("w_advection")) w_advection=get_component_field_value(current_state, "w_advection")
    if (is_component_field_available("w_viscosity")) w_viscosity=get_component_field_value(current_state, "w_viscosity")
    if (is_component_field_available("w_buoyancy")) w_buoyancy=get_component_field_value(current_state, "w_buoyancy")

    do k=1, current_state%local_grid%size(Z_INDEX)
      upr(k)=current_state%u%data(k,current_state%column_local_y,current_state%column_local_x)
      uprm1(k)=current_state%u%data(k,current_state%column_local_y,current_state%column_local_x-1)
      if (allocated(current_state%global_grid%configuration%vertical%olubar)) then
        upr(k)=upr(k)-(current_state%global_grid%configuration%vertical%olubar(k)-current_state%ugal)
        uprm1(k)=uprm1(k)-(current_state%global_grid%configuration%vertical%olubar(k)-current_state%ugal)
      end if
      vpr(k)=current_state%v%data(k,current_state%column_local_y,current_state%column_local_x)
      vprm1(k)=current_state%v%data(k,current_state%column_local_y-1,current_state%column_local_x)
      if (allocated(current_state%global_grid%configuration%vertical%olvbar)) then
        vpr(k)=vpr(k)-(current_state%global_grid%configuration%vertical%olvbar(k)-current_state%vgal)
        vprm1(k)=vprm1(k)-(current_state%global_grid%configuration%vertical%olvbar(k)-current_state%vgal)
      end if
    end do
    do k=2, current_state%local_grid%size(Z_INDEX)-1
      if (allocated(tu_su)) tu_su(k)=tu_su(k)+2.0*upr(k)*current_state%su%data(k,current_state%column_local_y,&
           current_state%column_local_x)
      if (allocated(uu_advection)) uu_advection(k)=uu_advection(k)+2.0*upr(k)*u_advection%real_1d_array(k)
      if (allocated(uu_viscosity)) uu_viscosity(k)=uu_viscosity(k)+2.0*upr(k)*u_viscosity%real_1d_array(k)
      if (allocated(wu_u)) wu_u(k)=wu_u(k)+0.25*(upr(k)+upr(k+1)+uprm1(k)+uprm1(k+1))*0.25*&
           (upr(k)+upr(k+1)+uprm1(k)+uprm1(k+1))*current_state%w%data(k,current_state%column_local_y,&
           current_state%column_local_x)
      if (allocated(tv_sv)) tv_sv(k)=tv_sv(k)+2.0*vpr(k)*current_state%sv%data(k,current_state%column_local_y,&
           current_state%column_local_x)
      if (allocated(vv_advection)) vv_advection(k)=vv_advection(k)+2.0*vpr(k)*v_advection%real_1d_array(k)
      if (allocated(vv_viscosity)) uu_viscosity(k)=vv_viscosity(k)+2.0*vpr(k)*v_viscosity%real_1d_array(k)
      if (allocated(wv_v)) wv_v(k)=wv_v(k)+0.25*(vpr(k)+vpr(k+1)+vprm1(k)+vprm1(k+1))*0.25*&
           (vpr(k)+vpr(k+1)+vprm1(k)+vprm1(k+1))*current_state%w%data(k,current_state%column_local_y,&
           current_state%column_local_x)
      if (allocated(tw_sw)) tw_sw(k)=tw_sw(k)+2.0*current_state%w%data(k,current_state%column_local_y,&
           current_state%column_local_x)*current_state%sw%data(k,current_state%column_local_y,&
           current_state%column_local_x)
      if (allocated(ww_advection)) ww_advection(k)=ww_advection(k)+2.0*current_state%w%data(k,current_state%column_local_y,&
           current_state%column_local_x)*w_advection%real_1d_array(k)
      if (allocated(ww_viscosity)) ww_viscosity(k)=ww_viscosity(k)+2.0*current_state%w%data(k,current_state%column_local_y,&
           current_state%column_local_x)*w_viscosity%real_1d_array(k)
      if (allocated(ww_buoyancy)) ww_buoyancy(k)=ww_buoyancy(k)+2.0*current_state%w%data(k,current_state%column_local_y,&
           current_state%column_local_x)*w_buoyancy%real_1d_array(k)
    end do

    if (allocated(u_advection%real_1d_array)) deallocate(u_advection%real_1d_array)
    if (allocated(u_viscosity%real_1d_array)) deallocate(u_viscosity%real_1d_array)
    if (allocated(v_advection%real_1d_array)) deallocate(v_advection%real_1d_array)
    if (allocated(v_viscosity%real_1d_array)) deallocate(v_viscosity%real_1d_array)
    if (allocated(w_advection%real_1d_array)) deallocate(w_advection%real_1d_array)
    if (allocated(w_viscosity%real_1d_array)) deallocate(w_viscosity%real_1d_array)
    if (allocated(w_buoyancy%real_1d_array)) deallocate(w_buoyancy%real_1d_array)
  end subroutine compute_prognostic_budgets_for_column
  

  !> Clears the uw uv diagnostics
  subroutine clear_uw_vw()
    if (allocated(uw_advection)) uw_advection=0.0_DEFAULT_PRECISION
    if (allocated(vw_advection)) vw_advection=0.0_DEFAULT_PRECISION
    if (allocated(uw_viscosity)) uw_viscosity=0.0_DEFAULT_PRECISION
    if (allocated(vw_viscosity)) vw_viscosity=0.0_DEFAULT_PRECISION
    if (allocated(uw_buoyancy)) uw_buoyancy=0.0_DEFAULT_PRECISION
    if (allocated(vw_buoyancy)) vw_buoyancy=0.0_DEFAULT_PRECISION
    if (allocated(uw_tendency)) uw_tendency=0.0_DEFAULT_PRECISION
    if (allocated(vw_tendency)) vw_tendency=0.0_DEFAULT_PRECISION
    if (allocated(uw_w)) uw_w=0.0_DEFAULT_PRECISION
    if (allocated(vw_w)) vw_w=0.0_DEFAULT_PRECISION    
  end subroutine clear_uw_vw
  !> Clears the TKE diagnostics
  subroutine clear_TKE()
    if (allocated(sres)) sres=0.0_DEFAULT_PRECISION
    if (allocated(wke)) wke=0.0_DEFAULT_PRECISION
    if (allocated(wp)) wp=0.0_DEFAULT_PRECISION
    if (allocated(buoy)) buoy=0.0_DEFAULT_PRECISION
    if (allocated(tend)) tend=0.0_DEFAULT_PRECISION
  end subroutine clear_TKE

  !> Computes the uw uv diagnostics for a specific column
  !! @param current_state Current model state
  subroutine compute_uw_vw_for_column(current_state)
    type(model_state_type), target, intent(inout) :: current_state

    real(kind=DEFAULT_PRECISION), dimension(current_state%local_grid%size(Z_INDEX)) :: upr, vpr, uprm1, vprm1
    type(component_field_value_type) :: w_advection, v_advection, u_advection, w_viscosity, v_viscosity, u_viscosity, &
         w_buoyancy
    integer :: k

    if (is_component_field_available("w_advection")) w_advection=get_component_field_value(current_state, "w_advection")
    if (is_component_field_available("v_advection")) v_advection=get_component_field_value(current_state, "v_advection")
    if (is_component_field_available("u_advection")) u_advection=get_component_field_value(current_state, "u_advection")
    if (is_component_field_available("w_viscosity")) w_viscosity=get_component_field_value(current_state, "w_viscosity")
    if (is_component_field_available("v_viscosity")) v_viscosity=get_component_field_value(current_state, "v_viscosity")
    if (is_component_field_available("u_viscosity")) u_viscosity=get_component_field_value(current_state, "u_viscosity")
    if (is_component_field_available("w_buoyancy")) w_buoyancy=get_component_field_value(current_state, "w_buoyancy")

    do k=1, current_state%local_grid%size(Z_INDEX)
      upr(k)=current_state%u%data(k,current_state%column_local_y,current_state%column_local_x)
      uprm1(k)=current_state%u%data(k,current_state%column_local_y,current_state%column_local_x-1)
      if (allocated(current_state%global_grid%configuration%vertical%olubar)) then
        upr(k)=upr(k)-(current_state%global_grid%configuration%vertical%olubar(k)-current_state%ugal)
        uprm1(k)=uprm1(k)-(current_state%global_grid%configuration%vertical%olubar(k)-current_state%ugal)
      end if
      vpr(k)=current_state%v%data(k,current_state%column_local_y,current_state%column_local_x)
      vprm1(k)=current_state%v%data(k,current_state%column_local_y-1,current_state%column_local_x)
      if (allocated(current_state%global_grid%configuration%vertical%olvbar)) then
        vpr(k)=vpr(k)-(current_state%global_grid%configuration%vertical%olvbar(k)-current_state%vgal)
        vprm1(k)=vprm1(k)-(current_state%global_grid%configuration%vertical%olvbar(k)-current_state%vgal)
      end if
    end do
    do k=2, current_state%local_grid%size(Z_INDEX)-1
      if (allocated(uw_advection)) uw_advection(k)=uw_advection(k)+0.5*(upr(k)+uprm1(k))*0.5*&
           (w_advection%real_1d_array(k)+w_advection%real_1d_array(k-1))+0.25*(&
           current_state%w%data(k,current_state%column_local_y,current_state%column_local_x)+&
           current_state%w%data(k-1,current_state%column_local_y,current_state%column_local_x)+&
           current_state%w%data(k,current_state%column_local_y,current_state%column_local_x+1)+&
           current_state%w%data(k-1,current_state%column_local_y,current_state%column_local_x+1))*u_advection%real_1d_array(k)
      if (allocated(vw_advection)) vw_advection(k)=vw_advection(k)+0.5*(vpr(k)+vprm1(k))*0.5*&
           (w_advection%real_1d_array(k)+w_advection%real_1d_array(k-1))+0.25*(&
           current_state%w%data(k,current_state%column_local_y,current_state%column_local_x)+&
           current_state%w%data(k-1,current_state%column_local_y,current_state%column_local_x)+&
           current_state%w%data(k,current_state%column_local_y,current_state%column_local_x+1)+&
           current_state%w%data(k-1,current_state%column_local_y,current_state%column_local_x+1))*v_advection%real_1d_array(k)
      if (allocated(uw_viscosity)) uw_viscosity(k)=uw_viscosity(k)+0.5*(upr(k)+uprm1(k))*0.5*&
           (w_viscosity%real_1d_array(k)+w_viscosity%real_1d_array(k-1))+0.25*(&
           current_state%w%data(k,current_state%column_local_y,current_state%column_local_x)+&
           current_state%w%data(k-1,current_state%column_local_y,current_state%column_local_x)+&
           current_state%w%data(k,current_state%column_local_y,current_state%column_local_x+1)+&
           current_state%w%data(k-1,current_state%column_local_y,current_state%column_local_x+1))*u_viscosity%real_1d_array(k)
      if (allocated(vw_viscosity)) vw_viscosity(k)=vw_viscosity(k)+0.5*(vpr(k)+vprm1(k))*0.5*&
           (w_viscosity%real_1d_array(k)+w_viscosity%real_1d_array(k-1))+0.25*(&
           current_state%w%data(k,current_state%column_local_y,current_state%column_local_x)+&
           current_state%w%data(k-1,current_state%column_local_y,current_state%column_local_x)+&
           current_state%w%data(k,current_state%column_local_y,current_state%column_local_x+1)+&
           current_state%w%data(k-1,current_state%column_local_y,current_state%column_local_x+1))*v_viscosity%real_1d_array(k)
      if (allocated(uw_buoyancy)) uw_buoyancy(k)=uw_buoyancy(k)+0.5*(upr(k)+uprm1(k))*0.5*(&
           w_buoyancy%real_1d_array(k)+w_buoyancy%real_1d_array(k-1)) 
      if (allocated(vw_buoyancy)) vw_buoyancy(k)=vw_buoyancy(k)+0.5*(vpr(k)+vprm1(k))*0.5*(&
           w_buoyancy%real_1d_array(k)+w_buoyancy%real_1d_array(k-1))
      if (allocated(uw_tendency)) uw_tendency(k)=uw_tendency(k)+0.25*(&
           current_state%w%data(k,current_state%column_local_y,current_state%column_local_x)+&
           current_state%w%data(k-1,current_state%column_local_y,current_state%column_local_x)+&
           current_state%w%data(k,current_state%column_local_y,current_state%column_local_x+1)+&
           current_state%w%data(k-1,current_state%column_local_y,current_state%column_local_x+1))*&
           current_state%su%data(k,current_state%column_local_y,current_state%column_local_x)+0.5*(&
           current_state%sw%data(k,current_state%column_local_y,current_state%column_local_x)+&
           current_state%sw%data(k-1,current_state%column_local_y,current_state%column_local_x))*0.5*(upr(k)+uprm1(k))
      if (allocated(vw_tendency)) vw_tendency(k)=vw_tendency(k)+0.25*(&
           current_state%w%data(k,current_state%column_local_y,current_state%column_local_x)+&
           current_state%w%data(k-1,current_state%column_local_y,current_state%column_local_x)+&
           current_state%w%data(k,current_state%column_local_y,current_state%column_local_x+1)+&
           current_state%w%data(k-1,current_state%column_local_y,current_state%column_local_x+1))*&
           current_state%sv%data(k,current_state%column_local_y,current_state%column_local_x)+0.5*(&
           current_state%sw%data(k,current_state%column_local_y,current_state%column_local_x)+&
           current_state%sw%data(k-1,current_state%column_local_y,current_state%column_local_x))*0.5*(vpr(k)+vprm1(k))
      if (allocated(uw_w)) uw_w(k)=uw_w(k)+0.25*(upr(k)+upr(k+1)+uprm1(k)+uprm1(k+1))*&
           current_state%w%data(k,current_state%column_local_y,current_state%column_local_x)*&
           current_state%w%data(k,current_state%column_local_y,current_state%column_local_x)
      if (allocated(vw_w)) vw_w(k)=vw_w(k)+0.25*(vpr(k)+vpr(k+1)+vprm1(k)+vprm1(k+1))*&
           current_state%w%data(k,current_state%column_local_y,current_state%column_local_x)*&
           current_state%w%data(k,current_state%column_local_y,current_state%column_local_x)
    end do

    if (allocated(u_advection%real_1d_array)) deallocate(u_advection%real_1d_array)
    if (allocated(u_viscosity%real_1d_array)) deallocate(u_viscosity%real_1d_array)
    if (allocated(v_advection%real_1d_array)) deallocate(v_advection%real_1d_array)
    if (allocated(v_viscosity%real_1d_array)) deallocate(v_viscosity%real_1d_array)
    if (allocated(w_advection%real_1d_array)) deallocate(w_advection%real_1d_array)
    if (allocated(w_viscosity%real_1d_array)) deallocate(w_viscosity%real_1d_array)
    if (allocated(w_buoyancy%real_1d_array)) deallocate(w_buoyancy%real_1d_array)
  end subroutine compute_uw_vw_for_column
     

  !> Computes the TKE diagnostics for a specific column. 
  !! @param current_state The current model state
  subroutine compute_TKE_for_column(current_state)
  type(model_state_type), target, intent(inout) :: current_state

    real(kind=DEFAULT_PRECISION), dimension(current_state%local_grid%size(Z_INDEX)) :: upr, vpr, uprm1, vprm1, &
          uu_tendency,vv_tendency,ww_tendency
    real(kind=DEFAULT_PRECISION), dimension(current_state%local_grid%size(Z_INDEX)) :: umean, wu_umean, vmean, wv_vmean, &
          w_pprime_at_p, rke1, w_qvprime_at_w, w_qclprime_at_w, w_thprime_at_w, wq, rho, rec_rho, rhon, rec_rhon, &
          uw_tot, vw_tot,w_upr_at_w,w_vpr_at_w, w_buoyancy
    real(kind=DEFAULT_PRECISION) :: u_at_p, v_at_p, w_at_p
    real(kind=DEFAULT_PRECISION) ::  C_virtual
   
    integer :: k, n

    C_virtual = (ratio_mol_wts-1.0_DEFAULT_PRECISION)
    
    !Resolved diagnostics
! ***********************Buoyant production 1/2 ***************************
       
    do k=1, current_state%local_grid%size(Z_INDEX)
      rho(k)=current_state%global_grid%configuration%vertical%rho(k)
      rhon(k)=current_state%global_grid%configuration%vertical%rhon(k)
      rec_rho(k)=1.0_DEFAULT_PRECISION/rho(k)
      rec_rhon(k)=1.0_DEFAULT_PRECISION/rhon(k)
    enddo

   
! ***********************TKE Tendency ***************************  

    do k=2, current_state%local_grid%size(Z_INDEX)

      uu_tendency(k) = ((current_state%u%data(k,current_state%column_local_y,current_state%column_local_x) - &
                         current_state%global_grid%configuration%vertical%olubar(k) )**2 - &
                        (current_state%zu%data(k,current_state%column_local_y,current_state%column_local_x) - &
                         current_state%global_grid%configuration%vertical%olzubar(k) )**2 ) / &
                        current_state%dtm
                        
      vv_tendency(k) = ((current_state%v%data(k,current_state%column_local_y,current_state%column_local_x) - &
                         current_state%global_grid%configuration%vertical%olvbar(k) )**2 - &
                        (current_state%zv%data(k,current_state%column_local_y,current_state%column_local_x) - &
                         current_state%global_grid%configuration%vertical%olzvbar(k) )**2 ) / &
                        current_state%dtm                        

      ww_tendency(k) = (current_state%w%data(k,current_state%column_local_y,current_state%column_local_x) * &
                        current_state%w%data(k,current_state%column_local_y,current_state%column_local_x) - &
                        current_state%zw%data(k,current_state%column_local_y,current_state%column_local_x) * &
                        current_state%zw%data(k,current_state%column_local_y,current_state%column_local_x)) / &
                        current_state%dtm

    enddo

    uu_tendency(1) = -uu_tendency(2)
    vv_tendency(1) = -vv_tendency(2)
    
    if (allocated(tend)) then
      do k=2, current_state%local_grid%size(Z_INDEX)-1
        tend(k)=tend(k) + 0.5_DEFAULT_PRECISION * (& 
          0.5_DEFAULT_PRECISION * (uu_tendency(k)+uu_tendency(k+1)) + &
          0.5_DEFAULT_PRECISION * (vv_tendency(k)+vv_tendency(k+1)) + &
          ww_tendency(k) )
      enddo
    end if
        
! ***********************Shear production ***************************

    do k=2, current_state%local_grid%size(Z_INDEX)-1

      umean(k)=(current_state%global_grid%configuration%vertical%olubar(k+1) -&
                current_state%global_grid%configuration%vertical%olubar(k))* &
                current_state%global_grid%configuration%vertical%rdzn(k+1)

      w_upr_at_w(k) =current_state%w%data(k,current_state%column_local_y,current_state%column_local_x) * &
          (0.25_DEFAULT_PRECISION * ( &
          (current_state%u%data(k,current_state%column_local_y,current_state%column_local_x-1)  - &
           current_state%global_grid%configuration%vertical%olubar(k)) + &
          (current_state%u%data(k+1,current_state%column_local_y,current_state%column_local_x)  - &
           current_state%global_grid%configuration%vertical%olubar(k+1)) + &
          (current_state%u%data(k,current_state%column_local_y,current_state%column_local_x)    - &
           current_state%global_grid%configuration%vertical%olubar(k)) + &
          (current_state%u%data(k+1,current_state%column_local_y,current_state%column_local_x-1)- &
           current_state%global_grid%configuration%vertical%olubar(k+1)) ) )

      vmean(k)=(current_state%global_grid%configuration%vertical%olvbar(k+1) - &
                current_state%global_grid%configuration%vertical%olvbar(k)) * &
               current_state%global_grid%configuration%vertical%rdzn(k+1)

      w_vpr_at_w(k) =current_state%w%data(k,current_state%column_local_y,current_state%column_local_x) * &
          (0.25_DEFAULT_PRECISION * ( &
          (current_state%v%data(k,current_state%column_local_y-1,current_state%column_local_x)  - &
           current_state%global_grid%configuration%vertical%olvbar(k)) + &
          (current_state%v%data(k+1,current_state%column_local_y,current_state%column_local_x)  - &
           current_state%global_grid%configuration%vertical%olvbar(k+1)) + &
          (current_state%v%data(k,current_state%column_local_y,current_state%column_local_x)    - &
           current_state%global_grid%configuration%vertical%olvbar(k)) + &
          (current_state%v%data(k+1,current_state%column_local_y-1,current_state%column_local_x)- &
           current_state%global_grid%configuration%vertical%olvbar(k+1)) ) )         

      wu_umean(k)=(w_upr_at_w(k)*umean(k))
      wv_vmean(k)= (w_vpr_at_w(k)*vmean(k))

      if (allocated(sres)) then
        sres(k)=sres(k) - (wv_vmean(k) + wu_umean(k))
      end if
  
    end do
    sres(1)=0.0_DEFAULT_PRECISION
    sres(current_state%local_grid%size(Z_INDEX))=0.0_DEFAULT_PRECISION
    
    do k=2, current_state%local_grid%size(Z_INDEX)

! *********************** Pressure transport ***************************	
        !In current state - p=p/rho (rho here is on same levels as p)
        ! Note - calculating on z levels (i.e. w) 
        ! So need w'p' on p levels
        
      w_pprime_at_p(k) = 0.5_DEFAULT_PRECISION * & 
         (current_state%w%data(k,  current_state%column_local_y,current_state%column_local_x) + &
          current_state%w%data(k-1,current_state%column_local_y,current_state%column_local_x)) * &
         (current_state%global_grid%configuration%vertical%rhon(k) * &
          current_state%p%data(k,current_state%column_local_y,current_state%column_local_x))
        
      u_at_p = 0.5_DEFAULT_PRECISION * &
          ((current_state%u%data(k,current_state%column_local_y,current_state%column_local_x-1)- &
            current_state%global_grid%configuration%vertical%olubar(k)) + &
           (current_state%u%data(k,current_state%column_local_y,current_state%column_local_x)  - &
            current_state%global_grid%configuration%vertical%olubar(k))) 
   
      v_at_p = 0.5_DEFAULT_PRECISION * &
          ((current_state%v%data(k,current_state%column_local_y-1,current_state%column_local_x)- &
            current_state%global_grid%configuration%vertical%olvbar(k)) + &
           (current_state%v%data(k,current_state%column_local_y,current_state%column_local_x)  - &
            current_state%global_grid%configuration%vertical%olvbar(k)))

      w_at_p = 0.5_DEFAULT_PRECISION  * &
          (current_state%w%data(k,  current_state%column_local_y,current_state%column_local_x) + &
           current_state%w%data(k-1,current_state%column_local_y,current_state%column_local_x))
   
      rke1(k)= 0.5_DEFAULT_PRECISION * w_at_p * &
          ( u_at_p*u_at_p + v_at_p*v_at_p + w_at_p*w_at_p) * rec_rhon(k)

    end do
           
    w_pprime_at_p(current_state%local_grid%size(Z_INDEX)) = 0.0_DEFAULT_PRECISION 
    ! Zero gradient at surface     
    w_pprime_at_p(1)=w_pprime_at_p(2)
    if (allocated(wp)) then
      do k=1, current_state%local_grid%size(Z_INDEX)-1
        wp(k)= wp(k) - ((w_pprime_at_p(k+1) - w_pprime_at_p(k)) * &
          current_state%global_grid%configuration%vertical%rdzn(k+1) * rec_rho(k))
      end do
    end if
    
! ********************** Resolved turbulent transport ************************

    rke1(current_state%local_grid%size(Z_INDEX)) = 0.0_DEFAULT_PRECISION
    ! Zero gradient at surface    
    rke1(1)=rke1(2)
    
    if (allocated(wke)) then
      do k=1, current_state%local_grid%size(Z_INDEX)-1
        wke(k) = wke(k) -(rho(k) * (rke1(k+1) - rke1(k) ) * &
          current_state%global_grid%configuration%vertical%rdzn(k+1))
      end do
    end if
            
! *********************** Subgrid buoyant production*************************** 
!!!Using buoyancy.F90

#ifdef W_ACTIVE
    if (.not. current_state%passive_th .and. current_state%th%active) then
      do k=2,current_state%local_grid%size(Z_INDEX)-1    
        w_buoyancy(k)=(0.5_DEFAULT_PRECISION*current_state%global_grid%configuration%vertical%buoy_co(k)) * &
             ((current_state%th%data(k, current_state%column_local_y, current_state%column_local_x) -       &
               current_state%global_grid%configuration%vertical%olthbar(k)) +                               &
              (current_state%th%data(k+1, current_state%column_local_y, current_state%column_local_x) -     &
               current_state%global_grid%configuration%vertical%olthbar(k+1)))
       end do
    end if

    if (.not. current_state%passive_q .and. current_state%number_q_fields .gt. 0) then
      if (current_state%use_anelastic_equations) then                                                      
        do n=1,current_state%number_q_fields
          do k=2,current_state%local_grid%size(Z_INDEX)-1  
            w_buoyancy(k) = w_buoyancy(k) + &       
                 (0.5_DEFAULT_PRECISION*current_state%global_grid%configuration%vertical%buoy_co(k)) * &
                 current_state%cq(n) * &
                 (current_state%global_grid%configuration%vertical%thref(k)*&
                  (current_state%q(n)%data(k, current_state%column_local_y, current_state%column_local_x) - &
                   current_state%global_grid%configuration%vertical%olqbar(k,n)) + & 
                  current_state%global_grid%configuration%vertical%thref(k+1) * &
                  (current_state%q(n)%data(k+1, current_state%column_local_y, current_state%column_local_x) - &
                   current_state%global_grid%configuration%vertical%olqbar(k+1,n)))
          end do
        end do
      else                                                                     
        do n=1,current_state%number_q_fields
          do k=2,current_state%local_grid%size(Z_INDEX)-1
            w_buoyancy(k) = w_buoyancy(k) + & 
                  G*0.5_DEFAULT_PRECISION*current_state%cq(n)*&
                  (current_state%q(n)%data(k, current_state%column_local_y, current_state%column_local_x) -&
                   current_state%global_grid%configuration%vertical%olqbar(k,n) + & 
                   current_state%q(n)%data(k+1, current_state%column_local_y, current_state%column_local_x)- &
                   current_state%global_grid%configuration%vertical%olqbar(k,n))
          end do
        end do
      end if
    end if

    do k=2, current_state%local_grid%size(Z_INDEX)-1
      if (allocated(buoy)) then
        buoy(k) = buoy(k) + &
        current_state%w%data(k,current_state%column_local_y,current_state%column_local_x) * &
        w_buoyancy(k)
      end if
    end do
   
#endif
   
    buoy(1) = 0.0_DEFAULT_PRECISION
    buoy(current_state%local_grid%size(Z_INDEX)) = 0.0_DEFAULT_PRECISION
    
  end subroutine compute_TKE_for_column

  !> Clears the Q flux diagnostics, called at the start of a timestep
  subroutine clear_q_fluxes()
    if (allocated(q_flux_values)) q_flux_values=0.0_DEFAULT_PRECISION
    if (allocated(q_gradient)) q_gradient=0.0_DEFAULT_PRECISION
    if (allocated(q_diff)) q_diff=0.0_DEFAULT_PRECISION
    if (allocated(q_buoyancy)) q_buoyancy=0.0_DEFAULT_PRECISION
    if (allocated(q_tendency)) q_tendency=0.0_DEFAULT_PRECISION
  end subroutine clear_q_fluxes  

  !> Computes the Q flux diagnostics for a specific column
  !! @param current_state Current model state
  subroutine compute_q_flux_for_column(current_state)
    type(model_state_type), target, intent(inout) :: current_state

    integer :: k, n
    real(kind=DEFAULT_PRECISION), dimension(current_state%local_grid%size(Z_INDEX)) :: qpr
    type(component_field_value_type) :: w_advection_published_value, q_advection_published_value, w_viscosity_published_value, &
         q_diffusion_published_value, w_buoyancy_published_value

    if (allocated(q_gradient)) then
      w_advection_published_value=get_component_field_value(current_state, "w_advection")
      q_advection_published_value=get_component_field_value(current_state, "q_advection")
    end if
    if (allocated(q_diff)) then
      w_viscosity_published_value=get_component_field_value(current_state, "w_viscosity")
      q_diffusion_published_value=get_component_field_value(current_state, "q_diffusion")
    end if
    if (allocated(q_buoyancy)) then
      w_buoyancy_published_value=get_component_field_value(current_state, "w_buoyancy")
    end if

    do n=1, current_state%number_q_fields
      do k=1, current_state%local_grid%size(Z_INDEX)
        qpr(k)=current_state%q(n)%data(k,current_state%column_local_y,current_state%column_local_x)
        if (allocated(current_state%global_grid%configuration%vertical%olqbar)) then
          qpr(k)=qpr(k)-current_state%global_grid%configuration%vertical%olqbar(k,n)
        end if        
      end do
      do k=2, current_state%local_grid%size(Z_INDEX)-1
        if (allocated(q_flux_values)) q_flux_values(k,n)=q_flux_values(k,n)+&
             current_state%w%data(k,current_state%column_local_y,current_state%column_local_x)*&
           current_state%w%data(k,current_state%column_local_y,current_state%column_local_x)*0.5*(qpr(k)+qpr(k+1))
        if (allocated(q_tendency)) q_tendency(k,n)=q_tendency(k,n)+0.5*&
           (current_state%w%data(k,current_state%column_local_y,current_state%column_local_x)+&
           current_state%w%data(k-1,current_state%column_local_y,current_state%column_local_x))*&
           current_state%sq(n)%data(k,current_state%column_local_y,current_state%column_local_x)+0.5*&
           qpr(k)*(current_state%sw%data(k,current_state%column_local_y,current_state%column_local_x)+&
           current_state%sw%data(k-1,current_state%column_local_y,current_state%column_local_x))
        if (allocated(q_gradient)) then
          q_gradient(k,n)=q_gradient(k,n)+qpr(k)*0.5*(w_advection_published_value%real_1d_array(k)+&
               w_advection_published_value%real_1d_array(k-1))+0.5*&
               (current_state%w%data(k,current_state%column_local_y,current_state%column_local_x)+&
               current_state%w%data(k-1,current_state%column_local_y,current_state%column_local_x))*&
               q_advection_published_value%real_2d_array(k, n)
        end if
        if (allocated(q_diff)) then
          q_diff(k,n)=q_diff(k,n)+qpr(k)*0.5*(w_viscosity_published_value%real_1d_array(k)+&
               w_viscosity_published_value%real_1d_array(k-1))+&
               0.5*(current_state%w%data(k,current_state%column_local_y,current_state%column_local_x)+&
               current_state%w%data(k-1,current_state%column_local_y,current_state%column_local_x))*&
               q_diffusion_published_value%real_2d_array(k,n)
        end if
        if (allocated(q_buoyancy)) then
          q_buoyancy(k,n)=q_buoyancy(k,n)+qpr(k)*0.5*(w_buoyancy_published_value%real_1d_array(k)+&
               w_buoyancy_published_value%real_1d_array(k-1))
        end if
     end do
    end do
    if (allocated(w_advection_published_value%real_1d_array)) deallocate(w_advection_published_value%real_1d_array)
    if (allocated(q_advection_published_value%real_2d_array)) deallocate(q_advection_published_value%real_2d_array)
    if (allocated(w_viscosity_published_value%real_1d_array)) deallocate(w_viscosity_published_value%real_1d_array)
    if (allocated(q_diffusion_published_value%real_2d_array)) deallocate(q_diffusion_published_value%real_2d_array)
    if (allocated(w_buoyancy_published_value%real_1d_array)) deallocate(w_buoyancy_published_value%real_1d_array)
  end subroutine compute_q_flux_for_column

  !> Clears the heat flux diagnostics at the start of a timestep
  subroutine clear_theta_fluxes()
    if (allocated(th_flux_values)) th_flux_values=0.0_DEFAULT_PRECISION
    if (allocated(th_tendency)) th_tendency=0.0_DEFAULT_PRECISION
    if (allocated(th_gradient)) th_gradient=0.0_DEFAULT_PRECISION
    if (allocated(th_diff)) th_diff=0.0_DEFAULT_PRECISION
    if (allocated(th_buoyancy)) th_buoyancy=0.0_DEFAULT_PRECISION
  end subroutine clear_theta_fluxes  

  !> Computes the heat flux diagnostics for a specific column
  !! @param current_state Current model state
  subroutine compute_theta_flux_for_column(current_state)
    type(model_state_type), target, intent(inout) :: current_state

    integer :: k
    real(kind=DEFAULT_PRECISION), dimension(current_state%local_grid%size(Z_INDEX)) :: thpr
    type(component_field_value_type) :: w_advection_published_value, th_advection_published_value, w_viscosity_published_value, &
         th_diffusion_published_value, w_buoyancy_published_value

    do k=1, current_state%local_grid%size(Z_INDEX)
      thpr(k)=current_state%th%data(k,current_state%column_local_y,current_state%column_local_x)
      if (allocated(current_state%global_grid%configuration%vertical%olthbar)) then
        thpr(k)=thpr(k)-current_state%global_grid%configuration%vertical%olthbar(k)
      end if
    end do
    if (allocated(th_gradient)) then
      w_advection_published_value=get_component_field_value(current_state, "w_advection")
      th_advection_published_value=get_component_field_value(current_state, "th_advection")
    end if
    if (allocated(th_diff)) then
      w_viscosity_published_value=get_component_field_value(current_state, "w_viscosity")
      th_diffusion_published_value=get_component_field_value(current_state, "th_diffusion")
    end if
    if (allocated(th_buoyancy)) then
      w_buoyancy_published_value=get_component_field_value(current_state, "w_buoyancy")
    end if
    do k=2, current_state%local_grid%size(Z_INDEX)-1
      if (allocated(th_flux_values)) th_flux_values(k)=th_flux_values(k)+&
           current_state%w%data(k,current_state%column_local_y,current_state%column_local_x)*&
           current_state%w%data(k,current_state%column_local_y,current_state%column_local_x)*0.5*(thpr(k)+thpr(k+1))
      if (allocated(th_tendency)) th_tendency(k)=th_tendency(k)+0.5*&
           (current_state%w%data(k,current_state%column_local_y,current_state%column_local_x)+&
           current_state%w%data(k-1,current_state%column_local_y,current_state%column_local_x))*&
           current_state%sth%data(k,current_state%column_local_y,current_state%column_local_x)+0.5*&
           thpr(k)*(current_state%sw%data(k,current_state%column_local_y,current_state%column_local_x)+&
           current_state%sw%data(k-1,current_state%column_local_y,current_state%column_local_x))
      if (allocated(th_gradient)) then
        th_gradient(k)=th_gradient(k)+thpr(k)*0.5*(w_advection_published_value%real_1d_array(k)+&
             w_advection_published_value%real_1d_array(k-1))+0.5*&
             (current_state%w%data(k,current_state%column_local_y,current_state%column_local_x)+&
             current_state%w%data(k-1,current_state%column_local_y,current_state%column_local_x))*&
             th_advection_published_value%real_1d_array(k)
      end if
      if (allocated(th_diff)) then
        th_diff(k)=th_diff(k)+thpr(k)*0.5*(w_viscosity_published_value%real_1d_array(k)+&
             w_viscosity_published_value%real_1d_array(k-1))+&
             0.5*(current_state%w%data(k,current_state%column_local_y,current_state%column_local_x)+&
             current_state%w%data(k-1,current_state%column_local_y,current_state%column_local_x))*&
             th_diffusion_published_value%real_1d_array(k)
      end if
      if (allocated(th_buoyancy)) then
        th_buoyancy(k)=th_buoyancy(k)+thpr(k)*0.5*(w_buoyancy_published_value%real_1d_array(k)+&
             w_buoyancy_published_value%real_1d_array(k-1))
      end if
    end do
    if (allocated(w_advection_published_value%real_1d_array)) deallocate(w_advection_published_value%real_1d_array)
    if (allocated(th_advection_published_value%real_1d_array)) deallocate(th_advection_published_value%real_1d_array)
    if (allocated(w_viscosity_published_value%real_1d_array)) deallocate(w_viscosity_published_value%real_1d_array)
    if (allocated(th_diffusion_published_value%real_1d_array)) deallocate(th_diffusion_published_value%real_1d_array)
    if (allocated(w_buoyancy_published_value%real_1d_array)) deallocate(w_buoyancy_published_value%real_1d_array)
  end subroutine compute_theta_flux_for_column

  !> Determines whether a specific published field is a heat flux field
  !! @param name The name of the field to check
  !! @returns Whether the field name is a heat flux field
  logical function is_field_heat_flux(name)
    character(len=*), intent(in) :: name

    is_field_heat_flux=c_contains(heat_flux_fields, name)
  end function is_field_heat_flux
  
  !> Determines whether a specific published field is a TKE budget field
  !! @param name The name of the field to check
  !! @returns Whether the field name is a TKE budget field
  logical function is_field_tke_flux(name)
    character(len=*), intent(in) :: name

    is_field_tke_flux=c_contains(tke_fields, name)
  end function is_field_tke_flux

  !> Determines whether a specific published field is a q flux field
  !! @param name The name of the field to check
  !! @returns Whether the field name is a q flux field
  logical function is_field_q_flux(name)
    character(len=*), intent(in) :: name

    is_field_q_flux=c_contains(q_flux_fields, name)
  end function is_field_q_flux

  !> Determines whether a specific published field is a uw or uv field
  !! @param name The name of the field to check
  !! @returns Whether the field name is a uw or uv field
  logical function is_field_uw_vw(name)
    character(len=*), intent(in) :: name

    is_field_uw_vw=c_contains(uw_vw_fields, name)
  end function is_field_uw_vw  

  !> Determines whether a specific published field is a uu, vv or ww field
  !! @param name The name of the field to check
  !! @returns Whether the field name is a uu, vv or ww field
  logical function is_field_prognostic_budget(name)
    character(len=*), intent(in) :: name

    is_field_prognostic_budget=c_contains(prognostic_budget_fields, name)
  end function is_field_prognostic_budget

  !> Determines whether a specific published field is a thetal field
  !! @param name The name of the field to check
  !! @returns Whether the field name is a thetal field
  logical function is_field_thetal(name)
    character(len=*), intent(in) :: name

    is_field_thetal=c_contains(thetal_fields, name)
  end function is_field_thetal

  !> Determines whether a specific published field is a mse field
  !! @param name The name of the field to check
  !! @returns Whether the field name is a mse field
  logical function is_field_mse(name)
    character(len=*), intent(in) :: name

    is_field_mse=c_contains(mse_fields, name)
  end function is_field_mse

  !> Determines whether a specific published field is a mse field
  !! @param name The name of the field to check
  !! @returns Whether the field name is a mse field
  logical function is_field_qt(name)
    character(len=*), intent(in) :: name

    is_field_qt=c_contains(qt_fields, name)
  end function is_field_qt

  !> Determines whether a specific published field is a scalar field
  !! @param name The name of the field to check
  !! @returns Whether the field name is a scalar field
  logical function is_field_scalar(name)
    character(len=*), intent(in) :: name

    is_field_scalar=c_contains(scalar_fields, name)
  end function is_field_scalar
  
  !> Field value retrieval callback, this returns the value of a specific published field
  !! @param current_state Current model state
  !! @param name The name of the field to retrieve the value for
  !! @param field_value Populated with the value of the field
  subroutine field_value_retrieval_callback(current_state, name, field_value)
    type(model_state_type), target, intent(inout) :: current_state
    character(len=*), intent(in) :: name
    type(component_field_value_type), intent(out) :: field_value
    
    if (name .eq. "heat_flux_transport_local" .and. allocated(th_flux_values)) then
      call set_published_field_value(field_value, real_1d_field=th_flux_values)
    else if (name .eq. "heat_flux_gradient_local" .and. allocated(th_gradient)) then
      call set_published_field_value(field_value, real_1d_field=th_gradient)
    else if (name .eq. "heat_flux_dissipation_local" .and. allocated(th_diff)) then
      call set_published_field_value(field_value, real_1d_field=th_diff)
    else if (name .eq. "heat_flux_buoyancy_local" .and. allocated(th_buoyancy)) then
      call set_published_field_value(field_value, real_1d_field=th_buoyancy)
    else if (name .eq. "heat_flux_tendency_local" .and. allocated(th_tendency)) then
      call set_published_field_value(field_value, real_1d_field=th_tendency)
    else if (name .eq. "q_flux_transport_local" .and. allocated(q_flux_values)) then
      call set_published_field_value(field_value, real_2d_field=q_flux_values)
    else if (name .eq. "q_flux_gradient_local" .and. allocated(q_gradient)) then
      call set_published_field_value(field_value, real_2d_field=q_gradient)
    else if (name .eq. "q_flux_dissipation_local" .and. allocated(q_diff)) then
      call set_published_field_value(field_value, real_2d_field=q_diff)
    else if (name .eq. "q_flux_buoyancy_local" .and. allocated(q_buoyancy)) then
      call set_published_field_value(field_value, real_2d_field=q_buoyancy)
    else if (name .eq. "q_flux_tendency_local" .and. allocated(q_tendency)) then
      call set_published_field_value(field_value, real_2d_field=q_tendency)
    else if (name .eq. "uw_advection_local" .and. allocated(uw_advection)) then
      call set_published_field_value(field_value, real_1d_field=uw_advection)
    else if (name .eq. "vw_advection_local" .and. allocated(vw_advection)) then
      call set_published_field_value(field_value, real_1d_field=vw_advection)
    else if (name .eq. "uw_viscosity_local" .and. allocated(uw_viscosity)) then
      call set_published_field_value(field_value, real_1d_field=uw_viscosity)
    else if (name .eq. "vw_viscosity_local" .and. allocated(vw_viscosity)) then
      call set_published_field_value(field_value, real_1d_field=vw_viscosity)
    else if (name .eq. "uw_buoyancy_local" .and. allocated(uw_buoyancy)) then
      call set_published_field_value(field_value, real_1d_field=uw_buoyancy)
    else if (name .eq. "vw_buoyancy_local" .and. allocated(vw_buoyancy)) then
      call set_published_field_value(field_value, real_1d_field=vw_buoyancy)
    else if (name .eq. "uw_tendency_local" .and. allocated(uw_tendency)) then
      call set_published_field_value(field_value, real_1d_field=uw_tendency)
    else if (name .eq. "vw_tendency_local" .and. allocated(vw_tendency)) then
      call set_published_field_value(field_value, real_1d_field=vw_tendency)
    else if (name .eq. "uw_w_local" .and. allocated(uw_w)) then
      call set_published_field_value(field_value, real_1d_field=uw_w)
    else if (name .eq. "vw_w_local" .and. allocated(vw_w)) then
      call set_published_field_value(field_value, real_1d_field=vw_w)
    else if (name .eq. "resolved_pressure_transport_local" .and. allocated(wp)) then
      call set_published_field_value(field_value, real_1d_field=wp)
    else if (name .eq. "tke_tendency_local" .and. allocated(tend)) then
      call set_published_field_value(field_value, real_1d_field=tend)
    else if (name .eq. "resolved_shear_production_local" .and. allocated(sres)) then
      call set_published_field_value(field_value, real_1d_field=sres)
    else if (name .eq. "resolved_turbulent_transport_local" .and. allocated(wke)) then
      call set_published_field_value(field_value, real_1d_field=wke)
    else if (name .eq. "resolved_buoyant_production_local" .and. allocated(buoy)) then
      call set_published_field_value(field_value, real_1d_field=buoy)
    else if (name .eq. "tu_su_local" .and. allocated(tu_su)) then
      call set_published_field_value(field_value, real_1d_field=tu_su)
    else if (name .eq. "uu_advection_local" .and. allocated(uu_advection)) then
      call set_published_field_value(field_value, real_1d_field=uu_advection)
    else if (name .eq. "uu_viscosity_local" .and. allocated(uu_viscosity)) then
      call set_published_field_value(field_value, real_1d_field=uu_viscosity)
    else if (name .eq. "wu_u_local" .and. allocated(wu_u)) then
      call set_published_field_value(field_value, real_1d_field=wu_u)
    else if (name .eq. "tv_sv_local" .and. allocated(tv_sv)) then
      call set_published_field_value(field_value, real_1d_field=tv_sv)
    else if (name .eq. "vv_advection_local" .and. allocated(vv_advection)) then
      call set_published_field_value(field_value, real_1d_field=vv_advection)
    else if (name .eq. "vv_viscosity_local" .and. allocated(vv_viscosity)) then
      call set_published_field_value(field_value, real_1d_field=vv_viscosity)
    else if (name .eq. "wv_v_local" .and. allocated(wv_v)) then
      call set_published_field_value(field_value, real_1d_field=wv_v)
    else if (name .eq. "tw_sw_local" .and. allocated(tw_sw)) then
      call set_published_field_value(field_value, real_1d_field=tw_sw)
    else if (name .eq. "ww_advection_local" .and. allocated(ww_advection)) then
      call set_published_field_value(field_value, real_1d_field=ww_advection)
    else if (name .eq. "ww_viscosity_local" .and. allocated(ww_viscosity)) then
      call set_published_field_value(field_value, real_1d_field=ww_viscosity)
    else if (name .eq. "ww_buoyancy_local" .and. allocated(ww_buoyancy)) then
      call set_published_field_value(field_value, real_1d_field=ww_buoyancy)
    else if (name .eq. "u_thetal_local" .and. allocated(u_thetal)) then
      call set_published_field_value(field_value, real_1d_field=u_thetal)
    else if (name .eq. "us_thetal_local" .and. allocated(us_thetal)) then
      call set_published_field_value(field_value, real_1d_field=us_thetal)
    else if (name .eq. "u_thetal_advection_local" .and. allocated(u_thetal_advection)) then
      call set_published_field_value(field_value, real_1d_field=u_thetal_advection)
    else if (name .eq. "u_thetal_viscosity_diffusion_local" .and. allocated(u_thetal_viscosity_diffusion)) then
      call set_published_field_value(field_value, real_1d_field=u_thetal_viscosity_diffusion)
    else if (name .eq. "wu_thetal_local" .and. allocated(wu_thetal)) then
      call set_published_field_value(field_value, real_1d_field=wu_thetal)
    else if (name .eq. "v_thetal_local" .and. allocated(v_thetal)) then
      call set_published_field_value(field_value, real_1d_field=v_thetal)
    else if (name .eq. "vs_thetal_local" .and. allocated(vs_thetal)) then
      call set_published_field_value(field_value, real_1d_field=vs_thetal)
    else if (name .eq. "v_thetal_advection_local" .and. allocated(v_thetal_advection)) then
      call set_published_field_value(field_value, real_1d_field=v_thetal_advection)
    else if (name .eq. "v_thetal_viscosity_diffusion_local" .and. allocated(v_thetal_viscosity_diffusion)) then
      call set_published_field_value(field_value, real_1d_field=v_thetal_viscosity_diffusion)
    else if (name .eq. "wv_thetal_local" .and. allocated(wv_thetal)) then
      call set_published_field_value(field_value, real_1d_field=wv_thetal)
    else if (name .eq. "w_thetal_local" .and. allocated(w_thetal)) then
      call set_published_field_value(field_value, real_1d_field=w_thetal)
    else if (name .eq. "ws_thetal_local" .and. allocated(ws_thetal)) then
      call set_published_field_value(field_value, real_1d_field=ws_thetal)
    else if (name .eq. "w_thetal_advection_local" .and. allocated(w_thetal_advection)) then
      call set_published_field_value(field_value, real_1d_field=w_thetal_advection)
    else if (name .eq. "w_thetal_viscosity_diffusion_local" .and. allocated(w_thetal_viscosity_diffusion)) then
      call set_published_field_value(field_value, real_1d_field=w_thetal_viscosity_diffusion)
    else if (name .eq. "w_thetal_buoyancy_local" .and. allocated(w_thetal_buoyancy)) then
      call set_published_field_value(field_value, real_1d_field=w_thetal_buoyancy)
    else if (name .eq. "ww_thetal_local" .and. allocated(ww_thetal)) then
      call set_published_field_value(field_value, real_1d_field=ww_thetal)
    else if (name .eq. "thetal_thetal_local" .and. allocated(thetal_thetal)) then
      call set_published_field_value(field_value, real_1d_field=thetal_thetal)
    else if (name .eq. "sthetal_thetal_local" .and. allocated(sthetal_thetal)) then
      call set_published_field_value(field_value, real_1d_field=sthetal_thetal)
    else if (name .eq. "thetal_thetal_advection_local" .and. allocated(thetal_thetal_advection)) then
      call set_published_field_value(field_value, real_1d_field=thetal_thetal_advection)
    else if (name .eq. "thetal_thetal_diffusion_local" .and. allocated(thetal_thetal_diffusion)) then
      call set_published_field_value(field_value, real_1d_field=thetal_thetal_diffusion)
    else if (name .eq. "wthetal_thetal_local" .and. allocated(wthetal_thetal)) then
      call set_published_field_value(field_value, real_1d_field=wthetal_thetal)
    else if (name .eq. "u_mse_local" .and. allocated(u_mse)) then
      call set_published_field_value(field_value, real_1d_field=u_mse)
    else if (name .eq. "us_mse_local" .and. allocated(us_mse)) then
      call set_published_field_value(field_value, real_1d_field=us_mse)
    else if (name .eq. "u_mse_advection_local" .and. allocated(u_mse_advection)) then
      call set_published_field_value(field_value, real_1d_field=u_mse_advection)
    else if (name .eq. "u_mse_viscosity_diffusion_local" .and. allocated(u_mse_viscosity_diffusion)) then
      call set_published_field_value(field_value, real_1d_field=u_mse_viscosity_diffusion)
    else if (name .eq. "wu_mse_local" .and. allocated(wu_mse)) then
      call set_published_field_value(field_value, real_1d_field=wu_mse)
    else if (name .eq. "v_mse_local" .and. allocated(v_mse)) then
      call set_published_field_value(field_value, real_1d_field=v_mse)
    else if (name .eq. "vs_mse_local" .and. allocated(vs_mse)) then
      call set_published_field_value(field_value, real_1d_field=vs_mse)
    else if (name .eq. "v_mse_advection_local" .and. allocated(v_mse_advection)) then
      call set_published_field_value(field_value, real_1d_field=v_mse_advection)
    else if (name .eq. "v_mse_viscosity_diffusion_local" .and. allocated(v_mse_viscosity_diffusion)) then
      call set_published_field_value(field_value, real_1d_field=v_mse_viscosity_diffusion)
    else if (name .eq. "wv_mse_local" .and. allocated(wv_mse)) then
      call set_published_field_value(field_value, real_1d_field=wv_mse)
    else if (name .eq. "w_mse_local" .and. allocated(w_mse)) then
      call set_published_field_value(field_value, real_1d_field=w_mse)
    else if (name .eq. "ws_mse_local" .and. allocated(ws_mse)) then
      call set_published_field_value(field_value, real_1d_field=ws_mse)
    else if (name .eq. "w_mse_advection_local" .and. allocated(w_mse_advection)) then
      call set_published_field_value(field_value, real_1d_field=w_mse_advection)
    else if (name .eq. "w_mse_viscosity_diffusion_local" .and. allocated(w_mse_viscosity_diffusion)) then
      call set_published_field_value(field_value, real_1d_field=w_mse_viscosity_diffusion)
    else if (name .eq. "w_mse_buoyancy_local" .and. allocated(w_mse_buoyancy)) then
      call set_published_field_value(field_value, real_1d_field=w_mse_buoyancy)
    else if (name .eq. "ww_mse_local" .and. allocated(ww_mse)) then
      call set_published_field_value(field_value, real_1d_field=ww_mse)
    else if (name .eq. "mse_mse_local" .and. allocated(mse_mse)) then
      call set_published_field_value(field_value, real_1d_field=mse_mse)
    else if (name .eq. "smse_mse_local" .and. allocated(smse_mse)) then
      call set_published_field_value(field_value, real_1d_field=smse_mse)
    else if (name .eq. "mse_mse_advection_local" .and. allocated(mse_mse_advection)) then
      call set_published_field_value(field_value, real_1d_field=mse_mse_advection)
    else if (name .eq. "mse_mse_diffusion_local" .and. allocated(mse_mse_diffusion)) then
      call set_published_field_value(field_value, real_1d_field=mse_mse_diffusion)
    else if (name .eq. "wmse_mse_local" .and. allocated(wmse_mse)) then
      call set_published_field_value(field_value, real_1d_field=wmse_mse)
    else if (name .eq. "us_qt_local" .and. allocated(us_qt)) then
      call set_published_field_value(field_value, real_1d_field=us_qt)
    else if (name .eq. "u_qt_advection_local" .and. allocated(u_qt_advection)) then
      call set_published_field_value(field_value, real_1d_field=u_qt_advection)
    else if (name .eq. "u_qt_viscosity_diffusion_local" .and. allocated(u_qt_viscosity_diffusion)) then
      call set_published_field_value(field_value, real_1d_field=u_qt_viscosity_diffusion)
    else if (name .eq. "wu_qt_local" .and. allocated(wu_qt)) then
      call set_published_field_value(field_value, real_1d_field=wu_qt)
    else if (name .eq. "vs_qt_local" .and. allocated(vs_qt)) then
      call set_published_field_value(field_value, real_1d_field=vs_qt)
    else if (name .eq. "v_qt_advection_local" .and. allocated(v_qt_advection)) then
      call set_published_field_value(field_value, real_1d_field=v_qt_advection)
    else if (name .eq. "v_qt_viscosity_diffusion_local" .and. allocated(v_qt_viscosity_diffusion)) then
      call set_published_field_value(field_value, real_1d_field=v_qt_viscosity_diffusion)
    else if (name .eq. "wv_qt_local" .and. allocated(wv_qt)) then
      call set_published_field_value(field_value, real_1d_field=wv_qt)
    else if (name .eq. "w_qt_local" .and. allocated(w_qt)) then
      call set_published_field_value(field_value, real_1d_field=w_qt)
    else if (name .eq. "ws_qt_local" .and. allocated(ws_qt)) then
      call set_published_field_value(field_value, real_1d_field=ws_qt)
    else if (name .eq. "w_qt_advection_local" .and. allocated(w_qt_advection)) then
      call set_published_field_value(field_value, real_1d_field=w_qt_advection)
    else if (name .eq. "w_qt_viscosity_diffusion_local" .and. allocated(w_qt_viscosity_diffusion)) then
      call set_published_field_value(field_value, real_1d_field=w_qt_viscosity_diffusion)
    else if (name .eq. "w_qt_buoyancy_local" .and. allocated(w_qt_buoyancy)) then
      call set_published_field_value(field_value, real_1d_field=w_qt_buoyancy)
    else if (name .eq. "ww_qt_local" .and. allocated(ww_qt)) then
      call set_published_field_value(field_value, real_1d_field=ww_qt)
    else if (name .eq. "qt_qt_local" .and. allocated(qt_qt)) then
      call set_published_field_value(field_value, real_1d_field=qt_qt)
    else if (name .eq. "sqt_qt_local" .and. allocated(sqt_qt)) then
      call set_published_field_value(field_value, real_1d_field=sqt_qt)
    else if (name .eq. "qt_qt_advection_local" .and. allocated(qt_qt_advection)) then
      call set_published_field_value(field_value, real_1d_field=qt_qt_advection)
    else if (name .eq. "qt_qt_diffusion_local" .and. allocated(qt_qt_diffusion)) then
      call set_published_field_value(field_value, real_1d_field=qt_qt_diffusion)
    else if (name .eq. "wqt_qt_local" .and. allocated(wqt_qt)) then
      call set_published_field_value(field_value, real_1d_field=wqt_qt)
    else if (name .eq. "mflux_local") then
      field_value%scalar_real=mflux
    end if
  end subroutine field_value_retrieval_callback

  !> Sets the published field value from the temporary diagnostic values held by this component.
  !! @param field_value Populated with the value of the field
  !! @param real_1d_field Optional one dimensional real of values to publish
  !! @param real_2d_field Optional two dimensional real of values to publish
  subroutine set_published_field_value(field_value, real_1d_field, real_2d_field)
    type(component_field_value_type), intent(inout) :: field_value
    real(kind=DEFAULT_PRECISION), dimension(:), optional :: real_1d_field
    real(kind=DEFAULT_PRECISION), dimension(:,:), optional :: real_2d_field

    if (present(real_1d_field)) then
      allocate(field_value%real_1d_array(size(real_1d_field)), source=real_1d_field)
    else if (present(real_2d_field)) then
      allocate(field_value%real_2d_array(size(real_2d_field, 1), size(real_2d_field, 2)), source=real_2d_field)      
    end if    
  end subroutine set_published_field_value

  !> Sets the published value enabled state in the provided collection map
  !! @param collection The map to set this in
  !! @param field_name The name of the published field
  !! @param enabled_state Whether the field is enabled which is stored
  subroutine set_published_field_enabled_state(collection, field_name, enabled_state)
    type(hashmap_type), intent(inout) :: collection
    character(len=*), intent(in) :: field_name
    logical, intent(in) :: enabled_state
    
    call c_put_logical(collection, field_name, enabled_state)
  end subroutine set_published_field_enabled_state

  !> Retrieves whether a published field is enabled or not
  !! @param collection The map to look up the published field in
  !! @param field_name The name of the field to look up
  !! @returns Whether this published field is enabled or not
  logical function get_published_field_enabled_state(collection, field_name)
    type(hashmap_type), intent(inout) :: collection
    character(len=*), intent(in) :: field_name

    get_published_field_enabled_state=c_get_logical(collection, field_name)
  end function get_published_field_enabled_state
end module flux_budget_mod
