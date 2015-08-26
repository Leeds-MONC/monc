!> Specific theta advection, which involves the vertical advection of reference state and advection of mean baroclinicity
module thadvection_mod
  use datadefn_mod, only : DEFAULT_PRECISION
  use monc_component_mod, only : component_descriptor_type
  use state_mod, only : model_state_type
  use grids_mod, only : Z_INDEX
  use science_constants_mod, only : G
  use logging_mod, only : LOG_ERROR, log_master_log
  use optionsdatabase_mod, only : options_get_real, options_get_logical
implicit none

#ifndef TEST_MODE
  private
#endif

  logical :: baroclinicity_use_geostrophic_shear
  real(kind=DEFAULT_PRECISION) :: fcoriol, fcoriol_over_G, rate_change_geostrophic_wind_x, rate_change_geostrophic_wind_y, &
       multiplicative_factor_x, multiplicative_factor_y

  public thadvection_get_descriptor
contains

  !> Provides a description of this component for the core to register
  !! @returns The descriptor containing registration information for this component
  type(component_descriptor_type) function thadvection_get_descriptor()
    thadvection_get_descriptor%name="th_advection"
    thadvection_get_descriptor%version=0.1
    thadvection_get_descriptor%initialisation=>initialisation_callback
    thadvection_get_descriptor%timestep=>timestep_callback
  end function thadvection_get_descriptor

  !> Initialisation callback to set up the variables and data needed by the component
  !! @param current_state The current model state
  subroutine initialisation_callback(current_state)
    type(model_state_type), target, intent(inout) :: current_state

    baroclinicity_use_geostrophic_shear=options_get_logical(current_state%options_database, "baroclinicity_use_geostrophic_shear")
    fcoriol=options_get_real(current_state%options_database, "fcoriol")
    rate_change_geostrophic_wind_x=options_get_real(current_state%options_database, "rate_change_geostrophic_wind_x")
    rate_change_geostrophic_wind_y=options_get_real(current_state%options_database, "rate_change_geostrophic_wind_y")
    fcoriol_over_G = fcoriol/G
    multiplicative_factor_x=rate_change_geostrophic_wind_x*current_state%thref0*fcoriol_over_G
    multiplicative_factor_y=rate_change_geostrophic_wind_y*current_state%thref0*fcoriol_over_G
  end subroutine initialisation_callback  

  !> Timestep callback, will call the two separate procedures to do their advection if needed
  !! @param current_state The current model state
  subroutine timestep_callback(current_state)
    type(model_state_type), target, intent(inout) :: current_state

    call vertical_advection_of_reference_state(current_state, current_state%column_local_y, current_state%column_local_x)
    call advection_of_mean_baroclinicity(current_state, current_state%column_local_y, current_state%column_local_x)
  end subroutine timestep_callback

  !> Vertical advection of the reference state.  It doesn't seem consistent to do the advection in this way if
  !! TVD advection of the deviation from the reference state has been selected. Separate vertical advection of the reference
  !! state was introduced to improve energy conservation when carrying out idealized gravity wave simulations in a deep, dry
  !! isothermal layer, for which the difference in potential temp between top and bottom was of order 100K. In less extreme cases
  !! the benefits are unlikely to be significant and with TVD advection energy conservation has been compromised so the best
  !! way forward might be to recombine the reference state into l_th
  !! @param current_state The current model state
  !! @param local_y The local Y of the column
  !! @param local_x The local X of the column
  subroutine vertical_advection_of_reference_state(current_state, local_y, local_x)
    type(model_state_type), target, intent(inout) :: current_state
    integer, intent(in) :: local_y, local_x

    integer :: k
    real(kind=DEFAULT_PRECISION) :: sctmp1, sctmp2

    if (current_state%use_anelastic_equations) then
      ! This code only needs to be executed if anelastic, otherwise THREF is constant and the increment calculated here is zero
      do k=2, current_state%local_grid%size(Z_INDEX)
        sctmp1=current_state%global_grid%configuration%vertical%tzc1(k)*2.0_DEFAULT_PRECISION*&
             current_state%global_grid%configuration%vertical%dthref(k-1)
        sctmp2=current_state%global_grid%configuration%vertical%tzc2(k)*2.0_DEFAULT_PRECISION*&
             current_state%global_grid%configuration%vertical%dthref(k)
        current_state%sth%data(k, local_y, local_x)=current_state%sth%data(k, local_y, local_x)-(sctmp1*&
             current_state%w%data(k-1, local_y, local_x) + sctmp2*current_state%w%data(k, local_y, local_x))
      end do
    end if
  end subroutine vertical_advection_of_reference_state

  !> Performs advection of the mean baroclinicity if appropriate
  !! @param current_state The current model state
  !! @param local_y The local Y of the column
  !! @param local_x The local X of the column
  subroutine advection_of_mean_baroclinicity(current_state, local_y, local_x)
    type(model_state_type), target, intent(inout) :: current_state
    integer, intent(in) :: local_y, local_x

    integer :: k

    if (baroclinicity_use_geostrophic_shear) then
        if (current_state%passive_q) then
            if (current_state%use_anelastic_equations) then
                do k=2, current_state%local_grid%size(Z_INDEX)
                  current_state%sth%data(k, local_y, local_x)=current_state%sth%data(k, local_y, local_x)+&
                       current_state%global_grid%configuration%vertical%thref(k)*fcoriol_over_G*&
                       ((current_state%v%data(k, local_y, local_x) + current_state%vgal) * rate_change_geostrophic_wind_x-&
                       (current_state%u%data(k, local_y, local_x) + current_state%ugal) * rate_change_geostrophic_wind_y)
                end do
              else
                do k=2, current_state%local_grid%size(Z_INDEX)
                  current_state%sth%data(k, local_y, local_x)=current_state%sth%data(k, local_y, local_x)+&
                       ((current_state%v%data(k, local_y, local_x) + current_state%vgal) * multiplicative_factor_x-&
                       (current_state%u%data(k, local_y, local_x) + current_state%ugal) * multiplicative_factor_y)                    
                end do
              end if
        else
          call log_master_log(LOG_ERROR, "The combination if baroclinicity and active q is not yet allowed")
        end if
      end if
  end subroutine advection_of_mean_baroclinicity  
end module thadvection_mod
