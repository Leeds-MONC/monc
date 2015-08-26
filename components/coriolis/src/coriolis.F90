!> This calculates the coriolis and mean pressure gradient terms which impact su and sv fields
module coriolis_mod
  use datadefn_mod, only : DEFAULT_PRECISION
  use monc_component_mod, only : component_descriptor_type
  use optionsdatabase_mod, only : options_get_logical, options_get_real
  use state_mod, only : model_state_type
  use grids_mod, only : Z_INDEX, Y_INDEX, X_INDEX
  implicit none

#ifndef TEST_MODE
  private
#endif

  logical :: baroclinicity_use_geostrophic_shear
  real(kind=DEFAULT_PRECISION), dimension(:), allocatable :: geostrophic_wind_x, geostrophic_wind_y
  real(kind=DEFAULT_PRECISION) :: fcoriol

  public coriolis_get_descriptor

  contains

  !> Provides the descriptor back to the caller and is used in component registration
  !! @returns The termination check component descriptor
  type(component_descriptor_type) function coriolis_get_descriptor()
    coriolis_get_descriptor%name="coriolis"
    coriolis_get_descriptor%version=0.1
    coriolis_get_descriptor%initialisation=>initialisation_callback
    coriolis_get_descriptor%timestep=>timestep_callback
  end function coriolis_get_descriptor

  !> Initialisation call back which will read in the coriolis configuration and set up the geostrophic winds
  !! @param current_state The current model state
  subroutine initialisation_callback(current_state)
    type(model_state_type), target, intent(inout) :: current_state

    integer :: k
    
    baroclinicity_use_geostrophic_shear=options_get_logical(current_state%options_database, "baroclinicity_use_geostrophic_shear")
    fcoriol=options_get_real(current_state%options_database, "fcoriol")
    current_state%geostrophic_wind_rate_of_change_in_x=options_get_real(current_state%options_database, &
         "geostrophic_wind_rate_of_change_in_x")
    current_state%geostrophic_wind_rate_of_change_in_y=options_get_real(current_state%options_database, &
         "geostrophic_wind_rate_of_change_in_y")
    current_state%surface_geostrophic_wind_x=options_get_real(current_state%options_database, "surface_geostrophic_wind_x")
    current_state%surface_geostrophic_wind_y=options_get_real(current_state%options_database, "surface_geostrophic_wind_y")

    allocate(geostrophic_wind_x(current_state%local_grid%size(Z_INDEX)), &
         geostrophic_wind_y(current_state%local_grid%size(Z_INDEX)))

    do k=1,current_state%local_grid%size(Z_INDEX)
      geostrophic_wind_x(k)=current_state%surface_geostrophic_wind_x
      geostrophic_wind_y(k)=current_state%surface_geostrophic_wind_y
      if (baroclinicity_use_geostrophic_shear) then
        geostrophic_wind_x(k)=geostrophic_wind_x(k)+current_state%geostrophic_wind_rate_of_change_in_x*&
             current_state%global_grid%configuration%vertical%zn(k)
        geostrophic_wind_y(k)=geostrophic_wind_y(k)+current_state%geostrophic_wind_rate_of_change_in_y*&
             current_state%global_grid%configuration%vertical%zn(k)
      end if
    end do
  end subroutine initialisation_callback

  !> For each none halo cell this will calculate the coriolis terms for su and sv fields
  !! @param current_state The current model state
  subroutine timestep_callback(current_state)
    type(model_state_type), target, intent(inout) :: current_state

    integer :: k

    if (current_state%halo_column) then
      if (.not. ((current_state%column_local_y == current_state%local_grid%halo_size(Y_INDEX) .and. &
           current_state%column_local_x .le. current_state%local_grid%local_domain_end_index(X_INDEX) .and. &
           current_state%column_local_x .ge. current_state%local_grid%local_domain_start_index(X_INDEX)-1) .or. &
           (current_state%column_local_x == current_state%local_grid%halo_size(X_INDEX) .and. &
           current_state%column_local_y .ge. current_state%local_grid%local_domain_start_index(Y_INDEX) &
           .and. current_state%column_local_y .le. current_state%local_grid%local_domain_end_index(Y_INDEX)) )) return
    end if
    
    do k=2,current_state%local_grid%size(Z_INDEX)
#if defined(U_ACTIVE) && defined(V_ACTIVE)
      current_state%su%data(k, current_state%column_local_y, current_state%column_local_x)=&
           current_state%su%data(k, current_state%column_local_y, current_state%column_local_x)+fcoriol*&
           (0.25_DEFAULT_PRECISION*(current_state%v%data(k, current_state%column_local_y, current_state%column_local_x)+&
           current_state%v%data(k, current_state%column_local_y, current_state%column_local_x+1)+&
           current_state%v%data(k, current_state%column_local_y-1, current_state%column_local_x)+&
           current_state%v%data(k, current_state%column_local_y-1, current_state%column_local_x+1))+current_state%vgal-&
           geostrophic_wind_y(k))

      current_state%sv%data(k, current_state%column_local_y, current_state%column_local_x)=&
           current_state%sv%data(k, current_state%column_local_y, current_state%column_local_x)-fcoriol*&
           (0.25_DEFAULT_PRECISION*(current_state%u%data(k, current_state%column_local_y, current_state%column_local_x)+&
           current_state%u%data(k, current_state%column_local_y, current_state%column_local_x-1)+&
           current_state%u%data(k, current_state%column_local_y+1, current_state%column_local_x)+&
           current_state%u%data(k, current_state%column_local_y+1, current_state%column_local_x-1))+current_state%ugal-&
           geostrophic_wind_x(k))
#endif
    end do
  end subroutine timestep_callback
end module coriolis_mod
