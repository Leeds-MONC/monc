!> Derive 3-D diagnostic fields which are available in current_state 
module diagnostics_3d_mod
  !!
  !! module to derive 3-D fields which are not already stored
  !! in current_state and output as a 3-D diagnostic
  !!
  use monc_component_mod, only : COMPONENT_SCALAR_FIELD_TYPE, COMPONENT_DOUBLE_DATA_TYPE, &
       COMPONENT_ARRAY_FIELD_TYPE, COMPONENT_INTEGER_DATA_TYPE, &
       component_descriptor_type, component_field_value_type, component_field_information_type
  use grids_mod, only : Z_INDEX, Y_INDEX, X_INDEX
  use optionsdatabase_mod, only : options_get_real
  use state_mod, only : model_state_type
  use datadefn_mod, only : DEFAULT_PRECISION
  use q_indices_mod, only: get_q_index, standard_q_names
  use saturation_mod, only: qsaturation
  use science_constants_mod, only : rlvap_over_cp

  implicit none

#ifndef TEST_MODE
  private
#endif

  integer :: iqv, iql, iqr
  real(kind=DEFAULT_PRECISION), dimension(:,:,:), allocatable ::     &
       TdegK,               & ! absolute temperature in kelvin
       theta,               & ! potential temperature in kelvin (th + thref)
       liquid_ice_theta       ! liquid-ice potential temperature in kelvin
  real(kind=DEFAULT_PRECISION), dimension(:), allocatable ::         &    
       total_condensate
  real(kind=DEFAULT_PRECISION) :: qlcrit

  public diagnostics_3d_get_descriptor

contains

  !> Provides the component descriptor for the core to register
  !! @returns The descriptor describing this component
  type(component_descriptor_type) function diagnostics_3d_get_descriptor()
    diagnostics_3d_get_descriptor%name="diagnostics_3d"
    diagnostics_3d_get_descriptor%version=0.1
    diagnostics_3d_get_descriptor%initialisation=>initialisation_callback
    diagnostics_3d_get_descriptor%timestep=>timestep_callback

    diagnostics_3d_get_descriptor%field_value_retrieval=>field_value_retrieval_callback
    diagnostics_3d_get_descriptor%field_information_retrieval=>field_information_retrieval_callback
    allocate(diagnostics_3d_get_descriptor%published_fields(3))

    diagnostics_3d_get_descriptor%published_fields(1)="temperature"
    diagnostics_3d_get_descriptor%published_fields(2)="theta"
    diagnostics_3d_get_descriptor%published_fields(3)="thetali"
    
  end function diagnostics_3d_get_descriptor
  
  subroutine initialisation_callback(current_state)
    type(model_state_type), target, intent(inout) :: current_state

    if (current_state%th%active) then
       allocate(TdegK(current_state%local_grid%size(Z_INDEX),  &
            current_state%local_grid%size(Y_INDEX),            &
            current_state%local_grid%size(X_INDEX)))
       allocate(theta(current_state%local_grid%size(Z_INDEX),  &
            current_state%local_grid%size(Y_INDEX),            &
            current_state%local_grid%size(X_INDEX)))
       allocate(liquid_ice_theta(current_state%local_grid%size(Z_INDEX),  &
            current_state%local_grid%size(Y_INDEX),            &
            current_state%local_grid%size(X_INDEX)))
    endif

    if (.not. current_state%passive_q .and. current_state%number_q_fields .gt. 0) then
       allocate(total_condensate(current_state%local_grid%size(Z_INDEX)))
       iqv=get_q_index(standard_q_names%VAPOUR, 'diagnostics_3d') 
       iql=get_q_index(standard_q_names%CLOUD_LIQUID_MASS, 'diagnostics_3d')
       if (current_state%rain_water_mixing_ratio_index > 0) & 
            iqr = current_state%rain_water_mixing_ratio_index
    endif
    
  end subroutine initialisation_callback
  
  subroutine timestep_callback(current_state)
    type(model_state_type), target, intent(inout) :: current_state

    real(kind=DEFAULT_PRECISION) :: exner, pmb, qv, qs

    integer :: k
    integer :: current_y_index, current_x_index, target_x_index, target_y_index

    if (current_state%halo_column) return
       
    current_y_index=current_state%column_local_y
    current_x_index=current_state%column_local_x
    target_y_index=current_y_index-current_state%local_grid%halo_size(Y_INDEX)
    target_x_index=current_x_index-current_state%local_grid%halo_size(X_INDEX)

    if (current_state%th%active) then
       theta(:,target_y_index, target_x_index) =                            &
            (current_state%th%data(:,current_y_index,current_x_index)       &
            + current_state%global_grid%configuration%vertical%thref(:))
       liquid_ice_theta(:,target_y_index, target_x_index) =                 &
            (current_state%th%data(:,current_y_index,current_x_index)       &
            + current_state%global_grid%configuration%vertical%thref(:))  
            ! test for the qfields
       if (.not. current_state%passive_q .and. &
            current_state%number_q_fields .gt. 0) then
          total_condensate(:) =                                             &
               current_state%q(iql)%data(:,current_y_index,current_x_index) + &
               current_state%q(iqr)%data(:,current_y_index,current_x_index)
          liquid_ice_theta(:,target_y_index, target_x_index) =              &
               liquid_ice_theta(:,target_y_index, target_x_index) -         &
               ( total_condensate(:) * (rlvap_over_cp) )
       endif
       TdegK(:,target_y_index, target_x_index) =                            &
            (current_state%th%data(:,current_y_index,current_x_index)       &
            + current_state%global_grid%configuration%vertical%thref(:))    &
            * current_state%global_grid%configuration%vertical%rprefrcp(:)       
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
    field_information%number_dimensions=3
    field_information%dimension_sizes(1)=current_state%local_grid%size(Z_INDEX)
    field_information%dimension_sizes(2)=current_state%local_grid%size(Y_INDEX)
    field_information%dimension_sizes(3)=current_state%local_grid%size(X_INDEX)
    field_information%data_type=COMPONENT_DOUBLE_DATA_TYPE
    if (name .eq. "temperature" .or. name .eq. "theta" &
         .or. name .eq. "thetali") then
       field_information%enabled=current_state%th%active
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

    if (name .eq. "temperature") then
       allocate(field_value%real_3d_array(current_state%local_grid%size(Z_INDEX),  &
            current_state%local_grid%size(Y_INDEX),                                &
            current_state%local_grid%size(X_INDEX)))
       field_value%real_3d_array(:,:,:) = TdegK(:,:,:)
    elseif (name .eq. "theta") then
       allocate(field_value%real_3d_array(current_state%local_grid%size(Z_INDEX),  &
            current_state%local_grid%size(Y_INDEX),                                &
            current_state%local_grid%size(X_INDEX)))
       field_value%real_3d_array(:,:,:) = theta(:,:,:)
    elseif (name .eq. "thetali") then
       allocate(field_value%real_3d_array(current_state%local_grid%size(Z_INDEX),  &
            current_state%local_grid%size(Y_INDEX),                                &
            current_state%local_grid%size(X_INDEX)))
       field_value%real_3d_array(:,:,:) = liquid_ice_theta(:,:,:)
    end if
    
  end subroutine field_value_retrieval_callback
end module diagnostics_3d_mod
