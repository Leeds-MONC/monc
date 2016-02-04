module profile_diagnostics_mod
  use monc_component_mod, only : COMPONENT_SCALAR_FIELD_TYPE, COMPONENT_DOUBLE_DATA_TYPE, &
       COMPONENT_ARRAY_FIELD_TYPE, COMPONENT_INTEGER_DATA_TYPE, &
       component_descriptor_type, component_field_value_type, component_field_information_type
  use grids_mod, only : Z_INDEX, Y_INDEX, X_INDEX
  use optionsdatabase_mod, only : options_get_real
  use state_mod, only : model_state_type
  use datadefn_mod, only : DEFAULT_PRECISION
  use q_indices_mod, only: q_indices_add

  implicit none

#ifndef TEST_MODE
  private
#endif

  integer :: total_points, iqv, iql
  real(kind=DEFAULT_PRECISION), dimension(:), allocatable :: tempfac, theta_tot, theta_mean_loc, &
       qv_tot, qv_mean_loc, ql_tot, ql_mean_loc
  real(kind=DEFAULT_PRECISION) :: qlcrit

  public profile_diagnostics_get_descriptor

contains

  !> Provides the component descriptor for the core to register
  !! @returns The descriptor describing this component
  type(component_descriptor_type) function profile_diagnostics_get_descriptor()
    profile_diagnostics_get_descriptor%name="profile_diagnostics"
    profile_diagnostics_get_descriptor%version=0.1

    profile_diagnostics_get_descriptor%initialisation=>initialisation_callback
    profile_diagnostics_get_descriptor%timestep=>timestep_callback

    profile_diagnostics_get_descriptor%field_value_retrieval=>field_value_retrieval_callback
    profile_diagnostics_get_descriptor%field_information_retrieval=>field_information_retrieval_callback
    allocate(profile_diagnostics_get_descriptor%published_fields(3))

    profile_diagnostics_get_descriptor%published_fields(1)="theta_mean_local"
    profile_diagnostics_get_descriptor%published_fields(2)="vapour_mmr_mean_local"
    profile_diagnostics_get_descriptor%published_fields(3)="liquid_mmr_mean_local"

  end function profile_diagnostics_get_descriptor

  subroutine initialisation_callback(current_state)
    type(model_state_type), target, intent(inout) :: current_state

    integer :: k
    
    iqv=q_indices_add('vapour', 'simplecloud')
    iql=q_indices_add('cloud liquid mass', 'simplecloud')

    qlcrit=options_get_real(current_state%options_database, "qlcrit")
    total_points=(current_state%local_grid%size(Y_INDEX) * current_state%local_grid%size(X_INDEX))

    allocate(theta_tot(current_state%local_grid%size(Z_INDEX)), theta_mean_loc(current_state%local_grid%size(Z_INDEX)) &
         , qv_tot(current_state%local_grid%size(Z_INDEX)), qv_mean_loc(current_state%local_grid%size(Z_INDEX)) &
         , ql_tot(current_state%local_grid%size(Z_INDEX)), ql_mean_loc(current_state%local_grid%size(Z_INDEX)))
   
  end subroutine initialisation_callback  

  subroutine timestep_callback(current_state)
    type(model_state_type), target, intent(inout) :: current_state

    integer :: k, i
    real(kind=DEFAULT_PRECISION) :: cltop_col, clbas_col

    if (current_state%first_timestep_column) then
      theta_tot(:)=0.0_DEFAULT_PRECISION
      qv_tot(:)=0.0_DEFAULT_PRECISION
      ql_tot(:)=0.0_DEFAULT_PRECISION
    end if
    if (.not. current_state%halo_column) then
       do k=1, current_state%local_grid%size(Z_INDEX)
          theta_tot(k) = theta_tot(k) + (current_state%th%data(k,current_state%column_local_y,current_state%column_local_x) &
               + current_state%global_grid%configuration%vertical%thref(k))
          qv_tot(k) = qv_tot(k) + (current_state%q(iqv)%data(k,current_state%column_local_y,current_state%column_local_x)) 
          ql_tot(k) = ql_tot(k) + (current_state%q(iql)%data(k,current_state%column_local_y,current_state%column_local_x))
       enddo
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
    field_information%enabled=.true.

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

    if (name .eq. "theta_mean_local") then
       allocate(field_value%real_1d_array(current_state%local_grid%size(Z_INDEX)))
       do k = 1, current_state%local_grid%size(Z_INDEX)
          field_value%real_1d_array(k)=theta_tot(k)/total_points
       enddo
      ! deallocate(field_value%real_1d_array)
    else if (name .eq. "vapour_mmr_mean_local") then
       allocate(field_value%real_1d_array(current_state%local_grid%size(Z_INDEX)))
       do k = 1, current_state%local_grid%size(Z_INDEX)
          field_value%real_1d_array(k)=qv_tot(k)/total_points
       enddo
      ! deallocate(field_value%real_1d_array)
    else if (name .eq. "liquid_mmr_mean_local") then
       allocate(field_value%real_1d_array(current_state%local_grid%size(Z_INDEX)))
       do k = 1, current_state%local_grid%size(Z_INDEX)
          field_value%real_1d_array(k)=ql_tot(k)/total_points
       enddo
      ! deallocate(field_value%real_1d_array)  
    end if
  end subroutine field_value_retrieval_callback
end module profile_diagnostics_mod
