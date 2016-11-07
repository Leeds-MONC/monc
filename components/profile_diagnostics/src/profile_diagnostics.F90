module profile_diagnostics_mod
  use monc_component_mod, only : COMPONENT_SCALAR_FIELD_TYPE, COMPONENT_DOUBLE_DATA_TYPE, &
       COMPONENT_ARRAY_FIELD_TYPE, COMPONENT_INTEGER_DATA_TYPE, &
       component_descriptor_type, component_field_value_type, component_field_information_type
  use grids_mod, only : Z_INDEX, Y_INDEX, X_INDEX
  use optionsdatabase_mod, only : options_get_real
  use state_mod, only : model_state_type
  use datadefn_mod, only : DEFAULT_PRECISION
  use q_indices_mod, only: get_q_index, standard_q_names

  implicit none

#ifndef TEST_MODE
  private
#endif

  integer :: total_points, iqv, iql
  real(kind=DEFAULT_PRECISION), dimension(:), allocatable ::     &
       tempfac, u_wind_tot, uprime_tot, v_wind_tot, vprime_tot,  &
       ww_tot, theta_tot, qv_tot, ql_tot, w_wind_tot
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
    allocate(profile_diagnostics_get_descriptor%published_fields(9))

    profile_diagnostics_get_descriptor%published_fields(1)="theta_total_local"
    profile_diagnostics_get_descriptor%published_fields(2)="vapour_mmr_total_local"
    profile_diagnostics_get_descriptor%published_fields(3)="liquid_mmr_total_local"
    profile_diagnostics_get_descriptor%published_fields(4)="u_wind_total_local"
    profile_diagnostics_get_descriptor%published_fields(5)="uu_total_local"
    profile_diagnostics_get_descriptor%published_fields(6)="v_wind_total_local"
    profile_diagnostics_get_descriptor%published_fields(7)="vv_total_local"
    profile_diagnostics_get_descriptor%published_fields(8)="ww_total_local"
    profile_diagnostics_get_descriptor%published_fields(9)="w_wind_total_local"

  end function profile_diagnostics_get_descriptor

  subroutine initialisation_callback(current_state)
    type(model_state_type), target, intent(inout) :: current_state

    integer :: k

    total_points=(current_state%local_grid%size(Y_INDEX) * current_state%local_grid%size(X_INDEX))
    
    ! allocate local arrays for the horizontal wind averages
    allocate(u_wind_tot(current_state%local_grid%size(Z_INDEX)) &
         , v_wind_tot(current_state%local_grid%size(Z_INDEX))   &
         , ww_tot(current_state%local_grid%size(Z_INDEX)) &
         , w_wind_tot(current_state%local_grid%size(Z_INDEX)))

    if (allocated(current_state%global_grid%configuration%vertical%olubar)) then
      allocate(uprime_tot(current_state%local_grid%size(Z_INDEX)))
    end if
    
    if (allocated(current_state%global_grid%configuration%vertical%olvbar)) then
      allocate(vprime_tot(current_state%local_grid%size(Z_INDEX)))
    end if    
    
    if (current_state%th%active) then
        allocate(theta_tot(current_state%local_grid%size(Z_INDEX))) 
     endif
        
    if (.not. current_state%passive_q .and. current_state%number_q_fields .gt. 0) then                                     
       iqv=get_q_index(standard_q_names%VAPOUR, 'profile_diags')                                                                
       iql=get_q_index(standard_q_names%CLOUD_LIQUID_MASS, 'profile_diags')                                                     
       qlcrit=options_get_real(current_state%options_database, "qlcrit")                                                         
       allocate(qv_tot(current_state%local_grid%size(Z_INDEX))  &
         , ql_tot(current_state%local_grid%size(Z_INDEX)))
    endif   
  end subroutine initialisation_callback  

  subroutine timestep_callback(current_state)
    type(model_state_type), target, intent(inout) :: current_state

    integer :: k, i
    real(kind=DEFAULT_PRECISION) :: cltop_col, clbas_col

    if (current_state%first_timestep_column) then
       u_wind_tot(:) = 0.0
       if (allocated(uprime_tot)) uprime_tot(:) = 0.0
       v_wind_tot(:) = 0.0
       if (allocated(vprime_tot)) vprime_tot(:) = 0.0
       w_wind_tot(:) = 0.0
       ww_tot(:) = 0.0

       if (current_state%th%active) then 
          theta_tot(:)=0.0_DEFAULT_PRECISION
       endif
       if (.not. current_state%passive_q .and. &
            current_state%number_q_fields .gt. 0) then 
          qv_tot(:)=0.0_DEFAULT_PRECISION
          ql_tot(:)=0.0_DEFAULT_PRECISION
       endif
    end if
    if (.not. current_state%halo_column) then
       ! work out the sum of u and v wind over local domain
       do k=1, current_state%local_grid%size(Z_INDEX)
          u_wind_tot(k) = u_wind_tot(k) + & 
               (current_state%u%data(k,current_state%column_local_y,current_state%column_local_x)  &
                + current_state%ugal)
          if (allocated(uprime_tot)) then
            uprime_tot(k) = uprime_tot(k) + &
                 ((current_state%u%data(k,current_state%column_local_y,current_state%column_local_x) &
                 - (current_state%global_grid%configuration%vertical%olubar(k) - current_state%ugal))**2.)
          end if
          v_wind_tot(k) = v_wind_tot(k) + & 
               (current_state%v%data(k,current_state%column_local_y,current_state%column_local_x)  &
               + current_state%vgal)
          if (allocated(vprime_tot)) then
            vprime_tot(k) = vprime_tot(k) + &
                 ((current_state%v%data(k,current_state%column_local_y,current_state%column_local_x) &
                 - (current_state%global_grid%configuration%vertical%olvbar(k) - current_state%vgal))**2.)
          end if
          ww_tot(k) = ww_tot(k) + &
               (current_state%w%data(k,current_state%column_local_y,current_state%column_local_x)**2.)
          w_wind_tot(k) = w_wind_tot(k) + & 
               (current_state%w%data(k,current_state%column_local_y,current_state%column_local_x))
       enddo
       if (current_state%th%active) then
          do k=1, current_state%local_grid%size(Z_INDEX)
             theta_tot(k) = theta_tot(k) + & 
                  (current_state%th%data(k,current_state%column_local_y,current_state%column_local_x) &
                  + current_state%global_grid%configuration%vertical%thref(k))
          enddo
       endif
       if (.not. current_state%passive_q .and. current_state%number_q_fields .gt. 0) then
          do k=1, current_state%local_grid%size(Z_INDEX)
             qv_tot(k) = qv_tot(k) + (current_state%q(iqv)%data(k,current_state%column_local_y,current_state%column_local_x))   
             ql_tot(k) = ql_tot(k) + (current_state%q(iql)%data(k,current_state%column_local_y,current_state%column_local_x))
          enddo
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
    if (name .eq. "theta_total_local") then
      field_information%enabled=current_state%th%active
    else if (name .eq. "vapour_mmr_total_local" .or. name .eq. "liquid_mmr_total_local") then
      field_information%enabled=.not. current_state%passive_q .and. current_state%number_q_fields .gt. 0
    else if (name .eq. "uu_total_local") then
      field_information%enabled=allocated(uprime_tot)
    else if (name .eq. "vv_total_local") then
      field_information%enabled=allocated(vprime_tot)
    else 
      field_information%enabled=.true.
    end if
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

    if (name .eq. "u_wind_total_local") then
       allocate(field_value%real_1d_array(current_state%local_grid%size(Z_INDEX)))
       do k = 1, current_state%local_grid%size(Z_INDEX)
          field_value%real_1d_array(k)=u_wind_tot(k)
       enddo
    else if (name .eq. "uu_total_local") then
       allocate(field_value%real_1d_array(current_state%local_grid%size(Z_INDEX)))
       do k = 1, current_state%local_grid%size(Z_INDEX)
          field_value%real_1d_array(k)=uprime_tot(k)
       enddo
    else if (name .eq. "v_wind_total_local") then
       allocate(field_value%real_1d_array(current_state%local_grid%size(Z_INDEX)))
       do k = 1, current_state%local_grid%size(Z_INDEX)
          field_value%real_1d_array(k)=v_wind_tot(k)
       enddo
    else if (name .eq. "vv_total_local") then
       allocate(field_value%real_1d_array(current_state%local_grid%size(Z_INDEX)))
       do k = 1, current_state%local_grid%size(Z_INDEX)
          field_value%real_1d_array(k)=vprime_tot(k)
       enddo
    else if (name .eq. "ww_total_local") then
       allocate(field_value%real_1d_array(current_state%local_grid%size(Z_INDEX)))
       do k = 1, current_state%local_grid%size(Z_INDEX)
          field_value%real_1d_array(k)=ww_tot(k)
       enddo   
    else if (name .eq. "theta_total_local") then
       allocate(field_value%real_1d_array(current_state%local_grid%size(Z_INDEX)))
       do k = 1, current_state%local_grid%size(Z_INDEX)
          field_value%real_1d_array(k)=theta_tot(k)
       enddo
    else if (name .eq. "vapour_mmr_total_local") then
       allocate(field_value%real_1d_array(current_state%local_grid%size(Z_INDEX)))
       do k = 1, current_state%local_grid%size(Z_INDEX)
          field_value%real_1d_array(k)=qv_tot(k)
       enddo
    else if (name .eq. "liquid_mmr_total_local") then
       allocate(field_value%real_1d_array(current_state%local_grid%size(Z_INDEX)))
       do k = 1, current_state%local_grid%size(Z_INDEX)
          field_value%real_1d_array(k)=ql_tot(k)
       enddo
    else if (name .eq. "w_wind_total_local") then
       allocate(field_value%real_1d_array(current_state%local_grid%size(Z_INDEX)))
       do k = 1, current_state%local_grid%size(Z_INDEX)
          field_value%real_1d_array(k)=w_wind_tot(k)
       enddo
    end if
  end subroutine field_value_retrieval_callback
end module profile_diagnostics_mod
