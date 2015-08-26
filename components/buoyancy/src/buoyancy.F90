!> Calculates buoyancy terms for the SW field
module buoyancy_mod
  use datadefn_mod, only : DEFAULT_PRECISION
  use monc_component_mod, only : COMPONENT_ARRAY_FIELD_TYPE, COMPONENT_DOUBLE_DATA_TYPE, component_descriptor_type, &
       component_field_value_type, component_field_information_type
  use registry_mod, only : is_component_enabled
  use optionsdatabase_mod, only : options_has_key, options_get_real_array
  use state_mod, only : FORWARD_STEPPING, model_state_type
  use grids_mod, only : Z_INDEX
  use science_constants_mod
  use q_indices_mod, only: q_indices_add
implicit none

#ifndef TEST_MODE
  private
#endif
  real(kind=DEFAULT_PRECISION), dimension(:), allocatable :: w_buoyancy

  real(kind=DEFAULT_PRECISION) :: G_over_2

  integer :: iqv ! Index for water vapour

  public buoyancy_get_descriptor
contains

  !> Provides the descriptor back to the caller and is used in component registration
  !! @returns The termination check component descriptor
  type(component_descriptor_type) function buoyancy_get_descriptor()
    buoyancy_get_descriptor%name="buoyancy"
    buoyancy_get_descriptor%version=0.1
    buoyancy_get_descriptor%initialisation=>initialisation_callback
    buoyancy_get_descriptor%timestep=>timestep_callback
    buoyancy_get_descriptor%finalisation=>finalisation_callback

    buoyancy_get_descriptor%field_value_retrieval=>field_value_retrieval_callback
    buoyancy_get_descriptor%field_information_retrieval=>field_information_retrieval_callback
    allocate(buoyancy_get_descriptor%published_fields(1))
    buoyancy_get_descriptor%published_fields(1)="w_buoyancy"
  end function buoyancy_get_descriptor

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
    field_information%number_dimensions=1
    field_information%dimension_sizes(1)=current_state%local_grid%size(Z_INDEX)
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
    
    if (name .eq. "w_buoyancy") then
      allocate(field_value%real_1d_array(size(w_buoyancy)), source=w_buoyancy)   
    end if
  end subroutine field_value_retrieval_callback

  !> The initialisation callback sets up the buoyancy coefficient
  !! @param current_state The current model state
  subroutine initialisation_callback(current_state)
    type(model_state_type), target, intent(inout) :: current_state

    integer :: z_size

    if (.not. current_state%passive_q .and. current_state%number_q_fields > 0)then
      if (.not. allocated(current_state%cq))then
        allocate(current_state%cq(current_state%number_q_fields))
        current_state%cq=0.0_DEFAULT_PRECISION
      end if
      iqv = q_indices_add('vapour', 'buoyancy')
      current_state%cq(iqv) = ratio_mol_wts-1.0
    end if

    G_over_2 = 0.5_DEFAULT_PRECISION*G
    z_size=current_state%global_grid%size(Z_INDEX)
    allocate(w_buoyancy(z_size))
  end subroutine initialisation_callback  

  subroutine finalisation_callback(current_state)
    type(model_state_type), target, intent(inout) :: current_state

    if (allocated(w_buoyancy)) deallocate(w_buoyancy)
  end subroutine finalisation_callback

  !> Called for each column per timestep this will calculate the buoyancy terms for the SW field
  !! @param current_state The current model state
  subroutine timestep_callback(current_state)
    type(model_state_type), target, intent(inout) :: current_state

    integer :: k, n
    
#ifdef W_ACTIVE
    if (.not. current_state%passive_th .and. current_state%th%active) then
      do k=2,current_state%local_grid%size(Z_INDEX)-1    
        w_buoyancy(k)=(0.5_DEFAULT_PRECISION*current_state%global_grid%configuration%vertical%buoy_co(k))*&
             (current_state%th%data(k, current_state%column_local_y, current_state%column_local_x)&
             +current_state%th%data(k+1, current_state%column_local_y, current_state%column_local_x))
        current_state%sw%data(k, current_state%column_local_y, current_state%column_local_x)=&
             current_state%sw%data(k, current_state%column_local_y, current_state%column_local_x)+w_buoyancy(k)             
      end do
    end if
    if (.not. current_state%passive_q .and. current_state%number_q_fields .gt. 0) then
      if (current_state%use_anelastic_equations) then                                                      
        do n=1,current_state%number_q_fields
          do k=2,current_state%local_grid%size(Z_INDEX)-1            
            current_state%sw%data(k, current_state%column_local_y, current_state%column_local_x)=&
                 current_state%sw%data(k, current_state%column_local_y, current_state%column_local_x)+&
                 (0.5_DEFAULT_PRECISION*current_state%global_grid%configuration%vertical%buoy_co(k))*&
                 current_state%cq(n)* (current_state%global_grid%configuration%vertical%thref(k)*&
                 current_state%q(n)%data(k, current_state%column_local_y, current_state%column_local_x)+&
                 current_state%global_grid%configuration%vertical%thref(k+1)*&
                 current_state%q(n)%data(k+1, current_state%column_local_y, current_state%column_local_x))
          end do
        end do
      else                                                                     
        do n=1,current_state%number_q_fields
          do k=2,current_state%local_grid%size(Z_INDEX)-1
             current_state%sw%data(k, current_state%column_local_y, current_state%column_local_x)=&
                  current_state%sw%data(k, current_state%column_local_y, current_state%column_local_x)+&
                  G_over_2*current_state%cq(n)*&
                  (current_state%q(n)%data(k, current_state%column_local_y, current_state%column_local_x)+&
                  current_state%q(n)%data(k+1, current_state%column_local_y, current_state%column_local_x))
          end do
        end do
      end if
    end if
#endif
  end subroutine timestep_callback
end module buoyancy_mod
