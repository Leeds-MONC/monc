!> This component sets the source term for the lowest level (Level 1) so 
!> that, depending on surface consitionm, there is consistent lower boundary 
!> condition.
module set_consistent_lowbc_mod
  use datadefn_mod, only : DEFAULT_PRECISION
  use monc_component_mod, only : component_descriptor_type
  use optionsdatabase_mod, only : options_get_real
  use state_mod, only : PRESCRIBED_SURFACE_FLUXES, PRESCRIBED_SURFACE_VALUES, &
                        model_state_type
  use grids_mod, only : Z_INDEX
  use q_indices_mod, only: get_q_index, standard_q_names

  logical :: advect_flow, advect_th, advect_q

  integer :: iqv  ! index for vapour
  
  public set_consistent_lowbc_descriptor

contains
  !> Provides the descriptor back to the caller and is used in component registration
  !! @returns The termination check component descriptor
  type(component_descriptor_type) function set_consistent_lowbc_get_descriptor()
    set_consistent_lowbc_get_descriptor%name="set_consistent_lowbc"
    set_consistent_lowbc_get_descriptor%version=0.1
    set_consistent_lowbc_get_descriptor%initialisation=>initialisation_callback
    set_consistent_lowbc_get_descriptor%timestep=>timestep_callback
  end function set_consistent_lowbc_get_descriptor

  subroutine initialisation_callback(current_state)
    ! copy of the initialisation callback from pwadvection, used to 
    ! identify whether acting on flow, theta or q fields
    type(model_state_type), target, intent(inout) :: current_state

    ! Determine vapour index
    if (.not. current_state%passive_q) then 
       iqv = get_q_index(standard_q_names%VAPOUR, 'lowerbc')
    endif
    
  end subroutine initialisation_callback

    subroutine timestep_callback(current_state)
    type(model_state_type), target, intent(inout) :: current_state

    integer :: current_x_index, current_y_index

    current_x_index=current_state%column_local_x
    current_y_index=current_state%column_local_y

    if (current_state%halo_column) return

    call set_flow_lowbc(current_state, current_x_index, current_y_index)
    if (current_state%th%active) then 
       call set_th_lowbc(current_state, current_x_index, current_y_index)
    endif
    if (current_state%number_q_fields .gt. 0) then
       call set_q_lowbc(current_state, current_x_index, current_y_index)
    endif
  end subroutine timestep_callback

  subroutine set_flow_lowbc(current_state, current_x_index, current_y_index)
    type(model_state_type), target, intent(inout) :: current_state
    integer, intent(in) ::  current_x_index, current_y_index
    
    if (current_state%use_viscosity_and_diffusion .and. &
        current_state%use_surface_boundary_conditions) then 
#ifdef U_ACTIVE
        current_state%su%data(1, current_y_index, current_x_index)= &
               -current_state%su%data(2, current_y_index, current_x_index)
#endif
#ifdef V_ACTIVE
        current_state%sv%data(1, current_y_index, current_x_index)= &
               -current_state%sv%data(2, current_y_index, current_x_index)   
#endif
    else 
#ifdef U_ACTIVE
       current_state%su%data(1, current_y_index, current_x_index)= &               
                current_state%su%data(2, current_y_index, current_x_index)
#endif
#ifdef V_ACTIVE
       current_state%sv%data(1, current_y_index, current_x_index)= &
                current_state%sv%data(2, current_y_index, current_x_index)
#endif
    endif

  end subroutine set_flow_lowbc

  subroutine set_th_lowbc(current_state, current_x_index, current_y_index)
    type(model_state_type), target, intent(inout) :: current_state
    integer, intent(in) ::  current_x_index, current_y_index

    if (current_state%use_surface_boundary_conditions) then 
       if (current_state%type_of_surface_boundary_conditions == PRESCRIBED_SURFACE_FLUXES) then
          current_state%sth%data(1, current_y_index, current_x_index)= &
               current_state%sth%data(2, current_y_index, current_x_index)
       else
          current_state%sth%data(1, current_y_index, current_x_index)= &
               -current_state%sth%data(2, current_y_index, current_x_index)        
       endif
    endif

  end subroutine set_th_lowbc

  subroutine set_q_lowbc(current_state, current_x_index, current_y_index)
    type(model_state_type), target, intent(inout) :: current_state
    integer, intent(in) ::  current_x_index, current_y_index

    integer :: i  
    
    if (current_state%number_q_fields .gt. 0) then
       do n=1,current_state%number_q_fields
          current_state%sq(n)%data(1, current_y_index, current_x_index)= &
               current_state%sq(n)%data(2,current_y_index, current_x_index)
       enddo
    endif
    if (current_state%use_surface_boundary_conditions .and. &
         current_state%type_of_surface_boundary_conditions == PRESCRIBED_SURFACE_VALUES) then
       current_state%sq(iqv)%data(1, current_y_index, current_x_index)= &
            -(current_state%sq(iqv)%data(2,current_y_index, current_x_index))
    endif

  end subroutine set_q_lowbc

end module set_consistent_lowbc_mod
