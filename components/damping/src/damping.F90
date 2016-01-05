!> Damping applied to the W field at some point to stop stuff flying up and off
module damping_mod
  use monc_component_mod, only : component_descriptor_type
  use grids_mod, only : Z_INDEX
  use state_mod, only : model_state_type
  use collections_mod, only : map_type
  use optionsdatabase_mod, only : options_get_real
  use registry_mod, only : is_component_enabled
  use logging_mod, only : LOG_ERROR, log_master_log
  use datadefn_mod, only : DEFAULT_PRECISION
  implicit none

#ifndef TEST_MODE
  private
#endif

  real(kind=DEFAULT_PRECISION) :: dmptim,& !< Layer timescale
       zdmp,& !< The point (m) where the damping starts
       hdmp   !< The height (m) of the damping layer

  public damping_get_descriptor

contains

  !> Descriptor of this component for registration
  !! @returns The damping component descriptor
  type(component_descriptor_type) function damping_get_descriptor()
    damping_get_descriptor%name="damping"
    damping_get_descriptor%version=0.1
    damping_get_descriptor%initialisation=>init_callback
    damping_get_descriptor%timestep=>timestep_callback
  end function damping_get_descriptor

  !> On initialisation will set up data structures and field values
  !! @param current_state The current model state_mod
  subroutine init_callback(current_state)
    type(model_state_type), target, intent(inout) :: current_state

    integer :: k

    if (.not. is_component_enabled(current_state%options_database, "mean_profiles")) then
      call log_master_log(LOG_ERROR, "Damping requires the mean profiles component to be enabled")
    end if    

    dmptim=options_get_real(current_state%options_database, "dmptim")
    zdmp=options_get_real(current_state%options_database, "zdmp")
    hdmp=options_get_real(current_state%options_database, "hdmp")

    allocate(current_state%global_grid%configuration%vertical%dmpco(current_state%local_grid%size(Z_INDEX)), &
         current_state%global_grid%configuration%vertical%dmpcoz(current_state%local_grid%size(Z_INDEX)))
    current_state%global_grid%configuration%vertical%dmpco(:)=0.
    current_state%global_grid%configuration%vertical%dmpcoz(:)=0.
    do k=current_state%local_grid%size(Z_INDEX),1,-1
      current_state%global_grid%configuration%vertical%kdmpmin=k
      if (current_state%global_grid%configuration%vertical%zn(k) .ge. zdmp) then
        current_state%global_grid%configuration%vertical%dmpco(k)=dmptim*(exp((&
             current_state%global_grid%configuration%vertical%zn(k)-zdmp)/hdmp)-1.0)
      end if
      if (current_state%global_grid%configuration%vertical%z(k) .ge. zdmp) then
        current_state%global_grid%configuration%vertical%dmpcoz(K)=dmptim*(exp((&
             current_state%global_grid%configuration%vertical%z(K)-zdmp)/hdmp)-1.0)
      end if
      
      if(current_state%global_grid%configuration%vertical%zn(k).lt. zdmp) exit
    end do
  end subroutine init_callback

  !> For each data column will calculate the damping term and apply this to the source term for that field
  !! @param current_state The current model state_mod
  subroutine timestep_callback(current_state)
    type(model_state_type), target, intent(inout) :: current_state

    integer :: k, i

    if (current_state%halo_column .and. current_state%timestep <3) return
    
    do k=current_state%global_grid%configuration%vertical%kdmpmin,current_state%local_grid%size(Z_INDEX)
#ifdef U_ACTIVE
      current_state%su%data(k, current_state%column_local_y, current_state%column_local_x)=current_state%su%data(k, &
           current_state%column_local_y, current_state%column_local_x)-&
           current_state%global_grid%configuration%vertical%dmpco(k)*(current_state%zu%data(k, current_state%column_local_y, &
           current_state%column_local_x)- (current_state%global_grid%configuration%vertical%olzubar(k)-current_state%ugal))
#endif
#ifdef V_ACTIVE
      current_state%sv%data(k, current_state%column_local_y, current_state%column_local_x)=current_state%sv%data(k, &
           current_state%column_local_y, current_state%column_local_x)-&
           current_state%global_grid%configuration%vertical%dmpco(k)*(current_state%zv%data(k, current_state%column_local_y, &
           current_state%column_local_x)-(current_state%global_grid%configuration%vertical%olzvbar(k)-current_state%vgal))
#endif
      if (current_state%th%active) then
        current_state%sth%data(k, current_state%column_local_y, current_state%column_local_x)=current_state%sth%data(k, &
             current_state%column_local_y, current_state%column_local_x)-&
             current_state%global_grid%configuration%vertical%dmpco(k)*(current_state%zth%data(k, current_state%column_local_y, &
             current_state%column_local_x)-current_state%global_grid%configuration%vertical%olzthbar(k))
      end if
      do i=1,current_state%number_q_fields
        if (current_state%q(i)%active) then
          current_state%sq(i)%data(k, current_state%column_local_y, current_state%column_local_x)=current_state%sq(i)%data(k, &
               current_state%column_local_y, current_state%column_local_x)-&
               current_state%global_grid%configuration%vertical%dmpco(k)*&
               (current_state%zq(i)%data(k, current_state%column_local_y, current_state%column_local_x)-&
               current_state%global_grid%configuration%vertical%olzqbar(k,i))
        end if
      end do
    end do
#ifdef W_ACTIVE
    do k=current_state%global_grid%configuration%vertical%kdmpmin,current_state%local_grid%size(Z_INDEX)-1    
      current_state%sw%data(k, current_state%column_local_y, current_state%column_local_x)=current_state%sw%data(k, &
           current_state%column_local_y, current_state%column_local_x)-&
           current_state%global_grid%configuration%vertical%dmpcoz(k)*&
           current_state%zw%data(k, current_state%column_local_y, current_state%column_local_x)
    end do
#endif
  end subroutine timestep_callback
end module damping_mod
