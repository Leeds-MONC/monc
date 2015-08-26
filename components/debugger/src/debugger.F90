!> General purpose debugger. By changing the priority and other logic we can plug it in
!! whereever we want in the run to dump out information
module debugger_mod
  use grids_mod
  use monc_component_mod
  use state_mod
  use logging_mod
  use prognostics_mod
  use conversions_mod
  implicit none

#ifndef TEST_MODE
  private
#endif

public debugger_get_descriptor

contains

  !> Provides the component descriptor for the core to register
  !! @returns The descriptor describing this component
  type(component_descriptor_type) function debugger_get_descriptor()
    debugger_get_descriptor%name="debugger"
    debugger_get_descriptor%version=0.1
    debugger_get_descriptor%initialisation=>init_callback
    debugger_get_descriptor%timestep=>timestep_callback
  end function debugger_get_descriptor

  !> Called on MONC initialisation
  !! @param current_state The current model stat
  subroutine init_callback(current_state)
    type(model_state_type), target, intent(inout) :: current_state

    call log_log(LOG_WARN, "Debugger is active - disable this for production runs")
  end subroutine init_callback  

  !> Produces debugging information on each timestep
  subroutine timestep_callback(current_state)
    type(model_state_type), target, intent(inout) :: current_state

    integer :: i, j, k
    character :: halo_classifier

    if (.not. current_state%first_timestep_column) return
!!$    if (.not. current_state%last_timestep_column) return
!!$    do i=1,current_state%local_grid%size(X_INDEX) + current_state%local_grid%halo_size(X_INDEX)*2
!!$      do j=1,current_state%local_grid%size(Y_INDEX) + current_state%local_grid%halo_size(Y_INDEX)*2
!!$        do k=1,current_state%local_grid%size(Z_INDEX) + current_state%local_grid%halo_size(Z_INDEX)*2
!!$          haloClassifier = merge('H', 'D', i .le. current_state%local_grid%halo_size(X_INDEX) .or. &
!!$               i .gt. current_state%local_grid%size(X_INDEX) + current_state%local_grid%halo_size(X_INDEX) .or. &
!!$               j .le. current_state%local_grid%halo_size(Y_INDEX) .or. &
!!$               j .gt. current_state%local_grid%size(Y_INDEX) + current_state%local_grid%halo_size(Y_INDEX) .or. &
!!$               k .le. current_state%local_grid%halo_size(Z_INDEX) .or. &
!!$               k .gt. current_state%local_grid%size(Z_INDEX) + current_state%local_grid%halo_size(Z_INDEX))
!!$          call log_log(LOG_DEBUG, "Ts: "//trim(conv_to_string(current_state%timestep))//"("//trim(conv_to_string(k))//","//&
!!$               trim(conv_to_string(j))//"," //trim(conv_to_string(i))//") q="&
!!$               //trim(conv_to_string(current_state%q(1)%data(k,j,i)))//" zq="//trim(conv_to_string(current_state%zq(1)%data(k,j,i)))&
!!$               //"("//haloClassifier//")")
!!$        end do        
!!$      end do      
!!$    end do
    do i=1,3
      write(*,*) current_state%parallel%my_rank, size(current_state%u%data, i), lbound(current_state%u%data, i), &
           ubound(current_state%u%data, i)
    end do    
  end subroutine timestep_callback  
end module debugger_mod
