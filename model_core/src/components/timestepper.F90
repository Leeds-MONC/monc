!> Performs the actual time stepping over groups of components. Each group can be the whole (which is one call per
!! component per timestep) or column, which calls components for each column of the timestep. Groups are executed
!! sequentially in the order that they have been configured (which is already set up in the registry)
module timestepper_mod
  use state_mod, only : model_state_type
  use grids_mod, only : X_INDEX, Y_INDEX
  use registry_mod, only : GROUP_TYPE_WHOLE, GROUP_TYPE_COLUMN, group_descriptor_type, get_ordered_groups, &
       execute_timestep_callbacks
  implicit none

#ifndef TEST_MODE
  private
#endif

  type(group_descriptor_type), dimension(:), allocatable :: group_descriptors !< Prefetched ordered group descriptors

  public init_timestepper, timestep, finalise_timestepper
contains

  !> Initialises the timestepper by prefetching the groups in the order that they will be executed, this is for optimised
  !! execution in the timestep calls
  subroutine init_timestepper()
    call get_ordered_groups(group_descriptors)
  end subroutine init_timestepper  

  !> Performs a timestep, which is comprised of executing each group of components in the order that they have been configured
  !! in. The components in a group can be called, depending on the type, just once per timestep (WHOLE) or per column (COLUMN).
  !! @param current_state The current model state
  subroutine timestep(current_state)
    type(model_state_type), intent(inout) :: current_state

    integer :: i

    do i=1,size(group_descriptors)      
      if (group_descriptors(i)%type == GROUP_TYPE_WHOLE) then
        call timestep_whole(current_state, group_descriptors(i))
      else if (group_descriptors(i)%type == GROUP_TYPE_COLUMN) then
        call timestep_column(current_state, group_descriptors(i))
      end if
    end do    
  end subroutine timestep

  !> Finalises the timestepper by cleaning up allocated memory
  subroutine finalise_timestepper()
    deallocate(group_descriptors)
  end subroutine finalise_timestepper  

  !> Performs timestepping for a group of components on a per column basis. Each component in the group is executed 
  !! for every column.
  !! @param current_state The current model state
  !! @param group_descriptor Description of the group of components to execute
  subroutine timestep_column(current_state, group_descriptor)
    type(model_state_type), intent(inout) :: current_state
    type(group_descriptor_type), intent(in) :: group_descriptor

    current_state%column_global_x=current_state%local_grid%start(X_INDEX) - current_state%local_grid%halo_size(X_INDEX)
    current_state%column_local_x=1
    do while (current_state%column_global_x .le. &
         current_state%local_grid%end(X_INDEX)+current_state%local_grid%halo_size(X_INDEX))
      current_state%column_global_y = current_state%local_grid%start(Y_INDEX) - current_state%local_grid%halo_size(Y_INDEX)
      current_state%column_local_y=1
      do while (current_state%column_global_y .le. &
           current_state%local_grid%end(Y_INDEX)+current_state%local_grid%halo_size(Y_INDEX))
        call update_state_sitation_flags(current_state)
        call execute_timestep_callbacks(current_state, group_descriptor%id)
        current_state%column_global_y = current_state%column_global_y + 1
        current_state%column_local_y = current_state%column_local_y + 1
      end do
      current_state%column_global_x = current_state%column_global_x + 1
      current_state%column_local_x = current_state%column_local_x + 1
    end do    
  end subroutine timestep_column

  !> Executes a timestep for components in a group which are designed to be executed once per timestep
  !! @param current_state The current model state
  !! @param group_descriptor Description of the group of components to execute
  subroutine timestep_whole(current_state, group_descriptor)
    type(model_state_type), intent(inout) :: current_state
    type(group_descriptor_type), intent(in) :: group_descriptor

    call execute_timestep_callbacks(current_state, group_descriptor%id)
  end subroutine timestep_whole

  !> Updates the states situation flags for easy retrieval in the components that are
  !! run per timestep
  !! @param state The current model state
  subroutine update_state_sitation_flags(current_state)
    type(model_state_type), intent(inout) :: current_state

    current_state%first_timestep_column = (current_state%column_local_x == 1 .and. current_state%column_local_y == 1)
    current_state%last_timestep_column = (current_state%column_global_x == &
         current_state%local_grid%end(X_INDEX) + current_state%local_grid%halo_size(X_INDEX) .and. &
         current_state%column_global_y == current_state%local_grid%end(Y_INDEX) + current_state%local_grid%halo_size(Y_INDEX))   

    current_state%first_nonhalo_timestep_column = (current_state%column_local_x == current_state%local_grid%halo_size(X_INDEX)+1 &
       .and. current_state%column_local_y == current_state%local_grid%halo_size(Y_INDEX)+1)

    current_state%halo_column = current_state%column_local_y .lt. current_state%local_grid%local_domain_start_index(Y_INDEX) .or.&
         current_state%column_local_x .lt. current_state%local_grid%local_domain_start_index(X_INDEX) .or.&
         current_state%column_local_y .gt. current_state%local_grid%local_domain_end_index(Y_INDEX) .or.&
         current_state%column_local_x .gt. current_state%local_grid%local_domain_end_index(X_INDEX)
  end subroutine update_state_sitation_flags
end module timestepper_mod
