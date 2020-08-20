module xiosbridge_mod
  use monc_component_mod, only : component_descriptor_type
  use state_mod, only : model_state_type
  use grids_mod, only : global_grid_type, local_grid_type, &
                        X_INDEX, Y_INDEX, Z_INDEX
  use xios, only : xios_initialize, xios_finalize, xios_duration, xios_context, &
                   xios_context_initialize, xios_get_handle, &
                   xios_set_current_context, xios_context_finalize, &
                   xios_close_context_definition, xios_set_timestep, &
                   xios_update_calendar, xios_set_axis_attr, xios_set_domain_attr, &
                   xios_send_field
  use mpi
  implicit none

#ifndef TEST_MODE
  private
#endif

  public xiosbridge_get_descriptor

  integer :: comm

!!  character(len=*), parameter :: id  = "xios_monc"
!!  character(len=*), parameter :: ctx = "checkpoint"
!!  type(xios_duration) :: dtime
!  !type(xios_context)          :: ctx_hdl

!!~ real (kind=8), allocatable :: lat_val(:), lon_val(:), lval(:)



contains

  !> Provides the component descriptor for the core to register
  !! @returns The descriptor describing this component
  type(component_descriptor_type) function xiosbridge_get_descriptor()
    xiosbridge_get_descriptor%name="xiosbridge"
    xiosbridge_get_descriptor%version=0.1
    xiosbridge_get_descriptor%initialisation=>init_callback
    xiosbridge_get_descriptor%timestep=>timestep_callback
    xiosbridge_get_descriptor%finalisation=>finalisation_callback
  end function xiosbridge_get_descriptor


!=========================================================================================================================
!  !> Called on MONC initialisation
!  !! @param current_state The current model stat
  subroutine init_callback(current_state)
    type(model_state_type), target, intent(inout) :: current_state
!!    type(global_grid_type) :: gg
!!    type(local_grid_type)  :: lg

!!    integer :: i

!!    gg = current_state%global_grid
!!    lg = current_state%local_grid

!!   ! Initialisation goes here
    print*, "MMRR1 calling xios init"
!!    print*, "MMRR2 INDICES ", X_INDEX, Y_INDEX, Z_INDEX
!!    print*, "MMRR2 global grid ", gg%size(X_INDEX), gg%size(Y_INDEX), gg%size(Z_INDEX)
!!    print*, "MMRR3 local grid x", lg%local_domain_start_index(X_INDEX), lg%local_domain_end_index(X_INDEX)
!!    print*, "MMRR4 local grid y", lg%local_domain_start_index(Y_INDEX), lg%local_domain_end_index(Y_INDEX)
!!    print*, "MMRR5 local grid z", lg%local_domain_start_index(Z_INDEX), lg%local_domain_end_index(Z_INDEX)
!!    print*, "MMRR6 local grid x, start, end, size", lg%start(X_INDEX), lg%end(X_INDEX), lg%size(X_INDEX)
!!    print*, "MMRR7 local grid y, start, end, size", lg%start(Y_INDEX), lg%end(Y_INDEX), lg%size(Y_INDEX)
!!    print*, "MMRR8 local grid z, start, end, size", lg%start(Z_INDEX), lg%end(Z_INDEX), lg%size(Z_INDEX)

!    comm = current_state%parallel%monc_communicator
!    ! See pages 27-28 of XIOS Fortran reference guide
!    call xios_initialize(id, local_comm=comm)
!    call xios_context_initialize(ctx, comm)
!!    call xios_get_handle(ctx, ctx_hdl)
!!    call xios_set_current_context(ctx_hdl)     

!!    allocate(lval(lg%size(Z_INDEX)))
!!    do i = 1, lg%size(Z_INDEX)
!!      lval(i) = i - 1
!!    end do
!!    call xios_set_axis_attr("z_axis", size = gg%size(Z_INDEX), value = lval)

!!    allocate(lon_val(lg%size(X_INDEX)), lat_val(lg%size(Y_INDEX)))
!!    do i = 1, lg%size(X_INDEX)    
!!      lon_val(i) = lg%start(X_INDEX) + i - 2
!!    end do
!!    do i = 1, lg%size(Y_INDEX)
!!      lat_val(i) = lg%start(Y_INDEX) + i - 2
!!    end do

!!    call xios_set_domain_attr("domain_A", data_dim = 2, &
!!       ni_glo      = gg%size(X_INDEX), &
!!       ibegin      = lg%start(X_INDEX),      iend     = lg%end(X_INDEX), &
!!       data_ibegin = -lg%halo_size(X_INDEX), data_ni  = lg%size(X_INDEX) + 2 * lg%halo_size(X_INDEX), &
!!       nj_glo      = gg%size(Y_INDEX), &
!!       jbegin      = lg%start(Y_INDEX),      jend     = lg%end(Y_INDEX), &
!!       data_jbegin = -lg%halo_size(Y_INDEX), data_nj  = lg%size(Y_INDEX) + 2 * lg%halo_size(Y_INDEX), &
!!       lonvalue    = lon_val,                latvalue = lat_val &
!!       )

!    dtime%second = current_state%dtm
!    call xios_set_timestep(dtime)

!    call xios_close_context_definition()     

!!    call xios_update_calendar(current_state%timestep)        
  end subroutine init_callback

  !=====================================================================
  subroutine timestep_callback(current_state)
    type(model_state_type), target, intent(inout) :: current_state
!    !print*, "XIOS timestep callback"

!!~     type(global_grid_type) :: gg
!!~     type(local_grid_type)  :: lg

!!~     real (kind=8), allocatable :: tmp(:,:,:)
!!~     integer :: i, k, nx, ny, nz

!!~     gg = current_state%global_grid
!!~     lg = current_state%local_grid
!!~     nx = lg%size(X_INDEX) + 2 * lg%halo_size(X_INDEX)
!!~     ny = lg%size(Y_INDEX) + 2 * lg%halo_size(Y_INDEX)
!!~     nz = lg%size(Z_INDEX)
!!~     allocate(tmp(nx, ny, nz))

!!~     do i = 1, nx
!!~       do k = 1, nz
!!~         tmp(i, :, k) = current_state%u%data(k, :, i)
!!~       end do
!!~     end do

!!~    ! Timestep stuff goes here
!!~    print*, "MMRR called xios timestep", &
!!~            current_state%timestep, current_state%dtm, current_state%time
!     call xios_update_calendar(current_state%timestep)
!!~     call xios_send_field("u",  tmp)

!!~     deallocate(tmp)

  end subroutine timestep_callback

  !=====================================================================
  subroutine finalisation_callback(current_state)
    type(model_state_type), target, intent(inout) :: current_state

    print*, "XIOS Finalisation callback"

!!~     ! Finalisation stuff here
!!~     print*, "MMRR calling xios finalisation"
!    call xios_context_finalize()
!!~     deallocate(lat_val, lon_val, lval)
!    call xios_finalize()
  end subroutine finalisation_callback
end module xiosbridge_mod
