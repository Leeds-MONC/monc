!> Applies lateral boundary conditions. 
!! Note that these boundary conditions only make sense if bi-periodicity is switched off in the 
!! the halo swapping and solver.  The FFT solver should have bi-periodic conditions.
module lateral_bcs_mod
  use monc_component_mod, only : component_descriptor_type
  use state_mod, only : model_state_type
  use grids_mod, only : X_INDEX, Y_INDEX, Z_INDEX
  use logging_mod, only : LOG_DEBUG, log_get_logging_level, log_log
  use datadefn_mod, only : DEFAULT_PRECISION
  implicit none

#ifndef TEST_MODE
  private
#endif  

  INTEGER :: lbc_type_x, lbc_type_y ! probably need to add these to state

  INTEGER :: idepth, jdepth ! depth of lateral bcs in i and j directions

  INTEGER, PARAMETER :: LBC_PERIODIC = 1
  INTEGER, PARAMETER :: LBC_RIGID    = 2
  INTEGER, PARAMETER :: LBC_BLEND   = 3

  public lateral_bcs_get_descriptor

contains

  !> Provides a description of this component for the core to register
  !! @returns The descriptor containing registration information for this component
  type(component_descriptor_type) function lateral_bcs_get_descriptor()
    lateral_bcs_get_descriptor%name="lateral_bcs"
    lateral_bcs_get_descriptor%version=0.1
    lateral_bcs_get_descriptor%initialisation=>initialisation_callback
    lateral_bcs_get_descriptor%timestep=>timestep_callback
  end function lateral_bcs_get_descriptor  

  !> Initialisation callback.  Sets up boundary depth and other configurable options
  !! @param current_state The current model state
  subroutine initialisation_callback(current_state)
    type(model_state_type), target, intent(inout) :: current_state
    
    ! Hard wiring these for now
    idepth=3
    jdepth=2

    lbc_type_x=LBC_RIGID
    lbc_type_y=LBC_PERIODIC

  end subroutine initialisation_callback  

  !> Timestep callback which applies lateral boundary conditions 
  !!
  !! @param current_state The current model state_mod
  subroutine timestep_callback(current_state)
    type(model_state_type), target, intent(inout) :: current_state

    INTEGER :: istart, iend, jstart, jend, kstart, kend

    INTEGER :: i,j,k,iq

    LOGICAL :: is_east, is_west, is_north, is_south

    if (current_state%timestep==1) return

    is_west=current_state%parallel%wrapped_around(X_INDEX,1) 
    is_east=current_state%parallel%wrapped_around(X_INDEX,2) 
    is_south=current_state%parallel%wrapped_around(Y_INDEX,1) 
    is_north=current_state%parallel%wrapped_around(Y_INDEX,2) 

    ! k points
    if (is_west .or. is_east .or. is_north .or. is_south)then
      kstart=current_state%local_grid%local_domain_start_index(Z_INDEX)
      kend=current_state%local_grid%local_domain_end_index(Z_INDEX)
    end if

    ! Eastern boundary
    if (is_east)then

      iend=current_state%local_grid%local_domain_end_index(X_INDEX) 
      istart = iend - idepth + 1
      jstart=current_state%local_grid%local_domain_start_index(Y_INDEX)
      jend=current_state%local_grid%local_domain_end_index(Y_INDEX)

      if (lbc_type_x==LBC_PERIODIC)then
        ! how should this be done best Luis?
      else if (lbc_type_x==LBC_RIGID)then
        ! No communications needed assuming idepth and jdepth fit 
        ! on a single pe (what if they don't?)
        call apply_rigid(current_state,istart,iend,jstart,jend,kstart,kend)
        ! propagate winds to halos
        do i=iend+1,iend+current_state%local_grid%halo_size(X_INDEX)
          current_state%u%data(:,:,i) = current_state%u%data(:,:,iend)
          current_state%zu%data(:,:,i) = current_state%zu%data(:,:,iend)
          current_state%w%data(:,:,i) = current_state%w%data(:,:,iend)
          current_state%zw%data(:,:,i) = current_state%zw%data(:,:,iend)
        end do
        
      else if (lbc_type_x==LBC_BLEND)then
        ! To be coded....
      end if

    end if

    ! Western boundary
    if (is_west)then

      istart=current_state%local_grid%local_domain_start_index(X_INDEX)
      iend=istart + idepth - 1
      jstart=current_state%local_grid%local_domain_start_index(Y_INDEX)
      jend=current_state%local_grid%local_domain_end_index(Y_INDEX)


      if (lbc_type_x==LBC_PERIODIC)then
        ! how should this be done best Luis?
      else if (lbc_type_x==LBC_RIGID)then
        ! No communications needed assuming idepth and jdepth fit 
        ! on a single pe (what if they don't?)
        call apply_rigid(current_state,istart,iend,jstart,jend,kstart,kend)
        ! propagate winds to halos
        do i=istart-current_state%local_grid%halo_size(X_INDEX),istart-1
          current_state%u%data(:,:,i) = current_state%u%data(:,:,istart)
          current_state%zu%data(:,:,i) = current_state%zu%data(:,:,istart)
          current_state%w%data(:,:,i) = current_state%w%data(:,:,istart)
          current_state%zw%data(:,:,i) = current_state%zw%data(:,:,istart)
        end do
      else if (lbc_type_x==LBC_BLEND)then
        ! To be coded....
      end if
    end if
      
    
    ! Northern boundary
    if (is_north)then

      istart=current_state%local_grid%local_domain_start_index(X_INDEX)
      iend=current_state%local_grid%local_domain_end_index(X_INDEX)
      jend=current_state%local_grid%local_domain_end_index(Y_INDEX)
      jstart=jend - jdepth + 1

      if (lbc_type_y==LBC_PERIODIC)then
        ! how should this be done best Luis?
      else if (lbc_type_y==LBC_RIGID)then
        ! No communications needed assuming idepth and jdepth fit 
        ! on a single pe (what if they don't?)
        call apply_rigid(current_state,istart,iend,jstart,jend,kstart,kend)
        do j=jend+1,jend+current_state%local_grid%halo_size(Y_INDEX)
          current_state%v%data(:,:,j) = current_state%v%data(:,:,jend)
          current_state%zv%data(:,:,j) = current_state%zv%data(:,:,jend)
          current_state%w%data(:,:,j) = current_state%w%data(:,:,jend)
          current_state%zw%data(:,:,j) = current_state%zw%data(:,:,jend)
        end do
      else if (lbc_type_y==LBC_BLEND)then
        ! To be coded....
      end if

    end if
      

    ! Southern boundary
    if (is_south)then

      istart=current_state%local_grid%local_domain_start_index(X_INDEX)
      iend=current_state%local_grid%local_domain_end_index(X_INDEX)
      jstart=current_state%local_grid%local_domain_start_index(Y_INDEX)
      jend=jstart + jdepth - 1

      if (lbc_type_y==LBC_PERIODIC)then
        ! how should this be done best Luis?
      else if (lbc_type_y==LBC_RIGID)then
        ! No communications needed assuming idepth and jdepth fit 
        ! on a single pe (what if they don't?)
        call apply_rigid(current_state,istart,iend,jstart,jend,kstart,kend)
        do j=jstart-current_state%local_grid%halo_size(Y_INDEX),jstart-1
          current_state%v%data(:,:,j) = current_state%v%data(:,:,jstart)
          current_state%zv%data(:,:,j) = current_state%zv%data(:,:,jstart)
          current_state%w%data(:,:,j) = current_state%w%data(:,:,jstart)
          current_state%zw%data(:,:,j) = current_state%zw%data(:,:,jstart)
        end do
      else if (lbc_type_y==LBC_BLEND)then
        ! To be coded....
      end if

    end if
  end subroutine timestep_callback

  subroutine apply_rigid(current_state,istart,iend,jstart,jend,kstart,kend)
    type(model_state_type), target, intent(inout) :: current_state

    INTEGER :: istart, iend, jstart, jend, kstart, kend

    INTEGER :: i,j,k,iq,it

    do i=istart, iend
      do j=jstart, jend
        do k=kstart,kend
          ! Update to old values, i.e. fixed in time
#ifdef U_ACTIVE
!          current_state%u%data(k,j,i) = current_state%zu%data(k,j,i)
#endif
#ifdef V_ACTIVE
!          current_state%v%data(k,j,i) = current_state%zv%data(k,j,i)
#endif
#ifdef W_ACTIVE
!          current_state%w%data(k,j,i) = current_state%zw%data(k,j,i)
#endif
          current_state%th%data(k,j,i) = current_state%zth%data(k,j,i)
        end do
      end do
    end do
    do iq = 1, current_state%number_q_fields
      do i=istart, iend
        do j=jstart, jend
          do k=kstart,kend
            current_state%q(iq)%data(k,j,i) = current_state%zq(iq)%data(k,j,i) 
          end do
        end do
      end do
    end do
    do it = 1, current_state%n_tracers
      do i=istart, iend
        do j=jstart, jend
          do k=kstart,kend
            current_state%tracer(it)%data(k,j,i) = current_state%ztracer(it)%data(k,j,i) 
          end do
        end do
      end do
    end do

  end subroutine apply_rigid

end module lateral_bcs_mod
