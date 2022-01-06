!> Does the swapping and smoothing which is called for each column as part of the pressure-terms group of components
!! Note that this does not currently implement smoothing for mean profiles
module swapsmooth_mod
  use monc_component_mod, only : component_descriptor_type
  use state_mod, only : model_state_type, FORWARD_STEPPING
  use grids_mod, only : X_INDEX, Y_INDEX, Z_INDEX
  use datadefn_mod, only : DEFAULT_PRECISION
  implicit none

#ifndef TEST_MODE
  private
#endif

  logical :: mean_profiles_active=.false. !< Whether or not mean profiles need smoothing

  public swapsmooth_get_descriptor

contains

  !> Provides the descriptor back to the caller and is used in component registration
  !! @returns The SwapSmooth component descriptor
  type(component_descriptor_type) function swapsmooth_get_descriptor()
    swapsmooth_get_descriptor%name="swap_smooth"
    swapsmooth_get_descriptor%version=0.1
    swapsmooth_get_descriptor%initialisation=>initialisation_callback
    swapsmooth_get_descriptor%timestep=>timestep_callback
  end function swapsmooth_get_descriptor

  !> Initialises the swap and smooth component
  !! @param current_state The current model state
  subroutine initialisation_callback(current_state)
    type(model_state_type), target, intent(inout) :: current_state

#ifdef U_ACTIVE
    mean_profiles_active=allocated(current_state%global_grid%configuration%vertical%olubar)
    return
#endif
#ifdef V_ACTIVE
    mean_profiles_active=allocated(current_state%global_grid%configuration%vertical%olvbar)
    return
#endif
    if (current_state%th%active) then
      mean_profiles_active=allocated(current_state%global_grid%configuration%vertical%olthbar)   
      return
    end if
    if (current_state%number_q_fields .gt. 0) then
       mean_profiles_active=allocated(current_state%global_grid%configuration%vertical%olqbar)
       return
    end if
  end subroutine initialisation_callback  

  !> Called for each non halo timestep column and will perform swapping and smoothing as required on that column
  !! @param current_state The current model state_mod
  subroutine timestep_callback(current_state)
    type(model_state_type), target, intent(inout) :: current_state

    if (.not. current_state%halo_column) then
      if (current_state%field_stepping == FORWARD_STEPPING) then
        call swap_and_smooth_classic(current_state, .false.)
      else
        ! Centred stepping
        call swap_and_smooth_robert_filter(current_state)
      end if
    end if

    if (mean_profiles_active .and. current_state%last_timestep_column) then
      if (current_state%field_stepping == FORWARD_STEPPING) then
        call classic_for_average_profiles(current_state, .false.)
      else
        call robert_filter_for_average_profiles(current_state)
      end if
    end if    
  end subroutine timestep_callback

  !> Classic swap and smooth based upon the old or no smoothing
  !! @param current_state The current model state_mod
  !! @param old_smoother Whether to use the old smoother or not (not means no smoothing)
  subroutine swap_and_smooth_classic(current_state, old_smoother)
    type(model_state_type), intent(inout) :: current_state
    logical, intent(in) :: old_smoother

    integer :: y_index, x_index, k, n
    real(kind=DEFAULT_PRECISION) :: c1, c2, existing_value

    if (old_smoother) then
      c1 = 1.0_DEFAULT_PRECISION - current_state%tsmth                                               
      c2 = 2.0_DEFAULT_PRECISION * current_state%tsmth - 1.0_DEFAULT_PRECISION
    else
      c1 = 1.0_DEFAULT_PRECISION
      c2 = -1.0_DEFAULT_PRECISION
    end if

    x_index=current_state%column_local_x
    y_index=current_state%column_local_y
    
    do k=1,current_state%global_grid%size(Z_INDEX)
#ifdef U_ACTIVE
      existing_value = current_state%u%data(k,y_index,x_index) + current_state%zu%data(k,y_index,x_index)
      current_state%u%data(k,y_index,x_index)=existing_value * c1 + current_state%u%data(k,y_index,x_index) * c2
      current_state%zu%data(k,y_index,x_index)=existing_value - current_state%u%data(k,y_index,x_index)
#endif
#ifdef V_ACTIVE          
      existing_value = current_state%v%data(k,y_index,x_index) + current_state%zv%data(k,y_index,x_index)
      current_state%v%data(k,y_index,x_index)=existing_value * c1 + current_state%v%data(k,y_index,x_index) * c2
      current_state%zv%data(k,y_index,x_index)=existing_value - current_state%v%data(k,y_index,x_index)
#endif        
#ifdef W_ACTIVE
      existing_value = current_state%w%data(k,y_index,x_index) + current_state%zw%data(k,y_index,x_index)
      current_state%w%data(k,y_index,x_index)=existing_value * c1 + current_state%w%data(k,y_index,x_index) * c2
      current_state%zw%data(k,y_index,x_index)=existing_value - current_state%w%data(k,y_index,x_index)
#endif        
      if (current_state%th%active) then
        existing_value = current_state%th%data(k,y_index,x_index) + current_state%zth%data(k,y_index,x_index)
        current_state%th%data(k,y_index,x_index)=existing_value * c1 + current_state%th%data(k,y_index,x_index) * c2
        current_state%zth%data(k,y_index,x_index)=existing_value - current_state%th%data(k,y_index,x_index)
      end if
      do n=1,current_state%number_q_fields
        if (current_state%q(n)%active) then
          existing_value = current_state%q(n)%data(k,y_index,x_index) + current_state%zq(n)%data(k,y_index,x_index)
          current_state%q(n)%data(k,y_index,x_index)=existing_value * c1 + current_state%q(n)%data(k,y_index,x_index) * c2
          current_state%zq(n)%data(k,y_index,x_index)=existing_value - current_state%q(n)%data(k,y_index,x_index)
        end if
      end do
      do n=1,current_state%n_tracers
        existing_value = current_state%tracer(n)%data(k,y_index,x_index) + current_state%ztracer(n)%data(k,y_index,x_index)
        current_state%tracer(n)%data(k,y_index,x_index)=existing_value * c1 + current_state%tracer(n)%data(k,y_index,x_index) * c2
        current_state%ztracer(n)%data(k,y_index,x_index)=existing_value - current_state%tracer(n)%data(k,y_index,x_index)
      end do
    end do
  end subroutine swap_and_smooth_classic

  !> Swap and smooth with a Robert filter
  !! @param current_state The current model state
  subroutine swap_and_smooth_robert_filter(current_state)
    type(model_state_type), intent(inout) :: current_state

    integer :: y_index, x_index, k, n
    real(kind=DEFAULT_PRECISION) :: c1, c2, existing_value

    x_index=current_state%column_local_x
    y_index=current_state%column_local_y

    c1 = 1.0_DEFAULT_PRECISION - 2.0_DEFAULT_PRECISION*current_state%tsmth                                               
    c2 = current_state%tsmth
   
    do k=1,current_state%global_grid%size(Z_INDEX)
#ifdef U_ACTIVE
      existing_value = current_state%u%data(k,y_index,x_index)
      current_state%u%data(k,y_index,x_index)=current_state%zu%data(k,y_index,x_index)
      current_state%zu%data(k,y_index,x_index)=c1*existing_value+c2*(current_state%u%data(k, y_index, x_index)+&
           current_state%savu%data(k,y_index,x_index) -current_state%ugal)
#endif
#ifdef V_ACTIVE   
      existing_value = current_state%v%data(k,y_index,x_index)
      current_state%v%data(k,y_index,x_index)=current_state%zv%data(k,y_index,x_index)
      current_state%zv%data(k,y_index,x_index)=c1*existing_value+c2*(current_state%v%data(k, y_index, x_index)+&
           current_state%savv%data(k,y_index,x_index)-current_state%vgal)
#endif
#ifdef W_ACTIVE
      existing_value = current_state%w%data(k,y_index,x_index)
      current_state%w%data(k,y_index,x_index)=current_state%zw%data(k,y_index,x_index)
      current_state%zw%data(k,y_index,x_index)=c1*existing_value+c2*(current_state%w%data(k, y_index, x_index)+&
           current_state%savw%data(k,y_index,x_index))
#endif
      if (current_state%th%active) then
        ! Uses the partial smooth of theta from stepfields 
        existing_value = current_state%zth%data(k,y_index,x_index)
        current_state%zth%data(k,y_index,x_index)=current_state%th%data(k,y_index,x_index) + current_state%tsmth * existing_value
        current_state%th%data(k,y_index,x_index)=existing_value
      end if
      do n=1, current_state%number_q_fields
        if (current_state%q(n)%active) then
          ! Uses the partial smooth of q from stepfields 
          existing_value = current_state%zq(n)%data(k,y_index,x_index)
          current_state%zq(n)%data(k,y_index,x_index)=current_state%q(n)%data(k,y_index,x_index)+&
               current_state%tsmth * existing_value
          current_state%q(n)%data(k,y_index,x_index)=existing_value
        end if
      end do
      do n=1,current_state%n_tracers
        existing_value = current_state%ztracer(n)%data(k,y_index,x_index)
        current_state%ztracer(n)%data(k,y_index,x_index)=current_state%tracer(n)%data(k,y_index,x_index)+&
             current_state%tsmth * existing_value
        current_state%tracer(n)%data(k,y_index,x_index)=existing_value
      end do
    end do
  end subroutine swap_and_smooth_robert_filter

  !> Does swapping and smoothing (using classic algorithm) for the average profiles (the bars)
  !! @param current_state The current model state
  subroutine classic_for_average_profiles(current_state, old_smoother)
    type(model_state_type), intent(inout) :: current_state
    logical, intent(in) :: old_smoother

    integer :: k, n
    real(kind=DEFAULT_PRECISION) :: c1, c2

    if (old_smoother) then
      c1 = 1.0_DEFAULT_PRECISION - current_state%tsmth                                               
      c2 = 2.0_DEFAULT_PRECISION * current_state%tsmth - 1.0_DEFAULT_PRECISION
    else
      c1 = 1.0_DEFAULT_PRECISION
      c2 = -1.0_DEFAULT_PRECISION
    end if

    do k=1,current_state%global_grid%size(Z_INDEX)
#ifdef U_ACTIVE
      current_state%global_grid%configuration%vertical%olzubar(k)=current_state%global_grid%configuration%vertical%olubar(k) +&
           current_state%global_grid%configuration%vertical%olzubar(k)
      current_state%global_grid%configuration%vertical%olubar(k)=current_state%global_grid%configuration%vertical%olzubar(k) *&
           c1 + current_state%global_grid%configuration%vertical%olubar(k) * c2
      current_state%global_grid%configuration%vertical%olzubar(k)=current_state%global_grid%configuration%vertical%olzubar(k) -&
           current_state%global_grid%configuration%vertical%olubar(k)
#endif
#ifdef V_ACTIVE
      current_state%global_grid%configuration%vertical%olzvbar(k)=current_state%global_grid%configuration%vertical%olvbar(k) +&
           current_state%global_grid%configuration%vertical%olzvbar(k)
      current_state%global_grid%configuration%vertical%olvbar(k)=current_state%global_grid%configuration%vertical%olzvbar(k) *&
           c1 + current_state%global_grid%configuration%vertical%olvbar(k) * c2
      current_state%global_grid%configuration%vertical%olzvbar(k)=current_state%global_grid%configuration%vertical%olzvbar(k) -&
           current_state%global_grid%configuration%vertical%olvbar(k)
#endif
      if (current_state%th%active) then
        current_state%global_grid%configuration%vertical%olzthbar(k)=current_state%global_grid%configuration%vertical%olthbar(k)+&
             current_state%global_grid%configuration%vertical%olzthbar(k)
        current_state%global_grid%configuration%vertical%olthbar(k)=current_state%global_grid%configuration%vertical%olzthbar(k)*&
             c1 + current_state%global_grid%configuration%vertical%olthbar(k) * c2
        current_state%global_grid%configuration%vertical%olzthbar(k)=&
             current_state%global_grid%configuration%vertical%olzthbar(k)-&
             current_state%global_grid%configuration%vertical%olthbar(k)
      end if
      if (current_state%number_q_fields .gt. 0) then
        do n=1, current_state%number_q_fields
          current_state%global_grid%configuration%vertical%olzqbar(k,n)=&
               current_state%global_grid%configuration%vertical%olqbar(k,n)+&
               current_state%global_grid%configuration%vertical%olzqbar(k,n)
          current_state%global_grid%configuration%vertical%olqbar(k,n)=&
               current_state%global_grid%configuration%vertical%olzqbar(k,n)*&
               c1 + current_state%global_grid%configuration%vertical%olqbar(k,n) * c2
          current_state%global_grid%configuration%vertical%olzqbar(k,n)=&
               current_state%global_grid%configuration%vertical%olzqbar(k,n)-&
               current_state%global_grid%configuration%vertical%olqbar(k,n)
        end do
      end if      
    end do    
  end subroutine classic_for_average_profiles  

  !> Does swapping and smoothing (using robert filter algorithm) for the average profiles (the bars)
  !! @param current_state The current model state
  subroutine robert_filter_for_average_profiles(current_state)
    type(model_state_type), intent(inout) :: current_state

    integer :: k, n
    real(kind=DEFAULT_PRECISION) :: c1, c2, existing_value

    c1 = 1.0_DEFAULT_PRECISION - 2.0_DEFAULT_PRECISION*current_state%tsmth                                               
    c2 = current_state%tsmth

    do k=1,current_state%global_grid%size(Z_INDEX)
#ifdef U_ACTIVE
      existing_value=current_state%global_grid%configuration%vertical%olubar(k)
      current_state%global_grid%configuration%vertical%olubar(k)=current_state%global_grid%configuration%vertical%olzubar(k)
      current_state%global_grid%configuration%vertical%olzubar(k)=c1*existing_value+c2*&
           (current_state%global_grid%configuration%vertical%olubar(k) + &
           current_state%global_grid%configuration%vertical%savolubar(k))
#endif
#ifdef V_ACTIVE
      existing_value=current_state%global_grid%configuration%vertical%olvbar(k)
      current_state%global_grid%configuration%vertical%olvbar(k)=current_state%global_grid%configuration%vertical%olzvbar(k)
      current_state%global_grid%configuration%vertical%olzvbar(k)=c1*existing_value+c2*&
           (current_state%global_grid%configuration%vertical%olvbar(k) + &
           current_state%global_grid%configuration%vertical%savolvbar(k))
#endif
      if (current_state%th%active) then
        existing_value=current_state%global_grid%configuration%vertical%olzthbar(k)
        current_state%global_grid%configuration%vertical%olzthbar(k)=&
             current_state%global_grid%configuration%vertical%olthbar(k) + current_state%tsmth * existing_value
        current_state%global_grid%configuration%vertical%olthbar(k)=existing_value
      end if
      if (current_state%number_q_fields .gt. 0) then
        do n=1, current_state%number_q_fields
          existing_value=current_state%global_grid%configuration%vertical%olzqbar(k,n)
          current_state%global_grid%configuration%vertical%olzqbar(k,n)=&
               current_state%global_grid%configuration%vertical%olqbar(k,n) + current_state%tsmth * existing_value
          current_state%global_grid%configuration%vertical%olqbar(k,n)=existing_value
        end do
      end if
    end do
  end subroutine robert_filter_for_average_profiles  
end module swapsmooth_mod
