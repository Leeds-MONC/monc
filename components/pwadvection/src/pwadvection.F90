!> Piacsek-Williams advection scheme
module pwadvection_mod
  use datadefn_mod, only : DEFAULT_PRECISION
  use monc_component_mod, only : component_descriptor_type
  use optionsdatabase_mod, only : options_get_string
  use collections_mod, only : map_type
  use state_mod, only : model_state_type
  use grids_mod, only : Z_INDEX
implicit none

#ifndef TEST_MODE
  private
#endif

  logical :: advect_flow, advect_th, advect_q

  logical :: l_toplevel=.false.

 public pwadvection_get_descriptor
contains

  !> Provides a description of this component for the core to register
  !! @returns The descriptor containing registration information for this component
  type(component_descriptor_type) function pwadvection_get_descriptor()
    pwadvection_get_descriptor%name="pw_advection"
    pwadvection_get_descriptor%version=0.1
    pwadvection_get_descriptor%initialisation=>initialisation_callback
    pwadvection_get_descriptor%timestep=>timestep_callback
  end function pwadvection_get_descriptor

  !> Initialisation callback, will set up the configuration of this advection scheme
  !! @param current_state The current model state
  subroutine initialisation_callback(current_state)
    type(model_state_type), target, intent(inout) :: current_state

    advect_flow=determine_if_advection_here(options_get_string(current_state%options_database, "advection_flow_fields"))    
    advect_th=determine_if_advection_here(options_get_string(current_state%options_database, "advection_theta_field"))
    advect_q=determine_if_advection_here(options_get_string(current_state%options_database, "advection_q_fields"))    
  end subroutine initialisation_callback  

  !> Called per column of data, this will perform Piacsek-Williams advection on the applicable fields for non halo data
  !! @param current_state The current model state
  subroutine timestep_callback(current_state)
    type(model_state_type), target, intent(inout) :: current_state

    integer :: current_x_index, current_y_index

    current_x_index=current_state%column_local_x
    current_y_index=current_state%column_local_y

    if (current_state%halo_column) return

    if (advect_flow) call advect_flow_fields(current_state, current_x_index, current_y_index)
    if (advect_th) call advect_th_field(current_state, current_x_index, current_y_index)
    if (advect_q) call advect_q_field(current_state, current_x_index, current_y_index)
  end subroutine timestep_callback

  !> Advects the q fields in a column
  !! @param current_state The current model state
  !! @param current_x_index The current slice, x, index
  !! @param current_y_index The current column, y, index
  subroutine advect_q_field(current_state, current_x_index, current_y_index)
    type(model_state_type), target, intent(inout) :: current_state
    integer, intent(in) ::  current_x_index, current_y_index

    integer :: k, n

    do n=1,current_state%number_q_fields
      do k=2,current_state%local_grid%size(Z_INDEX)-1
        current_state%sq(n)%data(k, current_y_index, current_x_index)=0.0_DEFAULT_PRECISION
#ifdef U_ACTIVE
        current_state%sq(n)%data(k, current_y_index, current_x_index)=&
             current_state%sq(n)%data(k, current_y_index, current_x_index)+current_state%global_grid%configuration%horizontal%cx*&
             0.5_DEFAULT_PRECISION*(current_state%u%data(k, current_y_index, current_x_index-1)*&
             current_state%q(n)%data(k, current_y_index, current_x_index-1)-&
             current_state%u%data(k, current_y_index, current_x_index)*&
             current_state%q(n)%data(k, current_y_index, current_x_index+1))
#endif
#ifdef V_ACTIVE
        current_state%sq(n)%data(k, current_y_index, current_x_index)=&
             current_state%sq(n)%data(k, current_y_index, current_x_index)+&
             current_state%global_grid%configuration%horizontal%cy*0.5_DEFAULT_PRECISION*&
             (current_state%v%data(k, current_y_index-1, current_x_index)*&
             current_state%q(n)%data(k, current_y_index-1, current_x_index)-&
             current_state%v%data(k, current_y_index, current_x_index)*&
             current_state%q(n)%data(k, current_y_index+1, current_x_index))
#endif
#ifdef W_ACTIVE
        current_state%sq(n)%data(k, current_y_index, current_x_index)=&
             current_state%sq(n)%data(k, current_y_index, current_x_index)+&
             2.0_DEFAULT_PRECISION*(current_state%global_grid%configuration%vertical%tzc1(k)*&
             current_state%w%data(k-1, current_y_index, current_x_index)*&
             current_state%q(n)%data(k-1, current_y_index, current_x_index)-&
             current_state%global_grid%configuration%vertical%tzc2(k)*&
             current_state%w%data(k, current_y_index, current_x_index)*&
             current_state%q(n)%data(k+1, current_y_index, current_x_index))
#endif
      end do

      if (l_toplevel)then
      k=current_state%local_grid%size(Z_INDEX)
      current_state%sq(n)%data(k, current_y_index, current_x_index)=0.0_DEFAULT_PRECISION
#ifdef U_ACTIVE
      current_state%sq(n)%data(k, current_y_index, current_x_index)=&
           current_state%sq(n)%data(k, current_y_index, current_x_index)+current_state%global_grid%configuration%horizontal%cx*&
           0.5_DEFAULT_PRECISION*(current_state%u%data(k, current_y_index, current_x_index-1)*&
           current_state%q(n)%data(k, current_y_index, current_x_index-1)-&
           current_state%u%data(k, current_y_index, current_x_index)*&
           current_state%q(n)%data(k, current_y_index, current_x_index+1))
#endif
#ifdef V_ACTIVE
      current_state%sq(n)%data(k, current_y_index, current_x_index)=&
           current_state%sq(n)%data(k, current_y_index, current_x_index)+current_state%global_grid%configuration%horizontal%cy*&
           0.5_DEFAULT_PRECISION*(current_state%v%data(k, current_y_index-1, current_x_index)*&
           current_state%q(n)%data(k, current_y_index-1, current_x_index)-&
           current_state%v%data(k, current_y_index, current_x_index)*&
           current_state%q(n)%data(k, current_y_index+1, current_x_index))
#endif
#ifdef W_ACTIVE
      current_state%sq(n)%data(k, current_y_index, current_x_index)=&
           current_state%sq(n)%data(k, current_y_index, current_x_index)+&
           current_state%global_grid%configuration%vertical%tzc1(k)*2.0_DEFAULT_PRECISION*&
           current_state%w%data(k-1, current_y_index, current_x_index)*&
           current_state%q(n)%data(k-1, current_y_index, current_x_index)
#endif
    end if
    end do
  end subroutine advect_q_field  

  !> Advects the theta field in a column
  !! @param current_state The current model state
  !! @param current_x_index The current slice, x, index
  !! @param current_y_index The current column, y, index
  subroutine advect_th_field(current_state, current_x_index, current_y_index)
    type(model_state_type), target, intent(inout) :: current_state
    integer, intent(in) ::  current_x_index, current_y_index

    integer :: k

    if (current_state%th%active) then
      do k=2,current_state%local_grid%size(Z_INDEX)-1
        current_state%sth%data(k, current_y_index, current_x_index)=0.0_DEFAULT_PRECISION
#ifdef U_ACTIVE
        current_state%sth%data(k, current_y_index, current_x_index)=current_state%sth%data(k, current_y_index, current_x_index)+&
             current_state%global_grid%configuration%horizontal%cx*&
             0.5_DEFAULT_PRECISION*(current_state%u%data(k, current_y_index, current_x_index-1)*&
             current_state%th%data(k, current_y_index, current_x_index-1)-&
             current_state%u%data(k, current_y_index, current_x_index)*&
             current_state%th%data(k, current_y_index, current_x_index+1))
#endif
#ifdef V_ACTIVE
        current_state%sth%data(k, current_y_index, current_x_index)=current_state%sth%data(k, current_y_index, current_x_index)+&
             current_state%global_grid%configuration%horizontal%cy*0.5_DEFAULT_PRECISION*&
             (current_state%v%data(k, current_y_index-1, current_x_index)*&
             current_state%th%data(k, current_y_index-1, current_x_index)-&
             current_state%v%data(k, current_y_index, current_x_index)*&
             current_state%th%data(k, current_y_index+1, current_x_index))
#endif
#ifdef W_ACTIVE
        current_state%sth%data(k, current_y_index, current_x_index)=current_state%sth%data(k, current_y_index, current_x_index)+&
             2.0_DEFAULT_PRECISION*(current_state%global_grid%configuration%vertical%tzc1(k)*&
             current_state%w%data(k-1, current_y_index, current_x_index)*&
             current_state%th%data(k-1, current_y_index, current_x_index)-&
             current_state%global_grid%configuration%vertical%tzc2(k)*&
             current_state%w%data(k, current_y_index, current_x_index)*&
             current_state%th%data(k+1, current_y_index, current_x_index))
#endif
      end do

      if (l_toplevel)then
        k=current_state%local_grid%size(Z_INDEX)
        current_state%sth%data(k, current_y_index, current_x_index)=0.0_DEFAULT_PRECISION
#ifdef U_ACTIVE
        current_state%sth%data(k, current_y_index, current_x_index)=current_state%sth%data(k, current_y_index, current_x_index)+&
             current_state%global_grid%configuration%horizontal%cx*&
             0.5_DEFAULT_PRECISION*(current_state%u%data(k, current_y_index, current_x_index-1)*&
             current_state%th%data(k, current_y_index, current_x_index-1)-&
             current_state%u%data(k, current_y_index, current_x_index)*&
             current_state%th%data(k, current_y_index, current_x_index+1))
#endif
#ifdef V_ACTIVE
        current_state%sth%data(k, current_y_index, current_x_index)=current_state%sth%data(k, current_y_index, current_x_index)+&
             current_state%global_grid%configuration%horizontal%cy*0.5_DEFAULT_PRECISION*&
             (current_state%v%data(k, current_y_index-1, current_x_index)*&
             current_state%th%data(k, current_y_index-1, current_x_index)-&
             current_state%v%data(k, current_y_index, current_x_index)*&
             current_state%th%data(k, current_y_index+1, current_x_index))
#endif
#ifdef W_ACTIVE
        current_state%sth%data(k, current_y_index, current_x_index)=current_state%sth%data(k, current_y_index, current_x_index)+&
             current_state%global_grid%configuration%vertical%tzc1(k)*2.0_DEFAULT_PRECISION*&
             current_state%w%data(k-1, current_y_index, current_x_index)*current_state%th%data(k-1, current_y_index, &
             current_x_index)
#endif
      end if
    end if
  end subroutine advect_th_field  

  !> Advects the flow fields depending upon which fields are active in the model in a column
  !! @param current_state The current model state
  !! @param current_x_index The current slice, x, index
  !! @param current_y_index The current column, y, index
  subroutine advect_flow_fields(current_state, current_x_index, current_y_index)
    type(model_state_type), target, intent(inout) :: current_state
    integer, intent(in) ::  current_x_index, current_y_index

    integer :: k

    do k=2,current_state%local_grid%size(Z_INDEX)-1      
#ifdef U_ACTIVE
      current_state%su%data(k, current_y_index, current_x_index)=&
           current_state%global_grid%configuration%horizontal%tcx*(current_state%u%data(k, current_y_index, current_x_index-1)*&
           (current_state%u%data(k, current_y_index, current_x_index)+&
           current_state%u%data(k, current_y_index, current_x_index-1))-&
           current_state%u%data(k, current_y_index, current_x_index+1)*&
           (current_state%u%data(k, current_y_index, current_x_index)+&
           current_state%u%data(k, current_y_index, current_x_index+1)))
#ifdef V_ACTIVE
      current_state%su%data(k, current_y_index, current_x_index)=current_state%su%data(k, current_y_index, current_x_index)+&
           current_state%global_grid%configuration%horizontal%tcy*(current_state%u%data(k, current_y_index-1, current_x_index)*&
           (current_state%v%data(k, current_y_index-1, current_x_index)+&
           current_state%v%data(k, current_y_index-1, current_x_index+1))-&
           current_state%u%data(k, current_y_index+1, current_x_index)*&
           (current_state%v%data(k, current_y_index, current_x_index)+&
           current_state%v%data(k, current_y_index, current_x_index+1)))
#endif
#ifdef W_ACTIVE
      current_state%su%data(k, current_y_index, current_x_index)=current_state%su%data(k, current_y_index, current_x_index)+&
           (current_state%global_grid%configuration%vertical%tzc1(k)*current_state%u%data(k-1, current_y_index, current_x_index)*&
           (current_state%w%data(k-1, current_y_index, current_x_index)+&
           current_state%w%data(k-1, current_y_index, current_x_index+1))-&
           current_state%global_grid%configuration%vertical%tzc2(k)*current_state%u%data(k+1, current_y_index, current_x_index)*&
           (current_state%w%data(k, current_y_index, current_x_index)+&
           current_state%w%data(k, current_y_index, current_x_index+1)))
#endif
#endif

#ifdef V_ACTIVE
      current_state%sv%data(k, current_y_index, current_x_index)=current_state%global_grid%configuration%horizontal%tcy*(&
           current_state%v%data(k, current_y_index-1, current_x_index)*&
           (current_state%v%data(k, current_y_index, current_x_index)+&
           current_state%v%data(k, current_y_index-1, current_x_index))-&
           current_state%v%data(k, current_y_index+1, current_x_index)*&
           (current_state%v%data(k, current_y_index, current_x_index)+&
           current_state%v%data(k, current_y_index+1, current_x_index)))
#ifdef U_ACTIVE
      current_state%sv%data(k, current_y_index, current_x_index)=current_state%sv%data(k, current_y_index, current_x_index)+&
           current_state%global_grid%configuration%horizontal%tcx*(current_state%v%data(k, current_y_index, current_x_index-1)*&
           (current_state%u%data(k, current_y_index, current_x_index-1)+&
           current_state%u%data(k, current_y_index+1, current_x_index-1))-&
           current_state%v%data(k, current_y_index, current_x_index+1)*&
           (current_state%u%data(k, current_y_index, current_x_index)+&
           current_state%u%data(k, current_y_index+1, current_x_index)))
#endif
#ifdef W_ACTIVE
      current_state%sv%data(k, current_y_index, current_x_index)=current_state%sv%data(k, current_y_index, current_x_index)+&
           (current_state%global_grid%configuration%vertical%tzc1(k)*current_state%v%data(k-1, current_y_index, current_x_index)*&
           (current_state%w%data(k-1, current_y_index, current_x_index)+&
           current_state%w%data(k-1, current_y_index+1, current_x_index))-&
           current_state%global_grid%configuration%vertical%tzc2(k)*current_state%v%data(k+1, current_y_index, current_x_index)*&
           (current_state%w%data(k, current_y_index, current_x_index)+&
           current_state%w%data(k, current_y_index+1, current_x_index)))
#endif
#endif

#ifdef W_ACTIVE
      current_state%sw%data(k, current_y_index, current_x_index)=(current_state%global_grid%configuration%vertical%tzd1(k)*&
           current_state%w%data(k-1, current_y_index, current_x_index)*&
           (current_state%w%data(k, current_y_index, current_x_index)+&
           current_state%w%data(k-1, current_y_index, current_x_index))-&
           current_state%global_grid%configuration%vertical%tzd2(k)*current_state%w%data(k+1, current_y_index, current_x_index)*&
           (current_state%w%data(k, current_y_index, current_x_index)+&
           current_state%w%data(k+1, current_y_index, current_x_index)))
#ifdef U_ACTIVE
      current_state%sw%data(k, current_y_index, current_x_index)=current_state%sw%data(k, current_y_index, current_x_index)+&
           current_state%global_grid%configuration%horizontal%tcx*(current_state%w%data(k, current_y_index, current_x_index-1)*&
           (current_state%u%data(k, current_y_index, current_x_index-1)+&
           current_state%u%data(k+1, current_y_index, current_x_index-1))-&
           current_state%w%data(k, current_y_index, current_x_index+1)*&
           (current_state%u%data(k, current_y_index, current_x_index)+&
           current_state%u%data(k+1, current_y_index, current_x_index)))
#endif
#ifdef V_ACTIVE
      current_state%sw%data(k, current_y_index, current_x_index)=current_state%sw%data(k, current_y_index, current_x_index)+&
           current_state%global_grid%configuration%horizontal%tcy*(current_state%w%data(k, current_y_index-1, current_x_index)*&
           (current_state%v%data(k, current_y_index-1, current_x_index)+&
           current_state%v%data(k+1, current_y_index-1, current_x_index))-&
           current_state%w%data(k, current_y_index+1, current_x_index)*&
           (current_state%v%data(k, current_y_index, current_x_index)+&
           current_state%v%data(k+1, current_y_index, current_x_index)))
#endif
#endif
    end do

    if (l_toplevel)then
    k=current_state%local_grid%size(Z_INDEX)
#ifdef U_ACTIVE
    current_state%su%data(k, current_y_index, current_x_index)=current_state%global_grid%configuration%horizontal%tcx*&
         (current_state%u%data(k, current_y_index, current_x_index-1)*&
         (current_state%u%data(k, current_y_index, current_x_index)+&
         current_state%u%data(k, current_y_index, current_x_index-1))-&
         current_state%u%data(k, current_y_index, current_x_index+1)*&
         (current_state%u%data(k, current_y_index, current_x_index)+&
         current_state%u%data(k, current_y_index, current_x_index+1)))
#ifdef V_ACTIVE
    current_state%su%data(k, current_y_index, current_x_index)=current_state%su%data(k, current_y_index, current_x_index)+&
         current_state%global_grid%configuration%horizontal%tcy*(current_state%u%data(k, current_y_index-1, current_x_index)*&
         (current_state%v%data(k, current_y_index-1, current_x_index)+&
         current_state%v%data(k, current_y_index-1, current_x_index+1))-&
         current_state%u%data(k, current_y_index+1, current_x_index)*&
         (current_state%v%data(k, current_y_index, current_x_index)+&
         current_state%v%data(k, current_y_index, current_x_index+1)))
#endif
#ifdef W_ACTIVE
    current_state%su%data(k, current_y_index, current_x_index)=current_state%su%data(k, current_y_index, current_x_index)+&
         current_state%global_grid%configuration%vertical%tzc1(k)*current_state%u%data(k-1, current_y_index, current_x_index)*&
         (current_state%w%data(k-1, current_y_index, current_x_index)+&
         current_state%w%data(k-1, current_y_index, current_x_index+1))
#endif
#endif

#ifdef V_ACTIVE
    current_state%sv%data(k, current_y_index, current_x_index)=current_state%global_grid%configuration%horizontal%tcy*&
         (current_state%v%data(k, current_y_index-1, current_x_index)*&
         (current_state%v%data(k, current_y_index, current_x_index)+&
         current_state%v%data(k, current_y_index-1, current_x_index))-&
         current_state%v%data(k, current_y_index+1, current_x_index)*&
         (current_state%v%data(k, current_y_index, current_x_index)+&
         current_state%v%data(k, current_y_index+1, current_x_index)))
#ifdef U_ACTIVE
    current_state%sv%data(k, current_y_index, current_x_index)=current_state%sv%data(k, current_y_index, current_x_index)+&
         current_state%global_grid%configuration%horizontal%tcx*(current_state%v%data(k, current_y_index, current_x_index-1)*&
         (current_state%u%data(k, current_y_index, current_x_index-1)+&
         current_state%u%data(k, current_y_index+1, current_x_index-1))-&
         current_state%v%data(k, current_y_index, current_x_index+1)*&
         (current_state%u%data(k, current_y_index, current_x_index)+&
         current_state%u%data(k, current_y_index+1, current_x_index)))
#endif
#ifdef W_ACTIVE
    current_state%sv%data(k, current_y_index, current_x_index)=current_state%sv%data(k, current_y_index, current_x_index)+&
         current_state%global_grid%configuration%vertical%tzc1(k)*current_state%v%data(k-1, current_y_index, current_x_index)*&
         (current_state%w%data(k-1, current_y_index, current_x_index)+&
         current_state%w%data(k-1, current_y_index+1, current_x_index))
#endif
#endif
  end if
  end subroutine advect_flow_fields

  !> Parses a field string (read in from the configuration file) and determines whether this algorithm should be used
  !! for advecting that field
  !! @param field The string configuration of field advection
  !! @returns Whether or not the field is advected here
  logical function determine_if_advection_here(field)
    character(len=*), intent(in) :: field

    if (len_trim(field) .ne. 0) then
      if (trim(field) .eq. "pw" .or. trim(field) .eq. "any") then
        determine_if_advection_here=.true.
      else
        determine_if_advection_here=.false.
      end if
    else
      determine_if_advection_here=.true.
    end if
  end function determine_if_advection_here  
end module pwadvection_mod
