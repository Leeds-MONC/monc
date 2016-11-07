module scalar_diagnostics_mod
  use monc_component_mod, only : COMPONENT_SCALAR_FIELD_TYPE, COMPONENT_DOUBLE_DATA_TYPE, &
       COMPONENT_ARRAY_FIELD_TYPE, COMPONENT_INTEGER_DATA_TYPE, &
       component_descriptor_type, component_field_value_type, component_field_information_type
  use grids_mod, only : Z_INDEX, Y_INDEX, X_INDEX
  use optionsdatabase_mod, only : options_get_real
  use state_mod, only : model_state_type
  use datadefn_mod, only : DEFAULT_PRECISION
  use science_constants_mod, only : cp, rlvap

  implicit none

#ifndef TEST_MODE
  private
#endif

  integer :: total_points, ncl_col
  real(kind=DEFAULT_PRECISION), dimension(:), allocatable :: tempfac
  real(kind=DEFAULT_PRECISION) :: qlcrit
  real(kind=DEFAULT_PRECISION), dimension(:,:), allocatable :: vwp, lwp, wmax, wmin, &
       qlmax, hqlmax, cltop, clbas,  senhf, lathf

  public scalar_diagnostics_get_descriptor

contains

  !> Provides the component descriptor for the core to register
  !! @returns The descriptor describing this component
  type(component_descriptor_type) function scalar_diagnostics_get_descriptor()
    scalar_diagnostics_get_descriptor%name="scalar_diagnostics"
    scalar_diagnostics_get_descriptor%version=0.1

    scalar_diagnostics_get_descriptor%initialisation=>initialisation_callback
    scalar_diagnostics_get_descriptor%timestep=>timestep_callback

    scalar_diagnostics_get_descriptor%field_value_retrieval=>field_value_retrieval_callback
    scalar_diagnostics_get_descriptor%field_information_retrieval=>field_information_retrieval_callback
    allocate(scalar_diagnostics_get_descriptor%published_fields(10))

    scalar_diagnostics_get_descriptor%published_fields(1)="vwp_local"
    scalar_diagnostics_get_descriptor%published_fields(2)="lwp_local"
    scalar_diagnostics_get_descriptor%published_fields(3)="qlmax_local"
    scalar_diagnostics_get_descriptor%published_fields(4)="hqlmax_local"
    scalar_diagnostics_get_descriptor%published_fields(5)="cltop_local"
    scalar_diagnostics_get_descriptor%published_fields(6)="clbas_local"
    scalar_diagnostics_get_descriptor%published_fields(7)="wmax_local"
    scalar_diagnostics_get_descriptor%published_fields(8)="wmin_local"
    scalar_diagnostics_get_descriptor%published_fields(9)="senhf_local"
    scalar_diagnostics_get_descriptor%published_fields(10)="lathf_local"
    

  end function scalar_diagnostics_get_descriptor

  subroutine initialisation_callback(current_state)
    type(model_state_type), target, intent(inout) :: current_state

    integer :: k
    integer :: y_size_local, x_size_local
    
    ! qlcrit declared in the config file
    qlcrit=options_get_real(current_state%options_database, "qlcrit")
 
    y_size_local = current_state%local_grid%size(Y_INDEX)
    x_size_local = current_state%local_grid%size(X_INDEX)
    
    ! allocate scalar diagnostics as 2-D fields (horizontal slices) that 
    ! are published to the ioserver so that the ioserver can do the manipulation
    ! to obtain the scalar field
    allocate(vwp(y_size_local, x_size_local), lwp(y_size_local, x_size_local), &
         wmax(y_size_local, x_size_local), wmin(y_size_local, x_size_local), &
         qlmax(y_size_local, x_size_local), hqlmax(y_size_local, x_size_local), &
         cltop(y_size_local, x_size_local), clbas(y_size_local, x_size_local), &
         senhf(y_size_local, x_size_local),lathf(y_size_local, x_size_local))

    allocate(tempfac(current_state%local_grid%size(Z_INDEX)))
    
    do k=2, current_state%local_grid%size(Z_INDEX)
      tempfac(k)=current_state%global_grid%configuration%vertical%dz(k)*&
           current_state%global_grid%configuration%vertical%rhon(k)
    end do    

  end subroutine initialisation_callback  

  subroutine timestep_callback(current_state)
    type(model_state_type), target, intent(inout) :: current_state

    integer :: k, i
    integer :: current_y_index, current_x_index, target_x_index, target_y_index

    current_y_index=current_state%column_local_y
    current_x_index=current_state%column_local_x
    target_y_index=current_y_index-current_state%local_grid%halo_size(Y_INDEX)
    target_x_index=current_x_index-current_state%local_grid%halo_size(X_INDEX)

    if (current_state%first_timestep_column) then
       ! water vapour path for each column
       vwp(:,:)=0.0
       ! liquid water path for each column
       lwp(:,:)=0.0
       ! maximum vertical velocity for each column
       wmax(:,:)=0.0
       ! minimum vertical velocity for each column
       wmin(:,:)=0.0
       ! maximum liquid water content in a column
       qlmax(:,:)=0.0
       ! the height of the maximum liquid water content in a column
       hqlmax(:,:)=0.0
       ! cloud top height where liqud water content is greater than qlcrit
       cltop(:,:)=0.0
       ! minimum cloud base where liquid water content is greater than q;crit
       clbas(:,:)=0.0
       ! surface sensible heat flux
       senhf(:,:)=0.0
       ! surface latent heat flux
       lathf(:,:)=0.0
    end if
    
    if (.not. current_state%halo_column) then
       if (current_state%number_q_fields .gt. 0) then
          !
          ! calculate the lwc maximum and height of lwc max for each column
          !
          if (current_state%liquid_water_mixing_ratio_index .gt. 0 .and. &
               current_state%number_q_fields .ge. current_state%liquid_water_mixing_ratio_index) then
             qlmax(target_y_index, target_x_index) = &
                  maxval(current_state%q(current_state%liquid_water_mixing_ratio_index)%data &
                  (:,current_y_index, current_x_index))
             !hqlmax(current_y_index, current_x_index) = &
             !     current_state%global_grid%configuration%vertical% &
             !     zn(maxloc(current_state%q(current_state%liquid_water_mixing_ratio_index)%data &
             !     (:,current_y_index, current_x_index)))
          end if
          !
          ! calculate the cloud top maximum and minimum for each column 
          !
          do k = 2, current_state%local_grid%size(Z_INDEX)
             if (current_state%q(current_state%liquid_water_mixing_ratio_index)%data(k, &
                  current_y_index, current_x_index) .gt. qlcrit) then
                cltop(target_y_index, target_x_index) = &
                     current_state%global_grid%configuration%vertical%zn(k)
             endif
             if (current_state%q(current_state%liquid_water_mixing_ratio_index)%data( &
                  current_state%local_grid%size(Z_INDEX)+1-k, current_y_index, current_x_index) .gt. &
                qlcrit) then
                clbas(target_y_index, target_x_index)= &
                     current_state%global_grid%configuration%vertical%zn(current_state%local_grid%size(Z_INDEX)+1-k)
             end if
          enddo
          !
          ! calculate the vapour and liquid water path
          !
          if (current_state%water_vapour_mixing_ratio_index .gt. 0 .and. &
               current_state%number_q_fields .ge. current_state%water_vapour_mixing_ratio_index) then
             do k = 2, current_state%local_grid%size(Z_INDEX)
                vwp(target_y_index, target_x_index)=vwp(target_y_index, target_x_index) &
                     +tempfac(k)*current_state%q(current_state%water_vapour_mixing_ratio_index)%data(k, &
                     current_y_index, current_x_index)
                lwp(target_y_index, target_x_index)=lwp(target_y_index, target_x_index) &
                     +tempfac(k)*current_state%q(current_state%liquid_water_mixing_ratio_index)%data(k, &
                     current_y_index, current_x_index)
             enddo
          end if
        end if

       ! surface flux diagnostics
       if (current_state%use_surface_boundary_conditions) then
          if (current_state%water_vapour_mixing_ratio_index .gt. 0 .and. &
               current_state%number_q_fields .ge. current_state%water_vapour_mixing_ratio_index) then 
             lathf(target_y_index, target_x_index)= (current_state%diff_coefficient%data(1, current_y_index, current_x_index) *   &                   
                  current_state%global_grid%configuration%vertical%rdzn(2) *  & 
                  (current_state%q(current_state%water_vapour_mixing_ratio_index)%data(1,current_y_index,current_x_index) - & 
                  current_state%q(current_state%water_vapour_mixing_ratio_index)%data(2,current_y_index,current_x_index))) &
                  * rlvap * current_state%global_grid%configuration%vertical%rhon(1)
          endif

	  if (current_state%th%active) then	     
              senhf(target_y_index, target_x_index)=(current_state%diff_coefficient%data(1, current_y_index, current_x_index)  &            
                   * current_state%global_grid%configuration%vertical%rdzn(2)  &            
                   * (current_state%th%data(1, current_y_index, current_x_index) &            
                   - current_state%th%data(2, current_y_index, current_x_index) &            
                   - current_state%global_grid%configuration%vertical%dthref(1))) &
                   * current_state%global_grid%configuration%vertical%rhon(1)*cp
          endif
       endif      
       wmax(target_y_index, target_x_index)=maxval(current_state%w%data(:, current_y_index, current_x_index))
       wmin(target_y_index, target_x_index)=minval(current_state%w%data(:, current_y_index, current_x_index))
    end if
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
    field_information%number_dimensions=2
    field_information%dimension_sizes(1)=current_state%local_grid%size(Y_INDEX)
    field_information%dimension_sizes(2)=current_state%local_grid%size(X_INDEX)
    field_information%data_type=COMPONENT_DOUBLE_DATA_TYPE

    if (name .eq. "senhf_local") then
      field_information%enabled=current_state%use_surface_boundary_conditions .and. current_state%th%active
    else if (name .eq. "lathf_local") then
      field_information%enabled=current_state%use_surface_boundary_conditions .and. &
           current_state%water_vapour_mixing_ratio_index .gt. 0 .and. &
           current_state%number_q_fields .ge. current_state%water_vapour_mixing_ratio_index
    else if (name .eq. "cltop_local" .or. name .eq. "clbas_local") then
      field_information%enabled=current_state%number_q_fields .gt. 0
    else if (name .eq. "qlmax_local") then
      field_information%enabled=current_state%number_q_fields .gt. 0 .and. current_state%liquid_water_mixing_ratio_index .gt. 0 &
           .and. current_state%number_q_fields .ge. current_state%liquid_water_mixing_ratio_index
    else if (name .eq. "vwp_local" .or. name .eq. "lwp_local") then
      field_information%enabled=current_state%number_q_fields .gt. 0 .and. current_state%water_vapour_mixing_ratio_index .gt. 0 &
           .and. current_state%number_q_fields .ge. current_state%water_vapour_mixing_ratio_index
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
    
    integer :: i

    if (name .eq. "wmax_local") then
      allocate(field_value%real_2d_array(current_state%local_grid%size(Y_INDEX), &
           current_state%local_grid%size(X_INDEX)))
       field_value%real_2d_array(:,:)=wmax(:,:)
    else if (name .eq. "wmin_local") then
       allocate(field_value%real_2d_array(current_state%local_grid%size(Y_INDEX), &
           current_state%local_grid%size(X_INDEX))) 
       field_value%real_2d_array(:,:)=wmin(:,:)
    else if (name .eq. "qlmax_local") then
       allocate(field_value%real_2d_array(current_state%local_grid%size(Y_INDEX), &
           current_state%local_grid%size(X_INDEX))) 
       field_value%real_2d_array(:,:)=qlmax(:,:)
    !else if (name .eq. "hqlmax_local") then
    !   allocate(field_value%real_2d_array(current_state%local_grid%size(Y_INDEX), &
    !       current_state%local_grid%size(X_INDEX))) 
    !   field_value%real_2d_array(:,:)=hqlmax(:,:)
    else if (name .eq. "cltop_local") then
       allocate(field_value%real_2d_array(current_state%local_grid%size(Y_INDEX), &
           current_state%local_grid%size(X_INDEX))) 
       field_value%real_2d_array(:,:)=cltop(:,:)
    else if (name .eq. "clbas_local") then
       allocate(field_value%real_2d_array(current_state%local_grid%size(Y_INDEX), &
           current_state%local_grid%size(X_INDEX))) 
       field_value%real_2d_array(:,:)=clbas(:,:)
    else if (name .eq. "vwp_local") then
       allocate(field_value%real_2d_array(current_state%local_grid%size(Y_INDEX), &
           current_state%local_grid%size(X_INDEX))) 
       field_value%real_2d_array(:,:)=vwp(:,:)
    else if (name .eq. "lwp_local") then
       allocate(field_value%real_2d_array(current_state%local_grid%size(Y_INDEX), &
           current_state%local_grid%size(X_INDEX))) 
       field_value%real_2d_array(:,:)=lwp(:,:) 
    else if (name .eq. "senhf_local") then
       allocate(field_value%real_2d_array(current_state%local_grid%size(Y_INDEX), &
           current_state%local_grid%size(X_INDEX))) 
       field_value%real_2d_array(:,:)=senhf(:,:)
    else if (name .eq. "lathf_local") then
       allocate(field_value%real_2d_array(current_state%local_grid%size(Y_INDEX), &
           current_state%local_grid%size(X_INDEX))) 
       field_value%real_2d_array(:,:)=lathf(:,:) 
    end if
    
  end subroutine field_value_retrieval_callback
end module scalar_diagnostics_mod
