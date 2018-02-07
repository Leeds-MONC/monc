module scalar_diagnostics_mod
  use monc_component_mod, only : COMPONENT_SCALAR_FIELD_TYPE, COMPONENT_DOUBLE_DATA_TYPE, &
       COMPONENT_ARRAY_FIELD_TYPE, COMPONENT_INTEGER_DATA_TYPE, &
       component_descriptor_type, component_field_value_type, component_field_information_type
  use grids_mod, only : Z_INDEX, Y_INDEX, X_INDEX
  use optionsdatabase_mod, only : options_get_real
  use state_mod, only : model_state_type
  use datadefn_mod, only : DEFAULT_PRECISION
  use science_constants_mod, only : cp, rlvap
  use registry_mod, only : is_component_enabled

  implicit none

#ifndef TEST_MODE
  private
#endif

  integer :: total_points, ncl_col, iqv=0, iql=0, iqr=0, iqi=0, iqs=0,    &
       iqg=0
  real(kind=DEFAULT_PRECISION), dimension(:), allocatable :: dz_rhon_fac
  real(kind=DEFAULT_PRECISION), dimension(:), allocatable :: ww_prime_res, uu_prime_res, &
       vv_prime_res
  real(kind=DEFAULT_PRECISION) :: qlcrit
  real(kind=DEFAULT_PRECISION), dimension(:,:), allocatable :: vwp, lwp, wmax, wmin, &
       qlmax, hqlmax, cltop, clbas,  senhf, lathf, rwp, iwp, swp, gwp, tot_iwp,      &
       reske

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
    allocate(scalar_diagnostics_get_descriptor%published_fields(16))

    scalar_diagnostics_get_descriptor%published_fields(1)="vwp"
    scalar_diagnostics_get_descriptor%published_fields(2)="lwp"
    scalar_diagnostics_get_descriptor%published_fields(3)="qlmax"
    scalar_diagnostics_get_descriptor%published_fields(4)="hqlmax"
    scalar_diagnostics_get_descriptor%published_fields(5)="cltop"
    scalar_diagnostics_get_descriptor%published_fields(6)="clbas"
    scalar_diagnostics_get_descriptor%published_fields(7)="wmax"
    scalar_diagnostics_get_descriptor%published_fields(8)="wmin"
    scalar_diagnostics_get_descriptor%published_fields(9)="senhf"
    scalar_diagnostics_get_descriptor%published_fields(10)="lathf"
    scalar_diagnostics_get_descriptor%published_fields(11)="rwp"
    scalar_diagnostics_get_descriptor%published_fields(12)="iwp"
    scalar_diagnostics_get_descriptor%published_fields(13)="swp"
    scalar_diagnostics_get_descriptor%published_fields(14)="gwp"
    scalar_diagnostics_get_descriptor%published_fields(15)="tot_iwp"
    scalar_diagnostics_get_descriptor%published_fields(16)="reske"
    
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
    allocate(wmax(y_size_local, x_size_local), wmin(y_size_local, x_size_local), &
         reske(y_size_local, x_size_local),                                      &
         senhf(y_size_local, x_size_local),lathf(y_size_local, x_size_local))
    ! allocate the 1d velocity arrays for the kinetic energy calc
    allocate(ww_prime_res(current_state%local_grid%size(Z_INDEX)), &
         uu_prime_res(current_state%local_grid%size(Z_INDEX)),     &
         vv_prime_res(current_state%local_grid%size(Z_INDEX)))
    ww_prime_res(:) = 0.0
    uu_prime_res(:) = 0.0
    vv_prime_res(:) = 0.0

    if (.not. current_state%passive_q .and. current_state%number_q_fields .gt. 0) then
       iqv = current_state%water_vapour_mixing_ratio_index
       iql = current_state%liquid_water_mixing_ratio_index
       allocate(vwp(y_size_local, x_size_local), lwp(y_size_local, x_size_local), &
            qlmax(y_size_local, x_size_local), hqlmax(y_size_local, x_size_local), &
            cltop(y_size_local, x_size_local), clbas(y_size_local, x_size_local))
       ! allocate other hydrometeors. Allocation dependent on index being set in
       ! appropriate microphysics scheme (see casim component from example)
       if (current_state%rain_water_mixing_ratio_index > 0) then
          iqr = current_state%rain_water_mixing_ratio_index
          allocate(rwp(y_size_local, x_size_local))
       endif
       if (current_state%ice_water_mixing_ratio_index > 0) then
          iqi = current_state%ice_water_mixing_ratio_index 
          allocate(iwp(y_size_local, x_size_local))
          allocate(tot_iwp(y_size_local, x_size_local))
       endif
       if (current_state%snow_water_mixing_ratio_index > 0) then
          iqs = current_state%snow_water_mixing_ratio_index
          allocate(swp(y_size_local, x_size_local))
       endif
       if (current_state%graupel_water_mixing_ratio_index > 0) then
          iqg = current_state%graupel_water_mixing_ratio_index
          allocate(gwp(y_size_local, x_size_local))
       endif
    endif
    
    allocate(dz_rhon_fac(current_state%local_grid%size(Z_INDEX)))    
    do k=2, current_state%local_grid%size(Z_INDEX)
       ! used in the water path calculation
       dz_rhon_fac(k)=current_state%global_grid%configuration%vertical%dz(k)*&
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
       ! maximum vertical velocity for each column
       wmax(:,:)=0.0
       ! minimum vertical velocity for each column
       wmin(:,:)=0.0
       ! resolved ke
       reske(:,:) = 0.0
       
       if (.not. current_state%passive_q .and. current_state%number_q_fields .gt. 0) then
          ! maximum liquid water content in a column
          qlmax(:,:)=0.0
          ! the height of the maximum liquid water content in a column
          hqlmax(:,:)=0.0
          ! cloud top height where liqud water content is greater than qlcrit
          cltop(:,:)=0.0
          ! minimum cloud base where liquid water content is greater than q;crit
          clbas(:,:)=0.0
          ! water vapour path for each column
          vwp(:,:)=0.0
          ! liquid water path for each column
          lwp(:,:)=0.0
          ! rain water path for each column
          if (current_state%rain_water_mixing_ratio_index > 0) rwp(:,:)=0.0
          ! ice water path for each column
          if (current_state%ice_water_mixing_ratio_index > 0) then
             iwp(:,:)=0.0
             ! total ice water path (iwp + swp + gwp) for each column
             tot_iwp(:,:)=0.0
          endif
          ! snow water path for each column
          if (current_state%snow_water_mixing_ratio_index > 0) swp(:,:)=0.0
          ! graupel water path for each column
          if (current_state%graupel_water_mixing_ratio_index > 0) gwp(:,:)=0.0
       endif
     
       ! surface sensible heat flux
       senhf(:,:)=0.0
       ! surface latent heat flux
       lathf(:,:)=0.0
    end if
    
    if (.not. current_state%halo_column) then
       ! maximum and minimum vertical velocity in each column
       wmax(target_y_index, target_x_index)=maxval(current_state%w%data(:, current_y_index, current_x_index))
       wmin(target_y_index, target_x_index)=minval(current_state%w%data(:, current_y_index, current_x_index))
               
       ! work out the column resolved ww, uu, vv
       ww_prime_res(:) = &
            (current_state%w%data(:,current_state%column_local_y,current_state%column_local_x)**2.)
       if (allocated(current_state%global_grid%configuration%vertical%olubar)) then
          uu_prime_res(:) = &
               ((current_state%u%data(:,current_state%column_local_y,current_state%column_local_x) &
               - (current_state%global_grid%configuration%vertical%olubar(:) - current_state%ugal))**2.)
       else
          uu_prime_res(:) = 0.0
       endif
       if (allocated(current_state%global_grid%configuration%vertical%olvbar)) then
          vv_prime_res(:) = &
               ((current_state%v%data(:,current_state%column_local_y,current_state%column_local_x) &
               - (current_state%global_grid%configuration%vertical%olvbar(:) - current_state%vgal))**2.)
       else
          vv_prime_res(:) = 0.0
       endif
       ! use column resolved ww, uu, vv to derive total resolved KE for column
       do k=2, current_state%local_grid%size(Z_INDEX)-1
          reske(target_y_index, target_x_index) = reske(target_y_index, target_x_index) + &
               (uu_prime_res(k) + uu_prime_res(k+1)                                       &
               + vv_prime_res(k) + vv_prime_res(k+1)                                      &
               + 2_DEFAULT_PRECISION * ww_prime_res(k))*0.25_DEFAULT_PRECISION*           &
               current_state%global_grid%configuration%vertical%dzn(k+1)*                 &
               current_state%global_grid%configuration%vertical%rho(k)
       enddo
       ! divide reske by altitude to make it column mean (as in the LEM)
       reske(target_y_index, target_x_index) = reske(target_y_index, target_x_index)/ &
            current_state%global_grid%configuration%vertical%z(current_state%local_grid%size(Z_INDEX))
       ! subke is derive in the subgrid_profile_diagnostics component
       
       if (.not. current_state%passive_q .and. current_state%number_q_fields .gt. 0) then
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
          endif
          !
          ! calculate the vapour and liquid water path
          !
          if (current_state%water_vapour_mixing_ratio_index .gt. 0 .and. &
               current_state%number_q_fields .ge. current_state%water_vapour_mixing_ratio_index) then
             do k = 2, current_state%local_grid%size(Z_INDEX)
                vwp(target_y_index, target_x_index)=vwp(target_y_index, target_x_index) &
                     +dz_rhon_fac(k)*current_state%q(current_state%water_vapour_mixing_ratio_index)%data(k, &
                     current_y_index, current_x_index)
                lwp(target_y_index, target_x_index)=lwp(target_y_index, target_x_index) &
                     +dz_rhon_fac(k)*current_state%q(current_state%liquid_water_mixing_ratio_index)%data(k, &
                     current_y_index, current_x_index)
             enddo
             if (current_state%rain_water_mixing_ratio_index > 0) then
                do k = 2, current_state%local_grid%size(Z_INDEX)
                   rwp(target_y_index, target_x_index)=rwp(target_y_index, target_x_index) &
                        +dz_rhon_fac(k)*current_state%q(iqr)%data(k,current_y_index, current_x_index)
                enddo
             endif
             if (current_state%ice_water_mixing_ratio_index > 0) then
                do k = 2, current_state%local_grid%size(Z_INDEX)
                   iwp(target_y_index, target_x_index)=iwp(target_y_index, target_x_index) &
                        +dz_rhon_fac(k)*current_state%q(iqi)%data(k,current_y_index, current_x_index)
                enddo
                tot_iwp(target_y_index, target_x_index)=tot_iwp(target_y_index, target_x_index)+ &
                     iwp(target_y_index, target_x_index)
             endif
             if (current_state%snow_water_mixing_ratio_index > 0) then
                do k = 2, current_state%local_grid%size(Z_INDEX)
                   swp(target_y_index, target_x_index)=swp(target_y_index, target_x_index) &
                        +dz_rhon_fac(k)*current_state%q(iqs)%data(k,current_y_index, current_x_index)
                enddo
                tot_iwp(target_y_index, target_x_index)=tot_iwp(target_y_index, target_x_index)+ &
                     swp(target_y_index, target_x_index)
             endif
             if (current_state%graupel_water_mixing_ratio_index > 0) then
                do k = 2, current_state%local_grid%size(Z_INDEX)
                   gwp(target_y_index, target_x_index)=gwp(target_y_index, target_x_index) &
                        +dz_rhon_fac(k)*current_state%q(iqg)%data(k,current_y_index, current_x_index)
                enddo
                tot_iwp(target_y_index, target_x_index)=tot_iwp(target_y_index, target_x_index)+ &
                     gwp(target_y_index, target_x_index)
             endif
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

    if (name .eq. "senhf") then
      field_information%enabled=current_state%use_surface_boundary_conditions .and. current_state%th%active
    else if (name .eq. "lathf") then
      field_information%enabled=current_state%use_surface_boundary_conditions .and. &
           current_state%water_vapour_mixing_ratio_index .gt. 0 .and. &
           current_state%number_q_fields .ge. current_state%water_vapour_mixing_ratio_index
   else if (name .eq. "qlmax".or. name .eq. "cltop" .or. name .eq. "clbas") then
      field_information%enabled=.not. current_state%passive_q .and. current_state%liquid_water_mixing_ratio_index .gt. 0 &
           .and. current_state%number_q_fields .ge. current_state%liquid_water_mixing_ratio_index
    else if (name .eq. "vwp" .or. name .eq. "lwp") then
      field_information%enabled=current_state%number_q_fields .gt. 0 .and. current_state%water_vapour_mixing_ratio_index .gt. 0 &
           .and. current_state%number_q_fields .ge. current_state%water_vapour_mixing_ratio_index
   else if (name .eq. "rwp" ) then
      field_information%enabled= current_state%rain_water_mixing_ratio_index .gt. 0
   else if (name .eq. "iwp" .or. name .eq. 'tot_iwp') then
      field_information%enabled= current_state%ice_water_mixing_ratio_index .gt. 0
   else if (name .eq. "swp" ) then
      field_information%enabled= current_state%snow_water_mixing_ratio_index .gt. 0
   else if (name .eq. "gwp" ) then
      field_information%enabled= current_state%graupel_water_mixing_ratio_index .gt. 0
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

    if (name .eq. "wmax") then
      allocate(field_value%real_2d_array(current_state%local_grid%size(Y_INDEX), &
           current_state%local_grid%size(X_INDEX)))
       field_value%real_2d_array(:,:)=wmax(:,:)
    else if (name .eq. "wmin") then
       allocate(field_value%real_2d_array(current_state%local_grid%size(Y_INDEX), &
           current_state%local_grid%size(X_INDEX))) 
       field_value%real_2d_array(:,:)=wmin(:,:)
    else if (name .eq. 'reske') then
       allocate(field_value%real_2d_array(current_state%local_grid%size(Y_INDEX), &
           current_state%local_grid%size(X_INDEX)))
       field_value%real_2d_array(:,:)=reske(:,:)
    else if (name .eq. "qlmax") then
       allocate(field_value%real_2d_array(current_state%local_grid%size(Y_INDEX), &
           current_state%local_grid%size(X_INDEX))) 
       field_value%real_2d_array(:,:)=qlmax(:,:)
    !else if (name .eq. "hqlmax") then
    !   allocate(field_value%real_2d_array(current_state%local_grid%size(Y_INDEX), &
    !       current_state%local_grid%size(X_INDEX))) 
    !   field_value%real_2d_array(:,:)=hqlmax(:,:)
    else if (name .eq. "cltop") then
       allocate(field_value%real_2d_array(current_state%local_grid%size(Y_INDEX), &
           current_state%local_grid%size(X_INDEX))) 
       field_value%real_2d_array(:,:)=cltop(:,:)
    else if (name .eq. "clbas") then
       allocate(field_value%real_2d_array(current_state%local_grid%size(Y_INDEX), &
           current_state%local_grid%size(X_INDEX))) 
       field_value%real_2d_array(:,:)=clbas(:,:)
    else if (name .eq. "vwp") then
       allocate(field_value%real_2d_array(current_state%local_grid%size(Y_INDEX), &
           current_state%local_grid%size(X_INDEX))) 
       field_value%real_2d_array(:,:)=vwp(:,:)
    else if (name .eq. "lwp") then
       allocate(field_value%real_2d_array(current_state%local_grid%size(Y_INDEX), &
           current_state%local_grid%size(X_INDEX))) 
       field_value%real_2d_array(:,:)=lwp(:,:)
    else if (name .eq. "rwp") then
       allocate(field_value%real_2d_array(current_state%local_grid%size(Y_INDEX), &
           current_state%local_grid%size(X_INDEX))) 
       field_value%real_2d_array(:,:)=rwp(:,:)
    else if (name .eq. "iwp") then
       allocate(field_value%real_2d_array(current_state%local_grid%size(Y_INDEX), &
           current_state%local_grid%size(X_INDEX))) 
       field_value%real_2d_array(:,:)=iwp(:,:)
     else if (name .eq. "swp") then
       allocate(field_value%real_2d_array(current_state%local_grid%size(Y_INDEX), &
           current_state%local_grid%size(X_INDEX))) 
       field_value%real_2d_array(:,:)=swp(:,:)
    else if (name .eq. "gwp") then
       allocate(field_value%real_2d_array(current_state%local_grid%size(Y_INDEX), &
           current_state%local_grid%size(X_INDEX))) 
       field_value%real_2d_array(:,:)=gwp(:,:)
    else if (name .eq. "tot_iwp") then
       allocate(field_value%real_2d_array(current_state%local_grid%size(Y_INDEX), &
           current_state%local_grid%size(X_INDEX))) 
       field_value%real_2d_array(:,:)=tot_iwp(:,:)
    else if (name .eq. "senhf") then
       allocate(field_value%real_2d_array(current_state%local_grid%size(Y_INDEX), &
           current_state%local_grid%size(X_INDEX))) 
       field_value%real_2d_array(:,:)=senhf(:,:)
    else if (name .eq. "lathf") then
       allocate(field_value%real_2d_array(current_state%local_grid%size(Y_INDEX), &
           current_state%local_grid%size(X_INDEX))) 
       field_value%real_2d_array(:,:)=lathf(:,:) 
    end if
    
  end subroutine field_value_retrieval_callback
end module scalar_diagnostics_mod
