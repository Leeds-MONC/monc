!> Specific theta advection, which involves the vertical advection of reference state and advection of mean baroclinicity
module thadvection_mod
  use datadefn_mod, only : DEFAULT_PRECISION
  use monc_component_mod, only : COMPONENT_ARRAY_FIELD_TYPE, &
      COMPONENT_DOUBLE_DATA_TYPE, component_descriptor_type, &
      component_field_value_type, component_field_information_type
  use state_mod, only : model_state_type
  use grids_mod, only : Z_INDEX, Y_INDEX, X_INDEX
  use science_constants_mod, only : G
  use logging_mod, only : LOG_ERROR, log_master_log
  use optionsdatabase_mod, only : options_get_real, options_get_integer, options_get_logical
implicit none

#ifndef TEST_MODE
  private
#endif

  logical :: baroclinicity_use_geostrophic_shear
  real(kind=DEFAULT_PRECISION) :: fcoriol, fcoriol_over_G, rate_change_geostrophic_wind_x, rate_change_geostrophic_wind_y, &
       multiplicative_factor_x, multiplicative_factor_y

  ! Local tendency diagnostic variables for this component
  ! 3D tendency fields and logicals for their use
  real(kind=DEFAULT_PRECISION), dimension(:,:,:), allocatable :: tend_3d_th, tend_3d_tabs
  logical :: l_tend_3d_th, l_tend_3d_tabs
  ! Local mean tendency profile fields and logicals for their use
  real(kind=DEFAULT_PRECISION), dimension(:), allocatable ::     tend_pr_tot_th, tend_pr_tot_tabs
  logical :: l_tend_pr_tot_th, l_tend_pr_tot_tabs

  integer :: diagnostic_generation_frequency


  public thadvection_get_descriptor
contains

  !> Provides a description of this component for the core to register
  !! @returns The descriptor containing registration information for this component
  type(component_descriptor_type) function thadvection_get_descriptor()
    thadvection_get_descriptor%name="th_advection"
    thadvection_get_descriptor%version=0.1
    thadvection_get_descriptor%initialisation=>initialisation_callback
    thadvection_get_descriptor%timestep=>timestep_callback
    thadvection_get_descriptor%finalisation=>finalisation_callback

    thadvection_get_descriptor%field_value_retrieval=>field_value_retrieval_callback
    thadvection_get_descriptor%field_information_retrieval=>field_information_retrieval_callback
    allocate(thadvection_get_descriptor%published_fields(2+2))

    thadvection_get_descriptor%published_fields(1)= "tend_th_thadvection_3d_local"
    thadvection_get_descriptor%published_fields(2)= "tend_tabs_thadvection_3d_local"

    thadvection_get_descriptor%published_fields(2+1)= "tend_th_thadvection_profile_total_local"
    thadvection_get_descriptor%published_fields(2+2)= "tend_tabs_thadvection_profile_total_local"

  end function thadvection_get_descriptor

  !> Initialisation callback to set up the variables and data needed by the component
  !! @param current_state The current model state
  subroutine initialisation_callback(current_state)
    type(model_state_type), target, intent(inout) :: current_state

    baroclinicity_use_geostrophic_shear=options_get_logical(current_state%options_database, "baroclinicity_use_geostrophic_shear")
    fcoriol=options_get_real(current_state%options_database, "fcoriol")
    rate_change_geostrophic_wind_x=options_get_real(current_state%options_database, "rate_change_geostrophic_wind_x")
    rate_change_geostrophic_wind_y=options_get_real(current_state%options_database, "rate_change_geostrophic_wind_y")
    fcoriol_over_G = fcoriol/G
    multiplicative_factor_x=rate_change_geostrophic_wind_x*current_state%thref0*fcoriol_over_G
    multiplicative_factor_y=rate_change_geostrophic_wind_y*current_state%thref0*fcoriol_over_G

    ! Set tendency diagnostic logicals based on availability
    ! Need to use 3d tendencies to compute the profiles, so they will be allocated
    !      in the case where profiles are available
    l_tend_pr_tot_th  = current_state%th%active 
    l_tend_pr_tot_tabs  = l_tend_pr_tot_th

    l_tend_3d_th  = current_state%th%active .or. l_tend_pr_tot_th
    l_tend_3d_tabs = l_tend_3d_th

    ! Allocate 3d tendency fields upon availability
    if (l_tend_3d_th) then
      allocate( tend_3d_th(current_state%local_grid%size(Z_INDEX),  &
                           current_state%local_grid%size(Y_INDEX),  &
                           current_state%local_grid%size(X_INDEX)   )    )
    endif
    if (l_tend_3d_tabs) then
      allocate( tend_3d_tabs(current_state%local_grid%size(Z_INDEX),  &
                             current_state%local_grid%size(Y_INDEX),  &
                             current_state%local_grid%size(X_INDEX)   )    )
    endif

    ! Allocate profile tendency fields upon availability
    if (l_tend_pr_tot_th) then
      allocate( tend_pr_tot_th(current_state%local_grid%size(Z_INDEX)) )
    endif
    if (l_tend_pr_tot_tabs) then
      allocate( tend_pr_tot_tabs(current_state%local_grid%size(Z_INDEX)) )
    endif

    ! Save the sampling_frequency to force diagnostic calculation on select time steps
    diagnostic_generation_frequency=options_get_integer(current_state%options_database, "sampling_frequency")

  end subroutine initialisation_callback  


  subroutine finalisation_callback(current_state)
    type(model_state_type), target, intent(inout) :: current_state

    if (allocated(tend_3d_th)) deallocate(tend_3d_th)
    if (allocated(tend_3d_tabs)) deallocate(tend_3d_tabs)

    if (allocated(tend_pr_tot_th)) deallocate(tend_pr_tot_th)
    if (allocated(tend_pr_tot_tabs)) deallocate(tend_pr_tot_tabs)

  end subroutine finalisation_callback


  !> Timestep callback, will call the two separate procedures to do their advection if needed
  !! @param current_state The current model state
  !! @param target_(x/y)_index This is the index with the halos subtracted. This is needed so that diagnostic does
  !!                           not include halos and to prevent array out-of-bounds
  subroutine timestep_callback(current_state)
    type(model_state_type), target, intent(inout) :: current_state

    integer :: current_x_index, current_y_index, target_x_index, target_y_index

    current_x_index=current_state%column_local_x
    current_y_index=current_state%column_local_y
    target_y_index=current_y_index-current_state%local_grid%halo_size(Y_INDEX)
    target_x_index=current_x_index-current_state%local_grid%halo_size(X_INDEX)

    ! Zero profile tendency totals on first instance in the sum
    if (current_state%first_timestep_column) then
      if (l_tend_pr_tot_th) then
        tend_pr_tot_th(:)=0.0_DEFAULT_PRECISION
      endif
      if (l_tend_pr_tot_tabs) then
        tend_pr_tot_tabs(:)=0.0_DEFAULT_PRECISION
      endif
    endif  ! zero totals

    if (mod(current_state%timestep, diagnostic_generation_frequency) == 0 .and. .not. current_state%halo_column) then
      call save_precomponent_tendencies(current_state, current_x_index, current_y_index, target_x_index, target_y_index)
    end if

    call vertical_advection_of_reference_state(current_state, current_state%column_local_y, current_state%column_local_x)
    call advection_of_mean_baroclinicity(current_state, current_state%column_local_y, current_state%column_local_x)

    if (mod(current_state%timestep, diagnostic_generation_frequency) == 0 .and. .not. current_state%halo_column) then
      call compute_component_tendencies(current_state, current_x_index, current_y_index, target_x_index, target_y_index)
    end if

  end subroutine timestep_callback

  !> Vertical advection of the reference state.  It doesn't seem consistent to do the advection in this way if
  !! TVD advection of the deviation from the reference state has been selected. Separate vertical advection of the reference
  !! state was introduced to improve energy conservation when carrying out idealized gravity wave simulations in a deep, dry
  !! isothermal layer, for which the difference in potential temp between top and bottom was of order 100K. In less extreme cases
  !! the benefits are unlikely to be significant and with TVD advection energy conservation has been compromised so the best
  !! way forward might be to recombine the reference state into l_th
  !! @param current_state The current model state
  !! @param local_y The local Y of the column
  !! @param local_x The local X of the column
  subroutine vertical_advection_of_reference_state(current_state, local_y, local_x)
    type(model_state_type), target, intent(inout) :: current_state
    integer, intent(in) :: local_y, local_x

    integer :: k
    real(kind=DEFAULT_PRECISION) :: sctmp1, sctmp2

    if (current_state%use_anelastic_equations) then
      ! This code only needs to be executed if anelastic, otherwise THREF is constant and the increment calculated here is zero
      do k=2, current_state%local_grid%size(Z_INDEX)
        sctmp1=current_state%global_grid%configuration%vertical%tzc1(k)*2.0_DEFAULT_PRECISION*&
             current_state%global_grid%configuration%vertical%dthref(k-1)
        sctmp2=current_state%global_grid%configuration%vertical%tzc2(k)*2.0_DEFAULT_PRECISION*&
             current_state%global_grid%configuration%vertical%dthref(k)
        current_state%sth%data(k, local_y, local_x)=current_state%sth%data(k, local_y, local_x)-(sctmp1*&
             current_state%w%data(k-1, local_y, local_x) + sctmp2*current_state%w%data(k, local_y, local_x))
      end do
    end if
  end subroutine vertical_advection_of_reference_state

  !> Performs advection of the mean baroclinicity if appropriate
  !! @param current_state The current model state
  !! @param local_y The local Y of the column
  !! @param local_x The local X of the column
  subroutine advection_of_mean_baroclinicity(current_state, local_y, local_x)
    type(model_state_type), target, intent(inout) :: current_state
    integer, intent(in) :: local_y, local_x

    integer :: k

    if (baroclinicity_use_geostrophic_shear) then
        if (current_state%passive_q) then
            if (current_state%use_anelastic_equations) then
                do k=2, current_state%local_grid%size(Z_INDEX)
                  current_state%sth%data(k, local_y, local_x)=current_state%sth%data(k, local_y, local_x)+&
                       current_state%global_grid%configuration%vertical%thref(k)*fcoriol_over_G*&
                       ((current_state%v%data(k, local_y, local_x) + current_state%vgal) * rate_change_geostrophic_wind_x-&
                       (current_state%u%data(k, local_y, local_x) + current_state%ugal) * rate_change_geostrophic_wind_y)
                end do
              else
                do k=2, current_state%local_grid%size(Z_INDEX)
                  current_state%sth%data(k, local_y, local_x)=current_state%sth%data(k, local_y, local_x)+&
                       ((current_state%v%data(k, local_y, local_x) + current_state%vgal) * multiplicative_factor_x-&
                       (current_state%u%data(k, local_y, local_x) + current_state%ugal) * multiplicative_factor_y)                    
                end do
              end if
        else
          call log_master_log(LOG_ERROR, "The combination if baroclinicity and active q is not yet allowed")
        end if
      end if
  end subroutine advection_of_mean_baroclinicity  


   !> Save the 3d tendencies coming into the component.
  !! @param current_state Current model state
  !! @param cxn The current slice, x, index
  !! @param cyn The current column, y, index.
  !! @param txn target_x_index
  !! @param tyn target_y_index
  subroutine save_precomponent_tendencies(current_state, cxn, cyn, txn, tyn)
    type(model_state_type), target, intent(in) :: current_state
    integer, intent(in) ::  cxn, cyn, txn, tyn

    ! Save 3d tendency fields upon request (of 3d or profiles) and availability
    if (l_tend_3d_th) then
      tend_3d_th(:,tyn,txn)=current_state%sth%data(:,cyn,cxn)
    endif
    if (l_tend_3d_tabs) then
      tend_3d_tabs(:,tyn,txn)=current_state%sth%data(:,cyn,cxn) * current_state%global_grid%configuration%vertical%rprefrcp(:)
    endif

  end subroutine save_precomponent_tendencies


   !> Computation of component tendencies
  !! @param current_state Current model state
  !! @param cxn The current slice, x, index
  !! @param cyn The current column, y, index.
  !! @param txn target_x_index
  !! @param tyn target_y_index
  subroutine compute_component_tendencies(current_state, cxn, cyn, txn, tyn)
    type(model_state_type), target, intent(inout) :: current_state
    integer, intent(in) ::  cxn, cyn, txn, tyn

    ! Calculate change in tendency due to component
    if (l_tend_3d_th) then
      tend_3d_th(:,tyn,txn)=current_state%sth%data(:,cyn,cxn)     - tend_3d_th(:,tyn,txn)
    endif
    if (l_tend_3d_tabs) then
      tend_3d_tabs(:,tyn,txn)=   &
       current_state%sth%data(:,cyn,cxn) * current_state%global_grid%configuration%vertical%rprefrcp(:)   &
        - tend_3d_tabs(:,tyn,txn)
    endif

   ! Add local tendency fields to the profile total
    if (l_tend_pr_tot_th) then
      tend_pr_tot_th(:)=tend_pr_tot_th(:) + tend_3d_th(:,tyn,txn)
    endif
    if (l_tend_pr_tot_tabs) then
      tend_pr_tot_tabs(:)=tend_pr_tot_tabs(:) + tend_3d_tabs(:,tyn,txn)
    endif

  end subroutine compute_component_tendencies


   !> Field information retrieval callback, this returns information for a specific component's published field
  !! @param current_state Current model state
  !! @param name The name of the field to retrieve information for
  !! @param field_information Populated with information about the field
  !! @param strcomp Starting index within 1st argument string that matches substring (2nd argument); 0 if not a match.
  subroutine field_information_retrieval_callback(current_state, name, field_information)
    type(model_state_type), target, intent(inout) :: current_state
    character(len=*), intent(in) :: name
    type(component_field_information_type), intent(out) :: field_information
    integer :: strcomp

    ! Field information for 3d
    strcomp=INDEX(name, "thadvection_3d_local")
    if (strcomp .ne. 0) then
      field_information%field_type=COMPONENT_ARRAY_FIELD_TYPE
      field_information%number_dimensions=3
      field_information%dimension_sizes(1)=current_state%local_grid%size(Z_INDEX)
      field_information%dimension_sizes(2)=current_state%local_grid%size(Y_INDEX)
      field_information%dimension_sizes(3)=current_state%local_grid%size(X_INDEX)
      field_information%data_type=COMPONENT_DOUBLE_DATA_TYPE

      if (name .eq. "tend_th_thadvection_3d_local") then
        field_information%enabled=l_tend_3d_th
      else if (name .eq. "tend_tabs_thadvection_3d_local") then
        field_information%enabled=l_tend_3d_tabs
      else
        field_information%enabled=.true.
      end if

    end if !end 3d check

    ! Field information for profiles
    strcomp=INDEX(name, "thadvection_profile_total_local")
    if (strcomp .ne. 0) then
      field_information%field_type=COMPONENT_ARRAY_FIELD_TYPE
      field_information%number_dimensions=1
      field_information%dimension_sizes(1)=current_state%local_grid%size(Z_INDEX)
      field_information%data_type=COMPONENT_DOUBLE_DATA_TYPE

      if (name .eq. "tend_th_thadvection_profile_total_local") then
        field_information%enabled=l_tend_pr_tot_th
      else if (name .eq. "tend_tabs_thadvection_profile_total_local") then
        field_information%enabled=l_tend_pr_tot_tabs
      else
        field_information%enabled=.true.
      end if

    end if !end profile check

  end subroutine field_information_retrieval_callback


    !> Field value retrieval callback, this returns the value of a specific published field
  !! @param current_state Current model state
  !! @param name The name of the field to retrieve the value for
  !! @param field_value Populated with the value of the field
  subroutine field_value_retrieval_callback(current_state, name, field_value)
    type(model_state_type), target, intent(inout) :: current_state
    character(len=*), intent(in) :: name
    type(component_field_value_type), intent(out) :: field_value

    ! 3d Tendency Fields
    if (name .eq. "tend_th_thadvection_3d_local" .and. allocated(tend_3d_th)) then
      call set_published_field_value(field_value, real_3d_field=tend_3d_th)
    else if (name .eq. "tend_tabs_thadvection_3d_local" .and. allocated(tend_3d_tabs)) then
      call set_published_field_value(field_value, real_3d_field=tend_3d_tabs)

    ! Profile Tendency Fields
    else if (name .eq. "tend_th_thadvection_profile_total_local" .and. allocated(tend_pr_tot_th)) then
      call set_published_field_value(field_value, real_1d_field=tend_pr_tot_th)
    else if (name .eq. "tend_tabs_thadvection_profile_total_local" .and. allocated(tend_pr_tot_tabs)) then
      call set_published_field_value(field_value, real_1d_field=tend_pr_tot_tabs)
    end if

  end subroutine field_value_retrieval_callback


  !> Sets the published field value from the temporary diagnostic values held by this component.
  !! @param field_value Populated with the value of the field
  !! @param real_1d_field Optional one dimensional real of values to publish
  !! @param real_2d_field Optional two dimensional real of values to publish
  subroutine set_published_field_value(field_value, real_1d_field, real_2d_field, real_3d_field)
    type(component_field_value_type), intent(inout) :: field_value
    real(kind=DEFAULT_PRECISION), dimension(:), optional :: real_1d_field
    real(kind=DEFAULT_PRECISION), dimension(:,:), optional :: real_2d_field
    real(kind=DEFAULT_PRECISION), dimension(:,:,:), optional :: real_3d_field

    if (present(real_1d_field)) then
      allocate(field_value%real_1d_array(size(real_1d_field)), source=real_1d_field)
    else if (present(real_2d_field)) then
      allocate(field_value%real_2d_array(size(real_2d_field, 1), size(real_2d_field, 2)), source=real_2d_field)
    else if (present(real_3d_field)) then
      allocate(field_value%real_3d_array(size(real_3d_field, 1), size(real_3d_field, 2), size(real_3d_field, 3)), &
               source=real_3d_field)
    end if
  end subroutine set_published_field_value


end module thadvection_mod
