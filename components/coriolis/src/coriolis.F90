!> This calculates the coriolis and mean pressure gradient terms which impact su and sv fields
module coriolis_mod
  use datadefn_mod, only : DEFAULT_PRECISION
  use monc_component_mod, only : COMPONENT_ARRAY_FIELD_TYPE, &
      COMPONENT_DOUBLE_DATA_TYPE, component_descriptor_type, &
      component_field_value_type, component_field_information_type
  use optionsdatabase_mod, only : options_get_logical, options_get_real, options_get_integer
  use state_mod, only : model_state_type
  use grids_mod, only : Z_INDEX, Y_INDEX, X_INDEX
  implicit none

#ifndef TEST_MODE
  private
#endif

  logical :: baroclinicity_use_geostrophic_shear
  real(kind=DEFAULT_PRECISION), dimension(:), allocatable :: geostrophic_wind_x, geostrophic_wind_y
  real(kind=DEFAULT_PRECISION) :: fcoriol

  ! Local tendency diagnostic variables for this component
  ! 3D tendency fields and logicals for their use
  real(kind=DEFAULT_PRECISION), dimension(:,:,:), allocatable ::     &
       tend_3d_u, tend_3d_v
  logical :: l_tend_3d_u, l_tend_3d_v
  ! Local mean tendency profile fields and logicals for their use
  real(kind=DEFAULT_PRECISION), dimension(:), allocatable ::         &
       tend_pr_tot_u, tend_pr_tot_v
  logical :: l_tend_pr_tot_u, l_tend_pr_tot_v

  integer :: diagnostic_generation_frequency

  public coriolis_get_descriptor

  contains

  !> Provides the descriptor back to the caller and is used in component registration
  !! @returns The termination check component descriptor
  type(component_descriptor_type) function coriolis_get_descriptor()
    coriolis_get_descriptor%name="coriolis"
    coriolis_get_descriptor%version=0.1
    coriolis_get_descriptor%initialisation=>initialisation_callback
    coriolis_get_descriptor%timestep=>timestep_callback
    coriolis_get_descriptor%finalisation=>finalisation_callback

    coriolis_get_descriptor%field_value_retrieval=>field_value_retrieval_callback
    coriolis_get_descriptor%field_information_retrieval=>field_information_retrieval_callback
    allocate(coriolis_get_descriptor%published_fields(2+2))

    coriolis_get_descriptor%published_fields(1)= "tend_u_coriolis_3d_local"
    coriolis_get_descriptor%published_fields(2)= "tend_v_coriolis_3d_local"

    coriolis_get_descriptor%published_fields(2+1)= "tend_u_coriolis_profile_total_local"
    coriolis_get_descriptor%published_fields(2+2)= "tend_v_coriolis_profile_total_local"

  end function coriolis_get_descriptor


  !> Field information retrieval callback, this returns information for a specific components published field
  !! @param current_state Current model state
  !! @param name The name of the field to retrieve information for
  !! @param field_information Populated with information about the field
  subroutine field_information_retrieval_callback(current_state, name, field_information)
    type(model_state_type), target, intent(inout) :: current_state
    character(len=*), intent(in) :: name
    type(component_field_information_type), intent(out) :: field_information
    integer :: strcomp

    ! Field information for 3d
    strcomp=INDEX(name, "coriolis_3d_local")
    if (strcomp .ne. 0) then
      field_information%field_type=COMPONENT_ARRAY_FIELD_TYPE
      field_information%number_dimensions=3
      field_information%dimension_sizes(1)=current_state%local_grid%size(Z_INDEX)
      field_information%dimension_sizes(2)=current_state%local_grid%size(Y_INDEX)
      field_information%dimension_sizes(3)=current_state%local_grid%size(X_INDEX)
      field_information%data_type=COMPONENT_DOUBLE_DATA_TYPE

      if      (name .eq. "tend_u_coriolis_3d_local") then
        field_information%enabled=l_tend_3d_u
      else if (name .eq. "tend_v_coriolis_3d_local") then
        field_information%enabled=l_tend_3d_v
      else
        field_information%enabled=.true.
      end if

    end if !end 3d check

    ! Field information for profiles
    strcomp=INDEX(name, "coriolis_profile_total_local")
    if (strcomp .ne. 0) then
      field_information%field_type=COMPONENT_ARRAY_FIELD_TYPE
      field_information%number_dimensions=1
      field_information%dimension_sizes(1)=current_state%local_grid%size(Z_INDEX)
      field_information%data_type=COMPONENT_DOUBLE_DATA_TYPE

      if      (name .eq. "tend_u_coriolis_profile_total_local") then
        field_information%enabled=l_tend_pr_tot_u
      else if (name .eq. "tend_v_coriolis_profile_total_local") then
        field_information%enabled=l_tend_pr_tot_v
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
    if      (name .eq. "tend_u_coriolis_3d_local" .and. allocated(tend_3d_u)) then
      call set_published_field_value(field_value, real_3d_field=tend_3d_u)
    else if (name .eq. "tend_v_coriolis_3d_local" .and. allocated(tend_3d_v)) then
      call set_published_field_value(field_value, real_3d_field=tend_3d_v)

    ! Profile Tendency Fields
    else if (name .eq. "tend_u_coriolis_profile_total_local" .and. allocated(tend_pr_tot_u)) then
      call set_published_field_value(field_value, real_1d_field=tend_pr_tot_u)
    else if (name .eq. "tend_v_coriolis_profile_total_local" .and. allocated(tend_pr_tot_v)) then
      call set_published_field_value(field_value, real_1d_field=tend_pr_tot_v)
    end if


  end subroutine field_value_retrieval_callback


  !> Initialisation call back which will read in the coriolis configuration and set up the geostrophic winds
  !! @param current_state The current model state
  subroutine initialisation_callback(current_state)
    type(model_state_type), target, intent(inout) :: current_state

    integer :: k

    baroclinicity_use_geostrophic_shear=options_get_logical(current_state%options_database, "baroclinicity_use_geostrophic_shear")
    fcoriol=options_get_real(current_state%options_database, "fcoriol")
    current_state%geostrophic_wind_rate_of_change_in_x=options_get_real(current_state%options_database, &
         "geostrophic_wind_rate_of_change_in_x")
    current_state%geostrophic_wind_rate_of_change_in_y=options_get_real(current_state%options_database, &
         "geostrophic_wind_rate_of_change_in_y")
    current_state%surface_geostrophic_wind_x=options_get_real(current_state%options_database, "surface_geostrophic_wind_x")
    current_state%surface_geostrophic_wind_y=options_get_real(current_state%options_database, "surface_geostrophic_wind_y")

    allocate(geostrophic_wind_x(current_state%local_grid%size(Z_INDEX)), &
         geostrophic_wind_y(current_state%local_grid%size(Z_INDEX)))

    do k=1,current_state%local_grid%size(Z_INDEX)
      geostrophic_wind_x(k)=current_state%surface_geostrophic_wind_x
      geostrophic_wind_y(k)=current_state%surface_geostrophic_wind_y
      if (baroclinicity_use_geostrophic_shear) then
        geostrophic_wind_x(k)=geostrophic_wind_x(k)+current_state%geostrophic_wind_rate_of_change_in_x*&
             current_state%global_grid%configuration%vertical%zn(k)
        geostrophic_wind_y(k)=geostrophic_wind_y(k)+current_state%geostrophic_wind_rate_of_change_in_y*&
             current_state%global_grid%configuration%vertical%zn(k)
      end if
    end do

    ! Tendency Logicals
    l_tend_pr_tot_u   = current_state%u%active
    l_tend_pr_tot_v   = current_state%v%active

    l_tend_3d_u   = current_state%u%active .or. l_tend_pr_tot_u
    l_tend_3d_v   = current_state%v%active .or. l_tend_pr_tot_v

    ! Allocate 3d tendency fields upon availability
    if (l_tend_3d_u) then
      allocate( tend_3d_u(current_state%local_grid%size(Z_INDEX),  &
                          current_state%local_grid%size(Y_INDEX),  &
                          current_state%local_grid%size(X_INDEX)   )    )
    endif
    if (l_tend_3d_v) then
      allocate( tend_3d_v(current_state%local_grid%size(Z_INDEX),  &
                          current_state%local_grid%size(Y_INDEX),  &
                          current_state%local_grid%size(X_INDEX)   )    )
    endif

    ! Allocate profile tendency fields upon availability
    if (l_tend_pr_tot_u) then
      allocate( tend_pr_tot_u(current_state%local_grid%size(Z_INDEX)) )
    endif
    if (l_tend_pr_tot_v) then
      allocate( tend_pr_tot_v(current_state%local_grid%size(Z_INDEX)) )
    endif

    ! Save the sampling_frequency to force diagnostic calculation on select time steps
    diagnostic_generation_frequency=options_get_integer(current_state%options_database, "sampling_frequency")

  end subroutine initialisation_callback


  subroutine finalisation_callback(current_state)
    type(model_state_type), target, intent(inout) :: current_state

    if (allocated(tend_3d_u)) deallocate(tend_3d_u)
    if (allocated(tend_3d_v)) deallocate(tend_3d_v)

    if (allocated(tend_pr_tot_u)) deallocate(tend_pr_tot_u)
    if (allocated(tend_pr_tot_v)) deallocate(tend_pr_tot_v)

  end subroutine finalisation_callback

  !> For each none halo cell this will calculate the coriolis terms for su and sv fields
  !! @param current_state The current model state
  subroutine timestep_callback(current_state)
    type(model_state_type), target, intent(inout) :: current_state

    integer :: local_y, locaL_x, k, target_x_index, target_y_index

    local_y=current_state%column_local_y
    local_x=current_state%column_local_x
    target_y_index=local_y-current_state%local_grid%halo_size(Y_INDEX)
    target_x_index=local_x-current_state%local_grid%halo_size(X_INDEX)

    ! Zero profile tendency totals on first instance in the sum
    if (current_state%first_timestep_column) then
      if (l_tend_pr_tot_u) then
        tend_pr_tot_u(:)= 0.0_DEFAULT_PRECISION
      endif
      if (l_tend_pr_tot_v) then
        tend_pr_tot_v(:)= 0.0_DEFAULT_PRECISION
      endif
    endif  ! zero totals

    if (current_state%halo_column) then
      if (.not. ((current_state%column_local_y == current_state%local_grid%halo_size(Y_INDEX) .and. &
           current_state%column_local_x .le. current_state%local_grid%local_domain_end_index(X_INDEX) .and. &
           current_state%column_local_x .ge. current_state%local_grid%local_domain_start_index(X_INDEX)-1) .or. &
           (current_state%column_local_x == current_state%local_grid%halo_size(X_INDEX) .and. &
           current_state%column_local_y .ge. current_state%local_grid%local_domain_start_index(Y_INDEX) &
           .and. current_state%column_local_y .le. current_state%local_grid%local_domain_end_index(Y_INDEX)) )) return
    end if

    if (mod(current_state%timestep, diagnostic_generation_frequency) == 0 .and. .not. current_state%halo_column) then
      call save_precomponent_tendencies(current_state, local_x, local_y, target_x_index, target_y_index)
    end if
    
    do k=2,current_state%local_grid%size(Z_INDEX)
#if defined(U_ACTIVE) && defined(V_ACTIVE)
      current_state%su%data(k, current_state%column_local_y, current_state%column_local_x)=&
           current_state%su%data(k, current_state%column_local_y, current_state%column_local_x)+fcoriol*&
           (0.25_DEFAULT_PRECISION*(current_state%v%data(k, current_state%column_local_y, current_state%column_local_x)+&
           current_state%v%data(k, current_state%column_local_y, current_state%column_local_x+1)+&
           current_state%v%data(k, current_state%column_local_y-1, current_state%column_local_x)+&
           current_state%v%data(k, current_state%column_local_y-1, current_state%column_local_x+1))+current_state%vgal-&
           geostrophic_wind_y(k))

      current_state%sv%data(k, current_state%column_local_y, current_state%column_local_x)=&
           current_state%sv%data(k, current_state%column_local_y, current_state%column_local_x)-fcoriol*&
           (0.25_DEFAULT_PRECISION*(current_state%u%data(k, current_state%column_local_y, current_state%column_local_x)+&
           current_state%u%data(k, current_state%column_local_y, current_state%column_local_x-1)+&
           current_state%u%data(k, current_state%column_local_y+1, current_state%column_local_x)+&
           current_state%u%data(k, current_state%column_local_y+1, current_state%column_local_x-1))+current_state%ugal-&
           geostrophic_wind_x(k))
#endif
    end do

    if (mod(current_state%timestep, diagnostic_generation_frequency) == 0  .and. .not. current_state%halo_column) then
      call compute_component_tendencies(current_state, local_x, local_y, target_x_index, target_y_index)
    end if

  end subroutine timestep_callback


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
    if (l_tend_3d_u) then
      tend_3d_u(:,tyn,txn)=current_state%su%data(:,cyn,cxn)
    endif
    if (l_tend_3d_v) then
      tend_3d_v(:,tyn,txn)=current_state%sv%data(:,cyn,cxn)
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
    if (l_tend_3d_u) then
      tend_3d_u(:,tyn,txn)=current_state%su%data(:,cyn,cxn) - tend_3d_u(:,tyn,txn)
    endif
    if (l_tend_3d_v) then
      tend_3d_v(:,tyn,txn)=current_state%sv%data(:,cyn,cxn) - tend_3d_v(:,tyn,txn)
    endif

   ! Add local tendency fields to the profile total
    if (l_tend_pr_tot_u) then
      tend_pr_tot_u(:)=tend_pr_tot_u(:) + tend_3d_u(:,tyn,txn)
    endif
    if (l_tend_pr_tot_v) then
      tend_pr_tot_v(:)=tend_pr_tot_v(:) + tend_3d_v(:,tyn,txn)
    endif

  end subroutine compute_component_tendencies

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


end module coriolis_mod
