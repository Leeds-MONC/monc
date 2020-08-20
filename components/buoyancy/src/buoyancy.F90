!> Calculates buoyancy terms for the SW field
module buoyancy_mod
  use datadefn_mod, only : DEFAULT_PRECISION
  use monc_component_mod, only : COMPONENT_ARRAY_FIELD_TYPE, &
      COMPONENT_DOUBLE_DATA_TYPE, component_descriptor_type, &
      component_field_value_type, component_field_information_type
  use registry_mod, only : is_component_enabled
  use optionsdatabase_mod, only : options_has_key, options_get_real_array, options_get_integer
  use state_mod, only : FORWARD_STEPPING, model_state_type
  use grids_mod, only : Z_INDEX, X_INDEX, Y_INDEX
  use science_constants_mod
  use q_indices_mod, only: get_q_index, standard_q_names
implicit none

#ifndef TEST_MODE
  private
#endif
  real(kind=DEFAULT_PRECISION), dimension(:), allocatable :: w_buoyancy

  real(kind=DEFAULT_PRECISION) :: G_over_2

  integer :: iqv ! Index for water vapour

  ! Local tendency diagnostic variables for this component
  ! 3D tendency fields and logicals for their use
  real(kind=DEFAULT_PRECISION), dimension(:,:,:), allocatable :: tend_3d_w
  logical :: l_tend_3d_w
  ! Local mean tendency profile fields and logicals for their use
  real(kind=DEFAULT_PRECISION), dimension(:), allocatable :: tend_pr_tot_w
  logical :: l_tend_pr_tot_w

  public buoyancy_get_descriptor

contains

  !> Provides the descriptor back to the caller and is used in component registration
  !! @returns The termination check component descriptor
  type(component_descriptor_type) function buoyancy_get_descriptor()
    buoyancy_get_descriptor%name="buoyancy"
    buoyancy_get_descriptor%version=0.1
    buoyancy_get_descriptor%initialisation=>initialisation_callback
    buoyancy_get_descriptor%timestep=>timestep_callback
    buoyancy_get_descriptor%finalisation=>finalisation_callback

    buoyancy_get_descriptor%field_value_retrieval=>field_value_retrieval_callback
    buoyancy_get_descriptor%field_information_retrieval=>field_information_retrieval_callback
    allocate(buoyancy_get_descriptor%published_fields(1+1+1))

    buoyancy_get_descriptor%published_fields(1)="w_buoyancy"

    buoyancy_get_descriptor%published_fields(1+1)= "tend_w_buoyancy_3d_local"

    buoyancy_get_descriptor%published_fields(1+1+1)= "tend_w_buoyancy_profile_total_local"

  end function buoyancy_get_descriptor

  !> Field information retrieval callback, this returns information for a specific components published field
  !! @param current_state Current model state
  !! @param name The name of the field to retrieve information for
  !! @param field_information Populated with information about the field
  subroutine field_information_retrieval_callback(current_state, name, field_information)
    type(model_state_type), target, intent(inout) :: current_state
    character(len=*), intent(in) :: name
    type(component_field_information_type), intent(out) :: field_information
    integer :: strcomp

    ! Field description is the same regardless of the specific field being retrieved
    field_information%field_type=COMPONENT_ARRAY_FIELD_TYPE
    field_information%data_type=COMPONENT_DOUBLE_DATA_TYPE
    field_information%number_dimensions=1
    field_information%dimension_sizes(1)=current_state%local_grid%size(Z_INDEX)
    field_information%enabled=.true.

    ! Field information for 3d
    strcomp=INDEX(name, "buoyancy_3d_local")
    if (strcomp .ne. 0) then
      field_information%field_type=COMPONENT_ARRAY_FIELD_TYPE
      field_information%number_dimensions=3
      field_information%dimension_sizes(1)=current_state%local_grid%size(Z_INDEX)
      field_information%dimension_sizes(2)=current_state%local_grid%size(Y_INDEX)
      field_information%dimension_sizes(3)=current_state%local_grid%size(X_INDEX)
      field_information%data_type=COMPONENT_DOUBLE_DATA_TYPE

      if (name .eq. "tend_w_buoyancy_3d_local") then
        field_information%enabled=l_tend_3d_w
      else
        field_information%enabled=.true.
      end if

    end if !end 3d check

    ! Field information for profiles
    strcomp=INDEX(name, "buoyancy_profile_total_local")
    if (strcomp .ne. 0) then
      field_information%field_type=COMPONENT_ARRAY_FIELD_TYPE
      field_information%number_dimensions=1
      field_information%dimension_sizes(1)=current_state%local_grid%size(Z_INDEX)
      field_information%data_type=COMPONENT_DOUBLE_DATA_TYPE

      if (name .eq. "tend_w_buoyancy_profile_total_local") then
        field_information%enabled=l_tend_pr_tot_w
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
    
    if (name .eq. "w_buoyancy") then
      allocate(field_value%real_1d_array(size(w_buoyancy)), source=w_buoyancy)   

    ! 3d Tendency Fields
    else if (name .eq. "tend_w_buoyancy_3d_local" .and. allocated(tend_3d_w)) then
      call set_published_field_value(field_value, real_3d_field=tend_3d_w)

    ! Profile Tendency Fields
    else if (name .eq. "tend_w_buoyancy_profile_total_local" .and. allocated(tend_pr_tot_w)) then
      call set_published_field_value(field_value, real_1d_field=tend_pr_tot_w)
    end if

  end subroutine field_value_retrieval_callback


  !> The initialisation callback sets up the buoyancy coefficient
  !! @param current_state The current model state
  subroutine initialisation_callback(current_state)
    type(model_state_type), target, intent(inout) :: current_state

    integer :: z_size

    if (.not. current_state%passive_q .and. current_state%number_q_fields > 0)then
      if (.not. allocated(current_state%cq))then
        allocate(current_state%cq(current_state%number_q_fields))
        current_state%cq=0.0_DEFAULT_PRECISION
      end if
      iqv = get_q_index(standard_q_names%VAPOUR, 'buoyancy')
      current_state%cq(iqv) = ratio_mol_wts-1.0
    end if

    G_over_2 = 0.5_DEFAULT_PRECISION*G
    z_size=current_state%global_grid%size(Z_INDEX)
    allocate(w_buoyancy(z_size))

    ! Tendency logicals
    l_tend_pr_tot_w   = current_state%w%active
    l_tend_3d_w   = current_state%w%active .or. l_tend_pr_tot_w

    ! Allocate 3d tendency fields upon availability
    if (l_tend_3d_w) then
      allocate( tend_3d_w(current_state%local_grid%size(Z_INDEX),  &
                          current_state%local_grid%size(Y_INDEX),  &
                          current_state%local_grid%size(X_INDEX)   )    )
    endif

    ! Allocate profile tendency fields upon availability
    if (l_tend_pr_tot_w) then
      allocate( tend_pr_tot_w(current_state%local_grid%size(Z_INDEX)) )
    endif

  end subroutine initialisation_callback  


  subroutine finalisation_callback(current_state)
    type(model_state_type), target, intent(inout) :: current_state

    if (allocated(w_buoyancy)) deallocate(w_buoyancy)

    if (allocated(tend_3d_w)) deallocate(tend_3d_w)

    if (allocated(tend_pr_tot_w)) deallocate(tend_pr_tot_w)

  end subroutine finalisation_callback


  !> Called for each column per timestep this will calculate the buoyancy terms for the SW field
  !! @param current_state The current model state
  !! @param target_(x/y)_index This is the index with the halos subtracted. This is needed so that diagnostic does
  !!                           not include halos and to prevent array out-of-bounds
  subroutine timestep_callback(current_state)
    type(model_state_type), target, intent(inout) :: current_state

    integer :: k, n
    integer :: current_x_index, current_y_index, target_x_index, target_y_index
    logical :: calculate_diagnostics

    calculate_diagnostics = current_state%diagnostic_sample_timestep &
                            .and. .not. current_state%halo_column

    current_x_index=current_state%column_local_x
    current_y_index=current_state%column_local_y
    target_y_index=current_y_index-current_state%local_grid%halo_size(Y_INDEX)
    target_x_index=current_x_index-current_state%local_grid%halo_size(X_INDEX)

    ! Zero profile tendency totals on first instance in the sum
    if (current_state%first_timestep_column) then
      if (l_tend_pr_tot_w) then
        tend_pr_tot_w(:)=0.0_DEFAULT_PRECISION
      endif
    endif  ! zero totals

    if (calculate_diagnostics) &
        call save_precomponent_tendencies(current_state, current_x_index, current_y_index, target_x_index, target_y_index)
    
#ifdef W_ACTIVE
    if (.not. current_state%passive_th .and. current_state%th%active) then
      do k=2,current_state%local_grid%size(Z_INDEX)-1    
        w_buoyancy(k)=(0.5_DEFAULT_PRECISION*current_state%global_grid%configuration%vertical%buoy_co(k))*&
             (current_state%th%data(k, current_state%column_local_y, current_state%column_local_x)&
             +current_state%th%data(k+1, current_state%column_local_y, current_state%column_local_x))
        current_state%sw%data(k, current_state%column_local_y, current_state%column_local_x)=&
             current_state%sw%data(k, current_state%column_local_y, current_state%column_local_x)+w_buoyancy(k)             
      end do
    end if
    if (.not. current_state%passive_q .and. current_state%number_q_fields .gt. 0) then
      if (current_state%use_anelastic_equations) then                                                      
        do n=1,current_state%number_q_fields
          do k=2,current_state%local_grid%size(Z_INDEX)-1            
            current_state%sw%data(k, current_state%column_local_y, current_state%column_local_x)=&
                 current_state%sw%data(k, current_state%column_local_y, current_state%column_local_x)+&
                 (0.5_DEFAULT_PRECISION*current_state%global_grid%configuration%vertical%buoy_co(k))*&
                 current_state%cq(n)* (current_state%global_grid%configuration%vertical%thref(k)*&
                 current_state%q(n)%data(k, current_state%column_local_y, current_state%column_local_x)+&
                 current_state%global_grid%configuration%vertical%thref(k+1)*&
                 current_state%q(n)%data(k+1, current_state%column_local_y, current_state%column_local_x))
          end do
        end do
      else                                                                     
        do n=1,current_state%number_q_fields
          do k=2,current_state%local_grid%size(Z_INDEX)-1
             current_state%sw%data(k, current_state%column_local_y, current_state%column_local_x)=&
                  current_state%sw%data(k, current_state%column_local_y, current_state%column_local_x)+&
                  G_over_2*current_state%cq(n)*&
                  (current_state%q(n)%data(k, current_state%column_local_y, current_state%column_local_x)+&
                  current_state%q(n)%data(k+1, current_state%column_local_y, current_state%column_local_x))
          end do
        end do
      end if
    end if
#endif

    if (calculate_diagnostics) &
        call compute_component_tendencies(current_state, current_x_index, current_y_index, target_x_index, target_y_index)

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
    if (l_tend_3d_w) then
      tend_3d_w(:,tyn,txn)=current_state%sw%data(:,cyn,cxn)
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
    if (l_tend_3d_w) then
      tend_3d_w(:,tyn,txn)=current_state%sw%data(:,cyn,cxn) - tend_3d_w(:,tyn,txn)
    endif

   ! Add local tendency fields to the profile total
    if (l_tend_pr_tot_w) then
      tend_pr_tot_w(:)=tend_pr_tot_w(:) + tend_3d_w(:,tyn,txn)
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

end module buoyancy_mod
