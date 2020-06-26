!> Damping applied to the W field at some point to stop stuff flying up and off
module damping_mod
  use monc_component_mod, only : COMPONENT_ARRAY_FIELD_TYPE, &
      COMPONENT_DOUBLE_DATA_TYPE, component_descriptor_type, &
      component_field_value_type, component_field_information_type
  use grids_mod, only : Z_INDEX, Y_INDEX, X_INDEX
  use state_mod, only : model_state_type
  use collections_mod, only : map_type
  use optionsdatabase_mod, only : options_get_real, options_get_integer
  use registry_mod, only : is_component_enabled
  use logging_mod, only : LOG_ERROR, log_master_log
  use datadefn_mod, only : DEFAULT_PRECISION
  use q_indices_mod, only: get_q_index, standard_q_names

  implicit none

#ifndef TEST_MODE
  private
#endif

  real(kind=DEFAULT_PRECISION) :: dmptim,& !< Layer timescale
       zdmp,& !< The point (m) where the damping starts
       hdmp   !< The height (m) of the damping layer

  ! Local tendency diagnostic variables for this component
  ! 3D tendency fields and logicals for their use
  real(kind=DEFAULT_PRECISION), dimension(:,:,:), allocatable ::     &
       tend_3d_u, tend_3d_v, tend_3d_w, tend_3d_th,tend_3d_qv,       &
       tend_3d_ql,tend_3d_qi,tend_3d_qr,tend_3d_qs,tend_3d_qg,       &
       tend_3d_tabs
  logical :: l_tend_3d_u, l_tend_3d_v, l_tend_3d_w, l_tend_3d_th,l_tend_3d_qv,       &
             l_tend_3d_ql,l_tend_3d_qi,l_tend_3d_qr,l_tend_3d_qs,l_tend_3d_qg,       &
             l_tend_3d_tabs
  ! Local mean tendency profile fields and logicals for their use
  real(kind=DEFAULT_PRECISION), dimension(:), allocatable ::         &
       tend_pr_tot_u, tend_pr_tot_v, tend_pr_tot_w, tend_pr_tot_th,tend_pr_tot_qv,       &
       tend_pr_tot_ql,tend_pr_tot_qi,tend_pr_tot_qr,tend_pr_tot_qs,tend_pr_tot_qg,       &
       tend_pr_tot_tabs
  logical :: l_tend_pr_tot_u, l_tend_pr_tot_v, l_tend_pr_tot_w, l_tend_pr_tot_th,l_tend_pr_tot_qv,       &
             l_tend_pr_tot_ql,l_tend_pr_tot_qi,l_tend_pr_tot_qr,l_tend_pr_tot_qs,l_tend_pr_tot_qg,       &
             l_tend_pr_tot_tabs
  ! q indices
  integer :: iqv=0, iql=0, iqr=0, iqi=0, iqs=0, iqg=0

  ! tke tendency diagnostic
  real(kind=DEFAULT_PRECISION), dimension(:), allocatable :: tend_pr_tot_tke
  logical :: l_tend_pr_tot_tke

  public damping_get_descriptor

contains

  !> Descriptor of this component for registration
  !! @returns The damping component descriptor
  type(component_descriptor_type) function damping_get_descriptor()
    damping_get_descriptor%name="damping"
    damping_get_descriptor%version=0.1
    damping_get_descriptor%initialisation=>init_callback
    damping_get_descriptor%timestep=>timestep_callback
    damping_get_descriptor%finalisation=>finalisation_callback

    damping_get_descriptor%field_value_retrieval=>field_value_retrieval_callback
    damping_get_descriptor%field_information_retrieval=>field_information_retrieval_callback
    allocate(damping_get_descriptor%published_fields(11+11+1))

    damping_get_descriptor%published_fields(1)= "tend_u_damping_3d_local"
    damping_get_descriptor%published_fields(2)= "tend_v_damping_3d_local"
    damping_get_descriptor%published_fields(3)= "tend_w_damping_3d_local"
    damping_get_descriptor%published_fields(4)= "tend_th_damping_3d_local"
    damping_get_descriptor%published_fields(5)= "tend_qv_damping_3d_local"
    damping_get_descriptor%published_fields(6)= "tend_ql_damping_3d_local"
    damping_get_descriptor%published_fields(7)= "tend_qi_damping_3d_local"
    damping_get_descriptor%published_fields(8)= "tend_qr_damping_3d_local"
    damping_get_descriptor%published_fields(9)= "tend_qs_damping_3d_local"
    damping_get_descriptor%published_fields(10)="tend_qg_damping_3d_local"
    damping_get_descriptor%published_fields(11)="tend_tabs_damping_3d_local"

    damping_get_descriptor%published_fields(11+1)= "tend_u_damping_profile_total_local"
    damping_get_descriptor%published_fields(11+2)= "tend_v_damping_profile_total_local"
    damping_get_descriptor%published_fields(11+3)= "tend_w_damping_profile_total_local"
    damping_get_descriptor%published_fields(11+4)= "tend_th_damping_profile_total_local"
    damping_get_descriptor%published_fields(11+5)= "tend_qv_damping_profile_total_local"
    damping_get_descriptor%published_fields(11+6)= "tend_ql_damping_profile_total_local"
    damping_get_descriptor%published_fields(11+7)= "tend_qi_damping_profile_total_local"
    damping_get_descriptor%published_fields(11+8)= "tend_qr_damping_profile_total_local"
    damping_get_descriptor%published_fields(11+9)= "tend_qs_damping_profile_total_local"
    damping_get_descriptor%published_fields(11+10)="tend_qg_damping_profile_total_local"
    damping_get_descriptor%published_fields(11+11)="tend_tabs_damping_profile_total_local"

    damping_get_descriptor%published_fields(11+11+1)="tend_tke_damping_profile_total_local"


  end function damping_get_descriptor

  !> On initialisation will set up data structures and field values
  !! @param current_state The current model state_mod
  subroutine init_callback(current_state)
    type(model_state_type), target, intent(inout) :: current_state

    integer :: k
    logical :: l_qdiag

    if (.not. is_component_enabled(current_state%options_database, "mean_profiles")) then
      call log_master_log(LOG_ERROR, "Damping requires the mean profiles component to be enabled")
    end if    

    dmptim=options_get_real(current_state%options_database, "dmptim")
    zdmp=options_get_real(current_state%options_database, "zdmp")
    hdmp=options_get_real(current_state%options_database, "hdmp")

    allocate(current_state%global_grid%configuration%vertical%dmpco(current_state%local_grid%size(Z_INDEX)), &
         current_state%global_grid%configuration%vertical%dmpcoz(current_state%local_grid%size(Z_INDEX)))
    current_state%global_grid%configuration%vertical%dmpco(:)=0.
    current_state%global_grid%configuration%vertical%dmpcoz(:)=0.
    do k=current_state%local_grid%size(Z_INDEX),1,-1
      current_state%global_grid%configuration%vertical%kdmpmin=k
      if (current_state%global_grid%configuration%vertical%zn(k) .ge. zdmp) then
        current_state%global_grid%configuration%vertical%dmpco(k)=dmptim*(exp((&
             current_state%global_grid%configuration%vertical%zn(k)-zdmp)/hdmp)-1.0)
      end if
      if (current_state%global_grid%configuration%vertical%z(k) .ge. zdmp) then
        current_state%global_grid%configuration%vertical%dmpcoz(K)=dmptim*(exp((&
             current_state%global_grid%configuration%vertical%z(K)-zdmp)/hdmp)-1.0)
      end if
      if(current_state%global_grid%configuration%vertical%zn(k).lt. zdmp) exit
    end do

    ! Set tendency diagnostic logicals based on availability
    ! Need to use 3d tendencies to compute the profiles, so they will be allocated
    !      in the case where profiles are available
    l_qdiag =  (.not. current_state%passive_q .and. current_state%number_q_fields .gt. 0)

    l_tend_pr_tot_u   = current_state%u%active
    l_tend_pr_tot_v   = current_state%v%active
    l_tend_pr_tot_w   = current_state%w%active
    l_tend_pr_tot_th  = current_state%th%active
    l_tend_pr_tot_qv  = l_qdiag .and. current_state%number_q_fields .ge. 1
    l_tend_pr_tot_ql  = l_qdiag .and. current_state%number_q_fields .ge. 2
    l_tend_pr_tot_qi  = l_qdiag .and. current_state%number_q_fields .ge. 11
    l_tend_pr_tot_qr  = l_qdiag .and. current_state%number_q_fields .ge. 11
    l_tend_pr_tot_qs  = l_qdiag .and. current_state%number_q_fields .ge. 11
    l_tend_pr_tot_qg  = l_qdiag .and. current_state%number_q_fields .ge. 11
    l_tend_pr_tot_tabs  = l_tend_pr_tot_th

    l_tend_3d_u   = current_state%u%active .or. l_tend_pr_tot_u
    l_tend_3d_v   = current_state%v%active .or. l_tend_pr_tot_v
    l_tend_3d_w   = current_state%w%active .or. l_tend_pr_tot_w
    l_tend_3d_th  = current_state%th%active .or. l_tend_pr_tot_th
    l_tend_3d_qv  = (l_qdiag .and. current_state%number_q_fields .ge. 1) .or. l_tend_pr_tot_qv
    l_tend_3d_ql  = (l_qdiag .and. current_state%number_q_fields .ge. 2) .or. l_tend_pr_tot_ql
    l_tend_3d_qi  = (l_qdiag .and. current_state%number_q_fields .ge. 11) .or. l_tend_pr_tot_qi
    l_tend_3d_qr  = (l_qdiag .and. current_state%number_q_fields .ge. 11) .or. l_tend_pr_tot_qr
    l_tend_3d_qs  = (l_qdiag .and. current_state%number_q_fields .ge. 11) .or. l_tend_pr_tot_qs
    l_tend_3d_qg  = (l_qdiag .and. current_state%number_q_fields .ge. 11) .or. l_tend_pr_tot_qg
    l_tend_3d_tabs = l_tend_3d_th

    l_tend_pr_tot_tke = current_state%u%active .and. current_state%v%active .and. current_state%w%active
 
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
    if (l_tend_3d_w) then
      allocate( tend_3d_w(current_state%local_grid%size(Z_INDEX),  &
                          current_state%local_grid%size(Y_INDEX),  &
                          current_state%local_grid%size(X_INDEX)   )    )
    endif
    if (l_tend_3d_th) then
      allocate( tend_3d_th(current_state%local_grid%size(Z_INDEX),  &
                           current_state%local_grid%size(Y_INDEX),  &
                           current_state%local_grid%size(X_INDEX)   )    )
    endif
    if (l_tend_3d_qv) then
      iqv=get_q_index(standard_q_names%VAPOUR, 'damping')
      allocate( tend_3d_qv(current_state%local_grid%size(Z_INDEX),  &
                           current_state%local_grid%size(Y_INDEX),  &
                           current_state%local_grid%size(X_INDEX)   )    )
    endif
    if (l_tend_3d_ql) then
      iql=get_q_index(standard_q_names%CLOUD_LIQUID_MASS, 'damping')
      allocate( tend_3d_ql(current_state%local_grid%size(Z_INDEX),  &
                           current_state%local_grid%size(Y_INDEX),  &
                           current_state%local_grid%size(X_INDEX)   )    )
    endif
    if (l_tend_3d_qi) then
      iqi=get_q_index(standard_q_names%ICE_MASS, 'damping')
      allocate( tend_3d_qi(current_state%local_grid%size(Z_INDEX),  &
                           current_state%local_grid%size(Y_INDEX),  &
                           current_state%local_grid%size(X_INDEX)   )    )
    endif
    if (l_tend_3d_qr) then
      iqr=get_q_index(standard_q_names%RAIN_MASS, 'damping')
      allocate( tend_3d_qr(current_state%local_grid%size(Z_INDEX),  &
                           current_state%local_grid%size(Y_INDEX),  &
                           current_state%local_grid%size(X_INDEX)   )    )
    endif
    if (l_tend_3d_qs) then
      iqs=get_q_index(standard_q_names%SNOW_MASS, 'damping')
      allocate( tend_3d_qs(current_state%local_grid%size(Z_INDEX),  &
                           current_state%local_grid%size(Y_INDEX),  &
                           current_state%local_grid%size(X_INDEX)   )    )
    endif
    if (l_tend_3d_qg) then
      iqg=get_q_index(standard_q_names%GRAUPEL_MASS, 'damping')
      allocate( tend_3d_qg(current_state%local_grid%size(Z_INDEX),  &
                           current_state%local_grid%size(Y_INDEX),  &
                           current_state%local_grid%size(X_INDEX)   )    )
    endif
    if (l_tend_3d_tabs) then
      allocate( tend_3d_tabs(current_state%local_grid%size(Z_INDEX),  &
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
    if (l_tend_pr_tot_w) then
      allocate( tend_pr_tot_w(current_state%local_grid%size(Z_INDEX)) )
    endif
    if (l_tend_pr_tot_th) then
      allocate( tend_pr_tot_th(current_state%local_grid%size(Z_INDEX)) )
    endif
    if (l_tend_pr_tot_qv) then
      allocate( tend_pr_tot_qv(current_state%local_grid%size(Z_INDEX)) )
    endif
    if (l_tend_pr_tot_ql) then
      allocate( tend_pr_tot_ql(current_state%local_grid%size(Z_INDEX)) )
    endif
    if (l_tend_pr_tot_qi) then
      allocate( tend_pr_tot_qi(current_state%local_grid%size(Z_INDEX)) )
    endif
    if (l_tend_pr_tot_qr) then
      allocate( tend_pr_tot_qr(current_state%local_grid%size(Z_INDEX)) )
    endif
    if (l_tend_pr_tot_qs) then
      allocate( tend_pr_tot_qs(current_state%local_grid%size(Z_INDEX)) )
    endif
    if (l_tend_pr_tot_qg) then
      allocate( tend_pr_tot_qg(current_state%local_grid%size(Z_INDEX)) )
    endif
    if (l_tend_pr_tot_tabs) then
      allocate( tend_pr_tot_tabs(current_state%local_grid%size(Z_INDEX)) )
    endif

    ! Allocate profile tendency tke upon availability
    if (l_tend_pr_tot_tke) then
      allocate( tend_pr_tot_tke(current_state%local_grid%size(Z_INDEX)) )
    endif

  end subroutine init_callback


  subroutine finalisation_callback(current_state)
    type(model_state_type), target, intent(inout) :: current_state

    if (allocated(tend_3d_u)) deallocate(tend_3d_u)
    if (allocated(tend_3d_v)) deallocate(tend_3d_v)
    if (allocated(tend_3d_w)) deallocate(tend_3d_w)
    if (allocated(tend_3d_th)) deallocate(tend_3d_th)
    if (allocated(tend_3d_qv)) deallocate(tend_3d_qv)
    if (allocated(tend_3d_ql)) deallocate(tend_3d_ql)
    if (allocated(tend_3d_qi)) deallocate(tend_3d_qi)
    if (allocated(tend_3d_qr)) deallocate(tend_3d_qr)
    if (allocated(tend_3d_qs)) deallocate(tend_3d_qs)
    if (allocated(tend_3d_qg)) deallocate(tend_3d_qg)
    if (allocated(tend_3d_tabs)) deallocate(tend_3d_tabs)

    if (allocated(tend_pr_tot_u)) deallocate(tend_pr_tot_u)
    if (allocated(tend_pr_tot_v)) deallocate(tend_pr_tot_v)
    if (allocated(tend_pr_tot_w)) deallocate(tend_pr_tot_w)
    if (allocated(tend_pr_tot_th)) deallocate(tend_pr_tot_th)
    if (allocated(tend_pr_tot_qv)) deallocate(tend_pr_tot_qv)
    if (allocated(tend_pr_tot_ql)) deallocate(tend_pr_tot_ql)
    if (allocated(tend_pr_tot_qi)) deallocate(tend_pr_tot_qi)
    if (allocated(tend_pr_tot_qr)) deallocate(tend_pr_tot_qr)
    if (allocated(tend_pr_tot_qs)) deallocate(tend_pr_tot_qs)
    if (allocated(tend_pr_tot_qg)) deallocate(tend_pr_tot_qg)
    if (allocated(tend_pr_tot_tabs)) deallocate(tend_pr_tot_tabs)

    if (allocated(tend_pr_tot_tke)) deallocate(tend_pr_tot_tke)

  end subroutine finalisation_callback


  !> For each data column will calculate the damping term and apply this to the source term for that field
  !! @param current_state The current model state_mod
  !! @param target_(x/y)_index This is the index with the halos subtracted. This is needed so that diagnostic does
  !!                           not include halos and to prevent array out-of-bounds
  subroutine timestep_callback(current_state)
    type(model_state_type), target, intent(inout) :: current_state

    integer :: k, i
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
      if (l_tend_pr_tot_u) then
        tend_pr_tot_u(:)= 0.0_DEFAULT_PRECISION
      endif
      if (l_tend_pr_tot_v) then
        tend_pr_tot_v(:)= 0.0_DEFAULT_PRECISION
      endif
      if (l_tend_pr_tot_w) then
        tend_pr_tot_w(:)= 0.0_DEFAULT_PRECISION
      endif
      if (l_tend_pr_tot_th) then
        tend_pr_tot_th(:)=0.0_DEFAULT_PRECISION
      endif
      if (l_tend_pr_tot_qv) then
        tend_pr_tot_qv(:)=0.0_DEFAULT_PRECISION
      endif
      if (l_tend_pr_tot_ql) then
        tend_pr_tot_ql(:)=0.0_DEFAULT_PRECISION
      endif
      if (l_tend_pr_tot_qi) then
        tend_pr_tot_qi(:)=0.0_DEFAULT_PRECISION
      endif
      if (l_tend_pr_tot_qr) then
        tend_pr_tot_qr(:)=0.0_DEFAULT_PRECISION
      endif
      if (l_tend_pr_tot_qs) then
        tend_pr_tot_qs(:)=0.0_DEFAULT_PRECISION
      endif
      if (l_tend_pr_tot_qg) then
        tend_pr_tot_qg(:)=0.0_DEFAULT_PRECISION
      endif
      if (l_tend_pr_tot_tabs) then
        tend_pr_tot_tabs(:)=0.0_DEFAULT_PRECISION
      endif

      if (l_tend_pr_tot_tke) then
        tend_pr_tot_tke(:)=0.0_DEFAULT_PRECISION
      endif
    endif  ! zero totals


    if (current_state%halo_column .and. current_state%timestep <3) return

    if (calculate_diagnostics) &
        call save_precomponent_tendencies(current_state, current_x_index, current_y_index, target_x_index, target_y_index)
    
    do k=current_state%global_grid%configuration%vertical%kdmpmin,current_state%local_grid%size(Z_INDEX)
#ifdef U_ACTIVE
      current_state%su%data(k, current_state%column_local_y, current_state%column_local_x)=current_state%su%data(k, &
           current_state%column_local_y, current_state%column_local_x)-&
           current_state%global_grid%configuration%vertical%dmpco(k)*(current_state%zu%data(k, current_state%column_local_y, &
           current_state%column_local_x)- (current_state%global_grid%configuration%vertical%olzubar(k)-current_state%ugal))
#endif
#ifdef V_ACTIVE
      current_state%sv%data(k, current_state%column_local_y, current_state%column_local_x)=current_state%sv%data(k, &
           current_state%column_local_y, current_state%column_local_x)-&
           current_state%global_grid%configuration%vertical%dmpco(k)*(current_state%zv%data(k, current_state%column_local_y, &
           current_state%column_local_x)-(current_state%global_grid%configuration%vertical%olzvbar(k)-current_state%vgal))
#endif
      if (current_state%th%active) then
        current_state%sth%data(k, current_state%column_local_y, current_state%column_local_x)=current_state%sth%data(k, &
             current_state%column_local_y, current_state%column_local_x)-&
             current_state%global_grid%configuration%vertical%dmpco(k)*(current_state%zth%data(k, current_state%column_local_y, &
             current_state%column_local_x)-current_state%global_grid%configuration%vertical%olzthbar(k))
      end if
      
      do i=1,current_state%number_q_fields
        if (current_state%q(i)%active) then
          current_state%sq(i)%data(k, current_state%column_local_y, current_state%column_local_x)=current_state%sq(i)%data(k, &
               current_state%column_local_y, current_state%column_local_x)-&
               current_state%global_grid%configuration%vertical%dmpco(k)*&
               (current_state%zq(i)%data(k, current_state%column_local_y, current_state%column_local_x)-&
               current_state%global_grid%configuration%vertical%olzqbar(k,i))
        end if
      end do
    end do
#ifdef W_ACTIVE
    do k=current_state%global_grid%configuration%vertical%kdmpmin,current_state%local_grid%size(Z_INDEX)-1    
      current_state%sw%data(k, current_state%column_local_y, current_state%column_local_x)=current_state%sw%data(k, &
           current_state%column_local_y, current_state%column_local_x)-&
           current_state%global_grid%configuration%vertical%dmpcoz(k)*&
           current_state%zw%data(k, current_state%column_local_y, current_state%column_local_x)
    end do
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
    if (l_tend_3d_u) then
      tend_3d_u(:,tyn,txn)=current_state%su%data(:,cyn,cxn)
    endif
    if (l_tend_3d_v) then
      tend_3d_v(:,tyn,txn)=current_state%sv%data(:,cyn,cxn)
    endif
    if (l_tend_3d_w) then
      tend_3d_w(:,tyn,txn)=current_state%sw%data(:,cyn,cxn)
    endif
    if (l_tend_3d_th) then
      tend_3d_th(:,tyn,txn)=current_state%sth%data(:,cyn,cxn)
    endif
    if (l_tend_3d_qv) then
      tend_3d_qv(:,tyn,txn)=current_state%sq(iqv)%data(:,cyn,cxn)
    endif
    if (l_tend_3d_ql) then
      tend_3d_ql(:,tyn,txn)=current_state%sq(iql)%data(:,cyn,cxn)
    endif
    if (l_tend_3d_qi) then
      tend_3d_qi(:,tyn,txn)=current_state%sq(iqi)%data(:,cyn,cxn)
    endif
    if (l_tend_3d_qr) then
      tend_3d_qr(:,tyn,txn)=current_state%sq(iqr)%data(:,cyn,cxn)
    endif
    if (l_tend_3d_qs) then
      tend_3d_qs(:,tyn,txn)=current_state%sq(iqs)%data(:,cyn,cxn)
    endif
    if (l_tend_3d_qg) then
      tend_3d_qg(:,tyn,txn)=current_state%sq(iqg)%data(:,cyn,cxn)
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

    real(kind=DEFAULT_PRECISION), dimension(current_state%local_grid%size(Z_INDEX)) :: &
                            uu_tendency, vv_tendency, ww_tendency
    integer :: k

    ! Calculate change in tendency due to component
    if (l_tend_3d_u) then
      tend_3d_u(:,tyn,txn)=current_state%su%data(:,cyn,cxn)       - tend_3d_u(:,tyn,txn)
    endif
    if (l_tend_3d_v) then
      tend_3d_v(:,tyn,txn)=current_state%sv%data(:,cyn,cxn)       - tend_3d_v(:,tyn,txn)
    endif
    if (l_tend_3d_w) then
      tend_3d_w(:,tyn,txn)=current_state%sw%data(:,cyn,cxn)       - tend_3d_w(:,tyn,txn)
    endif
    if (l_tend_3d_th) then
      tend_3d_th(:,tyn,txn)=current_state%sth%data(:,cyn,cxn)     - tend_3d_th(:,tyn,txn)
    endif
    if (l_tend_3d_qv) then
      tend_3d_qv(:,tyn,txn)=current_state%sq(iqv)%data(:,cyn,cxn) - tend_3d_qv(:,tyn,txn)
    endif
    if (l_tend_3d_ql) then
      tend_3d_ql(:,tyn,txn)=current_state%sq(iql)%data(:,cyn,cxn) - tend_3d_ql(:,tyn,txn)
    endif
    if (l_tend_3d_qi) then
      tend_3d_qi(:,tyn,txn)=current_state%sq(iqi)%data(:,cyn,cxn) - tend_3d_qi(:,tyn,txn)
    endif
    if (l_tend_3d_qr) then
      tend_3d_qr(:,tyn,txn)=current_state%sq(iqr)%data(:,cyn,cxn) - tend_3d_qr(:,tyn,txn)
    endif
    if (l_tend_3d_qs) then
      tend_3d_qs(:,tyn,txn)=current_state%sq(iqs)%data(:,cyn,cxn) - tend_3d_qs(:,tyn,txn)
    endif
    if (l_tend_3d_qg) then
      tend_3d_qg(:,tyn,txn)=current_state%sq(iqg)%data(:,cyn,cxn) - tend_3d_qg(:,tyn,txn)
    endif
    if (l_tend_3d_tabs) then
      tend_3d_tabs(:,tyn,txn)=   &
       current_state%sth%data(:,cyn,cxn) * current_state%global_grid%configuration%vertical%rprefrcp(:)   &
        - tend_3d_tabs(:,tyn,txn)
    endif

   ! Add local tendency fields to the profile total
    if (l_tend_pr_tot_u) then
      tend_pr_tot_u(:)=tend_pr_tot_u(:)   + tend_3d_u(:,tyn,txn)
    endif
    if (l_tend_pr_tot_v) then
      tend_pr_tot_v(:)=tend_pr_tot_v(:)   + tend_3d_v(:,tyn,txn)
    endif
    if (l_tend_pr_tot_w) then
      tend_pr_tot_w(:)=tend_pr_tot_w(:)   + tend_3d_w(:,tyn,txn)
    endif
    if (l_tend_pr_tot_th) then
      tend_pr_tot_th(:)=tend_pr_tot_th(:) + tend_3d_th(:,tyn,txn)
    endif
    if (l_tend_pr_tot_qv) then
      tend_pr_tot_qv(:)=tend_pr_tot_qv(:) + tend_3d_qv(:,tyn,txn)
    endif
    if (l_tend_pr_tot_ql) then
      tend_pr_tot_ql(:)=tend_pr_tot_ql(:) + tend_3d_ql(:,tyn,txn)
    endif
    if (l_tend_pr_tot_qi) then
      tend_pr_tot_qi(:)=tend_pr_tot_qi(:) + tend_3d_qi(:,tyn,txn)
    endif
    if (l_tend_pr_tot_qr) then
      tend_pr_tot_qr(:)=tend_pr_tot_qr(:) + tend_3d_qr(:,tyn,txn)
    endif
    if (l_tend_pr_tot_qs) then
      tend_pr_tot_qs(:)=tend_pr_tot_qs(:) + tend_3d_qs(:,tyn,txn)
    endif
    if (l_tend_pr_tot_qg) then
      tend_pr_tot_qg(:)=tend_pr_tot_qg(:) + tend_3d_qg(:,tyn,txn)
    endif
    if (l_tend_pr_tot_tabs) then
      tend_pr_tot_tabs(:)=tend_pr_tot_tabs(:) + tend_3d_tabs(:,tyn,txn)
    endif

! Estimate contribution to TKE budget due to damping.
! Currently uses mean state at beginning of timestep where mean state at end is really required 
! (hence 'provisional' comments); assumption is mean state is very slowly varying. 

    if (l_tend_pr_tot_tke) then
      do k=2, current_state%local_grid%size(Z_INDEX)
  
        uu_tendency(k) = ((current_state%zu%data(k,cyn,cxn)                                          &
                              + tend_3d_u(k,tyn,txn)*current_state%dtm * 2.0_DEFAULT_PRECISION       &
                           - current_state%global_grid%configuration%vertical%olzubar(k) )**2        & ! provisional
                          -                                                                          &
                          (current_state%zu%data(k,cyn,cxn)                                          &
                           - current_state%global_grid%configuration%vertical%olzubar(k) )**2        & ! n-1
                         ) /  (current_state%dtm * 2.0_DEFAULT_PRECISION)
                        
        vv_tendency(k) = ((current_state%zv%data(k,cyn,cxn)                                          &
                              + tend_3d_v(k,tyn,txn)*current_state%dtm * 2.0_DEFAULT_PRECISION       &
                           - current_state%global_grid%configuration%vertical%olzvbar(k) )**2        & ! provisional
                          -                                                                          &
                          (current_state%zv%data(k,cyn,cxn)                                          &
                           - current_state%global_grid%configuration%vertical%olzvbar(k) )**2        & ! n-1
                         ) /  (current_state%dtm * 2.0_DEFAULT_PRECISION)

        ww_tendency(k) = ((current_state%zw%data(k,cyn,cxn)                                          & ! w_bar assumed zero
                              + tend_3d_w(k,tyn,txn)*current_state%dtm * 2.0_DEFAULT_PRECISION )**2  & ! provisional
                          -                                                                          &
                          current_state%zw%data(k,cyn,cxn)**2                                        & ! n-1
                         ) /  (current_state%dtm * 2.0_DEFAULT_PRECISION)


      enddo

      ! handle surface so that interpolation goes to zero
      uu_tendency(1) = -uu_tendency(2)
      vv_tendency(1) = -vv_tendency(2)

      ! interpolate to w-levels    
      do k=2, current_state%local_grid%size(Z_INDEX)-1
        tend_pr_tot_tke(k)=tend_pr_tot_tke(k) + 0.5_DEFAULT_PRECISION * (& 
          0.5_DEFAULT_PRECISION * (uu_tendency(k)+uu_tendency(k+1)) + &
          0.5_DEFAULT_PRECISION * (vv_tendency(k)+vv_tendency(k+1)) + &
          ww_tendency(k) )
      enddo

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
    strcomp=INDEX(name, "_damping_3d_local")
    if (strcomp .ne. 0) then
      field_information%field_type=COMPONENT_ARRAY_FIELD_TYPE
      field_information%number_dimensions=3
      field_information%dimension_sizes(1)=current_state%local_grid%size(Z_INDEX)
      field_information%dimension_sizes(2)=current_state%local_grid%size(Y_INDEX)
      field_information%dimension_sizes(3)=current_state%local_grid%size(X_INDEX)
      field_information%data_type=COMPONENT_DOUBLE_DATA_TYPE

      if      (name .eq. "tend_u_damping_3d_local") then
        field_information%enabled=l_tend_3d_u
      else if (name .eq. "tend_v_damping_3d_local") then
        field_information%enabled=l_tend_3d_v
      else if (name .eq. "tend_w_damping_3d_local") then
        field_information%enabled=l_tend_3d_w
      else if (name .eq. "tend_th_damping_3d_local") then
        field_information%enabled=l_tend_3d_th
      else if (name .eq. "tend_qv_damping_3d_local") then
        field_information%enabled=l_tend_3d_qv
      else if (name .eq. "tend_ql_damping_3d_local") then
        field_information%enabled=l_tend_3d_ql
      else if (name .eq. "tend_qi_damping_3d_local") then
        field_information%enabled=l_tend_3d_qi
      else if (name .eq. "tend_qr_damping_3d_local") then
        field_information%enabled=l_tend_3d_qr
      else if (name .eq. "tend_qs_damping_3d_local") then
        field_information%enabled=l_tend_3d_qs
      else if (name .eq. "tend_qg_damping_3d_local") then
        field_information%enabled=l_tend_3d_qg
      else if (name .eq. "tend_tabs_damping_3d_local") then
        field_information%enabled=l_tend_3d_tabs
      else
        field_information%enabled=.true.
      end if

    end if !end 3d check

    ! Field information for profiles
    strcomp=INDEX(name, "_damping_profile_total_local")
    if (strcomp .ne. 0) then
      field_information%field_type=COMPONENT_ARRAY_FIELD_TYPE
      field_information%number_dimensions=1
      field_information%dimension_sizes(1)=current_state%local_grid%size(Z_INDEX)
      field_information%data_type=COMPONENT_DOUBLE_DATA_TYPE

      if      (name .eq. "tend_u_damping_profile_total_local") then
        field_information%enabled=l_tend_pr_tot_u
      else if (name .eq. "tend_v_damping_profile_total_local") then
        field_information%enabled=l_tend_pr_tot_v
      else if (name .eq. "tend_w_damping_profile_total_local") then
        field_information%enabled=l_tend_pr_tot_w
      else if (name .eq. "tend_th_damping_profile_total_local") then
        field_information%enabled=l_tend_pr_tot_th
      else if (name .eq. "tend_qv_damping_profile_total_local") then
        field_information%enabled=l_tend_pr_tot_qv
      else if (name .eq. "tend_ql_damping_profile_total_local") then
        field_information%enabled=l_tend_pr_tot_ql
      else if (name .eq. "tend_qi_damping_profile_total_local") then
        field_information%enabled=l_tend_pr_tot_qi
      else if (name .eq. "tend_qr_damping_profile_total_local") then
        field_information%enabled=l_tend_pr_tot_qr
      else if (name .eq. "tend_qs_damping_profile_total_local") then
        field_information%enabled=l_tend_pr_tot_qs
      else if (name .eq. "tend_qg_damping_profile_total_local") then
        field_information%enabled=l_tend_pr_tot_qg
      else if (name .eq. "tend_tabs_damping_profile_total_local") then
        field_information%enabled=l_tend_pr_tot_tabs

      else if (name .eq. "tend_tke_damping_profile_total_local") then
        field_information%enabled=l_tend_pr_tot_tke

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
    if (name .eq. "tend_u_damping_3d_local" .and. allocated(tend_3d_u)) then
      call set_published_field_value(field_value, real_3d_field=tend_3d_u)
    else if (name .eq. "tend_v_damping_3d_local" .and. allocated(tend_3d_v)) then
      call set_published_field_value(field_value, real_3d_field=tend_3d_v)
    else if (name .eq. "tend_w_damping_3d_local" .and. allocated(tend_3d_w)) then
      call set_published_field_value(field_value, real_3d_field=tend_3d_w)
    else if (name .eq. "tend_th_damping_3d_local" .and. allocated(tend_3d_th)) then
      call set_published_field_value(field_value, real_3d_field=tend_3d_th)
    else if (name .eq. "tend_qv_damping_3d_local" .and. allocated(tend_3d_qv)) then
      call set_published_field_value(field_value, real_3d_field=tend_3d_qv)
    else if (name .eq. "tend_ql_damping_3d_local" .and. allocated(tend_3d_ql)) then
      call set_published_field_value(field_value, real_3d_field=tend_3d_ql)
    else if (name .eq. "tend_qi_damping_3d_local" .and. allocated(tend_3d_qi)) then
      call set_published_field_value(field_value, real_3d_field=tend_3d_qi)
    else if (name .eq. "tend_qr_damping_3d_local" .and. allocated(tend_3d_qr)) then
      call set_published_field_value(field_value, real_3d_field=tend_3d_qr)
    else if (name .eq. "tend_qs_damping_3d_local" .and. allocated(tend_3d_qs)) then
      call set_published_field_value(field_value, real_3d_field=tend_3d_qs)
    else if (name .eq. "tend_qg_damping_3d_local" .and. allocated(tend_3d_qg)) then
      call set_published_field_value(field_value, real_3d_field=tend_3d_qg)
    else if (name .eq. "tend_tabs_damping_3d_local" .and. allocated(tend_3d_tabs)) then
      call set_published_field_value(field_value, real_3d_field=tend_3d_tabs)

    ! Profile Tendency Fields
    else if (name .eq. "tend_u_damping_profile_total_local" .and. allocated(tend_pr_tot_u)) then
      call set_published_field_value(field_value, real_1d_field=tend_pr_tot_u)
    else if (name .eq. "tend_v_damping_profile_total_local" .and. allocated(tend_pr_tot_v)) then
      call set_published_field_value(field_value, real_1d_field=tend_pr_tot_v)
    else if (name .eq. "tend_w_damping_profile_total_local" .and. allocated(tend_pr_tot_w)) then
      call set_published_field_value(field_value, real_1d_field=tend_pr_tot_w)
    else if (name .eq. "tend_th_damping_profile_total_local" .and. allocated(tend_pr_tot_th)) then
      call set_published_field_value(field_value, real_1d_field=tend_pr_tot_th)
    else if (name .eq. "tend_qv_damping_profile_total_local" .and. allocated(tend_pr_tot_qv)) then
      call set_published_field_value(field_value, real_1d_field=tend_pr_tot_qv)
    else if (name .eq. "tend_ql_damping_profile_total_local" .and. allocated(tend_pr_tot_ql)) then
      call set_published_field_value(field_value, real_1d_field=tend_pr_tot_ql)
    else if (name .eq. "tend_qi_damping_profile_total_local" .and. allocated(tend_pr_tot_qi)) then
      call set_published_field_value(field_value, real_1d_field=tend_pr_tot_qi)
    else if (name .eq. "tend_qr_damping_profile_total_local" .and. allocated(tend_pr_tot_qr)) then
      call set_published_field_value(field_value, real_1d_field=tend_pr_tot_qr)
    else if (name .eq. "tend_qs_damping_profile_total_local" .and. allocated(tend_pr_tot_qs)) then
      call set_published_field_value(field_value, real_1d_field=tend_pr_tot_qs)
    else if (name .eq. "tend_qg_damping_profile_total_local" .and. allocated(tend_pr_tot_qg)) then
      call set_published_field_value(field_value, real_1d_field=tend_pr_tot_qg)
    else if (name .eq. "tend_tabs_damping_profile_total_local" .and. allocated(tend_pr_tot_tabs)) then
      call set_published_field_value(field_value, real_1d_field=tend_pr_tot_tabs)

    else if (name .eq. "tend_tke_damping_profile_total_local" .and. allocated(tend_pr_tot_tke)) then
      call set_published_field_value(field_value, real_1d_field=tend_pr_tot_tke)

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

end module damping_mod

