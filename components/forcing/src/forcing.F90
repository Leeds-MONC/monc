!> Forcing, both subsidence and large scale
module forcing_mod
  use monc_component_mod, only : COMPONENT_ARRAY_FIELD_TYPE, COMPONENT_DOUBLE_DATA_TYPE, component_descriptor_type, &
       component_field_value_type, component_field_information_type
  use grids_mod, only : Z_INDEX, Y_INDEX, X_INDEX
  use state_mod, only : model_state_type
  use datadefn_mod, only : DEFAULT_PRECISION, STRING_LENGTH
  use optionsdatabase_mod, only : options_get_integer, options_get_logical, options_get_real, options_get_array_size, &
     options_get_logical_array, options_get_real_array, options_get_string_array, options_get_string
  use interpolation_mod, only: piecewise_linear_1d
  use q_indices_mod, only: get_q_index, standard_q_names
  use science_constants_mod, only: seconds_in_a_day
  use naming_conventions_mod
  use registry_mod, only : is_component_enabled
  use logging_mod, only : LOG_ERROR, log_master_log

  implicit none

#ifndef TEST_MODE
  private
#endif

  integer, parameter :: DIVERGENCE=0 ! Input for subsidence forcing is a divergence profile
  integer, parameter :: SUBSIDENCE=1 ! Input for subsidence forcing is the subsidence velocity profile

  integer, parameter :: TENDENCY=0   ! Input for large-scale forcing: values are tendencies (time derivatives)
  integer, parameter :: RELAXATION=1 ! Input for large-scale forcing: values are target values to relax to over timescale
  integer, parameter :: INCREMENT=2  ! Input for large-scale forcing: values are increments (deltas) over timescale

  real(kind=DEFAULT_PRECISION), allocatable :: theta_profile(:) ! Local profile to be used in the subsidence calculation
  real(kind=DEFAULT_PRECISION), allocatable :: q_profile(:) ! Local profile to be used in the subsidence calculation
  real(kind=DEFAULT_PRECISION), allocatable :: u_profile(:) ! Local profile to be used in the subsidence calculation
  real(kind=DEFAULT_PRECISION), allocatable :: v_profile(:) ! Local profile to be used in the subsidence calculation

  real(kind=DEFAULT_PRECISION), allocatable :: dtheta_profile(:) ! Local profile to be used in time-indpendent forcing
  real(kind=DEFAULT_PRECISION), allocatable :: dq_profile(:) ! Local profile to be used in the time-indpendent forcing
  real(kind=DEFAULT_PRECISION), allocatable :: du_profile(:) ! Local profile to be used in the time-indpendent forcing
  real(kind=DEFAULT_PRECISION), allocatable :: dv_profile(:) ! Local profile to be used in the time-indpendent forcing
  real(kind=DEFAULT_PRECISION), allocatable :: du_profile_diag(:), dv_profile_diag(:), dtheta_profile_diag(:), &
       dq_profile_diag(:,:)

  real(kind=DEFAULT_PRECISION) :: forcing_timescale_theta ! Timescale for forcing of theta
  real(kind=DEFAULT_PRECISION) :: forcing_timescale_q     ! Timescale for forcing of q
  real(kind=DEFAULT_PRECISION) :: forcing_timescale_u     ! Timescale for forcing of u
  real(kind=DEFAULT_PRECISION) :: forcing_timescale_v     ! Timescale for forcing of v

  logical :: l_constant_forcing_theta ! Use a time-independent forcing for theta
  logical :: l_constant_forcing_q     ! Use a time-independent forcing for q
  logical :: l_constant_forcing_u     ! Use a time-independent forcing for u
  logical :: l_constant_forcing_v     ! Use a time-independent forcing for v

  integer :: constant_forcing_type_theta=TENDENCY ! Method for large-scale forcing of theta
  integer :: constant_forcing_type_q=TENDENCY     ! Method for large-scale forcing of q
  integer :: constant_forcing_type_u=RELAXATION   ! Method for large-scale forcing of u
  integer :: constant_forcing_type_v=RELAXATION   ! Method for large-scale forcing of v

  logical :: l_subs_pl_theta ! if .true. then subsidence applied to theta field
  logical :: l_subs_pl_q     ! if .true. then subsidence applied to q fields
  
  logical :: l_subs_local_theta ! if .true. then subsidence applied locally (i.e. not with mean fields) to theta field
  logical :: l_subs_local_q     ! if .true. then subsidence applied locally (i.e. not with mean fields) to q fields

  character(len=STRING_LENGTH), dimension(:), allocatable :: names_force_pl_q  ! names of q variables to force

  public forcing_get_descriptor

contains

  !> Provides the component descriptor for the core to register
  !! @returns The descriptor describing this component
  type(component_descriptor_type) function forcing_get_descriptor()
    forcing_get_descriptor%name="forcing"
    forcing_get_descriptor%version=0.1
    forcing_get_descriptor%initialisation=>init_callback
    forcing_get_descriptor%timestep=>timestep_callback
    forcing_get_descriptor%finalisation=>finalisation_callback

    forcing_get_descriptor%field_value_retrieval=>field_value_retrieval_callback
    forcing_get_descriptor%field_information_retrieval=>field_information_retrieval_callback
    allocate(forcing_get_descriptor%published_fields(8))

    forcing_get_descriptor%published_fields(1)="u_subsidence"
    forcing_get_descriptor%published_fields(2)="v_subsidence"
    forcing_get_descriptor%published_fields(3)="th_subsidence"
    forcing_get_descriptor%published_fields(4)="q_subsidence"
    forcing_get_descriptor%published_fields(5)="u_large_scale"
    forcing_get_descriptor%published_fields(6)="v_large_scale"
    forcing_get_descriptor%published_fields(7)="th_large_scale"
    forcing_get_descriptor%published_fields(8)="q_large_scale"
  end function forcing_get_descriptor

  !> Field information retrieval callback, this returns information for a specific components published field
  !! @param current_state Current model state
  !! @param name The name of the field to retrieve information for
  !! @param field_information Populated with information about the field
  subroutine field_information_retrieval_callback(current_state, name, field_information)
    type(model_state_type), target, intent(inout) :: current_state
    character(len=*), intent(in) :: name
    type(component_field_information_type), intent(out) :: field_information

    field_information%field_type=COMPONENT_ARRAY_FIELD_TYPE
    field_information%data_type=COMPONENT_DOUBLE_DATA_TYPE
    field_information%number_dimensions=1
    field_information%dimension_sizes(1)=current_state%local_grid%size(Z_INDEX)

    if (name .eq. "u_subsidence") then
      field_information%enabled=current_state%u%active .and. l_subs_pl_theta .and. &
           allocated(current_state%global_grid%configuration%vertical%olzubar)
    else if (name .eq. "v_subsidence") then
      field_information%enabled=current_state%v%active .and. l_subs_pl_theta .and. &
           allocated(current_state%global_grid%configuration%vertical%olzvbar)
    else if (name .eq. "th_subsidence") then
      field_information%enabled=current_state%th%active .and. l_subs_pl_theta .and. &
           allocated(current_state%global_grid%configuration%vertical%olzthbar)
    else if (name .eq. "q_subsidence") then
      field_information%number_dimensions=2
      field_information%dimension_sizes(2)=current_state%number_q_fields
      field_information%enabled=current_state%number_q_fields .gt. 0 .and. l_subs_pl_q .and. &
           allocated(current_state%global_grid%configuration%vertical%olzqbar)
    else if (name .eq. "u_large_scale") then
      field_information%enabled=current_state%u%active .and. l_constant_forcing_u
    else if (name .eq. "v_large_scale") then
      field_information%enabled=current_state%v%active .and. l_constant_forcing_v
    else if (name .eq. "th_large_scale") then
      field_information%enabled=current_state%th%active .and. l_constant_forcing_theta
    else if (name .eq. "q_large_scale") then
      field_information%number_dimensions=2
      field_information%dimension_sizes(2)=current_state%number_q_fields
      field_information%enabled=current_state%number_q_fields .gt. 0 .and. l_constant_forcing_q
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

    integer :: k, n, column_size

    column_size=current_state%local_grid%size(Z_INDEX)

    if (name .eq. "u_subsidence") then
      allocate(field_value%real_1d_array(column_size))
      do k=2, column_size-1
        field_value%real_1d_array(k)=0.0_DEFAULT_PRECISION-2.0*&
             current_state%global_grid%configuration%vertical%w_subs(k-1)*&
             current_state%global_grid%configuration%vertical%tzc1(k)*(&
             current_state%global_grid%configuration%vertical%olzubar(k)-&
             current_state%global_grid%configuration%vertical%olzubar(k-1))+&
             current_state%global_grid%configuration%vertical%w_subs(k)*&
             current_state%global_grid%configuration%vertical%tzc2(k)*&
             (current_state%global_grid%configuration%vertical%olzubar(k+1)-&
             current_state%global_grid%configuration%vertical%olzubar(k))
      end do
    else if (name .eq. "v_subsidence") then
      allocate(field_value%real_1d_array(column_size))
      do k=2, column_size-1
        field_value%real_1d_array(k)=0.0_DEFAULT_PRECISION-2.0*&
             current_state%global_grid%configuration%vertical%w_subs(k-1)*&
             current_state%global_grid%configuration%vertical%tzc1(k)*(&
             current_state%global_grid%configuration%vertical%olzvbar(k)-&
             current_state%global_grid%configuration%vertical%olzvbar(k-1))+&
             current_state%global_grid%configuration%vertical%w_subs(k)*&
             current_state%global_grid%configuration%vertical%tzc2(k)*(&
             current_state%global_grid%configuration%vertical%olzvbar(k+1)-&
             current_state%global_grid%configuration%vertical%olzvbar(k))
      end do
    else if (name .eq. "th_subsidence") then
      allocate(field_value%real_1d_array(column_size))
      do k=2, column_size-1
        field_value%real_1d_array(k)=0.0_DEFAULT_PRECISION-2.0*&
             current_state%global_grid%configuration%vertical%w_subs(k-1)*&
             current_state%global_grid%configuration%vertical%tzc1(k)*(&
             current_state%global_grid%configuration%vertical%olzthbar(k)-&
             current_state%global_grid%configuration%vertical%olzthbar(k-1)+&
             current_state%global_grid%configuration%vertical%dthref(k-1))+&
             current_state%global_grid%configuration%vertical%w_subs(k)*&
             current_state%global_grid%configuration%vertical%tzc2(k)*(&
             current_state%global_grid%configuration%vertical%olzthbar(k+1)-&
             current_state%global_grid%configuration%vertical%olzthbar(k)+&
             current_state%global_grid%configuration%vertical%dthref(k))
      end do
    else if (name .eq. "q_subsidence") then
      allocate(field_value%real_2d_array(column_size, current_state%number_q_fields))
      do n=1, current_state%number_q_fields
        do k=2, column_size-1
          field_value%real_2d_array(k,n)=0.0_DEFAULT_PRECISION-2.0*&
               current_state%global_grid%configuration%vertical%w_subs(k-1)*&
               current_state%global_grid%configuration%vertical%tzc1(k)*(&
               current_state%global_grid%configuration%vertical%olzqbar(k,n)-&
               current_state%global_grid%configuration%vertical%olzqbar(k-1,n))+&
               current_state%global_grid%configuration%vertical%w_subs(k)&
               *current_state%global_grid%configuration%vertical%tzc2(k)*(&
               current_state%global_grid%configuration%vertical%olzqbar(k+1,n)-&
               current_state%global_grid%configuration%vertical%olzqbar(k,n))
        end do
      end do
    else if (name .eq. "u_large_scale") then
      allocate(field_value%real_1d_array(column_size))
      field_value%real_1d_array=get_averaged_diagnostics(current_state, du_profile_diag)
    else if (name .eq. "v_large_scale") then
      allocate(field_value%real_1d_array(column_size))
      field_value%real_1d_array=get_averaged_diagnostics(current_state, dv_profile_diag)
    else if (name .eq. "th_large_scale") then
      allocate(field_value%real_1d_array(column_size))
      field_value%real_1d_array=get_averaged_diagnostics(current_state, dtheta_profile_diag)
    else if (name .eq. "q_large_scale") then
      allocate(field_value%real_2d_array(column_size, current_state%number_q_fields))
      do n=1, current_state%number_q_fields
        field_value%real_2d_array(:,n)=get_averaged_diagnostics(current_state, dq_profile_diag(:,n))
      end do
    end if
  end subroutine field_value_retrieval_callback  

  !> Initialises the forcing data structures
  subroutine init_callback(current_state)
    type(model_state_type), target, intent(inout) :: current_state

    integer :: nq_force ! The number of q fields apply large-scale time-independent forcing
    integer :: nzq      ! The number of input levels for subsidence/divergence profile
    integer :: i,n  ! loop counters
    integer :: iq       ! temporary q varible index
 
    ! Input arrays for subsidence profile
    real(kind=DEFAULT_PRECISION), dimension(:), allocatable :: f_subs_pl  ! subsidence node for q variables
    real(kind=DEFAULT_PRECISION), dimension(:), allocatable :: z_subs_pl  ! subsidence node height values for q variables

    ! Input arrays for time-independent forcing
    real(kind=DEFAULT_PRECISION), dimension(:, :), allocatable :: f_force_pl_q   ! Forcing values for q variables
    real(kind=DEFAULT_PRECISION), dimension(:), allocatable :: z_force_pl_q      ! Forcing height values for q variables
    real(kind=DEFAULT_PRECISION), dimension(:), allocatable :: f_force_pl_theta  ! Forcing values for theta variable
    real(kind=DEFAULT_PRECISION), dimension(:), allocatable :: z_force_pl_theta  ! Forcing height values for theta variable
    real(kind=DEFAULT_PRECISION), dimension(:), allocatable :: f_force_pl_u      ! Forcing values for u variable
    real(kind=DEFAULT_PRECISION), dimension(:), allocatable :: z_force_pl_u      ! Forcing height values for u variable
    real(kind=DEFAULT_PRECISION), dimension(:), allocatable :: f_force_pl_v      ! Forcing values for v variable
    real(kind=DEFAULT_PRECISION), dimension(:), allocatable :: z_force_pl_v      ! Forcing height values for v variabl

    integer :: subsidence_input_type=DIVERGENCE  ! Determines if we're reading in a subsidence velocity or divergence
    
    real(kind=DEFAULT_PRECISION), allocatable :: f_force_pl_q_tmp(:) !temporary 1D storage of forcing for q fields
    real(kind=DEFAULT_PRECISION), allocatable :: zgrid(:)            ! z grid to use in interpolation   

    character(len=STRING_LENGTH), dimension(:), allocatable :: units_q_force  ! units of q variable forcing
    character(len=STRING_LENGTH) :: units_theta_force='unset'  ! units of theta variable forcing
    character(len=STRING_LENGTH) :: units_u_force='unset'  ! units of theta variable forcing
    character(len=STRING_LENGTH) :: units_v_force='unset'  ! units of theta variable forcing

    logical :: convert_input_theta_from_temperature=.false. ! If .true. input forcing data is for temperature and should
                                                            ! be converted to theta (potential temerature).

    allocate(theta_profile(current_state%local_grid%size(Z_INDEX)), &
       q_profile(current_state%local_grid%size(Z_INDEX)),           &
       u_profile(current_state%local_grid%size(Z_INDEX)),           &
       v_profile(current_state%local_grid%size(Z_INDEX)))

    allocate(dtheta_profile(current_state%local_grid%size(Z_INDEX)), &
       dq_profile(current_state%local_grid%size(Z_INDEX)),           &
       du_profile(current_state%local_grid%size(Z_INDEX)),           &
       dv_profile(current_state%local_grid%size(Z_INDEX)))

    allocate(du_profile_diag(current_state%local_grid%size(Z_INDEX)), dv_profile_diag(current_state%local_grid%size(Z_INDEX)), &
         dtheta_profile_diag(current_state%local_grid%size(Z_INDEX)), &
         dq_profile_diag(current_state%local_grid%size(Z_INDEX), current_state%number_q_fields))

    allocate(zgrid(current_state%local_grid%size(Z_INDEX)))

    ! Subsidence forcing initialization
    l_subs_pl_theta=options_get_logical(current_state%options_database, "l_subs_pl_theta")
    l_subs_pl_q=options_get_logical(current_state%options_database, "l_subs_pl_q")
    subsidence_input_type=options_get_integer(current_state%options_database, "subsidence_input_type")
    l_subs_local_theta=options_get_logical(current_state%options_database, "subsidence_local_theta")
    l_subs_local_q=options_get_logical(current_state%options_database, "subsidence_local_q")


    if ((l_subs_pl_theta .and. .not. l_subs_local_theta) .or. &
       (l_subs_pl_q .and. .not. l_subs_local_q))then
      if (.not. is_component_enabled(current_state%options_database, "mean_profiles")) then
        call log_master_log(LOG_ERROR, "Damping requires the mean profiles component to be enabled")
      end if
    end if

    if (l_subs_pl_theta .or. l_subs_pl_q)then
      allocate(z_subs_pl(options_get_array_size(current_state%options_database, "z_subs_pl")), &
           f_subs_pl(options_get_array_size(current_state%options_database, "f_subs_pl")))
      call options_get_real_array(current_state%options_database, "z_subs_pl", z_subs_pl)
      call options_get_real_array(current_state%options_database, "f_subs_pl", f_subs_pl)      
      ! Get profiles
      zgrid=current_state%global_grid%configuration%vertical%z(:)
      call piecewise_linear_1d(z_subs_pl(1:size(z_subs_pl)), f_subs_pl(1:size(f_subs_pl)), zgrid, &
         current_state%global_grid%configuration%vertical%w_subs)
      if (subsidence_input_type==DIVERGENCE) then
        current_state%global_grid%configuration%vertical%w_subs(:) = &
           -1.0*current_state%global_grid%configuration%vertical%w_subs(:)*zgrid(:)
      end if
      deallocate(z_subs_pl, f_subs_pl)
    end if

    ! Time independent large-scale forcing (proxy for e.g. advection/radiation)

    ! This probably isn't the right place to be doing this
    if (.not. allocated(current_state%l_forceq))then
      allocate(current_state%l_forceq(current_state%number_q_fields))
      current_state%l_forceq=.false.
    end if

    l_constant_forcing_theta=options_get_logical(current_state%options_database, "l_constant_forcing_theta")
    l_constant_forcing_q=options_get_logical(current_state%options_database, "l_constant_forcing_q")
    l_constant_forcing_u=options_get_logical(current_state%options_database, "l_constant_forcing_u")
    l_constant_forcing_v=options_get_logical(current_state%options_database, "l_constant_forcing_v")

    if (l_constant_forcing_q) then
      allocate(names_force_pl_q(options_get_array_size(current_state%options_database, "names_constant_forcing_q")))
      call options_get_string_array(current_state%options_database, "names_constant_forcing_q", names_force_pl_q)
    end if
    
    if (l_constant_forcing_theta)then
      allocate(z_force_pl_theta(options_get_array_size(current_state%options_database, "z_force_pl_theta")), &
           f_force_pl_theta(options_get_array_size(current_state%options_database, "f_force_pl_theta")))
      call options_get_real_array(current_state%options_database, "z_force_pl_theta", z_force_pl_theta)
      call options_get_real_array(current_state%options_database, "f_force_pl_theta", f_force_pl_theta)
      ! Get profiles
      zgrid=current_state%global_grid%configuration%vertical%zn(:)
      call piecewise_linear_1d(z_force_pl_theta(1:size(z_force_pl_theta)), f_force_pl_theta(1:size(f_force_pl_theta)), zgrid, &
         current_state%global_grid%configuration%vertical%theta_force)

      ! Unit conversions...
      if (convert_input_theta_from_temperature)then ! Input is temperature not theta
        current_state%global_grid%configuration%vertical%theta_force(:) =   &
           current_state%global_grid%configuration%vertical%theta_force(:)* &
           current_state%global_grid%configuration%vertical%prefrcp(:)
      end if

      units_theta_force=options_get_string(current_state%options_database, "units_theta_force")
      select case(trim(units_theta_force))
      case(k_per_day)
        current_state%global_grid%configuration%vertical%theta_force(:) = &
           current_state%global_grid%configuration%vertical%theta_force(:)/seconds_in_a_day
      case default !(k_per_second)
      end select
      deallocate(z_force_pl_theta, f_force_pl_theta)
    end if
                                                           
#ifdef U_ACTIVE
    if (l_constant_forcing_u)then
      allocate(z_force_pl_u(options_get_array_size(current_state%options_database, "z_force_pl_u")), &
           f_force_pl_u(options_get_array_size(current_state%options_database, "f_force_pl_u")))
      call options_get_real_array(current_state%options_database, "z_force_pl_u", z_force_pl_u)
      call options_get_real_array(current_state%options_database, "f_force_pl_u", f_force_pl_u)
      ! Get profiles
      zgrid=current_state%global_grid%configuration%vertical%zn(:)
      call piecewise_linear_1d(z_force_pl_u(1:size(z_force_pl_u)), f_force_pl_u(1:size(f_force_pl_u)), zgrid, &
         current_state%global_grid%configuration%vertical%u_force)

      ! Unit conversions...
      units_u_force=options_get_string(current_state%options_database, "units_u_force")
      select case(trim(units_u_force))
      case(m_per_second_per_day)
        current_state%global_grid%configuration%vertical%u_force(:) = &
           current_state%global_grid%configuration%vertical%u_force(:)/seconds_in_a_day
      case default  !(m_per_second_per_second)
      end select
      deallocate(z_force_pl_u, f_force_pl_u)
    end if
#endif
                                                           
#ifdef V_ACTIVE
    if (l_constant_forcing_v)then
      allocate(z_force_pl_v(options_get_array_size(current_state%options_database, "z_force_pl_v")), &
           f_force_pl_v(options_get_array_size(current_state%options_database, "f_force_pl_v")))
      call options_get_real_array(current_state%options_database, "z_force_pl_v", z_force_pl_v)
      call options_get_real_array(current_state%options_database, "f_force_pl_v", f_force_pl_v)
      ! Get profiles
      zgrid=current_state%global_grid%configuration%vertical%zn(:)
      call piecewise_linear_1d(z_force_pl_v(1:size(z_force_pl_v)), f_force_pl_v(1:size(f_force_pl_v)), zgrid, &
         current_state%global_grid%configuration%vertical%u_force)

      ! Unit conversions...
      units_v_force=options_get_string(current_state%options_database, "units_v_force")
      select case(trim(units_v_force))
      case(m_per_second_per_day)
        current_state%global_grid%configuration%vertical%v_force(:) = &
           current_state%global_grid%configuration%vertical%v_force(:)/seconds_in_a_day
      case default !(m_per_second_per_second)
      end select
      deallocate(z_force_pl_v, f_force_pl_v)
    end if 
#endif   
        
    if (l_constant_forcing_q) then
      nq_force=size(names_force_pl_q)
      allocate(z_force_pl_q(options_get_array_size(current_state%options_database, "z_force_pl_q")))
      call options_get_real_array(current_state%options_database, "z_force_pl_q", z_force_pl_q)
      nzq=size(z_force_pl_q)
      zgrid=current_state%global_grid%configuration%vertical%zn(:)
      allocate(f_force_pl_q_tmp(nq_force*nzq))
      call options_get_real_array(current_state%options_database, "f_force_pl_q", f_force_pl_q_tmp)
      allocate(f_force_pl_q(nzq, nq_force))
      f_force_pl_q(1:nzq, 1:nq_force)=reshape(f_force_pl_q_tmp, (/nzq, nq_force/))

      allocate(units_q_force(options_get_array_size(current_state%options_database, "units_q_force")))
      call options_get_string_array(current_state%options_database, "units_q_force", units_q_force)
      do n=1, nq_force
        iq=get_q_index(trim(names_force_pl_q(n)), 'forcing:time-independent')
        call piecewise_linear_1d(z_force_pl_q(1:nzq), f_force_pl_q(1:nzq,n), zgrid, &
           current_state%global_grid%configuration%vertical%q_force(:,iq))
        
        current_state%l_forceq(iq)=.true.

        ! Unit conversions...
        select case(trim(units_q_force(n)))
        case(kg_per_kg_per_day)
          current_state%global_grid%configuration%vertical%q_force(:,iq) = &
             current_state%global_grid%configuration%vertical%q_force(:,iq)/seconds_in_a_day
        case(g_per_kg_per_day)
          current_state%global_grid%configuration%vertical%q_force(:,iq) = &
             0.001*current_state%global_grid%configuration%vertical%q_force(:,iq)/seconds_in_a_day
        case(g_per_kg_per_second)
          current_state%global_grid%configuration%vertical%q_force(:,iq) = &
             0.001*current_state%global_grid%configuration%vertical%q_force(:,iq)
        case default !(kg_per_kg_per_second)
        end select
      end do
      deallocate(f_force_pl_q_tmp, units_q_force, f_force_pl_q, z_force_pl_q)  
    end if
    deallocate(zgrid)
  end subroutine init_callBack

  !> Called for each data column and will determine the forcing values in x and y which are then applied to the field
  !! source terms
  !! @param current_state The current model state_mod
  subroutine timestep_callback(current_state)
    type(model_state_type), target, intent(inout) :: current_state

    if (current_state%first_timestep_column) then
      du_profile_diag=0.0_DEFAULT_PRECISION
      dv_profile_diag=0.0_DEFAULT_PRECISION
      dtheta_profile_diag=0.0_DEFAULT_PRECISION
      dq_profile_diag=0.0_DEFAULT_PRECISION
    end if    

    if (current_state%halo_column .or. current_state%timestep<3) return
    if (l_subs_pl_theta) then
      call apply_subsidence_to_flow_fields(current_state)
      call apply_subsidence_to_theta(current_state)
    end if
    if (l_subs_pl_q) call apply_subsidence_to_q_fields(current_state)

    if (l_constant_forcing_theta)call apply_time_independent_forcing_to_theta(current_state)                                                         
#ifdef U_ACTIVE
    if (l_constant_forcing_u)call apply_time_independent_forcing_to_u(current_state)
#endif                                                       
#ifdef V_ACTIVE
    if (l_constant_forcing_v)call apply_time_independent_forcing_to_v(current_state)
#endif
    if (l_constant_forcing_q)call apply_time_independent_forcing_to_q(current_state)
  end subroutine timestep_callback

  subroutine apply_subsidence_to_flow_fields(current_state)
    type(model_state_type), intent(inout) :: current_state

    integer :: k
    real(kind=DEFAULT_PRECISION) :: usub, vsub    


    if (l_subs_local_theta)then ! Use local gradients not global means
      u_profile(:)=current_state%zu%data(:,current_state%column_local_y,current_state%column_local_x)
      v_profile(:)=current_state%zv%data(:,current_state%column_local_y,current_state%column_local_x)
    else
      u_profile(:)=current_state%global_grid%configuration%vertical%olzubar(:)
      v_profile(:)=current_state%global_grid%configuration%vertical%olzvbar(:)
    end if

    do k=2,current_state%local_grid%size(Z_INDEX)-1                                                             
#ifdef U_ACTIVE
      usub =  2.0 * (current_state%global_grid%configuration%vertical%w_subs(k-1)*              &
         current_state%global_grid%configuration%vertical%tzc1(k)*(u_profile(k)-u_profile(k-1)) &
           + current_state%global_grid%configuration%vertical%w_subs(k)*                        &
           current_state%global_grid%configuration%vertical%tzc2(k)*                            &
           (u_profile(k+1) - u_profile(k)))
      current_state%su%data(k,current_state%column_local_y,current_state%column_local_x) =      &
           current_state%su%data(k,current_state%column_local_y,current_state%column_local_x) - usub
#endif
#ifdef V_ACTIVE     
      vsub =  2.0 * (current_state%global_grid%configuration%vertical%w_subs(k-1)*              &
         current_state%global_grid%configuration%vertical%tzc1(k)*(v_profile(k)-v_profile(k-1)) &
         + current_state%global_grid%configuration%vertical%w_subs(k)*                          &
         current_state%global_grid%configuration%vertical%tzc2(k)*                              &
         (v_profile(k+1) - v_profile(k)))
      current_state%sv%data(k,current_state%column_local_y,current_state%column_local_x) =      &
           current_state%sv%data(k,current_state%column_local_y,current_state%column_local_x) - vsub
#endif
    end do
    k=current_state%local_grid%size(Z_INDEX)
#ifdef U_ACTIVE
    usub =  2.0 * (current_state%global_grid%configuration%vertical%w_subs(k-1)*                &
         current_state%global_grid%configuration%vertical%tzc1(k)*(u_profile(k)-u_profile(k-1)))
    current_state%su%data(k,current_state%column_local_y,current_state%column_local_x) =        &
         current_state%su%data(k,current_state%column_local_y,current_state%column_local_x) - usub
#endif
#ifdef V_ACTIVE
    vsub =  2.0 * (current_state%global_grid%configuration%vertical%w_subs(k-1)*                &
         current_state%global_grid%configuration%vertical%tzc1(k)*(v_profile(k)-v_profile(k-1)))
    current_state%sv%data(k,current_state%column_local_y,current_state%column_local_x) =        &
         current_state%sv%data(k,current_state%column_local_y,current_state%column_local_x) - vsub
#endif
  end subroutine apply_subsidence_to_flow_fields

  subroutine apply_subsidence_to_theta(current_state)
    type(model_state_type), intent(inout) :: current_state

    integer :: k
    real(kind=DEFAULT_PRECISION) :: thsub

    if (l_subs_local_theta)then ! Use local gradients not global means
      theta_profile(:)=current_state%zth%data(:,current_state%column_local_y,current_state%column_local_x) &
         + current_state%global_grid%configuration%vertical%thref(:)
    else
      theta_profile(:)=current_state%global_grid%configuration%vertical%olzthbar(:) &
         + current_state%global_grid%configuration%vertical%thref(:)
    end if

    do k=2,current_state%local_grid%size(Z_INDEX)-1
      thsub = current_state%global_grid%configuration%vertical%w_subs(k)*   &
         (theta_profile(k+1) - theta_profile(k))*                           &
         current_state%global_grid%configuration%vertical%rdzn(k+1)
      
      current_state%sth%data(k,current_state%column_local_y,current_state%column_local_x) = &
           current_state%sth%data(k,current_state%column_local_y,current_state%column_local_x) - thsub
    end do
    k=current_state%local_grid%size(Z_INDEX)
    thsub = current_state%global_grid%configuration%vertical%w_subs(k)* &
       (theta_profile(k) - theta_profile(k-1))*                         &
         current_state%global_grid%configuration%vertical%rdzn(k)

    current_state%sth%data(k,current_state%column_local_y,current_state%column_local_x) = &
         current_state%sth%data(k,current_state%column_local_y,current_state%column_local_x) - thsub
  end subroutine apply_subsidence_to_theta

  subroutine apply_subsidence_to_q_fields(current_state)
    type(model_state_type), intent(inout) :: current_state

    integer :: k, n
    real(kind=DEFAULT_PRECISION) :: qsub


    do n=1,current_state%number_q_fields
      if (l_subs_local_q)then ! Use local gradients not global means
        q_profile(:)=current_state%zq(n)%data(:,current_state%column_local_y,current_state%column_local_x)
      else
        q_profile(:)=current_state%global_grid%configuration%vertical%olzqbar(:,n)
      end if
      do k=2,current_state%local_grid%size(Z_INDEX)-1
        qsub = current_state%global_grid%configuration%vertical%w_subs(k)*  &
           (q_profile(k+1) - q_profile(k))*                               &
           current_state%global_grid%configuration%vertical%rdzn(k+1)
        current_state%sq(n)%data(k,current_state%column_local_y,current_state%column_local_x) = &
           current_state%sq(n)%data(k,current_state%column_local_y,current_state%column_local_x) - qsub
      end do
      k=current_state%local_grid%size(Z_INDEX)
      qsub = current_state%global_grid%configuration%vertical%w_subs(k)*    &
         (q_profile(k) - q_profile(k-1))*                                   &
         current_state%global_grid%configuration%vertical%rdzn(k)
      current_state%sq(n)%data(k,current_state%column_local_y,current_state%column_local_x) = &
         current_state%sq(n)%data(k,current_state%column_local_y,current_state%column_local_x) - qsub
    end do
  end subroutine apply_subsidence_to_q_fields

  subroutine apply_time_independent_forcing_to_theta(current_state)
    type(model_state_type), intent(inout) :: current_state

    integer :: k
    real(kind=DEFAULT_PRECISION) :: dtm_scale

    if (constant_forcing_type_theta==TENDENCY)then
      dtm_scale=current_state%dtm
    else !  constant_forcing_type_theta==(RELAXATION or INCREMENT)
      dtm_scale=current_state%dtm/forcing_timescale_theta
    end if
    
    if (constant_forcing_type_theta==RELAXATION)then
      dtheta_profile(:)=dtm_scale * (current_state%global_grid%configuration%vertical%theta_force(:) - &
         current_state%zth%data(:,current_state%column_local_y,current_state%column_local_x) -       &
         current_state%global_grid%configuration%vertical%thref(:))
    else !  constant_forcing_type_theta==(TENDENCY or INCREMENT)
      dtheta_profile(:)=dtm_scale * current_state%global_grid%configuration%vertical%theta_force(:)
    end if

    dtheta_profile_diag=dtheta_profile_diag+dtheta_profile

    do k=2,current_state%local_grid%size(Z_INDEX)-1
      current_state%sth%data(k,current_state%column_local_y,current_state%column_local_x) = &
           current_state%sth%data(k,current_state%column_local_y,current_state%column_local_x) &
           + dtheta_profile(k)
    end do

  end subroutine apply_time_independent_forcing_to_theta

  subroutine apply_time_independent_forcing_to_q(current_state)
    type(model_state_type), intent(inout) :: current_state

    integer :: n, k
    real(kind=DEFAULT_PRECISION) :: dtm_scale

    do n=1,current_state%number_q_fields
      if (current_state%l_forceq(n))then
        if (constant_forcing_type_q==TENDENCY)then
          dtm_scale=current_state%dtm
        else !  constant_forcing_type_q==(RELAXATION or INCREMENT)
          dtm_scale=current_state%dtm/forcing_timescale_q
        end if
        
        if (constant_forcing_type_q==RELAXATION)then
          dq_profile(:)=dtm_scale * (current_state%global_grid%configuration%vertical%q_force(:,n) - &
             current_state%zq(n)%data(:,current_state%column_local_y,current_state%column_local_x))
        else !  constant_forcing_type_q==(TENDENCY or INCREMENT)
          dq_profile(:)=dtm_scale * current_state%global_grid%configuration%vertical%q_force(:,n)
        end if

        dq_profile_diag(:,n)=dq_profile_diag(:,n)+dq_profile

        do k=2,current_state%local_grid%size(Z_INDEX)-1
          current_state%sq(n)%data(k,current_state%column_local_y,current_state%column_local_x) = &
             current_state%sq(n)%data(k,current_state%column_local_y,current_state%column_local_x) &
             + dq_profile(k)
        end do
      end if
    end do

  end subroutine apply_time_independent_forcing_to_q

  subroutine apply_time_independent_forcing_to_u(current_state)
    type(model_state_type), intent(inout) :: current_state

    integer :: k
    real(kind=DEFAULT_PRECISION) :: dtm_scale

    if (constant_forcing_type_u==TENDENCY)then
      dtm_scale=current_state%dtm
    else !  constant_forcing_type_u==(RELAXATION or INCREMENT)
      dtm_scale=current_state%dtm/forcing_timescale_u
    end if
    
    if (constant_forcing_type_u==RELAXATION)then
      du_profile(:)=dtm_scale * (current_state%global_grid%configuration%vertical%u_force(:) - &
         current_state%zu%data(:,current_state%column_local_y,current_state%column_local_x))
    else !  constant_forcing_type_u==(TENDENCY or INCREMENT)
      du_profile(:)=dtm_scale * current_state%global_grid%configuration%vertical%u_force(:)
    end if

    du_profile_diag=du_profile_diag+du_profile

    do k=2,current_state%local_grid%size(Z_INDEX)-1
      current_state%su%data(k,current_state%column_local_y,current_state%column_local_x) = &
           current_state%su%data(k,current_state%column_local_y,current_state%column_local_x) &
           + du_profile(k)
    end do

  end subroutine apply_time_independent_forcing_to_u

  subroutine apply_time_independent_forcing_to_v(current_state)
    type(model_state_type), intent(inout) :: current_state

    integer :: k
    real(kind=DEFAULT_PRECISION) :: dtm_scale

    if (constant_forcing_type_v==TENDENCY)then
      dtm_scale=current_state%dtm
    else !  constant_forcing_type_v==(RELAXATION or INCREMENT)
      dtm_scale=current_state%dtm/forcing_timescale_v
    end if
    
    if (constant_forcing_type_v==RELAXATION)then
      dv_profile(:)=dtm_scale * (current_state%global_grid%configuration%vertical%v_force(:) - &
         current_state%zv%data(:,current_state%column_local_y,current_state%column_local_x))
    else !  constant_forcing_type_v==(TENDENCY or INCREMENT)
      dv_profile(:)=dtm_scale * current_state%global_grid%configuration%vertical%v_force(:)
    end if

    dv_profile_diag=dv_profile_diag+dv_profile

    do k=2,current_state%local_grid%size(Z_INDEX)-1
      current_state%sv%data(k,current_state%column_local_y,current_state%column_local_x) = &
           current_state%sv%data(k,current_state%column_local_y,current_state%column_local_x) &
           + dv_profile(k)
    end do
  end subroutine apply_time_independent_forcing_to_v

  !> Finalises the component 
  !! @current_state Current model state
  subroutine finalisation_callback(current_state)
    type(model_state_type), target, intent(inout) :: current_state
    deallocate(theta_profile, q_profile, u_profile, v_profile)
    deallocate(dtheta_profile, dq_profile, du_profile, dv_profile)
    if (allocated(du_profile_diag)) deallocate(du_profile_diag)
    if (allocated(dv_profile_diag)) deallocate(dv_profile_diag)
    if (allocated(dtheta_profile_diag)) deallocate(dtheta_profile_diag)
    if (allocated(dq_profile_diag)) deallocate(dq_profile_diag)
  end subroutine finalisation_callback

  !> Averages some diagnostic values across all local horizontal points
  !! @param current_state Current model state
  !! @param diagnostics_summed The diagnostics, which are summed up values from each local point, which need to be averaged
  !! @returns The averaged diagnostics across all local horizontal points
  function get_averaged_diagnostics(current_state, diagnostics_summed)
    type(model_state_type), target, intent(inout) :: current_state
    real(kind=DEFAULT_PRECISION), dimension(:), intent(in) :: diagnostics_summed
    real(kind=DEFAULT_PRECISION), dimension(size(diagnostics_summed)) :: get_averaged_diagnostics

    integer :: horizontal_points

    horizontal_points=current_state%local_grid%size(X_INDEX) * current_state%local_grid%size(Y_INDEX)

    get_averaged_diagnostics(:)=diagnostics_summed(:)/horizontal_points
  end function get_averaged_diagnostics
end module forcing_mod
