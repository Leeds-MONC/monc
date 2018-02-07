!> Manages the grid based upon the model state_mod
module gridmanager_mod
  use datadefn_mod, only : STRING_LENGTH
  use monc_component_mod, only : component_descriptor_type
  use state_mod, only : model_state_type
  use optionsdatabase_mod, only : options_get_integer, options_get_logical, options_get_real, &
     options_get_logical_array, options_get_real_array, options_get_string_array, options_get_array_size, options_get_string
  use grids_mod, only : vertical_grid_configuration_type, X_INDEX, Y_INDEX, Z_INDEX
  use logging_mod, only : LOG_INFO, LOG_ERROR, log_master_log, log_log
  use conversions_mod, only : conv_to_string
  use datadefn_mod, only : DEFAULT_PRECISION
  use science_constants_mod, only : von_karman_constant, z0, z0th, cp, r_over_cp, r, g, rlvap_over_cp
  use saturation_mod, only : qsaturation, dqwsatdt
  use q_indices_mod, only: get_q_index, standard_q_names
  use interpolation_mod, only: piecewise_linear_1d

  implicit none

#ifndef TEST_MODE
  private
#endif
  
  ! 1 = No adjustment
  ! 2 = ensure P0=PSF by adjusting THREF profile by constant factor
  ! 3 = ensure P0=PSF by adjusting PSF (not advised)
  ! 4 = ensure P0=PSF by adjusting PTOP
  integer, parameter :: ANELASTIC_PROFILE_MODE=4
  real, parameter :: DEFAULT_SPACING = 1.E9  !< The default spacing used if no grid is active in a specific dimension
  real(kind=DEFAULT_PRECISION), dimension(:,:), allocatable :: qinit

  public gridmanager_get_descriptor

contains

  !> Provides the descriptor back to the caller and is used in component registration
  !! @returns The GridManager component descriptor
  type(component_descriptor_type) function gridmanager_get_descriptor()
    gridmanager_get_descriptor%name="grid_manager"
    gridmanager_get_descriptor%version=0.1
    gridmanager_get_descriptor%initialisation=>initialise_callback
    gridmanager_get_descriptor%finalisation=>finalise_callback
  end function gridmanager_get_descriptor

  !> Called during initialisation and will initialise the horizontal and vertical grid configurations
  !! Note that the model state_mod (from a checkpoint or external file) must have been initialised already
  !! @param current_state The current model state_mod
  subroutine initialise_callback(current_state)
    type(model_state_type), target, intent(inout) :: current_state

    integer :: dimensions

    if (.not. current_state%initialised) then
      call log_log(LOG_ERROR, "Must initialise the model state_mod before constructing the grid properties")
    end if

    call initialise_horizontalgrid_configuration_types(current_state)
    call initialise_verticalgrid_configuration_type(current_state)
    dimensions=1
    if (current_state%global_grid%active(X_INDEX)) dimensions = dimensions+1
    if (current_state%global_grid%active(Y_INDEX)) dimensions = dimensions+1   
    call log_master_log(LOG_INFO, trim(conv_to_string(dimensions))//"D system; z="//&
         trim(conv_to_string(current_state%global_grid%size(Z_INDEX)))//", y="//&
         trim(conv_to_string(current_state%global_grid%size(Y_INDEX)))//", x="//&
         trim(conv_to_string(current_state%global_grid%size(X_INDEX))))
  end subroutine initialise_callback

  !> Called as MONC exits, for good practice this will deallocate the memory used for the grids
  !! @param current_state The current model state_mod
  subroutine finalise_callback(current_state)
    type(model_state_type), target, intent(inout) :: current_state

    type(vertical_grid_configuration_type) :: vertical_grid

    vertical_grid=current_state%global_grid%configuration%vertical

    deallocate(vertical_grid%z, vertical_grid%zn, vertical_grid%dz, vertical_grid%dzn, vertical_grid%czb, vertical_grid%cza, &
         vertical_grid%czg, vertical_grid%czh, vertical_grid%rdz, vertical_grid%rdzn, vertical_grid%tzc1, vertical_grid%tzc2,&
         vertical_grid%tzd1, vertical_grid%tzd2, vertical_grid%thref, vertical_grid%theta_init,  vertical_grid%temp_init, &
         vertical_grid%rh_init, vertical_grid%tref, &
         vertical_grid%prefn, vertical_grid%pdiff, vertical_grid%prefrcp, vertical_grid%rprefrcp, vertical_grid%rho, &
         vertical_grid%rhon, vertical_grid%tstarpr, vertical_grid%qsat, vertical_grid%dqsatdt, vertical_grid%qsatfac, &
         vertical_grid%dthref, vertical_grid%rneutml, vertical_grid%rneutml_sq, vertical_grid%buoy_co, &
         vertical_grid%u_init, vertical_grid%v_init, vertical_grid%q_init,                              &
         vertical_grid%q_rand, vertical_grid%theta_rand, vertical_grid%w_subs, vertical_grid%w_rand, &
         vertical_grid%q_force, vertical_grid%theta_force, vertical_grid%u_force, vertical_grid%v_force &
         )
  end subroutine finalise_callback  

  !> Will initialise the vertical grid configuration
  !! @param current_state The current model state_mod
  subroutine initialise_verticalgrid_configuration_type(current_state)
    type(model_state_type), intent(inout) :: current_state

    call allocate_vertical_grid_data(current_state%global_grid%configuration%vertical, &
       current_state%global_grid%size(Z_INDEX), current_state%number_q_fields )
    call set_up_and_smooth_grid(current_state%global_grid%configuration%vertical, &
         current_state%global_grid%configuration%vertical%kgd, current_state%global_grid%configuration%vertical%hgd, &
         size(current_state%global_grid%configuration%vertical%kgd), current_state%global_grid%size(Z_INDEX), &
         current_state%global_grid%top(Z_INDEX), options_get_integer(current_state%options_database, "nsmth"), &
         current_state%origional_vertical_grid_setup, current_state%continuation_run)
    call set_vertical_reference_profile(current_state, current_state%global_grid%configuration%vertical, &
         current_state%global_grid%size(Z_INDEX))
    
  end subroutine initialise_verticalgrid_configuration_type

  !> Sets up the vertical grid reference profile at each point
  !! @param current_state The current model state_mod
  !! @param vertical_grid The vertical grid that we are working on
  !! @param kkp Number of grid points in a vertical column
  subroutine set_vertical_reference_profile(current_state, vertical_grid, kkp)
    type(model_state_type), intent(inout) :: current_state
    type(vertical_grid_configuration_type), intent(inout) :: vertical_grid
    integer, intent(in) :: kkp

    integer :: k

    call calculate_initial_profiles(current_state, vertical_grid)
    call set_up_vertical_reference_properties(current_state, vertical_grid, current_state%global_grid%size(Z_INDEX))
    call set_anelastic_pressure(current_state)
    ! 
    call set_qv_init_from_rh(current_state)

    do k=2,kkp-1
      ! for diffusion onto p-level from below
      vertical_grid%czb(k)=(vertical_grid%rho(k-1)/vertical_grid%rhon(k))/(vertical_grid%dz(k)*vertical_grid%dzn(k))
      ! for diffusion onto p-level from above
      vertical_grid%cza(k)=(vertical_grid%rho(k)/vertical_grid%rhon(k))/(vertical_grid%dz(k)*vertical_grid%dzn(k+1))
      vertical_grid%czg(k)=-vertical_grid%czb(k)-vertical_grid%cza(k)
      if (k .gt. 2) vertical_grid%czh(k)=vertical_grid%czb(k)*vertical_grid%cza(k-1)
    end do
    do k=2,kkp-1
      ! advection onto p-level from below
      vertical_grid%tzc1(k)=0.25_DEFAULT_PRECISION*vertical_grid%rdz(k)*vertical_grid%rho(k-1)/vertical_grid%rhon(k) 
      ! advection onto p-level from above
      vertical_grid%tzc2(k)=0.25_DEFAULT_PRECISION*vertical_grid%rdz(k)*vertical_grid%rho(k)/vertical_grid%rhon(k) 
    end do
    do k=2,kkp-1
      ! advection onto w-level (K) from below
      vertical_grid%tzd1(k)=0.25_DEFAULT_PRECISION*vertical_grid%rdzn(k+1)*vertical_grid%rhon(k)/vertical_grid%rho(k)
      ! advection onto w-level (K) from above
      vertical_grid%tzd2(k)=0.25_DEFAULT_PRECISION*vertical_grid%rdzn(k+1)*vertical_grid%rhon(k+1)/vertical_grid%rho(k)
    end do
    k=kkp
    vertical_grid%czb(k)=(vertical_grid%rho(k-1)/vertical_grid%rhon(k))/(vertical_grid%dz(k)*vertical_grid%dzn(k))
    vertical_grid%cza(k)=0.0_DEFAULT_PRECISION
    vertical_grid%czg(k)=-vertical_grid%czb(k)
    vertical_grid%czh(k)=vertical_grid%czb(k)*vertical_grid%cza(k-1)
    vertical_grid%tzc2(k)=0.25_DEFAULT_PRECISION*vertical_grid%rdz(k)*vertical_grid%rho(k)/vertical_grid%rhon(k)
    vertical_grid%tzc1(k)=0.25_DEFAULT_PRECISION*vertical_grid%rdz(k)*vertical_grid%rho(k-1)/vertical_grid%rhon(k) 
    vertical_grid%czn=vertical_grid%dzn(2)*0.5_DEFAULT_PRECISION
    vertical_grid%zlogm=log(1.0_DEFAULT_PRECISION+vertical_grid%zn(2)/z0)
    vertical_grid%zlogth=log((vertical_grid%zn(2)+z0)/z0th)
    vertical_grid%vk_on_zlogm=von_karman_constant/vertical_grid%zlogm
    call setup_reference_state_liquid_water_temperature_and_saturation(&
         current_state, vertical_grid, current_state%global_grid%size(Z_INDEX))
    call calculate_mixing_length_for_neutral_case(current_state, vertical_grid, current_state%global_grid%size(Z_INDEX))
    call set_buoyancy_coefficient(current_state, vertical_grid, current_state%global_grid%size(Z_INDEX))
  end subroutine set_vertical_reference_profile

  !> Calculates the initial profiles for U, V, TH & Q if required
  !! @param current_state The current model state_mod
  !! @param vertical_grid The vertical grid that we are working on
  !! @param kkp Number of grid points in a vertical column
  subroutine calculate_initial_profiles(current_state, vertical_grid)
    type(model_state_type), intent(inout) :: current_state
    type(vertical_grid_configuration_type), intent(inout) :: vertical_grid

    integer :: nq_init ! The number of q fields to initialize
    integer :: nzq     ! The number of input levels for q_init
    integer :: i,j,n, k ! loop counters
    integer :: iq  ! temporary q varible index

    real(kind=DEFAULT_PRECISION), dimension(:,:), allocatable :: f_init_pl_q       ! Initial node values for q variables
    real(kind=DEFAULT_PRECISION), dimension(:), allocatable :: z_init_pl_q      ! Initial node height values for q variables
    real(kind=DEFAULT_PRECISION), dimension(:), allocatable :: f_init_pl_theta  ! Initial node values for potential temperature variable
    real(kind=DEFAULT_PRECISION), dimension(:), allocatable :: z_init_pl_theta  ! Initial node height values for potential temperature variable
    real(kind=DEFAULT_PRECISION), dimension(:), allocatable :: f_init_pl_u      ! Initial node values for u variable
    real(kind=DEFAULT_PRECISION), dimension(:), allocatable :: z_init_pl_u      ! Initial node height values for u variable
    real(kind=DEFAULT_PRECISION), dimension(:), allocatable :: f_init_pl_v      ! Initial node values for v variable
    real(kind=DEFAULT_PRECISION), dimension(:), allocatable :: z_init_pl_v      ! Initial node height values for v variable

    real(kind=DEFAULT_PRECISION), dimension(:), allocatable :: f_thref   ! Initial node values for thref
    real(kind=DEFAULT_PRECISION), dimension(:), allocatable :: z_thref   ! Initial node height values for thref

    logical :: l_init_pl_u     ! if .true. then initialize u field
    logical :: l_init_pl_v     ! if .true. then initialize v field
    logical :: l_init_pl_theta ! if .true. then initialize potential temperature field
    logical :: l_init_pl_rh    ! if .true. then initialize relative humidity field
    logical :: l_init_pl_q     ! if .true. then initialize q fields
    logical :: l_thref         ! if .true. then initialize thref profile (overrides thref0)
    logical :: l_matchthref    ! if .true. then initialize thref to be the same as theta_init

    character(len=STRING_LENGTH), dimension(:), allocatable :: names_init_pl_q ! names of q variables to initialize
    
    real(kind=DEFAULT_PRECISION), allocatable :: f_init_pl_q_tmp(:) !temporary 1D storage of initial q field
    real(kind=DEFAULT_PRECISION), allocatable :: zgrid(:)  ! z grid to use in interpolation

    real(kind=DEFAULT_PRECISION) :: zztop ! top of the domain
    real(kind=DEFAULT_PRECISION) :: qsat

    allocate(zgrid(current_state%local_grid%local_domain_end_index(Z_INDEX)))
    
    zztop = current_state%global_grid%top(Z_INDEX)

    ! Initialize everything to zero.  This won't make sense for theta.
    vertical_grid%q_init = 0.0_DEFAULT_PRECISION
    vertical_grid%u_init = 0.0_DEFAULT_PRECISION
    vertical_grid%v_init = 0.0_DEFAULT_PRECISION
    vertical_grid%theta_init = 0.0_DEFAULT_PRECISION

    l_init_pl_theta=options_get_logical(current_state%options_database, "l_init_pl_theta")
    l_init_pl_rh=options_get_logical(current_state%options_database, "l_init_pl_rh") 
    l_init_pl_q=options_get_logical(current_state%options_database, "l_init_pl_q")
    if (l_init_pl_q) then
      allocate(names_init_pl_q(options_get_array_size(current_state%options_database, "names_init_pl_q")))
      call options_get_string_array(current_state%options_database, "names_init_pl_q", names_init_pl_q)
      do n = 1,size(names_init_pl_q)
         if (trim(names_init_pl_q(n)) .eq. 'vapour' .and. l_init_pl_rh) then 
            call log_master_log(LOG_ERROR, "Initialisation of vapour and RH - STOP")
         endif
      enddo
    end if
    l_init_pl_u=options_get_logical(current_state%options_database, "l_init_pl_u")
    l_init_pl_v=options_get_logical(current_state%options_database, "l_init_pl_v")

    l_thref=options_get_logical(current_state%options_database, "l_thref")
    l_matchthref=options_get_logical(current_state%options_database, "l_matchthref")

    if (l_thref)then
      if (.not. l_matchthref)then
        allocate(z_thref(options_get_array_size(current_state%options_database, "z_thref")), &
             f_thref(options_get_array_size(current_state%options_database, "f_thref")))
        call options_get_real_array(current_state%options_database, "z_thref", z_thref)
        call options_get_real_array(current_state%options_database, "f_thref", f_thref)
        call check_top(zztop, z_thref(size(z_thref)), 'z_thref')
        zgrid=current_state%global_grid%configuration%vertical%zn(:)
        call piecewise_linear_1d(z_thref(1:size(z_thref)), f_thref(1:size(f_thref)), zgrid, &
           current_state%global_grid%configuration%vertical%thref)
        deallocate(z_thref, f_thref)
      end if
    else
      current_state%global_grid%configuration%vertical%thref(:)=current_state%thref0
    end if

    if (l_init_pl_theta)then
      allocate(z_init_pl_theta(options_get_array_size(current_state%options_database, "z_init_pl_theta")), &
             f_init_pl_theta(options_get_array_size(current_state%options_database, "f_init_pl_theta")))
      call options_get_real_array(current_state%options_database, "z_init_pl_theta", z_init_pl_theta)
      call options_get_real_array(current_state%options_database, "f_init_pl_theta", f_init_pl_theta)
      call check_top(zztop, z_init_pl_theta(size(z_init_pl_theta)), 'z_init_pl_theta')
      call check_input_levels(size(z_init_pl_theta), size(f_init_pl_theta), "f_init_pl_theta")
      zgrid=current_state%global_grid%configuration%vertical%zn(:)
      call piecewise_linear_1d(z_init_pl_theta(1:size(z_init_pl_theta)), f_init_pl_theta(1:size(f_init_pl_theta)), zgrid, &
         current_state%global_grid%configuration%vertical%theta_init)
      if (l_matchthref) then
         if(.not. current_state%use_anelastic_equations) then
           call log_master_log(LOG_ERROR, "Non-anelastic equation set and l_maththref are incompatible")
         end if
         current_state%global_grid%configuration%vertical%thref = current_state%global_grid%configuration%vertical%theta_init
      end if
      if (.not. current_state%continuation_run) then
        do i=current_state%local_grid%local_domain_start_index(X_INDEX), current_state%local_grid%local_domain_end_index(X_INDEX)
          do j=current_state%local_grid%local_domain_start_index(Y_INDEX), current_state%local_grid%local_domain_end_index(Y_INDEX)
            current_state%th%data(:,j,i) = current_state%global_grid%configuration%vertical%theta_init(:) - &
                 current_state%global_grid%configuration%vertical%thref(:) 
          end do
        end do
      end if
      deallocate(z_init_pl_theta, f_init_pl_theta)
    end if

    if (l_init_pl_u)then
      allocate(z_init_pl_u(options_get_array_size(current_state%options_database, "z_init_pl_u")), &
             f_init_pl_u(options_get_array_size(current_state%options_database, "f_init_pl_u")))
      call options_get_real_array(current_state%options_database, "z_init_pl_u", z_init_pl_u)
      call options_get_real_array(current_state%options_database, "f_init_pl_u", f_init_pl_u)
      call check_top(zztop, z_init_pl_u(size(z_init_pl_u)), 'z_init_pl_u')
      call check_input_levels(size(z_init_pl_u), size(f_init_pl_u), "f_init_pl_u")
      zgrid=current_state%global_grid%configuration%vertical%zn(:)
      call piecewise_linear_1d(z_init_pl_u(1:size(z_init_pl_u)), f_init_pl_u(1:size(f_init_pl_u)), &
         zgrid, current_state%global_grid%configuration%vertical%u_init)
      if (.not. current_state%continuation_run) then
        do i=current_state%local_grid%local_domain_start_index(X_INDEX), current_state%local_grid%local_domain_end_index(X_INDEX)
          do j=current_state%local_grid%local_domain_start_index(Y_INDEX), current_state%local_grid%local_domain_end_index(Y_INDEX)
            current_state%u%data(:,j,i) = current_state%global_grid%configuration%vertical%u_init(:)
          end do
        end do
      end if
      deallocate(z_init_pl_u, f_init_pl_u)
    end if

    if (l_init_pl_v)then
      allocate(z_init_pl_v(options_get_array_size(current_state%options_database, "z_init_pl_v")), &
             f_init_pl_v(options_get_array_size(current_state%options_database, "f_init_pl_v")))
      call options_get_real_array(current_state%options_database, "z_init_pl_v", z_init_pl_v)
      call options_get_real_array(current_state%options_database, "f_init_pl_v", f_init_pl_v)
      call check_top(zztop, z_init_pl_v(size(z_init_pl_v)), 'z_init_pl_v')
      call check_input_levels(size(z_init_pl_v), size(f_init_pl_v), "f_init_pl_v")
      zgrid=current_state%global_grid%configuration%vertical%zn(:)
      call piecewise_linear_1d(z_init_pl_v(1:size(z_init_pl_v)), f_init_pl_v(1:size(f_init_pl_v)), &
         zgrid, current_state%global_grid%configuration%vertical%v_init)
      if (.not. current_state%continuation_run) then
        do i=current_state%local_grid%local_domain_start_index(X_INDEX), current_state%local_grid%local_domain_end_index(X_INDEX)
          do j=current_state%local_grid%local_domain_start_index(Y_INDEX), current_state%local_grid%local_domain_end_index(Y_INDEX)
            current_state%v%data(:,j,i) = current_state%global_grid%configuration%vertical%v_init(:)
          end do
        end do
      end if
      deallocate(z_init_pl_v, f_init_pl_v)
    end if

    if (l_init_pl_q)then
      nq_init=size(names_init_pl_q)
      allocate(z_init_pl_q(options_get_array_size(current_state%options_database, "z_init_pl_q")))
      call options_get_real_array(current_state%options_database, "z_init_pl_q", z_init_pl_q)
      nzq=size(z_init_pl_q)
      call check_top(zztop, z_init_pl_q(nzq), 'z_init_pl_q')
      zgrid=current_state%global_grid%configuration%vertical%zn(:)
      allocate(f_init_pl_q_tmp(nq_init*nzq))
      call options_get_real_array(current_state%options_database, "f_init_pl_q", f_init_pl_q_tmp)
      !call check_input_levels(size(z_init_pl_q), size(f_init_pl_q_tmp), "f_init_pl_q_tmp")
      allocate(f_init_pl_q(nzq, nq_init))
      f_init_pl_q(1:nzq, 1:nq_init)=reshape(f_init_pl_q_tmp, (/nzq, nq_init/))
      do n=1, nq_init
         iq=get_q_index(trim(names_init_pl_q(n)), 'piecewise_initialization')
         call check_input_levels(size(z_init_pl_q), size(f_init_pl_q(1:nzq,n)), "f_init_pl_q")
        call piecewise_linear_1d(z_init_pl_q(1:nzq), f_init_pl_q(1:nzq,n), zgrid, &
           current_state%global_grid%configuration%vertical%q_init(:,iq))
        if (.not. current_state%continuation_run) then
          do i=current_state%local_grid%local_domain_start_index(X_INDEX), &
               current_state%local_grid%local_domain_end_index(X_INDEX)
            do j=current_state%local_grid%local_domain_start_index(Y_INDEX), &
                 current_state%local_grid%local_domain_end_index(Y_INDEX)
              current_state%q(iq)%data(:,j,i) = current_state%global_grid%configuration%vertical%q_init(:, iq)
            end do
          end do
        end if
      end do
      deallocate(f_init_pl_q_tmp, z_init_pl_q, f_init_pl_q, names_init_pl_q)
   end if
    deallocate(zgrid)      
  end subroutine calculate_initial_profiles

  !> Calculates the mixing length for the neutral case
  !! @param current_state The current model state_mod
  !! @param vertical_grid The vertical grid that we are working on
  !! @param kkp Number of grid points in a vertical column
  subroutine calculate_mixing_length_for_neutral_case(current_state, vertical_grid, kkp)
    type(model_state_type), intent(inout) :: current_state
    type(vertical_grid_configuration_type), intent(inout) :: vertical_grid
    integer, intent(in) :: kkp

    integer :: k

    do k=2, kkp-1
      vertical_grid%rneutml(k)=sqrt(1.0_DEFAULT_PRECISION/(1.0_DEFAULT_PRECISION/(von_karman_constant*&
           (vertical_grid%z(k)+z0))**2+1.0_DEFAULT_PRECISION/current_state%rmlmax**2) )
      vertical_grid%rneutml_sq(k)=vertical_grid%rneutml(k)*vertical_grid%rneutml(k)
    end do    
  end subroutine calculate_mixing_length_for_neutral_case 

  !> Sets the buoyancy coefficient from the grid configuration and configuration
  !! @param current_state The current model state_mod
  !! @param vertical_grid The vertical grid that we are working on
  !! @param kkp Number of grid points in a vertical column
  subroutine set_buoyancy_coefficient(current_state, vertical_grid, kkp)
    type(model_state_type), intent(inout), target :: current_state
    type(vertical_grid_configuration_type), intent(inout) :: vertical_grid
    integer, intent(in) :: kkp

    integer :: k    

    if (.not. current_state%passive_th) then
      if(current_state%use_anelastic_equations)then                                                      
        do k=1, kkp-1          
          vertical_grid%buoy_co(k)=cp*(vertical_grid%prefn(k)**r_over_cp-vertical_grid%prefn(k+1)**r_over_cp)/&
               ((current_state%surface_reference_pressure**r_over_cp)*vertical_grid%dzn(k+1))
        end do
      else                                                                     
        vertical_grid%buoy_co(1:kkp-1)=G/current_state%thref0        ! _Boussinesq
      end if
      ! Dummy value at top level
      vertical_grid%buoy_co(kkp)=0.
    else
      vertical_grid%buoy_co(:)=0.
    end if
  end subroutine set_buoyancy_coefficient

  !> Setting up reference state liquid water temperature and saturation mixing ratio on main levels.
  !! @param current_state The current model state_mod
  !! @param vertical_grid The vertical grid that we are working on
  !! @param kkp Number of grid points in a vertical column
  subroutine setup_reference_state_liquid_water_temperature_and_saturation(current_state, vertical_grid, kkp)
    type(model_state_type), intent(inout) :: current_state
    type(vertical_grid_configuration_type), intent(inout) :: vertical_grid
    integer, intent(in) :: kkp

    real(kind=DEFAULT_PRECISION) :: delta_t=1.0_DEFAULT_PRECISION, qlinit, tinit, qsatin, dqsatdtin, dsatfacin
    integer :: iter, k
    
    do k=1,kkp
       vertical_grid%tref(k)=vertical_grid%thref(k)*(vertical_grid%prefn(k)/current_state%surface_reference_pressure)**r_over_cp
      vertical_grid%tstarpr(k)=0.0_DEFAULT_PRECISION
   end do
   if (current_state%th%active) then
      ! PREFRCP is used and hence calculated if theta is active
      do k = 1,kkp   
         vertical_grid%prefrcp(k)=(current_state%surface_reference_pressure/vertical_grid%prefn(k))**r_over_cp
         vertical_grid%rprefrcp(k)=1.0_DEFAULT_PRECISION/vertical_grid%prefrcp(k)
         ! Denotion between setup run and chain run in LEM - need to consider here too
         vertical_grid%qsat(k)=qsaturation(vertical_grid%tref(k), 0.01_DEFAULT_PRECISION*vertical_grid%prefn(k))
         vertical_grid%dqsatdt(k)=(qsaturation(vertical_grid%tref(k)+delta_t, 0.01_DEFAULT_PRECISION*vertical_grid%prefn(k)) -&
              qsaturation(vertical_grid%tref(k)-delta_t, 0.01_DEFAULT_PRECISION*vertical_grid%prefn(k)))/&
              (2.0_DEFAULT_PRECISION*delta_t)
         vertical_grid%qsatfac(k)=1.0_DEFAULT_PRECISION/(1.0_DEFAULT_PRECISION+rlvap_over_cp*vertical_grid%dqsatdt(k))
      end do
      if (current_state%calculate_th_and_q_init) then
         do k=1,kkp
            !       !Note that at this point THETA_INIT and QINIT(IQ=1) are still
            !       !theta_l and q_t, as read in from configuration.
            !       ! start from input QL profile
            qlinit=qinit(k, current_state%liquid_water_mixing_ratio_index)
            do iter=1,5
               !         ! calculate T and thence new q_l from Taylor expansion
               !         ! keeping theta_l and q_t fixed
               !         ! Note theta_l = theta - (L/c_p)*q_l here
               tinit     = vertical_grid%theta_init(k)*vertical_grid%rprefrcp(k) + rlvap_over_cp*qlinit
               qsatin    = qsaturation(tinit, 0.01_DEFAULT_PRECISION*vertical_grid%prefn(k))
               dqsatdtin = dqwsatdt(qsatin, tinit)
               dsatfacin=( 1.0_DEFAULT_PRECISION/(1.0_DEFAULT_PRECISION + rlvap_over_cp*dqsatdtin*vertical_grid%rprefrcp(k)))
               qlinit=max(0.0_DEFAULT_PRECISION, (qinit(k, current_state%water_vapour_mixing_ratio_index)-&
                    (qsatin+dqsatdtin*(vertical_grid%theta_init(k)*vertical_grid%rprefrcp(k)-tinit) ))*dsatfacin)
            end do
            qinit(k, current_state%liquid_water_mixing_ratio_index)=qlinit
            qinit(k, current_state%water_vapour_mixing_ratio_index)=qinit(k,current_state%water_vapour_mixing_ratio_index)-qlinit
            vertical_grid%theta_init(k)=vertical_grid%theta_init(k)+rlvap_over_cp*qlinit

            ! Denotion between setup run and chain run in LEM - need to consider here too
            vertical_grid%tstarpr(k)= tinit-vertical_grid%tref(k)
            vertical_grid%qsat(k)=qsatin
            vertical_grid%dqsatdt(k)=dqsatdtin        
            vertical_grid%qsatfac(k)= ( 1.0_DEFAULT_PRECISION/ ( 1.0_DEFAULT_PRECISION + rlvap_over_cp*dqsatdtin ) )
         end do
      endif
   end if
  end subroutine setup_reference_state_liquid_water_temperature_and_saturation  

  !> Sets up the reference properties for the vertical grid at each point
  !! @param current_state The current model state_mod
  !! @param vertical_grid The vertical grid that we are working on
  !! @param kkp The number of grid points in a vertical column
  subroutine set_up_vertical_reference_properties(current_state, vertical_grid, kkp)
    type(model_state_type), intent(inout) :: current_state
    type(vertical_grid_configuration_type), intent(inout) :: vertical_grid
    integer, intent(in) :: kkp
    
    integer :: k

    do k=1,kkp
!       vertical_grid%thref(k)=current_state%thref0
!       vertical_grid%theta_init(k)=vertical_grid%thref(k) ! In LEM this can also be set from configuration (TODO)       
       vertical_grid%prefn(k)=0.0_DEFAULT_PRECISION
       vertical_grid%pdiff(k)=0.0_DEFAULT_PRECISION
       vertical_grid%rho(k)=current_state%rhobous
       vertical_grid%rhon(k)=current_state%rhobous
    end do
    do k=1,kkp-1
      vertical_grid%dthref(k)=vertical_grid%thref(k+1)-vertical_grid%thref(k)
    end do
    vertical_grid%dthref(kkp)=0.0_DEFAULT_PRECISION
  end subroutine set_up_vertical_reference_properties  

  !> Sets up and smooths the vertical grid. This is based upon the grid configuration already read in
  !! @param vertical_grid The vertical grid
  !! @param kgd The grid heights per division
  !! @param hgd The real world (m) heights per division
  !! @param kkp Number of points in the vertical domain
  !! @param zztop The real world (m) height of the top of the column
  !! @param nsmth Number of smoothing iterations to run on the grid
  !! @param origional_setup To use the origional vertical grid setup routine or the new one
  subroutine set_up_and_smooth_grid(vertical_grid, kgd, hgd, ninitp, kkp, zztop, nsmth, origional_setup, continuation_run)
    type(vertical_grid_configuration_type), intent(inout) :: vertical_grid
    integer, dimension(:), intent(in) :: kgd
    real(kind=DEFAULT_PRECISION), dimension(:), intent(in) :: hgd
    integer, intent(in) :: ninitp, kkp, nsmth
    real(kind=DEFAULT_PRECISION),intent(in) :: zztop
    logical, intent(in) :: origional_setup, continuation_run

    integer :: k

    if (.not. continuation_run) then
      if (origional_setup) then
        call original_vertical_grid_setup(vertical_grid, kgd, hgd, ninitp, kkp, zztop, nsmth)
      else
        call new_vertical_grid_setup(vertical_grid, kgd, kkp, zztop)
      end if
    end if
    
    ! Regardless of the vertical grid computation method, set the level deltas
    do k=2,kkp
       vertical_grid%dz(k)=vertical_grid%z(k)-vertical_grid%z(k-1)
       vertical_grid%dzn(k)= vertical_grid%zn(k)-vertical_grid%zn(k-1)                                                     
       vertical_grid%rdz(k)=1./vertical_grid%dz(k)                                                          
       vertical_grid%rdzn(k)=1./vertical_grid%dzn(k)                                                        
    end do
    vertical_grid%dzn(1)=0.d0
  end subroutine set_up_and_smooth_grid

  !> The original vertical grid setup and smoothing as per the LEM
  !! @param vertical_grid The vertical grid
  !! @param kgd The grid heights per division
  !! @param hgd The real world (m) heights per division
  !! @param kkp Number of points in the vertical domain
  !! @param zztop The real world (m) height of the top of the column
  !! @param nsmth Number of smoothing iterations to run on the grid
  subroutine original_vertical_grid_setup(vertical_grid, kgd, hgd, ninitp, kkp, zztop, nsmth)
    type(vertical_grid_configuration_type), intent(inout) :: vertical_grid
    integer, dimension(:), intent(in) :: kgd
    real(kind=DEFAULT_PRECISION), dimension(:), intent(in) :: hgd
    integer, intent(in) :: ninitp, kkp, nsmth
    real(kind=DEFAULT_PRECISION), intent(in) :: zztop

    integer :: n, k
    
    call create_linear_grid(vertical_grid, kgd, hgd, ninitp, kkp, zztop)
    ! Smooth grid
    vertical_grid%z(1)=0.0_DEFAULT_PRECISION
    vertical_grid%z(kkp)=zztop
    do n=1,nsmth
      do k=2,kkp
        vertical_grid%zn(k)=0.5_DEFAULT_PRECISION*(vertical_grid%z(k)+vertical_grid%z(k-1))
      end do
      do k=2,kkp-1
        vertical_grid%z(k)=0.5_DEFAULT_PRECISION*(vertical_grid%zn(k)+vertical_grid%zn(k+1))
      end do
    end do
    ! Fourth order interpolation
    do k=3,kkp-1
      vertical_grid%zn(k)=0.0625_DEFAULT_PRECISION*(9.0_DEFAULT_PRECISION*&
           (vertical_grid%z(k-1)+vertical_grid%z(k))-vertical_grid%z(k+1)-vertical_grid%z(k-2))
    end do
    vertical_grid%zn(2)=0.5_DEFAULT_PRECISION*(vertical_grid%z(1)+vertical_grid%z(2))
    vertical_grid%zn(1)=-vertical_grid%zn(2)
    vertical_grid%zn(kkp)=0.5_DEFAULT_PRECISION*(vertical_grid%z(kkp-1)+vertical_grid%z(kkp))
  end subroutine original_vertical_grid_setup

  !> The newer vertical grid setup routine
  !! @param vertical_grid The vertical grid
  !! @param kgd The grid heights per division
  !! @param kkp Number of points in the vertical domain
  !! @param zztop The real world (m) height of the top of the column
  subroutine new_vertical_grid_setup(vertical_grid, kgd, kkp, zztop)
    type(vertical_grid_configuration_type), intent(inout) :: vertical_grid
    integer, dimension(:), intent(in) :: kgd
    integer, intent(in) :: kkp
    real(kind=DEFAULT_PRECISION), intent(in) :: zztop

    real(kind=DEFAULT_PRECISION) :: a(2*kkp), r1, d1, dd, d0, a0
    logical :: first_gt=.true.
    integer :: k, k0

    r1=1.10_DEFAULT_PRECISION
    d1=10.0_DEFAULT_PRECISION

    dd=0.5_DEFAULT_PRECISION*d1
    d0=dd
    a(1)=-dd
    a(2)=0.0_DEFAULT_PRECISION
    do k=3, kkp*2
      if (d0 .gt. dd .or. k==3) then
        a(k)=a(k-1)+dd
        if (.not. (dd .gt. 25.0_DEFAULT_PRECISION .and. a(k) .lt. 2000.0_DEFAULT_PRECISION)) then
          if (a(k) .lt. 2000.0_DEFAULT_PRECISION) then
            dd=dd*r1
          else
            dd=dd*(1.0_DEFAULT_PRECISION+(r1-1.0_DEFAULT_PRECISION)/1.5_DEFAULT_PRECISION)
          end if
        end if
        d0=(zztop-a(k))/real(kkp*2-k, kind=DEFAULT_PRECISION)
      else
        if (first_gt) then
          k0=k
          a0=a(k-1)+d0
          first_gt=.false.
        end if
        a(k)=a0+real(k-k0, kind=DEFAULT_PRECISION)*d0
      end if
    end do    
    
    do k=1, kkp
      vertical_grid%z(k)=a(k*2)
      vertical_grid%zn(k)=a(2*k-1)
    end do
  end subroutine new_vertical_grid_setup

  !> Creates the linear vertical grid based upon the configuration properties.
  !!
  !! This will correspond the vertical k grid points to the configuration which might be only
  !! a few reference properties.
  !! @param vertical_grid The vertical grid
  !! @param kgd The grid heights per division
  !! @param hgd The real world (m) heights per division
  !! @param kkp Number of points in the vertical domain
  !! @param zztop The real world (m) height of the top of the column
  !! @param nsmth Number of smoothing iterations to run on the grid
  subroutine create_linear_grid(vertical_grid, kgd, hgd, ninitp, kkp, zztop)
    type(vertical_grid_configuration_type), intent(inout) :: vertical_grid
    integer, dimension(:), intent(in) :: kgd
    real(kind=DEFAULT_PRECISION), dimension(:), intent(in) :: hgd
    integer, intent(in) :: ninitp, kkp
    real(kind=DEFAULT_PRECISION), intent(in) :: zztop

    integer :: kmax, k, i
    real(kind=DEFAULT_PRECISION) :: zmax

    vertical_grid%z(1) = 0.0_DEFAULT_PRECISION
    kmax=kgd(1)
    zmax=0.0_DEFAULT_PRECISION
    if (kgd(1) .gt. 1) then
      do k=1,kgd(1)
        ! Loop up to first division point
        vertical_grid%z(k)=real(k-1, kind=DEFAULT_PRECISION)*hgd(1)/real(kgd(1)-1, kind=DEFAULT_PRECISION)
      end do
      do i=2,ninitp
        if(kgd(i) .gt. 0) then
          kmax=kgd(i)
          zmax=hgd(i)
          
          do k=kgd(i-1)+1,kgd(i)
            vertical_grid%z(k)=hgd(i-1)+(hgd(i)-hgd(i-1))*real(k-kgd(i-1), kind=DEFAULT_PRECISION)&
                 /real(kgd(i)-kgd(i-1), kind=DEFAULT_PRECISION)
          end do
        end if
      end do
   end if
    if(kmax .lt. kkp)then
      do k=kmax,kkp
        ! Handle any points above the kth max division
        vertical_grid%z(k)=zmax+(zztop-zmax)*real(k-kmax, kind=DEFAULT_PRECISION)/real(kkp-kmax, kind=DEFAULT_PRECISION)
     end do
    end if
  end subroutine create_linear_grid

  !> Allocates the data required for the vertical grid configuration
  !! @param vertical_grid The vertical grid that we are working with
  !! @param n The number of grid points in the vertical
  !! @param nq The number of q fields
  subroutine allocate_vertical_grid_data(vertical_grid, n, nq)
    type(vertical_grid_configuration_type), intent(inout) :: vertical_grid
    integer, intent(in) :: n
    integer, intent(in) :: nq

    allocate(vertical_grid%dz(n), vertical_grid%dzn(n),&
         vertical_grid%czb(n), vertical_grid%cza(n), vertical_grid%czg(n), vertical_grid%czh(n),&
         vertical_grid%rdz(n), vertical_grid%rdzn(n), vertical_grid%tzc1(n), vertical_grid%tzc2(n),&
         vertical_grid%tzd1(n), vertical_grid%tzd2(n), vertical_grid%theta_init(n), vertical_grid%temp_init(n), &
         vertical_grid%rh_init(n), &
         vertical_grid%tref(n), vertical_grid%prefn(n), vertical_grid%pdiff(n), vertical_grid%prefrcp(n), &
         vertical_grid%rprefrcp(n), vertical_grid%rho(n), vertical_grid%rhon(n), vertical_grid%tstarpr(n), &
         vertical_grid%qsat(n), vertical_grid%dqsatdt(n), vertical_grid%qsatfac(n), vertical_grid%dthref(n), &
         vertical_grid%rneutml(n), vertical_grid%rneutml_sq(n), vertical_grid%buoy_co(n), &
         vertical_grid%u_init(n), vertical_grid%v_init(n), vertical_grid%theta_rand(n), vertical_grid%w_rand(n), &
         vertical_grid%w_subs(n), vertical_grid%u_force(n), vertical_grid%v_force(n), vertical_grid%theta_force(n))

    if (.not. allocated(vertical_grid%thref)) allocate(vertical_grid%thref(n))
    if (.not. allocated(vertical_grid%z)) allocate(vertical_grid%z(n))
    if (.not. allocated(vertical_grid%zn)) allocate(vertical_grid%zn(n))

    allocate(vertical_grid%q_rand(n,nq), vertical_grid%q_init(n,nq), vertical_grid%q_force(n,nq))
  end subroutine allocate_vertical_grid_data  

  !> Initialises the horizontal grid configurations
  !! @param current_state The current model state_mod
  subroutine initialise_horizontalgrid_configuration_types(current_state)
    type(model_state_type), intent(inout) :: current_state

    current_state%global_grid%configuration%horizontal%dx = merge(real(current_state%global_grid%resolution(X_INDEX)), &
         DEFAULT_SPACING, current_state%global_grid%active(X_INDEX))
    current_state%global_grid%configuration%horizontal%dy = merge(real(current_state%global_grid%resolution(Y_INDEX)), &
         DEFAULT_SPACING, current_state%global_grid%active(Y_INDEX))

    current_state%global_grid%configuration%horizontal%cx=1./current_state%global_grid%configuration%horizontal%dx
    current_state%global_grid%configuration%horizontal%cy=1./current_state%global_grid%configuration%horizontal%dy
    current_state%global_grid%configuration%horizontal%cx2=current_state%global_grid%configuration%horizontal%cx ** 2
    current_state%global_grid%configuration%horizontal%cy2=current_state%global_grid%configuration%horizontal%cy ** 2
    current_state%global_grid%configuration%horizontal%cxy=current_state%global_grid%configuration%horizontal%cx * &
         current_state%global_grid%configuration%horizontal%cy
    current_state%global_grid%configuration%horizontal%tcx=&
         0.25_DEFAULT_PRECISION/current_state%global_grid%configuration%horizontal%dx
    current_state%global_grid%configuration%horizontal%tcy=&
         0.25_DEFAULT_PRECISION/current_state%global_grid%configuration%horizontal%dy
  end subroutine initialise_horizontalgrid_configuration_types  

  !> Set reference profile of potential temperature for the Boussinesq/Anelastic approximation
  !! Note that this is not in general the same as the profile defining the initial vertical distribution of potential
  !! temperature for the integration. In particular, while the later may contain sharp changes in gradient representing
  !! a capping inversion or the tropopause, for example, the reference profile should be smooth
  !! @param current_state The current model state
  subroutine set_anelastic_pressure(current_state)
    type(model_state_type), intent(inout) :: current_state
    
    if (current_state%use_anelastic_equations) then
      call compute_anelastic_pressure_profile_and_density(current_state)
    else
      if (current_state%passive_th) then
        current_state%global_grid%configuration%vertical%prefn=0.0_DEFAULT_PRECISION
      else
        call compute_anelastic_pressure_profile_only(current_state)
      end if
        current_state%global_grid%configuration%vertical%rho=current_state%rhobous
        current_state%global_grid%configuration%vertical%rhon=current_state%rhobous
        current_state%global_grid%configuration%vertical%pdiff=0.0_DEFAULT_PRECISION
    end if
  end subroutine set_anelastic_pressure  

  !> Computes the anelastic pressure only - if we are using Boussinesq approximation
  !! @param current_state The current model state
  subroutine compute_anelastic_pressure_profile_only(current_state)
    type(model_state_type), intent(inout) :: current_state

    integer :: ipass, k
    real(kind=DEFAULT_PRECISION) ::    p0    &!pressure at z=0 adjustments made after 1st iteration so P0=PSF after 2nd iteration
        ,   ptop  &!pressure at z=ZN(KKP)
        , thprof(current_state%local_grid%size(Z_INDEX))

    
    ! TODO: NOTE - we are mocking in thprof at the moment, this should be read from a configuration and used here instead
    thprof=0.0_DEFAULT_PRECISION
    ptop=0.0_DEFAULT_PRECISION
    current_state%global_grid%configuration%vertical%pdiff(current_state%local_grid%size(Z_INDEX))=0.0_DEFAULT_PRECISION

    do ipass=1,2 ! _after first pass adjust PTOP
      current_state%global_grid%configuration%vertical%prefn(current_state%local_grid%size(Z_INDEX))=&
           (ptop/current_state%surface_reference_pressure)**r_over_cp
      do k=current_state%local_grid%size(Z_INDEX)-1,1,-1
        current_state%global_grid%configuration%vertical%pdiff(k)=G*&
                     current_state%global_grid%configuration%vertical%dzn(k+1)/(0.5_DEFAULT_PRECISION*cp*&
                     (current_state%global_grid%configuration%vertical%thref(k)+&
                     current_state%global_grid%configuration%vertical%thref(k+1)))
      end do
      do k=current_state%local_grid%size(Z_INDEX)-1,1,-1
        current_state%global_grid%configuration%vertical%prefn(k)=&
             current_state%global_grid%configuration%vertical%prefn(k+1)+&
             current_state%global_grid%configuration%vertical%pdiff(k)
      end do
      do k=current_state%local_grid%size(Z_INDEX),1,-1
        current_state%global_grid%configuration%vertical%prefn(k)=current_state%surface_reference_pressure*&
             current_state%global_grid%configuration%vertical%prefn(k)**(1.0_DEFAULT_PRECISION/r_over_cp)
      end do
      p0=0.5_DEFAULT_PRECISION*(current_state%global_grid%configuration%vertical%prefn(1)+&
             current_state%global_grid%configuration%vertical%prefn(2))
      if (ipass .eq. 1) then                                                       
        ptop=current_state%surface_pressure**r_over_cp+ptop**r_over_cp-p0**r_over_cp
        if (ptop .le. 0.0_DEFAULT_PRECISION .and. current_state%parallel%my_rank==0) then
          call log_log(LOG_ERROR, "Negative ptop in setup of anelastic. Need a warmer THREF or different setup options")
        end if
        ptop=ptop**(1.0_DEFAULT_PRECISION/r_over_cp)
      end if
    end do
  end subroutine compute_anelastic_pressure_profile_only  

  !> Computes the anelastic pressure and density - if we are using Anelastic approximation
  !! @param current_state The current model state
  subroutine compute_anelastic_pressure_profile_and_density(current_state)
    type(model_state_type), intent(inout) :: current_state

    integer :: ipass, k
    real(kind=DEFAULT_PRECISION) ::    p0    &!pressure at z=0 adjustments made after 1st iteration so P0=PSF after 2nd iteration
        ,   ptop  &!pressure at z=ZN(KKP)
        ,   thfactor !factor for multiplying TH profile (if IADJANELP=2)

    ptop=0.0_DEFAULT_PRECISION
    current_state%global_grid%configuration%vertical%pdiff(current_state%local_grid%size(Z_INDEX))=0.0_DEFAULT_PRECISION
    do ipass=1,2 ! _after first pass, may adjust basic states
        if (ipass .eq. 1 .or. ANELASTIC_PROFILE_MODE .gt. 1) then                                     
            current_state%global_grid%configuration%vertical%prefn(current_state%local_grid%size(Z_INDEX))=&
                 (ptop/current_state%surface_reference_pressure)**r_over_cp
            do k=current_state%local_grid%size(Z_INDEX)-1,1,-1
                current_state%global_grid%configuration%vertical%pdiff(k)=G*&
                     current_state%global_grid%configuration%vertical%dzn(k+1)/(0.5_DEFAULT_PRECISION*cp*&
                     (current_state%global_grid%configuration%vertical%thref(k)+&
                     current_state%global_grid%configuration%vertical%thref(k+1)))                
            end do
            do k=current_state%local_grid%size(Z_INDEX)-1,1,-1
                current_state%global_grid%configuration%vertical%prefn(k)=&
                     current_state%global_grid%configuration%vertical%prefn(k+1)+&
                     current_state%global_grid%configuration%vertical%pdiff(k)
            end do
            do k=current_state%local_grid%size(Z_INDEX),1,-1
                current_state%global_grid%configuration%vertical%prefn(k)=current_state%surface_reference_pressure*&
                     current_state%global_grid%configuration%vertical%prefn(k)**(1.0_DEFAULT_PRECISION/r_over_cp)
            end do           
        end if                                                                    
        p0=0.5_DEFAULT_PRECISION*(current_state%global_grid%configuration%vertical%prefn(1)+&
             current_state%global_grid%configuration%vertical%prefn(2))
        !       !-------------------------------------------------------------
        !       ! _If IADJANELP>1 we adjust the basic states to ensure P0=PSF,
        !       !                                as follows:
        !       !    IADJANELP=2     adjust THREF profile by constant factor
        !       !    IADJANELP=3     adjust PSF
        !       !    IADJANELP=4     adjust PTOP
        !       ! _Option 3 tends to give rather large changes in PSF, so
        !       !   I prefer 2 or 4 for most purposes
        !       !-------------------------------------------------------------
        if (ipass .eq. 1 .and. ANELASTIC_PROFILE_MODE .eq. 2) then                                    
            !             ! _adjust THREF profile by constant factor to enforce
            !             !    P0 = (fixed) PSF
            thfactor=((p0/current_state%surface_reference_pressure)**r_over_cp-(ptop/current_state%surface_reference_pressure)**&
                 r_over_cp)/((current_state%surface_pressure/current_state%surface_reference_pressure)**r_over_cp-&
                (ptop/current_state%surface_reference_pressure)**r_over_cp)
            do k=1,current_state%local_grid%size(Z_INDEX)
                current_state%global_grid%configuration%vertical%thref(k)=&
                     current_state%global_grid%configuration%vertical%thref(k)*thfactor
            end do
        end if
        if (ipass .eq. 1 .and. ANELASTIC_PROFILE_MODE .eq. 4) then                                    
            !             ! _adjust PTOP so that P0 = (fixed) PSF
            ptop=current_state%surface_pressure**r_over_cp+ptop**r_over_cp-p0**r_over_cp
            if (ptop .le. 0.0_DEFAULT_PRECISION .and. current_state%parallel%my_rank==0) then
                call log_log(LOG_ERROR, "Negative ptop in setup of anelastic. Need a warmer THREF or different setup options")
            end if
            ptop=ptop**(1.0_DEFAULT_PRECISION/r_over_cp)
        end if                                                                    
    end do
    !     !---------------------------------------
    !     ! _Finally compute density from pressure
    !     !---------------------------------------
    do k=1,current_state%local_grid%size(Z_INDEX)
        current_state%global_grid%configuration%vertical%rhon(k)=current_state%global_grid%configuration%vertical%prefn(k)&
            /(r*current_state%global_grid%configuration%vertical%thref(k)*&
            (current_state%global_grid%configuration%vertical%prefn(k)/current_state%surface_reference_pressure)**r_over_cp)
    end do
    do k=1,current_state%local_grid%size(Z_INDEX)-1
        current_state%global_grid%configuration%vertical%rho(k)=sqrt(current_state%global_grid%configuration%vertical%rhon(k)*&
             current_state%global_grid%configuration%vertical%rhon(k+1))                                           
    end do
    current_state%global_grid%configuration%vertical%rho(current_state%local_grid%size(Z_INDEX))=&
         current_state%global_grid%configuration%vertical%rhon(current_state%local_grid%size(Z_INDEX))**2/&
         current_state%global_grid%configuration%vertical%rhon(current_state%local_grid%size(Z_INDEX)-1)    
  end subroutine compute_anelastic_pressure_profile_and_density  


  subroutine check_top(zztop, z, info)
    real(kind=DEFAULT_PRECISION), intent(in) :: zztop
    real(kind=DEFAULT_PRECISION), intent(in) :: z
    character(*), intent(in) :: info

    if (z<zztop)then
      call log_master_log(LOG_ERROR, "Top of input profile is below the top of the domain:"//trim(info))
    end if

  end subroutine check_top

  subroutine check_input_levels(z_levels, field_levels, field)
    integer, intent(in) :: z_levels
    integer, intent(in) :: field_levels
    character(*), intent(in) :: field

    if (z_levels /= field_levels)then
       call log_master_log(LOG_ERROR, "Input levels not equal for "//trim(field)//", z_levels = "// &
            trim(conv_to_string(z_levels))//" field_levels = "//conv_to_string(field_levels))
    end if

  end subroutine check_input_levels
  
  subroutine set_qv_init_from_rh(current_state)

    type(model_state_type), intent(inout) :: current_state

    logical :: l_init_pl_rh    ! if .true. then initialize relative humidity field
    real(kind=DEFAULT_PRECISION) :: zztop ! top of the domain
    real(kind=DEFAULT_PRECISION) :: qsat
    real(kind=DEFAULT_PRECISION), dimension(:), allocatable :: f_init_pl_rh     ! Initial node values for relative humidity variable
    real(kind=DEFAULT_PRECISION), dimension(:), allocatable :: z_init_pl_rh     ! Initial node height values for relative humidity variable
    real(kind=DEFAULT_PRECISION), allocatable :: zgrid(:)  ! z grid to use in interpolation
    real(kind=DEFAULT_PRECISION), allocatable :: TdegK(:)  ! temperature in Kelvin
    integer :: i,j,n, k ! loop counters
    integer :: iq  ! temporary q varible index

    type(vertical_grid_configuration_type) :: vertical_grid

    vertical_grid=current_state%global_grid%configuration%vertical

    allocate(zgrid(current_state%local_grid%local_domain_end_index(Z_INDEX)))
    
    zztop = current_state%global_grid%top(Z_INDEX)

    l_init_pl_rh=options_get_logical(current_state%options_database, "l_init_pl_rh") 

    if (l_init_pl_rh)then
       allocate(z_init_pl_rh(options_get_array_size(current_state%options_database, "z_init_pl_rh")), &
             f_init_pl_rh(options_get_array_size(current_state%options_database, "f_init_pl_rh")))
       call options_get_real_array(current_state%options_database, "z_init_pl_rh", z_init_pl_rh)
       call options_get_real_array(current_state%options_database, "f_init_pl_rh", f_init_pl_rh)
       call check_top(zztop, z_init_pl_rh(size(z_init_pl_rh)), 'z_init_pl_rh')
       call check_input_levels(size(z_init_pl_rh), size(f_init_pl_rh), "f_init_pl_rh")
       zgrid=current_state%global_grid%configuration%vertical%zn(:)
      call piecewise_linear_1d(z_init_pl_rh(1:size(z_init_pl_rh)), f_init_pl_rh(1:size(f_init_pl_rh)), zgrid, &
           current_state%global_grid%configuration%vertical%rh_init)
      
      if (.not. current_state%passive_q .and. current_state%th%active) then
         iq=get_q_index('vapour', 'piecewise_initialization')
         allocate(TdegK(current_state%local_grid%local_domain_end_index(Z_INDEX)))
         TdegK(:) = current_state%global_grid%configuration%vertical%theta_init(:)* &
              (vertical_grid%prefn(:)/current_state%surface_reference_pressure)**r_over_cp
         do k = current_state%local_grid%local_domain_start_index(Z_INDEX), &
              current_state%local_grid%local_domain_end_index(Z_INDEX)
            qsat=qsaturation(TdegK(k), current_state%global_grid%configuration%vertical%prefn(k)/100.)    
            current_state%global_grid%configuration%vertical%q_init(k, iq) = & 
                 (current_state%global_grid%configuration%vertical%rh_init(k)/100.0)*qsat
            !print *,  current_state%global_grid%configuration%vertical%rh_init(k), &
            !     current_state%global_grid%configuration%vertical%q_init(k, iq), &
            !     TdegK(k)
         enddo
         if (.not. current_state%continuation_run) then
            do i=current_state%local_grid%local_domain_start_index(X_INDEX), &
                 current_state%local_grid%local_domain_end_index(X_INDEX)
               do j=current_state%local_grid%local_domain_start_index(Y_INDEX), &
                    current_state%local_grid%local_domain_end_index(Y_INDEX)
                  current_state%q(iq)%data(:,j,i) = current_state%global_grid%configuration%vertical%q_init(:, iq)
               end do
            end do
         end if

         deallocate(TdegK)
      else
         call log_master_log(LOG_ERROR, "Initialising with RH but q and/or theta passive")
      end if
      
      deallocate(z_init_pl_rh, f_init_pl_rh)
    end if

  end subroutine set_qv_init_from_rh

end module gridmanager_mod
