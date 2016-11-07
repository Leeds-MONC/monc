module kidtestcase_mod
  use monc_component_mod, only : component_descriptor_type
  use grids_mod, only : global_grid_type, X_INDEX, Y_INDEX, Z_INDEX
  use state_mod, only : model_state_type
  use logging_mod, only : LOG_ERROR, log_log, log_master_log
  use datadefn_mod, only : DEFAULT_PRECISION
  use science_constants_mod, only : pi, r_over_cp, G, Ru=>r
  use q_indices_mod, only: get_q_index, standard_q_names
  use interpolation_mod, only: piecewise_linear_1d
  use saturation_mod, only: qsaturation
  use optionsdatabase_mod, only :  options_get_integer
  implicit none

#ifndef TEST_MODE
  private
#endif

  public kidtestcase_get_descriptor


  !!! NOTE THE RENORMALIZATIONS BELOW WON@T WORK WITH MULTIPLE PROCESSES

  integer :: case_number

  integer, parameter :: CASE_CU = 1
  integer, parameter :: CASE_EGG = 2
  integer, parameter :: CASE_SC = 3
  integer, parameter :: CASE_SQUALL = 4
  integer, parameter :: CASE_HILL_FLOW = 5

contains

  !> Provides the component descriptor for the core to register
  !! @returns The descriptor describing this component
  type(component_descriptor_type) function kidtestcase_get_descriptor()
    kidtestcase_get_descriptor%name="kid_testcase"
    kidtestcase_get_descriptor%version=0.1
    
    kidtestcase_get_descriptor%initialisation=>initialise_callback
    kidtestcase_get_descriptor%timestep=>timestep_callback
  end function kidtestcase_get_descriptor

  subroutine initialise_callback(current_state)
    type(model_state_type), target, intent(inout) :: current_state

    case_number=options_get_integer(current_state%options_database, "case_number") 

    if (current_state%continuation_run) return

    select case(case_number)
    case(CASE_CU)
      call set_Cu_profiles(current_state)     
    case (CASE_EGG)
      call set_2D_Eggs(current_state)
    case (CASE_SC)
      call set_2D_Sc(current_state)
    case (CASE_SQUALL)
      call set_2D_squall(current_state)
    case (CASE_HILL_FLOW)
      call set_2D_hills(current_state)
    end select

  end subroutine initialise_callback

  subroutine set_Cu_profiles(current_state)
    type(model_state_type), target, intent(inout) :: current_state

    ! Set up the initial thermodynamic profile
    
    ! Levels for interpolation
    real(DEFAULT_PRECISION), allocatable :: &
       pHeight(:)         & ! height
       ,pTheta(:)          & ! theta
       ,pqv(:)             & ! qv
       ,pRH(:)               ! RH

    ! local allocatable arrays for temperature and presssure
    real(DEFAULT_PRECISION), allocatable :: &
         press_cu(:)  & ! pressure for Cu case
         ,temp_cu(:)  & ! temperature for Cu case
         ,rh_cu(:)     ! RH for Cu case

    real(DEFAULT_PRECISION) :: tempk, tempkm, delz, delt, tavi

    real(kind=DEFAULT_PRECISION), allocatable :: zgrid(:)  ! z grid to use in interpolation

    integer :: nlevs, km1, iq, i, j, k

    logical :: l_cu_cold=.false.


    real(kind=DEFAULT_PRECISION) :: zztop ! top of the domain
    
    zztop = current_state%global_grid%top(Z_INDEX)

    nlevs = 26

    allocate(zgrid(current_state%local_grid%local_domain_end_index(Z_INDEX)))

    allocate(pHeight(nlevs))
    allocate(pTheta(nlevs))
    allocate(pqv(nlevs))

    allocate(press_cu(nlevs))
    allocate(temp_cu(nlevs))
    allocate(rh_cu(nlevs))

    
    ! pqv in g/kg
    pqv=(/14.5,  14.5,  14.5,  14.0,  13.7,  13.9,  13.9,   &
         10.3,  10.3,  10.0,   9.9,   8.9,   7.9,   4.0,   2.3, &
         1.2,   1.2,   0.9,   0.6,   2.0,   1.6,   0.4,   1.5, &
         0.9,   0.5,   0.4/)

    press_cu=(/1014., 1010., 1000.,  990.,  975.,  960.,  950., &
          925.,  900.,  875.,  850.,  825.,  800.,  790.,  775., & 
          765.,  755.,  745.,  730.,  715.,  700.,  650.,  600., &
          550.,  500.,  450./)

    temp_cu=(/298.4, 298.0, 296.8, 295.7, 295.0, 293.7, 293.1, &
         291.4, 290.0, 288.0, 286.5, 285.1, 284.2, 284.5, 284.1, &
         284.4, 283.4, 284.2, 284.4, 283.2, 282.0, 280.0, 275.7, &
         272.5, 269.5, 266./)

    ! set qv kg/kg
    pqv(:) = pqv(:)*1.e-3

    if (l_cu_cold)then
      do k=1,nlevs
        rh_cu(k) =  pqv(k)/qsaturation(temp_cu(k),press_cu(k))
      end do
      
      temp_cu=temp_cu-20 ! reduce temperature for testing ice
      
      do k=1,nlevs
        pqv(k) =  rh_cu(k)*qsaturation(temp_cu(k),press_cu(k))
      end do
    end if

    
    ! calculate theta from temp_cu
    do k = 1, nlevs
       ptheta(k) = temp_cu(k)*(1.e3/press_cu(k))**r_over_cp
    enddo
    
    
    ! calculate approximate height from pressure 
    pheight(1) = 0.0
    do k = 2, nlevs
       km1 = k-1
       tempk = ptheta(k) * (1.e3/press_cu(k))**(-r_over_cp) &
            * (1. + .6*pqv(k))
       tempkm = ptheta(km1) * (1.e3/press_cu(km1))**(-r_over_cp) &
            * (1. + .6*pqv(km1))
       
       delt=tempk-tempkm
       if(delt.gt.1.e-4) then
          tavi=log(tempk/tempkm)/delt
       else
          tavi=1./tempk
       endif
       
       delz=-ru/(tavi*g)*log(press_cu(k)/press_cu(km1))
       pheight(k) = pheight(km1) + delz
    enddo
    
    call check_top(zztop, pheight(nlevs), 'kid_case:setup cu')

    zgrid=current_state%global_grid%configuration%vertical%zn(:)
    call piecewise_linear_1d(pheight, ptheta, zgrid, &
       current_state%global_grid%configuration%vertical%thref)

    current_state%global_grid%configuration%vertical%theta_init = current_state%global_grid%configuration%vertical%thref
    current_state%th%data=0.0

    iq=get_q_index(standard_q_names%VAPOUR, 'kidtestcase-Cu')
    
    call piecewise_linear_1d(pheight, pqv, zgrid, &
       current_state%global_grid%configuration%vertical%q_init(:,iq))    

    do i=current_state%local_grid%local_domain_start_index(X_INDEX), current_state%local_grid%local_domain_end_index(X_INDEX)
      do j=current_state%local_grid%local_domain_start_index(Y_INDEX), current_state%local_grid%local_domain_end_index(Y_INDEX)
        current_state%q(iq)%data(:,j,i) = current_state%global_grid%configuration%vertical%q_init(:, iq)
      end do
    end do

    deallocate(pqv)
    deallocate(rh_cu)
    deallocate(pTheta)
    deallocate(pHeight)

    deallocate(press_cu)
    deallocate(temp_cu)

  end subroutine set_Cu_profiles

  subroutine timestep_callback(current_state)
    type(model_state_type), target, intent(inout) :: current_state

    select case(case_number)
    case(CASE_CU)
      call set_2D_Cu_wind_field(current_state)
    end select

  end subroutine timestep_callback

  subroutine set_2D_Eggs(current_state)
    type(model_state_type), intent(inout), target :: current_state
    
    !
    ! Set up Chris' Egg shells
    !
    real(DEFAULT_PRECISION), dimension(:), allocatable :: dzp, dxp

    real(DEFAULT_PRECISION):: maxW, eta_1, xi_1, eta_2, xi_2, zeta_1, zeta_2
    real(DEFAULT_PRECISION):: x, z

    real(DEFAULT_PRECISION):: x1_1,x2_1,z1_1,z2_1, x1_2,x2_2,z1_2,z2_2

    real(DEFAULT_PRECISION), dimension(:,:), allocatable :: phi

    real(DEFAULT_PRECISION) :: renorm

    integer :: i,j,k,iq

    maxW=10.

    allocate(phi(0:current_state%local_grid%local_domain_end_index(X_INDEX), &
       1:current_state%local_grid%local_domain_end_index(Z_INDEX)),  &
       dzp(current_state%local_grid%local_domain_end_index(Z_INDEX)), &
       dxp(current_state%local_grid%local_domain_end_index(X_INDEX)))
   
    dzp=current_state%global_grid%configuration%vertical%dz    
    dxp = current_state%global_grid%configuration%horizontal%dx
    dzp(1)=dzp(2)

    x1_1=20.e3
    x2_1=110.e3
    z1_1=0.
    z2_1=10.e3

    x1_2=110.e3
    x2_2=200.e3
    z1_2=0.
    z2_2=10.e3


    do i=current_state%local_grid%local_domain_start_index(X_INDEX)-1, current_state%local_grid%local_domain_end_index(X_INDEX)
      x = (0.5 + current_state%local_grid%start(X_INDEX) + (i-current_state%local_grid%local_domain_start_index(X_INDEX))) * &
         current_state%global_grid%configuration%horizontal%dx
      do k=1, current_state%local_grid%local_domain_end_index(Z_INDEX) 
        z = current_state%global_grid%configuration%vertical%z(k)
        xi_1 = min(max(1.- (x-x1_1)/(x2_1-x1_1),0.0),1.0)
        eta_1=min(max(2*(z-z1_1)/(z2_1-z1_1) - 1.0,-1.0),1.0)
        xi_2 = min(max((x-x1_2)/(x2_2-x1_2),0.0),1.0)
        eta_2=min(max(2*(z-z1_2)/(z2_2-z1_2) - 1.0,-1.0),1.0)
        
        zeta_1 = 4.*(xi_1 - sqrt(xi_1)) + 1.0 + eta_1*eta_1
        zeta_2 = 4.*(xi_2 - sqrt(xi_2)) + 1.0 + eta_2*eta_2
        
        phi(i,k) = (1.0 - min(zeta_2, 1.0)) - (1.0 - min(zeta_1, 1.0)) 

      end do
    end do
    
    !  calculate rho*vel by derivation of streamfunction and normalize
    !  rho*ux velocity:
    do i=current_state%local_grid%local_domain_start_index(X_INDEX), current_state%local_grid%local_domain_end_index(X_INDEX)
      do k=2, current_state%local_grid%size(Z_INDEX)
        current_state%zu%data(k,:,i)=-(phi(i,k)-phi(i,k-1))/dzp(k) &
             /current_state%global_grid%configuration%vertical%rhon(k)
        current_state%zw%data(k,:,i)=(phi(i,k)-phi(i-1,k))/dxp(i)  &
             /current_state%global_grid%configuration%vertical%rho(k)
      enddo
      current_state%zu%data(1,:,i)=-current_state%zu%data(1,:,i)
      current_state%zw%data(1,:,i)=-current_state%zw%data(1,:,i)
    enddo


    ! Renormalize winds
    renorm = maxW/maxval(current_state%zw%data)
    
    current_state%zw%data = renorm*current_state%zw%data
    current_state%zu%data = renorm*current_state%zu%data

    current_state%u%data = current_state%zu%data
    current_state%w%data = current_state%zw%data

    ! Set some tracers
    
    iq=get_q_index('tracer', 'kidtestcase-Cu')
    do i=current_state%local_grid%local_domain_start_index(X_INDEX)-1, current_state%local_grid%local_domain_end_index(X_INDEX)
      x = (0.5 + current_state%local_grid%start(X_INDEX) + (i-current_state%local_grid%local_domain_start_index(X_INDEX))) * &
         current_state%global_grid%configuration%horizontal%dx
      if (abs(x - 0.5*(x2_1+x2_1)) < current_state%global_grid%configuration%horizontal%dx)then
        do k=1, current_state%local_grid%local_domain_end_index(Z_INDEX) 
        z = current_state%global_grid%configuration%vertical%z(k)
          current_state%q(iq)%data(k,:,i) = 1.0
        end do
      end if
    end do
    

    deallocate(phi,  dxp, dzp)
  end subroutine set_2D_Eggs

  subroutine set_2D_Cu_wind_field(current_state)
    type(model_state_type), intent(inout), target :: current_state

    !
    ! Set up the 2D wind field for cumulus based on 
    ! Morrison and Grabowski (2007) 
    ! 
    ! Using formulation described in Appendix of
    ! MG07

    real(DEFAULT_PRECISION):: maxW

    !local variables
    real(DEFAULT_PRECISION) :: t     ! temporary local time variable
    real(DEFAULT_PRECISION) :: Lz
    real(DEFAULT_PRECISION), dimension(:), allocatable :: zp, dzp, dxp!(nz+1)
    !real(DEFAULT_PRECISION) :: xp(nx+1), dxp(nx+1)
    real(DEFAULT_PRECISION), dimension(:,:), allocatable :: phi!(0:nx+1,nz+1) 
                      ! streamfunction for cumulus 
                      ! convection
    !real(DEFAULT_PRECISION) :: ux(nx+1,nz+1), uz(nx+1,nz+1)
    real(DEFAULT_PRECISION) :: zscale1, zscale2
                      ! depth of the inflow and 
                      ! outflow (respectively)
    real(DEFAULT_PRECISION) :: ampl0, ampl20, xscale0, ampa, ampb &
         , amp2a, amp2b, tscale1, tscale2, ampl, ampl2 &
         , xscale, t1, ztop, x0, zz1, zz2 &
         , xl, xx, zl, xcen, zsh
    real(DEFAULT_PRECISION) :: scal_fac & ! factor to scale w and v to 0 above 
                           ! zscale2
         , dw0 & ! change in w from zscale2 to 0
         , dv0   ! change in v from zscale2 to 0
    
    integer :: nxmid, i, k, x

    real(DEFAULT_PRECISION) :: alpha, beta, hx, hz, z0

    allocate(phi(0:current_state%local_grid%local_domain_end_index(X_INDEX), &
         1:current_state%local_grid%local_domain_end_index(Z_INDEX)),  &
         zp(current_state%local_grid%local_domain_end_index(Z_INDEX)), &
         dzp(current_state%local_grid%local_domain_end_index(Z_INDEX)), &
         dxp(current_state%local_grid%local_domain_end_index(X_INDEX)))

    ! DETERMINE THE DEPTH OF THE INFLOW (ZSCALE1) 
    ! AND OUTFLOW (ZSCALE2)(in METERS)
    zscale1=1.7*1.e3
    zscale2=2.7*1.e3

    maxW=1.0
 
    ! CALCULATE X AND Z DISTANCES (IN METERS)    
    dzp=current_state%global_grid%configuration%vertical%dz    
    dxp = current_state%global_grid%configuration%horizontal%dx
    zp=current_state%global_grid%configuration%vertical%z

    dzp(1)=dzp(2)
    
    !
    ! INITIAL DATA FOR THE STREAMFUNCTION. AMPL IS WMAX IN M/S,
    ! AMPL2 IS THE VALUE OF LINEAR SHEAR OVER THE ZSCALE2 DEPTH (M/S)
    ! XSCALE IS WIDTH OF THE UPDRAFT IN M
    ! AMPA AND AMPB CONTROL THE TEMPORAL FLUCTUATIONS OF WMAX 
    ! AMP2A AND AMP2B CONTROL THE TEMPORAL FLUCTUATIONS OF SHEAR (TILT)
    ! TSCALE1 AND TSCALE2 ARE PERIODS OF COSINE FLUCTUATIONS OF THE ABOVE

    ! AMPL0=1.0 original value
    AMPL0=maxW
    AMPL20=0.
    XSCALE=1.8*1.e3
    AMPA=3.5
    AMPB=3.0
    AMP2A=0.6*1.e3
    AMP2B=0.5*1.e3
    tscale1=600.
    tscale2=900.

    t=current_state%time
    ! set parameters for each time period
    !
    !  AMPL AND AMPL2 VARYING IN TIME
    !
    if (t.lt.300.) then
      AMPL=AMPL0
      AMPL2=AMPL20
    elseif (t.ge.300..and.t.lt.900.) then
      AMPL=AMPL0
      AMPL2=AMP2A*(cos(pi*((t-300.)/tscale1 - 1.)) + 1.)
    elseif (t.ge.900..and.t.le.1500.) then
      AMPL=AMPA*(cos(pi*((t-900.)/tscale1 + 1.)) +1.) + AMPL0
      AMPL2=AMP2A*(cos(pi*((t-300.)/tscale1 - 1.)) + 1.)
    elseif (t.ge.1500..and.t.lt.2100.) then
      AMPL=AMPB*(cos(pi*(t-1500.)/tscale2) +1.) + 2.*AMPL0
      AMPL2=AMP2B*(cos(pi*((t-1500.)/tscale2 - 1.)) + 1.)
    elseif (t.ge.2100..and.t.lt.2400.) then
      AMPL=AMPB*(cos(pi*(t-1500.)/tscale2) +1.) + 2.*AMPL0
      AMPL2=AMP2B*(cos(pi*((t-1500.)/tscale2 - 1.)) + 1.)
    else
      t1=2400.
      AMPL=AMPB*(cos(pi*(t1-1500.)/tscale2) +1.) + 2.*AMPL0
      AMPL2=AMP2B*(cos(pi*((t1-1500.)/tscale2 - 1.)) + 1.)
    endif
    AMPL=AMPL/pi*XSCALE
    !
    ! DEFINE STREAMFUNCTION AS A FUNCTION OF HEIGHT
    ! ALSO, CENTRALIZE THE UPDRAFT TO OCCUPY ONLY THE INNER XSCALE
    ! OF THE DOMAIN
    !
    !ZTOP=zp(nz+1)/ZSCALE1
    XCEN=.5*(current_state%global_grid%configuration%horizontal%dx * current_state%global_grid%size(X_INDEX))
    !NXMID=current_state%global_grid%configuration%horizontal%nx/2+1
    X0=(current_state%global_grid%configuration%horizontal%dx * current_state%global_grid%size(X_INDEX)-XSCALE)/2.

    do i=current_state%local_grid%local_domain_start_index(X_INDEX)-1, current_state%local_grid%local_domain_end_index(X_INDEX)
      x = (0.5 + current_state%local_grid%start(X_INDEX) + (i-current_state%local_grid%local_domain_start_index(X_INDEX))) * &
         current_state%global_grid%configuration%horizontal%dx
      do k=1, current_state%local_grid%local_domain_end_index(Z_INDEX) 
#
         hx=xscale
          IF (abs(x - xcen) <= 0.5*xscale) THEN
             alpha=1.0
          ELSE
             alpha=0.0
          ENDIF
          IF (x <= xcen + 0.5*xscale)THEN
             beta=1.0
          ELSE
             beta=-1.0
          END IF
          IF (zp(k) <= zscale1)THEN
             z0=0.0
             hz=2.*zscale1
          ELSE
             z0=700.
             hz=2000.
          END IF
          PHI(i,k)=-cos(alpha*pi*(x-X0)/hx)*sin(beta*pi*(zp(k)-z0)/hz)
          !     ADD LINEAR SHEAR TO PRODUCE A WEAK TILT OF THE UPDRAFT
          ZSH=zp(k)/ZSCALE2
          PHI(i,k)=AMPL*PHI(i,k) - AMPL2*.5*ZSH**2.
      end do
    end do
    
    !  calculate rho*vel by derivation of streamfunction and normalize
    !  rho*ux velocity:
    do i=current_state%local_grid%local_domain_start_index(X_INDEX), current_state%local_grid%local_domain_end_index(X_INDEX)
      do k=2, current_state%local_grid%size(Z_INDEX)
        current_state%zu%data(k,:,i)=-(phi(i,k)-phi(i,k-1))/dzp(k) &
             /current_state%global_grid%configuration%vertical%rhon(k)
        current_state%zw%data(k,:,i)=(phi(i,k)-phi(i-1,k))/dxp(i)  &
             /current_state%global_grid%configuration%vertical%rho(k)
      enddo
      current_state%zu%data(1,:,i)=-current_state%zu%data(1,:,i)
      current_state%zw%data(1,:,i)=-current_state%zw%data(1,:,i)
    enddo

    current_state%u%data = current_state%zu%data
    current_state%w%data = current_state%zw%data

    deallocate(phi, zp, dxp, dzp)

  end subroutine set_2D_Cu_wind_field 

  subroutine set_2D_Sc(current_state)
    type(model_state_type), intent(inout), target :: current_state
    
    real(DEFAULT_PRECISION), dimension(:), allocatable :: dzp, dxp

    real(DEFAULT_PRECISION):: x, z

    real(DEFAULT_PRECISION), dimension(:,:), allocatable :: phi

    real(DEFAULT_PRECISION) :: renorm, maxW

    real(DEFAULT_PRECISION) :: zscale, xscale, x0

    integer :: i,j,k,iq

    maxW=1.

    allocate(phi(0:current_state%local_grid%local_domain_end_index(X_INDEX), &
       1:current_state%local_grid%local_domain_end_index(Z_INDEX)),  &
       dzp(current_state%local_grid%local_domain_end_index(Z_INDEX)), &
       dxp(current_state%local_grid%local_domain_end_index(X_INDEX)))
   
    dzp=current_state%global_grid%configuration%vertical%dz    
    dxp = current_state%global_grid%configuration%horizontal%dx
    dzp(1)=dzp(2)

    ! parameters - these might be read in through config
    
    xscale=current_state%global_grid%configuration%horizontal%dx * current_state%global_grid%size(X_INDEX) 
    zscale=current_state%global_grid%top(Z_INDEX) 

    x0=(current_state%global_grid%configuration%horizontal%dx * current_state%global_grid%size(X_INDEX)-xscale)/2.

    do i=current_state%local_grid%local_domain_start_index(X_INDEX)-1, current_state%local_grid%local_domain_end_index(X_INDEX)
      x = (0.5 + current_state%local_grid%start(X_INDEX) + (i-current_state%local_grid%local_domain_start_index(X_INDEX))) * &
         current_state%global_grid%configuration%horizontal%dx
      do k=1, current_state%local_grid%local_domain_end_index(Z_INDEX) 
        z = current_state%global_grid%configuration%vertical%z(k)
        
        phi(i,k) = -cos(2*pi*(x - x0)/xscale) * sin(pi*z/zscale)

      end do
    end do
    
    !  calculate rho*vel by derivation of streamfunction and normalize
    !  rho*ux velocity:
    do i=current_state%local_grid%local_domain_start_index(X_INDEX), current_state%local_grid%local_domain_end_index(X_INDEX)
      do k=2, current_state%local_grid%size(Z_INDEX)
        current_state%zu%data(k,:,i)=-(phi(i,k)-phi(i,k-1))/dzp(k) &
             /current_state%global_grid%configuration%vertical%rhon(k)
        current_state%zw%data(k,:,i)=(phi(i,k)-phi(i-1,k))/dxp(i)  &
             /current_state%global_grid%configuration%vertical%rho(k)
      enddo
      current_state%zu%data(1,:,i)=-current_state%zu%data(1,:,i)
      current_state%zw%data(1,:,i)=-current_state%zw%data(1,:,i)
    enddo

    ! Renormalize winds
    renorm = maxW/maxval(current_state%zw%data)
    
    current_state%zw%data = renorm*current_state%zw%data
    current_state%zu%data = renorm*current_state%zu%data

    current_state%u%data = current_state%zu%data
    current_state%w%data = current_state%zw%data


    deallocate(phi,  dxp, dzp)
  end subroutine set_2D_Sc


  subroutine set_2D_Squall(current_state)
    type(model_state_type), intent(inout), target :: current_state
    
    real(DEFAULT_PRECISION), dimension(:), allocatable :: dzp, dxp

    real(DEFAULT_PRECISION):: x, z

    real(DEFAULT_PRECISION), dimension(:,:), allocatable :: phi

    real(DEFAULT_PRECISION) :: x1, x2, z1, z2, z3, xc1, xc2

    real(DEFAULT_PRECISION) :: w_conv, w_strat_dn, w_strat_up

    real(DEFAULT_PRECISION) :: u_S, u_Q

    real(DEFAULT_PRECISION) :: xm1, xm2, zm1, zm2, zm3

    real(DEFAULT_PRECISION) :: renorm

    integer :: i,j,k,iq

    allocate(phi(0:current_state%local_grid%local_domain_end_index(X_INDEX), &
       1:current_state%local_grid%local_domain_end_index(Z_INDEX)),  &
       dzp(current_state%local_grid%local_domain_end_index(Z_INDEX)), &
       dxp(current_state%local_grid%local_domain_end_index(X_INDEX)))
   
    dzp=current_state%global_grid%configuration%vertical%dz    
    dxp = current_state%global_grid%configuration%horizontal%dx
    dzp(1)=dzp(2)

    ! parameters - these might be read in through config
    x1 = 10000.0  ! width of updraught
    x2 = 40000.0  ! width of stratiform region
    z1 = 10000.0  ! Top of convective region
    z2 = 4000.0   ! Top of stratiform region
    z3 = 500.     ! Bottom of stratiform region
    
    xc1 = 40000.  ! Location of center of convective region
    xc2 = 125000. ! Location of center of stratiform region
    
    w_conv     = 30.0 ! peak updraught in convective region
    w_strat_dn = 1.3  ! peak downdraught in stratus region
    w_strat_up = 1.5  ! peak updraught in stratus region
    
    u_S = 0.001       ! wind shear
    u_Q = 4.0         ! mean wind
    
    
    do i=current_state%local_grid%local_domain_start_index(X_INDEX)-1, current_state%local_grid%local_domain_end_index(X_INDEX)
      x = (0.5 + current_state%local_grid%start(X_INDEX) + (i-current_state%local_grid%local_domain_start_index(X_INDEX))) * &
         current_state%global_grid%configuration%horizontal%dx
      do k=1, current_state%local_grid%local_domain_end_index(Z_INDEX) 
        z = current_state%global_grid%configuration%vertical%z(k)
        
        ! convective part        
        xm1 = max(-x1, min(x1, x - xc1))
        xm2 = max(-x2, min(x2, x - xc2))
        zm1 = min(z, z1)
        zm2 = min(z2, max(0., z - z3))
        zm3 = min(z2, max(0., z - z3 - z2))
        phi(i,k) = (w_conv*2*x1/pi)*sin(0.5*pi*xm1/x1)*sin(pi*zm1/z1)

        ! large scale constant + shear
        phi(i,k) = phi(i,k) + (-0.5*u_S*z*z - u_Q*z)

        ! lower stratiform downdraught
        phi(i,k) = phi(i,k) - (w_strat_dn*2.0*x2/pi)*sin(0.5*pi*xm2/x2)*sin(pi*zm2/z2)**2 

        ! upper stratiform updraught
        phi(i,k) = phi(i,k) + (w_strat_up*2.0*x2/pi)*sin(0.5*pi*xm2/x2)*sin(pi*zm3/z2)**2 
      end do
    end do
    
    !  calculate rho*vel by derivation of streamfunction and normalize
    !  rho*ux velocity:
    do i=current_state%local_grid%local_domain_start_index(X_INDEX), current_state%local_grid%local_domain_end_index(X_INDEX)
      do k=2, current_state%local_grid%size(Z_INDEX)
        current_state%zu%data(k,:,i)=-(phi(i,k)-phi(i,k-1))/dzp(k) &
             /current_state%global_grid%configuration%vertical%rhon(k)
        current_state%zw%data(k,:,i)=(phi(i,k)-phi(i-1,k))/dxp(i)  &
             /current_state%global_grid%configuration%vertical%rho(k)
      enddo
      current_state%zu%data(1,:,i)=-current_state%zu%data(1,:,i)
      current_state%zw%data(1,:,i)=-current_state%zw%data(1,:,i)
    enddo

    ! Renormalize winds
    renorm = w_conv/maxval(current_state%zw%data)
    
    current_state%zw%data = renorm*current_state%zw%data
    current_state%zu%data = renorm*current_state%zu%data

    current_state%u%data = current_state%zu%data
    current_state%w%data = current_state%zw%data

    deallocate(phi,  dxp, dzp)
  end subroutine set_2D_Squall


  subroutine set_2D_Hills(current_state)
    type(model_state_type), intent(inout), target :: current_state
    
    real(DEFAULT_PRECISION), dimension(:), allocatable :: dzp, dxp

    real(DEFAULT_PRECISION):: x, z

    real(DEFAULT_PRECISION), dimension(:,:), allocatable :: phi

    real(DEFAULT_PRECISION) :: xa, xb, pratio, ha, hb, ca, cb, wa, wb, wmax, xcentre

    real(DEFAULT_PRECISION) :: Ra, Rb, thresh, zoff

    real(DEFAULT_PRECISION) :: znd, xnd

    real(DEFAULT_PRECISION) :: renorm

    integer :: i,j,k,iq

    allocate(phi(0:current_state%local_grid%local_domain_end_index(X_INDEX), &
       1:current_state%local_grid%local_domain_end_index(Z_INDEX)),  &
       dzp(current_state%local_grid%local_domain_end_index(Z_INDEX)), &
       dxp(current_state%local_grid%local_domain_end_index(X_INDEX)))
   
    dzp=current_state%global_grid%configuration%vertical%dz    
    dxp = current_state%global_grid%configuration%horizontal%dx
    dzp(1)=dzp(2)

    xcentre=0.5*(current_state%global_grid%configuration%horizontal%dx * current_state%global_grid%size(X_INDEX))

    !=========================================
    ! Set parameters defining the hills
    !=========================================
    pratio = 1.08  ! ratio of circleA radius to hillA height

    xa = -2000.  ! location of hill a 
    xb = 2000. ! location of hill b 
    xa = xa + xcentre
    xb = xb + xcentre

    ha = 2800.  ! height of hill a (+zoff) 
    hb = 2300.  ! height of hill b (+zoff)

    ca = ha/pratio  ! height of circle a
    cb = hb/pratio  ! height of circle b

    wa = 2000.  ! width of hill a
    wb = 2000.  ! width of hill a

    !=========================================
    ! Determine streamline to use as hills
    !=========================================

    Ra=pratio
    Rb=sqrt((xa-xb)*(xa-xb)/(wb*wb)+ha*ha/(cb*cb))
    thresh = -(ha*(1.-1./(Ra*Ra)) +  ha*(1.-1./(Rb*Rb)))

    zoff=-0.5*thresh
    !=========================================
    ! Set amplitude  
    !=========================================
    wmax = 4.
    
    do i=current_state%local_grid%local_domain_start_index(X_INDEX)-1, current_state%local_grid%local_domain_end_index(X_INDEX)
      x = (0.5 + current_state%local_grid%start(X_INDEX) + (i-current_state%local_grid%local_domain_start_index(X_INDEX))) * &
         current_state%global_grid%configuration%horizontal%dx
      xnd=x
      do k=1, current_state%local_grid%local_domain_end_index(Z_INDEX) 
        z = current_state%global_grid%configuration%vertical%z(k)

        znd=(z+zoff)
        Ra=sqrt((xnd-xa)*(xnd-xa)/(wa*wa)+znd*znd/(ca*ca))
        Rb=sqrt((xnd-xb)*(xnd-xb)/(wb*wb)+znd*znd/(cb*cb))
        phi(i,k) = -(znd*(1.-1./(Ra*Ra)) +  znd*(1.-1./(Rb*Rb)))


        if (phi(i,k) > .5*thresh)then
          phi(i,k) = .5*thresh
        end if

      end do
    end do
    
    !  calculate rho*vel by derivation of streamfunction and normalize
    !  rho*ux velocity:
    do i=current_state%local_grid%local_domain_start_index(X_INDEX), current_state%local_grid%local_domain_end_index(X_INDEX)
      do k=2, current_state%local_grid%size(Z_INDEX)
        current_state%zu%data(k,:,i)=-(phi(i,k)-phi(i,k-1))/dzp(k) &
             /current_state%global_grid%configuration%vertical%rhon(k)
        current_state%zw%data(k,:,i)=(phi(i,k)-phi(i-1,k))/dxp(i)  &
             /current_state%global_grid%configuration%vertical%rho(k)
      enddo
      current_state%zu%data(1,:,i)=-current_state%zu%data(1,:,i)
      current_state%zw%data(1,:,i)=-current_state%zw%data(1,:,i)
    enddo

    ! Renormalize winds
    renorm = wmax/maxval(current_state%zw%data)
    
    current_state%zw%data = renorm*current_state%zw%data
    current_state%zu%data = renorm*current_state%zu%data

    current_state%u%data = current_state%zu%data
    current_state%w%data = current_state%zw%data

    deallocate(phi,  dxp, dzp)
  end subroutine set_2D_Hills

  subroutine check_top(zztop, z, info)
    real(kind=DEFAULT_PRECISION), intent(in) :: zztop
    real(kind=DEFAULT_PRECISION), intent(in) :: z
    character(*), intent(in) :: info

    if (z<zztop)then
      call log_master_log(LOG_ERROR, "Top of input profile is below the top of the domain:"//trim(info))
    end if

  end subroutine check_top


end module kidtestcase_mod
