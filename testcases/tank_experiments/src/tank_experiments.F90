module tank_experiments_mod

  ! This is essentially the same as cold tank_experiments, but we also fix the environmental 
  ! relative humidity by adjusting the vapour content.
       
  use datadefn_mod, only : DEFAULT_PRECISION, STRING_LENGTH
  use monc_component_mod, only : component_descriptor_type
  use state_mod, only : model_state_type
  use grids_mod, only : X_INDEX, Y_INDEX, Z_INDEX
  use science_constants_mod, only : pi, r_over_cp
  use maths_mod, only : random
  use optionsdatabase_mod, only :  options_get_real_array, options_get_real, options_get_logical, &
     options_get_array_size, options_get_string_array
  use saturation_mod, only: qsaturation
  use q_indices_mod, only: get_q_index, standard_q_names
  implicit none

#ifndef TEST_MODE
  private
#endif
  
  integer, parameter :: UNSET_REAL=-999.0
  integer, parameter :: MAXBUBBLES=10 !maximum number of bubbles
  integer, parameter :: MAXONOFF=10   !maximum number of times for switching on or off
  integer, parameter :: MAXQIN=10     ! maximum number of q variables to initialize

  real(kind=DEFAULT_PRECISION), dimension(MAXBUBBLES) :: x_cen=UNSET_REAL, y_cen, z_cen   ! coordinates of tank_experiments centre
  real(kind=DEFAULT_PRECISION), dimension(MAXBUBBLES) :: x_rad, y_rad, z_rad   ! radial parameters  
  real(kind=DEFAULT_PRECISION), dimension(MAXBUBBLES) :: th_pert               ! maximum theta perturbation
  real(kind=DEFAULT_PRECISION), dimension(MAXBUBBLES) :: t_pert                ! maximum temperature perturbation (can be used instead of theta)
  real(kind=DEFAULT_PRECISION), dimension(MAXBUBBLES) :: edge                  ! fraction of radius for edge smoothing

  ! humidity
  real(kind=DEFAULT_PRECISION) :: RH                                  ! initial relative humidity of environment
  real(kind=DEFAULT_PRECISION), dimension(MAXBUBBLES) :: bubble_RH    ! initial relative humidity of bubble

  ! q variables
  real(kind=DEFAULT_PRECISION), dimension(MAXBUBBLES*MAXQIN) :: bubble_q_values_tmp ! initial mixing ratios of bubble q variables
  real(kind=DEFAULT_PRECISION), allocatable :: bubble_q_values(:,:) ! initial mixing ratios of bubble q variables
  character(len=STRING_LENGTH) :: bubble_q_names(MAXQIN)='unset'  ! names of q variables to initialize

  ! Source descriptions
  real(kind=DEFAULT_PRECISION), dimension(MAXBUBBLES) :: x_cen_src=UNSET_REAL, y_cen_src, z_cen_src   ! coordinates of tank_experiments centre
  real(kind=DEFAULT_PRECISION), dimension(MAXBUBBLES) :: x_rad_src, y_rad_src, z_rad_src   ! radial parameters  
  real(kind=DEFAULT_PRECISION), dimension(MAXBUBBLES) :: edge_src                          ! fraction of radius for edge smoothing
  real(kind=DEFAULT_PRECISION), dimension(MAXBUBBLES) :: dth_src                           ! theta source
  real(kind=DEFAULT_PRECISION), allocatable :: onoff_src(:,:)   ! times for switching sources
  real(kind=DEFAULT_PRECISION), dimension(MAXBUBBLES*MAXONOFF):: onoff_src_tmp=UNSET_REAL    ! times for switching sources

  real(kind=DEFAULT_PRECISION) :: barrier_x ! x Location of  "barrier"
  real(kind=DEFAULT_PRECISION) :: barrier_y ! x Location of  "barrier"
  real(kind=DEFAULT_PRECISION) :: left_tank_delta_theta  ! temperature difference in left side of the tank
  real(kind=DEFAULT_PRECISION) :: right_tank_delta_theta ! temperature difference in right side of the tank
  real(kind=DEFAULT_PRECISION) :: back_tank_delta_theta  ! temperature difference in back side of the tank
  real(kind=DEFAULT_PRECISION) :: front_tank_delta_theta ! temperature difference in front side of the tank
  real(kind=DEFAULT_PRECISION) :: left_tank_delta_rh  ! rh difference in left side of the tank
  real(kind=DEFAULT_PRECISION) :: right_tank_delta_rh ! rh difference in right side of the tank
  real(kind=DEFAULT_PRECISION) :: back_tank_delta_rh  ! rh difference in back side of the tank
  real(kind=DEFAULT_PRECISION) :: front_tank_delta_rh ! rh difference in front side of the tank

  integer :: number_of_bubbles ! The number of bubbles
  integer :: number_of_sources ! The number of sources
  integer :: number_of_onoff   ! The number of times to switch sources on or off

  integer :: nq_bubbles ! Number of q variables to initialize in the bubbles

  logical, dimension(MAXBUBBLES) :: source_on=.false. 
  
  logical :: l_bubbles   ! Add bubbles
  logical :: l_splittank ! Split the tank
  logical :: l_sources   ! Add sources

  ! Use this value of PI here (slightly different from the scientific constants to match the LEM tank_experiments setup)
  real(kind=DEFAULT_PRECISION), parameter :: my_pi = 4.0_DEFAULT_PRECISION*atan(1.0_DEFAULT_PRECISION)

  logical :: l_straka   ! Use the cosine form as in the Straka test
  logical :: l_use_t    ! Use temperature perturbation instead of theta

  logical :: l_random     ! set this to .false. to turn off random swapping
  logical :: l_moist      ! if .true. then a moist case will be used, initialized with RH values above
  integer, parameter :: n_swap = 100

  integer :: iqv ! index for vapour

  public tank_experiments_get_descriptor
contains

  type(component_descriptor_type) function tank_experiments_get_descriptor()
    tank_experiments_get_descriptor%name="tank_experiments"
    tank_experiments_get_descriptor%version=0.1
    tank_experiments_get_descriptor%initialisation=>initialisation_callback
    tank_experiments_get_descriptor%timestep=>timestep_callback
    tank_experiments_get_descriptor%finalisation=>finalisation_callback
  end function tank_experiments_get_descriptor

  !! Note that this is not the most efficient way to iterate through theta (j heavy), but it is the same as the LEM set up 
  !! so directly comparable and probably doesn't matter too much as it is just called onec in the initialisation
  subroutine initialisation_callback(current_state)
    type(model_state_type), intent(inout), target :: current_state

    integer :: i ! loop counter

    ! Read in parameters from options database
    l_bubbles=options_get_logical(current_state%options_database, "l_bubbles") 
    if (l_bubbles)then
      l_straka=options_get_logical(current_state%options_database, "l_bubble_straka") 
      l_use_t=options_get_logical(current_state%options_database, "l_bubble_use_t") 
    
      call options_get_real_array(current_state%options_database, "bubble_x_cen", x_cen) 
      call options_get_real_array(current_state%options_database, "bubble_y_cen", y_cen) 
      call options_get_real_array(current_state%options_database, "bubble_z_cen", z_cen) 
      call options_get_real_array(current_state%options_database, "bubble_x_rad", x_rad) 
      call options_get_real_array(current_state%options_database, "bubble_y_rad", y_rad) 
      call options_get_real_array(current_state%options_database, "bubble_z_rad", z_rad) 
      if (l_use_t)then
        call options_get_real_array(current_state%options_database, "bubble_t_pert", t_pert) 
      else
        call options_get_real_array(current_state%options_database, "bubble_th_pert", th_pert) 
      end if

      call options_get_real_array(current_state%options_database, "bubble_edge", edge) 
      call options_get_real_array(current_state%options_database, "bubble_RH", bubble_RH)
      l_random=options_get_logical(current_state%options_database, "bubble_lrandom")

      do i=1,size(x_cen)
        if (x_cen(i) == UNSET_REAL) exit 
      end do
      number_of_bubbles=i-1

      nq_bubbles=options_get_array_size(current_state%options_database, "bubble_q_names")
      if (nq_bubbles >  0)then
        call options_get_string_array(current_state%options_database, "bubble_q_names", bubble_q_names)
        
        call options_get_real_array(current_state%options_database, "bubble_q_values", bubble_q_values_tmp)
        allocate(bubble_q_values(nq_bubbles, number_of_bubbles))
        bubble_q_values=reshape(bubble_q_values_tmp, (/ nq_bubbles, number_of_bubbles/))
      end if

    end if

    l_moist=options_get_logical(current_state%options_database, "tank_lmoist")
    if (l_moist)RH=options_get_real(current_state%options_database, "tank_RH")

    l_splittank=options_get_logical(current_state%options_database, "l_splittank") 
    if (l_splittank)then
      left_tank_delta_theta=options_get_real(current_state%options_database, "left_tank_delta_theta") 
      right_tank_delta_theta=options_get_real(current_state%options_database, "right_tank_delta_theta") 
      front_tank_delta_theta=options_get_real(current_state%options_database, "front_tank_delta_theta") 
      back_tank_delta_theta=options_get_real(current_state%options_database, "back_tank_delta_theta") 
      left_tank_delta_rh=options_get_real(current_state%options_database, "left_tank_delta_rh") 
      right_tank_delta_rh=options_get_real(current_state%options_database, "right_tank_delta_rh") 
      front_tank_delta_rh=options_get_real(current_state%options_database, "front_tank_delta_rh") 
      back_tank_delta_rh=options_get_real(current_state%options_database, "back_tank_delta_rh") 
      barrier_x=options_get_real(current_state%options_database, "barrier_x") 
      barrier_y=options_get_real(current_state%options_database, "barrier_y") 
    end if

    ! Add in the q variable index for vapour
    if (l_moist) iqv=get_q_index(standard_q_names%VAPOUR, 'tank_experiments')

    if (.not. current_state%continuation_run) call generate_bubbles(current_state)    

    ! Sources
    l_sources=options_get_logical(current_state%options_database, "l_sources") 
    if (l_sources)then
      call options_get_real_array(current_state%options_database, "source_x_cen", x_cen_src) 
      call options_get_real_array(current_state%options_database, "source_y_cen", y_cen_src) 
      call options_get_real_array(current_state%options_database, "source_z_cen", z_cen_src) 
      call options_get_real_array(current_state%options_database, "source_x_rad", x_rad_src) 
      call options_get_real_array(current_state%options_database, "source_y_rad", y_rad_src) 
      call options_get_real_array(current_state%options_database, "source_z_rad", z_rad_src) 
      call options_get_real_array(current_state%options_database, "source_edge", edge_src) 
      call options_get_real_array(current_state%options_database, "source_dth", dth_src) 
      call options_get_real_array(current_state%options_database, "source_onoff", onoff_src_tmp) 
      ! call options_get_string_array(current_state%options_database, "source_q_names", names_q_src)
      ! call options_get_real_array(current_state%options_database, "source_dq", dq_src) 
      do i=1,size(x_cen_src)
        if (x_cen_src(i) == UNSET_REAL) exit 
      end do
      number_of_sources=i-1
      
      do i=1,size(onoff_src_tmp)
        if (onoff_src_tmp(i) == UNSET_REAL) exit 
      end do
      number_of_onoff=(i-1)/number_of_sources
      allocate(onoff_src(number_of_onoff, number_of_sources))
      onoff_src=reshape(onoff_src_tmp(1:i-1), (/ number_of_onoff, number_of_sources/))
    endif 

  end subroutine initialisation_callback

  subroutine generate_bubbles(current_state)
    type(model_state_type), intent(inout), target :: current_state

    real(kind=DEFAULT_PRECISION) :: x, y              ! x and y coordinate
    real(kind=DEFAULT_PRECISION) :: rad            ! radial distance    
    real(kind=DEFAULT_PRECISION) :: xi, delx, dely, delz  ! useful things
    
    integer :: i,j,k,ibub,iq,n ! loop counters

    real(kind=DEFAULT_PRECISION) :: qsat  ! temporary qsat value
    real(kind=DEFAULT_PRECISION) :: exner ! temporary exner value
    real(kind=DEFAULT_PRECISION) :: TdegK ! temporary temperature value in Kelvin

    integer :: i_dum, j_first, j_second
    logical, allocatable :: l_start(:)
    integer, allocatable :: j_start(:), j_end(:)
    real(kind=DEFAULT_PRECISION) :: th_temp, q_temp, RH_tank, th_pert_loc

    ! random number bits
    allocate(j_start(number_of_bubbles))
    allocate(j_end(number_of_bubbles))
    allocate(l_start(number_of_bubbles))
    l_start(:) = .true.
    i_dum   = -10.0_DEFAULT_PRECISION * current_state%parallel%my_rank
    j_start(:) = 0
    j_end(:)   = 0

    do i=current_state%local_grid%local_domain_start_index(X_INDEX), current_state%local_grid%local_domain_end_index(X_INDEX)
      x = (current_state%local_grid%start(X_INDEX) + (i-current_state%local_grid%local_domain_start_index(X_INDEX))) * &
         current_state%global_grid%configuration%horizontal%dx
      do k=current_state%local_grid%local_domain_start_index(Z_INDEX), current_state%local_grid%local_domain_end_index(Z_INDEX)
        exner=(current_state%global_grid%configuration%vertical%prefn(k)/100000.)**(r_over_cp)
        do j=current_state%local_grid%local_domain_start_index(Y_INDEX), current_state%local_grid%local_domain_end_index(Y_INDEX)          
          y = (current_state%local_grid%start(Y_INDEX) + (j-current_state%local_grid%local_domain_start_index(Y_INDEX))) * &
               current_state%global_grid%configuration%horizontal%dy       
          
          if (l_splittank)then
            if (x < barrier_x)current_state%th%data(k,j,i) = left_tank_delta_theta
            if (x >= barrier_x)current_state%th%data(k,j,i) = right_tank_delta_theta
            
            if (y < barrier_y)current_state%th%data(k,j,i) = current_state%th%data(k,j,i)+back_tank_delta_theta
            if (y >= barrier_y)current_state%th%data(k,j,i) = current_state%th%data(k,j,i)+front_tank_delta_theta
          end if

          if (l_moist)then
            TdegK=(current_state%global_grid%configuration%vertical%thref(k)+current_state%th%data(k,j,i))*exner	
            qsat=qsaturation(TdegK, current_state%global_grid%configuration%vertical%prefn(k)/100.)
            RH_tank=RH
            if (l_splittank)then
              if (x < barrier_x)RH_tank = RH_tank+left_tank_delta_rh
              if (x >= barrier_x)RH_tank = RH_tank+right_tank_delta_rh
              
              if (y < barrier_y)RH_tank = RH_tank+back_tank_delta_rh
              if (y >= barrier_y)RH_tank = RH_tank+front_tank_delta_rh
            end if
            current_state%q(iqv)%data(k,j,i)  = qsat*RH_tank/100.
          end if

          if (l_bubbles)then
            do ibub=1,number_of_bubbles
              delx = ( x - x_cen(ibub) ) / x_rad(ibub)
              delz = (current_state%global_grid%configuration%vertical%zn(k) - z_cen(ibub)) / z_rad(ibub)   
              dely = ( y - y_cen(ibub) ) / y_rad(ibub)          
              rad = sqrt(delx*delx + dely*dely + delz*delz)

              if (l_use_t)then
                th_pert_loc = t_pert(ibub)/exner
              else
                th_pert_loc = th_pert(ibub)
              end if
              
              if ( rad .le. 1.0_DEFAULT_PRECISION ) then
                ! record start and end values of j index
                j_end(ibub)=j
                if (l_start(ibub)) then
                  j_start(ibub)=j
                  l_start(ibub)=.false.
                end if
                if (rad .le. 1.0_DEFAULT_PRECISION - edge(ibub)) then
                  if (l_straka)then
                    current_state%th%data(k,j,i) = current_state%th%data(k,j,i)+th_pert_loc*0.5*(1.0 + cos(pi*rad))
                  else
                    current_state%th%data(k,j,i) = current_state%th%data(k,j,i)+th_pert_loc
                  end if
                else
                  xi = cos(0.5_DEFAULT_PRECISION * my_pi * (rad-(1.0_DEFAULT_PRECISION-edge(ibub)))/edge(ibub))              
                  xi = min(xi,1.0_DEFAULT_PRECISION)
                  xi = max(xi,0.0_DEFAULT_PRECISION)
                  current_state%th%data(k,j,i) = current_state%th%data(k,j,i)+th_pert_loc * xi * xi
                end if
                if (l_moist)then
                  TdegK=(current_state%global_grid%configuration%vertical%thref(k)+current_state%th%data(k,j,i))*exner	
                  qsat=qsaturation(TdegK, current_state%global_grid%configuration%vertical%prefn(k)/100.)
                  current_state%q(iqv)%data(k,j,i)  = qsat*bubble_RH(ibub)/100.
                  do n=1,nq_bubbles
                    iq=get_q_index(trim(bubble_q_names(n)), 'bubble_initialization')
                    current_state%q(iq)%data(k,j,i) = bubble_q_values(n, ibub)
                  end do

                end if
              end if
              
            end do
          end if
        end do
        if (l_bubbles)then
          do ibub=1,number_of_bubbles
            ! reset logical after loop
            l_start=.true.
            ! randomly swap around the values within the warm bubble
            if (l_random) then
              if ((j_end(ibub)-j_start(ibub)) .gt. 0) then
                do j=1, n_swap
                  j_first  = j_start(ibub) + int((j_end(ibub)-j_start(ibub)) * random(i_dum))
                  j_second = j_start(ibub) + int((j_end(ibub)-j_start(ibub)) * random(i_dum))
                  th_temp = current_state%th%data(k,j_first,i)
                  current_state%th%data(k,j_first,i) = current_state%th%data(k,j_second,i)
                  current_state%th%data(k,j_second,i) = th_temp
                  if (l_moist)then ! This maintains RH
                    q_temp = current_state%q(iqv)%data(k,j_first,i)
                    current_state%q(iqv)%data(k,j_first,i) = current_state%q(iqv)%data(k,j_second,i)
                    current_state%q(iqv)%data(k,j_second,i) = q_temp
                  end if
                end do
              end if
            end if
          end do
        end if
      end do
    end do

    ! Tidy up
    deallocate(j_start)
    deallocate(j_end)
    deallocate(l_start)

    if (nq_bubbles >  0)deallocate(bubble_q_values)


  end subroutine generate_bubbles

  subroutine timestep_callback(current_state)
    type(model_state_type), intent(inout), target :: current_state

    real(kind=DEFAULT_PRECISION) :: x, y              ! x and y coordinate
    real(kind=DEFAULT_PRECISION) :: rad            ! radial distance    
    real(kind=DEFAULT_PRECISION) :: xi, delx, dely, delz  ! useful things

    integer :: isrc, ion, i,j,k ! look counters


    if (current_state%first_timestep_column)then
      do isrc=1,number_of_sources
        if (mod(count(onoff_src(1:number_of_onoff, isrc) < current_state%time),2)==1)then
          source_on(isrc)=.true.
        else
          source_on(isrc)=.false.
        end if
      end do
    end if

    if (current_state%halo_column .or. .not. l_sources) return


    do isrc=1,number_of_sources
      if (.not. source_on(isrc)) cycle
      do i=current_state%local_grid%local_domain_start_index(X_INDEX), current_state%local_grid%local_domain_end_index(X_INDEX)
        x = (current_state%local_grid%start(X_INDEX) + (i-current_state%local_grid%local_domain_start_index(X_INDEX))) * &
           current_state%global_grid%configuration%horizontal%dx
        do k=current_state%local_grid%local_domain_start_index(Z_INDEX), current_state%local_grid%local_domain_end_index(Z_INDEX)
          do j=current_state%local_grid%local_domain_start_index(Y_INDEX), current_state%local_grid%local_domain_end_index(Y_INDEX)          
            y = (current_state%local_grid%start(Y_INDEX) + (j-current_state%local_grid%local_domain_start_index(Y_INDEX))) * &
               current_state%global_grid%configuration%horizontal%dy       
          
            delx = ( x - x_cen_src(isrc) ) / x_rad_src(isrc)
            delz = (current_state%global_grid%configuration%vertical%zn(k) - z_cen_src(isrc)) / z_rad_src(isrc)   
            dely = ( y - y_cen_src(isrc) ) / y_rad_src(isrc)          
            rad = sqrt(delx*delx + dely*dely + delz*delz)
            
            if ( rad .le. 1.0_DEFAULT_PRECISION ) then
              if (rad .le. 1.0_DEFAULT_PRECISION - edge_src(isrc)) then
                current_state%sth%data(k,j,i) = current_state%sth%data(k,j,i)+dth_src(isrc)
              else
                xi = cos(0.5_DEFAULT_PRECISION * my_pi * (rad-(1.0_DEFAULT_PRECISION-edge(isrc)))/edge(isrc))              
                xi = min(xi,1.0_DEFAULT_PRECISION)
                xi = max(xi,0.0_DEFAULT_PRECISION)
                current_state%sth%data(k,j,i) = current_state%sth%data(k,j,i)+dth_src(isrc) * xi * xi
              end if
            end if

          end do
        end do
      end do
    end do
    

  end subroutine timestep_callback

  subroutine finalisation_callback(current_state)
    type(model_state_type), intent(inout), target :: current_state
    
    if (allocated(onoff_src)) deallocate(onoff_src)
  end subroutine finalisation_callback
  
end module tank_experiments_mod
