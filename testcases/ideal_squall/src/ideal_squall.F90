module ideal_squall_mod

  ! This module provides forcing for a squall simulation with sponge layers at the 
  ! limits of the bi-periodic domain (to prevent wrapping-round).  This should be 
  ! turned off if/when the model is developed to run with a stretched grid in the 
  ! horizontal.
       
  use datadefn_mod, only : DEFAULT_PRECISION
  use monc_component_mod, only : component_descriptor_type
  use state_mod, only : FORWARD_STEPPING, model_state_type
  use grids_mod, only : vertical_grid_configuration_type, X_INDEX, Y_INDEX, Z_INDEX
  use science_constants_mod, only : pi
  use maths_mod, only : random
  use optionsdatabase_mod, only :  options_get_real_array, options_get_real, &
     options_get_logical, options_get_array_size, options_get_string_array

  use q_indices_mod, only: get_q_index, standard_q_names

  implicit none

#ifndef TEST_MODE
  private
#endif
  
  real(kind=DEFAULT_PRECISION) :: shear_factor 
  real(kind=DEFAULT_PRECISION) :: shear_top 
  real(kind=DEFAULT_PRECISION) :: u_trans 
  real(kind=DEFAULT_PRECISION) :: evaporation_factor 
  real(kind=DEFAULT_PRECISION) :: alpha_george   
  real(kind=DEFAULT_PRECISION) :: force_width 
  real(kind=DEFAULT_PRECISION) :: force_z_centre 
  real(kind=DEFAULT_PRECISION) :: spinup_time 
  real(kind=DEFAULT_PRECISION) :: transition_time
  real(kind=DEFAULT_PRECISION) :: xdmp



  public ideal_squall_get_descriptor
contains

  type(component_descriptor_type) function ideal_squall_get_descriptor()
    ideal_squall_get_descriptor%name="ideal_squall"
    ideal_squall_get_descriptor%version=0.1
    ideal_squall_get_descriptor%initialisation=>initialisation_callback
    ideal_squall_get_descriptor%timestep=>timestep_callback
  end function ideal_squall_get_descriptor

  !! Note that this is not the most efficient way to iterate through theta (j heavy), but it is the same as the LEM set up 
  !! so directly comparable and probably doesn't matter too much as it is just called onec in the initialisation
  subroutine initialisation_callback(current_state)
    type(model_state_type), intent(inout), target :: current_state

    type(vertical_grid_configuration_type) :: vertical

    integer :: i,j,k ! loop counters

    vertical=current_state%global_grid%configuration%vertical

    ! Read in parameters from options database
    shear_factor=options_get_real(current_state%options_database, "shear_factor") 
    shear_top=options_get_real(current_state%options_database, "shear_top") 
    evaporation_factor=options_get_real(current_state%options_database, "evaporation_factor") 
    alpha_george=options_get_real(current_state%options_database, "alpha_george") 
    force_width=options_get_real(current_state%options_database, "force_width") 
    force_z_centre=options_get_real(current_state%options_database, "force_z_centre") 
    spinup_time=options_get_real(current_state%options_database, "spinup_time") 
    transition_time=options_get_real(current_state%options_database, "transition_time") 
    xdmp=options_get_real(current_state%options_database, "xdmp") 
    
    
    ! Set up initial winds (component initialisation should be called after
    ! Gridmanager which also sets up the winds.
    ! Squall is set up in to vary in the x direction and be quasi-uniform in y

    vertical%v_init(:) = 0.0
    do k=current_state%local_grid%local_domain_start_index(Z_INDEX), &
             current_state%local_grid%local_domain_end_index(Z_INDEX)  
      if (vertical%z(k) < shear_top)then
        vertical%u_init(k) = shear_factor*vertical%z(k)/shear_top
      else
        vertical%u_init(k) = shear_factor
      end if
!      current_state%u%data(k,:,:) = vertical%u_init(k)
    end do

!    do i=current_state%local_grid%local_domain_start_index(X_INDEX), current_state%local_grid%local_domain_end_index(X_INDEX)
!      do j=current_state%local_grid%local_domain_start_index(Y_INDEX), current_state%local_grid%local_domain_end_index(Y_INDEX)
!        current_state%u%data(:,j,i) = vertical%u_init(:)
!      end do
!    end do
  end subroutine initialisation_callback

  subroutine timestep_callback(current_state)
    type(model_state_type), intent(inout), target :: current_state

    type(vertical_grid_configuration_type):: vertical

    real(kind=DEFAULT_PRECISION) :: x, y       ! x, y position
    real(kind=DEFAULT_PRECISION) :: xoff, zoff ! offsets from centre of forcing region
    real(kind=DEFAULT_PRECISION) :: xsize, xc  ! size,centre of x domain
    real(kind=DEFAULT_PRECISION) :: dtm        ! current timestep


    real(kind=DEFAULT_PRECISION) :: gamma_george   
    real(kind=DEFAULT_PRECISION) :: damping_coefficient
  
    integer :: iq,k ! loop counters
    INTEGER :: icol, jcol
    if (current_state%halo_column) return

    icol=current_state%column_local_x
    jcol=current_state%column_local_y
      
    x = (current_state%local_grid%start(X_INDEX) + (icol-current_state%local_grid%local_domain_start_index(X_INDEX))) * &
       current_state%global_grid%configuration%horizontal%dx
    y = (current_state%local_grid%start(Y_INDEX) + (jcol-current_state%local_grid%local_domain_start_index(Y_INDEX))) * &
       current_state%global_grid%configuration%horizontal%dy 

    dtm = current_state%dtm*2.0
    if (current_state%field_stepping == FORWARD_STEPPING) dtm=current_state%dtm

    vertical=current_state%global_grid%configuration%vertical

    xsize = current_state%global_grid%size(X_INDEX)*current_state%global_grid%configuration%horizontal%dx
    xc = 0.5*xsize

    if (current_state%time < spinup_time)then
      gamma_george=1.0
      if (current_state%time > spinup_time-transition_time)then
        gamma_george= 1.0 - (current_state%time + transition_time - spinup_time)/transition_time 
      end if
      
      xoff=(x - xc)/force_width
      do k=current_state%local_grid%local_domain_start_index(Z_INDEX), &
         current_state%local_grid%local_domain_end_index(Z_INDEX)  
        if (vertical%z(k) < 2.0*force_z_centre)then
          zoff=vertical%z(k)/force_z_centre
          if (abs(xoff) < 1.0 .and. abs(zoff) < 1.0)then
            current_state%su%data(k,jcol,icol) = current_state%su%data(k,jcol,icol) + &
               alpha_george*gamma_george*cos(0.5*pi*xoff)**2 &
               / ((cosh(2.5*zoff)**2)**2)
          end if
        end if
      end do
    end if

!    Damping at the sides - we're currently stuck with bi-periodic boundaries, so let's do a little fudge
    
    
    if (x < xdmp .or. x > xsize - xdmp)then
      if (x < xdmp)damping_coefficient = exp(-8.0*((2.*x-xdmp)/xdmp)**2)
      if (x > xsize - xdmp)damping_coefficient = exp(-8.0*((2.*(xsize-x)-xdmp)/xdmp)**2)
      
      do k=current_state%local_grid%local_domain_start_index(Z_INDEX), &
         current_state%local_grid%local_domain_end_index(Z_INDEX)  
        
        do iq=1,current_state%number_q_fields
          current_state%sq(iq)%data(k,jcol,icol) = &
             current_state%sq(iq)%data(k,jcol,icol)*(1.-damping_coefficient) &
             - damping_coefficient*(current_state%zq(iq)%data(k,jcol,icol)  &
                   - vertical%q_init(k,iq))/dtm
        end do
        current_state%sth%data(k,jcol,icol) = &
           current_state%sth%data(k,jcol,icol)*(1.-damping_coefficient) &
           - damping_coefficient*(current_state%zth%data(k,jcol,icol) + vertical%thref(k)  &
           - vertical%theta_init(k))/dtm
        current_state%su%data(k,jcol,icol) = &
           current_state%su%data(k,jcol,icol)*(1.-damping_coefficient) &
           - damping_coefficient*(current_state%zu%data(k,jcol,icol)  &
           - vertical%u_init(k) + current_state%ugal)/dtm
        ! NB not relaxing v winds
        current_state%sw%data(k,jcol,icol) = &
           current_state%sw%data(k,jcol,icol)*(1.-damping_coefficient) &
           - damping_coefficient*(current_state%zw%data(k,jcol,icol))/dtm
      end do
    end if

  end subroutine timestep_callback

end module ideal_squall_mod
