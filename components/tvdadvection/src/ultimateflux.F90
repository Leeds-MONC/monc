!> Calculates the effective face values for advection using Leonard's ultimate quickest scheme with
!! first multi-dimension limiter and gradient terms added to the quickest scheme.
!!
!! Note that all face courant numbers (e.g. vface) have a factor of 0.5 included. This routine calculates
!! the right/top face flux although it saves it in FFLXL(J+1) etc...
module ultimateflux_mod
  use prognostics_mod, only : prognostic_field_type
  use grids_mod, only : local_grid_type, grid_configuration_type, X_INDEX, Y_INDEX, Z_INDEX
  use state_mod, only : parallel_state_type
  use datadefn_mod, only : DEFAULT_PRECISION
  implicit none

#ifndef TEST_MODE
  private
#endif

  public ultflx

  real(kind=DEFAULT_PRECISION) :: r6 !< A sixth used when calculating fluxes in specific directions

contains

  subroutine ultflx(y_flow_index, x_flow_index, u, v, w, y_scalar_index, x_scalar_index, zf, &
       local_grid, grid_config, parallel, kdof, dt, flux_y, flux_z, flux_x, flux_previous_x, rdz, rdzn, dzn, kmin, kmax)

    integer, intent(in) ::y_flow_index, x_flow_index, y_scalar_index, x_scalar_index   &! loop counter
         ,kdof,&                  ! =1 for advection of W, 0 otherwise to shift vertical grid index
         kmin, kmax
    real(kind=DEFAULT_PRECISION), intent(in) ::dt   ! timestep (s)
    type(prognostic_field_type), intent(inout) :: u, w, v, zf
    type(grid_configuration_type), intent(inout) :: grid_config
    type(local_grid_type), intent(inout) :: local_grid
    type(parallel_state_type), intent(inout) :: parallel
    real(kind=DEFAULT_PRECISION),  intent(in), dimension(:) :: rdz, rdzn, dzn
    real(kind=DEFAULT_PRECISION),  intent(inout), dimension(:) ::&
         flux_z        &! flux through bottom cell face
         ,flux_y         &! flux through left (y-dirn) cell face
         ,flux_x       &   ! flux through left (x-dirn) cell face
         ,flux_previous_x
    r6 = 1.0_DEFAULT_PRECISION/6.0_DEFAULT_PRECISION

#ifdef W_ACTIVE
    call handle_vertical_fluxes(y_flow_index, x_flow_index, u, v, w, y_scalar_index, x_scalar_index, zf, &
         local_grid, grid_config, kdof, dt, flux_z, rdz, rdzn, dzn, kmin, kmax)
#endif

#ifdef V_ACTIVE
    if (y_scalar_index .ne. local_grid%local_domain_end_index(Y_INDEX) .or. parallel%my_coords(Y_INDEX) .ne. &
         parallel%dim_sizes(Y_INDEX)-1) then 
      call handle_y_direction_fluxes(y_flow_index, x_flow_index, u, v, w, y_scalar_index, x_scalar_index, zf, local_grid, &
           grid_config, kdof, dt, flux_y, rdz, kmin, kmax)
    end if
#endif

#ifdef U_ACTIVE  
    flux_x(:) = flux_previous_x(:)  
    call handle_x_direction_fluxes(y_flow_index, x_flow_index, u, v, w, y_scalar_index, x_scalar_index, zf, &
         local_grid, grid_config, kdof, dt, flux_previous_x, rdz, kmin, kmax)
#endif
  end subroutine ultflx

  subroutine handle_vertical_fluxes(y_flow_index, x_flow_index, u, v, w, y_scalar_index, x_scalar_index, zf, &
       local_grid, grid_config, kdof, dt, flux_z, rdz, rdzn, dzn, kmin, kmax)
    integer, intent(in) :: kdof, y_flow_index, x_flow_index, y_scalar_index, x_scalar_index, kmin, kmax
    real(kind=DEFAULT_PRECISION), intent(in) :: dt
    type(prognostic_field_type), intent(inout) :: u, w, v, zf
    type(grid_configuration_type), intent(inout) :: grid_config
    type(local_grid_type), intent(inout) :: local_grid
    real(kind=DEFAULT_PRECISION), intent(in), dimension(:) :: rdz, rdzn, dzn
    real(kind=DEFAULT_PRECISION), intent(inout), dimension(:) :: flux_z

    if (kmin==1) call handle_vertical_fluxes_bottomofcolumn(y_flow_index, x_flow_index, u, v, w, y_scalar_index, x_scalar_index, &
         zf, grid_config, kdof, dt, flux_z, rdz, rdzn, dzn)
    call handle_vertical_fluxes_topofcolumn(y_flow_index, x_flow_index, u, v, w, y_scalar_index, x_scalar_index, zf, &
         local_grid, grid_config, kdof, dt, flux_z, rdz, rdzn, dzn)
    call handle_vertical_fluxes_middleofcolumn(y_flow_index, x_flow_index, u, v, w, y_scalar_index, x_scalar_index, zf, &
         local_grid, grid_config, kdof, dt, flux_z, rdz, rdzn, dzn)
    if (kmin .gt. 1) flux_z(2)=0.0_DEFAULT_PRECISION
  end subroutine handle_vertical_fluxes

  subroutine handle_x_direction_fluxes(y_flow_index, x_flow_index, u, v, w, y_scalar_index, x_scalar_index, zf, local_grid, &
       grid_config, kdof, dt, flux_x, rdz, kmin, kmax)
    integer, intent(in) :: y_flow_index, x_flow_index, y_scalar_index, x_scalar_index, kdof, kmin, kmax
    real(kind=DEFAULT_PRECISION), intent(in) :: dt
    type(prognostic_field_type), intent(inout) :: u, w, v, zf
    type(grid_configuration_type), intent(inout) :: grid_config
    type(local_grid_type), intent(inout) :: local_grid
    real(kind=DEFAULT_PRECISION), intent(in), dimension(:) :: rdz
    real(kind=DEFAULT_PRECISION), intent(out), dimension(:) :: flux_x

    integer :: k, kp1, kp0

    real(kind=DEFAULT_PRECISION) :: advneg,advpos,& !< Sign indicators (0 or 1)
         fgt1, fgt2,&  !< Gradient terms
         fc,& !< Upwinded nodal points in z
         fcurvs,& !< Curvature of F
         fd,& !< Downwind nodal points in z
         fdels,&  !< Downwind-upwind difference of F
         fu,& !< Upwind nodal point
         sum_cfl_out !< Sum of absolute value of out flowing CFL numbers

    do k=2, kmax
      kp1=min(k+1,local_grid%size(Z_INDEX))
      kp0=kp1-1
      call calculate_stencil_for_u(k, advneg, advpos, fc, fcurvs, fd, fdels, fu, &
           u%data(:, y_flow_index, x_flow_index), zf%data(:, y_scalar_index, x_scalar_index-1), &
           zf%data(:, y_scalar_index, x_scalar_index+2), zf%data(:, y_scalar_index, x_scalar_index), &
           zf%data(:, y_scalar_index, x_scalar_index+1))      
      if(abs(fcurvs) .ge. abs(fdels))then
        flux_x(k) = fc
      else        
        if (u%data(k, y_flow_index, x_flow_index) .gt. 0.0_DEFAULT_PRECISION) then
          sum_cfl_out=0.0_DEFAULT_PRECISION
#ifdef W_ACTIVE
          sum_cfl_out = sum_cfl_out + rdz(k)*(max(0.0_DEFAULT_PRECISION,&
               w%data(k, y_flow_index, x_flow_index))+abs(min(0.0_DEFAULT_PRECISION,w%data(k-1, y_flow_index, x_flow_index))))
#endif
#ifdef V_ACTIVE
          sum_cfl_out = sum_cfl_out + grid_config%horizontal%cy*(max(0.0_DEFAULT_PRECISION,&
               v%data(k, y_flow_index, x_flow_index))+abs(min(0.0_DEFAULT_PRECISION,v%data(k, y_flow_index-1, x_flow_index))))
#endif
#ifdef U_ACTIVE
          sum_cfl_out = sum_cfl_out + grid_config%horizontal%cx*(u%data(k, y_flow_index, x_flow_index)+&
               abs(min(0.0_DEFAULT_PRECISION,u%data(k, y_flow_index, x_flow_index-1))))
#endif
          sum_cfl_out = sum_cfl_out * dt          
        else
          sum_cfl_out=0.0_DEFAULT_PRECISION
#ifdef W_ACTIVE
          sum_cfl_out = sum_cfl_out + rdz(k)*&
               (max(0.0_DEFAULT_PRECISION,w%data(k, y_flow_index, x_flow_index+1))+&
               abs(min(0.0_DEFAULT_PRECISION,w%data(k-1, y_flow_index, x_flow_index+1))))
#endif
#ifdef V_ACTIVE
          sum_cfl_out = sum_cfl_out + grid_config%horizontal%cy*(max(0.0_DEFAULT_PRECISION,&
               v%data(k, y_flow_index, x_flow_index+1))+&
               abs(min(0.0_DEFAULT_PRECISION,v%data(k, y_flow_index-1, x_flow_index+1))))
#endif
#ifdef U_ACTIVE
          sum_cfl_out = sum_cfl_out + grid_config%horizontal%cx*(max(0.0_DEFAULT_PRECISION,&
               u%data(k, y_flow_index, x_flow_index+1))-u%data(k, y_flow_index, x_flow_index))
#endif
          sum_cfl_out = sum_cfl_out * dt
        end if
        ! First gradient term in Z direction
        fgt1=0.0_DEFAULT_PRECISION
#ifdef W_ACTIVE
        fgt1=calculate_gradient_term_in_z(k, w%data(:, y_flow_index, x_flow_index), &
             w%data(:, y_flow_index, x_flow_index+1), dt, advneg, advpos, kdof, kp0, kp1, &
             zf%data(:, y_scalar_index, x_scalar_index), zf%data(:, y_scalar_index, x_scalar_index+1), rdz)
#endif

        ! Second gradient term in Y direction
        fgt2=0.0_DEFAULT_PRECISION
#ifdef V_ACTIVE
        fgt2=calculate_gradient_term_in_y(k, v, y_flow_index, x_flow_index+1, &
             grid_config%horizontal%cy, dt, fc, zf, y_scalar_index, x_scalar_index, multiplication_factor_one=advpos, &
             term_data_two=zf, term_data_two_y_index=y_scalar_index, term_data_two_x_index=x_scalar_index+1, &
             multiplication_factor_two=advneg, data_two=v, data_two_y_index=y_flow_index, data_two_x_index=x_flow_index)
#endif
        ! Now calculate the fluxes in X
        flux_x(k) = calculate_flux_in_x(u%data(k, y_flow_index, x_flow_index), dt, sum_cfl_out, fdels, fd, fc, fu, fgt1, fgt2, &
             grid_config%horizontal%cx, fcurvs)
      end if
    end do
  end subroutine handle_x_direction_fluxes

  subroutine handle_y_direction_fluxes(y_flow_index, x_flow_index, u, v, w, y_scalar_index, x_scalar_index, zf, &
       local_grid, grid_config, kdof, dt, flux_y, rdz, kmin, kmax)
    integer, intent(in) :: y_flow_index, x_flow_index, y_scalar_index, x_scalar_index, kdof, kmin, kmax
    real(kind=DEFAULT_PRECISION), intent(in) :: dt
    type(grid_configuration_type), intent(inout) :: grid_config
    type(prognostic_field_type), intent(inout) :: u, w, v, zf
    type(local_grid_type), intent(inout) :: local_grid
    real(kind=DEFAULT_PRECISION), intent(in), dimension(:) :: rdz
    real(kind=DEFAULT_PRECISION), intent(inout), dimension(:) :: flux_y

    integer :: k, kp1, kp0, jofset, jperiod

    real(kind=DEFAULT_PRECISION) :: advneg,advpos,& !< Sign indicators (0 or 1)
         fgt1, fgt2,&  !< Gradient terms
         fc,& !< Upwinded nodal points in z
         fcurvs,& !< Curvature of F
         fd,& !< Downwind nodal points in z
         fdels,&  !< Downwind-upwind difference of F
         fu,& !< Upwind nodal point
         sum_cfl_out !< Sum of absolute value of out flowing CFL numbers

    do k=2, kmax
      kp1=min(k+1,local_grid%size(Z_INDEX))
      kp0=kp1-1
      ! Set up the stencil_mod
      call calculate_stencil_for_y(y_flow_index, x_flow_index, y_scalar_index, x_scalar_index, k, advneg, &
           advpos, fc, fcurvs, fd, fdels, fu, jofset, jperiod, v, zf)
      if (abs(fcurvs) .ge. abs(fdels) )then
        flux_y(k) = fc
      else
        sum_cfl_out=0.0_DEFAULT_PRECISION
#ifdef W_ACTIVE
        sum_cfl_out = sum_cfl_out + rdz(k)* (max(0.0_DEFAULT_PRECISION,w%data(k, y_flow_index+jofset, &
             x_flow_index))+ abs(min(0.0_DEFAULT_PRECISION,w%data(k-1, y_flow_index+jofset, x_flow_index))))
#endif
#ifdef V_ACTIVE
        sum_cfl_out = sum_cfl_out + grid_config%horizontal%cy*(max(0.0_DEFAULT_PRECISION,v%data(&
             k, y_flow_index+jofset, x_flow_index)) +&
             abs(min(0.0_DEFAULT_PRECISION,v%data(k, y_flow_index-1 +jofset+jperiod, x_flow_index))))
#endif
#ifdef U_ACTIVE
        sum_cfl_out = sum_cfl_out + grid_config%horizontal%cx*(max(0.0_DEFAULT_PRECISION,u%data(&
             k, y_flow_index+jofset, x_flow_index)) + &
             abs(min(0.0_DEFAULT_PRECISION,u%data(k, y_flow_index+jofset, x_flow_index-1))))
#endif
        sum_cfl_out = sum_cfl_out * dt
        ! First gradient term in Z direction
        fgt1=0.0_DEFAULT_PRECISION
#ifdef W_ACTIVE
        fgt1= calculate_gradient_term_in_z(k, w%data(:, y_flow_index, x_flow_index), w%data(:, y_flow_index+1, &
             x_flow_index), dt, advneg, advpos, kdof, kp0, kp1, &
             zf%data(:, y_scalar_index+jofset, x_scalar_index), rdz=rdz)
#endif
        ! Second gradient term in X direction
        fgt2=0.0_DEFAULT_PRECISION
#ifdef U_ACTIVE
        fgt2=calculate_gradient_term_in_x(k, u%data(:, y_flow_index, x_flow_index-1), u%data(:, y_flow_index, &
             x_flow_index), u%data(:, y_flow_index+1, x_flow_index-1), u%data(:, y_flow_index+1, x_flow_index), &
             grid_config%horizontal%cx, dt, fc, zf%data(k, y_scalar_index+jofset, x_scalar_index+1), &
             zf%data(k, y_scalar_index+jofset, x_scalar_index-1), 0)
#endif
        ! Now calculate the flux in Y
        flux_y(k)=calculate_flux_in_y(v%data(k, y_flow_index, x_flow_index), dt, sum_cfl_out, &
             fdels, fd, fc, fu, fgt1, fgt2, grid_config%horizontal%cy, fcurvs) 
      end if
    end do
  end subroutine handle_y_direction_fluxes

  subroutine handle_vertical_fluxes_middleofcolumn(y_flow_index, x_flow_index, u, v, w, y_scalar_index, x_scalar_index, zf, &
       local_grid, grid_config, kdof, dt, flux_z, rdz, rdzn, dzn)
    integer, intent(in) :: kdof, y_flow_index, x_flow_index, y_scalar_index, x_scalar_index
    real(kind=DEFAULT_PRECISION), intent(in) :: dt
    type(grid_configuration_type), intent(inout) :: grid_config
    type(prognostic_field_type), intent(inout) :: u, w, v, zf
    type(local_grid_type), intent(inout) :: local_grid
    real(kind=DEFAULT_PRECISION), intent(in), dimension(:) :: rdz, rdzn, dzn
    real(kind=DEFAULT_PRECISION), intent(inout), dimension(:) :: flux_z

    integer :: kofset, k

    real(kind=DEFAULT_PRECISION) :: advneg,advpos,& !< Sign indicators (0 or 1)
         fgt1, fgt2,&  !< Gradient terms
         fc,& !< Upwinded nodal points in z
         fcurvs,& !< Curvature of F
         fd,& !< Downwind nodal points in z
         fdels,&  !< Downwind-upwind difference of F
         fu,& !< Upwind nodal point
         rdc,&  !< Central 1/grid size
         rdu,&  !< Upwinded 1/grid size
         sum_cfl_out !< Sum of absolute value of out flowing CFL numbers

    do k=2,local_grid%size(Z_INDEX)-2
      call calculate_stencil_for_w(y_flow_index, x_flow_index, y_scalar_index, x_scalar_index, k, advneg, &
           advpos, fc, fcurvs, fd, fdels, fu, kofset, w, zf)

      if(abs(fcurvs) .ge. abs(fdels))then
        flux_z(k+1) = fc
      else        
        sum_cfl_out=0.0_DEFAULT_PRECISION
#ifdef W_ACTIVE
        sum_cfl_out=sum_cfl_out+rdz(k+kofset)*(max(0.0_DEFAULT_PRECISION,w%data(k+kofset, y_flow_index, &
             x_flow_index))+ abs(min(0.0_DEFAULT_PRECISION,w%data(k-1+kofset, y_flow_index, x_flow_index))))
#endif
#ifdef V_ACTIVE
        sum_cfl_out=sum_cfl_out+grid_config%horizontal%cy*(max(0.0_DEFAULT_PRECISION,&
             v%data(k+kofset, y_flow_index, x_flow_index))+&
             abs(min(0.0_DEFAULT_PRECISION,v%data(k+kofset, y_flow_index-1, x_flow_index))))
#endif
#ifdef U_ACTIVE
        sum_cfl_out=sum_cfl_out+grid_config%horizontal%cx*(max(0.0_DEFAULT_PRECISION,&
             u%data(k+kofset, y_flow_index, x_flow_index))+ &
             abs(min(0.0_DEFAULT_PRECISION,u%data(k+kofset, y_flow_index, x_flow_index-1))))
#endif
        sum_cfl_out = sum_cfl_out * dt

        ! First gradient term, in Y direction
        fgt1=0.0_DEFAULT_PRECISION
#ifdef V_ACTIVE
        fgt1=calculate_gradient_term_in_y(k, v, y_flow_index, x_flow_index, grid_config%horizontal%cy, &
             dt, fc, zf, y_scalar_index, x_scalar_index, kofset)
#endif
        ! Second gradient TERM, in X direction
        fgt2=0.0_DEFAULT_PRECISION
#ifdef U_ACTIVE
        fgt2=calculate_gradient_term_in_x(k, u%data(:, y_flow_index, x_flow_index-1), u%data(:, y_flow_index, &
             x_flow_index), u%data(:, y_flow_index, x_flow_index-1), u%data(:, y_flow_index, &
             x_flow_index), grid_config%horizontal%cx, dt, fc, zf%data(k+kofset, y_scalar_index, x_scalar_index+1), &
             zf%data(k+kofset, y_scalar_index, x_scalar_index-1), 1)
#endif
        rdu=rdzn(k) *advpos + rdzn(k+2)*advneg
        rdc=rdz(k+kdof)*advpos + rdz(k+1+kdof)*advneg

        ! Calculate the fluxes in W
        flux_z(k+1) = calculate_flux_in_w(k, w%data(k, y_flow_index, x_flow_index), &
             dt, sum_cfl_out, rdc, fdels, fd, fc, fu, rdu, fgt1, fgt2, dzn, rdzn)
      end if
    end do
  end subroutine handle_vertical_fluxes_middleofcolumn

  subroutine handle_vertical_fluxes_bottomofcolumn(y_flow_index, x_flow_index, u, v, w, y_scalar_index, x_scalar_index, zf, &
       grid_config, kdof, dt, flux_z, rdz, rdzn, dzn)
    integer, intent(in) :: y_flow_index, x_flow_index, y_scalar_index, x_scalar_index, kdof
    real(kind=DEFAULT_PRECISION), intent(in) :: dt
    type(grid_configuration_type), intent(inout) :: grid_config
    type(prognostic_field_type), intent(inout) :: u, w, v, zf
    real(kind=DEFAULT_PRECISION), intent(in), dimension(:) :: rdz, rdzn, dzn
    real(kind=DEFAULT_PRECISION), intent(out), dimension(:) :: flux_z

    integer :: k

    real(kind=DEFAULT_PRECISION) :: fgt1, fgt2,&  !< Gradient terms
         fc,& !< Upwinded nodal points in z
         fcurvs,& !< Curvature of F
         fd,& !< Downwind nodal points in z
         fdels,&  !< Downwind-upwind difference of F
         fu,& !< Upwind nodal point
         rdc,&  !< Central 1/grid size
         rdu,&  !< Upwinded 1/grid size
         sum_cfl_out !< Sum of absolute value of out flowing CFL numbers

    k=1
    if(w%data(k, y_flow_index, x_flow_index) .ge. 0.0_DEFAULT_PRECISION) then
      fu = zf%data(1, y_scalar_index, x_scalar_index)
      fc = zf%data(k, y_scalar_index, x_scalar_index)
      fd = zf%data(k+1, y_scalar_index, x_scalar_index)
    else
      fu = zf%data(k+2, y_scalar_index, x_scalar_index)
      fc = zf%data(k+1, y_scalar_index, x_scalar_index)
      fd = zf%data(k, y_scalar_index, x_scalar_index)
    end if
    fcurvs = fd - 2.0_DEFAULT_PRECISION*fc + fu
    fdels  = fd - fu
    if(abs(fcurvs) .ge. abs(fdels))then
      flux_z(k+1) = fc
    else
      if(w%data(k, y_flow_index, x_flow_index) .ge. 0.0_DEFAULT_PRECISION) then
        sum_cfl_out=0.0_DEFAULT_PRECISION
#ifdef W_ACTIVE
        sum_cfl_out = sum_cfl_out + w%data(k, y_flow_index, x_flow_index)*rdz(k+kdof)
#endif
#ifdef V_ACTIVE
        sum_cfl_out = sum_cfl_out + grid_config%horizontal%cy*(max(0.0_DEFAULT_PRECISION,v%data(k, y_flow_index, x_flow_index))&
             +abs(min(0.0_DEFAULT_PRECISION,v%data(k, y_flow_index-1, x_flow_index))))
#endif
#ifdef U_ACTIVE
        sum_cfl_out = sum_cfl_out + grid_config%horizontal%cx*(max(0.0_DEFAULT_PRECISION,u%data(k, y_flow_index, x_flow_index))&
             +abs(min(0.0_DEFAULT_PRECISION,u%data(k, y_flow_index, x_flow_index-1))))
#endif
        sum_cfl_out = sum_cfl_out * dt

        ! First gradient term in Y direction
        fgt1=0.0_DEFAULT_PRECISION
#ifdef V_ACTIVE
        fgt1=calculate_gradient_term_in_y(k, v, y_flow_index, x_flow_index, &
             grid_config%horizontal%cy, dt, fc, zf, y_scalar_index, x_scalar_index)
#endif
        ! Second gradient term in X direction
        fgt2=0.0_DEFAULT_PRECISION
#ifdef U_ACTIVE
        fgt2=calculate_gradient_term_in_x(k, u%data(:, y_flow_index, x_flow_index-1), u%data(:, y_flow_index, &
             x_flow_index), u%data(:, y_flow_index, x_flow_index-1), u%data(:, y_flow_index, &
             x_flow_index), grid_config%horizontal%cx, dt, fc, zf%data(k, y_scalar_index, x_scalar_index+1), &
             zf%data(k, y_scalar_index, x_scalar_index-1), 1)
#endif
        rdu=rdzn(2)  ! k map_typeped onto 2
        rdc=rdz(k+kdof)
      else
        sum_cfl_out=0.0_DEFAULT_PRECISION
#ifdef W_ACTIVE
        sum_cfl_out=sum_cfl_out+rdz(k+1+kdof)*(max(0.0_DEFAULT_PRECISION,w%data(k+1, y_flow_index, &
             x_flow_index)) -w%data(k, y_flow_index, x_flow_index))
#endif
#ifdef V_ACTIVE
        sum_cfl_out=sum_cfl_out+grid_config%horizontal%cy*(max(0.0_DEFAULT_PRECISION,v%data(k+1, y_flow_index, x_flow_index))&
             +abs(min(0.0_DEFAULT_PRECISION,v%data(k+1, y_flow_index-1, x_flow_index))))
#endif
#ifdef U_ACTIVE
        sum_cfl_out=sum_cfl_out + grid_config%horizontal%cx*(max(0.0_DEFAULT_PRECISION,u%data(k+1, y_flow_index, x_flow_index))&
             +abs(min(0.0_DEFAULT_PRECISION,u%data(k+1, y_flow_index, x_flow_index-1))))
#endif
        sum_cfl_out = sum_cfl_out * dt

        ! First gradient term in Y direction
        fgt1=0.0_DEFAULT_PRECISION
#ifdef V_ACTIVE
        fgt1=calculate_gradient_term_in_y(k, v, y_flow_index, x_flow_index, &
             grid_config%horizontal%cy, dt, fc, zf, y_scalar_index, x_scalar_index, 1)
#endif
        ! Second gradient term in X direction
        fgt2=0.0_DEFAULT_PRECISION
#ifdef U_ACTIVE
        fgt2=calculate_gradient_term_in_x(k, u%data(:, y_flow_index, x_flow_index-1), u%data(:, y_flow_index, &
             x_flow_index), u%data(:, y_flow_index, x_flow_index-1), u%data(:, y_flow_index, &
             x_flow_index), grid_config%horizontal%cx, dt, fc, zf%data(k+1, y_scalar_index, x_scalar_index+1), &
             zf%data(k+1, y_scalar_index, x_scalar_index-1), 1)
#endif
        rdu=rdzn(k+2)
        rdc=rdz(k+1+kdof)
      end if
      ! Calculate fluxes in W
      flux_z(k+1) = calculate_flux_in_w(k, w%data(k, y_flow_index, x_flow_index), &
           dt, sum_cfl_out, rdc, fdels, fd, fc, fu, rdu, fgt1, fgt2, dzn, rdzn)
    end if
  end subroutine handle_vertical_fluxes_bottomofcolumn

  subroutine handle_vertical_fluxes_topofcolumn(y_flow_index, x_flow_index, u, v, w, y_scalar_index, x_scalar_index, zf, &
       local_grid, grid_config, kdof, dt, flux_z, rdz, rdzn, dzn)
    integer, intent(in) :: y_flow_index, x_flow_index, y_scalar_index, x_scalar_index, kdof
    real(kind=DEFAULT_PRECISION), intent(in) :: dt
    type(grid_configuration_type), intent(inout) :: grid_config
    type(prognostic_field_type), intent(inout) :: u, w, v, zf
    type(local_grid_type), intent(inout) :: local_grid
    real(kind=DEFAULT_PRECISION), intent(in), dimension(:) :: rdz, rdzn, dzn
    real(kind=DEFAULT_PRECISION), intent(inout), dimension(:) :: flux_z

    integer :: k

    real(kind=DEFAULT_PRECISION) :: fgt1, fgt2,&  !< Gradient terms
         fc,& !< Upwinded nodal points in z
         fcurvs,& !< Curvature of F
         fd,& !< Downwind nodal points in z
         fdels,&  !< Downwind-upwind difference of F
         fu,& !< Upwind nodal point
         rdc,&  !< Central 1/grid size
         rdu,&  !< Upwinded 1/grid size
         sum_cfl_out !< Sum of absolute value of out flowing CFL numbers

    k=local_grid%size(Z_INDEX)-1

    if(w%data(k, y_flow_index, x_flow_index) .ge. 0.0_DEFAULT_PRECISION) then
      fu = zf%data(k-1, y_scalar_index, x_scalar_index)
      fc = zf%data(k, y_scalar_index, x_scalar_index)
      fd = zf%data(k+1, y_scalar_index, x_scalar_index)
    else
      ! Extrapolate linearly to level equidistant above upper boundary
      fu = 2.0_DEFAULT_PRECISION*zf%data(k+1, y_scalar_index, x_scalar_index) - zf%data(k, y_scalar_index, x_scalar_index)
      fc = zf%data(k+1, y_scalar_index, x_scalar_index)
      fd = zf%data(k, y_scalar_index, x_scalar_index)
    end if
    fcurvs = fd - 2.0_DEFAULT_PRECISION*fc + fu
    fdels  = fd - fu
    if(abs(fcurvs) .ge. abs(fdels))then
      flux_z(k+1) = fc
    else
      if(w%data(k, y_flow_index, x_flow_index) .ge. 0.0_DEFAULT_PRECISION) then
        sum_cfl_out=0.0_DEFAULT_PRECISION
#ifdef W_ACTIVE
        sum_cfl_out = sum_cfl_out + rdz(k+kdof)*(w%data(k, y_flow_index, x_flow_index)&
             +abs(min(0.0_DEFAULT_PRECISION,w%data(k-1, y_flow_index, x_flow_index))))
#endif
#ifdef V_ACTIVE
        sum_cfl_out = sum_cfl_out + grid_config%horizontal%cy*(max(0.0_DEFAULT_PRECISION,v%data(k, y_flow_index, x_flow_index))+&
             abs(min(0.0_DEFAULT_PRECISION,v%data(k, y_flow_index-1, x_flow_index))))
#endif
#ifdef U_ACTIVE
        sum_cfl_out = sum_cfl_out + grid_config%horizontal%cx*(max(0.0_DEFAULT_PRECISION,u%data(k, y_flow_index, x_flow_index))&
             +abs(min(0.0_DEFAULT_PRECISION,u%data(k, y_flow_index, x_flow_index-1))))
#endif
        sum_cfl_out = sum_cfl_out * dt

        ! First gradient term in Y direction
        fgt1=0.0_DEFAULT_PRECISION
#ifdef V_ACTIVE
        fgt1=calculate_gradient_term_in_y(k, v, y_flow_index, x_flow_index, grid_config%horizontal%cy, dt, fc, &
             zf, y_scalar_index, x_scalar_index)
#endif
        ! Second gradient term in X direction
        fgt2=0.0_DEFAULT_PRECISION
#ifdef U_ACTIVE
        fgt2=calculate_gradient_term_in_x(k, u%data(:, y_flow_index, x_flow_index-1), u%data(:, y_flow_index, &
             x_flow_index), u%data(:, y_flow_index, x_flow_index-1), u%data(:, y_flow_index, &
             x_flow_index), grid_config%horizontal%cx, dt, fc, zf%data(k, y_scalar_index, x_scalar_index+1), &
             zf%data(k, y_scalar_index, x_scalar_index-1), 1)
#endif
        rdu=rdzn(k)
        rdc=rdz(k+kdof)
      else
        sum_cfl_out=0.0_DEFAULT_PRECISION
#ifdef W_ACTIVE
        sum_cfl_out = sum_cfl_out + (- w%data(k, y_flow_index, x_flow_index)&
             *rdz(min(local_grid%size(Z_INDEX),k+1+kdof)))
#endif
#ifdef V_ACTIVE
        sum_cfl_out = sum_cfl_out + grid_config%horizontal%cy*(max(0.0_DEFAULT_PRECISION,&
             v%data(k+1, y_flow_index, x_flow_index))+&
             abs(min(0.0_DEFAULT_PRECISION,v%data(k+1, y_flow_index-1, x_flow_index))))
#endif
#ifdef U_ACTIVE
        sum_cfl_out = sum_cfl_out + grid_config%horizontal%cx*(max(0.0_DEFAULT_PRECISION,&
             u%data(k+1, y_flow_index, x_flow_index))&
             +abs(min(0.0_DEFAULT_PRECISION,u%data(k+1, y_flow_index, x_flow_index-1))))
#endif
        sum_cfl_out = sum_cfl_out * dt

        ! First gradient term in Y direction
        fgt1=0.0_DEFAULT_PRECISION
#ifdef V_ACTIVE
        fgt1=calculate_gradient_term_in_y(k, v, y_flow_index, x_flow_index, grid_config%horizontal%cy, &
             dt, fc, zf, y_scalar_index, x_scalar_index, 1)
#endif
        ! Second gradient term in X direction
        fgt2=0.0_DEFAULT_PRECISION
#ifdef U_ACTIVE
        fgt2=calculate_gradient_term_in_x(k, u%data(:, y_flow_index, x_flow_index-1), u%data(:, y_flow_index, &
             x_flow_index), u%data(:, y_flow_index, x_flow_index-1), u%data(:, y_flow_index, &
             x_flow_index), grid_config%horizontal%cx, dt, fc, zf%data(k+1, y_scalar_index, x_scalar_index+1), &
             zf%data(k+1, y_scalar_index, x_scalar_index-1), 1)
#endif
        rdu=rdzn(local_grid%size(Z_INDEX))
        rdc=rdz(min(local_grid%size(Z_INDEX),k+1+kdof))
      end if
      ! Now calculate fluxes in W
      flux_z(k+1) = calculate_flux_in_w(k, w%data(k, y_flow_index, x_flow_index), &
           dt, sum_cfl_out, rdc, fdels, fd, fc, fu, rdu, fgt1, fgt2, dzn, rdzn)
    end if
  end subroutine handle_vertical_fluxes_topofcolumn

  real(kind=DEFAULT_PRECISION) function calculate_flux_in_w(k, data_value, dt, sum_cfl_out, rdc, fdels, fd, fc, fu, rdu, fgt1, &
       fgt2, dzn, rdzn)
    integer, intent(in) :: k
    real(kind=DEFAULT_PRECISION), intent(in), dimension(:) :: dzn, rdzn
    real(kind=DEFAULT_PRECISION), intent(in) :: dt, sum_cfl_out, rdc, fdels, fd, fc, fu, rdu, fgt1, fgt2
    real(kind=DEFAULT_PRECISION), intent(in) :: data_value

    real(kind=DEFAULT_PRECISION) :: dd,&         !< Downwind grid spacing
         rdd,&        !< Downwinded 1/grid size
         cflmod,&     !< Absolute value of face CFL number
         hrcfl,&      !< 1/SUM_CFL_OUT
         hcfl,&       !< 0.5*CFLMOD
         cflcoef,&    !< Coefficient in the QUICKEST calculation
         rfdels,&     !< 1/FZDELS
         ffs,&        !< 1d QUICKEST face value of F
         fnfs,&       !< Normalised face value
         fncs,&       !< Normalised central node value
         fnrefs,&     !< Normalised reference value
         ftemp        !< Limited face value

    dd=dzn(k+1)
    rdd=rdzn(k+1)
    cflmod=abs(data_value*rdzn(k+1)*dt)
    hrcfl=1.0_DEFAULT_PRECISION/(sum_cfl_out+1.e-30_DEFAULT_PRECISION)
    hcfl=0.5_DEFAULT_PRECISION*cflmod
    cflcoef=r6*(1.0_DEFAULT_PRECISION-cflmod*cflmod)*(dd*dd*rdc)
    rfdels=1.0_DEFAULT_PRECISION/(fdels+1.e-30_DEFAULT_PRECISION)
    ffs=0.5_DEFAULT_PRECISION*(fd+fc)-hcfl*(fd-fc)-cflcoef*((fd-fc)*rdd-(fc-fu)*rdu) - fgt1 - fgt2
    fnfs=(ffs-fu)*rfdels
    fncs=(fc-fu)*rfdels
    fnrefs=fncs*hrcfl
    ftemp=max(fncs,min(fnfs,fnrefs,1.0_DEFAULT_PRECISION))
    calculate_flux_in_w = ftemp*fdels + fu
  end function calculate_flux_in_w

  real(kind=DEFAULT_PRECISION) function calculate_flux_in_x(data_value, dt, sum_cfl_out, fdels, fd, fc, fu, fgt1, fgt2, cx, fcurvs)
    real(kind=DEFAULT_PRECISION), intent(in) :: sum_cfl_out, fdels, fd, fc, fu, fgt1, fgt2, cx, dt, fcurvs
    real(kind=DEFAULT_PRECISION), intent(in) :: data_value

    real(kind=DEFAULT_PRECISION) :: cflmod,&     !< Absolute value of face CFL number
         hrcfl,&      !< 1/SUM_CFL_OUT
         hcfl,&       !< 0.5*CFLMOD
         cflcoef,&    !< Coefficient in the QUICKEST calculation
         rfdels,&     !< 1/FZDELS
         ffs,&        !< 1d QUICKEST face value of F
         fnfs,&       !< Normalised face value
         fncs,&       !< Normalised central node value
         fnrefs,&     !< Normalised reference value
         ftemp        !< Limited face value

    cflmod=abs(data_value*cx*dt)
    hrcfl=1.0_DEFAULT_PRECISION/(sum_cfl_out+1.e-30_DEFAULT_PRECISION)
    hcfl=0.5_DEFAULT_PRECISION*cflmod
    cflcoef=r6*(1.0_DEFAULT_PRECISION-cflmod*cflmod)
    rfdels=1.0_DEFAULT_PRECISION/(fdels+1.e-30_DEFAULT_PRECISION)
    ffs=0.5_DEFAULT_PRECISION*(fd+fc)-hcfl*(fd-fc)-cflcoef*fcurvs - fgt1 - fgt2
    fnfs=(ffs-fu)*rfdels
    fncs=(fc-fu)*rfdels
    fnrefs=fncs*hrcfl
    ftemp=max(fncs,min(fnfs,fnrefs,1.0_DEFAULT_PRECISION))
    calculate_flux_in_x = ftemp*fdels + fu    
  end function calculate_flux_in_x

  real(kind=DEFAULT_PRECISION) function calculate_flux_in_y(data_value, dt, sum_cfl_out, fdels, fd, fc, fu, fgt1, fgt2, cy, fcurvs)
    real(kind=DEFAULT_PRECISION), intent(in) :: dt, sum_cfl_out, fdels, fd, fc, fu, fgt1, fgt2, cy, fcurvs
    real(kind=DEFAULT_PRECISION), intent(in) :: data_value

    real(kind=DEFAULT_PRECISION) :: cflmod,&     !< Absolute value of face CFL number
         hrcfl,&      !< 1/SUM_CFL_OUT
         hcfl,&       !< 0.5*CFLMOD
         cflcoef,&    !< Coefficient in the QUICKEST calculation
         rfdels,&     !< 1/FZDELS
         ffs,&        !< 1d QUICKEST face value of F
         fnfs,&       !< Normalised face value
         fncs,&       !< Normalised central node value
         fnrefs,&     !< Normalised reference value
         ftemp        !< Limited face value

    cflmod=abs(data_value*cy*dt)
    hrcfl=1.0_DEFAULT_PRECISION/(sum_cfl_out+1.e-30_DEFAULT_PRECISION)
    hcfl=0.5_DEFAULT_PRECISION*cflmod
    cflcoef=r6*(1.0_DEFAULT_PRECISION-cflmod*cflmod)
    rfdels=1.0_DEFAULT_PRECISION/(fdels+1.e-30_DEFAULT_PRECISION)
    ffs=0.5_DEFAULT_PRECISION*(fd+fc)-hcfl*(fd-fc)-cflcoef*fcurvs - fgt1 - fgt2
    fnfs=(ffs-fu)*rfdels
    fncs=(fc-fu)*rfdels
    fnrefs=fncs*hrcfl
    ftemp=max(fncs,min(fnfs,fnrefs,1.0_DEFAULT_PRECISION))
    calculate_flux_in_y = ftemp*fdels + fu
  end function calculate_flux_in_y

  subroutine calculate_stencil_for_u(k, advneg, advpos, fc, fcurvs, fd, fdels, fu, &
       flow_field, advection_field_one, advection_field_two, advection_field_three, advection_field_four)
    integer, intent(in) :: k
    real(kind=DEFAULT_PRECISION), intent(in), dimension(:) :: flow_field, advection_field_one, advection_field_two, &
         advection_field_three, advection_field_four
    real(kind=DEFAULT_PRECISION), intent(out) :: advneg, advpos, fc, fcurvs, fd, fdels, fu

    real(kind=DEFAULT_PRECISION) :: sign_modified !< 0.5*sign of u, v or w

    sign_modified = sign(real(.5, kind=DEFAULT_PRECISION),flow_field(k))
    advpos = 0.5_DEFAULT_PRECISION+sign_modified
    advneg = 0.5_DEFAULT_PRECISION-sign_modified
    fu = advection_field_one(k)*advpos + advection_field_two(k)*advneg
    fc = advection_field_three(k)*advpos + advection_field_four(k)*advneg
    fd = advection_field_four(k)*advpos + advection_field_three(k)*advneg
    fdels = fd - fu
    fcurvs = fd - 2.0_DEFAULT_PRECISION*fc + fu
  end subroutine calculate_stencil_for_u

  subroutine calculate_stencil_for_y(y_flow_index, x_flow_index, y_scalar_index, x_scalar_index, k, advneg, advpos, &
       fc, fcurvs, fd, fdels, fu, jofset, jperiod, flow_field, zf)
    integer, intent(in) :: y_flow_index, x_flow_index, y_scalar_index, x_scalar_index, k
    type(prognostic_field_type), intent(inout) :: flow_field, zf
    integer, intent(out) :: jofset, jperiod
    real(kind=DEFAULT_PRECISION), intent(out) :: advneg, advpos, fc, fcurvs, fd, fdels, fu

    real(kind=DEFAULT_PRECISION) :: ajeq0, sign_modified !< 0.5*sign of u, v or w

    sign_modified=sign(0.5_DEFAULT_PRECISION, flow_field%data(k, y_flow_index, x_flow_index))
    advpos=0.5_DEFAULT_PRECISION+sign_modified
    advneg=0.5_DEFAULT_PRECISION-sign_modified
    jofset=nint(advneg)
    ajeq0=0.5_DEFAULT_PRECISION-sign(0.5_DEFAULT_PRECISION, (y_scalar_index-1)-0.5_DEFAULT_PRECISION)
    jperiod=merge(0,1, y_scalar_index .gt. 1) !nint(jjp*ajeq0)  ! =JJP OR 0

    fu = zf%data(k, y_scalar_index-1+jperiod, x_scalar_index)*advpos + zf%data(k, y_scalar_index+2, x_scalar_index)*advneg
    fc = zf%data(k, y_scalar_index, x_scalar_index)  *advpos + zf%data(k, y_scalar_index+1, x_scalar_index)*advneg
    fd = zf%data(k, y_scalar_index+1, x_scalar_index)*advpos + zf%data(k, y_scalar_index, x_scalar_index)*advneg
    fdels = fd - fu
    fcurvs = fd - 2.0_DEFAULT_PRECISION*fc + fu
  end subroutine calculate_stencil_for_y

  subroutine calculate_stencil_for_w(y_flow_index, x_flow_index, y_scalar_index, x_scalar_index, k, advneg, advpos, &
       fc, fcurvs, fd, fdels, fu, kofset, flow_field, zf)
    integer, intent(in) :: y_flow_index, x_flow_index, y_scalar_index, x_scalar_index, k
    type(prognostic_field_type), intent(inout) :: flow_field, zf
    integer, intent(out) :: kofset
    real(kind=DEFAULT_PRECISION), intent(out) :: advneg, advpos, fc, fcurvs, fd, fdels, fu

    real(kind=DEFAULT_PRECISION) :: sign_modified !< 0.5*sign of u, v or w

    sign_modified=sign(real(.5, kind=DEFAULT_PRECISION),flow_field%data(k, y_flow_index, x_flow_index))
    advpos=0.5_DEFAULT_PRECISION+sign_modified
    advneg=0.5_DEFAULT_PRECISION-sign_modified
    kofset=nint(advneg)
    fu = zf%data(k-1, y_scalar_index, x_scalar_index)*advpos + zf%data(k+2, y_scalar_index, x_scalar_index)*advneg
    fc = zf%data(k+kofset, y_scalar_index, x_scalar_index)
    fd = zf%data(k+1, y_scalar_index, x_scalar_index)*advpos + zf%data(k, y_scalar_index, x_scalar_index)  *advneg
    fcurvs = fd - 2.0_DEFAULT_PRECISION*fc + fu
    fdels  = fd - fu
  end subroutine calculate_stencil_for_w

  real(kind=DEFAULT_PRECISION) function calculate_gradient_term_in_z(k, column_one, column_two, dt, advneg, &
       advpos, kdof, kp0, kp1, source_column_one, source_column_two, rdz)
    integer, intent(in) :: k, kdof, kp0, kp1
    real(kind=DEFAULT_PRECISION), intent(in) :: dt, advneg, advpos
    real(kind=DEFAULT_PRECISION), intent(in), dimension(:) :: column_one, column_two, source_column_one
    real(kind=DEFAULT_PRECISION), intent(in), dimension(:), optional :: source_column_two
    real(kind=DEFAULT_PRECISION), intent(in), dimension(:) :: rdz

    integer :: ksu, ksl, kposwf
    real(kind=DEFAULT_PRECISION) :: specific_face, temp_src

    specific_face=0.125_DEFAULT_PRECISION*(column_one(k)+column_one(k-1)+column_two(k)+column_two(k-1))
    kposwf=nint(0.5_DEFAULT_PRECISION+sign(real(0.5, kind=DEFAULT_PRECISION), specific_face))
    ksu=kp1*(1-kposwf)+(k-1)*kposwf
    ksl=kp0*(1-kposwf)+k*kposwf
    temp_src=source_column_one(ksl)-source_column_one(ksu)
    if (present(source_column_two)) temp_src = temp_src * advpos+(source_column_two(ksl)-source_column_two(ksu))*advneg
    calculate_gradient_term_in_z=abs(specific_face)*rdz(k+kdof)*dt*temp_src
  end function calculate_gradient_term_in_z

  real(kind=DEFAULT_PRECISION) function calculate_gradient_term_in_x(k, column_one, column_two, column_three, column_four, &
       cx, dt, fc, source_value_one, source_value_two, kpv)
    integer, intent(in) :: k, kpv
    real(kind=DEFAULT_PRECISION), intent(in) :: cx, dt, fc
    real(kind=DEFAULT_PRECISION), intent(in) :: source_value_one, source_value_two
    real(kind=DEFAULT_PRECISION), intent(in), dimension(:) :: column_one, column_two, column_three, column_four

    real(kind=DEFAULT_PRECISION) :: specific_face, posuf

    specific_face=0.125_DEFAULT_PRECISION*(column_one(k)+column_two(k)+column_three(k+kpv)+column_four(k+kpv))
    posuf=0.5_DEFAULT_PRECISION+sign(real(.5, kind=DEFAULT_PRECISION), specific_face)
    calculate_gradient_term_in_x=ABS(specific_face)*cx*dt*(fc-(source_value_one*(1.-posuf)+source_value_two*posuf))
  end function calculate_gradient_term_in_x

  real(kind=DEFAULT_PRECISION) function calculate_gradient_term_in_y(k, data, data_y_index, data_x_index, cy, &
       dt, fc, term_data_one, term_data_one_y_index, term_data_one_x_index, koffset_one, multiplication_factor_one, &
       term_data_two, term_data_two_y_index, term_data_two_x_index, &
       koffset_two, multiplication_factor_two, data_two, data_two_y_index, data_two_x_index)

    integer, intent(in) :: k, data_y_index, data_x_index, term_data_one_y_index, term_data_one_x_index
    integer, intent(in), optional :: term_data_two_y_index, term_data_two_x_index, data_two_y_index, data_two_x_index
    integer, intent(in), optional :: koffset_one, koffset_two
    real(kind=DEFAULT_PRECISION), intent(in), optional :: multiplication_factor_one, multiplication_factor_two
    real(kind=DEFAULT_PRECISION), intent(in) :: cy, dt, fc
    type(prognostic_field_type), intent(in) :: data, term_data_one
    type(prognostic_field_type), intent(in), optional :: term_data_two, data_two

    real(kind=DEFAULT_PRECISION) :: specific_face, calc_term_one, calc_term_two
    integer :: jsu, k_index

    calc_term_two=0.

    if (present(data_two)) then
      specific_face=0.125_DEFAULT_PRECISION*(data%data(k, data_y_index, data_x_index)+&
           data_two%data(k, data_two_y_index, data_two_x_index)+&
           data%data(k, data_y_index-1, data_x_index)+ data_two%data(k, data_two_y_index-1, data_two_x_index))
    else
      specific_face=0.125_DEFAULT_PRECISION*(data%data(k, data_y_index-1, data_x_index)+&
           data%data(k, data_y_index, data_x_index)+&
           data%data(k+1, data_y_index-1, data_x_index)+data%data(k+1, data_y_index, data_x_index))
    end if    

    if (present(koffset_one)) then
      k_index=k+koffset_one
    else
      k_index=k
    end if

    jsu=term_data_one_y_index-nint(sign(real(1.0, kind=DEFAULT_PRECISION), specific_face))
    calc_term_one = term_data_one%data(k_index, jsu, term_data_one_x_index)
    if (present(multiplication_factor_one)) calc_term_one = calc_term_one * multiplication_factor_one
    if (present(term_data_two)) then
      if (present(koffset_two)) then
        k_index=k+koffset_two
      else
        k_index=k
      end if
      jsu=term_data_two_y_index-nint(sign(real(1.0, kind=DEFAULT_PRECISION), specific_face))
      calc_term_two = term_data_two%data(k_index, jsu, term_data_two_x_index)
      if (present(multiplication_factor_two)) calc_term_two = calc_term_two * multiplication_factor_two
      calc_term_one = calc_term_one + calc_term_two
    end if
    calculate_gradient_term_in_y=abs(specific_face)*cy*dt*(fc-calc_term_one)
  end function calculate_gradient_term_in_y
end module ultimateflux_mod
