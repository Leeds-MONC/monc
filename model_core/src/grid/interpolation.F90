!< Utility interpolation routines
module interpolation_mod
  use datadefn_mod, only : DEFAULT_PRECISION
  use logging_mod, only : LOG_ERROR, log_master_log, LOG_DEBUG, log_get_logging_level, log_log

  implicit none

contains

  !> Does a simple 1d piecewise linear interpolation
  !! @param zvals input z nodes
  !! @param vals  input nodal values
  !! @param zgrid grid to interpolate onto
  !! @param field output interpolated values
  subroutine piecewise_linear_1d(zvals, vals, zgrid, field)

    real(kind=DEFAULT_PRECISION), intent(in) :: zvals(:), vals(:)
    real(kind=DEFAULT_PRECISION), intent(inout) :: zgrid(:)
    real(kind=DEFAULT_PRECISION), intent(inout) :: field(:)

    real(kind=DEFAULT_PRECISION) :: initgd_lem(size(zvals)+1)
    real(kind=DEFAULT_PRECISION) :: zngd_lem(size(zvals)+1)
    real(kind=DEFAULT_PRECISION) :: field_lem(size(field))

    real(kind=DEFAULT_PRECISION) :: verylow, veryhigh
    
    integer :: nn,k  ! loop counters
    integer :: nnodes    ! number of input values

    nnodes=size(zvals)


    ! Code replicated from the LEM. This duplicates the interpolation of the 
    ! LEM exactly
    verylow = -1.0e5          
    
    initgd_lem(1) = vals(1) 
    zngd_lem(1) = verylow 
    DO nn=2,nnodes+1              
       initgd_lem(nn) = vals(nn-1)
       zngd_lem(nn) = zvals(nn-1) 
    ENDDO

    do nn=1,nnodes
       DO k=1,size(field)
          IF( zngd_lem(nn).LE.zgrid(k) .AND. zgrid(k).LT.zngd_lem(nn+1)) then 
             field(k) = initgd_lem(nn) + ( initgd_lem(nn+1) - initgd_lem(nn))*(zgrid(k)- zngd_lem(nn))/ &
                  (zngd_lem(nn+1) - zngd_lem(nn))
          end IF
       end DO
    end do

    ! Add check in for when zgrid(k_top) equals zngd_lem(nn+1), which is the case for subsidence
    if (zgrid(size(field)) == zngd_lem(nnodes+1)) Then
       field(size(field)) = initgd_lem(nnodes+1)
    endif
    
  end subroutine piecewise_linear_1d

  !> Does a simple 1d linear interpolation to a point
  !! @param zvals input z nodes
  !! @param vals  input nodal values
  !! @param z location to interpolate onto
  !! @param f output interpolated value
  subroutine interpolate_point_linear_1d(zvals, vals, z, f, extrapolate)

    real(kind=DEFAULT_PRECISION), intent(in) :: zvals(:), vals(:)
    real(kind=DEFAULT_PRECISION), intent(in) :: z
    real(kind=DEFAULT_PRECISION), intent(out) :: f
    character(*), intent(in), optional :: extrapolate
    
    integer :: nn  ! loop counter
    integer :: nnodes    ! number of input values
    
    integer, parameter :: MAXCHARS=20
    character(MAXCHARS) :: ext_type

    nnodes=size(zvals)
    ! suggested fix from JME
    do
      if (nnodes == 1) exit
      if (zvals(nnodes-1) < zvals(nnodes)) exit
      nnodes = nnodes - 1
    enddo
    !
    ext_type='linear'
    if (present(extrapolate))ext_type=trim(extrapolate)

    if (z < zvals(1))then 
      nn=1
      select case (trim(ext_type))
      case ('linear')
        f = vals(nn) + (vals(nn+1) - vals(nn))/(zvals(nn+1) - zvals(nn)) &
           *(z - zvals(nn))
      case ('constant')
        f = vals(nn)
      end select
      return
    end if
      
    if (present(extrapolate))ext_type=trim(extrapolate)

    if (z >= zvals(nnodes))then 
      nn=nnodes
      select case (trim(ext_type))
      case ('linear')
        f = vals(nn) + (vals(nn-1) - vals(nn))/(zvals(nn-1) - zvals(nn)) &
           *(z - zvals(nn))
      case ('constant')
        f = vals(nn)
      end select
      return
    end if

    do nn=1,nnodes-1
      if (zvals(nn) <= z .and. z < zvals(nn+1))then
        f = vals(nn) + (vals(nn+1) - vals(nn))/(zvals(nn+1) - zvals(nn)) &
           *(z - zvals(nn))
        exit
      end if
    end do

  end subroutine interpolate_point_linear_1d

  !> Does a simple 1d linear interpolation to a point
  !! @param zvals_in input z nodes
  !! @param vals_in input nodal values
  !! @param z_out location to interpolate onto
  !! @param f output interpolated value
  subroutine piecewise_linear_2d(zvals_in, time_vals, vals_in, z_out, field)

    ! Assumes input variables (vals) are 2-D, with dims (z, time) 

    real(kind=DEFAULT_PRECISION), intent(in) :: zvals_in(:), time_vals(:)
    real(kind=DEFAULT_PRECISION), intent(in) :: vals_in(:,:)
    real(kind=DEFAULT_PRECISION), intent(in) :: z_out(:)
    real(kind=DEFAULT_PRECISION), intent(out) :: field(:,:)

    real(kind=DEFAULT_PRECISION) :: scale_tmp
    
    integer :: nn, k_monc, k_force                     ! loop counter
    integer :: nz_force, nt_force, nz_monc, nt_monc    ! time and height array sizes for forcing and monc grids
    integer :: nnodes                                  ! number of input values
    real(kind=DEFAULT_PRECISION), dimension(:), allocatable :: zvals, z
    real(kind=DEFAULT_PRECISION), dimension(:,:), allocatable :: vals
    

    nz_force = size(zvals_in)
    nt_force = size(time_vals)
    nz_monc  = size(z_out)
    nt_monc  = size(time_vals) ! time is intepolated in the timestep callback

    allocate(zvals(nz_force),z(nz_monc),vals(nz_force,nt_force))

    zvals=zvals_in

    if ( zvals(1) .GT. zvals(nz_force) ) then   ! pressure
      zvals=log10(zvals_in(nz_force:1:-1))
      z=log10(z_out(nz_monc:1:-1))
      vals=vals_in(nz_force:1:-1,:)
    else
      zvals=zvals_in
      z=z_out
      vals=vals_in
    end if

       do k_monc=1,nz_monc                                                     
          do k_force=1,nz_force-1                                          
             if( z(k_monc) >= zvals(k_force) .AND. z(k_monc) < zvals(k_force+1) ) then                  
                scale_tmp = ( z(k_monc) - zvals(k_force) ) /              &
                     ( zvals(k_force+1) - zvals(k_force) )               
                do nn=1, nt_force                                       
                   field(k_monc,nn) = vals(k_force,nn) +                  &           
                        (  vals(k_force+1,nn) - vals(k_force,nn) )        &          
                        * scale_tmp                           
                enddo
             endif
          enddo
       enddo
       ! now examine the cases below and above forlevs(1) and forlevs(ktmfor
       ! uses the local vertical gradient in the forcing to determine the   
       ! new values                                                         
       do k_monc=1,nz_monc                                                
          if ( z(k_monc) >= zvals(nz_force) ) then                    
             scale_tmp = ( z(k_monc) - zvals(nz_force) )                   &                
                  / ( zvals(nz_force) - zvals(nz_force-1) )           
             do nn=1,nt_force                                            
                field(k_monc,nn) = vals(nz_force,nn) +                  &           
                     (  vals(nz_force,nn) - vals(nz_force-1,nn) )        &          
                     * scale_tmp   
             enddo
          elseif ( z(k_monc) < zvals(1) )THEN                     
             scale_tmp = ( z(k_monc) - zvals(1) )                        &                      
                  / ( zvals(1) - zvals(2) )                       
             do nn=1,nt_force  
                field(k_monc,nn) = vals(1,nn) +                  &           
                     (  vals(1,nn) - vals(2,nn) )        &          
                     * scale_tmp   
             enddo
          endif
       enddo
       !                                                                    
    if ( zvals(nz_force) .GT. zvals(1) ) then   ! pressure (flipped coordinates)
      field=field(nz_monc:1:-1,:)
    endif

  end subroutine piecewise_linear_2d



  ! "extrapolate" is a bad name for this option, as extrapolate is only used to determine how to 
  ! do the interpolation so that it is either linear or constant (not actually interpolation)
  ! Actual extrapolation is only constant and automatic, it's just constantly replicated. Oy.
  ! Pleasantly, the default is to do internal linear interpolation.
  !  See interpolate_point_linear_1d for what looks to be the intended implementation.
  subroutine interpolate_point_linear_2d(zvals, vals, z, f, extrapolate)

    ! 2-d because the "vals" array is 2d, probably height and time

    real(kind=DEFAULT_PRECISION), intent(in) :: zvals(:), vals(:,:)  ! forcing_times(:),forcing_values(height,time)
    real(kind=DEFAULT_PRECISION), intent(in) :: z                    ! time to interpolate to
    real(kind=DEFAULT_PRECISION), intent(out) :: f(:)                ! output profile(height)
    character(*), intent(in), optional :: extrapolate
    
    integer :: nn  ! loop counter
    integer :: nnodes    ! number of input values
    
    integer, parameter :: MAXCHARS=20
    character(MAXCHARS) :: ext_type

    nnodes=size(zvals)
    ! suggested fix from JME
    do
      if (nnodes == 1) exit
      if (zvals(nnodes-1) < zvals(nnodes)) exit
      nnodes = nnodes - 1
    enddo
    !
    ext_type='linear'
    if (present(extrapolate))ext_type=trim(extrapolate)
    
    ! test where model time is outside the bounds of the input forcing times

    ! 1) less than the lowest time level in forcing times
    if (z < zvals(1))then
       nn=1
      f(:) = vals(:,nn)
      return
    end if
      
    !if (present(extrapolate))ext_type=trim(extrapolate)

    ! 2) greater than the last time level in the forcing times, make forcing constant
    if (z >= zvals(nnodes))then 
       nn=nnodes
       f(:) = vals(:,nn)
       return
    end if

    ! Time is within the bounds of the forcing input time

    do nn = 1, nnodes-1       
       if (zvals(nn) <= z .and. z < zvals(nn+1)) then
          select case (trim(ext_type))
          case ('linear')
             f(:) = vals(:,nn) + (vals(:,nn+1) - vals(:,nn))/(zvals(nn+1) - zvals(nn)) &
                  *(z - zvals(nn))
          case ('constant')
             f(:) = vals(:,nn)
          end select
          exit
       endif
    enddo

    
  end subroutine interpolate_point_linear_2d

end module interpolation_mod
