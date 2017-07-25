!< Utility interpolation routines
module interpolation_mod
  use datadefn_mod, only : DEFAULT_PRECISION

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

!!$    do k=1, size(field)
!!$      do nn=2,nnodes
!!$        if (zgrid(k) < zvals(nn))then
!!$          field(k) = vals(nn-1) + (vals(nn) - vals(nn-1))/(zvals(nn) - zvals(nn-1)) &
!!$             *(zgrid(k) - zvals(nn-1))
!!$          exit
!!$        end if
!!$      end do
!!$    end do

    ! now repeat the code of the LEM as a check
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

    !do k=1, size(field)
    !   print *, 'field_monc, field_lem =', field(k), field_lem(k), k
    !enddo

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

end module interpolation_mod
