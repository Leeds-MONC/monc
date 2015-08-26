module maths_mod
  use datadefn_mod, only : DEFAULT_PRECISION
  implicit none

#ifndef TEST_MODE
  private
#endif

  integer, parameter :: k4b=selected_int_kind(9)
  integer(k4b), parameter :: ia=16807,im=2147483647
  integer(k4b), parameter :: iq=127773,ir=2836

  public random
contains

  !> returns a scalar random number, the initial seed idum must be negative
  !! usage: idum = -k !(set idum < 0 to initialise)
  !! do i=1,nvals
  !!   my_rand(i) = random(idum)
  !! end do
  !! @param idum Initial seed which must be negative
  real(kind=DEFAULT_PRECISION) function random(idum)    
    integer(k4b),intent(inout) :: idum
        
    real,save :: am
    integer(k4b), save :: ix=-1,iy=-1,k

    if (idum <=0 .or. iy < 0) then
      am = nearest(1.0,-1.0)/im
      iy=ior(ieor(888889999,abs(idum)),1)
      ix=ieor(777755555,abs(idum))
      idum=abs(idum)+1
    end if
    ix = ieor(ix,ishft(ix,13))
    ix = ieor(ix,ishft(ix,-17))
    ix = ieor(ix,ishft(ix,5))
    k=iy/iq
    iy=ia*(iy-k*iq)-ir*k
    if(iy < 0) iy = iy + im
    random = am*ior(iand(im,ieor(ix,iy)),1)
  end function random
end module maths_mod
