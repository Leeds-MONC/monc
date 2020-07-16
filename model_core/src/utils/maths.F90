module maths_mod
  use datadefn_mod, only : DEFAULT_PRECISION
  implicit none

#ifndef TEST_MODE
  private
#endif

  integer, parameter :: k4b=selected_int_kind(9)
  integer(k4b), parameter :: ia=16807,im=2147483647
  integer(k4b), parameter :: iq=127773,ir=2836

  public random, sort_1d
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


  !> Combines with MergeSortMerge sorting algorithm taken from:
  !  https://rosettacode.org/wiki/Sorting_algorithms/Merge_sort#Fortran
  !  and modified to match local type
  !! @A array of values to be sorted, returned sorted
  !! @N size of A
  !! @T I don't really understand T
  recursive subroutine sort_1d(A,N,T)

    integer, intent(in) :: N
    real(kind=DEFAULT_PRECISION), dimension(N), intent(in out) :: A
    real(kind=DEFAULT_PRECISION), dimension((N+1)/2), intent (out) :: T

    integer :: NA,NB
    real(kind=DEFAULT_PRECISION) :: V

    if (N < 2) return
    if (N == 2) then
      if (A(1) > A(2)) then
        V = A(1)
        A(1) = A(2)
        A(2) = V
      endif
      return
    endif
    NA=(N+1)/2
    NB=N-NA

    call sort_1d(A,NA,T)
    call sort_1d(A(NA+1),NB,T)

    if (A(NA) > A(NA+1)) then
      T(1:NA)=A(1:NA)
      call MergeSortMerge(T,NA,A(NA+1),NB,A,N)
    endif
    return

  end subroutine sort_1d


  !> Combines with sort_1d sorting algorithm taken from:
  !  https://rosettacode.org/wiki/Sorting_algorithms/Merge_sort#Fortran
  !  and modified to match local type and renamed to avoid confusion with intrinsic merge
  !  All parameters based on sort_1d.  No need to modify anything.
  subroutine MergeSortMerge(A,NA,B,NB,C,NC)

    integer, intent(in) :: NA,NB,NC                              ! Normal usage: NA+NB = NC
    real(kind=DEFAULT_PRECISION), intent(in out) :: A(NA)        ! B overlays C(NA+1:NC)
    real(kind=DEFAULT_PRECISION), intent(in)     :: B(NB)
    real(kind=DEFAULT_PRECISION), intent(in out) :: C(NC)

    integer :: I,J,K

    I = 1; J = 1; K = 1;
    do while(I <= NA .and. J <= NB)
      if (A(I) <= B(J)) then
         C(K) = A(I)
         I = I+1
      else
         C(K) = B(J)
         J = J+1
      endif
      K = K + 1
    enddo
    do while (I <= NA)
      C(K) = A(I)
      I = I + 1
      K = K + 1
    enddo
    return

  end subroutine MergeSortMerge
end module maths_mod
