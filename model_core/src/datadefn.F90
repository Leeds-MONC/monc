!> Contains common definitions for the data and datatypes used by MONC
module datadefn_mod
  use mpi
  implicit none

#ifndef TEST_MODE
  private
#endif

  integer, public, parameter :: STRING_LENGTH=150 !< Default length of strings
  integer, public, parameter :: LONG_STRING_LENGTH=STRING_LENGTH + 50!< Length of longer strings

  integer, public, parameter :: SINGLE_PRECISION = selected_real_kind(6,30)   !< Single precision (32 bit) kind
  integer, public, parameter :: DOUBLE_PRECISION = selected_real_kind(15,307) !< Double precision (64 bit) kind

  !< Default precision which is used for prognostic data and calculations
  integer, public, parameter :: DEFAULT_PRECISION = DOUBLE_PRECISION 
  !< Solver precision is used in the interative solver
  integer, public, parameter :: SOLVER_PRECISION = DOUBLE_PRECISION 
  !< MPI communication type which we use for the prognostic and calculation data
  integer, public :: PRECISION_TYPE, SINGLE_PRECISION_TYPE, DOUBLE_PRECISION_TYPE

  public init_data_defn

contains

  !> Will initialise the data definitions. This should be called as soon as MONC starts up
  subroutine init_data_defn()
    if (DEFAULT_PRECISION .eq. DOUBLE_PRECISION) then
      PRECISION_TYPE = MPI_DOUBLE_PRECISION
    else
      PRECISION_TYPE = MPI_REAL
    endif
  !> Initialise single and double precision type for use of mixed precision, e.g
  !> in the iterative solver
    SINGLE_PRECISION_TYPE = MPI_REAL
    DOUBLE_PRECISION_TYPE = MPI_DOUBLE_PRECISION

  end subroutine init_data_defn
end module datadefn_mod
