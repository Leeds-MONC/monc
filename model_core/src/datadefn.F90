!> Contains common definitions for the data and datatypes used by MONC
module datadefn_mod
  use mpi
  implicit none

#ifndef TEST_MODE
  private
#endif

  integer, public, parameter :: STRING_LENGTH=150 !< Default length of strings
  integer, public, parameter :: LONG_STRING_LENGTH=STRING_LENGTH + 50!< Length of longer strings

  !< Single precision (32 bit) kind
  integer, public, parameter :: SINGLE_P = 6    !< precision (significant figures)
  integer, public, parameter :: SINGLE_R = 30   !< range (exponent)
  integer, public, parameter :: SINGLE_PRECISION = selected_real_kind(SINGLE_P, SINGLE_R)

  !< Double precision (64 bit) kind
  integer, public, parameter :: DOUBLE_P = 15   !< precision (significant figures)
  integer, public, parameter :: DOUBLE_R = 307  !< range (exponent)
  integer, public, parameter :: DOUBLE_PRECISION = selected_real_kind(DOUBLE_P, DOUBLE_R)

  !< Default precision which is used for prognostic data and calculations
  integer, public, parameter :: DEFAULT_PRECISION = DOUBLE_PRECISION 
  !< Solver precision is used in the interative solver
  integer, public, parameter :: SOLVER_PRECISION = DOUBLE_PRECISION 
  !< MPI communication type which we use for the prognostic and calculation data
  integer, public :: PRECISION_TYPE, SINGLE_PRECISION_TYPE, DOUBLE_PRECISION_TYPE


  ! Configuration precision toggle
  !   ---- HARD-CODED ONLY ----
  logical, public, parameter :: l_config_double = .true.  ! read config reals as DOUBLE_PRECISION
  !logical, public, parameter :: l_config_double = .false.  ! read config reals as fortran reals (old)
  integer, public :: config_precision, config_range

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

  !> Record the expected configuration precision for verification
    if (l_config_double) then
      config_precision = DOUBLE_P
      config_range = DOUBLE_R
    else
      config_precision = SINGLE_P
      config_range = SINGLE_R
    end if
  end subroutine init_data_defn
end module datadefn_mod
