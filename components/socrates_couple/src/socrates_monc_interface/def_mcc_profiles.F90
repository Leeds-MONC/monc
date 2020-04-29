module def_mcc_profiles

  use datadefn_mod, only : DEFAULT_PRECISION
  
  implicit none

  type str_mcc_profiles

     ! all read from configuration
     integer :: levs
     ! number of McClatchey levels (read from file)
     integer :: irad_levs            ! the total number of levels for radiation
                                     ! monc levels + mcc to the top of atmosphere
     integer :: n_levels             ! irad_levs -1, not sure if this is needed in
                                     ! socrates. It is used in old LEM ES code.
     integer :: cut                  ! pressure level in the McClatchey profile that
                                     ! is closest to prefn(Z_SIZE) - pmindiff

     ! McClatchey profile arrays, which are read from file 
     ! Note: assume the all the McClatchey profiles exist on the grid interface
     real(kind=DEFAULT_PRECISION), allocatable :: p_level(:)
     real(kind=DEFAULT_PRECISION), allocatable :: t_level(:)
     real(kind=DEFAULT_PRECISION), allocatable :: q_level(:)
     real(kind=DEFAULT_PRECISION), allocatable :: o3_level(:)
     ! Derive Mcc profiles for the centre of grid. These profiles are only used to 
     ! work out the weightings for merge data
     real(kind=DEFAULT_PRECISION), allocatable :: p_n(:)
     real(kind=DEFAULT_PRECISION), allocatable :: t_n(:)
     real(kind=DEFAULT_PRECISION), allocatable :: q_n(:)
     real(kind=DEFAULT_PRECISION), allocatable :: o3_n(:)

  end type str_mcc_profiles
end module def_mcc_profiles
