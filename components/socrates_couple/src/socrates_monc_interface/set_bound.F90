! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Subroutine to set the grid used by the core radiation code
!
!------------------------------------------------------------------------------
SUBROUTINE set_bound(                                                       &

! Structures for the core radiation code interface
  control, atm, dimen, spectrum, bound,  socrates_derived_fields)

USE rad_pcf
USE def_control, ONLY: StrCtrl
USE def_atm,     ONLY: StrAtm 
USE def_spectrum, ONLY: StrSpecData
USE def_dimen,    ONLY: StrDim
USE def_bound,    ONLY: StrBound, allocate_bound
use def_socrates_derived_fields, only: str_socrates_derived_fields

IMPLICIT NONE

! Control options:
TYPE(StrCtrl),      INTENT(IN)    :: control

! Atmospheric properties:
TYPE(StrAtm),       INTENT(IN) :: atm

! Dimensions:
TYPE(StrDim),       INTENT(IN)  :: dimen

! Spectral data:
TYPE (StrSpecData), INTENT(IN)  :: spectrum

! Boundary conditions:
TYPE(StrBound),     INTENT(OUT) :: bound

! solar variables from MONC solar_position_angle_mod and timestep_callback
type (str_socrates_derived_fields), intent(in) :: socrates_derived_fields

call allocate_bound(bound, dimen, spectrum)

!     Set the surface basis functions for a Lambertian surface.
bound%n_brdf_basis_fnc = 1
!     By setting F_{1,0,0,0} equal to 4 we can set rho_alb equal to
!     the diffuse albedo
bound%f_brdf(1,0,0,0)  = 4.0 
!     Need to set this initially using input from config
bound%rho_alb(1:atm%n_profile, ip_surf_alb_diff, 1:spectrum%basic%n_band) = &
     socrates_derived_fields%albedoin1
bound%rho_alb(1:atm%n_profile, ip_surf_alb_dir,  1:spectrum%basic%n_band) = & 
     socrates_derived_fields%albedoin1
!     Holds the secant of the zenith angle for the fluxes
bound%zen_0(1:atm%n_profile) = socrates_derived_fields%sec_out
! sol_const = default_solar_const * scs (calced in timestep_callback)
! fractional_lit = lit in UM (is 0 or 1)
bound%solar_irrad(:) = socrates_derived_fields%sol_const * socrates_derived_fields%fraction_lit
!     pass from the MONC surface temperature
bound%t_ground(1:atm%n_profile) = socrates_derived_fields%srf_temperature ! Needs sorting

END SUBROUTINE set_bound
