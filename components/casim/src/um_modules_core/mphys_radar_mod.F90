! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
MODULE mphys_radar_mod

! Description:
! Holds reflectivity constants required by the microphysics
! scheme

! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Large Scale Precipitation

! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 programming standards.

!  USE mphys_constants_mod, ONLY: mprog_min
  USE science_constants_mod, ONLY: mprog_min

IMPLICIT NONE

REAL, PARAMETER :: kliq  = 0.93
REAL, PARAMETER :: kice  = 0.174
REAL, PARAMETER :: mm6m3 = 1.0e18

! Define reflectivity limit
REAL, PARAMETER :: ref_lim = -35.0 ! dBZ
! Convert this to linear units (mm6 m-3) using 10.0**p
! Where p = -35 dBZ / 10.0 .
REAL, PARAMETER :: ref_lim_lin = 3.1623e-4

! Mixing ratio limit (below which we ignore the species)
! Set this to be the same as the absolute value used in the
! rest of the microphysics.
REAL, PARAMETER :: mr_lim = 1.e-8
!REAL :: mr_lim = mprog_min

! Cloud fraction limit (below which we ignore to avoid massive
! reflectivity values)- currently set as 1% of the grid box
REAL, PARAMETER :: cf_lim = 0.01

! Cloud drop number concentration limit (below which we ignore liquid cloud)
! Equivalent to 5 per cm3, used by most of the UM routines
REAL, PARAMETER :: nd_lim  = 5.0e6

REAL, PARAMETER :: rho_g   = 500.0 ! 900.0 in Mark's HWT code.
REAL, PARAMETER :: rho_i   = 900.0
REAL, PARAMETER :: rho_i2  = 900.0

REAL, PARAMETER :: ref_mom = 4.0 ! Radar reflectivity moment


END MODULE mphys_radar_mod
