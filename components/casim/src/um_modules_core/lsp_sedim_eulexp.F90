! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!  Microphysics hydrometeor Eulerian sedimentation scheme
MODULE lsp_sedim_eulexp_mod

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='LSP_SEDIM_EULEXP_MOD'

CONTAINS

SUBROUTINE lsp_sedim_eulexp(                                                   &
  points,m0,dhi,dhir,rho,rhor,                                                 &
  flux_fromabove, fallspeed_thislayer,                                         &
  mixratio_thislayer, fallspeed_fromabove,                                     &
  total_flux_out)

!USE lsprec_mod, ONLY: zero, one

! Use in KIND for large scale precip, used for compressed variables passed down
! from here
!USE um_types,             ONLY: real_lsprec

use variable_precision, only: wp, iwp, defp

USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim
IMPLICIT NONE

! Description:
!  Dummy replacement for lsp_sedim_eulexp, to permit compile of
!  CASIM sedimentation in MONC.

! Method:
!   Based on method described in Rotstayn (1997)(QJRMS, 123, 1227-1282)
!
! Code Owner: Please refer to the UM file CodeOwners.txt

! Subroutine arguments

      ! Intent (In)
INTEGER :: points       ! number of points to process

REAL (KIND=wp) ::                                                     &
  m0,                                                                          &
                         ! Small mass (kg/kg) defined in c_lspmic
  dhi(points),                                                                 &
                         ! CFL limit (s m-1)
  dhir(points),                                                                &
                         ! 1.0/DHI (m s-1)
  rho(points),                                                                 &
                         ! Air density (kg m-3)
  rhor(points),                                                                &
                         ! 1.0/Rho
  flux_fromabove(points),                                                      &
  fallspeed_thislayer(points)

    ! Intent (InOut)
REAL (KIND=wp) ::                                                     &
  mixratio_thislayer(points),                                                  &
  fallspeed_fromabove(points)

    ! Intent (Out)
REAL (KIND=wp) ::                                                     &
  total_flux_out(points)

! Local variables

REAL (KIND=wp) ::                                                     &
  mixratio_fromabove,                                                          &
                                 ! Mixing Ratio from above
  flux_out,                                                                    &
                                 ! Temporary flux out of layer
  expfactor                  ! Exponential Factor

INTEGER :: i                    ! Loop counter

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='LSP_SEDIM_EULEXP'


!-----------------------------------------------------------------------

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)


IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE lsp_sedim_eulexp
END MODULE lsp_sedim_eulexp_mod
