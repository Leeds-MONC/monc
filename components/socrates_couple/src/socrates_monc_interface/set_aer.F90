! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Subroutine to set the grid used by the core radiation code
!
!------------------------------------------------------------------------------
SUBROUTINE set_aer(                                                       &

! Structures for the core radiation code interface
  control, atm, dimen, spectrum, aer )

USE rad_pcf
USE def_control, ONLY: StrCtrl
USE def_atm,     ONLY: StrAtm 
USE def_spectrum, ONLY: StrSpecData
USE def_dimen,    ONLY: StrDim
USE def_aer,    ONLY: StrAer, allocate_aer, allocate_aer_prsc

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
TYPE(StrAer),     INTENT(OUT) :: aer

call allocate_aer(aer, dimen, spectrum)
CALL allocate_aer_prsc(aer, dimen, spectrum)


END SUBROUTINE set_aer
