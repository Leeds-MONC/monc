! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

! Description:
!   This module defines the controlling structure for LW calculations.

! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Radiation Control

MODULE lw_control_mod

USE def_control, ONLY: StrCtrl

IMPLICIT NONE

TYPE (StrCtrl), SAVE :: lw_control

END MODULE lw_control_mod
