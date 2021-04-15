! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

! Dummy version of um_types so CASIM will build in both MONC and UM
! lsp_sed_eulexp (from the UM) requires the precision to be set using
! real_lsprec. In the UM and rose-stem, this can be either single
! or double. The following code ensures CASIM will build in MONC with
! lsp_sed_eulexp call and single precision microphysics
!

MODULE um_types

use variable_precision, only: wp

IMPLICIT NONE

!Large scale precipitation scheme
INTEGER, PARAMETER :: real_lsprec = wp


END MODULE um_types
