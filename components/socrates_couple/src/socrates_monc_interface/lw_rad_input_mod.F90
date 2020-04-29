! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
! Input control for lw radiation.

! Description:
!   Module containing control as used by the lw radiation code.

! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Radiation Control

! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP

MODULE lw_rad_input_mod

!USE lw_control_mod, ONLY: lw_control
USE lw_control_default_mod, ONLY: lw_control_default
USE missing_data_mod, ONLY: imdi
USE yomhook,  ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim
USE errormessagelength_mod, ONLY: errormessagelength

IMPLICIT NONE

!DrHook-related parameters
INTEGER(KIND=jpim), PARAMETER, PRIVATE :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER, PRIVATE :: zhook_out = 1

! ----------------------------------------------------------------
CONTAINS

! Subroutine to set the input values of the lw control structure.

SUBROUTINE lw_input(current_state, lw_control)

USE rad_pcf, ONLY: ip_cloud_mix_max, ip_cloud_mix_random,               &
  ip_cloud_triple, ip_cloud_part_corr, ip_cloud_part_corr_cnv,          &
  ip_cloud_mcica, ip_cloud_clear,                                       &
  ip_solver_homogen_direct, ip_solver_mix_direct_hogan,                 &
  ip_solver_triple_hogan, ip_solver_triple_app_scat,                    &
  ip_solver_mix_app_scat, ip_solver_no_scat,                            &
  ip_cloud_homogen, ip_cloud_ice_water, ip_cloud_conv_strat,            &
  ip_cloud_csiw,                                                        &
  ip_scaling, ip_mcica,                                                 &
  ip_max_rand, ip_exponential_rand, ip_rand,                            &
  ip_scatter_approx, ip_no_scatter_abs, ip_no_scatter_ext
USE ereport_mod, ONLY: ereport
USE def_control,  ONLY: StrCtrl

use state_mod, only : model_state_type
use optionsdatabase_mod, only : options_get_real, options_get_string,    & 
     options_get_integer, options_get_logical

IMPLICIT NONE

type(model_state_type), target, intent(inout) :: current_state

TYPE(StrCtrl),      INTENT(INOUT) :: lw_control

INTEGER :: j,                        &
     errorstatus      ! Return code : 0 Normal Exit : >0 Error
INTEGER :: namelist_unit

REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*) :: RoutineName
PARAMETER (   RoutineName='lw_rad_input_mod')
CHARACTER(LEN=errormessagelength) :: CMessage  ! Error message if Errorstatus >0

IF (lhook) CALL dr_hook('LW_INPUT',zhook_in,zhook_handle)

! Set default values of control variables.

CALL lw_control_default(lw_control)

  lw_control%spectral_file          =  &
       ADJUSTL(options_get_string(current_state%options_database, "spectral_file_lw"))
  lw_control%i_gas_overlap          =  &
       options_get_integer(current_state%options_database, "i_gas_overlap_lw")
  lw_control%i_cloud_representation =  &
       options_get_integer(current_state%options_database, "i_cloud_representation")
  lw_control%i_st_water             =  &
       options_get_integer(current_state%options_database, "i_water_lw")
  lw_control%i_st_ice               =  &
       options_get_integer(current_state%options_database, "i_ice_lw")
  lw_control%i_scatter_method       =  & 
       options_get_integer(current_state%options_database, "i_scatter_method_lw")

! Set i_cloud from combination of i_overlap and i_representation

  IF (lw_control%i_cloud_representation == ip_cloud_homogen) THEN

    lw_control%i_solver=ip_solver_homogen_direct
    IF (lw_control%i_inhom == ip_mcica) THEN
       ! The code below builds but if options is selected, MONC WILL NOT 
       ! run as other modules required - AH, 26/06/15
      lw_control%i_cloud=ip_cloud_mcica
    ELSE
      IF (lw_control%i_overlap == ip_max_rand) THEN
        lw_control%i_cloud=ip_cloud_mix_max
      ELSE IF (lw_control%i_overlap == ip_exponential_rand) THEN
        lw_control%i_cloud=ip_cloud_part_corr
      ELSE IF (lw_control%i_overlap == ip_rand) THEN
        lw_control%i_cloud=ip_cloud_mix_random
      ELSE
        ErrorStatus = 100
        CMessage = 'The selected cloud overlap is not available'
        CALL ereport(RoutineName, ErrorStatus, CMessage)
      END IF
    END IF

  ELSE IF (lw_control%i_cloud_representation == ip_cloud_ice_water) THEN

    IF (lw_control%i_inhom == ip_mcica) THEN
      lw_control%i_cloud=ip_cloud_mcica
      IF ( (lw_control%i_scatter_method == ip_no_scatter_abs) .OR.   &
           (lw_control%i_scatter_method == ip_no_scatter_ext) ) THEN
         ! The code below builds but if options is selected, MONC WILL NOT 
         ! run as other modules required - AH, 26/06/15 
        lw_control%i_solver      =ip_solver_no_scat
        lw_control%i_solver_clear=ip_solver_no_scat
      ELSE
        lw_control%i_solver=ip_solver_homogen_direct
      END IF
    ELSE
      IF (lw_control%i_scatter_method == ip_scatter_approx) THEN
        lw_control%i_solver=ip_solver_mix_app_scat
      ELSE
        lw_control%i_solver=ip_solver_mix_direct_hogan
      END IF
      IF (lw_control%i_overlap == ip_max_rand) THEN
        lw_control%i_cloud=ip_cloud_mix_max
      ELSE IF (lw_control%i_overlap == ip_exponential_rand) THEN
        lw_control%i_cloud=ip_cloud_part_corr
      ELSE IF (lw_control%i_overlap == ip_rand) THEN
        lw_control%i_cloud=ip_cloud_mix_random
      ELSE
        ErrorStatus = 100
        CMessage = 'The selected cloud overlap is not available'
        CALL ereport(RoutineName, ErrorStatus, CMessage)
      END IF
    END IF

  ELSE IF ((lw_control%i_cloud_representation == ip_cloud_conv_strat) .OR.  &
           (lw_control%i_cloud_representation == ip_cloud_csiw)) THEN

    IF (lw_control%i_inhom == ip_mcica) THEN
      ErrorStatus = 100
      CMessage = 'McICA is not compatible with the selected'//                 &
                 ' cloud representation'
      CALL ereport(RoutineName, ErrorStatus, CMessage)
    ELSE
      IF (lw_control%i_scatter_method == ip_scatter_approx) THEN
        lw_control%i_solver=ip_solver_triple_app_scat
      ELSE
        lw_control%i_solver=ip_solver_triple_hogan
      END IF
      IF (lw_control%i_overlap == ip_max_rand) THEN
        lw_control%i_cloud=ip_cloud_triple
      ELSE IF (lw_control%i_overlap == ip_exponential_rand) THEN
        lw_control%i_cloud=ip_cloud_part_corr_cnv
      ELSE IF (lw_control%i_overlap == ip_rand) THEN
        lw_control%i_cloud=ip_cloud_mix_random
      ELSE
        ErrorStatus = 100
        CMessage = 'The selected cloud overlap is not available'
        CALL ereport(RoutineName, ErrorStatus, CMessage)
      END IF
    END IF

  ELSE

    ! No treatment of cloud for LW radiation
    lw_control%l_cloud=.FALSE.
    lw_control%i_cloud=ip_cloud_clear
    IF ( (lw_control%i_scatter_method == ip_no_scatter_abs) .OR.            &
         (lw_control%i_scatter_method == ip_no_scatter_ext) ) THEN
      lw_control%i_solver      =ip_solver_no_scat
      lw_control%i_solver_clear=ip_solver_no_scat
    ELSE
      lw_control%i_solver=ip_solver_homogen_direct
    END IF
    lw_control%l_microphysics=.FALSE.
  END IF

IF (lhook) CALL dr_hook('LW_INPUT',zhook_out,zhook_handle)

END SUBROUTINE lw_input

END MODULE lw_rad_input_mod
