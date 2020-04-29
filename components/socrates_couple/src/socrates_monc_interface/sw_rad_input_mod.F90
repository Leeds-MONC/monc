! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
! Input control for sw radiation.

! Description:
!   Module containing control as used by the sw radiation code.

! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Radiation Control

! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP

MODULE sw_rad_input_mod

USE sw_control_mod, ONLY: sw_control
USE sw_control_default_mod, ONLY: sw_control_default
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

! Subroutine to set the input values of the sw control structure

SUBROUTINE sw_input(current_state)

USE rad_pcf, ONLY: ip_cloud_mix_max, ip_cloud_mix_random,               &
  ip_cloud_triple, ip_cloud_part_corr, ip_cloud_part_corr_cnv,          &
  ip_cloud_mcica, ip_cloud_clear,                                       &
  ip_solver_homogen_direct, ip_solver_mix_direct_hogan,                 &
  ip_solver_triple_hogan,                                               &
  ip_cloud_homogen, ip_cloud_ice_water, ip_cloud_conv_strat,            &
  ip_cloud_csiw,                                                        &
  ip_scaling, ip_mcica,                                                 &
  ip_max_rand, ip_exponential_rand, ip_rand
USE ereport_mod, ONLY: ereport

use state_mod, only : model_state_type
use optionsdatabase_mod, only : options_get_real,options_get_string,    & 
     options_get_integer, options_get_logical

IMPLICIT NONE

type(model_state_type), target, intent(inout) :: current_state
INTEGER :: j,                                                     &
           errorstatus      ! Return code : 0 Normal Exit : >0 Error
INTEGER :: namelist_unit

REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*) :: RoutineName
PARAMETER (   RoutineName='sw_rad_input_mod')
CHARACTER(LEN=errormessagelength) :: CMessage ! Error message if Errorstatus >0

IF (lhook) CALL dr_hook('SW_INPUT',zhook_in,zhook_handle)

  ! Set default values of control variables.

CALL sw_control_default(sw_control)

! Transfer namelist items to the data structure.
sw_control%spectral_file          =  &
     ADJUSTL(options_get_string(current_state%options_database, "spectral_file_sw"))
sw_control%i_gas_overlap          =  &
     options_get_integer(current_state%options_database, "i_gas_overlap_sw")
sw_control%i_st_water             =  & 
     options_get_integer(current_state%options_database, "i_water_sw")
sw_control%i_st_ice               =  &
     options_get_integer(current_state%options_database, "i_ice_sw")
sw_control%i_cloud_representation =  & 
     options_get_integer(current_state%options_database, "i_cloud_representation")

! Set i_cloud from combination of i_overlap and 

IF (sw_control%i_cloud_representation == ip_cloud_homogen) THEN

   sw_control%i_solver=ip_solver_homogen_direct
   IF (sw_control%i_inhom == ip_mcica) THEN
       ! The code below builds but if options is selected, MONC WILL NOT 
       ! run as other modules required - AH, 26/06/15
      sw_control%i_cloud=ip_cloud_mcica
   ELSE
      IF (sw_control%i_overlap == ip_max_rand) THEN
         sw_control%i_cloud=ip_cloud_mix_max
      ELSE IF (sw_control%i_overlap == ip_exponential_rand) THEN
         sw_control%i_cloud=ip_cloud_part_corr
      ELSE IF (sw_control%i_overlap == ip_rand) THEN
         sw_control%i_cloud=ip_cloud_mix_random
      ELSE
         ErrorStatus = 100
         CMessage = 'The selected cloud overlap is not available'
         CALL ereport(RoutineName, ErrorStatus, CMessage)
      END IF
   END IF
   
ELSE IF (sw_control%i_cloud_representation == ip_cloud_ice_water) THEN
   
   IF (sw_control%i_inhom == ip_mcica) THEN
       ! The code below builds but if options is selected, MONC WILL NOT 
       ! run as other modules required - AH, 26/06/15
      sw_control%i_cloud=ip_cloud_mcica
      sw_control%i_solver=ip_solver_homogen_direct
   ELSE
      sw_control%i_solver=ip_solver_mix_direct_hogan
      IF (sw_control%i_overlap == ip_max_rand) THEN
         sw_control%i_cloud=ip_cloud_mix_max
      ELSE IF (sw_control%i_overlap == ip_exponential_rand) THEN
         sw_control%i_cloud=ip_cloud_part_corr
      ELSE IF (sw_control%i_overlap == ip_rand) THEN
         sw_control%i_cloud=ip_cloud_mix_random
      ELSE
         ErrorStatus = 100
         CMessage = 'The selected cloud overlap is not available'
         CALL ereport(RoutineName, ErrorStatus, CMessage)
      END IF
   END IF

ELSE IF ((sw_control%i_cloud_representation == ip_cloud_conv_strat) .OR.  &
     (sw_control%i_cloud_representation == ip_cloud_csiw)) THEN
   
   IF (sw_control%i_inhom == ip_mcica) THEN
      ErrorStatus = 100
      CMessage = 'McICA is not compatible with the selected'//                 &
           ' cloud representation'
      CALL ereport(RoutineName, ErrorStatus, CMessage)
   ELSE
      sw_control%i_solver=ip_solver_triple_hogan
      IF (sw_control%i_overlap == ip_max_rand) THEN
         sw_control%i_cloud=ip_cloud_triple
      ELSE IF (sw_control%i_overlap == ip_exponential_rand) THEN
         sw_control%i_cloud=ip_cloud_part_corr_cnv
      ELSE IF (sw_control%i_overlap == ip_rand) THEN
         sw_control%i_cloud=ip_cloud_mix_random
      ELSE
         ErrorStatus = 100
         CMessage = 'The selected cloud overlap is not available'
         CALL ereport(RoutineName, ErrorStatus, CMessage)
      END IF
   END IF
   
ELSE
   
   ! No treatment of cloud for SW radiation
   sw_control%l_cloud        = .FALSE.
   sw_control%i_cloud        = ip_cloud_clear
   sw_control%i_solver       = ip_solver_homogen_direct
   sw_control%l_microphysics = .FALSE.
   
END IF

IF (lhook) CALL dr_hook('SW_INPUT',zhook_out,zhook_handle)

END SUBROUTINE sw_input

END MODULE sw_rad_input_mod
