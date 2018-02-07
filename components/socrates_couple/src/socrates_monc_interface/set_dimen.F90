! *****************************copyright*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************copyright*******************************
!
! Subroutine to set dimensions for the radiation code.
!
! Purpose:
!   To set the dimensions of arrays for use in the radiation code
!   depending on the options chosen.
!
!------------------------------------------------------------------------------
SUBROUTINE set_dimen(control, dimen, spectrum,                                 &
  n_points, model_levels, cloud_levels)

!!!NEED TO WORK OUT THE GRID FROM MERGE DATA + MODEL LEVELS in CURRENT_STATE

USE rad_pcf
USE def_control,  ONLY: StrCtrl
USE def_dimen,    ONLY: StrDim
USE def_spectrum, ONLY: StrSpecData
USE ereport_mod, ONLY: ereport
USE errormessagelength_mod, ONLY: errormessagelength


IMPLICIT NONE


! Control options:
  TYPE(StrCtrl),      INTENT(IN)    :: control

! Dimensions:
  TYPE(StrDim),       INTENT(INOUT) :: dimen

! Spectral data:
  TYPE (StrSpecData), INTENT(IN)    :: spectrum

INTEGER, INTENT(IN) :: n_points
!   Number of grid points to operate on
!! MONC: both variables below are irad_levs from the
!!       def_mcc_profiles.F90 (type)
INTEGER, INTENT(IN) :: model_levels
!   Number of theta levels in the model
INTEGER, INTENT(IN) :: cloud_levels
!   Number of potentially cloudy layers in the model

INTEGER                      :: ierr = i_normal
CHARACTER (LEN=*), PARAMETER :: RoutineName = 'set_dimen'
CHARACTER (LEN=errormessagelength) :: cmessage

! Grid dimensions
! (On some architectures an odd size is preferred to avoid memory
! bank conflicts).
dimen%nd_profile = 2*(n_points/2)+1
dimen%nd_layer = model_levels

! Cloud
dimen%nd_cloud_type      = 4
dimen%nd_cloud_component = 4
dimen%nd_cloud_representation = 4
SELECT CASE(control%i_cloud)
CASE (ip_cloud_column_max)
  dimen%nd_column = 3 * cloud_levels + 2
CASE DEFAULT
  dimen%nd_column = 1
END SELECT
dimen%nd_subcol_gen = 1
dimen%nd_subcol_req = 1

dimen%id_cloud_top = dimen%nd_layer + 1 - cloud_levels
IF (control%l_cloud) THEN
  dimen%nd_layer_clr = dimen%id_cloud_top - 1
ELSE
  dimen%nd_layer_clr = dimen%nd_layer
END IF

! Aerosol
dimen%nd_aerosol_mode = 1

! Arrays for prescribed optical properties (not used here)
dimen%nd_profile_aerosol_prsc   = 1
dimen%nd_profile_cloud_prsc     = 1
dimen%nd_opt_level_aerosol_prsc = 1
dimen%nd_opt_level_cloud_prsc   = 1


! Tiled surface.
dimen%nd_point_tile = MAX(1,n_points)
dimen%nd_tile     = 2


! Solver dependent dimensions
IF (control%i_angular_integration == ip_two_stream) THEN

  dimen%nd_viewing_level    = 1
  dimen%nd_radiance_profile = 1
  dimen%nd_j_profile        = 1
  dimen%nd_direction        = 1
  dimen%nd_brdf_basis_fnc   = 2
  dimen%nd_brdf_trunc       = 1
  dimen%nd_flux_profile     = dimen%nd_profile
  dimen%nd_channel          = 1

  dimen%nd_2sg_profile      = dimen%nd_profile
  dimen%nd_source_coeff     = 2
  dimen%nd_max_order        = 1
  dimen%nd_sph_coeff        = 1

  SELECT CASE(control%i_solver)

  CASE(ip_solver_mix_app_scat, ip_solver_mix_direct, ip_solver_mix_direct_hogan)
    dimen%nd_overlap_coeff=8
    dimen%nd_region=2

  CASE(ip_solver_triple_app_scat, ip_solver_triple, ip_solver_triple_hogan)
    dimen%nd_overlap_coeff=18
    dimen%nd_region=3

  CASE DEFAULT
    dimen%nd_overlap_coeff=1
    dimen%nd_region=2

  END SELECT

ELSE IF (control%i_angular_integration == ip_spherical_harmonic) THEN

  dimen%nd_direction        = 1
  dimen%nd_brdf_basis_fnc   = 2
  dimen%nd_radiance_profile = dimen%nd_profile
  dimen%nd_j_profile        = 1

  IF (control%i_sph_mode == ip_sph_mode_flux) THEN
    dimen%nd_viewing_level  = model_levels+1
    dimen%nd_brdf_trunc     = 1
    dimen%nd_flux_profile   = dimen%nd_profile
    
    ! When calculating fluxes all spectral bands are mapped into a
    ! single output channel.
    dimen%nd_channel        = 1
  ELSE
    dimen%nd_viewing_level  = 1
    dimen%nd_brdf_trunc     = 1
    dimen%nd_flux_profile   = 1
        
    ! We assume for now that each band in the spectral file
    ! corresponds to a particular channel.
    dimen%nd_channel        = spectrum%basic%n_band
  END IF

  dimen%nd_2sg_profile      = 1
  dimen%nd_source_coeff     = 1
  dimen%nd_overlap_coeff    = 1
  dimen%nd_region           = 1
  dimen%nd_max_order        = control%ls_global_trunc + 2
  IF (control%i_truncation == ip_trunc_triangular) THEN
    dimen%nd_sph_coeff =                                                       &
      (control%ls_global_trunc+3)*(control%ls_global_trunc+4)/2
  ELSE IF (control%i_truncation == ip_trunc_azim_sym) THEN
    dimen%nd_sph_coeff = control%ls_global_trunc+2
  ELSE
    cmessage = 'Illegal truncation'
    ierr=i_err_fatal
    GO TO 9999
  END IF

END IF


9999 CONTINUE
! Check error condition
IF (ierr /= i_normal) THEN
  CALL ereport(RoutineName, ierr, cmessage)
END IF

END SUBROUTINE set_dimen
