! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Subroutine to set algorithmic options for the core radiation code
!
! Purpose:
!   Algorithmic options and array sizes to be set interactively
!   are determined.
!
!------------------------------------------------------------------------------
SUBROUTINE set_control(control, spectrum)


USE rad_pcf
USE def_control,  ONLY: StrCtrl, allocate_control
USE def_spectrum, ONLY: StrSpecData
USE ereport_mod,  ONLY: ereport
USE errormessagelength_mod, ONLY: errormessagelength

IMPLICIT NONE


! Control options:
TYPE(StrCtrl),      INTENT(INOUT) :: control

! Spectral data:
TYPE (StrSpecData), INTENT(IN)    :: spectrum

! Local variables.
INTEGER :: i
!   Loop variable

INTEGER                      :: ierr = i_normal
CHARACTER (LEN=*), PARAMETER :: RoutineName = 'set_control'
CHARACTER (LEN=errormessagelength) :: cmessage


! By default in MONC we assume mixing ratios so set following 
! to true
LOGICAL :: l_mixing_ratio                                       = .True.
!   True if mixing ratios are with respect to dry mass

! Set the last band to use as the last band in the spectral file
  control%last_band = spectrum%basic%n_band

! Allocate band-by-band control options
  CALL allocate_control(control, spectrum)

! Set diagnostic flags
SELECT CASE (control%isolir)

CASE (ip_solar)
  control%l_clear = .TRUE.
  control%l_cloud_extinction     = .TRUE.
  control%l_ls_cloud_extinction  = .TRUE.
  control%l_cnv_cloud_extinction = .FALSE.
  control%l_flux_direct_diag     = .FALSE.
  control%l_flux_down_diag       = .FALSE.
  control%l_flux_up_diag         = .FALSE.
  control%l_flux_down_diag_surf  = .FALSE.
  control%l_flux_down_clear_diag_surf = .FALSE.

CASE (ip_infra_red)
  control%l_clear = .TRUE.
  control%l_cloud_absorptivity     = .TRUE.
  control%l_ls_cloud_absorptivity  = .TRUE.
  control%l_cnv_cloud_absorptivity = .FALSE.

END SELECT

! Control flag for corrections to the direct solar flux at the surface
! for sloping terrain
control%l_orog = .FALSE.

! Decide on the final options for aerosols:
control%l_aerosol_mode = .FALSE.
control%l_aerosol = .FALSE.
control%l_aerosol_ccn=.FALSE.


! Set properties for individual bands.
DO i = 1, spectrum%basic%n_band
  control%map_channel(i)           = 1
  control%weight_band(i)           = 1.0
  control%i_scatter_method_band(i) = control%i_scatter_method
  control%i_gas_overlap_band(i)    = control%i_gas_overlap
  IF (ANY(spectrum%gas%i_scale_fnc(i,:) == ip_scale_ses2)) THEN
    ! If SES2 scaling is used in this band then the overlap must also use SES2:
    control%i_gas_overlap_band(i)  = ip_overlap_mix_ses2
  END IF
END DO


IF (control%i_angular_integration == ip_two_stream) THEN

  IF (control%l_rescale) control%n_order_forward=2
  control%l_tile=.FALSE.

ELSE IF (control%i_angular_integration == ip_spherical_harmonic) THEN

  IF (control%i_sph_mode == ip_sph_mode_flux) THEN
    ! Map all bands to a single output channel
    control%n_channel = 1
    control%map_channel(1:spectrum%basic%n_band) = 1
  ELSE
    ! Map each band to a separate output channel
    control%n_channel = spectrum%basic%n_band
    DO i = 1, spectrum%basic%n_band
      control%map_channel(i) = i
    END DO
  END IF

  IF (control%l_rescale) control%n_order_forward=control%ls_global_trunc+1

  ! As currently implemented, Euler's transformation is applied
  ! only in its most basic form, adding just half of the last
  ! term in an alternating series.
  IF (control%l_euler_trnf) THEN
    control%euler_factor=0.5
  ELSE
    control%euler_factor=1.0
  END IF

  ! Clear-sky fluxes are not available from the spherical harmonic
  ! code in the same call as cloudy fluxes yet. If required, they
  ! should be diagnosed by using a separate call to the code with
  ! clouds switched off.
  IF ( control%l_clear ) THEN
    cmessage = 'Clear-sky fluxes not directly available in harmonics'
    ierr=i_err_fatal
    GO TO 9999
  END IF

  ! We permit tiling of sea-ice points only with the two-stream
  ! option at present.
  control%l_tile=.FALSE.

END IF


9999 CONTINUE
! Check error condition
IF (ierr /= i_normal) THEN
  CALL ereport(RoutineName, ierr, cmessage)
END IF

END SUBROUTINE set_control
