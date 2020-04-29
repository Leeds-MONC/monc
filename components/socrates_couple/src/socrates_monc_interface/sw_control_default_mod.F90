! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

! Subroutine to set the default values of the control structure.

! Description:
!   Performs suitable default initializations of the controlling 
!   structure for SW radiation.

! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Radiation Control

MODULE sw_control_default_mod

CONTAINS

SUBROUTINE sw_control_default(sw_control)

USE def_control, ONLY: StrCtrl
USE yomhook,     ONLY: lhook, dr_hook
USE parkind1,    ONLY: jprb, jpim
USE rad_pcf,     ONLY: ip_solar, ip_two_stream, ip_pifm80, &
                       ip_solver_homogen_direct, ip_trunc_triangular, &
                       ip_sph_direct, ip_sph_mode_rad, ip_scatter_full, &
                       ip_max_rand, &
                       ip_homogeneous, ip_overlap_k_eqv_scl, &
                       ip_cloud_ice_water

IMPLICIT NONE

TYPE (StrCtrl), INTENT(INOUT) :: sw_control
!   The block of controlling options for the code

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle


IF (lhook) CALL dr_hook('SW_CONTROL_DEFAULT',zhook_in,zhook_handle)

! Spectral region and bands
sw_control%isolir     = ip_solar
sw_control%first_band = 1

! Physical processes
sw_control%l_microphysics = .TRUE.
sw_control%l_gas          = .TRUE.
sw_control%l_rayleigh     = .TRUE.
sw_control%l_continuum    = .TRUE.
sw_control%l_cloud        = .TRUE.
sw_control%l_drop         = .TRUE.
sw_control%l_ice          = .TRUE.
sw_control%l_aerosol      = .FALSE.
sw_control%l_aerosol_ccn  = .FALSE.

! Properties of clouds
sw_control%l_local_cnv_partition = .FALSE.
sw_control%l_global_cloud_top    = .TRUE.

! Angular integration (including algorithmic options)
sw_control%n_channel              = 1
sw_control%i_angular_integration  = ip_two_stream
sw_control%i_2stream              = ip_pifm80
sw_control%i_solver_clear         = ip_solver_homogen_direct
sw_control%n_order_gauss          = 0
sw_control%i_truncation           = ip_trunc_triangular
sw_control%i_sph_algorithm        = ip_sph_direct
sw_control%n_order_phase_solar    = 1
sw_control%ls_global_trunc        = 9
sw_control%ms_min                 = 0
sw_control%ms_max                 = 0
sw_control%ls_brdf_trunc          = 0
sw_control%accuracy_adaptive      = 1.0e-04
sw_control%l_rescale              = .TRUE.
sw_control%l_henyey_greenstein_pf = .TRUE.
sw_control%i_sph_mode             = ip_sph_mode_rad
sw_control%i_solar_src            = 3
sw_control%i_scatter_method       = ip_scatter_full

! Switches for diagnostic output
sw_control%l_blue_flux_surf = .TRUE.

! Satellite data
sw_control%sat_hgt      = 0.0
sw_control%sat_lon      = 0.0
sw_control%sat_lat      = 0.0
sw_control%max_view_lon = 0.0
sw_control%min_view_lon = 0.0
sw_control%max_view_lat = 0.0
sw_control%min_view_lat = 0.0

! Originally in sw_rad_iput, defaults do not need to change for MONC
sw_control%i_gas_overlap          = ip_overlap_k_eqv_scl
sw_control%l_o2                   = .TRUE.
sw_control%l_o3                   = .TRUE.
sw_control%l_h2o                  = .TRUE.
sw_control%l_co2                  = .TRUE.
sw_control%l_n2o                  = .TRUE.
sw_control%l_ch4                  = .TRUE.
sw_control%l_co                   = .FALSE.
sw_control%l_nh3                  = .FALSE.
sw_control%l_tio                  = .FALSE.
sw_control%l_vo                   = .FALSE.
sw_control%l_h2                   = .FALSE.
sw_control%l_he                   = .FALSE.
sw_control%l_na                   = .FALSE.
sw_control%l_k                    = .FALSE.
sw_control%i_st_water             = 5
sw_control%i_cnv_water            = 5
sw_control%i_st_ice               = 8
sw_control%i_cnv_ice              = 8
sw_control%i_cloud_representation = ip_cloud_ice_water
sw_control%i_overlap              = ip_max_rand
sw_control%i_inhom                = ip_homogeneous

IF (lhook) CALL dr_hook('SW_CONTROL_DEFAULT',zhook_out,zhook_handle)

END SUBROUTINE sw_control_default

END MODULE sw_control_default_mod
