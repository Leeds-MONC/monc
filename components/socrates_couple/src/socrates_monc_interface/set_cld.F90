! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Subroutine to set the grid used by the core radiation code
!
!------------------------------------------------------------------------------
SUBROUTINE set_cld(                                                       &

! Structures for the core radiation code interface
  control, atm, dimen, spectrum, cld,                                     &
! Grid
  n_profile, n_layer, nclds,                                              & 
! monc fields
  socrates_opt, merge_fields) 

USE rad_pcf
USE def_control, ONLY: StrCtrl
USE def_atm,     ONLY: StrAtm 
USE def_spectrum, ONLY: StrSpecData
USE def_dimen,    ONLY: StrDim
USE def_cld,    ONLY: StrCld, allocate_cld, allocate_cld_prsc
! monc-socrates couple structures
use def_merge_atm, only: str_merge_atm
use def_socrates_options, only: str_socrates_options

! Monc specific modules
use datadefn_mod, only : DEFAULT_PRECISION
use logging_mod, only : LOG_INFO, LOG_WARN, LOG_ERROR, &
       LOG_DEBUG, log_master_log, log_log, log_get_logging_level, &
       log_master_log
use science_constants_mod, only : pi

IMPLICIT NONE

! Control options:
TYPE(StrCtrl),      INTENT(IN)    :: control

! Atmospheric properties:
TYPE(StrAtm),       INTENT(IN) :: atm

! Dimensions:
TYPE(StrDim),       INTENT(IN)  :: dimen

! Spectral data:
TYPE (StrSpecData), INTENT(IN)  :: spectrum

! Cloud data:
TYPE(StrCld),     INTENT(OUT) :: cld

! McClatchey profiles plus monc profiles, flipped
type (str_merge_atm), intent(in) :: merge_fields

! MONC  options read from configuration
type (str_socrates_options), intent(in) :: socrates_opt

! declaration same as that in socrates_init on UM trunk
REAL ::                                                           &
    condensed_min_dim(dimen%nd_cloud_component)                   &
!       Minimum dimensions of condensed components
    , condensed_max_dim(dimen%nd_cloud_component)

!     local variables:
INTEGER ::                                                        &
     i, j, l                                                      &
!             loop variable
     , i_scheme
!             parametrization scheme

INTEGER, intent(in) ::                                            &
     n_profile                                                    &
     !             NUMBER OF PROFILES
     , n_layer                                                      &
     !             Number of layers seen by the radiation code
     , nclds
     !             NUMBER OF CLOUDY LEVELS

!     Parameters for the aggregate parametrization.
REAL, PARAMETER :: a0_agg_cold = 7.5094588e-04
REAL, PARAMETER :: b0_agg_cold = 5.0830326e-07
REAL, PARAMETER :: a0_agg_warm = 1.3505403e-04
REAL, PARAMETER :: b0_agg_warm = 2.6517429e-05
REAL, PARAMETER :: t_switch    = 216.208
REAL, PARAMETER :: t0_agg      = 279.5
REAL, PARAMETER :: s0_agg      = 0.05

!     Coefficients for effective size of ice crystal from Sun scheme
REAL, PARAMETER :: a1_s=45.9866
REAL, PARAMETER :: a2_s=0.2214
REAL, PARAMETER :: a3_s=0.7957
REAL, PARAMETER :: a4_s=0.2535
REAL, PARAMETER :: a5_s=83.15
REAL, PARAMETER :: a6_s=1.2351
REAL, PARAMETER :: a7_s=0.0105
!     Latitude threshold for division of ice crystal correction
REAL, PARAMETER :: a8_s=30.0

!     Parameters for Baran's diagnostic parameterisation of effective size.
!     m and n are the slope and intercept (positive sign), and max and min the 
!     maximum and minimum values. Values for parameterisation expressed in
!     units of metres.
REAL, PARAMETER :: m_ice_baran = 1.868e-6
REAL, PARAMETER :: n_ice_baran = 353.613e-6
REAL, PARAMETER :: min_ice_baran = 7.0e-6
REAL, PARAMETER :: max_ice_baran = 156.631e-6

! Variables for liquid re calculation
REAL :: beta

!     Functions called:
INTEGER, EXTERNAL :: set_n_cloud_parameter
!       Function to find number of parameters for clouds

call allocate_cld(cld, dimen, spectrum)
call allocate_cld_prsc(cld, dimen, spectrum)

!     Select parametrization for water in stratiform clouds:
IF ( (control%i_st_water <= spectrum%dim%nd_drop_type) .AND.      &
     (spectrum%drop%l_drop_type(control%i_st_water)) ) THEN
   i_scheme=spectrum%drop%i_drop_parm(control%i_st_water)
   cld%i_condensed_param(ip_clcmp_st_water)=i_scheme
   cld%condensed_n_phf(ip_clcmp_st_water)=                        &
        spectrum%drop%n_phf(control%i_st_water)
   condensed_min_dim(ip_clcmp_st_water)                           &
        =spectrum%drop%parm_min_dim(control%i_st_water)
   condensed_max_dim(ip_clcmp_st_water)                           &
        =spectrum%drop%parm_max_dim(control%i_st_water)
ELSE
   call log_master_log &
        (LOG_ERROR, "Socrates error: no data exist for type of selected droplet - STOP")
END IF

DO i=1, spectrum%basic%n_band
   ! DEPENDS ON: set_n_cloud_parameter
   DO j=1, set_n_cloud_parameter(i_scheme                         &
        , ip_clcmp_st_water, cld%condensed_n_phf(ip_clcmp_st_water))
      cld%condensed_param_list(j, ip_clcmp_st_water, i)               &
           =spectrum%drop%parm_list(j, i, control%i_st_water)
   END DO
END DO

!     Select parametrization for ice in stratiform clouds:
IF ( (control%i_st_ice <= spectrum%dim%nd_ice_type) .AND.                               &
     (spectrum%ice%l_ice_type(control%i_st_ice)) ) THEN
   i_scheme=spectrum%ice%i_ice_parm(control%i_st_ice)
   cld%i_condensed_param(ip_clcmp_st_ice)=i_scheme
   cld%condensed_n_phf(ip_clcmp_st_ice)=                                   &
        spectrum%ice%n_phf(control%i_st_ice)
   condensed_min_dim(ip_clcmp_st_ice)                             &
        =spectrum%ice%parm_min_dim(control%i_st_ice)
   condensed_max_dim(ip_clcmp_st_ice)                             &
        =spectrum%ice%parm_max_dim(control%i_st_ice)
ELSE
  call log_master_log &
        (LOG_ERROR, "Socrates error: no data exist for type of ice crystal - STOP")
END IF

DO i=1, spectrum%basic%n_band
   ! DEPENDS ON: set_n_cloud_parameter
   DO j=1, set_n_cloud_parameter(i_scheme                         &
        , ip_clcmp_st_ice, cld%condensed_n_phf(ip_clcmp_st_ice))
      cld%condensed_param_list(j, ip_clcmp_st_ice, i)                 &
           = spectrum%ice%parm_list(j, i, control%i_st_ice)
   END DO
END DO

IF (socrates_opt%cloud_representation == ip_cloud_ice_water) THEN

  cld%n_cloud_type=2

  ! Here the clouds are split into two separate types.
  ! The partitioning between ice and water is regarded as
  ! determining the areas within the grid_box covered by
  ! ice or water cloud, rather than as determining the in-cloud
  ! mixing ratios. The grid-box mean ice water contents may
  ! be predicted by the ice microphysics scheme or may be
  ! determined as a function of the temperature (LSP_FOCWWIL).

  ! Set the components within the clouds. Here we have two
  ! components: stratiform ice and water.
  cld%n_condensed=2
  cld%type_condensed(1)=ip_clcmp_st_water
  cld%type_condensed(2)=ip_clcmp_st_ice
  
  !DO j=1, n_layer
  !  DO l=1, n_profile
  !    cld%condensed_mix_ratio(l, j, ip_clcmp_cnv_water)=0.0e+00
  !    cld%condensed_mix_ratio(l, j, ip_clcmp_cnv_ice)=0.0e+00
  !  END DO
  !END DO
 
  ! Stratiform clouds:
  ! Note: n_layer is the merge_field index and using merge_field so cloud arrays are
  !       the correct way for the radiation, i.e. k=0 is TOA k=n_layer is surface
  DO j = 1, n_layer
     DO l=1, n_profile
        ! calculate the total_water (ice + water)
        
       ! Add the merged fields here 
        IF (merge_fields%total_cloud_fraction(j) > &
             EPSILON(merge_fields%total_cloud_fraction(j))) THEN
          cld%condensed_mix_ratio(l, j, ip_clcmp_st_water) =                    &
               ( merge_fields%ql_n(j) +  merge_fields%qi_n(j) )/ &
               merge_fields%total_cloud_fraction(j)
       ELSE
          cld%condensed_mix_ratio(l, j, ip_clcmp_st_water) = 0.0
       END IF
       cld%condensed_mix_ratio(l, j, ip_clcmp_st_ice) =                         &
            cld%condensed_mix_ratio(l, j, ip_clcmp_st_water)
    END DO
 END DO
 !
 ! Cloud fractions:
 !
 DO j = 1, n_layer
    DO l=1, n_profile
       cld%w_cloud(l, j) = merge_fields%total_cloud_fraction(j)
      ! Partition clouds using the ratio of cloud water contents done
      ! in merge fields, so set using the merged fields here.
      cld%frac_cloud(l, j, ip_cloud_type_sw) = merge_fields%liquid_cloud_fraction(j)
      cld%frac_cloud(l, j, ip_cloud_type_si) = merge_fields%ice_cloud_fraction(j)
    END DO      ! L (N_PROFILE)
 END DO      ! I (N_LAYER)

  if (socrates_opt%l_fix_re) then
     DO j = 1, n_layer
        DO l=1, n_profile 
           cld%condensed_dim_char(l, j, ip_clcmp_cnv_water) = 0.0_DEFAULT_PRECISION
           ! AH - convert fixed re from microns to metres
           cld%condensed_dim_char(l, j, ip_clcmp_st_water)  = socrates_opt%fixed_cloud_re * 1.e-6
        enddo
     enddo
  endif

  if (socrates_opt%l_use_ndrop) then
     IF (socrates_opt%l_use_liu_spec) THEN
!!$        DO j = 1, n_layer
!!$           ! Find the total mixing ratio of water substance in the cloud.
!!$           DO l=1, n_profile
!!$         
!!$              ! Apply Liu spectral dispersion
!!$              beta = aparam                                               & 
!!$             *((MAX(eps,cld%condensed_mix_ratio(l, j, ip_clcmp_st_water)*atm%density(l,j)*1.0e-3 &
!!$             /(merge_fields%cloudnumber_n(j)*1e-06))**(bparam))
!!$
!!$              cld%condensed_dim_char(l, j, ip_clcmp_cnv_water) = 0.0_DEFAULT_PRECISION
!!$              cld%condensed_dim_char(l, j, ip_clcmp_st_water)  = MAX(0.0_DEFAULT_PRECISION,       &
!!$                   3.0_DEFAULT_PRECISION*cld%condensed_mix_ratio(l, j, ip_clcmp_st_water)*atm%density(l, j) &
!!$                   /(4.0_DEFAULT_PRECISION*pi*socrates_opt%rho_water* &
!!$                   (beta**(-3.0))*merge_fields%cloudnumber_n(j))) &
!!$                   **(1.0_DEFAULT_PRECISION/3.0_DEFAULT_PRECISION)
!!$           END DO
!!$        END DO
     ELSE
        ! We only need to find land across n_profile
        DO j=1, n_layer
           ! Find the total mixing ratio of water substance in the cloud.
           DO l=1, n_profile
              cld%condensed_dim_char(l, j, ip_clcmp_cnv_water) = 0.0_DEFAULT_PRECISION
              cld%condensed_dim_char(l, j, ip_clcmp_st_water)  = MAX(0.0_DEFAULT_PRECISION,       &
                   3.0_DEFAULT_PRECISION*cld%condensed_mix_ratio(l, j, ip_clcmp_st_water)         &
                   *atm%density(l, j)/(4.0_DEFAULT_PRECISION*pi*socrates_opt%rho_water*      &
                   socrates_opt%kparam*merge_fields%cloudnumber_n(j))      &
                   **(1.0_DEFAULT_PRECISION/3.0_DEFAULT_PRECISION))
            END DO
        END DO
     END IF
  END if
endif
 
  ! Constrain the sizes of droplets to lie within the range of
  ! validity of the parametrization scheme.
DO j=1, n_layer
   DO l=1, n_profile
      cld%condensed_dim_char(l, j, ip_clcmp_st_water)              &
           =MAX(condensed_min_dim(ip_clcmp_st_water)             &
           , MIN(condensed_max_dim(ip_clcmp_st_water)            &
           , cld%condensed_dim_char(l, j, ip_clcmp_st_water)))
      cld%condensed_dim_char(l, j, ip_clcmp_cnv_water)             &
           =MAX(condensed_min_dim(ip_clcmp_cnv_water)            &
           , MIN(condensed_max_dim(ip_clcmp_cnv_water)           &
           , cld%condensed_dim_char(l, j, ip_clcmp_cnv_water)))

      cld%condensed_dim_char(l, j, ip_clcmp_st_ice) = 0.0
   END DO
   
END DO

!     SET THE CHARACTERISTIC DIMENSIONS OF ICE CRYSTALS:

!     ICE CRYSTALS IN STRATIFORM CLOUDS:

SELECT CASE (cld%i_condensed_param(ip_clcmp_st_ice))

CASE (ip_ice_agg_de, ip_ice_agg_de_sun)

  !      Aggregate parametrization based on effective dimension.
  !      In the initial form, the same approach is used for stratiform
  !      and convective cloud.

  !      The fit provided here is based on Stephan Havemann's fit of
  !      Dge with temperature, consistent with David Mitchell's treatment
  !      of the variation of the size distribution with temperature. The
  !      parametrization of the optical properties is based on De
  !      (=(3/2)volume/projected area), whereas Stephan's fit gives Dge
  !      (=(2*SQRT(3)/3)*volume/projected area), which explains the
  !      conversion factor. The fit to Dge is in two sections, because
  !      Mitchell's relationship predicts a cusp at 216.208 K. Limits
  !      of 8 and 124 microns are imposed on Dge: these are based on this
  !      relationship and should be reviewed if it is changed. Note also
  !      that the relationship given here is for polycrystals only.
  DO i=1, n_layer
    DO l=1, n_profile
      !          Preliminary calculation of Dge.
      IF (atm%t(l, i) < t_switch) THEN
        cld%condensed_dim_char(l, i, ip_clcmp_st_ice)                  &
          = a0_agg_cold*EXP(s0_agg*(atm%t(l, i)-t0_agg))+b0_agg_cold
      ELSE
        cld%condensed_dim_char(l, i, ip_clcmp_st_ice)                  &
          = a0_agg_warm*EXP(s0_agg*(atm%t(l, i)-t0_agg))+b0_agg_warm
      END IF
      !          Limit and convert to De.
      cld%condensed_dim_char(l, i, ip_clcmp_st_ice)                    &
        = (3.0/2.0)*(3.0/(2.0*SQRT(3.0)))*                         &
          MIN(1.24e-04, MAX(8.0e-06,                               &
          cld%condensed_dim_char(l, i, ip_clcmp_st_ice)))
    END DO
  END DO

END SELECT ! I_CONDENSED_PARAM(IP_CLCMP_ST_ICE)
!
!     CONSTRAIN THE SIZES OF ICE CRYSTALS TO LIE WITHIN THE RANGE
!     OF VALIDITY OF THE PARAMETRIZATION SCHEME.
DO i=1, n_layer
  DO l=1, n_profile
    cld%condensed_dim_char(l, i, ip_clcmp_st_ice)                   &
       =MAX(condensed_min_dim(ip_clcmp_st_ice)                  &
       , MIN(condensed_max_dim(ip_clcmp_st_ice)                 &
       , cld%condensed_dim_char(l, i, ip_clcmp_st_ice)))
    cld%condensed_dim_char(l, i, ip_clcmp_cnv_ice)                  &
       =MAX(condensed_min_dim(ip_clcmp_cnv_ice)                 &
       , MIN(condensed_max_dim(ip_clcmp_cnv_ice)                &
       , cld%condensed_dim_char(l, i, ip_clcmp_cnv_ice)))
  END DO
END DO


END SUBROUTINE set_cld
