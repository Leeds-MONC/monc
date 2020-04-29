! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Subroutine to set the grid used by the core radiation code
!
!------------------------------------------------------------------------------
SUBROUTINE set_atm(                                                       &

! Structures for the core radiation code interface
  control, atm, dimen, spectrum,                                           &

! Grid
  n_profile, n_layer,                                                     &
! monc fields
  socrates_opt, merge_fields)

USE rad_pcf
USE def_control, ONLY: StrCtrl
USE def_atm,     ONLY: StrAtm, allocate_atm
USE def_spectrum, ONLY: StrSpecData
USE def_dimen,    ONLY: StrDim
use def_socrates_derived_fields, only: str_socrates_derived_fields

USE gas_list_pcf, ONLY: ip_h2o, ip_co2, ip_o3, ip_n2o,     &
                        ip_ch4, ip_o2, ip_cfc11, ip_cfc12, &
                        ip_cfc113, ip_cfc114, ip_hcfc22,   &
                        ip_hfc125, ip_hfc134a
! monc-socrates couple structures
use def_merge_atm, only: str_merge_atm
use def_socrates_options, only: str_socrates_options
! From MONC model_core use science constants for consitency
use datadefn_mod, only : DEFAULT_PRECISION
use science_constants_mod, only: r, ratio_mol_wts

IMPLICIT NONE

! Control options:
TYPE(StrCtrl),      INTENT(IN)    :: control

! Atmospheric properties:
TYPE(StrAtm),       INTENT(OUT) :: atm

! Dimensions:
TYPE(StrDim),       INTENT(IN)  :: dimen

! Spectral data:
TYPE (StrSpecData), INTENT(IN)  :: spectrum

! McClatchey profiles plus monc profiles, flipped
type (str_merge_atm), intent(in) :: merge_fields
!
! MONC  options read from configuration
type (str_socrates_options), intent(in) :: socrates_opt

INTEGER, INTENT(IN) :: n_profile
!   Number of atmospheric profiles for radiation calculations
INTEGER, INTENT(IN) :: n_layer
!   Number of atmospheric layers for radiation calculations

integer :: i, k, l ! loop counters

call allocate_atm(atm, dimen, spectrum)

! Setup atmosphere for radiation (upside down!)
atm%n_profile = n_profile
atm%n_layer   = n_layer

do l=1, atm%n_profile

   atm%p_level(l,:) = merge_fields%pres_level(:) ! pref + mcclatchy going from the top 1:n_layer+1 
   atm%t_level(l,:) = merge_fields%t_level(:)! equivalent to t_bdy and tac = atm%t
   atm%mass(l,:)   = merge_fields%mass(:)
   atm%p(l,:)   = merge_fields%pres_n(:)
   atm%t(l,:)   =  merge_fields%t_n(:)  
   if (Spectrum%Cont%index_water > 0) THEN
      ! The following lines come from socrates:src/modules_core/rad_ccf. 
      ! and accounts for the impact of water vapour of air density
      DO k=1, atm%n_layer
         atm%density(l, k)=atm%p(l, k)/(r*atm%t(l, k)*(1.0e+00_DEFAULT_PRECISION     &
              + (ratio_mol_wts-1.0_DEFAULT_PRECISION)                               &
              *merge_fields%qv_n(k)))
      enddo
   else
      atm%density(l, :)=atm%p(l, :)/(r*atm%t(l, :))
   endif
   
enddo 

DO i=1, spectrum%gas%n_absorb
   DO k=1, atm%n_layer
      DO l=1, atm%n_profile

         IF (spectrum%gas%type_absorb(i) == ip_h2o) THEN
            atm%gas_mix_ratio(l, k, i) = merge_fields%qv_n(k)
         ELSE IF (spectrum%gas%type_absorb(i) == ip_co2) THEN
            atm%gas_mix_ratio(l, k, i) = socrates_opt%co2_mmr
         ELSE IF (spectrum%gas%type_absorb(i) == ip_o3) THEN
            atm%gas_mix_ratio(l, k, i) = merge_fields%o3_n(k)
         ELSE IF (spectrum%gas%type_absorb(i) == ip_n2o) THEN
            atm%gas_mix_ratio(l, k, i) = socrates_opt%n2o_mmr
         ELSE IF (spectrum%gas%type_absorb(i) == ip_ch4) THEN
            atm%gas_mix_ratio(l, k, i) = socrates_opt%ch4_mmr
         ELSE IF (spectrum%gas%type_absorb(i) == ip_o2) THEN
            atm%gas_mix_ratio(l, k, i) = socrates_opt%o2_mmr
         ELSE IF (spectrum%gas%type_absorb(i) == ip_cfc11) THEN
            atm%gas_mix_ratio(l, k, i) = socrates_opt%cfc11_mmr
         ELSE IF (spectrum%gas%type_absorb(i) == ip_cfc12) THEN
            atm%gas_mix_ratio(l, k, i) = socrates_opt%cfc12_mmr 
         ELSE IF (spectrum%gas%type_absorb(i) == ip_cfc113) THEN
            atm%gas_mix_ratio(l, k, i) = socrates_opt%cfc113_mmr 
         ELSE IF (spectrum%gas%type_absorb(i) == ip_cfc114) THEN
            atm%gas_mix_ratio(l, k, i) = socrates_opt%cfc114_mmr
         ELSE IF (spectrum%gas%type_absorb(i) == ip_hcfc22) THEN
            atm%gas_mix_ratio(l, k, i) = socrates_opt%hcfc22_mmr
         ELSE IF (spectrum%gas%type_absorb(i) == ip_hfc125) THEN
            atm%gas_mix_ratio(l, k, i) = socrates_opt%hfc125_mmr
         ELSE IF (spectrum%gas%type_absorb(i) == ip_hfc134a) THEN
            atm%gas_mix_ratio(l, k, i) = socrates_opt%hfc134a_mmr
         ELSE
           atm%gas_mix_ratio(l, k, i) = 0.0 
         END IF
      END DO
   END DO
END DO

!print *, 'vapour, co2, o3, n2o, ch4'

!DO k=1, atm%n_layer
!   print *, k, atm%gas_mix_ratio(1, k, 1),atm%gas_mix_ratio(1, k, 2), atm%gas_mix_ratio(1, k, 3), &
!        atm%gas_mix_ratio(1, k, 4), atm%gas_mix_ratio(1, k, 5)
!enddo

END SUBROUTINE set_atm
