!******************************COPYRIGHT****************************************
! (c) British Antarctic Survey and University of Bristol, 2023.
! All rights reserved.
!
! This routine has been licensed to the other JULES partners for use and
! distribution under the JULES collaboration agreement, subject to the terms and
! conditions set out therein.
!
! [Met Office Ref SC0237]
!******************************COPYRIGHT****************************************

MODULE sf_ssi_imp_wtrac_mod

!-----------------------------------------------------------------------------
! Description:
!   Calculate final water tracer surface fluxes over sea and sea ice.
!
! Method:
!   Update the explicit water tracer surface fluxes to account for
!   changes in the normal water flux due to:
!   a) the implicit calculation (in sf_im_pt2)
!   b) any change in temperature caused by surface melt (sf_melt) (only
!      applicable for sea ice).
!
!   The water tracer surface flux is updated by an increment which is simply
!   the change in the normal water flux scaled by the source water tracer
!   ratio.
!
! Code Owner: Please refer to ModuleLeaders.txt
! This file belongs in Surface
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to JULES coding standards v1.
!
!-----------------------------------------------------------------------------

USE um_types, ONLY: real_jlslsm

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='SF_SSI_IMP_WTRAC_MOD'

CONTAINS

SUBROUTINE sf_ssi_imp_wtrac (i_wt,                                             &
 ssi_pts, sea_pts, ssi_index, sea_index, std_ratio, tile_frac, fssi, fqw,      &
 qw_1, q_sfc_ratio, qw_1_wtrac,                                                &
 fqw_wtrac, fqw_norm_wtrac, fqw_1_wtrac, dfqw_1_wtrac, e_wtrac)

USE atm_fields_bounds_mod,    ONLY: tdims
USE theta_field_sizes,        ONLY: t_i_length
USE jules_water_tracers_mod,  ONLY: wtrac_calc_ratio_fn_jules

USE parkind1, ONLY: jprb, jpim
USE yomhook, ONLY: lhook, dr_hook

IMPLICIT NONE

INTEGER, INTENT(IN) ::                                                         &
 i_wt                                                                          &
                       ! IN Water tracer number
,ssi_pts                                                                       &
                       ! IN Total number of sea or sea-ice points
,sea_pts                                                                       &
                       ! IN Number of sea points.
,ssi_index(ssi_pts)                                                            &
                       ! IN Index of sea and sea-ice points.
,sea_index(ssi_pts)
                       ! IN Index of sea points.

REAL(KIND=real_jlslsm), INTENT(IN) ::                                          &
 std_ratio             ! IN Standard water tracer to water ratio

REAL(KIND=real_jlslsm), INTENT(IN) ::                                          &
 tile_frac(ssi_pts)                                                            &
                       ! IN 'Tile' fraction so leads fraction if ocean call
                       !    or sea ice category fraction if sea ice call
,fssi(ssi_pts)                                                                 &
                       ! IN Non-land fraction of grid box

,fqw(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)                      &
                       ! IN Final normal water surface flux at end of implicit
                       !    scheme
,qw_1(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)                     &
                       ! IN Atmospheric water content

,q_sfc_ratio(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)              &
                       ! IN Water tracer ratio at surface (sea or sea ice)

,qw_1_wtrac(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)               &
                       ! IN Atmospheric water tracer content

,fqw_wtrac(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)                &
                       ! IN Water tracer surface flux for 'tile'
                       !    (i.e. open ocean or sea ice category)
                       !    Calculated in explicit scheme
                       
,fqw_norm_wtrac(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)
                       ! IN NORMAL water tracer surface flux for
                       !   'tile' (i.e. open ocean or sea ice category)
                       !   Calculated in explicit scheme

REAL(KIND=real_jlslsm), INTENT(IN OUT) ::                                      &
 fqw_1_wtrac(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)              &
                       ! IN OUT GBM water tracer surface flux
,dfqw_1_wtrac(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)
                       ! IN OUT GBM increment to water tracer surface flux

REAL(KIND=real_jlslsm), INTENT(OUT) ::                                      &
 e_wtrac(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)
                       ! OUT Water tracer evaporation for 'tile' (i.e.
                       !        open ocean or sea ice category)
                       ! Flux corrected for implicit correction to
                       !       normal water flux


! Local variables

INTEGER ::                                                                     &
 i,j                                                                           &
             ! Horizontal field index
,k                                                                             &
             ! Point counter
,l
             ! Sea or sea ice field index

REAL(KIND=real_jlslsm) ::                                                      &
 e_diff                                                                        &
             ! Difference between the implicit and explicit surface flux
,qw_atm_ratio
             ! Water tracer ratio of the atmospheric water content

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='SF_SSI_IMP_WTRAC'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

!$OMP PARALLEL DO IF(sea_pts > 1)                                              &
!$OMP SCHEDULE(STATIC) DEFAULT(NONE)                                           &
!$OMP PRIVATE(l,k,i,j,e_diff,qw_atm_ratio)                                     &
!$OMP SHARED(sea_pts,sea_index,ssi_index,t_i_length,fqw_1_wtrac,fqw,           &
!$OMP        e_wtrac,tile_frac,fssi,std_ratio,dfqw_1_wtrac,                    &
!$OMP        i_wt,qw_1,qw_1_wtrac,q_sfc_ratio,fqw_wtrac,fqw_norm_wtrac)
DO k = 1,sea_pts
  l = sea_index(k)
  j=(ssi_index(l) - 1) / t_i_length + 1
  i = ssi_index(l) - (j-1) * t_i_length

  ! Calculate difference between normal water evaporation (calculated in the
  ! explicit scheme) and the final actual water evaporation field 
  e_diff = fqw(i,j) - fqw_norm_wtrac(i,j)

  ! Correct water tracer flux for the additional water flux (assume no
  ! fractionation here)  
  IF (e_diff < 0.0) THEN
    ! Source is the atmosphere
    qw_atm_ratio = wtrac_calc_ratio_fn_jules(i_wt,                             &
                                             qw_1_wtrac(i,j),qw_1(i,j))
    e_wtrac(i,j) = fqw_wtrac(i,j) + qw_atm_ratio * e_diff

  ELSE IF (e_diff > 0.0) THEN
    ! Source is the sea or sea ice
    e_wtrac(i,j) = fqw_wtrac(i,j) + q_sfc_ratio(i,j) * e_diff

  END IF

  ! Add new flux to GBM fields
  ! (where tile_frac = leads fraction for ocean call
  !                  = sea ice fraction for sea ice call)
  dfqw_1_wtrac(i,j) = dfqw_1_wtrac(i,j)                                        &
                           + (e_wtrac(i,j) - fqw_wtrac(i,j))                   &
                              * tile_frac(l) * fssi(l)
  fqw_1_wtrac(i,j)  = fqw_1_wtrac(i,j)                                         &
                           +  e_wtrac(i,j) * tile_frac(l) * fssi(l)

END DO
!$OMP END PARALLEL DO

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

RETURN
END SUBROUTINE sf_ssi_imp_wtrac

END MODULE sf_ssi_imp_wtrac_mod
