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

MODULE sf_resist_wtrac_mod

!-----------------------------------------------------------------------------
! Description:
!   Set resistance factors for water tracers (which could include the kinetic
!   fractionation factors for water isotopes in the future).
!
!   Resistance factors are set for the four different evaporative sources:
!   snow, canopy intercepted water, bare soil/transpiration from vegetation
!   and open water (i.e. lakes).
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

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='SF_RESIST_WTRAC_MOD'

CONTAINS

SUBROUTINE sf_resist_wtrac (i_wt, n, land_pts, surft_pts, surft_index,         &
 snowdep_surft, snow_surft, resfs, dq, catch, evap_srce_frac, cansnowtile,     &
 resf_wtrac)

USE jules_water_tracers_mod, ONLY: n_evap_srce, evap_srce_snow,                &
                                   evap_srce_cano, evap_srce_soil,             &
                                   evap_srce_lake
USE jules_vegetation_mod,    ONLY: can_model
USE jules_surface_types_mod, ONLY: lake
USE jules_science_fixes_mod, ONLY: l_fix_neg_snow
USE water_constants_mod,     ONLY: tm, rho_ice

USE parkind1, ONLY: jprb, jpim
USE yomhook, ONLY: lhook, dr_hook

IMPLICIT NONE

INTEGER, INTENT(IN) ::                                                         &
 i_wt                                                                          &
                       ! IN Water tracer index
,n                                                                             &
                       ! IN Tile index
,land_pts                                                                      &
                       ! IN Total number of land points.
,surft_pts                                                                     &
                       ! IN Number of tile points.
,surft_index(land_pts)
                       ! IN Index of tile points.

REAL(KIND=real_jlslsm), INTENT(IN) ::                                          &
 snowdep_surft(land_pts)                                                       &
                       ! IN Snow depth (m)
,snow_surft(land_pts)                                                          &
                       ! IN Lying snow on tiles (kg/m2)
,resfs(land_pts)                                                               &
                       ! IN Combined soil, stomatal and aerodynamic
                       !     resistance factor for fraction 1-fracaero.
,dq(land_pts)                                                                  &
                       ! IN Specific humidity difference between lowest
                       !    atmospheric level and the surface (Q1 - Q*)
,catch(land_pts)                                                               &
                       ! IN Surface/canopy water capacity (kg/m2)
,evap_srce_frac(land_pts,n_evap_srce)
                       ! IN Evaporative source fractional contribution to
                       !     total evaporation from tile

LOGICAL, INTENT (IN) ::                                                        &
 cansnowtile           ! IN Switch for pft canopy snow model

REAL(KIND=real_jlslsm), INTENT(OUT) ::                                         &
 resf_wtrac(land_pts,n_evap_srce)
                       ! OUT Resitance factor for water tracers for individual
                       !     evaporative source types

! Local variables

INTEGER ::                                                                     &
 k                                                                             &
             ! Loop counter (tile field index)
,l                                                                             &
             ! Loop counter (land point field index)
,n_evap
             ! Loop counter (evaporative sources index)

REAL(KIND=real_jlslsm) ::                                                      &
 snow_test
             ! Quantity used to test for presence of snow

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='SF_RESIST_WTRAC'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

!$OMP PARALLEL DO IF(surft_pts > 1)                                            &
!$OMP SCHEDULE(STATIC)                                                         &
!$OMP DEFAULT(NONE)                                                            &
!$OMP PRIVATE(l,k,n_evap,snow_test)                                            &
!$OMP SHARED(surft_pts,surft_index,i_wt,evap_srce_frac,dq,resf_wtrac,          &
!$OMP        snowdep_surft,snow_surft,catch,can_model,cansnowtile,resfs,n,     &
!$OMP        lake,l_fix_neg_snow)
DO k = 1,surft_pts
  l = surft_index(k)

  ! Calculate snow quantity used in scheme to test for presence of snow
  IF (l_fix_neg_snow) THEN
    snow_test = MAX(snowdep_surft(l), snow_surft(l) / rho_ice)
  ELSE
    snow_test = snowdep_surft(l)
  END IF

  ! Set resistance factors for water tracers

  DO n_evap = 1, n_evap_srce
    ! Initialise to 1.0 (i.e. aerodynamic resistance only)
    resf_wtrac(l,n_evap) = 1.0
  END DO

  ! Snow sublimation
  IF (snow_test > 0.0) THEN
    IF (can_model == 4 .AND. cansnowtile) THEN
      ! Canopy snow model and snow sublimation has resistance factor < 1
      resf_wtrac(l,evap_srce_snow) = resfs(l)
    ELSE
      ! For all non-vegetative tiles + vegetative tiles when no canopy snow
      !  model
      resf_wtrac(l,evap_srce_snow) = 1.0
    END IF
  END IF

  ! Soil evapotranspiration (Bare soil evaporation + Transpiration)
  ! Evaporation is limited by soil moisture availability and stomatal
  ! resistance.
  IF (dq(l) < 0.0) THEN
    ! Evaporation (with or without snow)
    resf_wtrac(l,evap_srce_soil) = resfs(l) * 1.0
  END IF

  ! If using bug version, then there is additional resistance if there is snow
  ! regardless of the direction of the flux
  IF (.NOT. l_fix_neg_snow) THEN
    IF ( .NOT. (can_model == 4 .AND. cansnowtile)) THEN
      IF (snow_test > 0.0 .AND. dq(l) > 0.0) THEN
        resf_wtrac(l,evap_srce_soil) = resfs(l) * 1.0
      END IF
    END IF
  END IF

  ! Canopy intercepted water
  resf_wtrac(l,evap_srce_cano) = 1.0

  ! Lake evaporation
  resf_wtrac(l,evap_srce_lake) = 1.0

END DO
!$OMP END PARALLEL DO

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

RETURN
END SUBROUTINE sf_resist_wtrac
END MODULE sf_resist_wtrac_mod
