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

MODULE sf_evap_srce_frac_mod

!-----------------------------------------------------------------------------
! Description:
!   Calculate the fractional contribution of each evaporative source to the
!   total tile evaporative flux.
!
!   The different evaporative sources are: snow, canopy intercepted water,
!   bare soil/transpiration from vegetation and open water (i.e. lakes).
!
!   This code is based on how the total tile resistance factors are
!   set in sf_resist.
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

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='SF_EVAP_SRCE_FRAC_MOD'

CONTAINS

SUBROUTINE sf_evap_srce_frac (land_pts, surft_pts, surft_index, flake,         &
 fracaero_t, fracaero_s, snowdep_surft, snow_surft, tstar_surft, dq, catch,    &
 cansnowtile, evap_srce_frac)

USE jules_water_tracers_mod, ONLY: n_evap_srce, evap_srce_snow,                &
                                   evap_srce_cano, evap_srce_soil,             &
                                   evap_srce_lake
USE jules_vegetation_mod,    ONLY: can_model
USE jules_science_fixes_mod, ONLY: l_fix_neg_snow
USE water_constants_mod,     ONLY: tm, rho_ice

USE parkind1, ONLY: jprb, jpim
USE yomhook, ONLY: lhook, dr_hook

IMPLICIT NONE

INTEGER, INTENT(IN) ::                                                         &
 land_pts                                                                      &
                       ! IN Total number of land points.
,surft_pts                                                                     &
                       ! IN Number of tile points.
,surft_index(land_pts)
                       ! IN Index of tile points.

REAL(KIND=real_jlslsm), INTENT(IN) ::                                          &
 flake(land_pts)                                                               &
                       ! IN Lake fraction
,fracaero_t(land_pts)                                                          &
                       ! IN Fraction of surface moisture flux with only
                       !    aerodynamic resistance
,fracaero_s(land_pts)                                                          &
                       ! IN Fraction of surface moisture flux with only
                       !    aerodynamic resistance over the frozen part of the
                       !    surface only
,snowdep_surft(land_pts)                                                       &
                       ! IN Snow depth (m)
,snow_surft(land_pts)                                                          &
                       ! IN Lying snow on tiles (kg/m2)

,tstar_surft(land_pts)                                                         &
                       ! IN Surface temperature (K)
,dq(land_pts)                                                                  &
                       ! IN Specific humidity difference between lowest
                       !    atmospheric level and the surface (Q1 - Q*)
,catch(land_pts)
                       ! IN Surface/canopy water capacity (km/m2)

LOGICAL, INTENT (IN) ::                                                        &
 cansnowtile           ! IN Switch for pft canopy snow model

REAL(KIND=real_jlslsm), INTENT(OUT) ::                                         &
 evap_srce_frac(land_pts,n_evap_srce)
                       ! OUT Evaporative source fractional contribution to
                       !     total evaporation from tile

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

CHARACTER(LEN=*), PARAMETER :: RoutineName='SF_EVAP_SRCE_FRAC'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

!$OMP PARALLEL DO IF(surft_pts > 1)                                            &
!$OMP SCHEDULE(STATIC)                                                         &
!$OMP DEFAULT(NONE)                                                            &
!$OMP PRIVATE(l,k,n_evap,snow_test)                                            &
!$OMP SHARED(surft_pts,surft_index,evap_srce_frac,                             &
!$OMP        flake,fracaero_t,fracaero_s,dq,snowdep_surft,snow_surft,          &
!$OMP        can_model,tstar_surft,cansnowtile,catch,l_fix_neg_snow)
DO k = 1,surft_pts
  l = surft_index(k)

  ! Calculate snow quantity used in scheme to test for presence of snow
  IF (l_fix_neg_snow) THEN
    snow_test = MAX(snowdep_surft(l), snow_surft(l) / rho_ice)
  ELSE
    snow_test = snowdep_surft(l)
  END IF

  ! Set tile fractions for each evaporative source

  DO n_evap = 1, n_evap_srce
    evap_srce_frac(l,n_evap) = 0.0
  END DO

  IF (dq(l) >= 0.0) THEN
    ! Condensation - only aerodynamic resistance for downward fluxes

    IF (tstar_surft(l) > tm) THEN
        
      IF (catch(l) > 0.0) THEN
        ! Downward moisture flux is added to canopy water and lakes
        evap_srce_frac(l,evap_srce_cano) = 1.0 - flake(l)
      ELSE
        ! Downward moisture flux added to soil and lakes  
        evap_srce_frac(l,evap_srce_soil) = 1.0 - flake(l)
      END IF
      evap_srce_frac(l,evap_srce_lake) = flake(l)

    ELSE
      ! All downward moisture flux is added to snow
      evap_srce_frac(l,evap_srce_snow) = 1.0

    END IF

    ! If using bug version, then there is additional resistance if there is snow
    ! regardless of the direction of the flux
    ! Not clear how to split this up.....
    IF (.NOT. l_fix_neg_snow) THEN
      IF ( .NOT. (can_model == 4 .AND. cansnowtile)) THEN
        IF (snow_test > 0.0) THEN
          evap_srce_frac(l,evap_srce_snow) = fracaero_t(l) * (1.0 - flake(l))
          evap_srce_frac(l,evap_srce_soil) = (1.0-fracaero_t(l))               &
                                                       * (1.0 - flake(l))
          evap_srce_frac(l,evap_srce_lake) = flake(l)
          evap_srce_frac(l,evap_srce_cano) = 0.0
        END IF
      END IF
    END IF

  ELSE
    ! Evaporation or sublimation

    IF (l_fix_neg_snow) THEN

      ! This doesn't work for l_fix_neg_snow = F as fracareao_s not set!!!!

      evap_srce_frac(l,evap_srce_snow) = fracaero_s(l) * (1.0 - flake(l))
      evap_srce_frac(l,evap_srce_cano) = (fracaero_t(l)-fracaero_s(l))         &
                                                * (1.0 - flake(l))
      evap_srce_frac(l,evap_srce_soil) = (1.0-fracaero_t(l)) * (1.0 - flake(l))
      evap_srce_frac(l,evap_srce_lake) = flake(l)

    ELSE
      IF (snow_test > 0.0) THEN
        evap_srce_frac(l,evap_srce_snow) = fracaero_t(l) * (1.0 - flake(l))
        evap_srce_frac(l,evap_srce_soil) = (1.0-fracaero_t(l))                 &
                                                * (1.0 - flake(l))
        evap_srce_frac(l,evap_srce_lake) = flake(l)
        evap_srce_frac(l,evap_srce_cano) = 0.0
      ELSE
        evap_srce_frac(l,evap_srce_cano) = fracaero_t(l) * (1.0 - flake(l))
        evap_srce_frac(l,evap_srce_soil) = (1.0-fracaero_t(l))                &
                                                 * (1.0 - flake(l))
        evap_srce_frac(l,evap_srce_lake) = flake(l)
        evap_srce_frac(l,evap_srce_snow) = 0.0
      END IF
    END IF

  END IF

  ! There are different settings if using the canopy snow model
  IF (can_model == 4 .AND. cansnowtile .AND. snow_test > 0.0) THEN

    ! In this case, fracaero_t = 0.0
    evap_srce_frac(l,evap_srce_snow) = (1.0 - flake(l))
    evap_srce_frac(l,evap_srce_lake) = flake(l)
    evap_srce_frac(l,evap_srce_cano) = 0.0
    evap_srce_frac(l,evap_srce_soil) = 0.0

  END IF
END DO
!$OMP END PARALLEL DO

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

RETURN
END SUBROUTINE sf_evap_srce_frac
END MODULE sf_evap_srce_frac_mod
