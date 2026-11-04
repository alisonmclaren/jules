!******************************COPYRIGHT****************************************
! (c) British Antarctic Survey 2025.
! All rights reserved.
!
! This routine has been licensed to the other JULES partners for use and
! distribution under the JULES collaboration agreement, subject to the terms and
! conditions set out therein.
!
! [Met Office Ref SC0237]
!******************************COPYRIGHT****************************************

MODULE sf_flux_land_nofrac_wtrac_mod

!-----------------------------------------------------------------------------
! Description:
!   Calculate water tracer surface fluxes over land assuming no fractionation
!   occurs for isotopes.
!
!   Unlike normal water, these are calculated separately for the four
!   different evaporative sources:
!   snow, canopy intercepted water, bare soil/transpiration from veg and
!   open water (i.e. lakes).
!
!   This is necessary due to evaporative sources having different water
!   tracer concentrations.
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

CHARACTER(LEN=*), PARAMETER, PRIVATE ::                                        &
                                 ModuleName='SF_FLUX_LAND_NOFRAC_WTRAC_MOD'

CONTAINS

SUBROUTINE sf_flux_land_nofrac_wtrac (i_wt,                                    &
 land_pts, surft_pts, surft_index, evap_srce_frac, q_elev, q_elev_wtrac,       &
 snow_surft, snow_surft_wtrac, smc_soilt, smc_soilt_wtrac,                     &
 canopy, canopy_wtrac, fqw_evapsrce_norm_wtrac, fqw_evapsrce_wtrac,            &
 fqw_1_wtrac)

USE jules_water_tracers_mod, ONLY: n_evap_srce, wtrac_lake_fixed_ratio,        &
                                    evap_srce_snow, evap_srce_cano,            &
                                    evap_srce_soil, evap_srce_lake,            &
                                    wtrac_calc_ratio_fn_jules,                &
                                    standard_ratio_wtrac
                                    
USE jules_water_iso_mod, ONLY: wtrac_eqm_factor_liqvap, wtrac_eqm_factor_icevap

USE parkind1, ONLY: jprb, jpim
USE yomhook,  ONLY: lhook, dr_hook

IMPLICIT NONE

INTEGER, INTENT(IN) ::                                                         &
 i_wt                                                                          &
                       ! IN Water tracer index
,land_pts                                                                      &
                       ! IN Total number of land points.
,surft_pts                                                                     &
                       ! IN Number of tile points.
,surft_index(land_pts)
                       ! IN Index of tile points.

REAL(KIND=real_jlslsm), INTENT(IN) ::                                          &
 evap_srce_frac(land_pts,n_evap_srce)                                          &
                       ! IN Evaporative source fractional contribution to
                       !    total evaporation from tile
,q_elev(land_pts)                                                              &
                       ! IN Specific humidity at elevated height
,q_elev_wtrac(land_pts)                                                        &
                       ! IN Water tracer specific humidity at elevated
                       !    height
,snow_surft(land_pts)                                                          &
                       ! IN Lying snow
,snow_surft_wtrac(land_pts)                                                    &
                       ! IN Water tracer content of lying snow
,smc_soilt(land_pts)                                                           &
                       ! IN Available water content in soil layers
,smc_soilt_wtrac(land_pts)                                                     &
                       ! IN Available water tracer content in soil layers
,canopy(land_pts)                                                              &
                       ! IN Surface/canopy water for snow-free land tiles
,canopy_wtrac(land_pts)
                       ! IN Water tracer content in canopy

REAL(KIND=real_jlslsm), INTENT(IN) ::                                          &
 fqw_evapsrce_norm_wtrac(land_pts,n_evap_srce)
                       ! IN Normal water tracer surface flux per evaporative
                       !     source

REAL(KIND=real_jlslsm), INTENT(IN OUT) ::                                      &
 fqw_evapsrce_wtrac(land_pts,n_evap_srce)                                      &
                       ! IN OUT Water tracer surface flux per evaporative
                       !        source
,fqw_1_wtrac(land_pts)
                       ! IN OUT Water tracer surface flux per tile

! Local variables

INTEGER ::                                                                     &
 k                                                                             &
             ! Loop counter (tile field index).
,l                                                                             &
             ! Loop counter (land point field index).
,n_evap
             ! Loop counter for evaporative sources

REAL(KIND=real_jlslsm) ::                                                      &
 r_wtrac,                                                                      &
             ! Ratio of water tracer to water for evaporative source
 atm_r_wtrac
             ! Ratio of water tracer to water for atmosphere

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='SF_FLUX_LAND_NOFRAC_WTRAC'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

!$OMP PARALLEL DO IF(surft_pts > 1)                                            &
!$OMP SCHEDULE(STATIC)                                                         &
!$OMP DEFAULT(NONE)                                                            &
!$OMP PRIVATE(l,k,n_evap,r_wtrac,atm_r_wtrac)                                  &
!$OMP SHARED(surft_pts,surft_index,i_wt,snow_surft_wtrac,snow_surft,           &
!$OMP        smc_soilt_wtrac,smc_soilt,canopy_wtrac,canopy,                    &
!$OMP        q_elev,q_elev_wtrac,fqw_evapsrce_wtrac,fqw_evapsrce_norm_wtrac,   &
!$OMP        fqw_1_wtrac,evap_srce_frac,wtrac_lake_fixed_ratio)
DO k = 1,surft_pts
  l = surft_index(k)

  ! Calculate atmosperic ratio
  atm_r_wtrac = wtrac_calc_ratio_fn_jules(i_wt,q_elev_wtrac(l),q_elev(l))
  
  DO n_evap = 1, n_evap_srce
    SELECT CASE(n_evap)

    CASE (evap_srce_snow)
      ! Snow sublimation

      ! Set ratio of water tracer to water for source
      r_wtrac = wtrac_calc_ratio_fn_jules(i_wt,snow_surft_wtrac(l),            &
                                               snow_surft(l))
      
    CASE (evap_srce_soil)
      ! Soil evapotranspiration

      ! Set ratio of water tracer to water for source
      r_wtrac = wtrac_calc_ratio_fn_jules(i_wt,smc_soilt_wtrac(l),smc_soilt(l))

    CASE (evap_srce_cano)
      ! Canopy intercepted water evaporation

      ! Set ratio of water tracer to water for source
      r_wtrac = wtrac_calc_ratio_fn_jules(i_wt, canopy_wtrac(l), canopy(l))

    CASE (evap_srce_lake)
      ! Lake evaporation

      ! Set ratio of water tracer to water for source
      r_wtrac = wtrac_lake_fixed_ratio(i_wt)

    END SELECT

    ! Set surface flux to scaled normal water flux

    IF (fqw_evapsrce_norm_wtrac(l,n_evap) > 0.0) THEN
      ! Evaporation - Scale flux using surface source
      fqw_evapsrce_wtrac(l,n_evap) =                                           &
                            r_wtrac * fqw_evapsrce_norm_wtrac(l,n_evap)

    ELSE 
      ! Condensation - Scale flux using atmosphere ratio
      fqw_evapsrce_wtrac(l,n_evap) =                                           &
                          atm_r_wtrac * fqw_evapsrce_norm_wtrac(l,n_evap)
    END IF

    ! Calculate total tile flux
    fqw_1_wtrac(l) = fqw_1_wtrac(l) + fqw_evapsrce_wtrac(l,n_evap)

  END DO  ! evaporation sources

END DO  ! tile points
!$OMP END PARALLEL DO

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

RETURN
END SUBROUTINE sf_flux_land_nofrac_wtrac

END MODULE sf_flux_land_nofrac_wtrac_mod
