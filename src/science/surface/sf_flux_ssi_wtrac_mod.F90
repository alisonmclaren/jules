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

MODULE sf_flux_ssi_wtrac_mod

!-----------------------------------------------------------------------------
! Description:
!   Calculate water tracer surface fluxes over sea and sea ice.
!   This includes fractionation for water isotopes.
!
! Code Owner: Please refer to ModuleLeaders.txt
! This file belongs in Surface
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to JULES coding standards v1.
!
!-----------------------------------------------------------------------------

USE um_types, ONLY: real_jlslsm, real_32

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='SF_FLUX_SSI_WTRAC_MOD'

CONTAINS

SUBROUTINE sf_flux_ssi_wtrac (i_wt, ssi_pts, sea_pts, sea_index,               &
 sea_point,  salinityfactor, alpha1, rhokh_1, qstar, tstar,                    &
 dtstar, resf, v_s, z1, r_sfc_wtrac, ssi_type, q_elev_wtrac,                   &
 fqw_1_wtrac)

USE water_tracers_mod,   ONLY: wtrac_info
USE jules_water_iso_mod, ONLY: wtrac_eqm_factor_liqvap,                        &
                               wtrac_eqm_factor_icevap,                        &
                               wtrac_kin_factor_method,                        &
                               wtrac_kin_factor_kj,wtrac_kin_factor_const
USE Ereport_mod,             ONLY: ereport

USE parkind1, ONLY: jprb, jpim
USE yomhook, ONLY: lhook, dr_hook


IMPLICIT NONE

INTEGER, INTENT(IN) ::                                                         &
 i_wt                                                                          &
                       ! IN Water tracer index
,ssi_pts                                                                       &
                       ! IN Total number of sea or sea-ice points
,sea_pts                                                                       &
                       ! IN Number of sea points.
,sea_index(ssi_pts)
                       ! IN Index of sea points.


REAL(KIND=real_jlslsm), INTENT(IN) ::                                          &
 sea_point                                                                     &
                       ! IN If sea_point = 0.0, then flux is corrected for
                       ! change in sfc temperature. If = 1.0, it is not.
,salinityfactor
                       ! IN Factor allowing for the effect of the salinity
                       ! of sea water on the evaporative flux.

REAL(KIND=real_jlslsm), INTENT(IN) ::                                          &
 alpha1(ssi_pts)                                                               &
                       ! IN Mean gradient of saturated specific humidity with
                       !    respect to temperature between the bottom model
                       !    layer and tile surfaces
,rhokh_1(ssi_pts)                                                              &
                       ! IN Surface exchange coefficients for land tiles
,qstar(ssi_pts)                                                                &
                       ! IN Surface saturated specific humidity
,tstar(ssi_pts)                                                                &
                       ! IN Surface temperature

,dtstar(ssi_pts)                                                               &
                       ! IN Change in TSTAR over timestep for land tiles
,resf(ssi_pts)                                                                 &
                       ! IN Resistance factor for normal water evaporation
,v_s(ssi_pts)                                                                  &
                       ! IN u_star
,z1(ssi_pts)                                                                   &
                       ! IN Height of lowest atmospheric level
,q_elev_wtrac(ssi_pts)                                                         &
                       ! IN Water tracer specific humidity at elevated
                       !    height
,r_sfc_wtrac(ssi_pts)
                       ! IN Water tracer to water ratio of sea or sea ice
                       !    surface

CHARACTER(LEN=3), INTENT(IN) :: ssi_type   ! IN 'sea' or 'ice'

REAL(KIND=real_jlslsm), INTENT(IN OUT) ::                                      &
 fqw_1_wtrac(ssi_pts)
                       ! IN OUT Water tracer surface flux per tile

! Local variables

INTEGER ::                                                                     &
 k                                                                             &
             ! Point counter
,l
             ! Sea or sea ice field index

REAL(KIND=real_jlslsm) ::                                                      &
 eqm_factor                                                                    &
             ! Equilibrium fractionation factor
,eqm_factor_32                                                                 &
             ! Equilibrium fractionation factor at real32 precision
,kin_factor                                                                    &
             ! Kinetic fractionation factor
,alpha1_wtrac                                                                  &
             ! alpha1 for water tracers
,epot                                                                          &
             ! Potential evaporation for water tracers
,resf_wtrac
             ! Resistance factor (including kinetic
             !  fractionation for istopes)

INTEGER :: ErrorStatus        ! Error code
CHARACTER(LEN=256) :: message

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='SF_FLUX_SSI_WTRAC'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

ErrorStatus=0

!$OMP PARALLEL DO IF(sea_pts > 1)                                              &
!$OMP SCHEDULE(STATIC)                                                         &
!$OMP DEFAULT(NONE)                                                            &
!$OMP PRIVATE(l,k,eqm_factor,eqm_factor_32,alpha1_wtrac,epot,kin_factor,       &
!$OMP         resf_wtrac,ErrorStatus,message)                                  &
!$OMP SHARED(sea_pts,sea_index,ssi_type,resf,r_sfc_wtrac, wtrac_info,          &
!$OMP        alpha1,rhokh_1,salinityfactor,qstar,q_elev_wtrac,fqw_1_wtrac,     &
!$OMP        tstar,dtstar,sea_point,i_wt,v_s,z1)
DO k = 1,sea_pts
  l = sea_index(k)

  IF (wtrac_info(i_wt)%is_iso) THEN

    SELECT CASE(ssi_type)

    CASE ('sea')
      ! Open ocean evaporation

      ! Get isotopic equilibrium fractionation factor
      ! Calculate at real32 precision and convert to jlslsm
      ! (This saves having to have different versions of the function at
      !  different precisions)
      eqm_factor_32 = wtrac_eqm_factor_liqvap(i_wt,                            &
                                      REAL(tstar(l),KIND=real_32),'sfc_evap')
      eqm_factor = REAL(eqm_factor_32, KIND=real_jlslsm) 

      ! Get isotopic kinetic fractionation factor
      SELECT CASE(wtrac_kin_factor_method)
      CASE('K&J79')
         kin_factor = wtrac_kin_factor_kj(i_wt,v_s(l),z1(l))
      CASE('const')
         kin_factor = wtrac_kin_factor_const(i_wt)
      CASE DEFAULT
        ! Abort model
        ErrorStatus = 100
        message = 'Water isotope kinetic frac coeff not set'
        CALL Ereport(RoutineName, ErrorStatus, message)
      END SELECT
      
    CASE ('ice')
      ! Sea ice sublimation

      ! Set isotopic equilibrium fractionation factor
      ! (Temporarily set to 1 until isotopes are coded)
      eqm_factor = 1.0

      ! Set resistance factor
      ! (This will include kinetic fractionation for isotopes)
      kin_factor = 1.0

    END SELECT

  ELSE
    ! No fractionation
    eqm_factor = 1.0
    kin_factor = 1.0
  END IF

  ! Set total resistance factor 
  resf_wtrac = resf(l)*kin_factor

  ! Calculate water tracer equivalent of alpha1
  alpha1_wtrac = (r_sfc_wtrac(l)/eqm_factor) * alpha1(l)

  ! Calculate evaporation potential
  epot = rhokh_1(l) *                                                          &
             ( (r_sfc_wtrac(l)/eqm_factor) * salinityfactor * qstar(l)         &
                                - q_elev_wtrac(l) )

  ! Calculate explicit water tracer surface flux
  fqw_1_wtrac(l) = resf_wtrac * epot

  ! Correction to surface isotopic moisture fluxes due to change in
  ! surface temperature
  fqw_1_wtrac(l) = fqw_1_wtrac(l) + resf_wtrac                                 &
                          * rhokh_1(l) * alpha1_wtrac * dtstar(l)              &
                          * (1.0 - sea_point)

END DO
!$OMP END PARALLEL DO

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

RETURN
END SUBROUTINE sf_flux_ssi_wtrac
END MODULE sf_flux_ssi_wtrac_mod
