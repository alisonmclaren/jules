!******************************COPYRIGHT****************************************
! (c) British Antarctic Survey, 2025.
! All rights reserved.
!
! This routine has been licensed to the other JULES partners for use and
! distribution under the JULES collaboration agreement, subject to the terms and
! conditions set out therein.
!
! [Met Office Ref SC0237]
!******************************COPYRIGHT****************************************

MODULE sf_flux_land_norm_mod

!-----------------------------------------------------------------------------
! Description:
!   Calculate water surface fluxes over land for use by water tracers.
!
!   Water tracers require surface fluxes to be calculated separately for the
!   four different evaporative sources:
!   snow, canopy intercepted water, bare soil/transpiration from veg and
!   open water (i.e. lakes).
!
!   This is necessary due to evaporative sources having different water
!   tracer concentrations.  (And potentially different kinetic fractionation
!   factors for isotopes, although currently no fractionation is included.)
!
!   The module also contains a subroutine to check the calculated fluxes
!   against the water fluxes.
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

! Maximum difference between water and water tracer fqw_1 which trigger
! warning or abort
REAL(KIND=real_jlslsm), PARAMETER :: max_diff_warning = 1.0e-16
REAL(KIND=real_jlslsm), PARAMETER :: max_diff_abort   = 1.0e-12

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='SF_FLUX_LAND_NORM_MOD'

CONTAINS

SUBROUTINE sf_flux_land_norm (                                                 &
 land_pts, surft_pts, surft_index, alpha1, rhokh_1, qstar,                     &
 tstar, dtstar, resf_wtrac, evap_srce_frac, q_elev, fqw_evapsrce)

USE jules_water_tracers_mod, ONLY: n_evap_srce
                                    
USE parkind1, ONLY: jprb, jpim
USE yomhook,  ONLY: lhook, dr_hook

IMPLICIT NONE

INTEGER, INTENT(IN) ::                                                         &
 land_pts                                                                      &
                       ! IN Total number of land points.
,surft_pts                                                                     &
                       ! IN Number of tile points.
,surft_index(land_pts)
                       ! IN Index of tile points.

REAL(KIND=real_jlslsm), INTENT(IN) ::                                          &
 alpha1(land_pts)                                                              &
                       ! IN Mean gradient of saturated specific humidity with
                       !    respect to temperature between the bottom model
                       !    layer and tile surfaces
,rhokh_1(land_pts)                                                             &
                       ! IN Surface exchange coefficients for land tiles
,qstar(land_pts)                                                               &
                       ! IN Surface saturated specific humidity
,tstar(land_pts)                                                               &
                       ! IN Surface temperature
,dtstar(land_pts)                                                              &
                       ! IN Change in TSTAR over timestep for land tiles
,resf_wtrac(land_pts,n_evap_srce)                                              &
                       ! IN Resitance factor for water tracers for individual
                       !     evaporative source types
                       !     The normal water tracer values are used here
,evap_srce_frac(land_pts,n_evap_srce)                                          &
                       ! IN Evaporative source fractional contribution to
                       !    total evaporation from tile
,q_elev(land_pts)
                       ! IN Specific humidity at elevated height

REAL(KIND=real_jlslsm), INTENT(IN OUT) ::                                      &
 fqw_evapsrce(land_pts,n_evap_srce)
                       ! IN OUT Water surface flux per evaporative
                       !        source

! Local variables

INTEGER ::                                                                     &
 k                                                                             &
             ! Loop counter (tile field index).
,l                                                                             &
             ! Loop counter (land point field index).
,n_evap
             ! Loop counter for evaporative sources

REAL(KIND=real_jlslsm) ::                                                      &
 epot
             ! Potential evaporation for water tracers

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='SF_FLUX_LAND_NORM'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

!$OMP PARALLEL DO IF(surft_pts > 1)                                            &
!$OMP SCHEDULE(STATIC)                                                         &
!$OMP DEFAULT(NONE)                                                            &
!$OMP PRIVATE(l,k,n_evap,epot)                                                 &
!$OMP SHARED(tstar,surft_pts,surft_index,alpha1,rhokh_1,qstar,q_elev,          &
!$OMP        resf_wtrac,fqw_evapsrce,dtstar,evap_srce_frac)
DO k = 1,surft_pts
  l = surft_index(k)

  DO n_evap = 1, n_evap_srce

    ! Calculate evaporation potential
    epot = rhokh_1(l) * ( qstar(l) - q_elev(l) )

    ! Calculate explicit water surface flux
    fqw_evapsrce(l,n_evap) = resf_wtrac(l,n_evap) * epot

    ! Correction to surface flux due to change in surface temperature
    fqw_evapsrce(l,n_evap) =  fqw_evapsrce(l,n_evap)                           &
                              + resf_wtrac(l,n_evap)                           &
                                      * rhokh_1(l) * alpha1(l) * dtstar(l)

    ! Scale surface flux by evaporative source fraction to get tile mean value
    fqw_evapsrce(l,n_evap) =  fqw_evapsrce(l,n_evap)                           &
                                        * evap_srce_frac(l,n_evap)

  END DO  ! evaporation sources

END DO  ! tile points
!$OMP END PARALLEL DO

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

RETURN
END SUBROUTINE sf_flux_land_norm

!-----------------------------------------------------------------------------

SUBROUTINE sf_flux_land_check_wtrac (                                          &
 land_pts, nsurft, land_index, tile_frac, flandg,                              &
 snowdep_surft, canopy, flake, dq, fracaero, resfs, resft,                     &
 resf_wtrac, evap_srce_frac,                                                   &
 fqw_surft, fqw_1, fqw_surft_wtrac, fqw_1_wtrac, fqw_evapsrce_wtrac )

!-----------------------------------------------------------------------------
! Description:
!   Compare the 'normal' water tracer explicit surface flux with the
!   calculated water flux.  If the differences are large, then either
!   write a warning or abort model depending on size of error.
!-----------------------------------------------------------------------------

USE atm_fields_bounds_mod, ONLY: tdims
USE theta_field_sizes,     ONLY: t_i_length

USE jules_water_tracers_mod, ONLY: n_evap_srce

USE jules_print_mgr, ONLY: jules_message, jules_print, newline

USE ereport_mod, ONLY: ereport

USE parkind1, ONLY: jprb, jpim
USE yomhook,  ONLY: lhook, dr_hook

IMPLICIT NONE

INTEGER, INTENT(IN) ::                                                         &
 land_pts                                                                      &
                       ! IN Total number of land points.
,nsurft                                                                        &
                       ! IN No. of land-surface tiles
,land_index(land_pts)
                       ! IN Index of land points.

REAL(KIND=real_jlslsm), INTENT(IN) ::                                          &
 tile_frac(land_pts,nsurft)                                                    &
                       ! IN Tile fractions
,flandg(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)                   &
                       ! IN Land fraction on all tiles.
,snowdep_surft(land_pts,nsurft)                                                &
                       ! IN Snow depth
,canopy(land_pts, nsurft)                                                      &
                       ! IN Surface/canopy water for snow-free land tiles
,flake(land_pts,nsurft)                                                        &
                       ! IN Lake fraction.
,dq(land_pts,nsurft)                                                           &
                       ! IN Specific humidity difference between
                       !    surface and lowest atmospheric lev
,fracaero(land_pts,nsurft)                                                     &
                       ! IN Fraction of surface moisture flux with only
                       !    aerodynamic resistance for snow-free land tiles.
,resfs(land_pts,nsurft)                                                        &
                       ! IN Combined soil, stomatal and aerodynamic resistance
                       !     factor for fraction (1-fracaero) of snow-free land
                       !     tiles.
,resft(land_pts,nsurft)                                                        &
                       ! IN Total resistance factor.
                       !  fracaero+(1-fracaero)*resfs for snow-free land,
                       !  1 for snow.
,resf_wtrac(land_pts,n_evap_srce,nsurft)                                       &
                       ! IN Resitance factor for water tracers for individual
                       !     evaporative source types
,evap_srce_frac(land_pts,n_evap_srce,nsurft)
                       ! IN Evaporative source fractional contribution to
                       !    total evaporation from tile

REAL(KIND=real_jlslsm), INTENT(IN) ::                                          &
 fqw_surft(land_pts,nsurft)                                                    &
                       ! IN Surface FQW for land tiles
,fqw_1(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)                    &
                       ! IN GBM surface FQW
,fqw_surft_wtrac(land_pts,nsurft)                                              &
                       ! IN Water tracer surface flux per tile
,fqw_1_wtrac(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)              &
                       ! IN GBM Water tracer surface FQW
,fqw_evapsrce_wtrac(land_pts,n_evap_srce,nsurft)
                       ! IN Water tracer surface flux per evaporative
                       !        source

! Local variables

INTEGER ::                                                                     &
 i,j                                                                           &
             ! Loop counter (horizontal field index).
,l                                                                             &
             ! Loop counter (land point field index).
,n                                                                             &
             ! Loop counter (tile index)
,n_evap
             ! Loop counter for evaporative sources

INTEGER :: imax, jmax, lmax

INTEGER :: ErrorStatus

REAL(KIND=real_jlslsm) :: diff, largest_error

CHARACTER(LEN=256)     :: message                 ! Error reporting

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='SF_FLUX_LAND_CHECK_WTRAC'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! Compare water tracer flux to normal water flux and find largest error.

largest_error = 0.0

DO l = 1,land_pts
  j=(land_index(l) - 1) / t_i_length + 1
  i = land_index(l) - (j-1) * t_i_length
  diff = fqw_1(i,j) - fqw_1_wtrac(i,j)
  IF ( ABS(diff) > ABS(largest_error) ) THEN
    largest_error = diff
    imax          = i
    jmax          = j
    lmax          = l
  END IF
END DO

! If error is greater than max_diff_warning then write warning.

IF (ABS(largest_error) > max_diff_warning) THEN

  WRITE(jules_message,'(a,3g20.12,a,2(i0,1x),e14.6)')                          &
  newline// '*********WATER TRACER PROBLEM*********' //newline//               &
  'Large difference between tracer and water for explicit fqw_1 in JULES'      &
  //newline// 'after land calculations' //newline//                            &
  ' fqw_1, fqw_1_wtrac, max diff  = ' //newline,                               &
  fqw_1(imax,jmax),fqw_1_wtrac(imax,jmax), largest_error,                      &
  newline// ' at i,j,flandg = ', imax,jmax,flandg(imax,jmax)

  CALL jules_print(RoutineName, jules_message)

  ! Print out tile information

  DO n = 1, nsurft

    IF (tile_frac(lmax,n) > 0.0) THEN

      WRITE(jules_message,'(A,i3,e20.12)')                                     &
        newline// 'Individual file fluxes for n, tile_frac =',                 &
        n,tile_frac(lmax,n)
      CALL jules_print(RoutineName, jules_message)

      WRITE(jules_message,'(A,3e20.12)')                                       &
        newline// 'fqw_surft: water, wtrac, diff = ', fqw_surft(lmax,n),       &
      fqw_surft_wtrac(lmax,n), fqw_surft(lmax,n)-fqw_surft_wtrac(lmax,n)
      CALL jules_print(RoutineName, jules_message)

      WRITE(jules_message,'(A,7e20.12)')                                       &
       'snow,canopy,flake,fracaero,dq,resft,resfs = ',                         &
        snowdep_surft(lmax,n), canopy(lmax,n),                                 &
        flake(lmax,n),fracaero(lmax,n),dq(lmax,n),resfs(lmax,n),resft(lmax,n)
      CALL jules_print(RoutineName, jules_message)

      WRITE(jules_message,'(A,4e20.12)')                                       &
       'fqw_evapsrce_wtrac=', (fqw_evapsrce_wtrac(lmax,n_evap,n), n_evap=1,4)
      CALL jules_print(RoutineName, jules_message)

      WRITE(jules_message,'(A,4e20.12)')                                       &
        'resf_wtrac=', (resf_wtrac(lmax,n_evap,n), n_evap=1,4)
      CALL jules_print(RoutineName, jules_message)

      WRITE(jules_message,'(A,4e20.12)')                                       &
        'evap_srce_frac=', (evap_srce_frac(lmax,n_evap,n),n_evap=1,4)
      CALL jules_print(RoutineName, jules_message)

    END IF ! tile_frac > 0
  END DO   ! n

  ! Either write warning or abort model
  IF (ABS(largest_error) > max_diff_abort) THEN
    ! Abort model
    ErrorStatus = 100
    message = 'Problem with water tracers.'
  ELSE
    ! Output warning
    ErrorStatus = -100
    message = 'Problem with water tracers.'
  END IF

  CALL ereport(RoutineName, ErrorStatus, message)

END IF  ! Test on largest error size

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN

END SUBROUTINE sf_flux_land_check_wtrac

END MODULE sf_flux_land_norm_mod
