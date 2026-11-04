!******************************COPYRIGHT****************************************
! (c) British Antarctic Survey and University of Bristol, 2024.
! All rights reserved.
!
! This routine has been licensed to the other JULES partners for use and
! distribution under the JULES collaboration agreement, subject to the terms and
! conditions set out therein.
!
! [Met Office Ref SC0237]
!******************************COPYRIGHT****************************************

MODULE wtrac_checks_jls_iso_mod

!-------------------------------------------------------------------------------
! Description:
!   Check the water isotope prognostics in JULES and output warning or abort
!   model if ratio is over or under specified values.
!
! Code Owner: Please refer to ModuleLeaders.txt
!
! Code Description:
!   Language: Fortran 90.
!-------------------------------------------------------------------------------

USE um_types, ONLY: real_jlslsm
USE jules_water_tracers_mod,  ONLY: n_wtrac_jls
USE wtrac_checks_iso_mod, ONLY: max_r_warn, min_r_warn

#if defined(UM_JULES)
USE timestep_mod, ONLY: timestep_number
#else
USE model_time_mod, ONLY: timestep_number
#endif

IMPLICIT NONE

INTEGER, PARAMETER :: iwt=3     ! CHECK DHO

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='WTRAC_CHECKS_JLS_MOD'

CONTAINS

SUBROUTINE wtrac_checks_sno_iso(land_pts, nsurft, nsmax, surft_pts,            &
                            surft_index, land_index,                           &
                            snow_surft, snow_grnd_surft, sice_surft,           &
                            sliq_surft,                                        &
                            snow_surft_wtrac, snow_grnd_surft_wtrac,           &
                            sice_surft_wtrac, sliq_surft_wtrac)

! Test snow prognostics

USE parkind1,    ONLY: jprb, jpim
USE yomhook,     ONLY: lhook, dr_hook


IMPLICIT NONE

INTEGER, INTENT(IN) :: land_pts
INTEGER, INTENT(IN) :: nsurft
INTEGER, INTENT(IN) :: nsmax

INTEGER, INTENT(IN) :: surft_pts(nsurft)
INTEGER, INTENT(IN) :: surft_index(land_pts,nsurft)   ! Index of tile points
INTEGER, INTENT(IN) :: land_index(land_pts)

! Water prognostics
REAL(KIND=real_jlslsm), INTENT(IN) :: snow_surft(land_pts,nsurft)
REAL(KIND=real_jlslsm), INTENT(IN) :: snow_grnd_surft(land_pts,nsurft)
REAL(KIND=real_jlslsm), INTENT(IN) :: sice_surft(land_pts,nsurft,nsmax)
REAL(KIND=real_jlslsm), INTENT(IN) :: sliq_surft(land_pts,nsurft,nsmax)

! Equivalent water tracer prognostics
REAL(KIND=real_jlslsm), INTENT(IN) ::                                          &
              snow_surft_wtrac(land_pts,nsurft,n_wtrac_jls)
REAL(KIND=real_jlslsm), INTENT(IN) ::                                          &
              snow_grnd_surft_wtrac(land_pts,nsurft,n_wtrac_jls)
REAL(KIND=real_jlslsm), INTENT(IN) ::                                          &
              sice_surft_wtrac(land_pts,nsurft,nsmax,n_wtrac_jls)
REAL(KIND=real_jlslsm), INTENT(IN) ::                                          &
              sliq_surft_wtrac(land_pts,nsurft,nsmax,n_wtrac_jls)

CHARACTER(LEN=*), PARAMETER :: RoutineName='WTRAC_CHECKS_SNO'

!Dr Hook variables
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

!End of header

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! Canopy snow
CALL wtrac_checks_tiles_iso(land_pts, nsurft, 1, surft_pts, surft_index,       &
                            land_index,                                        &
                            snow_surft, snow_surft_wtrac, 'snow_surft')

! Ground snow
CALL wtrac_checks_tiles_iso(land_pts, nsurft, 1, surft_pts, surft_index,       &
                            land_index,                                        &
                   snow_grnd_surft, snow_grnd_surft_wtrac, 'snow_grnd_surft')

! Snow layer ice mass
CALL wtrac_checks_tiles_iso(land_pts, nsurft, nsmax, surft_pts, surft_index,   &
                            land_index,                                        &
                             sice_surft, sice_surft_wtrac, 'sice_surft')

! Snow layer liquid mass
CALL wtrac_checks_tiles_iso(land_pts, nsurft, nsmax, surft_pts, surft_index,   &
                            land_index,                                        &
                             sliq_surft, sliq_surft_wtrac, 'sliq_surft')

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE wtrac_checks_sno_iso

! ---------------------------------------------------------------------------

SUBROUTINE wtrac_checks_hyd_iso(land_pts, nsurft, nsoilt, sm_levels, soil_pts, &
                            soil_index, surft_pts, surft_index, land_index,    &
                            canopy_surft, sthu_soilt,                          &
                            sthf_soilt, smcl_soilt, sthzw_soilt,               &
                            canopy_surft_wtrac, sthu_soilt_wtrac,              &
                            sthf_soilt_wtrac, smcl_soilt_wtrac,                &
                            sthzw_soilt_wtrac)


! Test hydrology prognostics
! (Note, this code is assuming that nsoilt=1.)

USE parkind1,    ONLY: jprb, jpim
USE yomhook,     ONLY: lhook, dr_hook


IMPLICIT NONE

INTEGER, INTENT(IN) :: land_pts
INTEGER, INTENT(IN) :: nsurft
INTEGER, INTENT(IN) :: nsoilt
INTEGER, INTENT(IN) :: sm_levels

INTEGER, INTENT(IN) :: soil_pts
INTEGER, INTENT(IN) :: soil_index(land_pts)
INTEGER, INTENT(IN) :: surft_pts(nsurft)
INTEGER, INTENT(IN) :: surft_index(land_pts,nsurft)   ! Index of tile points
INTEGER, INTENT(IN) :: land_index(land_pts)

! Water prognostics
REAL(KIND=real_jlslsm), INTENT(IN) :: canopy_surft(land_pts,nsurft)
REAL(KIND=real_jlslsm), INTENT(IN) :: sthu_soilt(land_pts,nsoilt,sm_levels)
REAL(KIND=real_jlslsm), INTENT(IN) :: sthf_soilt(land_pts,nsoilt,sm_levels)
REAL(KIND=real_jlslsm), INTENT(IN) :: smcl_soilt(land_pts,nsoilt,sm_levels)
REAL(KIND=real_jlslsm), INTENT(IN) :: sthzw_soilt(land_pts,nsoilt)

! Equivalent water tracer prognostics
REAL(KIND=real_jlslsm), INTENT(IN) :: canopy_surft_wtrac                       &
                                      (land_pts,nsurft,n_wtrac_jls)
REAL(KIND=real_jlslsm), INTENT(IN) :: sthu_soilt_wtrac                         &
                                      (land_pts,nsoilt,sm_levels,n_wtrac_jls)
REAL(KIND=real_jlslsm), INTENT(IN) :: sthf_soilt_wtrac                         &
                                      (land_pts,nsoilt,sm_levels,n_wtrac_jls)
REAL(KIND=real_jlslsm), INTENT(IN) :: smcl_soilt_wtrac                         &
                                      (land_pts,nsoilt,sm_levels,n_wtrac_jls)
REAL(KIND=real_jlslsm), INTENT(IN) :: sthzw_soilt_wtrac                        &
                                      (land_pts,nsoilt,n_wtrac_jls)


INTEGER :: soilt_pts(nsoilt)
INTEGER :: soilt_index(land_pts,nsoilt)


CHARACTER(LEN=*), PARAMETER :: RoutineName='WTRAC_CHECKS_HYD'

!Dr Hook variables
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

!End of header

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! Set up 'tiled' soil points index (assuming nsoilt=1) so that the routine
! wtrac_check_tiles can be used for testing soil fields
! (To be able to use nsoilt>1 in this routine, then these fields need to be
! passed in.)

soilt_pts(1)     = soil_pts
soilt_index(:,1) = soil_index(:)

! Canopy
CALL wtrac_checks_tiles_iso(land_pts, nsurft, 1, surft_pts, surft_index,       &
                            land_index,                                        &
                            canopy_surft, canopy_surft_wtrac, 'canopy_surft')

! Unfrozen soil as fraction of saturation
CALL wtrac_checks_tiles_iso(land_pts, nsoilt, sm_levels, soilt_pts,            &
                            soilt_index, land_index,                           &
                            sthu_soilt, sthu_soilt_wtrac, 'sthu_soilt')

! Frozen soil as fraction of saturation
CALL wtrac_checks_tiles_iso(land_pts, nsoilt, sm_levels, soilt_pts,            &
                            soilt_index, land_index,                           &
                            sthf_soilt, sthf_soilt_wtrac, 'sthf_soilt')

! Soil moisture
CALL wtrac_checks_tiles_iso(land_pts, nsoilt, sm_levels, soilt_pts,            &
                            soilt_index, land_index,                           &
                            smcl_soilt, smcl_soilt_wtrac, 'smcl_soilt')

! Soil moisture fraction in deep-zw layer
CALL wtrac_checks_tiles_iso(land_pts, nsoilt, 1, soilt_pts, soilt_index,       &
                            land_index,                                        &
                            sthzw_soilt, sthzw_soilt_wtrac, 'sthzw_soilt')


IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE wtrac_checks_hyd_iso

! ---------------------------------------------------------------------------

SUBROUTINE wtrac_checks_tiles_iso(land_pts, nsurft, nlayer, surft_pts,         &
                              surft_index, land_index, q, q_wtrac, info_txt)

!
! Description:
!   Compare the normal water tracer to the water field for a single source
!   and output the largest error if larger than a specified relative value.
!   Also abort model if error is too large.
!

USE Ereport_mod,             ONLY: ereport
USE jules_print_mgr,         ONLY: jules_message, jules_print, newline
USE theta_field_sizes,          ONLY: t_i_length,t_j_length

USE yomhook,  ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim

IMPLICIT NONE

INTEGER, INTENT(IN) :: land_pts
INTEGER, INTENT(IN) :: nsurft
INTEGER, INTENT(IN) :: nlayer

INTEGER, INTENT(IN) :: surft_pts(nsurft)
INTEGER, INTENT(IN) :: surft_index(land_pts,nsurft)   ! Index of tile points

INTEGER, INTENT(IN) ::                                                         &
 land_index(land_pts)        ! IN LAND_INDEX(I)=J => the Jth
                             !    point in ROW_LENGTH,ROWS is the
                             !    land point.

! Water field:
REAL(KIND=real_jlslsm), INTENT(IN) :: q(land_pts,nsurft,nlayer)

! Water tracer field:
REAL(KIND=real_jlslsm), INTENT(IN) ::                                          &
                           q_wtrac(land_pts,nsurft,nlayer,n_wtrac_jls)

CHARACTER(LEN=*), INTENT(IN) :: info_txt      ! Text for error message


! Local

INTEGER :: i, j, l, n, k, ilyr                      ! Loop counters
INTEGER :: nmax_q, lmax_q, kmax_q, ilyrmax_q      
INTEGER :: nmin_q, lmin_q, kmin_q, ilyrmin_q      

INTEGER :: ErrorStatus        ! Error code

REAL(KIND=real_jlslsm) :: ratio, max_ratio, min_ratio

CHARACTER(LEN=256) :: message

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='WTRAC_CHECKS_TILES'

! End of header
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

min_ratio = 0.85
max_ratio = 0.85
nmax_q = 0
kmax_q = 0
lmax_q = 0
ilyrmax_q = 0
nmin_q = 0
kmin_q = 0
lmin_q = 0
ilyrmin_q = 0


! --------------------------------------------------------------------------
! Find smallest and largest ratio
! --------------------------------------------------------------------------

DO n = 1,nsurft
  DO k = 1,surft_pts(n)
    IF (surft_pts(n) > 0) THEN
      l = surft_index(k,n)

      DO ilyr = 1, nlayer
        IF (q(l,n,ilyr) > 1.e-10) THEN 
           ratio = q_wtrac(l,n,ilyr,iwt)/q(l,n,ilyr)
        ELSE
           ratio = -100.0
        END IF

        IF ( ratio /= -100.0 .AND. ratio > max_ratio) THEN
          max_ratio = ratio
          nmax_q = n
          kmax_q = k
          lmax_q = l
          ilyrmax_q = ilyr
        END IF

        IF ( ratio /= -100.0 .AND. ratio < min_ratio) THEN
          min_ratio = ratio
          nmin_q = n
          kmin_q = k
          lmin_q = l
          ilyrmin_q = ilyr
        END IF
      END DO
    END IF
  END DO
END DO

! --------------------------------------------------------------------------
! Print out large relative errors if greater than specified values
! --------------------------------------------------------------------------

IF (max_ratio > max_r_warn .AND. lmax_q > 0) THEN

 j = (land_index(lmax_q) - 1) / t_i_length + 1
 i = land_index(lmax_q) - (j-1) * t_i_length
   
  WRITE(jules_message,'(a,i0,1x,a,1x,5(i0,1x),3g20.12)')                       &
   'HIGH RATIO (HDO) END OF TS ', timestep_number,info_txt, nmax_q,lmax_q,     &
    i,j,ilyrmax_q,                                                             &
  q(lmax_q,nmax_q,ilyrmax_q),q_wtrac(lmax_q,nmax_q,ilyrmax_q,iwt), max_ratio
  CALL jules_print(RoutineName, jules_message)

END IF

IF (min_ratio < min_r_warn .AND. lmin_q > 0) THEN

 j = (land_index(lmin_q) - 1) / t_i_length + 1
 i = land_index(lmin_q) - (j-1) * t_i_length

  WRITE(jules_message,'(a,i0,1x,a,1x,5(i0,1x),3g20.12)')                       &
  'LOW RATIO (HDO) END OF TS ', timestep_number,info_txt, nmin_q,lmin_q,       &
  i,j,ilyrmin_q,                                                               &
  q(lmin_q,nmin_q,ilyrmin_q),q_wtrac(lmin_q,nmin_q,ilyrmin_q,iwt), min_ratio
  CALL jules_print(RoutineName, jules_message)

END IF


IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

RETURN
END SUBROUTINE wtrac_checks_tiles_iso

! ---------------------------------------------------------------------------

SUBROUTINE wtrac_checks_grid_iso(row_length, rows, q, q_wtrac, info_txt)
!
! Description:
!   Compare the normal water tracer to the water field for river storage
!   and output the largest error if larger than a specified relative value.
!   Also abort model if error is too large.
!

USE Ereport_mod,             ONLY: ereport
USE jules_print_mgr,         ONLY: jules_message, jules_print, newline

USE yomhook,  ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim

IMPLICIT NONE

INTEGER, INTENT(IN) :: row_length
INTEGER, INTENT(IN) :: rows

! Water field:
REAL(KIND=real_jlslsm), INTENT(IN) :: q(row_length,rows)

! Water tracer field:
REAL(KIND=real_jlslsm), INTENT(IN) :: q_wtrac(row_length,rows,n_wtrac_jls)

CHARACTER(LEN=*), INTENT(IN) :: info_txt      ! Text for error message


! Local

INTEGER :: i, j                               ! Loop counters
INTEGER :: imax_q, jmax_q
INTEGER :: imin_q, jmin_q

INTEGER :: ErrorStatus        ! Error code

REAL(KIND=real_jlslsm) :: ratio, max_ratio, min_ratio

CHARACTER(LEN=256) :: message

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='WTRAC_CHECKS_GRID'

! End of header
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

min_ratio = 0.85
max_ratio = 0.85
imin_q = 0
jmin_q = 0
imax_q = 0
jmax_q = 0


! --------------------------------------------------------------------------
! Find largest error
! --------------------------------------------------------------------------

DO j = 1, rows
  DO i = 1, row_length
    IF ( q(i,j) > 1.e-16) THEN
      ratio = q_wtrac(i,j,iwt) / q(i,j)
    ELSE
      ratio = -100.0
    END IF

    IF ( ratio /= -100.0 .AND. ratio > max_ratio) THEN
      max_ratio = ratio
       imax_q = i
       jmax_q = j
    END IF

    IF ( ratio /= -100.0 .AND. ratio < min_ratio) THEN
      min_ratio = ratio
      imin_q = i
      jmin_q = j
    END IF
  END DO
END DO


IF (max_ratio > max_r_warn .AND. imax_q > 0) THEN

  WRITE(jules_message,'(a,i0,1x,a,1x,2(i0,1x),3g20.12)')                       &
  'HIGH RATIO (HDO) END OF TS ', timestep_number,info_txt, imax_q,jmax_q,   &
  q(imax_q,jmax_q),q_wtrac(imax_q,jmax_q,iwt), max_ratio
  CALL jules_print(RoutineName, jules_message)

END IF

IF (min_ratio < min_r_warn  .AND. imin_q > 0) THEN

  WRITE(jules_message,'(a,i0,1x,a,1x,2(i0,1x),3g20.12)')                       &
  'LOW RATIO (HDO) END OF TS ', timestep_number,info_txt, imin_q,jmin_q,    &
  q(imin_q,jmin_q),q_wtrac(imin_q,jmin_q,iwt), min_ratio
  CALL jules_print(RoutineName, jules_message)

END IF

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

RETURN
END SUBROUTINE wtrac_checks_grid_iso

END MODULE wtrac_checks_jls_iso_mod




