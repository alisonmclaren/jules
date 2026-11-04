#if !defined(UM_JULES)
!******************************COPYRIGHT****************************************
! (c) British Antarctic Survey, 2026.
! All rights reserved.
!
! This routine has been licensed to the other JULES partners for use and
! distribution under the JULES collaboration agreement, subject to the terms and
! conditions set out therein.
!
! [Met Office Ref SC0237]
!******************************COPYRIGHT****************************************

MODULE water_tracers_set_mod

!-------------------------------------------------------------------------------
! Description:
!   This is a dummy module for water tracers in standalone JULES.
!   Note, water tracers/istopes do not work in standalone JULES.  
!
! Code Owner: Please refer to ModuleLeaders.txt
! This file belongs in Control
!
! Code Description:
!   Language: Fortran 90.
!-------------------------------------------------------------------------------

IMPLICIT NONE

! Central control on whether to ensure there are no negative water tracer
! ratios calculated in wtrac_calc_ratio
LOGICAL, PARAMETER :: l_wtrac_no_neg_ratio = .TRUE.

! Min value to determine whether there should be any water tracer
REAL(KIND=real_umphys), PARAMETER :: min_q_tiny_wt = 1.0e-20

! Min value to determine whether to do water tracer/isotope calculations
!REAL(KIND=real_umphys), PARAMETER :: min_q_wt_calc = 1.0e-18
REAL(KIND=real_umphys), PARAMETER :: min_q_wt_calc = 1.0e-12

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='WATER_TRACERS_SET_MOD'

CONTAINS


END MODULE water_tracers_set_mod
#endif
