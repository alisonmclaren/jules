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

MODULE water_isotopes_mod

!-------------------------------------------------------------------------------
! Description:
!   This is a dummy module for water isotopes in standalone JULES.
!   Note, water tracers/istopes do not work in standalone JULES.  
!
! Code Owner: Please refer to ModuleLeaders.txt
! This file belongs in Control
!
! Code Description:
!   Language: Fortran 90.
!-------------------------------------------------------------------------------

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='WATER_ISOTOPES_MOD'

CONTAINS

!-------------------------------------------------------------------------

FUNCTION wtrac_is_there_frac(process_txt)
!
! Is fractionation being modelled for this process?
!
IMPLICIT NONE

CHARACTER(len=8) :: process_txt         ! Phase change process

 Output
LOGICAL :: wtrac_is_there_frac          ! Fractionation?

wtrac_is_there_frac = .FALSE.

RETURN
END FUNCTION wtrac_is_there_frac

END MODULE water_isotopes_mod
#endif
