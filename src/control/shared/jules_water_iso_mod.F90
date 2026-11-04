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

MODULE jules_water_iso_mod

!-------------------------------------------------------------------------------
! Description:
!   This module contains water isotope parameters and common functions.
!
! Code Owner: Please refer to ModuleLeaders.txt
! This file belongs in Control
!
! Code Description:
!   Language: Fortran 90.
!-------------------------------------------------------------------------------

USE um_types,                ONLY: real_jlslsm, real_32
USE water_tracers_mod,       ONLY : wtrac_info
USE free_tracers_inputs_mod, ONLY : hdo_class_wtrac, h218o_class_wtrac
USE jules_print_mgr,         ONLY: jules_message, jules_print


IMPLICIT NONE

! Diffusivity ratio for H218O and HDO relative to H216O from Merlivat (1978) 
REAL(KIND=real_32), PARAMETER :: d_ratio_h218o = 0.9727
REAL(KIND=real_32), PARAMETER :: d_ratio_hdo   = 0.9757

! Which kinetic fractionation parameterisation to use for sea surface evap?
CHARACTER(LEN=5), PARAMETER :: wtrac_kin_factor_method = 'const'
!CHARACTER(LEN=5), PARAMETER :: wtrac_kin_factor_method = 'K&J79'
                                                   ! 'K&J79' or 'const'

! Constant kinetic fractionation values following Cauquoin and Werner (2021)
! These are mean values inbetween rough and smooth regimes (Benetti et al.
! 2014) where their alpha_k= 1.00475 and 1.00418. Converted using
! kin_factor = 1/alpha_k
REAL(KIND=real_jlslsm), PARAMETER :: kin_factor_h218o = 0.9953
REAL(KIND=real_jlslsm), PARAMETER :: kin_factor_hdo   = 0.9958

CONTAINS

!-------------------------------------------------------------------------

FUNCTION wtrac_eqm_factor_liqvap(i_wt,t,process_txt)

!
! Calculate liquid/vapor equilibrium fractionation factors based on
! Horita and Mesolowski (1994).
!
! Reference: Horita, J. and Wesolowski, D. J. (1994), Liquid-vapor 
! fractionation of oxygen and hydrogen isotopes of water from the 
! freezing to the critical temperature. Geochim. Cosmochim. Acta, 58(16), 
! 3425–3437.
!
! Calculated using real_32 option only. Doesn't warrant creating a second
! routine at real_64.

USE water_isotopes_mod, ONLY : wtrac_is_there_frac

IMPLICIT NONE

INTEGER, INTENT(IN)  :: i_wt                      ! Water tracer number

REAL(KIND=real_32), INTENT(IN) :: t               ! Temperature (in K)

CHARACTER(len=8) :: process_txt                   ! Phase process being
                                                  ! modelled

! Output
REAL(KIND=real_32) :: wtrac_eqm_factor_liqvap     ! Equilibrium fractionation
                                                  ! factor
! Local

LOGICAL :: l_iso_frac                             ! Fractionation?


! Are we modelling fractionation for this process?
! (This is only set to FALSE for testing purposes)
l_iso_frac = wtrac_is_there_frac(process_txt)

IF (l_iso_frac) THEN
  ! Set fractionation factors for isotopes
  SELECT CASE (wtrac_info(i_wt)%wt_class)

    CASE(h218o_class_wtrac)

      wtrac_eqm_factor_liqvap = exp((-7.685e-3)                                &
                                + (6.7123 / t)                                 &
                                - (1.6664e+03 / (t**2))                        &
                                + (0.35041e+06 / (t**3)))

    CASE(hdo_class_wtrac)

      wtrac_eqm_factor_liqvap = exp((1158.8e-12 * (t**3))                      &
                                - (1620.1e-09 * (t**2))                        &
                                + (794.84e-06 * t)                             &
                                - (161.04e-03)                                 &
                                + (2.9992e+06 / (t**3)))

    CASE DEFAULT

      ! For non-isotopic water tracers - no fractionation
      wtrac_eqm_factor_liqvap = 1.0 

  END SELECT
ELSE
  ! No fractionation for all water tracers
  wtrac_eqm_factor_liqvap = 1.0 
END IF

RETURN
   
END FUNCTION wtrac_eqm_factor_liqvap

!---------------------------------------------------------------------------

FUNCTION wtrac_eqm_factor_icevap(i_wt,t,process_txt)

!
! Calculate ice/vapor equilibrium fractionation factors based on
! Merlivat and Nief (1967) and Majoube (1971). 
!
! Reference: 
! For H218O:  Majoube (1971). Fractionnement en oxygene 18 entre la 
!             glace et la vapeur d'eau, J. Chem. Phys., 68, 625-636
! For HDO:    Merlivat, L. & Nief, G.(1967). Isotopic fractionation
!             during change of state solid-vapour and liquid-vapour of
!             water at temperatures below 08C, Tellus, 19, 122–126. 
!
! Calculated using real_32 option only. Doesn't warrant creating a second
! routine at real_64.

USE water_isotopes_mod, ONLY : wtrac_is_there_frac

IMPLICIT NONE

INTEGER, INTENT(IN)  :: i_wt                      ! Water tracer number

REAL(KIND=real_32), INTENT(IN) :: t               ! Temperature (in K)

CHARACTER(len=8), INTENT(IN)  :: process_txt      ! Phase process being
                                                  ! modelled
! Output
REAL(KIND=real_32) :: wtrac_eqm_factor_icevap     ! Equilibrium fractionation
                                                  ! factor

LOGICAL :: l_iso_frac                             ! Fractionation?


! Are we modelling fractionation for this process?
! (This is only set to FALSE for testing purposes)
l_iso_frac = wtrac_is_there_frac(process_txt)

IF (l_iso_frac) THEN
  ! Set fractionation factors for isotopes

  SELECT CASE (wtrac_info(i_wt)%wt_class)

    CASE(h218o_class_wtrac)

      wtrac_eqm_factor_icevap = exp(-28.224e-3 + (11.839 / t))

    CASE(hdo_class_wtrac)

      wtrac_eqm_factor_icevap = exp( -9.45e-2 + (16289.0 / (t**2.0)))

    CASE DEFAULT
      ! For non-isotopic water tracers - no fractionation
      wtrac_eqm_factor_icevap = 1.0

  END SELECT
ELSE
  ! No fractionation for all water tracers
  wtrac_eqm_factor_icevap = 1.0 
END IF

RETURN

END FUNCTION wtrac_eqm_factor_icevap

!-------------------------------------------------------------------------

FUNCTION wtrac_kin_factor_kj(i_wt,u_star,z1)

!
! Calculate a kinetic fractionation factor for surface evaporation
! using Merlivat and Jouzel (1979)
!
! Note, this routine is using u_star as input and calculating z0
! using Charnock's equation.  
  
USE planet_constants_mod, ONLY: g, vkman
USE jules_sea_seaice_mod, ONLY: charnock
USE water_isotopes_mod,   ONLY: l_wtrac_frac_sfc_evap

IMPLICIT NONE

INTEGER, INTENT(IN)  :: i_wt               
                                     ! Water tracer number

REAL(KIND=real_jlslsm), INTENT(IN) :: u_star
                                     ! Friction velocity (m/s)

REAL(KIND=real_jlslsm), INTENT(IN) :: z1
                                     ! Height of lowest atmos level


! Return
REAL(KIND=real_jlslsm) :: wtrac_kin_factor_kj

REAL(KIND=real_jlslsm), PARAMETER :: visc_air = 1.4e-5
                                     ! Kinematic viscosity of air
                                     ! Copied from sea_rough_int routine (m2/s)

REAL(KIND=real_jlslsm), PARAMETER :: d_air = 2.36e-5
                                     ! Molecular diffusivity of air

REAL(KIND=real_jlslsm) :: res_ratio ! turbulent resistance / molecular 
                                    ! resistance

REAL(KIND=real_jlslsm) :: d_ratio   ! Diffusivity ratio = H216O / isotope
REAL(KIND=real_jlslsm) :: n         ! factor
REAL(KIND=real_jlslsm) :: d_ratio_n ! d_ratio ** n
REAL(KIND=real_jlslsm) :: k_mj      ! k from Merlivat and Jouzel
REAL(KIND=real_jlslsm) :: Re        ! Reynolds number
REAL(KIND=real_jlslsm) :: z0        ! Roughness length

IF (l_wtrac_frac_sfc_evap) THEN

  SELECT CASE (wtrac_info(i_wt)%wt_class)

    CASE(h218o_class_wtrac)

      d_ratio = 1.0 / REAL(d_ratio_h218o,KIND=real_jlslsm)

    CASE(hdo_class_wtrac)

      d_ratio = 1.0 / REAL(d_ratio_hdo,KIND=real_jlslsm)

    CASE DEFAULT
      ! Should not be here! For non-isotopic water tracers - no fractionation

  END SELECT

  ! Charnock's equation
  z0 = charnock * u_star**2.0 / g     ! Read in z0hsea or z0h_ice

  ! Calculate Reynolds number
  Re = u_star * z0 / visc_air

  IF (Re < 1.0) THEN
    ! Smooth regime
    n = 2.0/3.0
    res_ratio = (1.0/vkman) * log (u_star* z1 / (30.0*visc_air))               &
                                / (13.6 * (visc_air/d_air)**n)
  
  ELSE
    ! Rough regime
    n = 1.0/2.0
    res_ratio = ( (1.0/vkman) * log(z1/z0) - 5.0 )                             &
                                 / ( 7.3*(Re**0.25) * (visc_air/d_air)**n )
  END IF

  d_ratio_n = d_ratio**n

  ! Get k as in MJ79
  k_mj = ( d_ratio_n - 1.0) / (d_ratio_n + res_ratio) 

  wtrac_kin_factor_kj = 1.0 - k_mj

ELSE
  ! No fractionation
  wtrac_kin_factor_kj = 1.0
END IF

RETURN

END FUNCTION wtrac_kin_factor_kj

!-------------------------------------------------------------------------

FUNCTION wtrac_kin_factor_const(i_wt)

!
! Use constant value for kinetic fractionation factor based on
! mean values imbetween rough and smooth K&J 79 regimes (Benetti et al. 2014)
!

USE water_isotopes_mod,   ONLY: l_wtrac_frac_sfc_evap

IMPLICIT NONE

INTEGER, INTENT(IN)  :: i_wt               
                                     ! Water tracer number


! Return
REAL(KIND=real_jlslsm) :: wtrac_kin_factor_const

IF (l_wtrac_frac_sfc_evap) THEN

  SELECT CASE (wtrac_info(i_wt)%wt_class)

    CASE(h218o_class_wtrac)

      wtrac_kin_factor_const = kin_factor_h218o

    CASE(hdo_class_wtrac)

      wtrac_kin_factor_const= kin_factor_hdo

    CASE DEFAULT
      ! Should not be here! For non-isotopic water tracers - no fractionation

  END SELECT
ELSE
  ! No fractionation
  wtrac_kin_factor_const = 1.0
END IF

END FUNCTION wtrac_kin_factor_const

END MODULE jules_water_iso_mod
