!C-----------------------------------------------------------------------
!C     Copyright (C) 2019 CEREA (ENPC) - INERIS
!C     SSH-aerosol is distributed under the GNU General Public License v3
!C-----------------------------------------------------------------------

!C------------------------------------------------------------------------
!C     
!C     -- DESCRIPTION 
!C     
!C     This file provides physical constants.
!C     
!C------------------------------------------------------------------------
!C     
!C     -- INPUT VARIABLES
!C     
!C     -- INPUT/OUTPUT VARIABLES
!C     
!C     -- OUTPUT VARIABLES
!C     
!C------------------------------------------------------------------------
!C     
!C     -- REMARKS
!C     
!C------------------------------------------------------------------------
!C     
!C     -- MODIFICATIONS
!C     
!C------------------------------------------------------------------------
!C     
!C     -- AUTHOR(S)
!C     
!C------------------------------------------------------------------------

!C     Pi.
      DOUBLE PRECISION PI
      PARAMETER(PI=3.14159265358979323846D0)

!C     Earth radius ([m]).
      DOUBLE PRECISION Ray
      PARAMETER(RAY=6371229.D0)

!C     Gravity acceleration.
      DOUBLE PRECISION PG
      PARAMETER(PG=9.8D0)

!C     Van Karman constant (diffusion).
      DOUBLE PRECISION PKa
      PARAMETER(PKa=0.4D0)

!C     Constants for diffusion.
      DOUBLE PRECISION Plinf
      PARAMETER(PLinf=80.D0)
      DOUBLE PRECISION Pz0
      PARAMETER(Pz0=0.6D0)

!C     For Richardson number.
      DOUBLE PRECISION PB,PC,PD
      PARAMETER(PB=5.D0,PC=5.D0,PD=5.D0)

!C     Avogadro number.
      DOUBLE PRECISION Navog
      PARAMETER (Navog = 6.02213D+23)

!C     Perfect gaz constant ([atm.L/mol/K]).
      DOUBLE PRECISION Pr
      PARAMETER (Pr = 8.20d-2)

!C     specific heat to constant pressure.
      DOUBLE PRECISION PCp
      PARAMETER (PCp = 1.005d3)

!C     Kelvin correction.
      DOUBLE PRECISION Ptkelvin
      PARAMETER (Ptkelvin = 273.15d0)

!C     Boltzmann constant
      DOUBLE PRECISION Pkboltz
      PARAMETER (Pkboltz = 1.38d-23)

!C     Maximum number of subcycle for advection
      INTEGER Ncycle_adv_max
      PARAMETER(Ncycle_adv_max = 100)

