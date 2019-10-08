!C-----------------------------------------------------------------------
!C     Copyright (C) 2019 CEREA (ENPC) - INERIS
!C     SSH-aerosol is distributed under the GNU General Public License v3
!C-----------------------------------------------------------------------

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!C     !C
!C     PI       pi number    
!C     cst_pi6  constant pi/6
!C     cst_FRAC3    infinite fraction 1./3 
!C     
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      DOUBLE PRECISION cst_FRAC3

      PARAMETER (cst_FRAC3=0.333333333333333333D0)

      DOUBLE PRECISION cst_PI6,cst_SQRT2
      PARAMETER(cst_PI6=0.52359877559829887307D0)
      PARAMETER (cst_SQRT2 = 1.41421356237309504880d0)

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!C     RGAS       perfect gas constant    
!C     ATM        atmospheric pressure    
!C     MUair      dynamic viscosity of air
!C     MMair      molar mass of air       
!C     RHOair     specific mass of air   
!C     SIGMAIR    surface tension	    
!C     PARMAIR    collision factor
!C     RHOwater   water density 
!C     MMwater    molar mass of water	  
!C     
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      DOUBLE PRECISION SIGMair,PARMair
      DOUBLE PRECISION RGAS,ATM
      DOUBLE PRECISION MMair,RHOair,MUair
      DOUBLE PRECISION RHOwater, MMwater
      
      PARAMETER (RGAS  = 8.31400D0) ! J.mol-1.K-1
      parameter (ATM      = 1.01325D05) ! Pascals       
      parameter (MUair    = 1.72500D-05) ! kg.m-1.s-1 
      parameter (MMair    = 2.89700D-02) ! kg.mol-1
      parameter (RHOair   = 1.29280D0)  ! kg.m-3
      parameter (SIGMair  = 3.617D0) ! N.m-1
      parameter (PARMair  = 9.700D01)   
      parameter (RHOwater = 1.D03) ! kg.m-3
      parameter (MMwater  = 1.801528d-02) ! kg.mol-1
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

!C     LWCmin = threshold lwc (g/m3) for aq chem
!C      DOUBLE PRECISION LWCmin  
!C      PARAMETER (LWCmin = 0.05d0)

!C     Threshold for Relative Humidity
      DOUBLE PRECISION Threshold_RH_inf, Threshold_RH_sup
      PARAMETER (Threshold_RH_inf=1.d-2)
      parameter (Threshold_RH_sup=97.d-2)

