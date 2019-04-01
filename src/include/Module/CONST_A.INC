!C-----------------------------------------------------------------------
!C     Copyright (!C) 2003-2007, ENPC - INRIA - EDF R&D
!C     Author(s): Kathleen Fahey, Edouard Debry and Bruno Sportisse
!C
!C     This file is part of the Size Resolved Aerosol Model (SIREAM),
!C     which is a component of the air quality modeling system 
!C     Polyphemus.
!C
!C     Polyphemus is developed in the INRIA - ENPC joint project-team
!C     CLIME and in the ENPC - EDF R&D joint laboratory CEREA.
!C    
!C     Polyphemus is free software; you can redistribute it and/or modify
!C     it under the terms of the GNU General Public License as published
!C     by the Free Software Foundation; either version 2 of the License,
!C     or (at your option) any later version.
!C     
!C     Polyphemus is distributed in the hope that it will be useful, but
!C     WITHOUT ANY WARRANTY; without even the implied warranty of
!C     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
!C     General Public License for more details.
!C     
!C     For more information, visit the Polyphemus web site:
!C     http://cerea.enpc.fr/polyphemus/
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
      
      PARAMETER (RGAS  = 8.31400D0,& ! J.mol-1.K-1
      &    ATM      = 1.01325D05,& ! Pascals       
      &    MUair    = 1.72500D-05,& ! kg.m-1.s-1 
      &    MMair    = 2.89700D-02,& ! kg.mol-1
      &    RHOair   = 1.29280D0,& ! kg.m-3
      &    SIGMair  = 3.617D0,& ! N.m-1
      &    PARMair  = 9.700D01,&  
      &    RHOwater = 1.D03,& ! kg.m-3
      &    MMwater  = 1.801528d-02) ! kg.mol-1
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

!C     LWCmin = threshold lwc (g/m3) for aq chem
!C      DOUBLE PRECISION LWCmin  
!C      PARAMETER (LWCmin = 0.05d0)

!C     Threshold for Relative Humidity
      DOUBLE PRECISION Threshold_RH_inf, Threshold_RH_sup
      PARAMETER (Threshold_RH_inf=1.d-2,&
      &    Threshold_RH_sup=97.d-2)

