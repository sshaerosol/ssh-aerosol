C-----------------------------------------------------------------------
C     Copyright (C) 2003-2007, ENPC - INRIA - EDF R&D
C     Author(s): Kathleen Fahey
C
C     This file is part of the Variable Size Resolved Model (VSRM),
C     based on the VSRM model of Carnegie Melon University. It is a
C     component of the air quality modeling system Polyphemus.
C
C     Polyphemus is developed in the INRIA - ENPC joint project-team
C     CLIME and in the ENPC - EDF R&D joint laboratory CEREA.
C    
C     Polyphemus is free software; you can redistribute it and/or modify
C     it under the terms of the GNU General Public License as published
C     by the Free Software Foundation; either version 2 of the License,
C     or (at your option) any later version.
C     
C     Polyphemus is distributed in the hope that it will be useful, but
C     WITHOUT ANY WARRANTY; without even the implied warranty of
C     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
C     General Public License for more details.
C     
C     For more information, visit the Polyphemus web site:
C     http://cerea.enpc.fr/polyphemus/
C-----------------------------------------------------------------------

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     Data for heterogeneous reactions:
C ----------------------------------------------------
C     HO2  -->  0.5 H2O2
C     NO2  -->  0.5 HONO + 0.5 HNO3
C     NO3  -->  HNO3
C     N2O5 -->  2 HNO3
C ----------------------------------------------------
C     Gamma: reaction probabilities                    
C     Ref: Heterogeneous chemistry and tropospheric     
C     ozone, D.J. Jacob, Atm. Env., 2000, 34,          
C     pp 2131-2159                                     
C     range for Gamma1:[0.1;1]       recommended: 0.2        
C     range for Gamma2:[1.d-6;1.d-3] recommended: 1.d-4
C     range for Gamma3:[2.d-4;1.d-2] recommended: 1.d-3
C     range for Gamma4:[0.01;1]      recommended: 0.03      
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      DOUBLE PRECISION Gamma(4)
      DATA Gamma / 0.2d0,
     &     0.0,
     &     1.d-3,
     &     0.03d0 /
      
      DOUBLE PRECISION SIGM_NO2,PARM_NO2
      PARAMETER(SIGM_NO2=3.765d0)
      PARAMETER(PARM_NO2=210.d0)
