C-----------------------------------------------------------------------
C     Copyright (C) 2019 CEREA (ENPC) - INERIS
C     SSH-aerosol is distributed under the GNU General Public License v3
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
