C-----------------------------------------------------------------------
C     Copyright (C) 2019 CEREA (ENPC) - INERIS
C     SSH-aerosol is distributed under the GNU General Public License v3
C-----------------------------------------------------------------------

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     
C     EPSJC     jacobian epsilon                 
C     EPSER     relative error_precision         
C     EPSEQ     relative equilibrium_precision   
C     TOLEQ     relative equilibrium_tolerance   
C     MAXIT     maximum number of iterations     
C     
C     MTSBL     metastble option of isorropia     
C     DMIN      minimum diameter_size (µm)        
C     ALFHP     H+ limiting rate (sec-1)         
C     
C     TINYM     threshold mass conc (µg.m-3)     
C     TINYN     threshold num conc (#part.m-3)   
C     
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      INTEGER MAXIT,NITEBI

	DOUBLE PRECISION DMIN
      DOUBLE PRECISION MTSBL,
     &     ALFHP,EPSER,EPSJC,
     &     TINYM,TINYN,EPSEQ,TOLEQ

	PARAMETER(DMIN = 0.001D0)
      PARAMETER(EPSJC = 1.D-8)
      PARAMETER(MTSBL = 1.0D0)
      PARAMETER(ALFHP = 0.1D0)
      PARAMETER(TINYM = 1.D-18)
      PARAMETER(TINYN = 1.D-10) 
      PARAMETER(EPSEQ = 1.D-02)
      PARAMETER(TOLEQ =2.D-03)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
