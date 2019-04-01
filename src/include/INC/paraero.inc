C-----------------------------------------------------------------------
C     Copyright (C) 2003-2007, ENPC - INRIA - EDF R&D
C     Author(s): Edouard Debry
C
C     This file is part of the Size Resolved Aerosol Model (SIREAM),
C     which is a component of the air quality modeling system 
C     Polyphemus.
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
C     ALFORG    organics limiting rate (sec-1)   
C     
C     TINYM     threshold mass conc (µg.m-3)     
C     TINYN     threshold num conc (#part.m-3)   
C     
C     DTAEROMIN minimal timestep for aerosols (s)
C     NITEBI number internal iterations EBI      
C
C     NITER_AEC_AQ number of iterations between
C                  isorropia and aec A module
C     NITER_AEC_DRY number of iterations between
C                   aec B module and pankow module
C     NITER_PKW     number of iterations
C                   in pankow module
C	  
C     NITER_POA     number of iterations
C                   in poa module
C     diam_pass		cut diameter = 100nm
C	                diam_pass is define in parameuler.inc 
C
C	  NITAERO		number of subcycling time step
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      INTEGER MAXIT,NITEBI

	DOUBLE PRECISION DMIN
      DOUBLE PRECISION MTSBL,
     &     ALFHP,ALFORG,EPSER,EPSJC,
     &     TINYM,TINYN,EPSEQ,TOLEQ,
     &     DTMAX,TINYA,TINYV,DTAEROMIN
      INTEGER NITER_PKW,NITER_AEC_AQ,NITER_AEC_DRY
	  INTEGER NITER_POA
	  double precision diam_pass
	  integer NITAERO

	PARAMETER(DMIN = 0.01D0)
      PARAMETER(EPSER = 0.01D0)                  
      PARAMETER(EPSJC = 1.D-10)
      PARAMETER(MTSBL = 1.0D0)
      PARAMETER(ALFHP = 0.1D0)
      PARAMETER(TINYM = 0.D0) ! 1.D-19) !!1.D-15)
      PARAMETER(TINYN = 1.D0) 
      PARAMETER(ALFORG = 0.4D0)
      PARAMETER(EPSEQ =0.D0)
      PARAMETER(TOLEQ =1.D-02)
      PARAMETER(MAXIT =100)
      PARAMETER(DTMAX =10.D0)
      PARAMETER(DTAEROMIN =10.D0) !30.D0)
      PARAMETER(NITEBI=5)
      PARAMETER(TINYA =1.D-20)
      PARAMETER(TINYV =1.D-30)
      PARAMETER(NITER_AEC_AQ = 1)
      PARAMETER(NITER_AEC_DRY = 1)
      PARAMETER(NITER_PKW = 5)
	  PARAMETER(NITER_POA = 10)

	  parameter(diam_pass = 1.D-01)
	  parameter(NITAERO = 10)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
