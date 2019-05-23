!C-----------------------------------------------------------------------
!C     Copyright (C) 2003-2007, ENPC - INRIA - EDF R&D
!C     Author(s): Edouard Debry
!C     
!C     This file is part of the Size Resolved Aerosol Model (SIREAM), a
!C     component of the air quality modeling system Polyphemus.
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

!**************************************************
!*       *
!*       LIST OF MOLAR WEIGHT EXPRESSED IN µg.mol-1 *
!*       FOR INTERNAL AND EXTERNAL SPECIES          *
!*       *
!*       SMW(*)   molar weight of solids            *
!*       *
!**************************************************
	DOUBLE PRECISION, dimension(SNaNO3:SLC) :: SMW = (/&
!******  
!	DATA SMW / 
     &      	   85.0D06,&	!NaNO3
     &             80.0D06,&	!NH4NO3
     &             58.5D06,&	!NaCl
     &             53.5D06,&	!NH4Cl
     &             142.0D06,&	!Na2SO4
     &             132.0D06,&	!NH42S4
     &             120.0D06,&	!NaHSO4
     &             115.0D06,&	!NH4HS4
     &             247.0D06  /)	!LC
!**************************************************
