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
!*       LIST OF SOLIDS MASS DENSITIES              *
!*       EXPRESSED IN µg.µm-3                       *
!*       *
!**************************************************
	DOUBLE PRECISION, dimension(SNaNO3:SLC) ::  SMD = (/&
!******  
!	DATA SMD / 2.260D-06,	!NaNO3
     &	           2.260D-06,&	!NaNO3
     &             1.725D-06,&	!NH4NO3
     &             2.165D-06,&	!NACl
     &             1.530D-06,&	!NH4Cl
     &             2.700D-06,&	!Na2SO4
     &             1.770D-06,&	!NH42S4
     &             2.740D-06,&	!NaHSO4
     &             1.780D-06,&	!NH4HS4
     &             1.770D-06  /)!LC
!**************************************************
