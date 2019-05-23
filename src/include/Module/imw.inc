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
!*       FOR INTERNAL AND IONIC SPECIES             *
!*       *
!*       IMW(*)   molar weight of inorganic species *
!*       in aqueous_phase                           *
!*       *
!**************************************************
	DOUBLE PRECISION, dimension(IH:IOH) :: IMW = (/&
!******  
!	DATA IMW /  1.0D06,	!H+
     &	 	    1.0D06,&	!H+
     &             23.0D06,&	!Na+
     &             18.0D06,&	!NH4+
     &             35.5D06,&	!Cl-
     &             96.0D06,&	!SO4 2-
     &             97.0D06,&	!HSO4-
     &             63.0D06,&	!NO3-
     &             18.0D06,&	!H2O
     &             17.0D06,&	!NH3
     &             36.5D06,&	!HCl
     &             63.0D06,&	!HNO3
     &             17.0D06  /)	!OH-
!**************************************************
