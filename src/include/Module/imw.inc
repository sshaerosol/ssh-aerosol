!!-----------------------------------------------------------------------
!!     Copyright (C) 2019 CEREA (ENPC) - INERIS
!!     SSH-aerosol is distributed under the GNU General Public License v3
!!-----------------------------------------------------------------------


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
