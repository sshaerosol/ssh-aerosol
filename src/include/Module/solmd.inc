!!-----------------------------------------------------------------------
!!     Copyright (C) 2019 CEREA (ENPC) - INERIS
!!     SSH-aerosol is distributed under the GNU General Public License v3
!!-----------------------------------------------------------------------

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
