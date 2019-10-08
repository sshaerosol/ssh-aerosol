!!-----------------------------------------------------------------------
!!     Copyright (C) 2019 CEREA (ENPC) - INERIS
!!     SSH-aerosol is distributed under the GNU General Public License v3
!!-----------------------------------------------------------------------

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
