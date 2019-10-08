!-----------------------------------------------------------------------
!     Copyright (C) 2019 CEREA (ENPC) - INERIS
!     SSH-aerosol is distributed under the GNU General Public License v3
!-----------------------------------------------------------------------

    double precision ::  ALFHP! percentage of H+ allowed to c/e(0.1)
!    double precision ::  EPSER
    double precision ::  TINYM,TINYN,MTSBL
!    double precision ::  DTAEROMIN
!    parameter(EPSER = 1.D-2) !! Relative error for time step adjustment
    parameter(ALFHP = 0.1D0)
    parameter(TINYM = 1.d-18)
    parameter(TINYN = 1.d-10)
    parameter(MTSBL = 1.0D0)
