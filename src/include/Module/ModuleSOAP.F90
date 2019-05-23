!-----------------------------------------------------------------------
!!     Copyright (C) 2012-2018, ENPC - EDF R&D - INERIS
!!     Author(s): Youngseob Kim
!!
!!     This file is part of the Size Composition Resolved Aerosol Model (SCRAM), a
!!     component of the SSH-aerosol model.
!!
!!     SSH-aerosol is a free software; you can redistribute it and/or modify
!!     it under the terms of the GNU General Public License as published
!!     by the Free Software Foundation; either version 2 of the License,
!!     or (at your option) any later version.
!!
!!     SSH-aerosol is distributed in the hope that it will be useful, but
!!     WITHOUT ANY WARRANTY; without even the implied warranty of
!!     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
!!     General Public License for more details.
!!
!!-----------------------------------------------------------------------
Module kSOAP

  use aInitialization
  implicit none

contains

  SUBROUTINE SOAP_EQ(lwcorg, lwc, rh, ionic, proton, &
       temp, aero, gas, liquid)

!------------------------------------------------------------------------
!
!     -- DESCRIPTION
!
!     This subroutine computes SOA concentration using 
!     the equilibrium approach.
!
!------------------------------------------------------------------------
!
!     -- INPUT VARIABLES
!
!     PROTON: hydronium ion concentration ([\mu g.m^-3]).
!     LWC: total liquid water content ([\mu g.m^-3]).
!     RH: relative humidity 0< <1 ([]).
!     TEMP: temperature ([Kelvin]).
!
!     -- INPUT/OUTPUT VARIABLES
!
!     AERO: aerosol bulk concentration ([\mu g.m^-3]).
!     GAS: gas concentration ([\mu g.m^-3]).
!
!     -- OUTPUT VARIABLES
!
!     ORGANION: organic ions ([\mu mol.m^-3]).
!     LWCORG: organic liquid water content ([\mu g.m^-3]).
!     CHP: hydronium ion concentration in water ([mol.L^-1]).
!
!------------------------------------------------------------------------
!
!     -- REMARKS
!
!------------------------------------------------------------------------
!
!     -- MODIFICATIONS
!
!
!------------------------------------------------------------------------
!
!     -- AUTHOR(S)
!
!     2018: Youngseob Kim, CEREA.
!
!------------------------------------------------------------------------

      IMPLICIT NONE

      INTEGER neq
      double precision q(N_size*(1+N_aerosol)+N_aerosol)
      DOUBLE PRECISION aero(N_aerosol),gas(N_aerosol)
      DOUBLE PRECISION lwc, lwcorg, rh, temp
      DOUBLE PRECISION ionic, proton, chp
!      DOUBLE PRECISION ionic_Nsize(N_size)
!      DOUBLE PRECISION chp_Nsize(N_size),lwc_Nsize(N_size)
      DOUBLE PRECISION chp_Nsize(N_size)
!      DOUBLE PRECISION liquid(12),liquid_Nsize(12,N_size)
      DOUBLE PRECISION liquid(12) !,liquid_Nsize(12,N_size)

      INTEGER i,j

      double precision DSD(N_size),csol(N_size), dt2

      neq = N_size * (1 + N_aerosol) + N_aerosol

      csol = 0.D0
      DSD = 0.D0    
      dt2 = 0.0

!    Calculate the concentration of hydronium ion in water
!    microg/m3(=micromol/m3) / microg/m3 (H+ molar mass: 1 g/mol) 
!     = micromol/microg * 1000 
!     = mol/kg = mol/L (Water density: 1 kg/L) 
      chp = proton / lwc * 1.0e3

      CALL soap_main(lwc, rh, temp, ionic, chp, lwcorg, &
           DT2, DSD, csol, liquid,&
           N_aerosol, neq, q, aero, gas, &
           lwc_Nsize, ionic_Nsize, chp_Nsize,liquid_Nsize,N_size)

!     In case there is no gas-phase species.
!     For instance, CB05 mechanism doesn't have GLY for PGLY.
!     If gaseoues species don't exist, gas(j) can't be a gas-phase
!     concentration of the species and it must be set to zero.
      DO i = 1,nesp_aec
         j = aec_species(i)
         IF (aerosol_species_interact(j).LT.0) THEN
            aero(j) = aero(j) + gas(j)
            gas(j) = 0.0
         ENDIF
      ENDDO

    END SUBROUTINE SOAP_EQ


    ! SUBROUTINE SOAP_DYN(rh, &
    !      ionic, proton, lwc,lwcorg, &
    !      temp, deltat, &
    !      DSD, neq, q, iq, liquid, &
    !      lwc_Nsize,ionic_Nsize,proton_Nsize,liquid_Nsize)

    ! SUBROUTINE SOAP_DYN(rh, &
    !      ionic, proton, lwc,lwcorg, &
    !      temp, deltat, &
    !      DSD, neq, q, iq, liquid, &
    !      ionic_Nsize,proton_Nsize,liquid_Nsize)

    SUBROUTINE SOAP_DYN(rh, &
         ionic, proton, lwc,lwcorg, &
         temp, deltat, &
         DSD, neq, q, iq, liquid)


!------------------------------------------------------------------------
!
!     -- DESCRIPTION
!
!     This subroutine computes SOA concentration using 
!     the dynamic approach.
!
!------------------------------------------------------------------------
!
!     -- INPUT VARIABLES
!
!     PROTON: hydronium ion concentration ([\mu g.m^-3]).
!     LWC: total liquid water content ([\mu g.m^-3]).7
!     RH: relative humidity 0< <1 ([]).
!     TEMP: temperature ([Kelvin]).
!
!     -- INPUT/OUTPUT VARIABLES
!
!     AERO: aerosol bulk concentration ([\mu g.m^-3]).
!     GAS: gas concentration ([\mu g.m^-3]).
!
!     -- OUTPUT VARIABLES
!
!     LWCORG: organic liquid water content ([\mu g.m^-3]).
!     CHP: hydronium ion concentration in water ([mol.L^-1]).
!
!------------------------------------------------------------------------
!
!     -- REMARKS
!
!------------------------------------------------------------------------
!
!     -- MODIFICATIONS
!
!
!------------------------------------------------------------------------
!
!     -- AUTHOR(S)
!
!     2018: Youngseob Kim, CEREA.
!
!------------------------------------------------------------------------

      IMPLICIT NONE
    
      INTEGER neq
      double precision q(N_size*(1+N_aerosol)+N_aerosol)
      INTEGER iq(N_aerosol, N_size)
      DOUBLE PRECISION lwc, lwcorg, rh, temp
      DOUBLE PRECISION ionic, proton, chp
      DOUBLE PRECISION liquid(12)

!      DOUBLE PRECISION lwc_Nsize(N_size),proton_Nsize(N_size)
      DOUBLE PRECISION proton_Nsize(N_size)
      DOUBLE PRECISION ionic_Nsize(N_size),chp_Nsize(N_size)
      DOUBLE PRECISION liquid_Nsize(12,N_size)

      INTEGER jj,jesp,js
      double precision qaero(N_aerosol), qgas(N_aerosol)
      double precision deltat
      double precision DSD(N_size)
      double precision csol(N_size) !! Concentration of solid particles (PBC + PMD)

      integer icpt, i
      integer :: ig(N_aerosol)

      neq = N_size * (1 + N_aerosol) + N_aerosol

      !   Pointer for aerosol concentrations
      icpt=N_size
      DO jesp=1,N_aerosol
         DO js=1,N_size
            icpt=icpt+1
            IQ(jesp,js)=icpt
         END DO
      END DO
      
      !   Pointer for gas-phase concentration
      DO jesp=1,N_aerosol
         icpt=icpt+1
         IG(jesp)=icpt
      END DO

      do js = 1, N_size
         q(js) = concentration_number(js)
      enddo

      i = 0 
      do jesp = 1,N_aerosol
         do js = 1, N_size
            i = i + 1
            q(N_size + i) = concentration_mass(js, jesp)
         enddo
      enddo
    
      do jesp = 1,N_aerosol
         q(N_size*(1+N_aerosol) + jesp) = concentration_gas(jesp)
      enddo

      lwcorg = 0.D0

!     ******compute total aerosol mass
      DO jesp=1,N_aerosol-1
         qaero(jesp)=0.D0
!         DO js=(ICUT2+1),nbin_aer
         DO js=1,N_size
            jj=IQ(jesp,js)
            qaero(jesp)=qaero(jesp)+q(jj)
         END DO
         qgas(jesp) = q(N_size * (1 + N_aerosol) + jesp)
      enddo
      qgas(N_aerosol)=0.
      qaero(N_aerosol)=0.

      if (lwc.GT.1.d-19) then
        chp = proton / lwc * 1.0e3
      else
        chp = 1.d-19
        lwc = 1.d-19
      endif

      do js = 1, N_size
        if (lwc_Nsize(js).GT.1.d-19) then
          chp_Nsize(js) = proton_Nsize(js) / lwc_Nsize(js) * 1.0e3
        else
          chp_Nsize(js) = 1.d-19
          lwc_Nsize(js) = 1.d-19
        endif
      enddo

      do js=1,N_size
         csol(js) = q(IQ(EMD,js)) + q(IQ(EBC,js))
      enddo

      !write(*,*) 'SOAPDYN', deltat,N_aerosol,N_size
      lwcorg=0.
      CALL soap_main(lwc, rh, temp, ionic, chp, lwcorg,&
           deltat,DSD,csol,liquid,&
           N_aerosol, neq, q, qaero, qgas, &
           lwc_Nsize, ionic_Nsize, chp_Nsize, liquid_Nsize, N_size)

      ! Get the calculated values from SOAP
      do js = 1, N_size
         concentration_number(js) = q(js) 
      enddo
      i = 0 
      do jesp = 1,N_aerosol
         do js = 1, N_size
            i = i + 1
            concentration_mass(js, jesp) = q(N_size + i) 
         enddo
      enddo
    
      do jesp = 1,N_aerosol
         concentration_gas(jesp) = q(N_size*(1+N_aerosol) + jesp) 
      enddo


   END SUBROUTINE SOAP_DYN
   
end Module kSOAP
