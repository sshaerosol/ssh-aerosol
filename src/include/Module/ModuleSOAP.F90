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
!     2019: Take into account number of layers of particle - Karine Sartelet 
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
      DOUBLE PRECISION liquid(12)

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
           N_aerosol, N_aerosol_layers, neq, q, aero, gas, &
           lwc_Nsize, ionic_Nsize, chp_Nsize,liquid_Nsize,N_size,isoapdyn,&
           aerosol_species_name, spec_name_len, molecular_weight_aer, &
           accomodation_coefficient,&
           nlayer, with_kelvin_effect, tequilibrium, dtaeromin, dorg,&
           coupled_phases, activity_model, epser_soap)

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


    SUBROUTINE SOAP_DYN(rh, &
         ionic, proton, lwc,lwcorg, &
         temp, deltat, &
         DSD, neq, liquid)


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
!     2019: Take into account number of layers of particle - Karine Sartelet 
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
      double precision q_soap(N_size*(1+N_aerosol_layers)+N_aerosol)
      INTEGER IQ(N_aerosol_layers, N_size)
      DOUBLE PRECISION lwc, lwcorg, rh, temp
      DOUBLE PRECISION ionic, proton, chp
      DOUBLE PRECISION liquid(12)

      INTEGER jj,jesp,js,s
      double precision qaero(N_aerosol), qgas(N_aerosol)
      double precision deltat
      double precision DSD(N_size)
      double precision csol(N_size) !! Concentration of solid particles (PBC + PMD)

      integer icpt, i

      neq = N_size * (1 + N_aerosol_layers) + N_aerosol

      !   Pointer for aerosol concentrations
      icpt=N_size
      DO s = 1,N_aerosol_layers
         DO js=1,N_size
            icpt=icpt+1
            IQ(s,js)=icpt
         END DO
      END DO
      
      do js = 1, N_size
         q_soap(js) = concentration_number(js)
      enddo
      
      i = 0 
      do jesp = 1,N_aerosol_layers
         do js = 1, N_size
            i = i + 1
            q_soap(N_size + i) = concentration_mass(js, jesp)
         enddo
      enddo
    
      do jesp = 1,N_aerosol
         q_soap(N_size*(1+N_aerosol_layers) + jesp) = concentration_gas(jesp)
      enddo

      lwcorg = 0.D0

!     ******compute total aerosol mass
      DO s=1,N_aerosol
         qaero(s)=0.D0
      ENDDO
      DO s=1,N_aerosol_layers
        jesp = List_species(s)
        if(aerosol_species_name(jesp).NE.'PH2O') then
!         DO js=(ICUT2+1),nbin_aer
          DO js=1,N_size
             jj=IQ(s,js)
             qaero(jesp)=qaero(jesp)+q_soap(jj)
          END DO
          qgas(jesp) = q_soap(N_size * (1 + N_aerosol_layers) + jesp)
        endif
      enddo
      qgas(N_aerosol)=0.

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
         csol(js) = q_soap(IQ(EMD,js)) + q_soap(IQ(EBC,js))
      enddo

      lwcorg=0.


      CALL soap_main(lwc, rh, temp, ionic, chp, lwcorg,&
           deltat,DSD,csol,liquid,&
           N_aerosol, N_aerosol_layers, neq, q_soap, qaero, qgas, &
           lwc_Nsize, ionic_Nsize, chp_Nsize, liquid_Nsize, N_size, isoapdyn,&
           aerosol_species_name, spec_name_len, molecular_weight_aer, accomodation_coefficient,&
           nlayer, with_kelvin_effect, tequilibrium, dtaeromin, dorg,&
           coupled_phases, activity_model, epser_soap)

      ! Get the calculated values from SOAP
      do js = 1, N_size
         concentration_number(js) = q_soap(js) 
      enddo
      i = 0 
      do jesp = 1,N_aerosol_layers
         do js = 1, N_size
            i = i + 1
            concentration_mass(js, jesp) = q_soap(N_size + i) 
         enddo
      enddo
    
      do jesp = 1,N_aerosol
         concentration_gas(jesp) = q_soap(N_size*(1+N_aerosol_layers) + jesp) 
      enddo

   END SUBROUTINE SOAP_DYN
   
end Module kSOAP
