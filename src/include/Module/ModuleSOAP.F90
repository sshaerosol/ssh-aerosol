!!-----------------------------------------------------------------------
!!     Copyright (C) 2019 CEREA (ENPC) - INERIS
!!     SSH-aerosol is distributed under the GNU General Public License v3
!!-----------------------------------------------------------------------

Module kSOAP

  use aInitialization
  implicit none

contains

  SUBROUTINE SSH_SOAP_EQ(lwcorg, lwc, rh, ionic, proton, &
       temp, aero, gas, liquid, delta_t, qaq)

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
!     QAQ: aerosol bulk concentration in the aqueous phase
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
      DOUBLE PRECISION delta_t
      DOUBLE PRECISION qaq(N_aerosol)

      INTEGER i,j,isoapdyn2

      double precision DSD(N_size),csol(N_size), dt2

      isoapdyn2=0

      neq = N_size * (1 + N_aerosol) + N_aerosol

      csol = 0.D0
      DSD = 0.D0    
      dt2 = delta_t

!    Calculate the concentration of hydronium ion in water
!    microg/m3(=micromol/m3) / microg/m3 (H+ molar mass: 1 g/mol) 
!     = micromol/microg * 1000 
      !     = mol/kg = mol/L (Water density: 1 kg/L)

      ! proton must not be negative,
      ! 'if' for proton  is added so that
      ! the program is not stopped.
      if (lwc>0.d0 .and. proton > 0.d0) then
         chp = proton / lwc * 1.0e3
      else
         chp= 1.0e-7
      endif
      
      CALL soap_main_ssh(lwc, rh, temp, ionic, chp, lwcorg, &
           DT2, DSD, csol, liquid,&
           N_aerosol, N_aerosol_layers, neq, q, aero, qaq, gas, &
           lwc_Nsize, ionic_Nsize, chp_Nsize,liquid_Nsize,N_size,isoapdyn2, &
           imethod, soap_inorg_loc,&
           aerosol_species_name, spec_name_len, molecular_weight_aer, &
           accomodation_coefficient, aerosol_type, &
           partitioning, smiles, saturation_vapor_pressure, enthalpy_vaporization, diffusion_coef,&
           nlayer, with_kelvin_effect, tequilibrium, dtaeromin, dorg,&
           coupled_phases, activity_model, epser_soap, i_hydrophilic, N_inert, N_inorganic,&
           with_oligomerization)

!     In case there is no gas-phase species.
!     For instance, CB05 mechanism doesn't have GLY for PGLY.
!     If gaseoues species don't exist, gas(j) can't be a gas-phase
      !     concentration of the species and it must be set to zero.
      frac_oligo(:)=0.d0
      DO i = 1,nesp_aec
         j = aec_species(i)        
         IF (aerosol_species_interact(j).LT.0) THEN
            aero(j) = aero(j) + gas(j)
            gas(j) = 0.0
         ENDIF

         IF (oligo_index(j)>0) THEN
             
            if (aero(oligo_index(j))+aero(j)>0.d0) then
               frac_oligo(j)=aero(oligo_index(j))/(aero(oligo_index(j))+aero(j))
               aero(j)=aero(j)+aero(oligo_index(j))
               aero(oligo_index(j))=0.               
            endif
         ENDIF        
      ENDDO
     
      if (soap_inorg==1) then
         lwc=aero(EH2O)
         proton = chp * lwc / 1.0e3
      endif
      
      
    END SUBROUTINE SSH_SOAP_EQ


    SUBROUTINE SSH_SOAP_DYN(rh, &
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
      DOUBLE PRECISION qaq(N_aerosol)

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

      q_soap = 0.0
      qgas = 0.0

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
        if (inon_volatile(jesp).EQ.0..and. aerosol_species_interact(jesp).GT.0) then
           q_soap(N_size*(1+N_aerosol_layers) + jesp) = concentration_gas(jesp)
        endif
      enddo

      lwcorg = 0.D0

!     ******compute total aerosol mass
      DO s=1,N_aerosol
         qaero(s)=0.D0
      ENDDO
      DO s=1,N_aerosol_layers
        jesp = List_species(s)
        if(s.NE.EH2O_layers) then
          DO js=1,N_size
             jj=IQ(s,js)
             qaero(jesp)=qaero(jesp)+q_soap(jj)
          END DO
          if (inon_volatile(jesp).EQ.0.and. aerosol_species_interact(jesp).GT.0) then
             qgas(jesp) = q_soap(N_size * (1 + N_aerosol_layers) + jesp)
          endif
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
         csol(js) = 0.d0
         do jesp=1,N_inert 
            csol(js) = csol(js) + q_soap(IQ(jesp,js))
         enddo
      enddo
      
      lwcorg=0.
      CALL soap_main_ssh(lwc, rh, temp, ionic, chp, lwcorg,&
           deltat,DSD,csol,liquid,&
           N_aerosol, N_aerosol_layers, neq, q_soap, qaero, qaq, qgas, &
           lwc_Nsize, ionic_Nsize, chp_Nsize, liquid_Nsize, N_size, isoapdyn, &
           imethod, soap_inorg_loc, &
           aerosol_species_name, spec_name_len, molecular_weight_aer, &
           accomodation_coefficient, aerosol_type, &
           partitioning, smiles, saturation_vapor_pressure, enthalpy_vaporization, diffusion_coef,&
           nlayer, with_kelvin_effect, tequilibrium, dtaeromin, dorg,&
           coupled_phases, activity_model, epser_soap, i_hydrophilic, N_inert, N_inorganic,&
           with_oligomerization)

      ! Get the calculated values from SOAP
      do js = 1, N_size
         concentration_number(js) = q_soap(js) 
      enddo
      i = 0 
      do jesp = 1,N_aerosol_layers
         s = List_species(jesp)
         do js = 1, N_size
            i = i + 1
            ! Need to redistribute the mass even for non volatile particles because layers may have changed
            if ((aerosol_type(s)==4.and.soap_inorg_loc>=0).or. &
                 ((aerosol_type(s)==3.or.s==EH2O).and.(soap_inorg_loc==1.or.soap_inorg_loc==-1))) then
               concentration_mass(js, jesp) = q_soap(N_size + i)
            endif         
         enddo
      enddo
    
      do jesp = 1,N_aerosol
         if (inon_volatile(jesp).EQ.0) then
            if ((aerosol_type(jesp)==4.and.soap_inorg_loc>=0).or. &
                 (aerosol_type(jesp)==3.and.(soap_inorg_loc==1.or.soap_inorg_loc==-1))) then
               concentration_gas(jesp) = q_soap(N_size*(1+N_aerosol_layers) + jesp)
            endif
         else
            if ( aerosol_species_interact(jesp).GT.0) then
               ! Add initial gas concentration that were not added to soap
               if ((aerosol_type(jesp)==4.and.soap_inorg_loc>=0).or. &
                    (aerosol_type(jesp)==3.and.(soap_inorg_loc==1.or.soap_inorg_loc==-1))) then
                  concentration_gas(jesp) = concentration_gas(jesp) + q_soap(N_size*(1+N_aerosol_layers) + jesp)
               endif
            endif
         endif
      enddo

   END SUBROUTINE SSH_SOAP_DYN

   SUBROUTINE SSH_SOAP_DYN_ICUT(rh, &
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
      double precision q_soap((N_size-ICUT_org)*(1+N_aerosol_layers)+N_aerosol)
      INTEGER IQ(N_aerosol_layers, (N_size-ICUT_org))
      DOUBLE PRECISION lwc, lwcorg, rh, temp
      DOUBLE PRECISION ionic, proton, chp
      DOUBLE PRECISION liquid(12)

      INTEGER jj,jesp,js,s,isoapdyn2
      double precision qaero(N_aerosol), qgas(N_aerosol)
      double precision deltat
      double precision DSD(N_size-ICUT_org)
      double precision csol(N_size-ICUT_org) !! Concentration of solid particles (PBC + PMD)
      DOUBLE PRECISION qaq(N_aerosol)

      integer icpt, i

      isoapdyn2=1

      neq = (N_size-ICUT_org) * (1 + N_aerosol_layers) + N_aerosol

      !   Pointer for aerosol concentrations
      icpt=N_size-ICUT_org
      DO s = 1,N_aerosol_layers
         DO js=ICUT_org+1,N_size
            icpt=icpt+1
            IQ(s,js-ICUT_org)=icpt
         END DO
      END DO

      q_soap = 0.0
      qgas = 0.0

      do js = ICUT_org+1, N_size
         q_soap(js-ICUT_org) = concentration_number(js)
      enddo
      
      i = 0 
      do jesp = 1,N_aerosol_layers
         do js = ICUT_org+1, N_size
            i = i + 1
            q_soap(N_size - ICUT_org + i) = concentration_mass(js, jesp)
         enddo
      enddo
      do jesp = 1,N_aerosol
        if (inon_volatile(jesp).EQ.0) then
           q_soap((N_size-ICUT_org)*(1+N_aerosol_layers) + jesp) = concentration_gas(jesp)
        endif
      enddo

      lwcorg = 0.D0

!     ******compute total aerosol mass
      DO s=1,N_aerosol
         qaero(s)=0.D0
      ENDDO
      DO s=1,N_aerosol_layers
        jesp = List_species(s)
        if(aerosol_species_name(jesp).NE.'PH2O') then
          DO js=ICUT_org+1,N_size
             jj=IQ(s,js-ICUT_org)
             qaero(jesp)=qaero(jesp)+q_soap(jj)
          END DO
          if (inon_volatile(jesp).EQ.0) then
             qgas(jesp) = q_soap((N_size-ICUT_org) * (1 + N_aerosol_layers) + jesp)
          endif
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

      do js=ICUT_org+1,N_size
         csol(js-ICUT_org) = 0.d0 
         do jesp=1,N_inert
            csol(js-ICUT_org) = csol(js-ICUT_org) + q_soap(IQ(jesp,js-ICUT_org))
         enddo
      enddo

      lwcorg=0.

           
      CALL soap_main_ssh(lwc, rh, temp, ionic, chp, lwcorg,&
           deltat,DSD,csol,liquid,&
           N_aerosol, N_aerosol_layers, neq, q_soap, qaero, qaq, qgas, &
           lwc_Nsize(ICUT_org+1:N_size), ionic_Nsize(ICUT_org+1:N_size), chp_Nsize(ICUT_org+1:N_size), &
           liquid_Nsize(:,ICUT_org+1:N_size), N_size-ICUT_org, isoapdyn2, &
           imethod, soap_inorg_loc, &
           aerosol_species_name, spec_name_len, molecular_weight_aer, accomodation_coefficient, &
           aerosol_type, partitioning, smiles, saturation_vapor_pressure, enthalpy_vaporization, diffusion_coef,&
           nlayer, with_kelvin_effect, tequilibrium, dtaeromin, dorg,&
           coupled_phases, activity_model, epser_soap, i_hydrophilic, N_inert, N_inorganic,&
           with_oligomerization)

      ! Get the calculated values from SOAP
      do js = ICUT_org+1, N_size
         concentration_number(js) = q_soap(js-ICUT_org) 
      enddo
      i = 0 
      do jesp = 1,N_aerosol_layers
         s = List_species(jesp)
         do js = ICUT_org+1, N_size
            i = i + 1
            ! Need to redistribute the mass even for non volatile particles because layers may have changed
            if ((aerosol_type(s)==4.and.soap_inorg_loc>=0).or. &
                 ((aerosol_type(s)==3.or.s==EH2O).and.(soap_inorg_loc==1.or.soap_inorg_loc==-1))) then
               concentration_mass(js, jesp) = q_soap(N_size + i - ICUT_org)
            endif           
         enddo        
      enddo
    
      do jesp = 1,N_aerosol
         if (inon_volatile(jesp).EQ.0) then
            if ((aerosol_type(jesp)==4.and.soap_inorg_loc>=0).or. &
                 (aerosol_type(jesp)==3.and.(soap_inorg_loc==1.or.soap_inorg_loc==-1))) then
               concentration_gas(jesp) = q_soap((N_size-ICUT_org)*(1+N_aerosol_layers) + jesp)
            endif
         else
            !Add initial gas concentration that were not added to soap
            if ((aerosol_type(jesp)==4.and.soap_inorg_loc>=0).or. &
                 (aerosol_type(jesp)==3.and.(soap_inorg_loc==1.or.soap_inorg_loc==-1))) then
               concentration_gas(jesp) = concentration_gas(jesp) + q_soap((N_size-ICUT_org)*(1+N_aerosol_layers) + jesp)
            endif
         endif
      enddo

      do js=ICUT_org+1,N_size
         proton_Nsize(js) = chp_Nsize(js) * lwc_Nsize(js) / 1.0e3
      enddo       

    END SUBROUTINE SSH_SOAP_DYN_ICUT
end Module kSOAP
