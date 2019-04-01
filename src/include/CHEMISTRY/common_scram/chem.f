C-----------------------------------------------------------------------
C     Copyright (C) 2001-2007, ENPC - INRIA - EDF R&D
C     Author(s): Denis Quélo
C     
C     This file is part of the Size Resolved Aerosol Model (SIREAM), a
C     component of the air quality modeling system Polyphemus.
C     
C     Polyphemus is developed in the INRIA - ENPC joint project-team
C     CLIME and in the ENPC - EDF R&D joint laboratory CEREA.
C     
C     Polyphemus is free software; you can redistribute it and/or modify
C     it under the terms of the GNU General Public License as published
C     by the Free Software Foundation; either version 2 of the License,
C     or (at your option) any later version.
C     
C     Polyphemus is distributed in the hope that it will be useful, but
C     WITHOUT ANY WARRANTY; without even the implied warranty of
C     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
C     General Public License for more details.
C     
C     For more information, visit the Polyphemus web site:
C     http://cerea.enpc.fr/polyphemus/
C-----------------------------------------------------------------------

      SUBROUTINE chem (ns,nr,nrphot,nreactphot
     $     ,nemis,nemisspecies,convers_factor,convers_factor_jac
     $     ,nzemis,LWCmin,Wmol,ts,DLattenuation
     $     ,DLhumid,DLtemp,DLpress,DLCsourc,DLCphotolysis_rates,delta_t
     $     ,DLattenuationf,DLhumidf,DLtempf, DLpressf,DLCsourcf
     $     ,DLCphotolysis_ratesf,ncycle,dlon,dlat,DLconc,icld,iheter
     $     ,ns_aer,nbin_aer,ncomp_aer,DLLWC,bin_bound_aer
     $     ,fixed_density_aer,wet_diameter_aer,hetero_species_index
     $     ,DLconc_aer,option_adaptive_time_step, ATOL, tstep_min
     $     ,option_photolysis,jBiPER,kBiPER,INUM,IDENS
     $     ,DLnumconc_aer,mass_density_aer)

C------------------------------------------------------------------------
C     
C     -- DESCRIPTION 
C     
C     This routine computes one timestep for gas-phase chemistry RADM.
C     Chemical kinetics is solved in each grid cell.
C     
C------------------------------------------------------------------------
C     
C     -- INPUT VARIABLES
C
C     NS: number of gas species.
C     NR: number of gas reactions.
C     NRPHOT: number of photolytic reactions.
C     NREACTPHOT: index of photolysis reactions.
C     NEMIS: number of gas with emissions.
C     NEMISSPECIES: index of gas species with emissions.
C     NZEMIS: highest vertical level for gas volumic emissions.
C     LWCmin: air liquid water content threshold.
C     Wmol: molar mass of gas species ([g/mol]).
C     TS: initial time (GMT, computed from January 1st, [s]).
C     DLATTENUATION: cloud attenuation field at initial time.
C     DLHUMID: specific humidity field at initial time ([%]).
C     DLTEMP: temperature field at initial time ([K]).
C     DLPRESS: pressure field at initial time ([Pa]).
C     DLCSOURC: array of chemical volumic emissions at initial time
C     # ([\mu.g/m^3/s]).
C     DLCPHOTOLYSIS_RATES: photochemical kinetic rates
C     # at initial time ([s^{-1}]).
C     DLATTENUATIONf: cloud attenuation field at final time.
C     DLHUMIDf: specific humidity field at final time ([%]).
C     DLTEMPf: temperature field at final time ([K]).
C     DLPRESSf: pressure field at final time ([Pa]).
C     DLCSOURCf: array of chemical volumic emissions at final time
C     # ([\mu.g/m^3/s]).
C     DLCPHOTOLYSIS_RATESf: photochemical kinetic rates
C     # at final time ([s^{-1}]).
C     DELTA_T: time step ([s]).
C     NCYCLE: number of cycle in chemistry computation.
C     DLON: longitude discretization.
C     DLAT: latitude discretization.
C     ICLD: take cloud into account (1).
C     IHETER: take heterogeneous reactions into account (1).
C     NS_AER: number of aerosol species.
C     NBIN_AER: number of aerosol bins.
C     DLLWC: air liquid water content ([fraction]).
C     BIN_BOUND_AER: aerosol diameters at bin bounds.
C     FIXED_DENSITY_AER: fixed aerosol density ([g/m^3]).
C     WET_DIAMETER_AER: aerosol wet diameters array ([\mu.m]).
C     HETERO_SPECIES_INDEX: index of gas species with heterogeneous reactions.
C     option_adaptive_time_step: 1 if adaptive time step.
C     ATOL = relative tolerance for deciding if the time step is kept.
C     tstep_min = minimum time step.
C     INUM = 1 if the number concentration is followed in the CTM
C     = 0 if the number concentration is not followed.
C     IDENS = 1 for varying aerosol density
C     = 0 for fixed aerosol density.
C     
C     -- INPUT/OUTPUT VARIABLES
C
C     DLCONC: gas concentrations ([\mu.g/m^3]).
C     DLCONC_AER: aerosol concentrations ([\mu.g/m^3]).
C     # Before entry, it is given at initial time of the timestep.
C     # On exit, it is computed at final time of the timestep.
C     
C     -- OUTPUT VARIABLES
C     
C------------------------------------------------------------------------
C     
C     -- REMARKS
C     
C------------------------------------------------------------------------
C     
C     -- MODIFICATIONS
C     
C     2002/02/26: new treatment of sources (Jaouad Boutahar, CEREA).
C     2006/09/29: updated header (Edouard Debry).
C     
C     2009/01/22: added adaptatif time stepping (K. Sartelet, CEREA)
C     2010/03/05: added option_photolysis index (Youngseob KIM)
C     
C------------------------------------------------------------------------
C     
C     -- AUTHOR(S)
C     
C     Denis Quélo, CEREA, June 2001.
C     
C------------------------------------------------------------------------

      IMPLICIT NONE

      INCLUDE 'CONST.INC'
      INCLUDE 'CONST_A.INC'
      INCLUDE 'paraero.inc'

      DOUBLE PRECISION ts,delta_t
      DOUBLE PRECISION tschem,tfchem

      integer i1,i2,j1,j2,nx,ny,nz,ns,nr,nrphot,nemis,nzemis

      DOUBLE PRECISION DLconc(ns),ZC(NS)
      DOUBLE PRECISION DLtemp,DLtempf
      DOUBLE PRECISION DLattenuation
      DOUBLE PRECISION DLattenuationf
      DOUBLE PRECISION DLhumid,DLhumidf
      DOUBLE PRECISION DLCsourc(Nemis)
      DOUBLE PRECISION DLCsourcf(Nemis)
      DOUBLE PRECISION ZCsourc(NS)
      DOUBLE PRECISION ZCsourcf(NS)
      DOUBLE PRECISION DLRki(Nr),DLRkf(Nr)
      DOUBLE PRECISION DLpress,DLpressf
      DOUBLE PRECISION DLCphotolysis_rates(NRphot)
      DOUBLE PRECISION DLCphotolysis_ratesf(NRphot)

      double precision dlon,dlat

      integer ncycle
      double precision convers_factor(ns)
      double precision convers_factor_jac(ns,ns)
      double precision Wmol(ns)

      DOUBLE PRECISION Zangzen,Zangzenf
      DOUBLE PRECISION Zatt,Zattf

      DOUBLE PRECISION muzero,DLmuzero
      EXTERNAL muzero

      integer nreactphot(nrphot)
      integer nemisspecies(nemis)

      INTEGER Jt,Jsp,i,Jb,Jemis

      INTEGER ns_aer,nbin_aer,ncomp_aer
      DOUBLE PRECISION DLLWC
      DOUBLE PRECISION bin_bound_aer(nbin_aer + 1)
      DOUBLE PRECISION wet_diameter_aer(nbin_aer)
      INTEGER hetero_species_index(4)
      DOUBLE PRECISION DLconc_aer(nbin_aer,ns_aer)
      DOUBLE PRECISION DLnumconc_aer(nbin_aer)

      INTEGER ICLD,IHETER,INUM,IDENS

c     DOUBLE PRECISION wet_diameter_aer_loc(Nbin_aer)
      DOUBLE PRECISION granulo_aer(Nbin_aer), conc_tot
      DOUBLE PRECISION lwca
      DOUBLE PRECISION LWCmin

      DOUBLE PRECISION MSF(nbin_aer)
      DOUBLE PRECISION DSF(nbin_aer)
      DOUBLE PRECISION DBF((nbin_aer/ncomp_aer)+1)!SZ
      INTEGER idx_bs(nbin_aer)

      DOUBLE PRECISION fixed_density_aer,RHOA
      double precision rho_dry(nbin_aer)
      double precision mass_density_aer(ns_aer)

      DOUBLE PRECISION DLk1(NS), DLk2(NS),ZC_old(NS)
      DOUBLE PRECISION EPSDLK
      PARAMETER (EPSDLK = 1.D-15)
      DOUBLE PRECISION supEdtstep, Edtstep(NS),ATOL
      DOUBLE PRECISION tstep,tstep_new,tstep_min,tfchem_tmp  

      INTEGER option_adaptive_time_step
      INTEGER option_photolysis

C     To take into account BiPER degradation inside the particle
      INTEGER jBiPER,kBiPER
      DOUBLE PRECISION saveDLBiPER(nbin_aer)


       double precision tot_before(Nbin_aer), tot_after(Nbin_aer)

C     Aerosol density converted in microg / microm^3.
      RHOA = fixed_density_aer * 1.D-09
C     Aerosol discretization converted in microm.
      DO Jb=1,(nbin_aer/ncomp_aer)+1!sz
         DBF(Jb) = bin_bound_aer(Jb) * 1.D06
      ENDDO

      DO Jb=1,nbin_aer
	  idx_bs(Jb)=(Jb-1)/ncomp_aer+1!relations between bin idx and size idx
      ENDDO
      
C     With real number concentration.
      IF (INUM.EQ.1) THEN
C     Compute aerosol density
         rho_dry = RHOA
         IF (IDENS.EQ.1) THEN   ! for varying density
            DO Jb=1,nbin_aer
               CALL COMPUTE_DENSITY(nbin_aer,ns_aer, ns_aer, TINYM,
     &              DLconc_aer,
     &              mass_density_aer,Jb,rho_dry(Jb))
            ENDDO
         ENDIF
         DO Jb = 1, nbin_aer
            conc_tot = 0.d0
            DO Jsp = 1, Ns_aer
               conc_tot = conc_tot + DLconc_aer(Jb,Jsp)
            ENDDO
C     Compute mass and diameter of each section
            IF (DLnumconc_aer(Jb) .GT. 0.d0) THEN
               MSF(Jb) = conc_tot/DLnumconc_aer(Jb)
            ELSE
               MSF(Jb) = 0.d0
            ENDIF

            if ((DLnumconc_aer(Jb).GT. TINYN .or. 
     s           conc_tot.GT.TINYM)
     s           .AND. IDENS .EQ. 1) then
               DSF(Jb) = (MSF(Jb)/cst_PI6/rho_dry(Jb))**cst_FRAC3

            else
               DSF(Jb) = DSQRT(DBF(idx_bs(Jb))* DBF(idx_bs(Jb)+1))!sz
            endif

            if (DSF(Jb) .LT. DBF(idx_bs(Jb)) .or. 
     &          DSF(Jb) .GT. DBF(idx_bs(Jb)+1)) THEN
               DSF(Jb) =  DSQRT(DBF(idx_bs(Jb)) * DBF(idx_bs(Jb)+1))
            endif

         ENDDO

      ELSE


         DO Jb = 1, nbin_aer
            DSF(Jb) = DSQRT(DBF(idx_bs(Jb)) * DBF(idx_bs(Jb)+1))!sz
            MSF(Jb) = RHOA * cst_pi6 * DSF(Jb)**3
         ENDDO
      ENDIF

C     Cloud liquid water content (g/m^3)
      lwca=DLLWC*1000.d0*
     &     DLpress/101325.d0*28.97d0/Pr/
     &     DLtemp

C     Loop on grid cells.
      DO Jsp=1,Ns
         ZCsourc(jsp)=0.D0
         ZCsourcf(jsp)=0.d0
      ENDDO
C     Spatial extraction for volumic sources.
      if (nemis.gt.0) then
c     if (jk.le.nzemis) then
         DO Jsp=1,Nemis
            Jemis = nemisspecies(jsp)+1
            ZCsourc(jemis) = DLCsourc(Jsp)
            ZCsourcf(jemis)= DLCsourcf(Jsp)
         ENDDO
c     else
c     DO Jsp=1,Nemis
c     ZCsourc(nemisspecies(jsp)+1)=0.D0
c     ZCsourcf(nemisspecies(jsp)+1)=0.D0
c     ENDDO
c     endif
      endif

C     Cloud attenuation.
      
      Zatt = DLattenuation
      Zattf = DLattenuationf
      
C     Projection.
      DO Jsp=1,Ns
         ZC(Jsp) = DLconc(Jsp)
      ENDDO

C     Initialization of granulo for kinetic cte of heterogeneous rxns
      
      DO Jb=1,Nbin_aer
         conc_tot=0.0D0
         tot_before(Jb) = 0.D0
         DO Jsp=1,Ns_aer-1
            conc_tot = conc_tot + DLconc_aer(Jb,Jsp) 
            tot_before(Jb) = tot_before(Jb) +  DLconc_aer(Jb,Jsp) 
         ENDDO
C     compute the particle number (mass/geometric mean diameter)
         IF (INUM.EQ.1) THEN
            granulo_aer(Jb) = DLnumconc_aer(Jb)
         ELSE
            if (MSF(Jb) .GT. 0.D0) THEN
               granulo_aer(Jb) = conc_tot/MSF(Jb)               
            else
               write(*,*) , "chem.f : error "
               stop
            endif
         ENDIF
c     wet_diameter_aer_loc(Jb) =
c     &                 wet_diameter_aer(jb)
      ENDDO

C     Integration of chemistry (eventually with subcycling).

      DO Jt=1,Ncycle
         tschem=ts+(Jt-1)*delta_t/Ncycle
         tfchem=tschem+delta_t/Ncycle
C     Integration of chemistry with adaptive time stepping
                  tfchem_tmp = tfchem
                  tstep = delta_t/Ncycle

         Do while (tschem.LT.tfchem) 
! Compute zenithal angles
            DLmuzero=muzero(tschem,Dlon,Dlat)
            Zangzen=dabs(DACOS(DLmuzero)*180.D0/PI)
            DLmuzero=muzero(tfchem_tmp,Dlon,Dlat)
            Zangzenf=dabs(DACOS(DLmuzero)*180.D0/PI)

C     If option_photolysis is 1,
C     photolytic reactions are calculated in kinetic.f

! Compute kinetic rates
            CALL Kinetic(Ns,Nbin_aer,Nr,
     s           IHETER,ICLD,DLRKi,DLtemp,
     s           DLhumid,DLpress,Zangzen,
     s           Zatt,lwca,granulo_aer,
     s           wet_diameter_aer,DSF,
     s           hetero_species_index,Wmol,LWCmin,
     s           option_photolysis)
            CALL Kinetic(Ns,Nbin_aer,Nr,
     s           IHETER,ICLD,DLRKf,DLtempf,
     s           DLhumidf,DLpressf,
     s           Zangzenf,Zattf,lwca,granulo_aer,
     s           wet_diameter_aer,DSF,
     s           hetero_species_index,Wmol,LWCmin,
     s           option_photolysis)

C     If option_photolysis is 2,
C     photolytic reactions may be read.
            IF (option_photolysis.eq.2) then
c     IF (Nrphot.GT.0) THEN
               DO i=1,Nrphot
                  DLRKi(Nreactphot(i)+1) = Zatt *
     $                 DLCphotolysis_rates(i)
                  DLRKf(Nreactphot(i)+1) = Zattf *
     $                 DLCphotolysis_ratesf(i)
               ENDDO
            ENDIF
            
            CALL roschem (NS,Nr,nemis,ZC,ZCsourc,ZCsourcf,
     s           convers_factor, convers_factor_jac,tschem,
     s           tfchem_tmp,DLRki,DLRkf,ZC_old,DLK1,DLK2)
            
            IF (jBiPER/=0) THEN
               DO Jb=1,nbin_aer
                  saveDLBiPER(Jb)=max(DLconc_aer(Jb,
     s                 jBiPER)
     s                 -DLRKi(kBiPER)*tstep*
     s                 DLconc_aer(Jb,jBiPER),0.0)
               ENDDO
            ENDIF
            
            IF(option_adaptive_time_step.EQ.1) then
C     Check that the time step was ok
               supEdtstep = 0.D0
               Do Jsp = 1,Ns
                  If((DLK1(Jsp).GT.EPSDLK
     &                 .OR.DLK2(Jsp).GT.EPSDLK)
     &                 .AND.ZC(Jsp).GT.EPSDLK) then
!     Estimate the relative error
                     Edtstep(Jsp) = 0.5D0 * 
     &                    dabs(DLk1(Jsp) + DLk2(Jsp))
     &                    / ZC(Jsp)
                     If(Edtstep(Jsp).GT.supEdtstep) then 
                        supEdtstep = Edtstep(Jsp)
                     Endif
                  Endif
               Enddo
               supEdtstep = supEdtstep/ATOL * tstep
c     supEdtstep = supEdtstep/ATOL
               If(supEdtstep.GT.1.D0
     &              .AND.tstep.GT.tstep_min) then
! The time step is rejected and the computation
! is redone with a smaller time step
                  tstep_new = tstep * 0.9d0 /dsqrt(supEdtstep)
                  tstep_new = DMAX1(tstep_new,tstep_min)
                  tfchem_tmp = tschem + tstep_new
                  If(tfchem_tmp.GT.tfchem) then
                     tfchem_tmp = tfchem
                     tstep_new = tfchem_tmp - tschem
                  Endif
                  tstep = tstep_new
                  Do Jsp=1,NS
                     ZC(Jsp) = ZC_old(Jsp)
                  Enddo
               Else
! The time step is accepted and the time is incremented
                  tschem = tfchem_tmp
                  tfchem_tmp = tfchem_tmp + tstep
                  If(tfchem_tmp.GT.tfchem) then
                     tfchem_tmp = tfchem
                     tstep = tfchem_tmp - tschem
                  Endif
                  IF(jBiPER/=0) THEN
                      DO Jb=1,nbin_aer
                     DLconc_aer(Jb,jBiPER)=saveDLBiPER(Jb)
                      ENDDO
                  ENDIF
               Endif
            ELSE
               tschem = tfchem
               IF(jBiPER/=0) THEN
                      DO Jb=1,nbin_aer
                  DLconc_aer(Jb,jBiPER)=saveDLBiPER(Jb)
                      ENDDO
               ENDIF  
            ENDIF

         Enddo                  !End loop Do while for time stepping

      ENDDO

C     Storage in the array of chemical concentrations.

               DO Jsp=1,NS
                  DLconc(Jsp) = ZC(Jsp)
                  ! NAN detection algorithm used
                  if (dlconc(Jsp) .ne. dlconc(Jsp)) then
                     write(*,*) "From chem 0D function:"
                     write(*,*) Jsp, dlconc(Jsp)
                     stop
                  endif
               ENDDO

      DO Jb=1,Nbin_aer
         DO Jsp=1,Ns_aer
            if (dlconc_aer(jb,jsp) .ne. dlconc_aer(jb,jsp)) then
               write(*,*) "From chem 0D function (aerosol conc.):"
               write(*,*) jb,jsp,dlconc_aer(jb,jsp)
               stop
            endif
         enddo
      enddo

      END
