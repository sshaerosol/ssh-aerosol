C-----------------------------------------------------------------------
C     Copyright (C) 2001-2007, ENPC - INRIA - EDF R&D
C     Author(s): Denis Qu�lo
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

      SUBROUTINE chem (i1,i2,j1,j2,Nx,Ny,Nz,ns,nr,nrphot,nreactphot
     $     ,nemis,nemisspecies,convers_factor,convers_factor_jac
     $     ,nzemis,LWCmin,Wmol,ts,DLattenuation
     $     ,DLhumid,DLtemp,DLpress,DLCsourc,DLCphotolysis_rates,delta_t
     $     ,DLattenuationf,DLhumidf,DLtempf, DLpressf,DLCsourcf
     $     ,DLCphotolysis_ratesf,ncycle,dlon,dlat,DLconc,icld,iheter
     $     ,ns_aer,nbin_aer,DLLWC,bin_bound_aer,fixed_density_aer
     $     ,wet_diameter_aer,hetero_species_index,DLconc_aer
     $     ,option_adaptive_time_step, ATOL, tstep_min
     $     ,option_photolysis, jBiPER, kBiPER)

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
C     Nx: number of cells in X direction (longitude).
C     Ny: number of cells in Y direction (latitude).
C     Nz: number of vertical levels.
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
C     DLATTENUATION: 3D cloud attenuation field at initial time.
C     DLHUMID: 3D specific humidity field at initial time ([%]).
C     DLTEMP: 3D temperature field at initial time ([K]).
C     DLPRESS: 3D pressure field at initial time ([Pa]).
C     DLCSOURC: array of chemical volumic emissions at initial time
C     # ([\mu.g/m^3/s]).
C     DLCPHOTOLYSIS_RATES: photochemical kinetic rates
C     # at initial time ([s^{-1}]).
C     DLATTENUATIONf: 3D cloud attenuation field at final time.
C     DLHUMIDf: 3D specific humidity field at final time ([%]).
C     DLTEMPf: 3D temperature field at final time ([K]).
C     DLPRESSf: 3D pressure field at final time ([Pa]).
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
C     DLLWC: 3D air liquid water content ([fraction]).
C     BIN_BOUND_AER: aerosol diameters at bin bounds.
C     FIXED_DENSITY_AER: fixed aerosol density ([g/m^3]).
C     WET_DIAMETER_AER: 4D aerosol wet diameters ([\mu.m]).
C     HETERO_SPECIES_INDEX: index of gas species with heterogeneous reactions.
C     option_adaptive_time_step: 1 if adaptive time step.
C     ATOL = relative tolerance for deciding if the time step is kept.
C     tstep_min = minimum time step.
C
C     -- INPUT/OUTPUT VARIABLES
C
C     DLCONC: array of 3D gas concentrations ([\mu.g/m^3]).
C     DLCONC_AER: array of 3D aerosol concentrations ([\mu.g/m^3]).
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
C     Denis Qu�lo, CEREA, June 2001.
C
C------------------------------------------------------------------------

      IMPLICIT NONE

      INCLUDE 'CONST.INC'
      INCLUDE 'CONST_A.INC'

      DOUBLE PRECISION ts,delta_t
      DOUBLE PRECISION tschem,tfchem

      integer i1,i2,j1,j2,nx,ny,nz,ns,nr,nrphot,nemis,nzemis

      DOUBLE PRECISION DLconc(NX,NY,NZ,ns),ZC(NS)
      DOUBLE PRECISION DLtemp(NX,NY,NZ),DLtempf(NX,NY,NZ)
      DOUBLE PRECISION DLattenuation(NX,NY,NZ)
      DOUBLE PRECISION DLattenuationf(NX,NY,NZ)
      DOUBLE PRECISION DLhumid(NX,NY,NZ),DLhumidf(NX,NY,NZ)
      DOUBLE PRECISION DLCsourc(Nx,Ny,Nzemis,Nemis)
      DOUBLE PRECISION DLCsourcf(Nx,Ny,Nzemis,Nemis)
      DOUBLE PRECISION ZCsourc(NS)
      DOUBLE PRECISION ZCsourcf(NS)
      DOUBLE PRECISION DLRki(Nr),DLRkf(Nr)
      DOUBLE PRECISION DLpress(NX,NY,NZ),DLpressf(NX,NY,NZ)
      DOUBLE PRECISION DLCphotolysis_rates(Nx,Ny,Nz,NRphot)
      DOUBLE PRECISION DLCphotolysis_ratesf(Nx,Ny,Nz,NRphot)

      double precision dlon(nx),dlat(ny)

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

      INTEGER Jt,Ji,Jj,Jk,Jsp,i,Jb,Jemis

      INTEGER ns_aer,nbin_aer
      DOUBLE PRECISION DLLWC(NX,NY,NZ)
      DOUBLE PRECISION bin_bound_aer(nbin_aer + 1)
      DOUBLE PRECISION wet_diameter_aer(Nx,Ny,Nz,nbin_aer)
      INTEGER hetero_species_index(4)
      DOUBLE PRECISION DLconc_aer(NX,NY,NZ,nbin_aer,ns_aer)

      INTEGER ICLD,IHETER

      DOUBLE PRECISION wet_diameter_aer_loc(Nbin_aer)
      DOUBLE PRECISION granulo_aer(Nbin_aer), conc_tot
      DOUBLE PRECISION lwca(NX,NY,NZ)
      DOUBLE PRECISION LWCmin

      DOUBLE PRECISION MSF(nbin_aer)
      DOUBLE PRECISION DSF(nbin_aer)
      DOUBLE PRECISION MBF(nbin_aer+1)
      DOUBLE PRECISION DBF(nbin_aer+1)

      DOUBLE PRECISION fixed_density_aer,RHOA

      DOUBLE PRECISION DLk1(NS), DLk2(NS),ZC_old(NS)
      DOUBLE PRECISION EPSDLK
      PARAMETER (EPSDLK = 1.D-15)
      DOUBLE PRECISION supEdtstep, Edtstep(NS),ATOL
      DOUBLE PRECISION tstep,tstep_new,tstep_min,tfchem_tmp

      INTEGER option_adaptive_time_step
      INTEGER option_photolysis

C     To take into account BiPER degradation inside the particle.
      INTEGER jBiPER, kBiPER
      DOUBLE PRECISION saveDLBiPER(nbin_aer)

C     Aerosol density converted in microg / microm^3.
      RHOA = fixed_density_aer * 1.D-09

C     Aerosol discretization converted in microm.
      DO Jb=1,nbin_aer+1
         DBF(Jb) = bin_bound_aer(Jb) * 1.D06
         MBF(Jb) = RHOA * cst_pi6 * DBF(Jb)**3
      ENDDO

      DO Jb=1,nbin_aer
         MSF(Jb) = DSQRT(MBF(Jb) * MBF(Jb+1))
         DSF(Jb) = DSQRT(DBF(Jb) * DBF(Jb+1))
      ENDDO

C     Cloud liquid water content (g/m^3)
      DO Jk=1,NZ
         DO Jj=j1,j2
            DO Ji=i1,i2
               lwca(Ji,Jj,Jk)=DLLWC(Ji,Jj,Jk)*1000.d0*
     &              DLpress(Ji,Jj,Jk)/101325.d0*28.97d0/Pr/
     &              DLtemp(Ji,Jj,Jk)
            ENDDO
         ENDDO
      ENDDO

C     Loop on grid cells.
      DO Jk=1,NZ
         DO Jj=j1,j2
            DO Ji=i1,i2

               DO Jsp=1,Ns
                  ZCsourc(jsp)=0.D0
                  ZCsourcf(jsp)=0.d0
               ENDDO
C     Spatial extraction for volumic sources.
               if (nemis.gt.0) then
                  if (jk.le.nzemis) then
                     DO Jsp=1,Nemis
                        Jemis = nemisspecies(jsp)+1
                        ZCsourc(jemis) = DLCsourc(Ji,Jj,Jk,Jsp)
                        ZCsourcf(jemis)= DLCsourcf(Ji,Jj,Jk,Jsp)
                     ENDDO
                  else
                     DO Jsp=1,Nemis
                        ZCsourc(nemisspecies(jsp)+1)=0.D0
                        ZCsourcf(nemisspecies(jsp)+1)=0.D0
                     ENDDO
                  endif
               endif

C     Cloud attenuation.

               Zatt = DLattenuation(Ji, Jj, Jk)
               Zattf = DLattenuationf(Ji, Jj, Jk)

C     Projection.

               DO Jsp=1,Ns
                  ZC(Jsp) = DLconc(Ji,Jj,Jk,Jsp)
               ENDDO

C     Initialization of granulo for kinetic cte of heterogeneous rxns

               DO Jb=1,Nbin_aer
                  conc_tot=0.0D0
                  DO Jsp=1,Ns_aer-1
                     conc_tot = conc_tot + DLconc_aer(Ji,Jj,Jk,Jb,Jsp)
                  ENDDO
C     compute the particle number (mass/geometric mean diameter)
                  granulo_aer(Jb)=conc_tot / MSF(Jb)
                  wet_diameter_aer_loc(Jb) =
     &                 wet_diameter_aer(Ji,Jj,Jk,jb)
               ENDDO

C     Integration of chemistry (eventually with subcycling).

               DO Jt=1,Ncycle
                  tschem=ts+(Jt-1)*delta_t/Ncycle
                  tfchem=tschem+delta_t/Ncycle
C     Integration of chemistry with adaptive time stepping
                  tfchem_tmp = tfchem
                  tstep = delta_t/Ncycle

                  Do while (tschem.LT.tfchem)

!     Compute zenithal angles
                     DLmuzero=muzero(tschem,Dlon(ji),Dlat(jj))
                     Zangzen=dabs(DACOS(DLmuzero)*180.D0/PI)
                     DLmuzero=muzero(tfchem_tmp,Dlon(ji),Dlat(jj))
                     Zangzenf=dabs(DACOS(DLmuzero)*180.D0/PI)

C     If option_photolysis is 1,
C     photolytic reactions are calculated in kinetic.f

!     Compute kinetic rates
                     CALL Kinetic(Ns,Nbin_aer,Nr,
     s                    IHETER,ICLD,DLRKi,DLtemp(Ji,Jj,Jk),
     s                    DLhumid(Ji,Jj,Jk),DLpress(Ji,Jj,Jk),Zangzen,
     s                    Zatt,lwca(Ji,Jj,Jk),granulo_aer,
     s                    wet_diameter_aer_loc,DSF,
     s                    hetero_species_index,Wmol,LWCmin,
     s                    option_photolysis)
                     CALL Kinetic(Ns,Nbin_aer,Nr,
     s                    IHETER,ICLD,DLRKf,DLtempf(Ji,Jj,Jk),
     s                    DLhumidf(Ji,Jj,Jk),DLpressf(Ji,Jj,Jk),
     s                    Zangzenf,Zattf,lwca(Ji,Jj,Jk),granulo_aer,
     s                    wet_diameter_aer_loc,DSF,
     s                    hetero_species_index,Wmol,LWCmin,
     s                    option_photolysis)

C     If option_photolysis is 2,
C     photolytic reactions may be read.
                     IF (option_photolysis.eq.2) then
c     IF (Nrphot.GT.0) THEN
                        DO i=1,Nrphot
                           DLRKi(Nreactphot(i)+1) = Zatt *
     $                          DLCphotolysis_rates(Ji, Jj, Jk, i)
                           DLRKf(Nreactphot(i)+1) = Zattf *
     $                          DLCphotolysis_ratesf(Ji, Jj, Jk, i)
                        ENDDO
                     ENDIF

!     Solve gas-phase chemistry for the time step
                     CALL roschem (NS,Nr,nemis,ZC,ZCsourc,ZCsourcf,
     s                    convers_factor, convers_factor_jac,tschem,
     s                    tfchem_tmp,DLRki,DLRkf,ZC_old,DLK1,DLK2)
                     IF (jBiPER /= 0) THEN
                        DO Jb = 1, nbin_aer
                           saveDLBiPER(Jb) =
     s                          max(DLconc_aer(Ji, Jj, Jk, Jb, jBiPER)
     s                          - DLRKi(kBiPER) * tstep
     s                          * DLconc_aer(Ji, Jj, Jk, Jb, jBiPER),
     s                          0.0)
                        ENDDO
                     ENDIF

                     IF(option_adaptive_time_step.EQ.1) then
C     Check that the time step was ok
                        supEdtstep = 0.D0
                        Do Jsp = 1,Ns
                           If((DLK1(Jsp).GT.EPSDLK
     &                          .OR.DLK2(Jsp).GT.EPSDLK)
     &                          .AND.ZC(Jsp).GT.EPSDLK) then
!     Estimate the relative error
                              Edtstep(Jsp) = 0.5D0 *
     &                             dabs(DLk1(Jsp) + DLk2(Jsp))
     &                             / ZC(Jsp)
                              If(Edtstep(Jsp).GT.supEdtstep) then
                                 supEdtstep = Edtstep(Jsp)
                              Endif
                           Endif
                        Enddo
                        supEdtstep = supEdtstep/ATOL * tstep
                        If(supEdtstep.GT.1.D0
     &                       .AND.tstep.GT.tstep_min) then
!     The time step is rejected and the computation
!     is redone with a smaller time step
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
!     The time step is accepted and the time is incremented
                           tschem = tfchem_tmp
                           tfchem_tmp = tfchem_tmp + tstep
                           If(tfchem_tmp.GT.tfchem) then
                              tfchem_tmp = tfchem
                              tstep = tfchem_tmp - tschem
                           Endif
                           IF (jBiPER /= 0 ) THEN
C    TODO_PIGA Jb was not enumerated (it was 'nbin_aer+1' because of the last loop)
                              DO Jb=1,nbin_aer
                                 DLconc_aer(Ji, Jj, Jk, Jb, jBiPER) =
     &                                saveDLBiPER(Jb)
                              ENDDO
                           ENDIF
                        Endif
                     ELSE
                        tschem = tfchem
                        IF (jBiPER /= 0) THEN
C    TODO_PIGA Jb was not enumerated (it was 'nbin_aer+1' because of the last loop)
                           DO Jb=1,nbin_aer
                              DLconc_aer(Ji, Jj, Jk, Jb, jBiPER) =
     &                             saveDLBiPER(Jb)
                           ENDDO
                        ENDIF
                     ENDIF

                  Enddo         !End loop Do while for time stepping

               ENDDO

C     Storage in the 3D array of chemical concentrations.

               DO Jsp=1,NS
                  DLconc(Ji,Jj,Jk,Jsp) = ZC(Jsp)
               ENDDO

            ENDDO               ! loop x.
         ENDDO                  ! loop y.
      ENDDO                     ! loop z.

      END
