      
      SUBROUTINE ssh_chem (ns,nr,nrphot,nreactphot
     $     ,nemis,nemisspecies,convers_factor,convers_factor_jac
     $     ,nzemis,LWCmin,Wmol,ts,DLattenuation
     $     ,DLhumid,DLtemp,DLpress,DLCsourc,DLCphotolysis_rates,delta_t
     $     ,DLattenuationf,DLhumidf,DLtempf, DLpressf,DLCsourcf
     $     ,DLCphotolysis_ratesf,dlon,dlat,DLconc,icld,iheter
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

      DOUBLE PRECISION ssh_muzero,DLmuzero
      EXTERNAL ssh_muzero

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
      DOUBLE PRECISION DBF((nbin_aer/ncomp_aer)+1) !SZ
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
      DO Jb=1,(nbin_aer/ncomp_aer)+1 !sz
         DBF(Jb) = bin_bound_aer(Jb) * 1.D06
      ENDDO

      DO Jb=1,nbin_aer
         idx_bs(Jb)=(Jb-1)/ncomp_aer+1 !relations between bin idx and size idx
      ENDDO

C     With real number concentration.
      IF (INUM.EQ.1) THEN
C     Compute aerosol density
         rho_dry = RHOA
         IF (IDENS.EQ.1) THEN   ! for varying density
            DO Jb=1,nbin_aer
               CALL SSH_COMPUTE_DENSITY(nbin_aer,ns_aer, ns_aer, TINYM,
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
               DSF(Jb) = DSQRT(DBF(idx_bs(Jb))* DBF(idx_bs(Jb)+1)) !sz
            endif

            if (DSF(Jb) .LT. DBF(idx_bs(Jb)) .or. 
     &           DSF(Jb) .GT. DBF(idx_bs(Jb)+1)) THEN
               DSF(Jb) =  DSQRT(DBF(idx_bs(Jb)) * DBF(idx_bs(Jb)+1))
            endif

         ENDDO

      ELSE


         DO Jb = 1, nbin_aer
            DSF(Jb) = DSQRT(DBF(idx_bs(Jb)) * DBF(idx_bs(Jb)+1)) !sz
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
               write(*,*) "chem.f : error "
               stop
            endif
         ENDIF
c     wet_diameter_aer_loc(Jb) =
c     &                 wet_diameter_aer(jb)
      ENDDO

C     Integration of chemistry (eventually with subcycling).

      
      tschem=ts
      tfchem=tschem+delta_t
C     Integration of chemistry with adaptive time stepping
      tfchem_tmp = tfchem
      tstep = delta_t

      Do while (tschem.LT.tfchem) 
!     Compute zenithal angles
         DLmuzero=ssh_muzero(tschem,Dlon,Dlat)
         Zangzen=dabs(DACOS(DLmuzero)*180.D0/PI)
         DLmuzero=ssh_muzero(tfchem_tmp,Dlon,Dlat)
         Zangzenf=dabs(DACOS(DLmuzero)*180.D0/PI)

C     If option_photolysis is 1,
C     photolytic reactions are calculated in kinetic.f

!     Compute kinetic rates
         CALL SSH_Kinetic(Ns,Nbin_aer,Nr,
     s        IHETER,ICLD,DLRKi,DLtemp,
     s        DLhumid,DLpress,Zangzen,
     s        Zatt,lwca,granulo_aer,
     s        wet_diameter_aer,DSF,
     s        hetero_species_index,Wmol,LWCmin,
     s        option_photolysis,0d0)
!     CALL Kinetic(Ns,Nbin_aer,Nr,
!     s           IHETER,ICLD,DLRKf,DLtempf,
!     s           DLhumidf,DLpressf,
!     s           Zangzenf,Zattf,lwca,granulo_aer,
!     s           wet_diameter_aer,DSF,
!     s           hetero_species_index,Wmol,LWCmin,
!     s           option_photolysis,RO2)
         DO i=1,Nr
            DLRKf(i) = DLRKi(i) 
         ENDDO

C     If option_photolysis is 2,
C     photolytic reactions may be read.
         IF (option_photolysis.eq.2) then
            DO i=1,Nrphot
               DLRKi(Nreactphot(i)) = Zatt *
     $              DLCphotolysis_rates(i)
               DLRKf(Nreactphot(i)) = DLRKi(Nreactphot(i)) 
            ENDDO
         ENDIF
         CALL SSH_roschem (NS,Nr,nemis,ZC,ZCsourc,ZCsourcf,
     s        convers_factor, convers_factor_jac,tschem,
     s        tfchem_tmp,DLRki,DLRkf,ZC_old,DLK1,DLK2)

!     IF (jBiPER/=0) THEN
!     DO Jb=1,nbin_aer
!     saveDLBiPER(Jb)=max(DLconc_aer(Jb,
!     s                 jBiPER)
!     s                 -DLRKi(kBiPER)*tstep*
!     s                 DLconc_aer(Jb,jBiPER),0.0)
!     ENDDO
!     ENDIF
         
         IF(option_adaptive_time_step.EQ.1) then
C     Check that the time step was ok
            supEdtstep = 0.D0
            Do Jsp = 1,Ns
               If((DLK1(Jsp).GT.EPSDLK
     &              .OR.DLK2(Jsp).GT.EPSDLK)
     &              .AND.ZC(Jsp).GT.EPSDLK) then
!     Estimate the relative error
                  Edtstep(Jsp) = 0.5D0 * 
     &                 dabs(DLk1(Jsp) + DLk2(Jsp))
     &                 / ZC(Jsp)
                  If(Edtstep(Jsp).GT.supEdtstep) then 
                     supEdtstep = Edtstep(Jsp)
                  Endif
               Endif
            Enddo
            supEdtstep = supEdtstep/ATOL * tstep
c     supEdtstep = supEdtstep/ATOL
            If(supEdtstep.GT.1.D0
     &           .AND.tstep.GT.tstep_min) then
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
!     IF(jBiPER/=0) THEN
!     DO Jb=1,nbin_aer
!     DLconc_aer(Jb,jBiPER)=saveDLBiPER(Jb)
!     ENDDO
!     ENDIF
            Endif
         ELSE
            tschem = tfchem
!     IF(jBiPER/=0) THEN
!     DO Jb=1,nbin_aer
!     DLconc_aer(Jb,jBiPER)=saveDLBiPER(Jb)
!     ENDDO
!     ENDIF  
         ENDIF

      Enddo                     !End loop Do while for time stepping

      

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


      
      SUBROUTINE ssh_chem_twostep (ns,nr,nrphot,nreactphot
     $     ,nemis,nemisspecies,convers_factor,convers_factor_jac
     $     ,nzemis,LWCmin,Wmol,ts,DLattenuation
     $     ,DLhumid,DLtemp,DLpress,DLCsourc,DLCphotolysis_rates,delta_t
     $     ,DLattenuationf,DLhumidf,DLtempf, DLpressf,DLCsourcf
     $     ,DLCphotolysis_ratesf,dlon,dlat,DLconc,icld,iheter
     $     ,ns_aer,nbin_aer,ncomp_aer,DLLWC,bin_bound_aer
     $     ,fixed_density_aer,wet_diameter_aer,hetero_species_index
     $     ,DLconc_aer,option_adaptive_time_step, TOL, tstep_min
     $     ,option_photolysis,jBiPER,kBiPER,INUM,IDENS
     $     ,DLnumconc_aer,mass_density_aer
     $     ,ncst_gas, cst_gas, cst_gas_index !genoa use constant gas conc.
     $     ,tag_RO2, nRO2_chem, iRO2, iRO2_cst, RO2index ! for RO2-RO2 reaction
     $     ,aerosol_species_interact,keep_gp, DLwall
     $     ,kwall_gas,kwall_particle,Cwall,aerosol_type,psat,dhvap,Tref)

C------------------------------------------------------------------------
C     
C     -- DESCRIPTION 
C     
C     This routine computes the chemistry with the 
C     TWO-STEP time numerical solver. It is based on
C     the application of a Gauss-Seidel iteration 
C     scheme to the two-step implicit backward 
C     differentiation (BDF2) formula.
C     
C     Reference: Verwer J. (1994). Gauss-Seidel 
C     iteration for stiff odes from chemical kinetics.
C     Journal on Scientific Computing, 15:1243-1250
C     
C------------------------------------------------------------------------


      IMPLICIT NONE

      INCLUDE 'CONST.INC'
      INCLUDE 'CONST_A.INC'
      INCLUDE 'paraero.inc'

      DOUBLE PRECISION ts,delta_t

      integer i1,i2,j1,j2,nx,ny,nz,ns,nr,nrphot,nemis,nzemis,Jsp2

      DOUBLE PRECISION DLconc(ns),ZC(NS),DLwall(ns),ZCwall(ns)
      DOUBLE PRECISION DLtemp,DLtempf
      DOUBLE PRECISION DLattenuation
      DOUBLE PRECISION DLattenuationf
      DOUBLE PRECISION DLhumid,DLhumidf
      DOUBLE PRECISION DLCsourc(Nemis)
      DOUBLE PRECISION DLCsourcf(Nemis)
      DOUBLE PRECISION ZCsourc(NS)
      DOUBLE PRECISION ZCsourcf(NS)
      DOUBLE PRECISION kwall_gas, kwall_particle, Cwall

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

      DOUBLE PRECISION ssh_muzero,DLmuzero
      EXTERNAL ssh_muzero

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

      DOUBLE PRECISION granulo_aer(Nbin_aer), conc_tot
      DOUBLE PRECISION lwca
      DOUBLE PRECISION LWCmin
      DOUBLE PRECISION wloss

      DOUBLE PRECISION MSF(nbin_aer)
      DOUBLE PRECISION DSF(nbin_aer)
      DOUBLE PRECISION DBF((nbin_aer/ncomp_aer)+1)
      INTEGER idx_bs(nbin_aer)

      DOUBLE PRECISION fixed_density_aer,RHOA
      double precision rho_dry(nbin_aer)
      double precision mass_density_aer(ns_aer)

C     Parameters initialized for the two-step solver
      integer m,j               !m nombre d'iteration
      DOUBLE PRECISION DLRki(Nr),DLRkf(Nr),dw(Nr,Ns)
      DOUBLE PRECISION tschem,tfchem,dun,dzero
      DOUBLE PRECISION alpha,TOL,rtol,atol
      PARAMETER (alpha=5.d0)
      PARAMETER (dun=1.d0)
      PARAMETER (dzero=0.d0)

      DOUBLE PRECISION tstep,tstep_min
      DOUBLE PRECISION chpr0(NS),chlo0(NS),conci(NS),concii(NS)
      DOUBLE PRECISION concz(NS),conczz(NS)
      DOUBLE PRECISION :: wk,dtnsave,error_max,c,gam,ratloss

      INTEGER option_adaptive_time_step
      INTEGER option_photolysis

C     To take into account BiPER degradation inside the particle
      INTEGER jBiPER,kBiPER
      DOUBLE PRECISION saveDLBiPER(nbin_aer)

C     genoa for RO2-RO2 reaction and constant profile
      INTEGER ncst_gas,nRO2_chem,iRO2,iRO2_cst,tag_RO2
      INTEGER cst_gas_index(ncst_gas), RO2index(nRO2_chem)
      DOUBLE PRECISION RO2,cst_gas(ncst_gas),ZC_cst(ncst_gas)

C     genoa keep_gp: gas-particle partitioning
      integer s, keep_gp        ! used for genoa with large timestep
      integer aerosol_species_interact(Ns),aerosol_type(ns_aer)
      DOUBLE PRECISION toadd, ZCtot_save(ns), psat(ns_aer)
      DOUBLE PRECISION kpart(ns_aer),psat_loc,dhvap(ns_aer),tref(ns_aer)
      DOUBLE PRECISION ratio_gas(ns),ZCtot(NS),DLconc_save(ns)

C     Aerosol density converted in microg / microm^3.
      RHOA = fixed_density_aer * 1.D-09
C     Aerosol discretization converted in microm.
      DO Jb=1,(nbin_aer/ncomp_aer)+1
         DBF(Jb) = bin_bound_aer(Jb) * 1.D06
      ENDDO

C     relations between bin idx and size idx
      DO Jb=1,nbin_aer
         idx_bs(Jb)=(Jb-1)/ncomp_aer+1
      ENDDO

      IF (kwall_gas>0.d0) THEN
         DO Jsp=1,ns_aer
            IF (aerosol_type(Jsp)==4 .and. psat(Jsp)>0) THEN
               if (Tref(Jsp)>0) then
                  psat_loc=psat(Jsp) !*
     &                 *exp(-1000.0*dhvap(Jsp)/8.314d0*
     &                 (1.0/DLtemp-1.0/Tref(Jsp)))
               else 
                  psat_loc=psat(Jsp) !*
     &                 *exp(-1000.0*dhvap(Jsp)/8.314d0*
     &                 (1.0/DLtemp-1.0/298.d0))
               endif
               kpart(Jsp)=760.d0*8.202d-5*DLtemp/(1.d6*200.d0*psat_loc)
               kpart(Jsp)=max(kpart(Jsp),kwall_gas/Cwall*tstep_min*5)
!     print*,Jsp,kpart(Jsp),psat_loc,dhvap(Jsp),psat(Jsp)
            ENDIF
         ENDDO
      ENDIF

C     With real number concentration.
      IF (INUM.EQ.1) THEN
C     Compute aerosol density
         rho_dry = RHOA
         IF (IDENS.EQ.1) THEN   ! for varying density
            DO Jb=1,nbin_aer
               CALL SSH_COMPUTE_DENSITY(nbin_aer,ns_aer, ns_aer, TINYM,
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
               DSF(Jb) = DSQRT(DBF(idx_bs(Jb))* DBF(idx_bs(Jb)+1)) !sz
            endif

            if (DSF(Jb) .LT. DBF(idx_bs(Jb)) .or. 
     &           DSF(Jb) .GT. DBF(idx_bs(Jb)+1)) THEN
               DSF(Jb) =  DSQRT(DBF(idx_bs(Jb)) * DBF(idx_bs(Jb)+1))
            endif

         ENDDO

      ELSE
         DO Jb = 1, nbin_aer
            DSF(Jb) = DSQRT(DBF(idx_bs(Jb)) * DBF(idx_bs(Jb)+1)) !sz
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
      endif

C     Cloud attenuation.
      Zatt = DLattenuation
      Zattf = DLattenuationf
      
C     Projection.
!     Conversion mug/m3 to molecules/cm3.
      DO Jsp=1,Ns
         ZC(Jsp) = DLconc(Jsp)* convers_factor(Jsp)
         ZCwall(Jsp) = DLwall(Jsp)* convers_factor(Jsp)
! genoa keep_gp
         DLconc_save(Jsp)=DLconc(Jsp)
         ZCtot(Jsp)=ZC(Jsp)
         ratio_gas(Jsp)=1.d0
      ENDDO        

! genoa keep_gp
      if (keep_gp==1) then
         do s = 1, Ns_aer
            Jsp=aerosol_species_interact(s)
!print*,s,Jsp
            if (Jsp>0.) then
               conc_tot=0.0D0
               DO Jb=1,Nbin_aer
                  conc_tot = conc_tot + DLconc_aer(Jb,s) 
               ENDDO
               
               ZCtot(Jsp) = (DLconc(Jsp)+conc_tot)*convers_factor(Jsp)
               if ( ZCtot(Jsp)>0.d0) then
                  ratio_gas(Jsp)= ZC(Jsp)/ZCtot(Jsp)
               else
                  ratio_gas(Jsp)= 1.d0
               endif
            endif
         ENDDO
      endif
!print*,maxval(ratio_gas)

C     !  input constant concentration genoa
      if (ncst_gas.gt.0) then
         do i1=1,ncst_gas
            ZC_cst(i1)=cst_gas(i1)* convers_factor(cst_gas_index(i1))
            ZC(cst_gas_index(i1))=ZC_cst(i1)
! genoa keep_gp
            ZCtot(cst_gas_index(i1))=ZC_cst(i1)
            ZCtot_save(cst_gas_index(i1))=ZCtot(cst_gas_index(i1))
            ratio_gas(i1)=1.d0
         enddo
      endif
      ZCtot_save(:)=ZCtot(:)    ! genoa keep_gp

C     Initialization of granulo for kinetic cte of heterogeneous rxns
      DO Jb=1,Nbin_aer
         conc_tot=0.0D0
         DO Jsp=1,Ns_aer-1
            conc_tot = conc_tot + DLconc_aer(Jb,Jsp) 
         ENDDO
C     compute the particle number (mass/geometric mean diameter)
         IF (INUM.EQ.1) THEN
            granulo_aer(Jb) = DLnumconc_aer(Jb)
         ELSE
            if (MSF(Jb) .GT. 0.D0) THEN
               granulo_aer(Jb) = conc_tot/MSF(Jb)               
            else
               write(*,*) "chem.f : error "
               stop
            endif
         ENDIF
      ENDDO

C     two-step solver starts
!tolerence
      rtol=TOL
      atol=1.d-6*TOL
      m=2                       !m nombre d'iteration
      tstep=tstep_min           ! initial step tstep: dtn
      tschem=ts                 ! initial timestep
      tfchem=tschem+delta_t     ! total chemistry step

! init conc.
      do Jsp=1,NS
         concii(Jsp)=ZCtot(Jsp) ! genoa keep_gp
         conci(Jsp)=ZCtot(Jsp)  ! genoa keep_gp
         concz(Jsp)=ZCwall(Jsp)
         conczz(Jsp)=ZCwall(Jsp)
      enddo

!premier calcul de l'ordre 1
      do j=1,m

!computes chpr and chlo which are the vectors with all the terms of production (chpr=P) and losses (chlo=LxC)

! compute RO2
         RO2=0.d0
         if (tag_RO2.ne.0) then
! add background if tag = 2,3,4 (from ro2 file)
            if (tag_RO2.ge.2
     &           .and.iRO2_cst.ne.0) then
               RO2 = RO2 + ZC_cst(iRO2_cst)
            endif
! add from primary VOCs (RO2 list) tag = 1,3
            if (tag_RO2 .eq. 1 .or. tag_RO2 .eq. 3 
     &           ) then
               do Jsp =1,nRO2_chem
                  RO2=RO2+ZC(RO2index(Jsp))
               enddo
            endif
         endif

! kinetic rate
! Compute zenithal angles
         DLmuzero=ssh_muzero(tschem,Dlon,Dlat)
         Zangzen=dabs(DACOS(DLmuzero)*180.D0/PI)
         DLmuzero=ssh_muzero(tfchem,Dlon,Dlat)
         Zangzenf=dabs(DACOS(DLmuzero)*180.D0/PI)
         CALL SSH_Kinetic(Ns,Nbin_aer,Nr,
     &        IHETER,ICLD,DLRKi,DLtemp,
     &        DLhumid,DLpress,Zangzen,
     &        Zatt,lwca,granulo_aer,
     &        wet_diameter_aer,DSF,
     &        hetero_species_index,Wmol,LWCmin,
     &        option_photolysis,RO2)

! keep inorganic constant
         if (ncst_gas.gt.0) then
            do i1=1,ncst_gas
               ZC(cst_gas_index(i1))=ZC_cst(i1)
               ZCtot(cst_gas_index(i1))=ZC_cst(i1) ! genoa keep_gp
            enddo
         endif

!C    !prod(Ns): chpr0
! W(Nr): reaction rates.: DLRkf
         call ssh_rates(Ns,Nr,DLRKi,ZC,DLRkf)
         call ssh_fexprod(Nr,Ns,DLRkf,chpr0)

!C    !loss(Ns): chlo0
!dw(Nr,Ns): derivative of reaction rates wrt Y.
         call ssh_dratedc(Ns,Nr,DLRKi,ZC,dw)
         call ssh_fexloss(Nr,Ns,dw,chlo0)

         IF (kwall_gas>0.d0) THEN 
            do Jsp=1,ns_aer
               if (aerosol_species_interact(Jsp)>0) then
                  if (aerosol_type(Jsp) == 4 .and. psat(Jsp)>0) then
                     Jsp2=aerosol_species_interact(Jsp)
                     chlo0(Jsp2)=chlo0(Jsp2)+kwall_gas*ratio_gas(Jsp2)
                     chpr0(Jsp2)=chpr0(Jsp2)+
     &                    kwall_gas/(kpart(Jsp)*Cwall)*ZCwall(Jsp2)
                     ZCwall(Jsp2)=(concz(Jsp2)+tstep*kwall_gas
     &                    *ZCtot(Jsp2)*ratio_gas(Jsp2))
     &                    /(dun+tstep*kwall_gas/(kpart(Jsp)*Cwall))
                     if(ZCwall(Jsp2)<dzero) ZCwall(Jsp2)=dzero
!     print*,ZCwall(Jsp2),kwall_gas,concz(Jsp2),ZCtot(Jsp2),
!     &                 kpart(Jsp),Cwall,kwall_gas/(kpart(Jsp)*Cwall)
                  endif
               endif
            enddo
         ENDIF

         do Jsp=1,NS
            

            
            if (chpr0(Jsp)>dzero.or.chlo0(Jsp)>dzero) then
!init
               ratloss=chlo0(Jsp)*ratio_gas(Jsp) ! genoa keep_gp
               ZCtot(Jsp)=(conci(Jsp)+tstep*chpr0(Jsp))
     &              /(dun+tstep*ratloss)
               ZC(Jsp)=ratio_gas(Jsp)*ZCtot(Jsp)
! clip
               if(ZC(Jsp)<dzero) ZC(Jsp)=dzero
               if(ZCtot(Jsp)<dzero) ZCtot(Jsp)=dzero
            endif
         enddo
      enddo

! keep inorganic constant
      if (ncst_gas.gt.0) then
         do i1=1,ncst_gas
            ZC(cst_gas_index(i1))=ZC_cst(i1)
            ZCtot(cst_gas_index(i1))=ZC_cst(i1) ! genoa keep_gp
         enddo
      endif

! calculate error
      error_max=0.d0 
      do Jsp=1,NS
         wk=atol+rtol*ZCtot(Jsp)*ratio_gas(Jsp) ! genoa keep_gp
         error_max=max(error_max,
     s        ratio_gas(Jsp)*abs(ZCtot(Jsp)-conci(Jsp))/wk)
      enddo

      IF (kwall_gas>0.d0) THEN
         do Jsp=1,NS         
            wk=atol+rtol*ZCwall(Jsp)
            error_max=max(error_max,
     s           abs(ZCwall(Jsp)-concz(Jsp))/wk)         
         enddo
      ENDIF

! set the current time
      tschem= tschem + tstep
! save timestep
      dtnsave=tstep

! assgin the conc. after 1st order calculation
      do Jsp=1,NS
         conci(Jsp)=ZCtot(Jsp)
         concz(Jsp)=ZCwall(Jsp)
      enddo

! pour evider diviser par zero
      if(error_max>dzero) then
         tstep=max(tstep_min,max(dun/alpha,min(alpha,0.8d0/
     &        (error_max**0.5d0)))*tstep)
      else
         tstep=alpha*tstep
      endif

      tstep=min(tstep,tfchem-tschem)
      if (tstep.gt.dzero) then
         c=dtnsave/tstep
      else
         c=1
      endif
      gam=(c+dun)/(c+2.d0)

!les calculs suivants de l'ordre 2
      do while (tschem<tfchem)
         do j=1,m
!computes chpr and chlo which are the vectors with all the terms of production (chpr=P) and losses (chlo=LxC)

! compute RO2 used in chem
            RO2=0.d0
            if (tag_RO2.ne.0) then
! add background tag = 2,3,4 (from ro2 file)
               if (tag_RO2.ge.2 
     &              .and.iRO2_cst.ne.0) then
                  RO2 = RO2 + ZC_cst(iRO2_cst)
               endif
! add from primary VOCs (RO2 list) tag = 1,3
               if (tag_RO2 .eq. 1 .or. tag_RO2 .eq. 3 
     &              ) then
                  do Jsp =1,nRO2_chem
                     RO2=RO2+ZC(RO2index(Jsp))
                  enddo
               endif
            endif

! kinetic rate
! Compute zenithal angles
            DLmuzero=ssh_muzero(tschem,Dlon,Dlat)
            Zangzen=dabs(DACOS(DLmuzero)*180.D0/PI)
            DLmuzero=ssh_muzero(tschem+tstep,Dlon,Dlat)
            Zangzenf=dabs(DACOS(DLmuzero)*180.D0/PI)
            CALL SSH_Kinetic(Ns,Nbin_aer,Nr,
     &           IHETER,ICLD,DLRKi,DLtemp,
     &           DLhumid,DLpress,Zangzen,
     &           Zatt,lwca,granulo_aer,
     &           wet_diameter_aer,DSF,
     &           hetero_species_index,Wmol,LWCmin,
     &           option_photolysis,RO2)

! keep inorganic constant
            if (ncst_gas.gt.0) then
               do i1=1,ncst_gas
                  ZC(cst_gas_index(i1))=ZC_cst(i1)
                  ZCtot(cst_gas_index(i1))=ZC_cst(i1) ! genoa keep_gp
               enddo
            endif

!C    !prod(Ns): chpr0
! W(Nr): reaction rates.: DLRkf
            call ssh_rates(Ns,Nr,DLRKi,ZC,DLRkf)
            call ssh_fexprod(Nr,Ns,DLRkf,chpr0)

!C    !loss(Ns): chlo0
!dw(Nr,Ns): derivative of reaction rates wrt Y.
            call ssh_dratedc(Ns,Nr,DLRKi,ZC,dw)
            call ssh_fexloss(Nr,Ns,dw,chlo0)

            IF (kwall_gas>0.d0) THEN
               do Jsp=1,ns_aer
                  if (aerosol_species_interact(Jsp)>0) then
                     if (aerosol_type(Jsp) == 4  .and. psat(Jsp)>0) then
                        Jsp2=aerosol_species_interact(Jsp)
                        chlo0(Jsp2)=chlo0(Jsp2)+kwall_gas
                        chpr0(Jsp2)=chpr0(Jsp2)+
     &                       kwall_gas/(kpart(Jsp)*Cwall)*ZCwall(Jsp2)
                        ZCwall(Jsp2)=(((c+dun)*(c+dun)*concz(Jsp2)-
     &                       conczz(Jsp2))/(c*c+2.d0*c)+
     &                       gam*tstep*kwall_gas*ZCtot(Jsp2)
     &                       * ratio_gas(Jsp2))/(dun+
     &                       gam*tstep*kwall_gas/(kpart(Jsp)*Cwall))
!print*,Jsp,tstep,kwall_gas/(kpart(Jsp)*Cwall)
                        if(ZCwall(Jsp2)<dzero) ZCwall(Jsp2)=dzero
                     endif
                  endif
               enddo
            ENDIF
            
            do Jsp=1,NS
               if (chpr0(Jsp)>dzero.or.chlo0(Jsp)>dzero) then
                  ratloss=chlo0(Jsp)*ratio_gas(Jsp)

! concc(Jsp)=((c+1)*(c+1)*conci(Jsp)-concii(Jsp))/(c*c+2*c)
! genoa keep_gp
                  ZCtot(Jsp)=(((c+dun)*(c+dun)*conci(Jsp)-
     &                 concii(Jsp))/(c*c+2.d0*c)+gam*tstep*
     &                 chpr0(Jsp))/(dun+gam*tstep*ratloss)
                  ZC(Jsp)=ratio_gas(Jsp)*ZCtot(Jsp)
! clip
                  if(ZC(Jsp)<dzero) ZC(Jsp)=dzero
                  if(ZCtot(Jsp)<dzero) ZCtot(Jsp)=dzero
               endif
            enddo
         enddo

! keep inorganic constant
         if (ncst_gas.gt.0) then
            do i1=1,ncst_gas
               ZC(cst_gas_index(i1))=ZC_cst(i1)
               ZCtot(cst_gas_index(i1))=ZC_cst(i1) ! genoa keep_gp
            enddo
         endif

         error_max=0.d0
         do Jsp=1,NS            ! genoa keep_gp
            wk=atol+rtol*abs(ZCtot(Jsp))*ratio_gas(Jsp)
            error_max=max(error_max,
     &           abs(2.0d0*(c*ratio_gas(Jsp)*ZCtot(Jsp)-
     &           (dun+c)*ratio_gas(Jsp)*conci(Jsp)+
     &           ratio_gas(Jsp)*concii(Jsp))/
     &           (c*(c+dun)*wk)))
         enddo

         IF (kwall_gas>0.d0) THEN
            do Jsp=1,NS
               wk=atol+rtol*ZCwall(Jsp)
               error_max=max(error_max,
     &              abs(2.0d0*(c*ZCwall(Jsp)-
     &              (dun+c)*concz(Jsp)+
     &              conczz(Jsp))/
     &              (c*(c+dun)*wk)))            
            enddo
         endif

         do while (error_max>10.0d0 .and. tstep>tstep_min)
            tstep=max(tstep_min,max(dun/alpha,min(alpha,
     &           0.8d0/(error_max**0.5d0)))*tstep)
            tstep=min(tstep,tfchem-tschem)
            if (tstep.gt.dzero) then
               c=dtnsave/tstep
            else
               c=1
            endif
            gam=(c+1)/(c+2.d0)
            do Jsp=1,NS
               ZCtot(Jsp)=conci(Jsp) !FCo ! genoa keep_gp
               ZC(Jsp)=ratio_gas(Jsp)*ZCtot(Jsp)
               ZCwall(Jsp)=concz(Jsp)
            enddo

            do j=1,m
!computes chpr and chlo which are the vectors with all the terms of production (chpr=P) and losses (chlo=LxC)

! compute RO2 used in chem
               RO2=0.d0
               if (tag_RO2.ne.0) then
! add background tag = 2,3,4 (from ro2 file)
                  if (tag_RO2.ge.2.
     &                 and.iRO2_cst.ne.0) then
                     RO2 = RO2 + ZC_cst(iRO2_cst)
                  endif
! add from primary VOCs (RO2 list) tag = 1,3
                  if (tag_RO2 .eq. 1 .or. tag_RO2 .eq. 3 
     &                 ) then
                     do Jsp =1,nRO2_chem
                        RO2=RO2+ZC(RO2index(Jsp))
                     enddo
                  endif
               endif

! kinetic rate
! Compute zenithal angles
               DLmuzero=ssh_muzero(tschem,Dlon,Dlat)
               Zangzen=dabs(DACOS(DLmuzero)*180.D0/PI)
               DLmuzero=ssh_muzero(tschem+tstep,Dlon,Dlat)
               Zangzenf=dabs(DACOS(DLmuzero)*180.D0/PI)
               CALL SSH_Kinetic(Ns,Nbin_aer,Nr,
     &              IHETER,ICLD,DLRKi,DLtemp,
     &              DLhumid,DLpress,Zangzen,
     &              Zatt,lwca,granulo_aer,
     &              wet_diameter_aer,DSF,
     &              hetero_species_index,Wmol,LWCmin,
     &              option_photolysis,RO2)

! keep inorganic constant
               if (ncst_gas.gt.0) then
                  do i1=1,ncst_gas
                     ZC(cst_gas_index(i1))=ZC_cst(i1)
                     ZCtot(cst_gas_index(i1))=ZC_cst(i1) ! genoa keep_gp
                  enddo
               endif

!C    !prod(Ns): chpr0
! W(Nr): reaction rates.: DLRkf
               call ssh_rates(Ns,Nr,DLRKi,ZC,DLRkf)
               call ssh_fexprod(Nr,Ns,DLRkf,chpr0)

!C    !loss(Ns): chlo0
!dw(Nr,Ns): derivative of reaction rates wrt Y.
               call ssh_dratedc(Ns,Nr,DLRKi,ZC,dw)
               call ssh_fexloss(Nr,Ns,dw,chlo0)

               IF (kwall_gas>0.d0) THEN
                  do Jsp=1,ns_aer
                     if (aerosol_species_interact(Jsp)>0) then
                        if (aerosol_type(Jsp)==4.and.psat(Jsp)>0) then
                           Jsp2=aerosol_species_interact(Jsp)
                           chlo0(Jsp2)=chlo0(Jsp2)+kwall_gas
                           chpr0(Jsp2)=chpr0(Jsp2)+
     &                          kwall_gas/(kpart(Jsp)*Cwall)
     &                          *ZCwall(Jsp2)
                           ZCwall(Jsp2)=(((c+dun)*(c+dun)*concz(Jsp2)-
     &                          conczz(Jsp2))/(c*c+2.d0*c)+
     &                          gam*tstep*kwall_gas*ZCtot(Jsp2)
     &                          * ratio_gas(Jsp2))
     &                          /(dun+gam*tstep*kwall_gas
     &                          /(kpart(Jsp)*Cwall))
                           if(ZCwall(Jsp2)<dzero) ZCwall(Jsp2)=dzero
                        endif
                     endif
                  enddo
               endif
            
               do Jsp=1,NS      ! genoa keep_gp
                  if (chpr0(Jsp)>dzero.or.chlo0(Jsp)>dzero) then
                     ratloss=chlo0(Jsp)*ratio_gas(Jsp)
                     
!     concc(Jsp)=((c+1)*(c+1)*conci(Jsp)-concii(Jsp))/(c*c+2*c)
                     ZCtot(Jsp)=(((c+dun)*(c+dun)*conci(Jsp)-
     &                    concii(Jsp))/(c*c+2.d0*c)+gam*tstep*
     &                    chpr0(Jsp))/(dun+gam*tstep*ratloss)
                     ZC(Jsp)=ratio_gas(Jsp)*ZCtot(Jsp)
!     clip
                     if(ZCtot(Jsp)<dzero) ZCtot(Jsp)=dzero
                     if(ZC(Jsp)<dzero) ZC(Jsp)=dzero
                  endif
               enddo
            enddo
            
!     keep inorganic constant
            if (ncst_gas.gt.0) then
               do i1=1,ncst_gas
                  ZC(cst_gas_index(i1))=ZC_cst(i1)
                  ZCtot(cst_gas_index(i1))=ZC_cst(i1) ! genoa keep_gp
               enddo
            endif
            
            error_max=0.d0
            do Jsp=1,NS         ! genoa keep_gp
               wk=atol+rtol*abs(ZCtot(Jsp))*ratio_gas(Jsp)
               error_max=max(error_max,
     &              abs(2.0d0*(c*ZCtot(Jsp)*ratio_gas(Jsp)-
     &              (dun+c)*ratio_gas(Jsp)*conci(Jsp)
     &              +ratio_gas(Jsp)*concii(Jsp))/
     &              (c*(c+dun)*wk)))
            enddo
            IF (kwall_gas>0.d0) THEN
               do Jsp=1,NS                     
                  wk=atol+rtol*ZCwall(Jsp)
                  error_max=max(error_max,
     &                 abs(2.0d0*(c*ZCwall(Jsp)-
     &                 (dun+c)*concz(Jsp)+
     &                 conczz(Jsp))/
     &                 (c*(c+dun)*wk)))
               enddo
            endif
         enddo
            
!save timestep
         dtnsave=tstep
! update current time
         tschem=tschem+tstep
! update timestep
         if(error_max>dzero) then
            tstep=max(tstep_min,max(dun/alpha,min(alpha,
     &           0.8d0/(error_max**0.5d0)))*tstep)
         else
            tstep=alpha*tstep
         endif

         tstep=min(tstep,tfchem-tschem)

         if (tstep.gt.dzero) then
            c=dtnsave/tstep
         else
            c=1
         endif
         gam=(c+dun)/(c+2.d0)

         do Jsp=1,NS
            conczz(Jsp)=concz(Jsp)
            concz(Jsp)=ZCwall(Jsp)
            concii(Jsp)=conci(Jsp)
            conci(Jsp)=ZCtot(Jsp) ! genoa keep_gp
! keep inorganic constant
            if (ncst_gas.gt.0) then
               do i1=1,ncst_gas
                  concii(cst_gas_index(i1))=ZC_cst(i1)
                  conci(cst_gas_index(i1))=ZC_cst(i1)
                  ZC(cst_gas_index(i1))=ZC_cst(i1)
                  ZCtot(cst_gas_index(i1))=ZC_cst(i1)
               enddo
            endif
         enddo
!print*,'      timestep',tstep,tschem
      enddo

C     two-step solver end

C     Storage in the array of chemical concentrations.

      DO Jsp=1,NS
! NAN detection algorithm used
         if (ZC(Jsp).ne.ZC(Jsp)) then
            write(*,*) "From chem 0D function:"
            write(*,*) Jsp, dlconc(Jsp),convers_factor(Jsp),
     &           ZC(Jsp)/convers_factor(Jsp)
            stop
         endif
! Conversion molecules/cm3 to mug/m3.
         DLconc(Jsp) = ZC(Jsp)/convers_factor(Jsp)
         DLwall(Jsp) = ZCwall(Jsp)/convers_factor(Jsp)
      ENDDO

      if (keep_gp==1) then
         do s = 1, Ns_aer
            Jsp=aerosol_species_interact(s)
            if (Jsp>0) then
!if(ratio_gas(Jsp)<1.d0) then
               conc_tot=0.0D0
               DO Jb=1,Nbin_aer
                  conc_tot=conc_tot+DLconc_aer(Jb,s) !*convers_factor(Jsp) 
               ENDDO
               
               toadd=(ZCtot(Jsp)-ZCtot_save(Jsp))/convers_factor(Jsp)
               if (DLconc_save(Jsp)+toadd>=0.d0) then
                  DLconc(Jsp)=DLconc_save(Jsp)+toadd
               else             !Not enough mass in the gas phase, have to be taken from the particle phase
!     print*,toadd,DLconc_save(Jsp)+conc_tot,
!     s                    ZCtot(Jsp)/convers_factor(Jsp),
!     s                    ZCtot_save(Jsp)/convers_factor(Jsp)
                  toadd=(conc_tot+toadd+DLconc_save(Jsp))/conc_tot
                  DLconc(Jsp)=0.d0
                  DO Jb=1,Nbin_aer
                     DLconc_aer(Jb,s) = DLconc_aer(Jb,s)*toadd 
                  ENDDO
               endif

               conc_tot=0.d0
               DO Jb=1,Nbin_aer
                  conc_tot=conc_tot+DLconc_aer(Jb,s) 
               ENDDO
!     print*,conc_tot+DLconc(Jsp)
!     s                 ,ZCtot(Jsp)/convers_factor(Jsp)
!endif
            endif
         ENDDO
      endif

! keep inorganic constant
      if (ncst_gas.gt.0) then
         do i1=1,ncst_gas
            DLconc(cst_gas_index(i1))=ZC(cst_gas_index(i1))/
     &           convers_factor(cst_gas_index(i1))
         enddo
      endif

! output RO2
      if (iRO2.ne.0.and.tag_RO2.ne.0) then
! output only ro2 from list
         DLconc(iRO2)=0.d0
         do Jsp =1,nRO2_chem
            DLconc(iRO2)= DLconc(iRO2)+
     &           ZC(RO2index(Jsp))/
     &           convers_factor(RO2index(Jsp))
         enddo
      endif

      IF (kwall_particle>0.d0) THEN
         wloss=exp(-kwall_particle*delta_t)
         DO Jb = 1, nbin_aer
            DLnumconc_aer(Jb)=DLnumconc_aer(Jb)*wloss
            DO Jsp = 1, Ns_aer
               DLconc_aer(Jb,Jsp)=DLconc_aer(Jb,Jsp)*wloss
            ENDDO
         ENDDO
      ENDIF

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
