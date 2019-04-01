C-----------------------------------------------------------------------
C     Copyright (C) 2003-2007, ENPC - INRIA - EDF R&D
C     Author(s): Kathleen Fahey
C
C     This file is part of the Regional Atmospheric Chemistry Modeling
C     (RACM), which is a component of the air quality modeling system
C     Polyphemus.
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

      SUBROUTINE HETRXN(Ns,Nbin_aer,temp,press,ICLD,lwctmp,
     &     WetDiam,granulo,rkin1,rkin2,rkin3,rkin4,
     &     dsf_aero,ispeclost,Wmol,LWCmin)

C------------------------------------------------------------------------
C
C     -- DESCRIPTION
C
C     This subroutine calculates the kinetic constants
C     of the heterogeneous reactions that happen on the surface of particles, as
C     described in Jacob (2000).  It is considered an irreversible
C     process with a first-order reaction rate.
C
C     HO2  -->  0.5 H2O2
C     NO2  -->  0.5 HONO + 0.5 HNO3
C     NO3  -->  HNO3
C     N2O5 -->  2 HNO3
C
C     REFERENCES:
C     Jacob, D.J. (2000) Heterogeneous chemistry and tropospheric
C     ozone.  Atmospheric Environment.  34, 2131-2159.
C
C     Jacob, D.J. (2004) Private communication to B. Sportisse.
C
C
C------------------------------------------------------------------------
C
C     -- INPUT VARIABLES
C     temp: Temerature ([K])
C     press: Pressure ([Pa])
C     LCL: logical for aqueous chemistry (from user choices).
C     lwctmp: Liquid Water content ([]g.m-3)
C     WetDiam: Aerosol size bins wet diameters ([\mu m])
C     granulo: Particles number in each bin.
C
C     -- INPUT/OUTPUT VARIABLES
C
C
C     -- OUTPUT VARIABLES
C     rkin1, rkin2, rkin3 and rkin4: kinetic constants of the 4 reactions
C     (in that order)
C
C------------------------------------------------------------------------
C
C     -- REMARKS
C
C------------------------------------------------------------------------
C
C     -- MODIFICATIONS
c
C     2005-05-23: kinetic constants in order to insert heterogeneous
C     reaction in the gaseous sheme (instead of in the aerosol module)
C     2005-11-28: add reactions for HO2 and NO3 on cloud droplet surfaces.
C     (Marilyne Tombette, CEREA)
C
C------------------------------------------------------------------------
C
C     -- AUTHOR(S)
C
C     Kathleen Fahey, CEREA
C
C------------------------------------------------------------------------



      IMPLICIT NONE

      INCLUDE 'CONST.INC'
      INCLUDE 'CONST_A.INC'
      INCLUDE 'hetrxn.inc'
      INCLUDE 'aerpar.inc'
      INCLUDE 'droppar.inc'

      INTEGER Ns,Nbin_aer


      INTEGER i,js
      INTEGER ispeclost(4)
      DOUBLE PRECISION Ndroplet

      DOUBLE PRECISION temp,press
      DOUBLE PRECISION rkin1,rkin2,rkin3,rkin4
      DOUBLE PRECISION vi
      DOUBLE PRECISION ktot(4),A
      DOUBLE PRECISION DIFF
      DOUBLE PRECISION avdiammeter

      DOUBLE PRECISION WetDiam(Nbin_aer)
      DOUBLE PRECISION WetDiammeter(Nbin_aer)
      DOUBLE PRECISION granulo(Nbin_aer)
      DOUBLE PRECISION dsf_aero(Nbin_aer)

      DOUBLE PRECISION lwctmp

      INTEGER ICLOUD
      INTEGER ICLD

      DOUBLE PRECISION sigm_tmp,parm_tmp,wmol_tmp
      DOUBLE PRECISION Wmol(Ns),LWCmin

C     Constants.
      avdiammeter = 1.d-6 * avdiam

      sigm_tmp = SIGM_NO2
      parm_tmp = PARM_NO2
      ICLOUD=0
      IF (ICLD.GE.1.AND.(lwctmp.GE.LWCmin)) THEN
         ICLOUD=1
      ENDIF

C     REACTION PROBABILITIES

      DO i=1,4
         ktot(i)=0.D0
         wmol_tmp = Wmol(ispeclost(i)+1)
         call COMPUTE_GAS_DIFFUSIVITY(TEMP,PRESS,
     &        sigm_tmp, wmol_tmp, parm_tmp, DIFF)

         call COMPUTE_QUADRATIC_MEAN_VELOCITY(TEMP,
     &        wmol_tmp,vi)

         do js=1,Nbin_aer
            IF ((ICLOUD.NE.1).OR.(dsf_aero(js).lt.dactiv*1.d6)) then
               WetDiammeter(js) = 1.d-6 * WetDiam(js)

C     CALCULATE k FOR EACH SECTION

               A = PI*(WetDiammeter(js)**2.d0)*granulo(js) !surface * nb_aero

C     CALCULATE k FOR AEROSOL DISTRIBUTION

               IF(Gamma(i).ne.0.d0) then
                  ktot(i) = ktot(i)
     &                      + (((WetDiammeter(js) / (2.d0 * DIFF))
     &                         + 4.d0 / (vi * Gamma(i)))**(-1.d0)) * A ! rate constant for each rxn
               ENDIF
            ENDIF
         enddo

C     If we are in a cloud, then reactions on droplet surface for N2O5.
         IF ((ICLOUD.EQ.1).and.(i.eq.4)) THEN
            Ndroplet = lwctmp
     &           / (RHOwater * 1.d3 !kg.m-3 -> g.m-3
     &           * pi / 6.d0 * avdiammeter**3.d0)
            A = PI * avdiammeter**2.d0 * Ndroplet !surface * nb_droplet
            IF(Gamma(i).ne.0.d0) then
               ktot(i) = ktot(i) + ((avdiammeter / (2.d0 * DIFF) +
     &           4.d0 / (vi*Gamma(i)))**(-1.d0))*A ! rate constant for each rxn
            ENDIF
         ENDIF
      enddo

      rkin1=ktot(1)
      rkin2=ktot(2)
      rkin3=ktot(3)
      rkin4=ktot(4)

      END
