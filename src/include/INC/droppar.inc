C-----------------------------------------------------------------------
C     Copyright (C) 2003-2007, ENPC - INRIA - EDF R&D
C     Author(s): Kathleen Fahey
C
C     This file is part of the Variable Size Resolved Model (VSRM),
C     based on the VSRM model of Carnegie Melon University. It is a
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

c     UPDATED DROPPAR.INC FOR THE SIZE RESOLVED MODEL
c     LAST UPDATE : 22 November 2005
c*************************************************************************
c     DROPPAR.INC
c*************************************************************************
c     
c     AQUEOUS-PHASE PARAMETERS AND VARIABLES

      INTEGER ksod,khyd,kamm,knit
      INTEGER kchl,ksvi,kwat,kec1,koc
      INTEGER kcru              !,khso5,khmsa,kform
      INTEGER ngso2,ngh2o2,nghcho
      INTEGER nghcooh,nghno2,ngno,ngno2
      INTEGER ngo3,ngpan,ngoh,ngho2,ngno3
      INTEGER ngch3o2,ngch3o2h,ngch3oh
      INTEGER ngch3co3h,meqn1,meqn2,kiron
      INTEGER iradical,ires
c     
c     NSP : NUMBER OF TOTAL PARTICLE SECTIONS
c     
c     
c     NBOUNDS : NUMBER OF SECTIONAL BOUNDARIES
c     
c     AQUEOUS-PHASE COMPONENTS
c     
c     IMPORTANT : ALL COMPONENTS HAVE THE SAME POSITIONS IN      
c     BOTH AEROSOL AND AQUEOUS MATRICES
c     NEVER CHANGE THIS CONVENTION BECAUSE AQMAIN
c     DEPENDS ON IT
c     
      parameter (ksod = nas)    ! Na(+)
      parameter (khyd = nah)    ! H(+)
      parameter (kamm = naa)    ! NH4(+)
      parameter (knit = nan)    ! NO3(-)
      parameter (kchl = nac)    ! Cl(-)
      parameter (ksvi = na4)    ! S(VI)
      parameter (kwat = naw)    ! H2O
      parameter (kec1 = nae)    ! EC
      parameter (koc  = nao)    ! OC
      parameter (kcru = nar)    ! Crustal
c     parameter (khso5 = 1)                ! HSO5-
c     parameter (khmsa = 2)                ! HMSA
c     parameter (kform = 3)                ! Formic Acid
c     
c     GASES IN LOCAL ARRAY
c     
c     INCORRECT  ******* TO BE FIXED *************      
      parameter (ngso2 = 11)
      parameter (ngh2o2 = 12)
      parameter (nghcho = 13)
      parameter (nghcooh = 14)
      parameter (nghno2 = 15)
      parameter (ngno = 16)
      parameter (ngno2 = 17)
      parameter (ngo3 = 18)
      parameter (ngpan = 19)
      parameter (ngoh = 20)
      parameter (ngho2 = 21)
      parameter (ngno3 = 22)
      parameter (ngch3o2 = 23)
      parameter (ngch3o2h = 24)
      parameter (ngch3oh = 25)
      parameter (ngch3co3h = 26)
c     
c     NUMBER OF EQUATIONS FOR AQUEOUS-PHASE CHEMISTRY SOLUTION
c     
      parameter (meqn1 = 13)
      parameter (meqn2 = 21)
c     
c     ACTIVATION DIAMETER (Dry)
c     
      double precision dactiv,dsep,avdiam
      parameter (dactiv = 0.2D0) ! in \mu m

c     
c     SEPARATION DIAMETER FOR TWO SECTIONS (DRY)
c     
      parameter (dsep = 1.5D0)  ! in \mu m
c     
c     WET DIAMETERS OF THE TWO SECTIONS

c     If bulk, diameter = avdiam
      parameter (avdiam = 20)   ! in \mu m     
c     
c     CHOICE OF EXPRESSION FOR IRON CHEMISTRY
c     = 0 (No iron/manganese chemistry)
c     KIRON = 1 (Phenomenological, Martin et al., 1991)
c     = 2 (Martin, 1984)
c     
      parameter (kiron = 2)     ! was 1
c     
c     CHOICE OF TURNING ON OR OFF RADICAL CHEMISTY
c     (IT IS BETTER TO TURN IT OFF DURING THE NIGHT)
c     
      parameter (iradical = 0)  ! kmf (07/24/02)
c     
c     CHOICE OF TURNING OFF CHLORINE CHEMISTRY
c     
      double precision chlorine
      parameter (chlorine = 0.d0)
c     
c     CHOICE OF TURNING OFF NO2 RXN (104) ! kmf (08/27/02)
c     
      double precision nitrogen
      parameter (nitrogen = 0.d0) ! uncertain rxn rate - producing too much SO4 here
c     
c     PARAMETER FOR SCALING OF PHOTOLYSIS RATES
c     
      double precision photo
      parameter (photo = 1.d0)  
c     
c     Fraction of crustal material that is alkaline 
c     
      double precision caratio  ! kmf (08/01/02)
      parameter (caratio = 0.d0) ! was 0.1  
c     
c     
c     Size resolution 1=bulk,
c     2=2-section sr
c     3=variable sr)
c     bulk treats all droplets as homogeneous
c     2-section size resolved divides droplets into 2 separate 
c     homogeneous sections
c     variable size-resolved - executes either the bulk or size-resolved
c     calculations based on initial conditions (NH3, HNO3, SO2, H2O2, 
c     lwc, and caratio)
c     
      parameter (ires = 1)
c     
c     fraction of liquid water content that goes to each s.r. section
c     
      double precision frac1,frac2
      parameter (frac1 = 0.8d0) ! fraction of lwc in sect. 1
      parameter (frac2 = 0.2d0) ! fraction of lwc in sect. 2
c     
c     
c     ASSUMPTION : Fe(3+) and Mn(2+) = 0.003%, 0.001% of Crustal mass
c     
      double precision firon,fman
      parameter (firon = 0.00003d0)
      parameter (fman = 0.00001d0)
C     
      double precision minso2
      parameter (minso2 = 1.5d-4)

C     LWCFOG: minimum threshold for fog
      
      double precision LWCFOG
      parameter (LWCFOG=0.05d0) ! in g.m-3 

C     Parameters for bimodal activated distribution: 
C     DM: mean diameter; SD: variance.

      double precision DM1_aq, DM2_aq, SD1_aq, SD2_aq
      
      parameter(DM1_aq = 0.4d0)
      parameter(DM2_aq = 2.5d0)
      parameter(SD1_aq = 1.8d0)
      parameter(SD2_aq = 2.15d0)

      double precision DG_AQ, ACC_AQ
      parameter(DG_AQ=1.d-5)
      parameter(ACC_AQ=0.1d0)
