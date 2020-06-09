C=======================================================================
C
C *** ISORROPIA CODE
C *** SUBROUTINE SSH_ISRP1F
C *** THIS SUBROUTINE SSH_IS THE DRIVER ROUTINE FOR THE FOREWARD PROBLEM OF 
C     AN AMMONIUM-SULFATE AEROSOL SYSTEM. 
C     THE COMPOSITION REGIME IS DETERMINED BY THE SULFATE RATIO AND BY 
C     THE AMBIENT RELATIVE HUMIDITY.
C
C *** COPYRIGHT 1996-2006, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY,
C *** GEORGIA INSTITUTE OF TECHNOLOGY
C *** WRITTEN BY ATHANASIOS NENES
C *** UPDATED BY CHRISTOS FOUNTOUKIS
C
C=======================================================================
C
      SUBROUTINE SSH_ISRP1F (WI, RHI, TEMPI)
      INCLUDE 'isrpia.inc'
      DIMENSION WI(NCOMP)
C
C *** INITIALIZE ALL VARIABLES IN COMMON BLOCK **************************
C
      CALL SSH_INIT1 (WI, RHI, TEMPI)
C
C *** CALCULATE SULFATE RATIO *******************************************
C
      SULRAT = W(3)/W(2)
C
C *** FIND CALCULATION REGIME FROM (SULRAT,RH) **************************
C
C *** SULFATE POOR 
C
      IF (2.0.LE.SULRAT) THEN 
      DC   = W(3) - 2.001D0*W(2)  ! For numerical stability
      W(3) = W(3) + MAX(-DC, ZERO)
C
      IF(METSTBL.EQ.1) THEN
         SCASE = 'A2'
         CALL SSH_CALCA2                 ! Only liquid (metastable)
      ELSE
C
         IF (RH.LT.DRNH42S4) THEN    
            SCASE = 'A1'
            CALL SSH_CALCA1              ! NH42SO4              ; case A1
C
         ELSEIF (DRNH42S4.LE.RH) THEN
            SCASE = 'A2'
            CALL SSH_CALCA2              ! Only liquid          ; case A2
         ENDIF
      ENDIF
C
C *** SULFATE RICH (NO ACID)
C
      ELSEIF (1.0.LE.SULRAT .AND. SULRAT.LT.2.0) THEN 
C
      IF(METSTBL.EQ.1) THEN
         SCASE = 'B4'
         CALL SSH_CALCB4                 ! Only liquid (metastable)
      ELSE
C
         IF (RH.LT.DRNH4HS4) THEN         
            SCASE = 'B1'
            CALL SSH_CALCB1              ! NH4HSO4,LC,NH42SO4   ; case B1
C
         ELSEIF (DRNH4HS4.LE.RH .AND. RH.LT.DRLC) THEN         
            SCASE = 'B2'
            CALL SSH_CALCB2              ! LC,NH42S4            ; case B2
C
         ELSEIF (DRLC.LE.RH .AND. RH.LT.DRNH42S4) THEN         
            SCASE = 'B3'
            CALL SSH_CALCB3              ! NH42S4               ; case B3
C
         ELSEIF (DRNH42S4.LE.RH) THEN         
            SCASE = 'B4'
            CALL SSH_CALCB4              ! Only liquid          ; case B4
         ENDIF
      ENDIF
      CALL SSH_CALCNH3
C
C *** SULFATE RICH (FREE ACID)
C
      ELSEIF (SULRAT.LT.1.0) THEN             
C
      IF(METSTBL.EQ.1) THEN
         SCASE = 'C2'
         CALL SSH_CALCC2                 ! Only liquid (metastable)
      ELSE
C
         IF (RH.LT.DRNH4HS4) THEN         
            SCASE = 'C1'
            CALL SSH_CALCC1              ! NH4HSO4              ; case C1
C
         ELSEIF (DRNH4HS4.LE.RH) THEN         
            SCASE = 'C2'
            CALL SSH_CALCC2              ! Only liquid          ; case C2
C
         ENDIF
      ENDIF
      CALL SSH_CALCNH3
      ENDIF
C
C *** RETURN POINT
C
      RETURN
C
C *** END OF SUBROUTINE SSH_ISRP1F *****************************************
C
      END
C=======================================================================
C
C *** ISORROPIA CODE
C *** SUBROUTINE SSH_ISRP2F
C *** THIS SUBROUTINE SSH_IS THE DRIVER ROUTINE FOR THE FOREWARD PROBLEM OF 
C     AN AMMONIUM-SULFATE-NITRATE AEROSOL SYSTEM. 
C     THE COMPOSITION REGIME IS DETERMINED BY THE SULFATE RATIO AND BY
C     THE AMBIENT RELATIVE HUMIDITY.
C
C *** COPYRIGHT 1996-2006, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY,
C *** GEORGIA INSTITUTE OF TECHNOLOGY
C *** WRITTEN BY ATHANASIOS NENES
C *** UPDATED BY CHRISTOS FOUNTOUKIS
C
C=======================================================================
C
      SUBROUTINE SSH_ISRP2F (WI, RHI, TEMPI)
      INCLUDE 'isrpia.inc'
      DIMENSION WI(NCOMP)
C
C *** INITIALIZE ALL VARIABLES IN COMMON BLOCK **************************
C
      CALL SSH_INIT2 (WI, RHI, TEMPI)
C
C *** CALCULATE SULFATE RATIO *******************************************
C
      SULRAT = W(3)/W(2)
C
C *** FIND CALCULATION REGIME FROM (SULRAT,RH) **************************
C
C *** SULFATE POOR 
C
      IF (2.0.LE.SULRAT) THEN                
C
      IF(METSTBL.EQ.1) THEN
         SCASE = 'D3'
         CALL SSH_CALCD3                 ! Only liquid (metastable)
      ELSE
C
         IF (RH.LT.DRNH4NO3) THEN    
            SCASE = 'D1'
            CALL SSH_CALCD1              ! NH42SO4,NH4NO3       ; case D1
C
         ELSEIF (DRNH4NO3.LE.RH .AND. RH.LT.DRNH42S4) THEN         
            SCASE = 'D2'
            CALL SSH_CALCD2              ! NH42S4               ; case D2
C
         ELSEIF (DRNH42S4.LE.RH) THEN
            SCASE = 'D3'
            CALL SSH_CALCD3              ! Only liquid          ; case D3
         ENDIF
      ENDIF
C
C *** SULFATE RICH (NO ACID)
C     FOR SOLVING THIS CASE, NITRIC ACID IS ASSUMED A MINOR SPECIES, 
C     THAT DOES NOT SIGNIFICANTLY PERTURB THE HSO4-SO4 EQUILIBRIUM.
C     SUBROUTINES CALCB? ARE CALLED, AND THEN THE NITRIC ACID IS DISSOLVED
C     FROM THE HNO3(G) -> (H+) + (NO3-) EQUILIBRIUM.
C
      ELSEIF (1.0.LE.SULRAT .AND. SULRAT.LT.2.0) THEN 
C
      IF(METSTBL.EQ.1) THEN
         SCASE = 'B4'
         CALL SSH_CALCB4                 ! Only liquid (metastable)
         SCASE = 'E4'
      ELSE
C
         IF (RH.LT.DRNH4HS4) THEN         
            SCASE = 'B1'
            CALL SSH_CALCB1              ! NH4HSO4,LC,NH42SO4   ; case E1
            SCASE = 'E1'
C
         ELSEIF (DRNH4HS4.LE.RH .AND. RH.LT.DRLC) THEN         
            SCASE = 'B2'
            CALL SSH_CALCB2              ! LC,NH42S4            ; case E2
            SCASE = 'E2'
C
         ELSEIF (DRLC.LE.RH .AND. RH.LT.DRNH42S4) THEN         
            SCASE = 'B3'
            CALL SSH_CALCB3              ! NH42S4               ; case E3
            SCASE = 'E3'
C
         ELSEIF (DRNH42S4.LE.RH) THEN         
            SCASE = 'B4'
            CALL SSH_CALCB4              ! Only liquid          ; case E4
            SCASE = 'E4'
         ENDIF
      ENDIF
C
      CALL SSH_CALCNA                 ! HNO3(g) DISSOLUTION
C
C *** SULFATE RICH (FREE ACID)
C     FOR SOLVING THIS CASE, NITRIC ACID IS ASSUMED A MINOR SPECIES, 
C     THAT DOES NOT SIGNIFICANTLY PERTURB THE HSO4-SO4 EQUILIBRIUM
C     SUBROUTINE SSH_CALCC? IS CALLED, AND THEN THE NITRIC ACID IS DISSOLVED
C     FROM THE HNO3(G) -> (H+) + (NO3-) EQUILIBRIUM.
C
      ELSEIF (SULRAT.LT.1.0) THEN             
C
      IF(METSTBL.EQ.1) THEN
         SCASE = 'C2'
         CALL SSH_CALCC2                 ! Only liquid (metastable)
         SCASE = 'F2'
      ELSE
C
         IF (RH.LT.DRNH4HS4) THEN         
            SCASE = 'C1'
            CALL SSH_CALCC1              ! NH4HSO4              ; case F1
            SCASE = 'F1'
C
         ELSEIF (DRNH4HS4.LE.RH) THEN         
            SCASE = 'C2'
            CALL SSH_CALCC2              ! Only liquid          ; case F2
            SCASE = 'F2'
         ENDIF
      ENDIF
C
      CALL SSH_CALCNA                 ! HNO3(g) DISSOLUTION
      ENDIF
C
C *** RETURN POINT
C
      RETURN
C
C *** END OF SUBROUTINE SSH_ISRP2F *****************************************
C
      END
C=======================================================================
C
C *** ISORROPIA CODE
C *** SUBROUTINE SSH_ISRP3F
C *** THIS SUBROUTINE SSH_IS THE DRIVER ROUTINE FOR THE FORWARD PROBLEM OF
C     AN AMMONIUM-SULFATE-NITRATE-CHLORIDE-SODIUM AEROSOL SYSTEM. 
C     THE COMPOSITION REGIME IS DETERMINED BY THE SULFATE & SODIUM 
C     RATIOS AND BY THE AMBIENT RELATIVE HUMIDITY.
C
C *** COPYRIGHT 1996-2006, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY,
C *** GEORGIA INSTITUTE OF TECHNOLOGY
C *** WRITTEN BY ATHANASIOS NENES
C *** UPDATED BY CHRISTOS FOUNTOUKIS
C
C=======================================================================
C
      SUBROUTINE SSH_ISRP3F (WI, RHI, TEMPI)
      INCLUDE 'isrpia.inc'
      DIMENSION WI(NCOMP)
C
C *** ADJUST FOR TOO LITTLE AMMONIUM AND CHLORIDE ***********************
C
      WI(3) = MAX (WI(3), 1.D-10)  ! NH4+ : 1e-4 umoles/m3
      WI(5) = MAX (WI(5), 1.D-10)  ! Cl-  : 1e-4 umoles/m3
C
C *** ADJUST FOR TOO LITTLE SODIUM, SULFATE AND NITRATE COMBINED ********
C
      IF (WI(1)+WI(2)+WI(4) .LE. 1d-10) THEN
         WI(1) = 1.D-10  ! Na+  : 1e-4 umoles/m3
         WI(2) = 1.D-10  ! SO4- : 1e-4 umoles/m3
      ENDIF
C
C *** INITIALIZE ALL VARIABLES IN COMMON BLOCK **************************
C
      CALL SSH_ISOINIT3 (WI, RHI, TEMPI)
C
C *** CHECK IF TOO MUCH SODIUM ; ADJUST AND ISSUE ERROR MESSAGE *********
C
      REST = 2.D0*W(2) + W(4) + W(5) 
      IF (W(1).GT.REST) THEN            ! NA > 2*SO4+CL+NO3 ?
         W(1) = (ONE-1D-6)*REST         ! Adjust Na amount
         CALL SSH_PUSHERR (0050, 'ISRP3F')  ! Warning error: Na adjusted
      ENDIF
C
C *** CALCULATE SULFATE & SODIUM RATIOS *********************************
C
      SULRAT = (W(1)+W(3))/W(2)
      SODRAT = W(1)/W(2)
C
C *** FIND CALCULATION REGIME FROM (SULRAT,RH) **************************
C
C *** SULFATE POOR ; SODIUM POOR
C
      IF (2.0.LE.SULRAT .AND. SODRAT.LT.2.0) THEN                
C
      IF(METSTBL.EQ.1) THEN
         SCASE = 'G5'
         CALL SSH_CALCG5                 ! Only liquid (metastable)
      ELSE
C
         IF (RH.LT.DRNH4NO3) THEN    
            SCASE = 'G1'
            CALL SSH_CALCG1              ! NH42SO4,NH4NO3,NH4CL,NA2SO4
C
         ELSEIF (DRNH4NO3.LE.RH .AND. RH.LT.DRNH4CL) THEN         
            SCASE = 'G2'
            CALL SSH_CALCG2              ! NH42SO4,NH4CL,NA2SO4
C
         ELSEIF (DRNH4CL.LE.RH  .AND. RH.LT.DRNH42S4) THEN         
            SCASE = 'G3'
            CALL SSH_CALCG3              ! NH42SO4,NA2SO4
C 
        ELSEIF (DRNH42S4.LE.RH  .AND. RH.LT.DRNA2SO4) THEN         
            SCASE = 'G4'
            CALL SSH_CALCG4              ! NA2SO4
C
         ELSEIF (DRNA2SO4.LE.RH) THEN         
            SCASE = 'G5'
            CALL SSH_CALCG5              ! Only liquid
         ENDIF
      ENDIF
C
C *** SULFATE POOR ; SODIUM RICH
C
      ELSE IF (SULRAT.GE.2.0 .AND. SODRAT.GE.2.0) THEN                
C
      IF(METSTBL.EQ.1) THEN
         SCASE = 'H6'
         CALL SSH_CALCH6                 ! Only liquid (metastable)
      ELSE
C
         IF (RH.LT.DRNH4NO3) THEN    
            SCASE = 'H1'
            CALL SSH_CALCH1              ! NH4NO3,NH4CL,NA2SO4,NACL,NANO3
C
         ELSEIF (DRNH4NO3.LE.RH .AND. RH.LT.DRNANO3) THEN         
            SCASE = 'H2'
            CALL SSH_CALCH2              ! NH4CL,NA2SO4,NACL,NANO3
C
         ELSEIF (DRNANO3.LE.RH  .AND. RH.LT.DRNACL) THEN         
            SCASE = 'H3'
            CALL SSH_CALCH3              ! NH4CL,NA2SO4,NACL
C
         ELSEIF (DRNACL.LE.RH   .AND. RH.LT.DRNH4Cl) THEN         
            SCASE = 'H4'
            CALL SSH_CALCH4              ! NH4CL,NA2SO4
C
         ELSEIF (DRNH4Cl.LE.RH .AND. RH.LT.DRNA2SO4) THEN         
            SCASE = 'H5'
            CALL SSH_CALCH5              ! NA2SO4
C
         ELSEIF (DRNA2SO4.LE.RH) THEN         
            SCASE = 'H6'
            CALL SSH_CALCH6              ! NO SOLID
         ENDIF
      ENDIF
C
C *** SULFATE RICH (NO ACID) 
C
      ELSEIF (1.0.LE.SULRAT .AND. SULRAT.LT.2.0) THEN 
C
      IF(METSTBL.EQ.1) THEN
         SCASE = 'I6'
         CALL SSH_CALCI6                 ! Only liquid (metastable)
      ELSE
C
         IF (RH.LT.DRNH4HS4) THEN         
            SCASE = 'I1'
            CALL SSH_CALCI1              ! NA2SO4,(NH4)2SO4,NAHSO4,NH4HSO4,LC
C
         ELSEIF (DRNH4HS4.LE.RH .AND. RH.LT.DRNAHSO4) THEN         
            SCASE = 'I2'
            CALL SSH_CALCI2              ! NA2SO4,(NH4)2SO4,NAHSO4,LC
C
         ELSEIF (DRNAHSO4.LE.RH .AND. RH.LT.DRLC) THEN         
            SCASE = 'I3'
            CALL SSH_CALCI3              ! NA2SO4,(NH4)2SO4,LC
C
         ELSEIF (DRLC.LE.RH     .AND. RH.LT.DRNH42S4) THEN         
            SCASE = 'I4'
            CALL SSH_CALCI4              ! NA2SO4,(NH4)2SO4
C
         ELSEIF (DRNH42S4.LE.RH .AND. RH.LT.DRNA2SO4) THEN         
            SCASE = 'I5'
            CALL SSH_CALCI5              ! NA2SO4
C
         ELSEIF (DRNA2SO4.LE.RH) THEN         
            SCASE = 'I6'
            CALL SSH_CALCI6              ! NO SOLIDS
         ENDIF
      ENDIF
C                                    
      CALL SSH_CALCNHA                ! MINOR SPECIES: HNO3, HCl       
      CALL SSH_CALCNH3                !                NH3 
C
C *** SULFATE RICH (FREE ACID)
C
      ELSEIF (SULRAT.LT.1.0) THEN             
C
      IF(METSTBL.EQ.1) THEN
         SCASE = 'J3'
         CALL SSH_CALCJ3                 ! Only liquid (metastable)
      ELSE
C
         IF (RH.LT.DRNH4HS4) THEN         
            SCASE = 'J1'
            CALL SSH_CALCJ1              ! NH4HSO4,NAHSO4
C
         ELSEIF (DRNH4HS4.LE.RH .AND. RH.LT.DRNAHSO4) THEN         
            SCASE = 'J2'
            CALL SSH_CALCJ2              ! NAHSO4
C
         ELSEIF (DRNAHSO4.LE.RH) THEN         
            SCASE = 'J3'
            CALL SSH_CALCJ3              
         ENDIF
      ENDIF
C                                    
      CALL SSH_CALCNHA                ! MINOR SPECIES: HNO3, HCl       
      CALL SSH_CALCNH3                !                NH3 
      ENDIF
C
C *** RETURN POINT
C
      RETURN
C
C *** END OF SUBROUTINE SSH_ISRP3F *****************************************
C
      END
C=======================================================================
C
C *** ISORROPIA CODE
C *** SUBROUTINE SSH_CALCA2
C *** CASE A2 
C
C     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
C     1. SULFATE POOR (SULRAT > 2.0)
C     2. LIQUID AEROSOL PHASE ONLY POSSIBLE
C
C     FOR CALCULATIONS, A BISECTION IS PERFORMED TOWARDS X, THE
C     AMOUNT OF HYDROGEN IONS (H+) FOUND IN THE LIQUID PHASE.
C     FOR EACH ESTIMATION OF H+, FUNCTION FUNCB2A CALCULATES THE
C     CONCENTRATION OF IONS FROM THE NH3(GAS) - NH4+(LIQ) EQUILIBRIUM.
C     ELECTRONEUTRALITY IS USED AS THE OBJECTIVE FUNCTION.
C
C *** COPYRIGHT 1996-2006, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY,
C *** GEORGIA INSTITUTE OF TECHNOLOGY
C *** WRITTEN BY ATHANASIOS NENES
C *** UPDATED BY CHRISTOS FOUNTOUKIS
C
C=======================================================================
C
      SUBROUTINE SSH_CALCA2
      INCLUDE 'isrpia.inc'
C
C *** SETUP PARAMETERS ************************************************
C
      CALAOU    =.TRUE.       ! Outer loop activity calculation flag
      OMELO     = TINY        ! Low  limit: SOLUTION IS VERY BASIC
      OMEHI     = 2.0D0*W(2)  ! High limit: FROM NH4+ -> NH3(g) + H+(aq)
C
C *** CALCULATE WATER CONTENT *****************************************
C
      MOLAL(5) = W(2)
      MOLAL(6) = ZERO
      CALL SSH_CALCMR
C
C *** INITIAL VALUES FOR BISECTION ************************************
C
      X1 = OMEHI
      Y1 = SSH_FUNCA2 (X1)
      IF (ABS(Y1).LE.EPS) RETURN
C
C *** ROOT TRACKING ; FOR THE RANGE OF HI AND LO **********************
C
      DX = (OMEHI-OMELO)/FLOAT(NDIV)
      DO 10 I=1,NDIV
         X2 = MAX(X1-DX, OMELO)
         Y2 = SSH_FUNCA2 (X2)
         IF (SIGN(1.d0,Y1)*SIGN(1.d0,Y2).LT.ZERO) GOTO 20  ! (Y1*Y2.LT.ZERO)
         X1 = X2
         Y1 = Y2
10    CONTINUE
      IF (ABS(Y2).LE.EPS) THEN
         RETURN
      ELSE
         CALL SSH_PUSHERR (0001, 'CALCA2')    ! WARNING ERROR: NO SOLUTION
         RETURN
      ENDIF
C
C *** PERFORM BISECTION ***********************************************
C
20    DO 30 I=1,MAXIT
         X3 = 0.5*(X1+X2)
         Y3 = SSH_FUNCA2 (X3)
         IF (SIGN(1.d0,Y1)*SIGN(1.d0,Y3) .LE. ZERO) THEN  ! (Y1*Y3 .LE. ZERO)
            Y2    = Y3
            X2    = X3
         ELSE
            Y1    = Y3
            X1    = X3
         ENDIF
         IF (ABS(X2-X1) .LE. EPS*X1) GOTO 40
30    CONTINUE
      CALL SSH_PUSHERR (0002, 'CALCA2')    ! WARNING ERROR: NO CONVERGENCE
C
C *** CONVERGED ; RETURN **********************************************
C
40    X3 = 0.5*(X1+X2)
      Y3 = SSH_FUNCA2 (X3)
      RETURN
C
C *** END OF SUBROUTINE SSH_CALCA2 ****************************************
C
      END



C=======================================================================
C
C *** ISORROPIA CODE
C *** FUNCTION FUNCA2
C *** CASE A2 
C     FUNCTION THAT SOLVES THE SYSTEM OF EQUATIONS FOR CASE A2 ; 
C     AND RETURNS THE VALUE OF THE ZEROED FUNCTION IN FUNCA2.
C
C=======================================================================
C
      DOUBLE PRECISION FUNCTION SSH_FUNCA2 (OMEGI)
      INCLUDE 'isrpia.inc'
      DOUBLE PRECISION LAMDA
C
C *** SETUP PARAMETERS ************************************************
C
      FRST   = .TRUE.
      CALAIN = .TRUE.
      PSI    = W(2)         ! INITIAL AMOUNT OF (NH4)2SO4 IN SOLUTION
C
C *** SOLVE EQUATIONS ; WITH ITERATIONS FOR ACTIVITY COEF. ************
C
      DO 10 I=1,NSWEEP
         A1    = XK1*WATER/GAMA(7)*(GAMA(8)/GAMA(7))**2.
         A2    = XK2*R*TEMP/XKW*(GAMA(8)/GAMA(9))**2.
         A3    = XKW*RH*WATER*WATER
C
         LAMDA = PSI/(A1/OMEGI+ONE)
         ZETA  = A3/OMEGI
C
C *** SPECIATION & WATER CONTENT ***************************************
C
         MOLAL (1) = OMEGI                        ! HI
         MOLAL (3) = MAX(W(3)/(ONE/A2/OMEGI + ONE), 2.*MOLAL(5))  ! NH4I
         MOLAL (5) = MAX(PSI-LAMDA,TINY)          ! SO4I
         MOLAL (6) = LAMDA                        ! HSO4I
         GNH3      = MAX (W(3)-MOLAL(3), TINY)    ! NH3GI
         COH       = ZETA                         ! OHI
C
C *** CALCULATE ACTIVITIES OR TERMINATE INTERNAL LOOP *****************
C
         IF (FRST.AND.CALAOU .OR. .NOT.FRST.AND.CALAIN) THEN
            CALL SSH_CALCACT     
         ELSE
            GOTO 20
         ENDIF
10    CONTINUE
C
C *** CALCULATE OBJECTIVE FUNCTION ************************************
C

C *** FOR INTERACTION WITH HYDROPHILIC ORGANIC COMPOUNDS
20    DENOM = (2.0*MOLAL(5)+MOLAL(6)+ORGANION)

      SSH_FUNCA2= (MOLAL(3)/DENOM - ONE) + MOLAL(1)/DENOM
      RETURN
C
C *** END OF FUNCTION FUNCA2********************************************
C
      END
C=======================================================================
C
C *** ISORROPIA CODE
C *** SUBROUTINE SSH_CALCA1
C *** CASE A1 
C
C     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
C     1. SULFATE POOR (SULRAT > 2.0)
C     2. SOLID AEROSOL ONLY
C     3. SOLIDS POSSIBLE : (NH4)2SO4
C
C     A SIMPLE MATERIAL BALANCE IS PERFORMED, AND THE SOLID (NH4)2SO4
C     IS CALCULATED FROM THE SULFATES. THE EXCESS AMMONIA REMAINS IN
C     THE GAS PHASE.
C
C *** COPYRIGHT 1996-2006, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY,
C *** GEORGIA INSTITUTE OF TECHNOLOGY
C *** WRITTEN BY ATHANASIOS NENES
C *** UPDATED BY CHRISTOS FOUNTOUKIS
C
C=======================================================================
C
      SUBROUTINE SSH_CALCA1
      INCLUDE 'isrpia.inc'
C
      CNH42S4 = W(2)
      GNH3    = MAX (W(3)-2.0*CNH42S4, ZERO)
      RETURN
C
C *** END OF SUBROUTINE SSH_CALCA1 ******************************************
C
      END



C=======================================================================
C
C *** ISORROPIA CODE
C *** SUBROUTINE SSH_CALCB4
C *** CASE B4 
C
C     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
C     1. SULFATE RICH, NO FREE ACID (1.0 <= SULRAT < 2.0)
C     2. LIQUID AEROSOL PHASE ONLY POSSIBLE
C
C     FOR CALCULATIONS, A BISECTION IS PERFORMED WITH RESPECT TO H+.
C     THE OBJECTIVE FUNCTION IS THE DIFFERENCE BETWEEN THE ESTIMATED H+
C     AND THAT CALCULATED FROM ELECTRONEUTRALITY.
C
C *** COPYRIGHT 1996-2006, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY,
C *** GEORGIA INSTITUTE OF TECHNOLOGY
C *** WRITTEN BY ATHANASIOS NENES
C *** UPDATED BY CHRISTOS FOUNTOUKIS
C
C=======================================================================
C
      SUBROUTINE SSH_CALCB4
      INCLUDE 'isrpia.inc'
C
C *** SOLVE EQUATIONS **************************************************
C
      FRST       = .TRUE.
      CALAIN     = .TRUE.
      CALAOU     = .TRUE.
C
C *** CALCULATE WATER CONTENT ******************************************
C
      CALL SSH_CALCB1A         ! GET DRY SALT CONTENT, AND USE FOR WATER.
      MOLALR(13) = CLC       
      MOLALR(9)  = CNH4HS4   
      MOLALR(4)  = CNH42S4   
      CLC        = ZERO
      CNH4HS4    = ZERO
      CNH42S4    = ZERO
      WATER      = MOLALR(13)/M0(13)+MOLALR(9)/M0(9)+MOLALR(4)/M0(4)
C *** FOR INTERACTION WITH HYDROPHILIC ORGANIC COMPOUNDS
C ADD ORGANIC LIQUID WATER CONTENT
      WATER = WATER + WATORG
C
      MOLAL(3)   = W(3)   ! NH4I
C
      DO 20 I=1,NSWEEP
         AK1   = XK1*((GAMA(8)/GAMA(7))**2.)*(WATER/GAMA(7))
         BET   = W(2)
         GAM   = MOLAL(3)
C
C *** FOR INTERACTION WITH HYDROPHILIC ORGANIC COMPOUNDS
         BB    = BET + AK1 - GAM + ORGANION
         CC    =-AK1*BET
         DD    = BB*BB - 4.D0*CC
C
C *** SPECIATION & WATER CONTENT ***************************************
C
         MOLAL (5) = MAX(TINY,MIN(0.5*(-BB + SQRT(DD)), W(2))) ! SO4I
         MOLAL (6) = MAX(TINY,MIN(W(2)-MOLAL(5),W(2)))         ! HSO4I
         MOLAL (1) = MAX(TINY,MIN(AK1*MOLAL(6)/MOLAL(5),W(2))) ! HI
         CALL SSH_CALCMR                                           ! Water content
C
C *** CALCULATE ACTIVITIES OR TERMINATE INTERNAL LOOP *****************
C
         IF (.NOT.CALAIN) GOTO 30
         CALL SSH_CALCACT
20    CONTINUE
C
30    RETURN
C
C *** END OF SUBROUTINE SSH_CALCB4 ******************************************
C
      END
C=======================================================================
C
C *** ISORROPIA CODE
C *** SUBROUTINE SSH_CALCB3
C *** CASE B3 
C
C     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
C     1. SULFATE RICH, NO FREE ACID (1.0 <= SULRAT < 2.0)
C     2. BOTH LIQUID & SOLID PHASE IS POSSIBLE
C     3. SOLIDS POSSIBLE: (NH4)2SO4
C
C *** COPYRIGHT 1996-2006, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY,
C *** GEORGIA INSTITUTE OF TECHNOLOGY
C *** WRITTEN BY ATHANASIOS NENES
C *** UPDATED BY CHRISTOS FOUNTOUKIS
C
C=======================================================================
C
      SUBROUTINE SSH_CALCB3
      INCLUDE 'isrpia.inc'
C    
C *** CALCULATE EQUIVALENT AMOUNT OF HSO4 AND SO4 ***********************
C
      X = MAX(2*W(2)-W(3), ZERO)   ! Equivalent NH4HSO4
      Y = MAX(W(3)  -W(2), ZERO)   ! Equivalent NH42SO4
C
C *** CALCULATE SPECIES ACCORDING TO RELATIVE ABUNDANCE OF HSO4 *********
C
      IF (X.LT.Y) THEN             ! LC is the MIN (x,y)
         SCASE   = 'B3 ; SUBCASE 1'
         TLC     = X
         TNH42S4 = Y-X
         CALL SSH_CALCB3A (TLC,TNH42S4)      ! LC + (NH4)2SO4 
      ELSE
         SCASE   = 'B3 ; SUBCASE 2'
         TLC     = Y
         TNH4HS4 = X-Y
         CALL SSH_CALCB3B (TLC,TNH4HS4)      ! LC + NH4HSO4
      ENDIF
C 
      RETURN
C
C *** END OF SUBROUTINE SSH_CALCB3 ******************************************
C
      END


C=======================================================================
C
C *** ISORROPIA CODE
C *** SUBROUTINE SSH_CALCB3A
C *** CASE B3 ; SUBCASE 1
C
C     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
C     1. SULFATE RICH (1.0 < SULRAT < 2.0)
C     2. BOTH LIQUID & SOLID PHASE IS POSSIBLE
C     3. SOLIDS POSSIBLE: (NH4)2SO4
C
C     FOR CALCULATIONS, A BISECTION IS PERFORMED TOWARDS ZETA, THE
C     AMOUNT OF SOLID (NH4)2SO4 DISSOLVED IN THE LIQUID PHASE.
C     FOR EACH ESTIMATION OF ZETA, FUNCTION FUNCB3A CALCULATES THE
C     AMOUNT OF H+ PRODUCED (BASED ON THE SO4 RELEASED INTO THE
C     SOLUTION). THE SOLUBILITY PRODUCT OF (NH4)2SO4 IS USED AS THE 
C     OBJECTIVE FUNCTION.
C
C *** COPYRIGHT 1996-2006, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY,
C *** GEORGIA INSTITUTE OF TECHNOLOGY
C *** WRITTEN BY ATHANASIOS NENES
C *** UPDATED BY CHRISTOS FOUNTOUKIS
C
C=======================================================================
C
      SUBROUTINE SSH_CALCB3A (TLC, TNH42S4)
      INCLUDE 'isrpia.inc'
C
      CALAOU = .TRUE.         ! Outer loop activity calculation flag
      ZLO    = ZERO           ! MIN DISSOLVED (NH4)2SO4
      ZHI    = TNH42S4        ! MAX DISSOLVED (NH4)2SO4
C
C *** INITIAL VALUES FOR BISECTION (DISSOLVED (NH4)2SO4 ****************
C
      Z1 = ZLO
      Y1 = SSH_FUNCB3A (Z1, TLC, TNH42S4)
      IF (ABS(Y1).LE.EPS) RETURN
      YLO= Y1
C
C *** ROOT TRACKING ; FOR THE RANGE OF HI AND LO ***********************
C
      DZ = (ZHI-ZLO)/FLOAT(NDIV)
      DO 10 I=1,NDIV
         Z2 = Z1+DZ
         Y2 = SSH_FUNCB3A (Z2, TLC, TNH42S4)
         IF (SIGN(1.d0,Y1)*SIGN(1.d0,Y2).LT.ZERO) GOTO 20  ! (Y1*Y2.LT.ZERO)
         Z1 = Z2
         Y1 = Y2
10    CONTINUE
C
C *** NO SUBDIVISION WITH SOLUTION FOUND 
C
      YHI= Y1                      ! Save Y-value at HI position
      IF (ABS(Y2) .LT. EPS) THEN   ! X2 IS A SOLUTION 
         RETURN
C
C *** { YLO, YHI } < 0.0 THE SOLUTION IS ALWAYS UNDERSATURATED WITH LC
C
      ELSE IF (YLO.LT.ZERO .AND. YHI.LT.ZERO) THEN
         Z1 = ZHI
         Z2 = ZHI
         GOTO 40
C
C *** { YLO, YHI } > 0.0 THE SOLUTION IS ALWAYS SUPERSATURATED WITH LC
C
      ELSE IF (YLO.GT.ZERO .AND. YHI.GT.ZERO) THEN
         Z1 = ZLO
         Z2 = ZLO
         GOTO 40
      ELSE
         CALL SSH_PUSHERR (0001, 'CALCB3A')    ! WARNING ERROR: NO SOLUTION
         RETURN
      ENDIF
C
C *** PERFORM BISECTION ***********************************************
C
20    DO 30 I=1,MAXIT
         Z3 = 0.5*(Z1+Z2)
         Y3 = SSH_FUNCB3A (Z3, TLC, TNH42S4)
         IF (SIGN(1.d0,Y1)*SIGN(1.d0,Y3) .LE. ZERO) THEN  ! (Y1*Y3 .LE. ZERO)
            Y2    = Y3
            Z2    = Z3
         ELSE
            Y1    = Y3
            Z1    = Z3
         ENDIF
         IF (ABS(Z2-Z1) .LE. EPS*Z1) GOTO 40
30    CONTINUE
      CALL SSH_PUSHERR (0002, 'CALCB3A')    ! WARNING ERROR: NO CONVERGENCE
C
C *** CONVERGED ; RETURN ************************************************
C
40    ZK = 0.5*(Z1+Z2)
      Y3 = SSH_FUNCB3A (ZK, TLC, TNH42S4)
C    
      RETURN
C
C *** END OF SUBROUTINE SSH_CALCB3A ******************************************
C
      END



C=======================================================================
C
C *** ISORROPIA CODE
C *** FUNCTION FUNCB3A
C *** CASE B3 ; SUBCASE 1
C     FUNCTION THAT SOLVES THE SYSTEM OF EQUATIONS FOR CASE B3
C     AND RETURNS THE VALUE OF THE ZEROED FUNCTION IN FUNCA3.
C
C=======================================================================
C
      DOUBLE PRECISION FUNCTION SSH_FUNCB3A (ZK, Y, X)
      INCLUDE 'isrpia.inc'
      DOUBLE PRECISION KK
C
C *** SOLVE EQUATIONS ; WITH ITERATIONS FOR ACTIVITY COEF. ************
C
      FRST   = .TRUE.
      CALAIN = .TRUE.
      DO 20 I=1,NSWEEP
C *** FOR INTERACTION WITH HYDROPHILIC ORGANIC COMPOUNDS
         GRAT1 = XK1*WATER/GAMA(7)*(GAMA(8)/GAMA(7))**2.
         DD    = SQRT( (ZK+GRAT1+Y+ORGANION)**2. 
     &        + 4.0*(Y*GRAT1 - ORGANION*(Y+ZK)) )
         KK    = 0.5*(-(ZK+GRAT1+Y+ORGANION) + DD)
C
C *** SPECIATION & WATER CONTENT ***************************************
C
C *** FOR INTERACTION WITH HYDROPHILIC ORGANIC COMPOUNDS
         MOLAL (1) = KK + ORGANION     ! HI
         MOLAL (5) = KK+ZK+Y           ! SO4I
         MOLAL (6) = MAX (Y-KK, TINY)  ! HSO4I
         MOLAL (3) = 3.0*Y+2*ZK        ! NH4I
         CNH42S4   = X-ZK              ! Solid (NH4)2SO4
         CALL SSH_CALCMR                   ! Water content
C
C *** CALCULATE ACTIVITIES OR TERMINATE INTERNAL LOOP *****************
C
         IF (FRST.AND.CALAOU .OR. .NOT.FRST.AND.CALAIN) THEN
            CALL SSH_CALCACT     
         ELSE
            GOTO 30
         ENDIF
20    CONTINUE
C
C *** CALCULATE OBJECTIVE FUNCTION ************************************
C
CCC30    FUNCB3A= ( SO4I*NH4I**2.0 )/( XK7*(WATER/GAMA(4))**3.0 )
30    SSH_FUNCB3A= MOLAL(5)*MOLAL(3)**2.0
      SSH_FUNCB3A= SSH_FUNCB3A/(XK7*(WATER/GAMA(4))**3.0) - ONE
      RETURN
C
C *** END OF FUNCTION FUNCB3A ********************************************
C
      END



C=======================================================================
C
C *** ISORROPIA CODE
C *** SUBROUTINE SSH_CALCB3B
C *** CASE B3 ; SUBCASE 2
C
C     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
C     1. SULFATE RICH (1.0 < SULRAT < 2.0)
C     2. LIQUID PHASE ONLY IS POSSIBLE
C
C     SPECIATION CALCULATIONS IS BASED ON THE HSO4 <--> SO4 EQUILIBRIUM. 
C
C *** COPYRIGHT 1996-2006, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY,
C *** GEORGIA INSTITUTE OF TECHNOLOGY
C *** WRITTEN BY ATHANASIOS NENES
C *** UPDATED BY CHRISTOS FOUNTOUKIS
C
C=======================================================================
C
      SUBROUTINE SSH_CALCB3B (Y, X)
      INCLUDE 'isrpia.inc'
      DOUBLE PRECISION KK
C
      CALAOU = .FALSE.        ! Outer loop activity calculation flag
      FRST   = .FALSE.
      CALAIN = .TRUE.
C
C *** SOLVE EQUATIONS ; WITH ITERATIONS FOR ACTIVITY COEF. ************
C
      DO 20 I=1,NSWEEP
         GRAT1 = XK1*WATER/GAMA(7)*(GAMA(8)/GAMA(7))**2.
C *** FOR INTERACTION WITH HYDROPHILIC ORGANIC COMPOUNDS
         DD    = SQRT( (GRAT1+Y+ORGANION)**2. + 4.0* (
     &        (X+Y)*GRAT1 - Y * ORGANION))
         KK    = 0.5*(-(GRAT1+Y+ORGANION) + DD )
C
C *** SPECIATION & WATER CONTENT ***************************************
C

C *** FOR INTERACTION WITH HYDROPHILIC ORGANIC COMPOUNDS
         MOLAL (1) = KK + ORGANION        ! HI
         MOLAL (5) = Y+KK                 ! SO4I
         MOLAL (6) = MAX (X+Y-KK, TINY)   ! HSO4I
         MOLAL (3) = 3.0*Y+X              ! NH4I
         CALL SSH_CALCMR                      ! Water content
C
C *** CALCULATE ACTIVITIES OR TERMINATE INTERNAL LOOP *****************
C
         IF (.NOT.CALAIN) GOTO 30
         CALL SSH_CALCACT     
20    CONTINUE
C    
30    RETURN
C
C *** END OF SUBROUTINE SSH_CALCB3B ******************************************
C
      END
C=======================================================================
C
C *** ISORROPIA CODE
C *** SUBROUTINE SSH_CALCB2
C *** CASE B2 
C
C     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
C     1. SULFATE RICH, NO FREE ACID (1.0 <= SULRAT < 2.0)
C     2. THERE IS BOTH A LIQUID & SOLID PHASE
C     3. SOLIDS POSSIBLE : LC, (NH4)2SO4
C
C     THERE ARE TWO POSSIBLE REGIMES HERE, DEPENDING ON THE SULFATE RATIO:
C     1. WHEN BOTH LC AND (NH4)2SO4 ARE POSSIBLE (SUBROUTINE SSH_CALCB2A)
C     2. WHEN ONLY LC IS POSSIBLE (SUBROUTINE SSH_CALCB2B).
C
C *** COPYRIGHT 1996-2006, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY,
C *** GEORGIA INSTITUTE OF TECHNOLOGY
C *** WRITTEN BY ATHANASIOS NENES
C *** UPDATED BY CHRISTOS FOUNTOUKIS
C
C=======================================================================
C
      SUBROUTINE SSH_CALCB2
      INCLUDE 'isrpia.inc'
C    
C *** CALCULATE EQUIVALENT AMOUNT OF HSO4 AND SO4 ***********************
C
      X = MAX(2*W(2)-W(3), TINY)   ! Equivalent NH4HSO4
      Y = MAX(W(3)  -W(2), TINY)   ! Equivalent NH42SO4
C
C *** CALCULATE SPECIES ACCORDING TO RELATIVE ABUNDANCE OF HSO4 *********
C
      IF (X.LE.Y) THEN             ! LC is the MIN (x,y)
         SCASE = 'B2 ; SUBCASE 1'
         CALL SSH_CALCB2A (X,Y-X)      ! LC + (NH4)2SO4 POSSIBLE
      ELSE
         SCASE = 'B2 ; SUBCASE 2'
         CALL SSH_CALCB2B (Y,X-Y)      ! LC ONLY POSSIBLE
      ENDIF
C 
      RETURN
C
C *** END OF SUBROUTINE SSH_CALCB2 ******************************************
C
      END



C=======================================================================
C
C *** ISORROPIA CODE
C *** SUBROUTINE SSH_CALCB2
C *** CASE B2 ; SUBCASE A. 
C
C     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
C     1. SULFATE RICH (1.0 < SULRAT < 2.0)
C     2. SOLID PHASE ONLY POSSIBLE
C     3. SOLIDS POSSIBLE: LC, (NH4)2SO4
C
C     THERE ARE TWO POSSIBLE REGIMES HERE, DEPENDING ON RELATIVE HUMIDITY:
C     1. WHEN RH >= MDRH ; LIQUID PHASE POSSIBLE (MDRH REGION)
C     2. WHEN RH < MDRH  ; ONLY SOLID PHASE POSSIBLE 
C
C     FOR SOLID CALCULATIONS, A MATERIAL BALANCE BASED ON THE STOICHIMETRIC
C     PROPORTION OF AMMONIUM AND SULFATE IS DONE TO CALCULATE THE AMOUNT 
C     OF LC AND (NH4)2SO4 IN THE SOLID PHASE.
C
C *** COPYRIGHT 1996-2006, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY,
C *** GEORGIA INSTITUTE OF TECHNOLOGY
C *** WRITTEN BY ATHANASIOS NENES
C *** UPDATED BY CHRISTOS FOUNTOUKIS
C
C=======================================================================
C
      SUBROUTINE SSH_CALCB2A (TLC, TNH42S4)
      INCLUDE 'isrpia.inc'
C
C *** REGIME DEPENDS UPON THE AMBIENT RELATIVE HUMIDITY *****************
C
      IF (RH.LT.DRMLCAS) THEN    
         SCASE   = 'B2 ; SUBCASE A1'    ! SOLIDS POSSIBLE ONLY
         CLC     = TLC
         CNH42S4 = TNH42S4
         SCASE   = 'B2 ; SUBCASE A1'
      ELSE
         SCASE = 'B2 ; SUBCASE A2'
         CALL SSH_CALCB2A2 (TLC, TNH42S4)   ! LIQUID & SOLID PHASE POSSIBLE
         SCASE = 'B2 ; SUBCASE A2'
      ENDIF
C
      RETURN
C
C *** END OF SUBROUTINE SSH_CALCB2A *****************************************
C
      END



C=======================================================================
C
C *** ISORROPIA CODE
C *** SUBROUTINE SSH_CALCB2A2
C *** CASE B2 ; SUBCASE A2. 
C
C     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
C     1. SULFATE RICH (1.0 < SULRAT < 2.0)
C     2. SOLID PHASE ONLY POSSIBLE
C     3. SOLIDS POSSIBLE: LC, (NH4)2SO4
C
C     THIS IS THE CASE WHERE THE RELATIVE HUMIDITY IS IN THE MUTUAL
C     DRH REGION. THE SOLUTION IS ASSUMED TO BE THE SUM OF TWO WEIGHTED
C     SOLUTIONS ; THE SOLID PHASE ONLY (SUBROUTINE SSH_CALCB2A1) AND THE
C     THE SOLID WITH LIQUID PHASE (SUBROUTINE SSH_CALCB3).
C
C *** COPYRIGHT 1996-2006, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY,
C *** GEORGIA INSTITUTE OF TECHNOLOGY
C *** WRITTEN BY ATHANASIOS NENES
C *** UPDATED BY CHRISTOS FOUNTOUKIS
C
C=======================================================================
C
      SUBROUTINE SSH_CALCB2A2 (TLC, TNH42S4)
      INCLUDE 'isrpia.inc'
C
C *** FIND WEIGHT FACTOR **********************************************
C
      IF (WFTYP.EQ.0) THEN
         WF = ZERO
      ELSEIF (WFTYP.EQ.1) THEN
         WF = 0.5D0
      ELSE
         WF = (DRLC-RH)/(DRLC-DRMLCAS)
      ENDIF
      ONEMWF  = ONE - WF
C
C *** FIND FIRST SECTION ; DRY ONE ************************************
C
      CLCO     = TLC                     ! FIRST (DRY) SOLUTION
      CNH42SO  = TNH42S4
C
C *** FIND SECOND SECTION ; DRY & LIQUID ******************************
C
      CLC     = ZERO
      CNH42S4 = ZERO
      CALL SSH_CALCB3                        ! SECOND (LIQUID) SOLUTION
C
C *** FIND SOLUTION AT MDRH BY WEIGHTING DRY & LIQUID SOLUTIONS.
C
      MOLAL(1)= ONEMWF*MOLAL(1)                                   ! H+
      MOLAL(3)= ONEMWF*(2.D0*(CNH42SO-CNH42S4) + 3.D0*(CLCO-CLC)) ! NH4+
      MOLAL(5)= ONEMWF*(CNH42SO-CNH42S4 + CLCO-CLC)               ! SO4--
      MOLAL(6)= ONEMWF*(CLCO-CLC)                                 ! HSO4-
C
      WATER   = ONEMWF*WATER
C
      CLC     = WF*CLCO    + ONEMWF*CLC
      CNH42S4 = WF*CNH42SO + ONEMWF*CNH42S4
C
      RETURN
C
C *** END OF SUBROUTINE SSH_CALCB2A2 ****************************************
C
      END



C=======================================================================
C
C *** ISORROPIA CODE
C *** SUBROUTINE SSH_CALCB2
C *** CASE B2 ; SUBCASE B 
C
C     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
C     1. SULFATE RICH (1.0 < SULRAT < 2.0)
C     2. BOTH LIQUID & SOLID PHASE IS POSSIBLE
C     3. SOLIDS POSSIBLE: LC
C
C     FOR CALCULATIONS, A BISECTION IS PERFORMED TOWARDS ZETA, THE
C     AMOUNT OF SOLID LC DISSOLVED IN THE LIQUID PHASE.
C     FOR EACH ESTIMATION OF ZETA, FUNCTION FUNCB2A CALCULATES THE
C     AMOUNT OF H+ PRODUCED (BASED ON THE HSO4, SO4 RELEASED INTO THE
C     SOLUTION). THE SOLUBILITY PRODUCT OF LC IS USED AS THE OBJECTIVE 
C     FUNCTION.
C
C *** COPYRIGHT 1996-2006, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY,
C *** GEORGIA INSTITUTE OF TECHNOLOGY
C *** WRITTEN BY ATHANASIOS NENES
C *** UPDATED BY CHRISTOS FOUNTOUKIS
C
C=======================================================================
C
      SUBROUTINE SSH_CALCB2B (TLC,TNH4HS4)
      INCLUDE 'isrpia.inc'
C
      CALAOU = .TRUE.       ! Outer loop activity calculation flag
      ZLO    = ZERO
      ZHI    = TLC          ! High limit: all of it in liquid phase
C
C *** INITIAL VALUES FOR BISECTION **************************************
C
      X1 = ZHI
      Y1 = SSH_FUNCB2B (X1,TNH4HS4,TLC)
      IF (ABS(Y1).LE.EPS) RETURN
      YHI= Y1                        ! Save Y-value at Hi position
C
C *** ROOT TRACKING ; FOR THE RANGE OF HI AND LO ************************
C
      DX = (ZHI-ZLO)/NDIV
      DO 10 I=1,NDIV
         X2 = X1-DX
         Y2 = SSH_FUNCB2B (X2,TNH4HS4,TLC)
         IF (SIGN(1.d0,Y1)*SIGN(1.d0,Y2).LT.ZERO) GOTO 20  ! (Y1*Y2.LT.ZERO)
         X1 = X2
         Y1 = Y2
10    CONTINUE
C
C *** NO SUBDIVISION WITH SOLUTION FOUND 
C
      YLO= Y1                      ! Save Y-value at LO position
      IF (ABS(Y2) .LT. EPS) THEN   ! X2 IS A SOLUTION 
         RETURN
C
C *** { YLO, YHI } < 0.0 THE SOLUTION IS ALWAYS UNDERSATURATED WITH LC
C
      ELSE IF (YLO.LT.ZERO .AND. YHI.LT.ZERO) THEN
         X1 = ZHI
         X2 = ZHI
         GOTO 40
C
C *** { YLO, YHI } > 0.0 THE SOLUTION IS ALWAYS SUPERSATURATED WITH LC
C
      ELSE IF (YLO.GT.ZERO .AND. YHI.GT.ZERO) THEN
         X1 = ZLO
         X2 = ZLO
         GOTO 40
      ELSE
         CALL SSH_PUSHERR (0001, 'CALCB2B')    ! WARNING ERROR: NO SOLUTION
         RETURN
      ENDIF
C
C *** PERFORM BISECTION *************************************************
C
20    DO 30 I=1,MAXIT
         X3 = 0.5*(X1+X2)
         Y3 = SSH_FUNCB2B (X3,TNH4HS4,TLC)
         IF (SIGN(1.d0,Y1)*SIGN(1.d0,Y3) .LE. ZERO) THEN  ! (Y1*Y3 .LE. ZERO)
            Y2    = Y3
            X2    = X3
         ELSE
            Y1    = Y3
            X1    = X3
         ENDIF
         IF (ABS(X2-X1) .LE. EPS*X1) GOTO 40
30    CONTINUE
      CALL SSH_PUSHERR (0002, 'CALCB2B')    ! WARNING ERROR: NO CONVERGENCE
C
C *** CONVERGED ; RETURN ************************************************
C
40    X3 = 0.5*(X1+X2)
      Y3 = SSH_FUNCB2B (X3,TNH4HS4,TLC)
C
      RETURN
C
C *** END OF SUBROUTINE SSH_CALCB2B *****************************************
C
      END



C=======================================================================
C
C *** ISORROPIA CODE
C *** FUNCTION SSH_FUNCB2B
C *** CASE B2 ; 
C     FUNCTION THAT SOLVES THE SYSTEM OF EQUATIONS FOR CASE B2 ; SUBCASE 2
C     AND RETURNS THE VALUE OF THE ZEROED FUNCTION IN FUNCB2B.
C
C=======================================================================
C
      DOUBLE PRECISION FUNCTION SSH_FUNCB2B (X,TNH4HS4,TLC)
      INCLUDE 'isrpia.inc'
C
C *** SOLVE EQUATIONS **************************************************
C
      FRST   = .TRUE.
      CALAIN = .TRUE.
      DO 20 I=1,NSWEEP
         GRAT2 = XK1*WATER*(GAMA(8)/GAMA(7))**2./GAMA(7)
C *** FOR INTERACTION WITH HYDROPHILIC ORGANIC COMPOUNDS
         PARM  = X+GRAT2+ORGANION
         DELTA = PARM*PARM + 4.0*((X+TNH4HS4)*GRAT2-X*ORGANION)
         OMEGA = 0.5*(-PARM + SQRT(DELTA))         ! Thetiki riza (ie:H+>0)
C
C *** SPECIATION & WATER CONTENT ***************************************
C
C *** FOR INTERACTION WITH HYDROPHILIC ORGANIC COMPOUNDS
         MOLAL (1) = OMEGA + ORGANION              ! HI
         MOLAL (3) = 3.0*X+TNH4HS4                 ! NH4I
         MOLAL (5) = X+OMEGA                       ! SO4I
         MOLAL (6) = MAX (X+TNH4HS4-OMEGA, TINY)   ! HSO4I
         CLC       = MAX(TLC-X,ZERO)               ! Solid LC
         CALL SSH_CALCMR                               ! Water content
C
C *** CALCULATE ACTIVITIES OR TERMINATE INTERNAL LOOP ******************
C
         IF (FRST.AND.CALAOU .OR. .NOT.FRST.AND.CALAIN) THEN
            CALL SSH_CALCACT     
         ELSE
            GOTO 30
         ENDIF
20    CONTINUE
C
C *** CALCULATE OBJECTIVE FUNCTION **************************************
C
CCC30    FUNCB2B= ( NH4I**3.*SO4I*HSO4I )/( XK13*(WATER/GAMA(13))**5. )
30    SSH_FUNCB2B= (MOLAL(3)**3.)*MOLAL(5)*MOLAL(6)
      SSH_FUNCB2B= SSH_FUNCB2B/(XK13*(WATER/GAMA(13))**5.) - ONE
      RETURN
C
C *** END OF FUNCTION FUNCB2B *******************************************
C
      END


C=======================================================================
C
C *** ISORROPIA CODE
C *** SUBROUTINE SSH_CALCB1
C *** CASE B1
C
C     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
C     1. SULFATE RICH, NO FREE ACID (1.0 <= SULRAT < 2.0)
C     2. THERE IS BOTH A LIQUID & SOLID PHASE
C     3. SOLIDS POSSIBLE : LC, (NH4)2SO4, NH4HSO4
C
C     THERE ARE TWO POSSIBLE REGIMES HERE, DEPENDING ON RELATIVE HUMIDITY:
C     1. WHEN RH >= MDRH ; LIQUID PHASE POSSIBLE (MDRH REGION)
C     2. WHEN RH < MDRH  ; ONLY SOLID PHASE POSSIBLE (SUBROUTINE SSH_CALCB1A)
C
C *** COPYRIGHT 1996-2006, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY,
C *** GEORGIA INSTITUTE OF TECHNOLOGY
C *** WRITTEN BY ATHANASIOS NENES
C *** UPDATED BY CHRISTOS FOUNTOUKIS
C
C=======================================================================
C
      SUBROUTINE SSH_CALCB1
      INCLUDE 'isrpia.inc'
C
C *** REGIME DEPENDS UPON THE AMBIENT RELATIVE HUMIDITY *****************
C
      IF (RH.LT.DRMLCAB) THEN    
         SCASE = 'B1 ; SUBCASE 1'  
         CALL SSH_CALCB1A              ! SOLID PHASE ONLY POSSIBLE
         SCASE = 'B1 ; SUBCASE 1'
      ELSE
         SCASE = 'B1 ; SUBCASE 2'
         CALL SSH_CALCB1B              ! LIQUID & SOLID PHASE POSSIBLE
         SCASE = 'B1 ; SUBCASE 2'
      ENDIF
C
      RETURN
C
C *** END OF SUBROUTINE SSH_CALCB1 ******************************************
C
      END



C=======================================================================
C
C *** ISORROPIA CODE
C *** SUBROUTINE SSH_CALCB1A
C *** CASE B1 ; SUBCASE 1
C
C     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
C     1. SULFATE RICH
C     2. THERE IS NO LIQUID PHASE
C     3. SOLIDS POSSIBLE: LC, { (NH4)2SO4  XOR  NH4HSO4 } (ONE OF TWO
C                         BUT NOT BOTH)
C
C     A SIMPLE MATERIAL BALANCE IS PERFORMED, AND THE AMOUNT OF LC
C     IS CALCULATED FROM THE (NH4)2SO4 AND NH4HSO4 WHICH IS LEAST
C     ABUNDANT (STOICHIMETRICALLY). THE REMAINING EXCESS OF SALT 
C     IS MIXED WITH THE LC.  
C
C *** COPYRIGHT 1996-2006, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY,
C *** GEORGIA INSTITUTE OF TECHNOLOGY
C *** WRITTEN BY ATHANASIOS NENES
C *** UPDATED BY CHRISTOS FOUNTOUKIS
C
C=======================================================================
C
      SUBROUTINE SSH_CALCB1A
      INCLUDE 'isrpia.inc'
C
C *** SETUP PARAMETERS ************************************************
C
      X = 2*W(2)-W(3)       ! Equivalent NH4HSO4
      Y = W(3)-W(2)         ! Equivalent (NH4)2SO4
C
C *** CALCULATE COMPOSITION *******************************************
C
      IF (X.LE.Y) THEN      ! LC is the MIN (x,y)
         CLC     = X        ! NH4HSO4 >= (NH4)2S04
         CNH4HS4 = ZERO
         CNH42S4 = Y-X
      ELSE
         CLC     = Y        ! NH4HSO4 <  (NH4)2S04
         CNH4HS4 = X-Y
         CNH42S4 = ZERO
      ENDIF
      RETURN
C
C *** END OF SUBROUTINE SSH_CALCB1 ******************************************
C
      END




C=======================================================================
C
C *** ISORROPIA CODE
C *** SUBROUTINE SSH_CALCB1B
C *** CASE B1 ; SUBCASE 2
C
C     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
C     1. SULFATE RICH
C     2. THERE IS BOTH A LIQUID & SOLID PHASE
C     3. SOLIDS POSSIBLE: LC, { (NH4)2SO4  XOR  NH4HSO4 } (ONE OF TWO
C                         BUT NOT BOTH)
C
C     THIS IS THE CASE WHERE THE RELATIVE HUMIDITY IS IN THE MUTUAL
C     DRH REGION. THE SOLUTION IS ASSUMED TO BE THE SUM OF TWO WEIGHTED
C     SOLUTIONS ; THE SOLID PHASE ONLY (SUBROUTINE SSH_CALCB1A) AND THE
C     THE SOLID WITH LIQUID PHASE (SUBROUTINE SSH_CALCB2).
C
C *** COPYRIGHT 1996-2006, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY,
C *** GEORGIA INSTITUTE OF TECHNOLOGY
C *** WRITTEN BY ATHANASIOS NENES
C *** UPDATED BY CHRISTOS FOUNTOUKIS
C
C=======================================================================
C
      SUBROUTINE SSH_CALCB1B
      INCLUDE 'isrpia.inc'
C
C *** FIND WEIGHT FACTOR **********************************************
C
      IF (WFTYP.EQ.0) THEN
         WF = ZERO
      ELSEIF (WFTYP.EQ.1) THEN
         WF = 0.5D0
      ELSE
         WF = (DRNH4HS4-RH)/(DRNH4HS4-DRMLCAB)
      ENDIF
      ONEMWF  = ONE - WF
C
C *** FIND FIRST SECTION ; DRY ONE ************************************
C
      CALL SSH_CALCB1A
      CLCO     = CLC               ! FIRST (DRY) SOLUTION
      CNH42SO  = CNH42S4
      CNH4HSO  = CNH4HS4
C
C *** FIND SECOND SECTION ; DRY & LIQUID ******************************
C
      CLC     = ZERO
      CNH42S4 = ZERO
      CNH4HS4 = ZERO
      CALL SSH_CALCB2                  ! SECOND (LIQUID) SOLUTION
C
C *** FIND SOLUTION AT MDRH BY WEIGHTING DRY & LIQUID SOLUTIONS.
C
      MOLAL(1)= ONEMWF*MOLAL(1)                                   ! H+
      MOLAL(3)= ONEMWF*(2.D0*(CNH42SO-CNH42S4) + (CNH4HSO-CNH4HS4)  
     &                + 3.D0*(CLCO-CLC))                          ! NH4+
      MOLAL(5)= ONEMWF*(CNH42SO-CNH42S4 + CLCO-CLC)               ! SO4--
      MOLAL(6)= ONEMWF*(CNH4HSO-CNH4HS4 + CLCO-CLC)               ! HSO4-
C
      WATER   = ONEMWF*WATER
C
      CLC     = WF*CLCO    + ONEMWF*CLC
      CNH42S4 = WF*CNH42SO + ONEMWF*CNH42S4
      CNH4HS4 = WF*CNH4HSO + ONEMWF*CNH4HS4
C
      RETURN
C
C *** END OF SUBROUTINE SSH_CALCB1B *****************************************
C
      END


C=======================================================================
C
C *** ISORROPIA CODE
C *** SUBROUTINE SSH_CALCC2
C *** CASE C2 
C
C     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
C     1. SULFATE RICH, FREE ACID (SULRAT < 1.0)
C     2. THERE IS ONLY A LIQUID PHASE
C
C *** COPYRIGHT 1996-2006, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY,
C *** GEORGIA INSTITUTE OF TECHNOLOGY
C *** WRITTEN BY ATHANASIOS NENES
C *** UPDATED BY CHRISTOS FOUNTOUKIS
C
C=======================================================================
C
      SUBROUTINE SSH_CALCC2
      INCLUDE 'isrpia.inc'
      DOUBLE PRECISION LAMDA, KAPA
C
      CALAOU =.TRUE.         ! Outer loop activity calculation flag
      FRST   =.TRUE.
      CALAIN =.TRUE.
C
C *** SOLVE EQUATIONS **************************************************
C
      LAMDA  = W(3)           ! NH4HSO4 INITIALLY IN SOLUTION
      PSI    = W(2)-W(3)      ! H2SO4 IN SOLUTION
      DO 20 I=1,NSWEEP
         PARM  = WATER*XK1/GAMA(7)*(GAMA(8)/GAMA(7))**2.
C *** FOR INTERACTION WITH HYDROPHILIC ORGANIC COMPOUNDS
         BB    = PSI+PARM+ORGANION
         CC    =-PARM*(LAMDA+PSI)
         KAPA  = 0.5*(-BB+SQRT(BB*BB-4.0*CC))
C
C *** SPECIATION & WATER CONTENT ***************************************
C
C *** FOR INTERACTION WITH HYDROPHILIC ORGANIC COMPOUNDS
         MOLAL(1) = PSI+KAPA + ORGANION                    ! HI
         MOLAL(3) = LAMDA                                  ! NH4I
         MOLAL(5) = KAPA                                   ! SO4I
         MOLAL(6) = MAX(LAMDA+PSI-KAPA, TINY)              ! HSO4I
         CH2SO4   = MAX(MOLAL(5)+MOLAL(6)-MOLAL(3), ZERO)  ! Free H2SO4
         CALL SSH_CALCMR                                       ! Water content
C
C *** CALCULATE ACTIVITIES OR TERMINATE INTERNAL LOOP *****************
C
         IF (.NOT.CALAIN) GOTO 30
         CALL SSH_CALCACT     
20    CONTINUE
C 
30    RETURN
C    
C *** END OF SUBROUTINE SSH_CALCC2 *****************************************
C
      END



C=======================================================================
C
C *** ISORROPIA CODE
C *** SUBROUTINE SSH_CALCC1
C *** CASE C1 
C
C     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
C     1. SULFATE RICH, FREE ACID (SULRAT < 1.0)
C     2. THERE IS BOTH A LIQUID & SOLID PHASE
C     3. SOLIDS POSSIBLE: NH4HSO4
C
C *** COPYRIGHT 1996-2006, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY,
C *** GEORGIA INSTITUTE OF TECHNOLOGY
C *** WRITTEN BY ATHANASIOS NENES
C *** UPDATED BY CHRISTOS FOUNTOUKIS
C
C=======================================================================
C
      SUBROUTINE SSH_CALCC1
      INCLUDE 'isrpia.inc'
      DOUBLE PRECISION KLO, KHI
C
      CALAOU = .TRUE.    ! Outer loop activity calculation flag
      KLO    = TINY    
      KHI    = W(3)
C
C *** INITIAL VALUES FOR BISECTION *************************************
C
      X1 = KLO
      Y1 = SSH_FUNCC1 (X1)
      IF (ABS(Y1).LE.EPS) GOTO 50
      YLO= Y1
C
C *** ROOT TRACKING ; FOR THE RANGE OF HI AND LO ***********************
C
      DX = (KHI-KLO)/FLOAT(NDIV)
      DO 10 I=1,NDIV
         X2 = X1+DX
         Y2 = SSH_FUNCC1 (X2)
         IF (SIGN(1.d0,Y1)*SIGN(1.d0,Y2) .LT. ZERO) GOTO 20 ! (Y1*Y2 .LT. ZERO)
         X1 = X2
         Y1 = Y2
10    CONTINUE
C
C *** NO SUBDIVISION WITH SOLUTION FOUND 
C
      YHI= Y2                 ! Save Y-value at HI position
      IF (ABS(Y2) .LT. EPS) THEN   ! X2 IS A SOLUTION 
         GOTO 50
C
C *** { YLO, YHI } < 0.0  SOLUTION IS ALWAYS UNDERSATURATED WITH NH4HS04
C
      ELSE IF (YLO.LT.ZERO .AND. YHI.LT.ZERO) THEN
         GOTO 50
C
C *** { YLO, YHI } > 0.0 SOLUTION IS ALWAYS SUPERSATURATED WITH NH4HS04
C
      ELSE IF (YLO.GT.ZERO .AND. YHI.GT.ZERO) THEN
         X1 = KLO
         X2 = KLO
         GOTO 40
      ELSE
         CALL SSH_PUSHERR (0001, 'CALCC1')    ! WARNING ERROR: NO SOLUTION
         GOTO 50
      ENDIF
C
C *** PERFORM BISECTION OF DISSOLVED NH4HSO4 **************************
C
20    DO 30 I=1,MAXIT
         X3 = 0.5*(X1+X2)
         Y3 = SSH_FUNCC1 (X3)
         IF (SIGN(1.d0,Y1)*SIGN(1.d0,Y3) .LE. ZERO) THEN  ! (Y1*Y3 .LE. ZERO)
            Y2    = Y3
            X2    = X3
         ELSE
            Y1    = Y3
            X1    = X3
         ENDIF
         IF (ABS(X2-X1) .LE. EPS*X1) GOTO 40
30    CONTINUE
      CALL SSH_PUSHERR (0002, 'CALCC1')    ! WARNING ERROR: NO CONVERGENCE
C
C *** CONVERGED ; RETURN ***********************************************
C
40    X3 = 0.5*(X1+X2)
      Y3 = SSH_FUNCC1 (X3)
C
50    RETURN
C
C *** END OF SUBROUTINE SSH_CALCC1 *****************************************
C
      END



C=======================================================================
C
C *** ISORROPIA CODE
C *** FUNCTION FUNCC1
C *** CASE C1 ; 
C     FUNCTION THAT SOLVES THE SYSTEM OF EQUATIONS FOR CASE C1
C     AND RETURNS THE VALUE OF THE ZEROED FUNCTION IN FUNCC1.
C
C *** COPYRIGHT 1996-2006, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY,
C *** GEORGIA INSTITUTE OF TECHNOLOGY
C *** WRITTEN BY ATHANASIOS NENES
C *** UPDATED BY CHRISTOS FOUNTOUKIS
C
C=======================================================================
C
      DOUBLE PRECISION FUNCTION SSH_FUNCC1 (KAPA)
      INCLUDE 'isrpia.inc'
      DOUBLE PRECISION KAPA, LAMDA
C
C *** SOLVE EQUATIONS **************************************************
C
      FRST   = .TRUE.
      CALAIN = .TRUE.
C
      PSI = W(2)-W(3)
      DO 20 I=1,NSWEEP
         PAR1  = XK1*WATER/GAMA(7)*(GAMA(8)/GAMA(7))**2.0
         PAR2  = XK12*(WATER/GAMA(9))**2.0
C *** FOR INTERACTION WITH HYDROPHILIC ORGANIC COMPOUNDS
         BB    = PSI + PAR1 + ORGANION
         CC    =-PAR1*(PSI+KAPA)
         LAMDA = 0.5*(-BB+SQRT(BB*BB-4*CC))
C
C *** SAVE CONCENTRATIONS IN MOLAL ARRAY *******************************
C
C *** FOR INTERACTION WITH HYDROPHILIC ORGANIC COMPOUNDS
         MOLAL(1) = PSI+LAMDA + ORGANION         ! HI
         MOLAL(3) = KAPA                         ! NH4I
         MOLAL(5) = LAMDA                        ! SO4I
         MOLAL(6) = MAX (ZERO, PSI+KAPA-LAMDA)   ! HSO4I
         CNH4HS4  = MAX(W(3)-MOLAL(3), ZERO)     ! Solid NH4HSO4
         CH2SO4   = MAX(PSI, ZERO)               ! Free H2SO4
         CALL SSH_CALCMR                             ! Water content
C
C *** CALCULATE ACTIVITIES OR TERMINATE INTERNAL LOOP *****************
C
         IF (FRST.AND.CALAOU .OR. .NOT.FRST.AND.CALAIN) THEN
            CALL SSH_CALCACT     
         ELSE
            GOTO 30
         ENDIF
20    CONTINUE
C
C *** CALCULATE ZERO FUNCTION *******************************************
C
CCC30    FUNCC1= (NH4I*HSO4I/PAR2) - ONE
30    SSH_FUNCC1= (MOLAL(3)*MOLAL(6)/PAR2) - ONE
      RETURN
C
C *** END OF FUNCTION FUNCC1 ********************************************
C
      END

C=======================================================================
C
C *** ISORROPIA CODE
C *** SUBROUTINE SSH_CALCD3
C *** CASE D3
C
C     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
C     1. SULFATE POOR (SULRAT > 2.0)
C     2. THERE IS OLNY A LIQUID PHASE
C
C *** COPYRIGHT 1996-2006, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY,
C *** GEORGIA INSTITUTE OF TECHNOLOGY
C *** WRITTEN BY ATHANASIOS NENES
C *** UPDATED BY CHRISTOS FOUNTOUKIS
C
C=======================================================================
C
      SUBROUTINE SSH_CALCD3
      INCLUDE 'isrpia.inc'
      COMMON /SSHSOLUT/ CHI1, CHI2, CHI3, CHI4, CHI5, CHI6, CHI7, CHI8,
     &               PSI1, PSI2, PSI3, PSI4, PSI5, PSI6, PSI7, PSI8,
     &               A1,   A2,   A3,   A4,   A5,   A6,   A7,   A8


C
C *** FIND DRY COMPOSITION **********************************************
C
      CALL SSH_CALCD1A
C
C *** SETUP PARAMETERS ************************************************
C
      CHI1 = CNH4NO3               ! Save from CALCD1 run
      CHI2 = CNH42S4
      CHI3 = GHNO3
      CHI4 = GNH3
C
      PSI1 = CNH4NO3               ! ASSIGN INITIAL PSI's
      PSI2 = CHI2
      PSI3 = ZERO   
      PSI4 = ZERO  
C
      MOLAL(5) = ZERO
      MOLAL(6) = ZERO
      MOLAL(3) = PSI1
      MOLAL(7) = PSI1
      CALL SSH_CALCMR                  ! Initial water
C
      CALAOU = .TRUE.              ! Outer loop activity calculation flag
      PSI4LO = TINY                ! Low  limit
      PSI4HI = CHI4                ! High limit
C
C *** INITIAL VALUES FOR BISECTION ************************************
C
60    X1 = PSI4LO
      Y1 = SSH_FUNCD3 (X1)
      IF (ABS(Y1).LE.EPS) RETURN
      YLO= Y1                 ! Save Y-value at HI position
C
C *** ROOT TRACKING ; FOR THE RANGE OF HI AND LO **********************
C
      DX = (PSI4HI-PSI4LO)/FLOAT(NDIV)
      DO 10 I=1,NDIV
         X2 = X1+DX
         Y2 = SSH_FUNCD3 (X2)
         IF (SIGN(1.d0,Y1)*SIGN(1.d0,Y2).LT.ZERO) GOTO 20  ! (Y1*Y2.LT.ZERO)
         X1 = X2
         Y1 = Y2
10    CONTINUE
C
C *** NO SUBDIVISION WITH SOLUTION FOUND 
C
      YHI= Y1                      ! Save Y-value at Hi position
      IF (ABS(Y2) .LT. EPS) THEN   ! X2 IS A SOLUTION 
         RETURN
C
C *** { YLO, YHI } < 0.0 THE SOLUTION IS ALWAYS UNDERSATURATED WITH NH3
C Physically I dont know when this might happen, but I have put this
C branch in for completeness. I assume there is no solution; all NO3 goes to the
C gas phase.
C
      ELSE IF (YLO.LT.ZERO .AND. YHI.LT.ZERO) THEN
         P4 = TINY ! PSI4LO ! CHI4
         YY = SSH_FUNCD3(P4)
         GOTO 50
C
C *** { YLO, YHI } > 0.0 THE SOLUTION IS ALWAYS SUPERSATURATED WITH NH3
C This happens when Sul.Rat. = 2.0, so some NH4+ from sulfate evaporates
C and goes to the gas phase ; so I redefine the LO and HI limits of PSI4
C and proceed again with root tracking.
C
      ELSE IF (YLO.GT.ZERO .AND. YHI.GT.ZERO) THEN
80     PSI4HI = PSI4LO
         PSI4LO = PSI4LO - 0.1*(PSI1+PSI2) ! No solution; some NH3 evaporates
         IF (PSI4LO.LT.-(PSI1+PSI2)) THEN
            CALL SSH_PUSHERR (0001, 'CALCD3')  ! WARNING ERROR: NO SOLUTION
            RETURN
         ELSE
            MOLAL(5) = ZERO
            MOLAL(6) = ZERO
            MOLAL(3) = PSI1
            MOLAL(7) = PSI1
            CALL SSH_CALCMR                  ! Initial water
            GOTO 60                        ! Redo root tracking
         ENDIF
      ENDIF
C
C *** PERFORM BISECTION ***********************************************
C
20    DO 30 I=1,MAXIT
         X3 = 0.5*(X1+X2)
         Y3 = SSH_FUNCD3 (X3)
         IF (SIGN(1.d0,Y1)*SIGN(1.d0,Y3) .LE. ZERO) THEN  ! (Y1*Y3 .LE. ZERO)
            Y2    = Y3
            X2    = X3
         ELSE
            Y1    = Y3
            X1    = X3
         ENDIF
         IF (ABS(X2-X1) .LE. EPS*ABS(X1)) GOTO 40
30    CONTINUE
      CALL SSH_PUSHERR (0002, 'CALCD3')    ! WARNING ERROR: NO CONVERGENCE
C
C *** CONVERGED ; RETURN **********************************************
C
40    X3 = 0.5*(X1+X2)
      Y3 = SSH_FUNCD3 (X3)
C 
C *** CALCULATE HSO4 SPECIATION AND RETURN *******************************
C
50    CONTINUE
      IF (MOLAL(1).GT.TINY) THEN
         CALL SSH_CALCHS4 (MOLAL(1), MOLAL(5), ZERO, DELTA)
c         MOLAL(1) = MOLAL(1) - DELTA                     ! H+   EFFECT
         MOLAL(1) = MAX(MOLAL(1) - DELTA,ZERO) ! YK
c      write(*,*) "isofwd.f: L1800 MOLAL(1): ",molal(1)
         MOLAL(5) = MAX(MOLAL(5) - DELTA,ZERO)                     ! SO4  EFFECT
         MOLAL(6) = DELTA                                ! HSO4 EFFECT
      ELSE
c	GOTO 80  !Shupeng ZHU 03122014
      ENDIF
      RETURN
C
C *** END OF SUBROUTINE SSH_CALCD3 ******************************************
C
      END



C=======================================================================
C
C *** ISORROPIA CODE
C *** FUNCTION SSH_FUNCD3
C *** CASE D3 
C     FUNCTION THAT SOLVES THE SYSTEM OF EQUATIONS FOR CASE D3 ; 
C     AND RETURNS THE VALUE OF THE ZEROED FUNCTION IN SSH_FUNCD3.
C
C=======================================================================
C
      DOUBLE PRECISION FUNCTION SSH_FUNCD3 (P4)
      INCLUDE 'isrpia.inc'
      COMMON /SSHSOLUT/ CHI1, CHI2, CHI3, CHI4, CHI5, CHI6, CHI7, CHI8,
     &               PSI1, PSI2, PSI3, PSI4, PSI5, PSI6, PSI7, PSI8,
     &               A1,   A2,   A3,   A4,   A5,   A6,   A7,   A8


C
C *** SETUP PARAMETERS ************************************************
C
      FRST   = .TRUE.
      CALAIN = .TRUE.
      PSI4   = P4
C
C *** SOLVE EQUATIONS ; WITH ITERATIONS FOR ACTIVITY COEF. ************
C
      DO 10 I=1,NSWEEP
         A2   = XK7*(WATER/GAMA(4))**3.0
         A3   = XK4*R*TEMP*(WATER/GAMA(10))**2.0
         A4   = (XK2/XKW)*R*TEMP*(GAMA(10)/GAMA(5))**2.0
         A7   = XKW *RH*WATER*WATER
C
         PSI3 = A3*A4*CHI3*(CHI4-PSI4) - PSI1*(2.D0*PSI2+PSI1+PSI4)
         PSI3 = PSI3/(A3*A4*(CHI4-PSI4) + 2.D0*PSI2+PSI1+PSI4) 
         PSI3 = MIN(MAX(PSI3, ZERO), CHI3)
C
C *** FOR INTERACTION WITH HYDROPHILIC ORGANIC COMPOUNDS
         BB   = PSI4 - PSI3 - ORGANION
CCCOLD         AHI  = 0.5*(-BB + SQRT(BB*BB + 4.d0*A7)) ! This is correct also
CCC         AHI  =2.0*A7/(BB+SQRT(BB*BB + 4.d0*A7)) ! Avoid overflow when HI->0
         DENM = BB+SQRT(BB*BB + 4.d0*A7)
         IF (DENM.LE.TINY) THEN       ! Avoid overflow when HI->0
            ABB  = ABS(BB)
            DENM = (BB+ABB) + 2.0*A7/ABB ! Taylor expansion of SQRT
         ENDIF
         AHI = 2.0*A7/DENM
C
C *** SPECIATION & WATER CONTENT ***************************************
C
         MOLAL (1) = AHI                             ! HI
         MOLAL (3) = PSI1 + PSI4 + 2.D0*PSI2         ! NH4I
         MOLAL (5) = PSI2                            ! SO4I
         MOLAL (6) = ZERO                            ! HSO4I
         MOLAL (7) = PSI3 + PSI1                     ! NO3I
         CNH42S4   = CHI2 - PSI2                     ! Solid (NH4)2SO4
         CNH4NO3   = ZERO                            ! Solid NH4NO3
         GHNO3     = CHI3 - PSI3                     ! Gas HNO3
         GNH3      = CHI4 - PSI4                     ! Gas NH3
         CALL SSH_CALCMR                                 ! Water content
C
C *** CALCULATE ACTIVITIES OR TERMINATE INTERNAL LOOP *****************
C
         IF (FRST.AND.CALAOU .OR. .NOT.FRST.AND.CALAIN) THEN
            CALL SSH_CALCACT     
         ELSE
            GOTO 20
         ENDIF
10    CONTINUE
C
C *** CALCULATE OBJECTIVE FUNCTION ************************************
C
20    CONTINUE
CCC      SSH_FUNCD3= NH4I/HI/MAX(GNH3,TINY)/A4 - ONE 
      SSH_FUNCD3= MOLAL(3)/MOLAL(1)/MAX(GNH3,TINY)/A4 - ONE 
      RETURN
C
C *** END OF FUNCTION SSH_FUNCD3 ********************************************
C
      END
C=======================================================================
C
C *** ISORROPIA CODE
C *** SUBROUTINE SSH_CALCD2
C *** CASE D2
C
C     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
C     1. SULFATE POOR (SULRAT > 2.0)
C     2. THERE IS BOTH A LIQUID & SOLID PHASE
C     3. SOLIDS POSSIBLE : (NH4)2SO4
C
C *** COPYRIGHT 1996-2006, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY,
C *** GEORGIA INSTITUTE OF TECHNOLOGY
C *** WRITTEN BY ATHANASIOS NENES
C *** UPDATED BY CHRISTOS FOUNTOUKIS
C
C=======================================================================
C
      SUBROUTINE SSH_CALCD2
      INCLUDE 'isrpia.inc'
      COMMON /SSHSOLUT/ CHI1, CHI2, CHI3, CHI4, CHI5, CHI6, CHI7, CHI8,
     &               PSI1, PSI2, PSI3, PSI4, PSI5, PSI6, PSI7, PSI8,
     &               A1,   A2,   A3,   A4,   A5,   A6,   A7,   A8


C
C *** FIND DRY COMPOSITION **********************************************
C
      CALL SSH_CALCD1A
C
C *** SETUP PARAMETERS ************************************************
C
      CHI1 = CNH4NO3               ! Save from CALCD1 run
      CHI2 = CNH42S4
      CHI3 = GHNO3
      CHI4 = GNH3
C
      PSI1 = CNH4NO3               ! ASSIGN INITIAL PSI's
      PSI2 = CNH42S4
      PSI3 = ZERO   
      PSI4 = ZERO  
C
      MOLAL(5) = ZERO
      MOLAL(6) = ZERO
      MOLAL(3) = PSI1
      MOLAL(7) = PSI1
      CALL SSH_CALCMR                  ! Initial water
C
      CALAOU = .TRUE.              ! Outer loop activity calculation flag
      PSI4LO = TINY                ! Low  limit
      PSI4HI = CHI4                ! High limit
C
C *** INITIAL VALUES FOR BISECTION ************************************
C
60    X1 = PSI4LO
      Y1 = SSH_FUNCD2 (X1)
      IF (ABS(Y1).LE.EPS) RETURN
      YLO= Y1                 ! Save Y-value at HI position
C
C *** ROOT TRACKING ; FOR THE RANGE OF HI AND LO **********************
C
      DX   = (PSI4HI-PSI4LO)/FLOAT(NDIV)
      DO 10 I=1,NDIV
         X2 = X1+DX
         Y2 = SSH_FUNCD2 (X2)
         IF (SIGN(1.d0,Y1)*SIGN(1.d0,Y2).LT.ZERO) THEN
C
C This is done, in case if Y(PSI4LO)>0, but Y(PSI4LO+DX) < 0 (i.e.undersat)
C
             IF (Y1 .LE. Y2) GOTO 20  ! (Y1*Y2.LT.ZERO)
         ENDIF
         X1 = X2
         Y1 = Y2
10    CONTINUE
C
C *** NO SUBDIVISION WITH SOLUTION FOUND 
C
      YHI= Y1                      ! Save Y-value at Hi position
      IF (ABS(Y2) .LT. EPS) THEN   ! X2 IS A SOLUTION 
         RETURN
C
C *** { YLO, YHI } < 0.0 THE SOLUTION IS ALWAYS UNDERSATURATED WITH NH3
C Physically I dont know when this might happen, but I have put this
C branch in for completeness. I assume there is no solution; all NO3 goes to the
C gas phase.
C
      ELSE IF (YLO.LT.ZERO .AND. YHI.LT.ZERO) THEN
         P4 = TINY ! PSI4LO ! CHI4
         YY = SSH_FUNCD2(P4)
         GOTO 50
C
C *** { YLO, YHI } > 0.0 THE SOLUTION IS ALWAYS SUPERSATURATED WITH NH3
C This happens when Sul.Rat. = 2.0, so some NH4+ from sulfate evaporates
C and goes to the gas phase ; so I redefine the LO and HI limits of PSI4
C and proceed again with root tracking.
C
      ELSE IF (YLO.GT.ZERO .AND. YHI.GT.ZERO) THEN
         PSI4HI = PSI4LO
         PSI4LO = PSI4LO - 0.1*(PSI1+PSI2) ! No solution; some NH3 evaporates
         IF (PSI4LO.LT.-(PSI1+PSI2)) THEN
            CALL SSH_PUSHERR (0001, 'CALCD2')  ! WARNING ERROR: NO SOLUTION
            RETURN
         ELSE
            MOLAL(5) = ZERO
            MOLAL(6) = ZERO
            MOLAL(3) = PSI1
            MOLAL(7) = PSI1
            CALL SSH_CALCMR                  ! Initial water
            GOTO 60                        ! Redo root tracking
         ENDIF
      ENDIF
C
C *** PERFORM BISECTION ***********************************************
C
20    DO 30 I=1,MAXIT
         X3 = 0.5*(X1+X2)
         Y3 = SSH_FUNCD2 (X3)
         IF (SIGN(1.d0,Y1)*SIGN(1.d0,Y3) .LE. ZERO) THEN  ! (Y1*Y3 .LE. ZERO)
            Y2    = Y3
            X2    = X3
         ELSE
            Y1    = Y3
            X1    = X3
         ENDIF
         IF (ABS(X2-X1) .LE. EPS*ABS(X1)) GOTO 40
30    CONTINUE
      CALL SSH_PUSHERR (0002, 'CALCD2')    ! WARNING ERROR: NO CONVERGENCE
C
C *** CONVERGED ; RETURN **********************************************
C
40    X3 = MIN(X1,X2)   ! 0.5*(X1+X2)  ! Get "low" side, it's acidic soln.
      Y3 = SSH_FUNCD2 (X3)
C 
C *** CALCULATE HSO4 SPECIATION AND RETURN *******************************
C
50    CONTINUE
      IF (MOLAL(1).GT.TINY) THEN
         CALL SSH_CALCHS4 (MOLAL(1), MOLAL(5), ZERO, DELTA)
c         MOLAL(1) = MOLAL(1) - DELTA                     ! H+   EFFECT
         MOLAL(1) = MAX(MOLAL(1) - DELTA,ZERO) ! YK 
c      write(*,*) "isofwd.f: L2027 MOLAL(1): ",molal(1)
c         MOLAL(5) = MOLAL(5) - DELTA                     ! SO4  EFFECT
         MOLAL(5) = MAX(MOLAL(5) - DELTA,ZERO) ! YK
         MOLAL(6) = DELTA                                ! HSO4 EFFECT
      ENDIF
      RETURN
C
C *** END OF SUBROUTINE SSH_CALCD2 ******************************************
C
      END



C=======================================================================
C
C *** ISORROPIA CODE
C *** FUNCTION SSH_FUNCD2
C *** CASE D2 
C     FUNCTION THAT SOLVES THE SYSTEM OF EQUATIONS FOR CASE D2 ; 
C     AND RETURNS THE VALUE OF THE ZEROED FUNCTION IN SSH_FUNCD2.
C
C=======================================================================
C
      DOUBLE PRECISION FUNCTION SSH_FUNCD2 (P4)
      INCLUDE 'isrpia.inc'
      COMMON /SSHSOLUT/ CHI1, CHI2, CHI3, CHI4, CHI5, CHI6, CHI7, CHI8,
     &               PSI1, PSI2, PSI3, PSI4, PSI5, PSI6, PSI7, PSI8,
     &               A1,   A2,   A3,   A4,   A5,   A6,   A7,   A8


C
C *** SETUP PARAMETERS ************************************************
C
      CALL SSH_RSTGAM       ! Reset activity coefficients to 0.1
      FRST   = .TRUE.
      CALAIN = .TRUE.
      PSI4   = P4
      PSI2   = CHI2
C
C *** SOLVE EQUATIONS ; WITH ITERATIONS FOR ACTIVITY COEF. ************
C
      DO 10 I=1,NSWEEP
         A2  = XK7*(WATER/GAMA(4))**3.0
         A3  = XK4*R*TEMP*(WATER/GAMA(10))**2.0
         A4  = (XK2/XKW)*R*TEMP*(GAMA(10)/GAMA(5))**2.0
         A7  = XKW *RH*WATER*WATER
C
         IF (CHI2.GT.TINY .AND. WATER.GT.TINY) THEN
            PSI14 = PSI1+PSI4
            CALL SSH_POLY3 (PSI14,0.25*PSI14**2.,-A2/4.D0, PSI2, ISLV)  ! PSI2
            IF (ISLV.EQ.0) THEN
                PSI2 = MIN (PSI2, CHI2)
            ELSE
                PSI2 = ZERO
            ENDIF
         ENDIF
C
         PSI3  = A3*A4*CHI3*(CHI4-PSI4) - PSI1*(2.D0*PSI2+PSI1+PSI4)
         PSI3  = PSI3/(A3*A4*(CHI4-PSI4) + 2.D0*PSI2+PSI1+PSI4) 
ccc         PSI3  = MIN(MAX(PSI3, ZERO), CHI3)
C
C *** FOR INTERACTION WITH HYDROPHILIC ORGANIC COMPOUNDS
         BB   = PSI4-PSI3 - ORGANION ! (BB > 0, acidic solution, <0 alkaline)
C
C Do not change computation scheme for H+, all others did not work well.
C
         DENM = BB+SQRT(BB*BB + 4.d0*A7)
         IF (DENM.LE.TINY) THEN       ! Avoid overflow when HI->0
            ABB  = ABS(BB)
            DENM = (BB+ABB) + 2.d0*A7/ABB ! Taylor expansion of SQRT
         ENDIF
         AHI = 2.d0*A7/DENM
C
C *** SPECIATION & WATER CONTENT ***************************************
C
         MOLAL (1) = AHI                              ! HI
         MOLAL (3) = PSI1 + PSI4 + 2.D0*PSI2          ! NH4
         MOLAL (5) = PSI2                             ! SO4
         MOLAL (6) = ZERO                             ! HSO4
         MOLAL (7) = PSI3 + PSI1                      ! NO3
         CNH42S4   = CHI2 - PSI2                      ! Solid (NH4)2SO4
         CNH4NO3   = ZERO                             ! Solid NH4NO3
         GHNO3     = CHI3 - PSI3                      ! Gas HNO3
         GNH3      = CHI4 - PSI4                      ! Gas NH3
         CALL SSH_CALCMR                                  ! Water content
C
C *** CALCULATE ACTIVITIES OR TERMINATE INTERNAL LOOP *****************
C
         IF (FRST.AND.CALAOU .OR. .NOT.FRST.AND.CALAIN) THEN
            CALL SSH_CALCACT     
         ELSE
            GOTO 20
         ENDIF
10    CONTINUE
C
C *** CALCULATE OBJECTIVE FUNCTION ************************************
C
20    CONTINUE
CCC      SSH_FUNCD2= NH4I/HI/MAX(GNH3,TINY)/A4 - ONE 
      SSH_FUNCD2= MOLAL(3)/MOLAL(1)/MAX(GNH3,TINY)/A4 - ONE 
      RETURN
C
C *** END OF FUNCTION SSH_FUNCD2 ********************************************
C
      END
C=======================================================================
C
C *** ISORROPIA CODE
C *** SUBROUTINE SSH_CALCD1
C *** CASE D1 
C
C     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
C     1. SULFATE POOR (SULRAT > 2.0)
C     2. SOLID AEROSOL ONLY
C     3. SOLIDS POSSIBLE : (NH4)2SO4, NH4NO3
C
C     THERE ARE TWO REGIMES DEFINED BY RELATIVE HUMIDITY:
C     1. RH < MDRH ; ONLY SOLID PHASE POSSIBLE (SUBROUTINE SSH_CALCD1A)
C     2. RH >= MDRH ; LIQUID PHASE POSSIBLE (MDRH REGION)
C
C *** COPYRIGHT 1996-2006, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY,
C *** GEORGIA INSTITUTE OF TECHNOLOGY
C *** WRITTEN BY ATHANASIOS NENES
C *** UPDATED BY CHRISTOS FOUNTOUKIS
C
C=======================================================================
C
      SUBROUTINE SSH_CALCD1
      INCLUDE 'isrpia.inc'
      EXTERNAL SSH_CALCD1A, SSH_CALCD2
C
C *** REGIME DEPENDS UPON THE AMBIENT RELATIVE HUMIDITY *****************
C
      IF (RH.LT.DRMASAN) THEN    
         SCASE = 'D1 ; SUBCASE 1'   ! SOLID PHASE ONLY POSSIBLE
         CALL SSH_CALCD1A            
         SCASE = 'D1 ; SUBCASE 1'
      ELSE
         SCASE = 'D1 ; SUBCASE 2'   ! LIQUID & SOLID PHASE POSSIBLE
         CALL SSH_CALCMDRH (RH, DRMASAN, DRNH4NO3, SSH_CALCD1A,
     c        SSH_CALCD2)
         SCASE = 'D1 ; SUBCASE 2'
      ENDIF
C 
      RETURN
C
C *** END OF SUBROUTINE SSH_CALCD1 ******************************************
C
      END


C=======================================================================
C
C *** ISORROPIA CODE
C *** SUBROUTINE SSH_CALCD1A
C *** CASE D1 ; SUBCASE 1
C
C     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
C     1. SULFATE POOR (SULRAT > 2.0)
C     2. SOLID AEROSOL ONLY
C     3. SOLIDS POSSIBLE : (NH4)2SO4, NH4NO3
C
C     THE SOLID (NH4)2SO4 IS CALCULATED FROM THE SULFATES, WHILE NH4NO3
C     IS CALCULATED FROM NH3-HNO3 EQUILIBRIUM. 'ZE' IS THE AMOUNT OF
C     NH4NO3 THAT VOLATIZES WHEN ALL POSSILBE NH4NO3 IS INITIALLY IN
C     THE SOLID PHASE.
C
C *** COPYRIGHT 1996-2006, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY,
C *** GEORGIA INSTITUTE OF TECHNOLOGY
C *** WRITTEN BY ATHANASIOS NENES
C *** UPDATED BY CHRISTOS FOUNTOUKIS
C
C=======================================================================
C
      SUBROUTINE SSH_CALCD1A
      INCLUDE 'isrpia.inc'
C
C *** SETUP PARAMETERS ************************************************
C
      PARM    = XK10/(R*TEMP)/(R*TEMP)
C
C *** CALCULATE NH4NO3 THAT VOLATIZES *********************************
C
      CNH42S4 = W(2)                                    
      X       = MAX(ZERO, MIN(W(3)-2.0*CNH42S4, W(4)))  ! MAX NH4NO3
      PS      = MAX(W(3) - X - 2.0*CNH42S4, ZERO)
      OM      = MAX(W(4) - X, ZERO)
C
      OMPS    = OM+PS
      DIAK    = SQRT(OMPS*OMPS + 4.0*PARM)              ! DIAKRINOUSA
      ZE      = MIN(X, 0.5*(-OMPS + DIAK))              ! THETIKI RIZA
C
C *** SPECIATION *******************************************************
C
      CNH4NO3 = X  - ZE    ! Solid NH4NO3
      GNH3    = PS + ZE    ! Gas NH3
      GHNO3   = OM + ZE    ! Gas HNO3
C
      RETURN
C
C *** END OF SUBROUTINE SSH_CALCD1A *****************************************
C
      END
C=======================================================================
C
C *** ISORROPIA CODE
C *** SUBROUTINE SSH_CALCG5
C *** CASE G5
C
C     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
C     1. SULFATE POOR (SULRAT > 2.0) ; SODIUM POOR (SODRAT < 2.0)
C     2. THERE IS BOTH A LIQUID & SOLID PHASE
C     3. SOLIDS POSSIBLE : (NH4)2SO4, NH4CL, NA2SO4
C
C *** COPYRIGHT 1996-2006, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY,
C *** GEORGIA INSTITUTE OF TECHNOLOGY
C *** WRITTEN BY ATHANASIOS NENES
C *** UPDATED BY CHRISTOS FOUNTOUKIS
C
C=======================================================================
C
      SUBROUTINE SSH_CALCG5
      INCLUDE 'isrpia.inc'
C
      DOUBLE PRECISION LAMDA
      COMMON /SSHCASEG/ CHI1, CHI2, CHI3, CHI4, CHI5, CHI6, LAMDA,
     &               PSI1, PSI2, PSI3, PSI4, PSI5, PSI6, PSI7,
     &               A1,   A2,   A3,   A4,   A5,   A6,   A7


C
C *** SETUP PARAMETERS ************************************************
C
      CALAOU = .TRUE.   
      CHI1   = 0.5*W(1)
      CHI2   = MAX (W(2)-CHI1, ZERO)
      CHI3   = ZERO
      CHI4   = MAX (W(3)-2.D0*CHI2, ZERO)
      CHI5   = W(4)
      CHI6   = W(5)
C 
      PSI1   = CHI1
      PSI2   = CHI2
      PSI6LO = TINY                  
      PSI6HI = CHI6-TINY    ! MIN(CHI6-TINY, CHI4)
C
      WATER  = CHI2/M0(4) + CHI1/M0(2)
C
C *** INITIAL VALUES FOR BISECTION ************************************
C
      X1 = PSI6LO
      Y1 = SSH_FUNCG5A (X1)
      IF (CHI6.LE.TINY) GOTO 50  
ccc      IF (ABS(Y1).LE.EPS .OR. CHI6.LE.TINY) GOTO 50  
ccc      IF (WATER .LE. TINY) RETURN                    ! No water
C
C *** ROOT TRACKING ; FOR THE RANGE OF HI AND LO **********************
C
      DX = (PSI6HI-PSI6LO)/FLOAT(NDIV)
      DO 10 I=1,NDIV
         X2 = X1+DX 
         Y2 = SSH_FUNCG5A (X2)
         IF (SIGN(1.d0,Y1)*SIGN(1.d0,Y2).LT.ZERO) GOTO 20  ! (Y1*Y2.LT.ZERO)
         X1 = X2
         Y1 = Y2
10    CONTINUE
C
C *** NO SUBDIVISION WITH SOLUTION; IF ABS(Y2)<EPS SOLUTION IS ASSUMED
C
      IF (ABS(Y2) .GT. EPS) Y2 = SSH_FUNCG5A (PSI6LO)
      GOTO 50
C
C *** PERFORM BISECTION ***********************************************
C
20    DO 30 I=1,MAXIT
         X3 = 0.5*(X1+X2)
         Y3 = SSH_FUNCG5A (X3)
         IF (SIGN(1.d0,Y1)*SIGN(1.d0,Y3) .LE. ZERO) THEN  ! (Y1*Y3 .LE. ZERO)
            Y2    = Y3
            X2    = X3
         ELSE
            Y1    = Y3
            X1    = X3
         ENDIF
         IF (ABS(X2-X1) .LE. EPS*X1) GOTO 40
30    CONTINUE
      CALL SSH_PUSHERR (0002, 'CALCG5')    ! WARNING ERROR: NO CONVERGENCE
C
C *** CONVERGED ; RETURN **********************************************
C
40    X3 = 0.5*(X1+X2)
      Y3 = SSH_FUNCG5A (X3)
C 
C *** CALCULATE HSO4 SPECIATION AND RETURN *******************************
C
50    CONTINUE
      IF (MOLAL(1).GT.TINY .AND. MOLAL(5).GT.TINY) THEN  ! If quadrat.called
         CALL SSH_CALCHS4 (MOLAL(1), MOLAL(5), ZERO, DELTA)
c         MOLAL(1) = MOLAL(1) - DELTA                    ! H+   EFFECT
         MOLAL(1) = MAX(MOLAL(1) - DELTA,ZERO) ! YK
c      write(*,*) "isofwd.f: L2321 MOLAL(1): ",molal(1)
c         MOLAL(5) = MOLAL(5) - DELTA                    ! SO4  EFFECT
         MOLAL(5) = MAX(MOLAL(5) - DELTA,ZERO) ! YK    
         MOLAL(6) = DELTA                               ! HSO4 EFFECT
      ENDIF
C
      RETURN
C
C *** END OF SUBROUTINE SSH_CALCG5 *******************************************
C
      END




C=======================================================================
C
C *** ISORROPIA CODE
C *** SUBROUTINE SSH_FUNCG5A
C *** CASE G5
C
C     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
C     1. SULFATE POOR (SULRAT > 2.0) ; SODIUM POOR (SODRAT < 2.0)
C     2. THERE IS BOTH A LIQUID & SOLID PHASE
C     3. SOLIDS POSSIBLE : (NH4)2SO4, NH4CL, NA2SO4
C
C *** COPYRIGHT 1996-2006, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY,
C *** GEORGIA INSTITUTE OF TECHNOLOGY
C *** WRITTEN BY ATHANASIOS NENES
C *** UPDATED BY CHRISTOS FOUNTOUKIS
C
C=======================================================================
C
      DOUBLE PRECISION FUNCTION SSH_FUNCG5A (X)
      INCLUDE 'isrpia.inc'
C
      DOUBLE PRECISION LAMDA
      COMMON /SSHCASEG/ CHI1, CHI2, CHI3, CHI4, CHI5, CHI6, LAMDA,
     &               PSI1, PSI2, PSI3, PSI4, PSI5, PSI6, PSI7,
     &               A1,   A2,   A3,   A4,   A5,   A6,   A7


C
C *** SETUP PARAMETERS ************************************************
C
      PSI6   = X
      FRST   = .TRUE.
      CALAIN = .TRUE. 
C
C *** SOLVE EQUATIONS ; WITH ITERATIONS FOR ACTIVITY COEF. ************
C
      DO 10 I=1,NSWEEP
C
      A1  = XK5 *(WATER/GAMA(2))**3.0
      A2  = XK7 *(WATER/GAMA(4))**3.0
      A4  = (XK2/XKW)*R*TEMP*(GAMA(10)/GAMA(5))**2.0
      A5  = XK4 *R*TEMP*(WATER/GAMA(10))**2.0
      A6  = XK3 *R*TEMP*(WATER/GAMA(11))**2.0
      AKK = A4*A6
C
C  CALCULATE DISSOCIATION QUANTITIES
C
      IF (CHI5.GE.TINY) THEN
         PSI5 = PSI6*CHI5/(A6/A5*(CHI6-PSI6) + PSI6)
      ELSE
         PSI5 = TINY
      ENDIF
C
CCC      IF(CHI4.GT.TINY) THEN
      IF(W(2).GT.TINY) THEN       ! Accounts for NH3 evaporation
         BB   =-(CHI4 + PSI6 + PSI5 + 1.d0/A4)
         CC   = CHI4*(PSI5+PSI6) - 2.d0*PSI2/A4
         DD   = MAX(BB*BB-4.d0*CC,ZERO)           ! Patch proposed by Uma Shankar, 19/11/01
         PSI4 =0.5d0*(-BB - SQRT(DD))
      ELSE
         PSI4 = TINY
      ENDIF
C
C *** CALCULATE SPECIATION ********************************************
C
      MOLAL (2) = 2.0D0*PSI1                          ! NAI
      MOLAL (3) = 2.0*PSI2 + PSI4                     ! NH4I
      MOLAL (4) = PSI6                                ! CLI
      MOLAL (5) = PSI2 + PSI1                         ! SO4I
      MOLAL (6) = ZERO
      MOLAL (7) = PSI5                                ! NO3I
C
C *** FOR INTERACTION WITH HYDROPHILIC ORGANIC COMPOUNDS
      SMIN      = 2.d0*MOLAL(5)+MOLAL(7)+MOLAL(4)-MOLAL(2)-MOLAL(3)
     &     + ORGANION
      CALL SSH_CALCPH (SMIN, HI, OHI)
      MOLAL (1) = HI
C 
      GNH3      = MAX(CHI4 - PSI4, TINY)              ! Gas NH3
      GHNO3     = MAX(CHI5 - PSI5, TINY)              ! Gas HNO3
      GHCL      = MAX(CHI6 - PSI6, TINY)              ! Gas HCl
C
      CNH42S4   = ZERO                                ! Solid (NH4)2SO4
      CNH4NO3   = ZERO                                ! Solid NH4NO3
      CNH4CL    = ZERO                                ! Solid NH4Cl
C
      CALL SSH_CALCMR                                     ! Water content
C
C *** CALCULATE ACTIVITIES OR TERMINATE INTERNAL LOOP *****************
C
      IF (FRST.AND.CALAOU .OR. .NOT.FRST.AND.CALAIN) THEN
         CALL SSH_CALCACT     
      ELSE
         GOTO 20
      ENDIF
10    CONTINUE
C
C *** CALCULATE FUNCTION VALUE FOR OUTER LOOP ***************************
C
20    SSH_FUNCG5A = MOLAL(1)*MOLAL(4)/GHCL/A6 - ONE
CCC         SSH_FUNCG5A = MOLAL(3)*MOLAL(4)/GHCL/GNH3/A6/A4 - ONE
C
      RETURN
C
C *** END OF FUNCTION SSH_FUNCG5A *******************************************
C
      END

C=======================================================================
C
C *** ISORROPIA CODE
C *** SUBROUTINE SSH_CALCG4
C *** CASE G4
C
C     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
C     1. SULFATE POOR (SULRAT > 2.0) ; SODIUM POOR (SODRAT < 2.0)
C     2. THERE IS BOTH A LIQUID & SOLID PHASE
C     3. SOLIDS POSSIBLE : (NH4)2SO4, NH4CL, NA2SO4
C
C *** COPYRIGHT 1996-2006, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY,
C *** GEORGIA INSTITUTE OF TECHNOLOGY
C *** WRITTEN BY ATHANASIOS NENES
C *** UPDATED BY CHRISTOS FOUNTOUKIS
C
C=======================================================================
C
      SUBROUTINE SSH_CALCG4
      INCLUDE 'isrpia.inc'
C
      DOUBLE PRECISION LAMDA
      COMMON /SSHCASEG/ CHI1, CHI2, CHI3, CHI4, CHI5, CHI6, LAMDA,
     &               PSI1, PSI2, PSI3, PSI4, PSI5, PSI6, PSI7,
     &               A1,   A2,   A3,   A4,   A5,   A6,   A7


C
C *** SETUP PARAMETERS ************************************************
C
      CALAOU = .TRUE.   
      CHI1   = 0.5*W(1)
      CHI2   = MAX (W(2)-CHI1, ZERO)
      CHI3   = ZERO
      CHI4   = MAX (W(3)-2.D0*CHI2, ZERO)
      CHI5   = W(4)
      CHI6   = W(5)
C 
      PSI2   = CHI2
      PSI6LO = TINY                  
      PSI6HI = CHI6-TINY    ! MIN(CHI6-TINY, CHI4)
C
      WATER  = CHI2/M0(4) + CHI1/M0(2)
C
C *** INITIAL VALUES FOR BISECTION ************************************
C
      X1 = PSI6LO
      Y1 = SSH_FUNCG4A (X1)
      IF (CHI6.LE.TINY) GOTO 50  
CCC      IF (ABS(Y1).LE.EPS .OR. CHI6.LE.TINY .OR. WATER .LE. TINY) GOTO 50
CCC      IF (WATER .LE. TINY) RETURN                    ! No water
C
C *** ROOT TRACKING ; FOR THE RANGE OF HI AND LO **********************
C
      DX = (PSI6HI-PSI6LO)/FLOAT(NDIV)
      DO 10 I=1,NDIV
         X2  = X1+DX
         Y2  = SSH_FUNCG4A (X2)
         IF (SIGN(1.d0,Y1)*SIGN(1.d0,Y2).LT.ZERO) GOTO 20  ! (Y1*Y2.LT.ZERO)
         X1  = X2
         Y1  = Y2
10    CONTINUE
C
C *** NO SUBDIVISION WITH SOLUTION; IF ABS(Y2)<EPS SOLUTION IS ASSUMED
C
      IF (ABS(Y2) .GT. EPS) Y2 = SSH_FUNCG4A (PSI6LO)
      GOTO 50
C
C *** PERFORM BISECTION ***********************************************
C
20    DO 30 I=1,MAXIT
         X3 = 0.5*(X1+X2)
         Y3 = SSH_FUNCG4A (X3)
         IF (SIGN(1.d0,Y1)*SIGN(1.d0,Y3) .LE. ZERO) THEN  ! (Y1*Y3 .LE. ZERO)
            Y2    = Y3
            X2    = X3
         ELSE
            Y1    = Y3
            X1    = X3
         ENDIF
         IF (ABS(X2-X1) .LE. EPS*X1) GOTO 40
30    CONTINUE
      CALL SSH_PUSHERR (0002, 'CALCG4')    ! WARNING ERROR: NO CONVERGENCE
C
C *** CONVERGED ; RETURN **********************************************
C
40    X3 = 0.5*(X1+X2)
      Y3 = SSH_FUNCG4A (X3)
C 
C *** CALCULATE HSO4 SPECIATION AND RETURN *******************************
C
50    CONTINUE
      IF (MOLAL(1).GT.TINY .AND. MOLAL(5).GT.TINY) THEN
         CALL SSH_CALCHS4 (MOLAL(1), MOLAL(5), ZERO, DELTA)
c         MOLAL(1) = MOLAL(1) - DELTA                     ! H+   EFFECT
         MOLAL(1) = MAX(MOLAL(1) - DELTA,ZERO) ! YK 
c      write(*,*) "isofwd.f: L2535 MOLAL(1): ",molal(1)
c         MOLAL(5) = MOLAL(5) - DELTA                     ! SO4  EFFECT
         MOLAL(5) = MAX(MOLAL(5) - DELTA,ZERO) ! YK
         MOLAL(6) = DELTA                                ! HSO4 EFFECT
      ENDIF
C
      RETURN
C
C *** END OF SUBROUTINE SSH_CALCG4 *******************************************
C
      END




C=======================================================================
C
C *** ISORROPIA CODE
C *** SUBROUTINE SSH_FUNCG4A
C *** CASE G4
C
C     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
C     1. SULFATE POOR (SULRAT > 2.0) ; SODIUM POOR (SODRAT < 2.0)
C     2. THERE IS BOTH A LIQUID & SOLID PHASE
C     3. SOLIDS POSSIBLE : (NH4)2SO4, NH4CL, NA2SO4
C
C *** COPYRIGHT 1996-2006, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY,
C *** GEORGIA INSTITUTE OF TECHNOLOGY
C *** WRITTEN BY ATHANASIOS NENES
C *** UPDATED BY CHRISTOS FOUNTOUKIS
C
C=======================================================================
C
      DOUBLE PRECISION FUNCTION SSH_FUNCG4A (X)
      INCLUDE 'isrpia.inc'
C
      DOUBLE PRECISION LAMDA, NAI, NH4I, NO3I
      COMMON /SSHCASEG/ CHI1, CHI2, CHI3, CHI4, CHI5, CHI6, LAMDA,
     &               PSI1, PSI2, PSI3, PSI4, PSI5, PSI6, PSI7,
     &               A1,   A2,   A3,   A4,   A5,   A6,   A7


C
C *** SETUP PARAMETERS ************************************************
C
      PSI6   = X
      PSI1   = CHI1
      FRST   = .TRUE.
      CALAIN = .TRUE. 
C
C *** SOLVE EQUATIONS ; WITH ITERATIONS FOR ACTIVITY COEF. ************
C
      DO 10 I=1,NSWEEP
C
      A1  = XK5 *(WATER/GAMA(2))**3.0
      A2  = XK7 *(WATER/GAMA(4))**3.0
      A4  = (XK2/XKW)*R*TEMP*(GAMA(10)/GAMA(5))**2.0
      A5  = XK4 *R*TEMP*(WATER/GAMA(10))**2.0
      A6  = XK3 *R*TEMP*(WATER/GAMA(11))**2.0
C
C  CALCULATE DISSOCIATION QUANTITIES
C
      IF (CHI5.GE.TINY) THEN
         PSI5 = PSI6*CHI5/(A6/A5*(CHI6-PSI6) + PSI6)
      ELSE
         PSI5 = TINY
      ENDIF
C
CCC      IF(CHI4.GT.TINY) THEN
      IF(W(2).GT.TINY) THEN       ! Accounts for NH3 evaporation
         BB   =-(CHI4 + PSI6 + PSI5 + 1.d0/A4)
         CC   = CHI4*(PSI5+PSI6) - 2.d0*PSI2/A4
         DD   = MAX(BB*BB-4.d0*CC,ZERO) ! Patch proposed by Uma shankar, 19/11/2001
         PSI4 =0.5d0*(-BB - SQRT(DD))
      ELSE
         PSI4 = TINY
      ENDIF
C
C  CALCULATE CONCENTRATIONS
C
      NH4I = 2.0*PSI2 + PSI4
      CLI  = PSI6
      SO4I = PSI2 + PSI1
      NO3I = PSI5
      NAI  = 2.0D0*PSI1  
C
      CALL SSH_CALCPH(2.d0*SO4I+NO3I+CLI-NAI-NH4I, HI, OHI)
C
C *** Na2SO4 DISSOLUTION
C
      IF (CHI1.GT.TINY .AND. WATER.GT.TINY) THEN        ! PSI1
         CALL SSH_POLY3 (PSI2, ZERO, -A1/4.D0, PSI1, ISLV)
         IF (ISLV.EQ.0) THEN
             PSI1 = MIN (PSI1, CHI1)
         ELSE
             PSI1 = ZERO
         ENDIF
      ELSE
         PSI1 = ZERO
      ENDIF
C
C *** SAVE CONCENTRATIONS IN MOLAL ARRAY ******************************
C
      MOLAL (1) = HI
      MOLAL (2) = NAI
      MOLAL (3) = NH4I
      MOLAL (4) = CLI
      MOLAL (5) = SO4I
      MOLAL (6) = ZERO
      MOLAL (7) = NO3I
C 
C *** CALCULATE GAS / SOLID SPECIES (LIQUID IN MOLAL ALREADY) *********
C
      GNH3      = MAX(CHI4 - PSI4, TINY)
      GHNO3     = MAX(CHI5 - PSI5, TINY)
      GHCL      = MAX(CHI6 - PSI6, TINY)
C
      CNH42S4   = ZERO
      CNH4NO3   = ZERO
      CNH4CL    = ZERO
      CNA2SO4   = MAX(CHI1-PSI1,ZERO)
C
C *** CALCULATE MOLALR ARRAY, WATER AND ACTIVITIES **********************
C
      CALL SSH_CALCMR
C
C *** CALCULATE ACTIVITIES OR TERMINATE INTERNAL LOOP *****************
C
      IF (FRST.AND.CALAOU .OR. .NOT.FRST.AND.CALAIN) THEN
         CALL SSH_CALCACT     
      ELSE
         GOTO 20
      ENDIF
10    CONTINUE
C
C *** CALCULATE FUNCTION VALUE FOR OUTER LOOP ***************************
C
20    SSH_FUNCG4A = MOLAL(1)*MOLAL(4)/GHCL/A6 - ONE
CCC         SSH_FUNCG4A = MOLAL(3)*MOLAL(4)/GHCL/GNH3/A6/A4 - ONE
C
      RETURN
C
C *** END OF FUNCTION SSH_FUNCG4A *******************************************
C
      END

C=======================================================================
C
C *** ISORROPIA CODE
C *** SUBROUTINE SSH_CALCG3
C *** CASE G3
C
C     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
C     1. SULFATE POOR (SULRAT > 2.0) ; SODIUM POOR (SODRAT < 2.0)
C     2. LIQUID & SOLID PHASE ARE BOTH POSSIBLE
C     3. SOLIDS POSSIBLE : (NH4)2SO4, NH4CL, NA2SO4
C
C *** COPYRIGHT 1996-2006, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY,
C *** GEORGIA INSTITUTE OF TECHNOLOGY
C *** WRITTEN BY ATHANASIOS NENES
C *** UPDATED BY CHRISTOS FOUNTOUKIS
C
C=======================================================================
C
      SUBROUTINE SSH_CALCG3
      INCLUDE 'isrpia.inc'
      EXTERNAL SSH_CALCG1A, SSH_CALCG4
C
C *** REGIME DEPENDS ON THE EXISTANCE OF WATER AND OF THE RH ************
C
      IF (W(4).GT.TINY .AND. W(5).GT.TINY) THEN ! NO3,CL EXIST, WATER POSSIBLE
         SCASE = 'G3 ; SUBCASE 1'  
         CALL SSH_CALCG3A
         SCASE = 'G3 ; SUBCASE 1' 
      ELSE                                      ! NO3, CL NON EXISTANT
         SCASE = 'G1 ; SUBCASE 1'  
         CALL SSH_CALCG1A
         SCASE = 'G1 ; SUBCASE 1'  
      ENDIF
C
      IF (WATER.LE.TINY) THEN
         IF (RH.LT.DRMG3) THEN        ! ONLY SOLIDS 
            WATER = TINY
            DO 10 I=1,NIONS
               MOLAL(I) = ZERO
10          CONTINUE
            CALL SSH_CALCG1A
            SCASE = 'G3 ; SUBCASE 2'  
            RETURN
         ELSE
            SCASE = 'G3 ; SUBCASE 3'  ! MDRH REGION (NA2SO4, NH42S4)  
            CALL SSH_CALCMDRH (RH, DRMG3, DRNH42S4, SSH_CALCG1A,
     &           SSH_CALCG4)
            SCASE = 'G3 ; SUBCASE 3'  
         ENDIF
      ENDIF
C 
      RETURN
C
C *** END OF SUBROUTINE SSH_CALCG3 ******************************************
C
      END


C=======================================================================
C
C *** ISORROPIA CODE
C *** SUBROUTINE SSH_CALCG3A
C *** CASE G3 ; SUBCASE 1
C
C     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
C     1. SULFATE POOR (SULRAT > 2.0) ; SODIUM POOR (SODRAT < 2.0)
C     2. THERE IS BOTH A LIQUID & SOLID PHASE
C     3. SOLIDS POSSIBLE : (NH4)2SO4, NH4CL, NA2SO4
C
C *** COPYRIGHT 1996-2006, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY,
C *** GEORGIA INSTITUTE OF TECHNOLOGY
C *** WRITTEN BY ATHANASIOS NENES
C *** UPDATED BY CHRISTOS FOUNTOUKIS
C
C=======================================================================
C
      SUBROUTINE SSH_CALCG3A
      INCLUDE 'isrpia.inc'
C
      DOUBLE PRECISION LAMDA
      COMMON /SSHCASEG/ CHI1, CHI2, CHI3, CHI4, CHI5, CHI6, LAMDA,
     &               PSI1, PSI2, PSI3, PSI4, PSI5, PSI6, PSI7,
     &               A1,   A2,   A3,   A4,   A5,   A6,   A7


C
C *** SETUP PARAMETERS ************************************************
C
      CALAOU = .TRUE.   
      CHI1   = 0.5*W(1)
      CHI2   = MAX (W(2)-CHI1, ZERO)
      CHI3   = ZERO
      CHI4   = MAX (W(3)-2.D0*CHI2, ZERO)
      CHI5   = W(4)
      CHI6   = W(5)
C 
      PSI6LO = TINY                  
      PSI6HI = CHI6-TINY    ! MIN(CHI6-TINY, CHI4)
C
      WATER  = TINY
C
C *** INITIAL VALUES FOR BISECTION ************************************
C
      X1 = PSI6LO
      Y1 = SSH_FUNCG3A (X1)
      IF (CHI6.LE.TINY) GOTO 50  
CCC      IF (ABS(Y1).LE.EPS .OR. CHI6.LE.TINY .OR. WATER .LE. TINY) GOTO 50
CCC      IF (WATER .LE. TINY) RETURN                    ! No water
C
C *** ROOT TRACKING ; FOR THE RANGE OF HI AND LO **********************
C
      DX = (PSI6HI-PSI6LO)/FLOAT(NDIV)
      DO 10 I=1,NDIV
         X2  = X1+DX 
         Y2  = SSH_FUNCG3A (X2)
C
         IF (SIGN(1.d0,Y1)*SIGN(1.d0,Y2).LT.ZERO) GOTO 20  ! (Y1*Y2.LT.ZERO)
         X1  = X2
         Y1  = Y2
10    CONTINUE
C
C *** NO SUBDIVISION WITH SOLUTION; IF ABS(Y2)<EPS SOLUTION IS ASSUMED
C
      IF (ABS(Y2) .GT. EPS) Y2 = SSH_FUNCG3A (PSI6LO)
      GOTO 50
C
C *** PERFORM BISECTION ***********************************************
C
20    DO 30 I=1,MAXIT
         X3 = 0.5*(X1+X2)
         Y3 = SSH_FUNCG3A (X3)
         IF (SIGN(1.d0,Y1)*SIGN(1.d0,Y3) .LE. ZERO) THEN  ! (Y1*Y3 .LE. ZERO)
            Y2    = Y3
            X2    = X3
         ELSE
            Y1    = Y3
            X1    = X3
         ENDIF
         IF (ABS(X2-X1) .LE. EPS*X1) GOTO 40
30    CONTINUE
      CALL SSH_PUSHERR (0002, 'CALCG3A')    ! WARNING ERROR: NO CONVERGENCE
C
C *** CONVERGED ; RETURN **********************************************
C
40    X3 = 0.5*(X1+X2)
      Y3 = SSH_FUNCG3A (X3)
C 
C *** FINAL CALCULATIONS *************************************************
C
50    CONTINUE
C
C *** Na2SO4 DISSOLUTION
C
      IF (CHI1.GT.TINY .AND. WATER.GT.TINY) THEN        ! PSI1
         CALL SSH_POLY3 (PSI2, ZERO, -A1/4.D0, PSI1, ISLV)
         IF (ISLV.EQ.0) THEN
             PSI1 = MIN (PSI1, CHI1)
         ELSE
             PSI1 = ZERO
         ENDIF
      ELSE
         PSI1 = ZERO
      ENDIF
      MOLAL(2) = 2.0D0*PSI1               ! Na+  EFFECT
      MOLAL(5) = MOLAL(5) + PSI1          ! SO4  EFFECT
      CNA2SO4  = MAX(CHI1 - PSI1, ZERO)   ! NA2SO4(s) depletion
C
C *** HSO4 equilibrium
C 
      IF (MOLAL(1).GT.TINY .AND. MOLAL(5).GT.TINY) THEN
         CALL SSH_CALCHS4 (MOLAL(1), MOLAL(5), ZERO, DELTA)
c         MOLAL(1) = MOLAL(1) - DELTA                     ! H+   EFFECT
         MOLAL(1) = MAX(MOLAL(1) - DELTA,ZERO) ! YK
c      write(*,*) "isofwd.f: L2848 MOLAL(1): ",molal(1)
c         MOLAL(5) = MOLAL(5) - DELTA                     ! SO4  EFFECT
         MOLAL(5) = MAX(MOLAL(5) - DELTA,ZERO) ! YK 
         MOLAL(6) = DELTA                                ! HSO4 EFFECT
      ENDIF
C
      RETURN
C
C *** END OF SUBROUTINE SSH_CALCG3A ******************************************
C
      END




C=======================================================================
C
C *** ISORROPIA CODE
C *** SUBROUTINE SSH_FUNCG3A
C *** CASE G3 ; SUBCASE 1
C
C     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
C     1. SULFATE POOR (SULRAT > 2.0) ; SODIUM POOR (SODRAT < 2.0)
C     2. THERE IS BOTH A LIQUID & SOLID PHASE
C     3. SOLIDS POSSIBLE : (NH4)2SO4, NH4CL, NA2SO4
C
C *** COPYRIGHT 1996-2006, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY,
C *** GEORGIA INSTITUTE OF TECHNOLOGY
C *** WRITTEN BY ATHANASIOS NENES
C *** UPDATED BY CHRISTOS FOUNTOUKIS
C
C=======================================================================
C
      DOUBLE PRECISION FUNCTION SSH_FUNCG3A (X)
      INCLUDE 'isrpia.inc'
C
      DOUBLE PRECISION LAMDA
      COMMON /SSHCASEG/ CHI1, CHI2, CHI3, CHI4, CHI5, CHI6, LAMDA,
     &               PSI1, PSI2, PSI3, PSI4, PSI5, PSI6, PSI7,
     &               A1,   A2,   A3,   A4,   A5,   A6,   A7


C
C *** SETUP PARAMETERS ************************************************
C
      PSI6   = X
      PSI2   = CHI2
      FRST   = .TRUE.
      CALAIN = .TRUE. 
C
C *** SOLVE EQUATIONS ; WITH ITERATIONS FOR ACTIVITY COEF. ************
C
      DO 10 I=1,NSWEEP
C
      A1  = XK5 *(WATER/GAMA(2))**3.0
      A2  = XK7 *(WATER/GAMA(4))**3.0
      A4  = (XK2/XKW)*R*TEMP*(GAMA(10)/GAMA(5))**2.0
      A5  = XK4 *R*TEMP*(WATER/GAMA(10))**2.0
      A6  = XK3 *R*TEMP*(WATER/GAMA(11))**2.0
C
C  CALCULATE DISSOCIATION QUANTITIES
C
      IF (CHI5.GE.TINY) THEN
         PSI5 = PSI6*CHI5/(A6/A5*(CHI6-PSI6) + PSI6)
      ELSE
         PSI5 = TINY
      ENDIF
C
CCC      IF(CHI4.GT.TINY) THEN
      IF(W(2).GT.TINY) THEN       ! Accounts for NH3 evaporation
         BB   =-(CHI4 + PSI6 + PSI5 + 1.d0/A4)
         CC   = CHI4*(PSI5+PSI6) - 2.d0*PSI2/A4
         DD   = MAX(BB*BB-4.d0*CC,ZERO)  ! Patch proposed by Uma Shankar, 19/11/01
         PSI4 =0.5d0*(-BB - SQRT(DD))
      ELSE
         PSI4 = TINY
      ENDIF
C
      IF (CHI2.GT.TINY .AND. WATER.GT.TINY) THEN     
         CALL SSH_POLY3 (PSI4, PSI4*PSI4/4.D0, -A2/4.D0, PSI20, ISLV)
         IF (ISLV.EQ.0) PSI2 = MIN (PSI20, CHI2)
      ENDIF
C 
C *** CALCULATE GAS / SOLID SPECIES (LIQUID IN MOLAL ALREADY) *********
C
      MOLAL (2) = ZERO                                ! Na
      MOLAL (3) = 2.0*PSI2 + PSI4                     ! NH4I
      MOLAL (4) = PSI6                                ! CLI
      MOLAL (5) = PSI2                                ! SO4I
      MOLAL (6) = ZERO                                ! HSO4
      MOLAL (7) = PSI5                                ! NO3I
C
C *** FOR INTERACTION WITH HYDROPHILIC ORGANIC COMPOUNDS
      SMIN      = 2.d0*MOLAL(5)+MOLAL(7)+MOLAL(4)-MOLAL(2)-MOLAL(3)
     &     + ORGANION
      CALL SSH_CALCPH (SMIN, HI, OHI)
      MOLAL (1) = HI
c
      GNH3      = MAX(CHI4 - PSI4, TINY)              ! Gas NH3
      GHNO3     = MAX(CHI5 - PSI5, TINY)              ! Gas HNO3
      GHCL      = MAX(CHI6 - PSI6, TINY)              ! Gas HCl
C
      CNH42S4   = CHI2 - PSI2                         ! Solid (NH4)2SO4
      CNH4NO3   = ZERO                                ! Solid NH4NO3
      CNH4CL    = ZERO                                ! Solid NH4Cl
C
      CALL SSH_CALCMR                                     ! Water content
C
C *** CALCULATE ACTIVITIES OR TERMINATE INTERNAL LOOP *****************
C
      IF (FRST.AND.CALAOU .OR. .NOT.FRST.AND.CALAIN) THEN
         CALL SSH_CALCACT     
      ELSE
         GOTO 20
      ENDIF
10    CONTINUE
C
C *** CALCULATE FUNCTION VALUE FOR OUTER LOOP ***************************
C
20    SSH_FUNCG3A = MOLAL(1)*MOLAL(4)/GHCL/A6 - ONE
CCC         SSH_FUNCG3A = MOLAL(3)*MOLAL(4)/GHCL/GNH3/A6/A4 - ONE
C
      RETURN
C
C *** END OF FUNCTION SSH_FUNCG3A *******************************************
C
      END

C=======================================================================
C
C *** ISORROPIA CODE
C *** SUBROUTINE SSH_CALCG2
C *** CASE G2
C
C     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
C     1. SULFATE POOR (SULRAT > 2.0) ; SODIUM POOR (SODRAT < 2.0)
C     2. LIQUID & SOLID PHASE ARE BOTH POSSIBLE
C     3. SOLIDS POSSIBLE : (NH4)2SO4, NH4CL, NA2SO4
C
C *** COPYRIGHT 1996-2006, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY,
C *** GEORGIA INSTITUTE OF TECHNOLOGY
C *** WRITTEN BY ATHANASIOS NENES
C *** UPDATED BY CHRISTOS FOUNTOUKIS
C
C=======================================================================
C
      SUBROUTINE SSH_CALCG2
      INCLUDE 'isrpia.inc'
      EXTERNAL SSH_CALCG1A, SSH_CALCG3A, SSH_CALCG4
C
C *** REGIME DEPENDS ON THE EXISTANCE OF NITRATES ***********************
C
      IF (W(4).GT.TINY) THEN        ! NO3 EXISTS, WATER POSSIBLE
         SCASE = 'G2 ; SUBCASE 1'  
         CALL SSH_CALCG2A
         SCASE = 'G2 ; SUBCASE 1' 
      ELSE                          ! NO3 NON EXISTANT, WATER NOT POSSIBLE
         SCASE = 'G1 ; SUBCASE 1'  
         CALL SSH_CALCG1A
         SCASE = 'G1 ; SUBCASE 1'  
      ENDIF
C
C *** REGIME DEPENDS ON THE EXISTANCE OF WATER AND OF THE RH ************
C
      IF (WATER.LE.TINY) THEN
         IF (RH.LT.DRMG2) THEN             ! ONLY SOLIDS 
            WATER = TINY
            DO 10 I=1,NIONS
               MOLAL(I) = ZERO
10          CONTINUE
            CALL SSH_CALCG1A
            SCASE = 'G2 ; SUBCASE 2'  
         ELSE
            IF (W(5).GT. TINY) THEN
               SCASE = 'G2 ; SUBCASE 3'    ! MDRH (NH4CL, NA2SO4, NH42S4)  
               CALL SSH_CALCMDRH (RH, DRMG2, DRNH4CL, SSH_CALCG1A,
     &              SSH_CALCG3A)
               SCASE = 'G2 ; SUBCASE 3'  
            ENDIF
            IF (WATER.LE.TINY .AND. RH.GE.DRMG3) THEN
               SCASE = 'G2 ; SUBCASE 4'    ! MDRH (NA2SO4, NH42S4)
               CALL SSH_CALCMDRH (RH, DRMG3, DRNH42S4, SSH_CALCG1A,
     &              SSH_CALCG4)
               SCASE = 'G2 ; SUBCASE 4'  
            ELSE
               WATER = TINY
               DO 20 I=1,NIONS
                  MOLAL(I) = ZERO
20             CONTINUE
               CALL SSH_CALCG1A
               SCASE = 'G2 ; SUBCASE 2'  
            ENDIF
         ENDIF
      ENDIF
C 
      RETURN
C
C *** END OF SUBROUTINE SSH_CALCG2 ******************************************
C
      END


C=======================================================================
C
C *** ISORROPIA CODE
C *** SUBROUTINE SSH_CALCG2A
C *** CASE G2 ; SUBCASE 1
C
C     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
C     1. SULFATE POOR (SULRAT > 2.0) ; SODIUM POOR (SODRAT < 2.0)
C     2. THERE IS BOTH A LIQUID & SOLID PHASE
C     3. SOLIDS POSSIBLE : (NH4)2SO4, NH4CL, NA2SO4
C
C *** COPYRIGHT 1996-2006, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY,
C *** GEORGIA INSTITUTE OF TECHNOLOGY
C *** WRITTEN BY ATHANASIOS NENES
C *** UPDATED BY CHRISTOS FOUNTOUKIS
C
C=======================================================================
C
      SUBROUTINE SSH_CALCG2A
      INCLUDE 'isrpia.inc'
C
      DOUBLE PRECISION LAMDA
      COMMON /SSHCASEG/ CHI1, CHI2, CHI3, CHI4, CHI5, CHI6, LAMDA,
     &               PSI1, PSI2, PSI3, PSI4, PSI5, PSI6, PSI7,
     &               A1,   A2,   A3,   A4,   A5,   A6,   A7


C
C *** SETUP PARAMETERS ************************************************
C
      CALAOU = .TRUE.   
      CHI1   = 0.5*W(1)
      CHI2   = MAX (W(2)-CHI1, ZERO)
      CHI3   = ZERO
      CHI4   = MAX (W(3)-2.D0*CHI2, ZERO)
      CHI5   = W(4)
      CHI6   = W(5)
C 
      PSI6LO = TINY                  
      PSI6HI = CHI6-TINY
C
      WATER  = TINY
C
C *** INITIAL VALUES FOR BISECTION ************************************
C
      X1 = PSI6LO
      Y1 = SSH_FUNCG2A (X1)
      IF (CHI6.LE.TINY) GOTO 50  
CCC      IF (ABS(Y1).LE.EPS .OR. CHI6.LE.TINY) GOTO 50  
CCC      IF (WATER .LE. TINY) GOTO 50               ! No water
C
C *** ROOT TRACKING ; FOR THE RANGE OF HI AND LO **********************
C
      DX = (PSI6HI-PSI6LO)/FLOAT(NDIV)
      DO 10 I=1,NDIV
         X2 = X1+DX 
         Y2 = SSH_FUNCG2A (X2)
         IF (SIGN(1.d0,Y1)*SIGN(1.d0,Y2).LT.ZERO) GOTO 20  ! (Y1*Y2.LT.ZERO)
         X1 = X2
         Y1 = Y2
10    CONTINUE
C
C *** NO SUBDIVISION WITH SOLUTION; IF ABS(Y2)<EPS SOLUTION IS ASSUMED
C
      IF (ABS(Y2) .GT. EPS) WATER = TINY
      GOTO 50
C
C *** PERFORM BISECTION ***********************************************
C
20    DO 30 I=1,MAXIT
         X3 = 0.5*(X1+X2)
         Y3 = SSH_FUNCG2A (X3)
         IF (SIGN(1.d0,Y1)*SIGN(1.d0,Y3) .LE. ZERO) THEN  ! (Y1*Y3 .LE. ZERO)
            Y2    = Y3
            X2    = X3
         ELSE
            Y1    = Y3
            X1    = X3
         ENDIF
         IF (ABS(X2-X1) .LE. EPS*X1) GOTO 40
30    CONTINUE
      CALL SSH_PUSHERR (0002, 'CALCG2A')    ! WARNING ERROR: NO CONVERGENCE
C
C *** CONVERGED ; RETURN **********************************************
C
40    X3 = 0.5*(X1+X2)
      IF (X3.LE.TINY2) THEN   ! PRACTICALLY NO NITRATES, SO DRY SOLUTION
         WATER = TINY
      ELSE
         Y3 = SSH_FUNCG2A (X3)
      ENDIF
C 
C *** CALCULATE HSO4 SPECIATION AND RETURN *******************************
C
50    CONTINUE
C
C *** Na2SO4 DISSOLUTION
C
      IF (CHI1.GT.TINY .AND. WATER.GT.TINY) THEN        ! PSI1
         CALL SSH_POLY3 (PSI2, ZERO, -A1/4.D0, PSI1, ISLV)
         IF (ISLV.EQ.0) THEN
             PSI1 = MIN (PSI1, CHI1)
         ELSE
             PSI1 = ZERO
         ENDIF
      ELSE
         PSI1 = ZERO
      ENDIF
      MOLAL(2) = 2.0D0*PSI1               ! Na+  EFFECT
      MOLAL(5) = MOLAL(5) + PSI1          ! SO4  EFFECT
      CNA2SO4  = MAX(CHI1 - PSI1, ZERO)   ! NA2SO4(s) depletion
C
C *** HSO4 equilibrium
C 
      IF (MOLAL(1).GT.TINY .AND. MOLAL(5).GT.TINY) THEN
         CALL SSH_CALCHS4 (MOLAL(1), MOLAL(5), ZERO, DELTA)
c         MOLAL(1) = MOLAL(1) - DELTA     ! H+   AFFECT
         MOLAL(1) = MAX(MOLAL(1) - DELTA,ZERO) ! YK
c      write(*,*) "isofwd.f: L3161 MOLAL(1): ",molal(1)
c         MOLAL(5) = MOLAL(5) - DELTA     ! SO4  AFFECT
         MOLAL(5) = MAX(MOLAL(5) - DELTA,ZERO) ! YK 
         MOLAL(6) = DELTA                ! HSO4 AFFECT
      ENDIF
C
      RETURN
C
C *** END OF SUBROUTINE SSH_CALCG2A ******************************************
C
      END




C=======================================================================
C
C *** ISORROPIA CODE
C *** SUBROUTINE SSH_FUNCG2A
C *** CASE G2
C
C     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
C     1. SULFATE POOR (SULRAT > 2.0) ; SODIUM POOR (SODRAT < 2.0)
C     2. THERE IS BOTH A LIQUID & SOLID PHASE
C     3. SOLIDS POSSIBLE : (NH4)2SO4, NH4CL, NA2SO4
C
C *** COPYRIGHT 1996-2006, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY,
C *** GEORGIA INSTITUTE OF TECHNOLOGY
C *** WRITTEN BY ATHANASIOS NENES
C *** UPDATED BY CHRISTOS FOUNTOUKIS
C
C=======================================================================
C
      DOUBLE PRECISION FUNCTION SSH_FUNCG2A (X)
      INCLUDE 'isrpia.inc'
C
      DOUBLE PRECISION LAMDA
      COMMON /SSHCASEG/ CHI1, CHI2, CHI3, CHI4, CHI5, CHI6, LAMDA,
     &               PSI1, PSI2, PSI3, PSI4, PSI5, PSI6, PSI7,
     &               A1,   A2,   A3,   A4,   A5,   A6,   A7


C
C *** SETUP PARAMETERS ************************************************
C
      PSI6   = X
      PSI2   = CHI2
      PSI3   = ZERO
      FRST   = .TRUE.
      CALAIN = .TRUE. 
C
C *** SOLVE EQUATIONS ; WITH ITERATIONS FOR ACTIVITY COEF. ************
C
      DO 10 I=1,NSWEEP
C
      A1  = XK5 *(WATER/GAMA(2))**3.0
      A2  = XK7 *(WATER/GAMA(4))**3.0
      A4  = (XK2/XKW)*R*TEMP*(GAMA(10)/GAMA(5))**2.0
      A5  = XK4 *R*TEMP*(WATER/GAMA(10))**2.0
      A6  = XK3 *R*TEMP*(WATER/GAMA(11))**2.0
C
      DENO = MAX(CHI6-PSI6-PSI3, ZERO)
      PSI5 = CHI5/((A6/A5)*(DENO/PSI6) + ONE)
C
      PSI4 = MIN(PSI5+PSI6,CHI4)
C
      IF (CHI2.GT.TINY .AND. WATER.GT.TINY) THEN     
         CALL SSH_POLY3 (PSI4, PSI4*PSI4/4.D0, -A2/4.D0, PSI20, ISLV)
         IF (ISLV.EQ.0) PSI2 = MIN (PSI20, CHI2)
      ENDIF
C
C *** SAVE CONCENTRATIONS IN MOLAL ARRAY ******************************
C
      MOLAL (2) = ZERO                             ! NA
      MOLAL (3) = 2.0*PSI2 + PSI4                  ! NH4I
      MOLAL (4) = PSI6                             ! CLI
      MOLAL (5) = PSI2                             ! SO4I
      MOLAL (6) = ZERO                             ! HSO4
      MOLAL (7) = PSI5                             ! NO3I
C
CCC      MOLAL (1) = MAX(CHI5 - PSI5, TINY)*A5/PSI5   ! HI
C *** FOR INTERACTION WITH HYDROPHILIC ORGANIC COMPOUNDS
      SMIN      = 2.d0*MOLAL(5)+MOLAL(7)+MOLAL(4)-MOLAL(2)-MOLAL(3)
     &     + ORGANION
      CALL SSH_CALCPH (SMIN, HI, OHI)
      MOLAL (1) = HI
C 
C *** CALCULATE GAS / SOLID SPECIES (LIQUID IN MOLAL ALREADY) *********
C
      GNH3      = MAX(CHI4 - PSI4, TINY)
      GHNO3     = MAX(CHI5 - PSI5, TINY)
      GHCL      = MAX(CHI6 - PSI6, TINY)
C
      CNH42S4   = MAX(CHI2 - PSI2, ZERO)
      CNH4NO3   = ZERO
C      
C *** NH4Cl(s) calculations
C
      A3   = XK6 /(R*TEMP*R*TEMP)
      IF (GNH3*GHCL.GT.A3) THEN
         DELT = MIN(GNH3, GHCL)
         BB = -(GNH3+GHCL)
         CC = GNH3*GHCL-A3
         DD = BB*BB - 4.D0*CC
         PSI31 = 0.5D0*(-BB + SQRT(DD))
         PSI32 = 0.5D0*(-BB - SQRT(DD))
         IF (DELT-PSI31.GT.ZERO .AND. PSI31.GT.ZERO) THEN
            PSI3 = PSI31
         ELSEIF (DELT-PSI32.GT.ZERO .AND. PSI32.GT.ZERO) THEN
            PSI3 = PSI32
         ELSE
            PSI3 = ZERO
         ENDIF
      ELSE
         PSI3 = ZERO
      ENDIF
C 
C *** CALCULATE GAS / SOLID SPECIES (LIQUID IN MOLAL ALREADY) *********
C
      GNH3    = MAX(GNH3 - PSI3, TINY)
      GHCL    = MAX(GHCL - PSI3, TINY)
      CNH4CL  = PSI3
C
C *** CALCULATE MOLALR ARRAY, WATER AND ACTIVITIES **********************
C
      CALL SSH_CALCMR
C
C *** CALCULATE ACTIVITIES OR TERMINATE INTERNAL LOOP *****************
C
      IF (FRST.AND.CALAOU .OR. .NOT.FRST.AND.CALAIN) THEN
         CALL SSH_CALCACT     
      ELSE
         GOTO 20
      ENDIF
10    CONTINUE
C
C *** CALCULATE FUNCTION VALUE FOR OUTER LOOP ***************************
C
20    IF (CHI4.LE.TINY) THEN
         SSH_FUNCG2A = MOLAL(1)*MOLAL(4)/GHCL/A6 - ONE
      ELSE
         SSH_FUNCG2A = MOLAL(3)*MOLAL(4)/GHCL/GNH3/A6/A4 - ONE
      ENDIF
C
      RETURN
C
C *** END OF FUNCTION SSH_FUNCG2A *******************************************
C
      END

C=======================================================================
C
C *** ISORROPIA CODE
C *** SUBROUTINE SSH_CALCG1
C *** CASE G1
C
C     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
C     1. SULFATE POOR (SULRAT > 2.0) ; SODIUM POOR (SODRAT < 2.0)
C     2. SOLID AEROSOL ONLY
C     3. SOLIDS POSSIBLE : (NH4)2SO4, NH4NO3, NH4CL, NA2SO4
C
C     THERE ARE TWO POSSIBLE REGIMES HERE, DEPENDING ON RELATIVE HUMIDITY:
C     1. WHEN RH >= MDRH ; LIQUID PHASE POSSIBLE (MDRH REGION)
C     2. WHEN RH < MDRH  ; ONLY SOLID PHASE POSSIBLE (SUBROUTINE SSH_CALCG1A)
C
C *** COPYRIGHT 1996-2006, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY,
C *** GEORGIA INSTITUTE OF TECHNOLOGY
C *** WRITTEN BY ATHANASIOS NENES
C *** UPDATED BY CHRISTOS FOUNTOUKIS
C
C=======================================================================
C
      SUBROUTINE SSH_CALCG1
      INCLUDE 'isrpia.inc'
      EXTERNAL SSH_CALCG1A, SSH_CALCG2A
C
C *** REGIME DEPENDS UPON THE AMBIENT RELATIVE HUMIDITY *****************
C
      IF (RH.LT.DRMG1) THEN    
         SCASE = 'G1 ; SUBCASE 1'  
         CALL SSH_CALCG1A              ! SOLID PHASE ONLY POSSIBLE
         SCASE = 'G1 ; SUBCASE 1'
      ELSE
         SCASE = 'G1 ; SUBCASE 2'  ! LIQUID & SOLID PHASE POSSIBLE
         CALL SSH_CALCMDRH (RH, DRMG1, DRNH4NO3, SSH_CALCG1A,
     c        SSH_CALCG2A)
         SCASE = 'G1 ; SUBCASE 2'
      ENDIF
C 
      RETURN
C
C *** END OF SUBROUTINE SSH_CALCG1 ******************************************
C
      END


C=======================================================================
C
C *** ISORROPIA CODE
C *** SUBROUTINE SSH_CALCG1A
C *** CASE G1 ; SUBCASE 1
C
C     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
C     1. SULFATE POOR (SULRAT > 2.0)
C     2. SOLID AEROSOL ONLY
C     3. SOLIDS POSSIBLE : (NH4)2SO4, NH4NO3
C
C     SOLID (NH4)2SO4 IS CALCULATED FROM THE SULFATES, WHILE NH4NO3
C     IS CALCULATED FROM NH3-HNO3 EQUILIBRIUM. 'ZE' IS THE AMOUNT OF
C     NH4NO3 THAT VOLATIZES WHEN ALL POSSILBE NH4NO3 IS INITIALLY IN
C     THE SOLID PHASE.
C
C *** COPYRIGHT 1996-2006, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY,
C *** GEORGIA INSTITUTE OF TECHNOLOGY
C *** WRITTEN BY ATHANASIOS NENES
C *** UPDATED BY CHRISTOS FOUNTOUKIS
C
C=======================================================================
C
      SUBROUTINE SSH_CALCG1A
      INCLUDE 'isrpia.inc'
      DOUBLE PRECISION LAMDA, LAMDA1, LAMDA2, KAPA, KAPA1, KAPA2
C
C *** CALCULATE NON VOLATILE SOLIDS ***********************************
C
      CNA2SO4 = 0.5*W(1)
      CNH42S4 = W(2) - CNA2SO4
C
C *** CALCULATE VOLATILE SPECIES **************************************
C
      ALF     = W(3) - 2.0*CNH42S4
      BET     = W(5)
      GAM     = W(4)
C
      RTSQ    = R*TEMP*R*TEMP
      A1      = XK6/RTSQ
      A2      = XK10/RTSQ
C
      THETA1  = GAM - BET*(A2/A1)
      THETA2  = A2/A1
C
C QUADRATIC EQUATION SOLUTION
C
      BB      = (THETA1-ALF-BET*(ONE+THETA2))/(ONE+THETA2)
      CC      = (ALF*BET-A1-BET*THETA1)/(ONE+THETA2)
      DD      = BB*BB - 4.0D0*CC
      IF (DD.LT.ZERO) GOTO 100   ! Solve each reaction seperately
C
C TWO ROOTS FOR KAPA, CHECK AND SEE IF ANY VALID
C
      SQDD    = SQRT(DD)
      KAPA1   = 0.5D0*(-BB+SQDD)
      KAPA2   = 0.5D0*(-BB-SQDD)
      LAMDA1  = THETA1 + THETA2*KAPA1
      LAMDA2  = THETA1 + THETA2*KAPA2
C
      IF (KAPA1.GE.ZERO .AND. LAMDA1.GE.ZERO) THEN
         IF (ALF-KAPA1-LAMDA1.GE.ZERO .AND.
     &       BET-KAPA1.GE.ZERO .AND. GAM-LAMDA1.GE.ZERO) THEN
             KAPA = KAPA1
             LAMDA= LAMDA1
             GOTO 200
         ENDIF
      ENDIF
C
      IF (KAPA2.GE.ZERO .AND. LAMDA2.GE.ZERO) THEN
         IF (ALF-KAPA2-LAMDA2.GE.ZERO .AND. 
     &       BET-KAPA2.GE.ZERO .AND. GAM-LAMDA2.GE.ZERO) THEN
             KAPA = KAPA2
             LAMDA= LAMDA2
             GOTO 200
         ENDIF
      ENDIF
C
C SEPERATE SOLUTION OF NH4CL & NH4NO3 EQUILIBRIA 
C 
100   KAPA  = ZERO
      LAMDA = ZERO
      DD1   = (ALF+BET)*(ALF+BET) - 4.0D0*(ALF*BET-A1)
      DD2   = (ALF+GAM)*(ALF+GAM) - 4.0D0*(ALF*GAM-A2)
C
C NH4CL EQUILIBRIUM
C
      IF (DD1.GE.ZERO) THEN
         SQDD1 = SQRT(DD1)
         KAPA1 = 0.5D0*(ALF+BET + SQDD1)
         KAPA2 = 0.5D0*(ALF+BET - SQDD1)
C
         IF (KAPA1.GE.ZERO .AND. KAPA1.LE.MIN(ALF,BET)) THEN
            KAPA = KAPA1 
         ELSE IF (KAPA2.GE.ZERO .AND. KAPA2.LE.MIN(ALF,BET)) THEN
            KAPA = KAPA2
         ELSE
            KAPA = ZERO
         ENDIF
      ENDIF
C
C NH4NO3 EQUILIBRIUM
C
      IF (DD2.GE.ZERO) THEN
         SQDD2 = SQRT(DD2)
         LAMDA1= 0.5D0*(ALF+GAM + SQDD2)
         LAMDA2= 0.5D0*(ALF+GAM - SQDD2)
C
         IF (LAMDA1.GE.ZERO .AND. LAMDA1.LE.MIN(ALF,GAM)) THEN
            LAMDA = LAMDA1 
         ELSE IF (LAMDA2.GE.ZERO .AND. LAMDA2.LE.MIN(ALF,GAM)) THEN
            LAMDA = LAMDA2
         ELSE
            LAMDA = ZERO
         ENDIF
      ENDIF
C
C IF BOTH KAPA, LAMDA ARE > 0, THEN APPLY EXISTANCE CRITERION
C
      IF (KAPA.GT.ZERO .AND. LAMDA.GT.ZERO) THEN
         IF (BET .LT. LAMDA/THETA1) THEN
            KAPA = ZERO
         ELSE
            LAMDA= ZERO
         ENDIF
      ENDIF
C
C *** CALCULATE COMPOSITION OF VOLATILE SPECIES ***********************
C
200   CONTINUE
      CNH4NO3 = LAMDA
      CNH4CL  = KAPA
C
      GNH3    = MAX(ALF - KAPA - LAMDA, ZERO)
      GHNO3   = MAX(GAM - LAMDA, ZERO)
      GHCL    = MAX(BET - KAPA, ZERO)
C
      RETURN
C
C *** END OF SUBROUTINE SSH_CALCG1A *****************************************
C
      END
C=======================================================================
C
C *** ISORROPIA CODE
C *** SUBROUTINE SSH_CALCH6
C *** CASE H6
C
C     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
C     1. SULFATE POOR (SULRAT > 2.0) ; SODIUM RICH (SODRAT >= 2.0)
C     2. THERE IS BOTH A LIQUID & SOLID PHASE
C     3. SOLIDS POSSIBLE : (NH4)2SO4, NH4CL, NA2SO4
C
C *** COPYRIGHT 1996-2006, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY,
C *** GEORGIA INSTITUTE OF TECHNOLOGY
C *** WRITTEN BY ATHANASIOS NENES
C *** UPDATED BY CHRISTOS FOUNTOUKIS
C
C=======================================================================
C
      SUBROUTINE SSH_CALCH6
      INCLUDE 'isrpia.inc'
C
      COMMON /SSHSOLUT/ CHI1, CHI2, CHI3, CHI4, CHI5, CHI6, CHI7, CHI8,
     &               PSI1, PSI2, PSI3, PSI4, PSI5, PSI6, PSI7, PSI8,
     &               A1,   A2,   A3,   A4,   A5,   A6,   A7,   A8


C
C *** SETUP PARAMETERS ************************************************
C
      CALAOU = .TRUE.   
      CHI1   = W(2)                                ! CNA2SO4
      CHI2   = ZERO                                ! CNH42S4
      CHI3   = ZERO                                ! CNH4CL
      FRNA   = MAX (W(1)-2.D0*CHI1, ZERO)       
      CHI8   = MIN (FRNA, W(4))                    ! CNANO3
      CHI4   = W(3)                                ! NH3(g)
      CHI5   = MAX (W(4)-CHI8, ZERO)               ! HNO3(g)
      CHI7   = MIN (MAX(FRNA-CHI8, ZERO), W(5))    ! CNACL
      CHI6   = MAX (W(5)-CHI7, ZERO)               ! HCL(g)
C
      PSI6LO = TINY                  
      PSI6HI = CHI6-TINY    ! MIN(CHI6-TINY, CHI4)
C
C *** INITIAL VALUES FOR BISECTION ************************************
C
      X1 = PSI6LO
      Y1 = SSH_FUNCH6A (X1)
      IF (ABS(Y1).LE.EPS .OR. CHI6.LE.TINY) GOTO 50  
C
C *** ROOT TRACKING ; FOR THE RANGE OF HI AND LO **********************
C
      DX = (PSI6HI-PSI6LO)/FLOAT(NDIV)
      DO 10 I=1,NDIV
         X2 = X1+DX 
         Y2 = SSH_FUNCH6A (X2)
         IF (SIGN(1.d0,Y1)*SIGN(1.d0,Y2).LT.ZERO) GOTO 20  ! (Y1*Y2.LT.ZERO)
         X1 = X2
         Y1 = Y2
10    CONTINUE
C
C *** NO SUBDIVISION WITH SOLUTION; IF ABS(Y2)<EPS SOLUTION IS ASSUMED
C
      IF (ABS(Y2) .GT. EPS) Y2 = SSH_FUNCH6A (PSI6LO)
      GOTO 50
C
C *** PERFORM BISECTION ***********************************************
C
20    DO 30 I=1,MAXIT
         X3 = 0.5*(X1+X2)
         Y3 = SSH_FUNCH6A (X3)
         IF (SIGN(1.d0,Y1)*SIGN(1.d0,Y3) .LE. ZERO) THEN  ! (Y1*Y3 .LE. ZERO)
            Y2    = Y3
            X2    = X3
         ELSE
            Y1    = Y3
            X1    = X3
         ENDIF
         IF (ABS(X2-X1) .LE. EPS*X1) GOTO 40
30    CONTINUE
      CALL SSH_PUSHERR (0002, 'CALCH6')    ! WARNING ERROR: NO CONVERGENCE
C
C *** CONVERGED ; RETURN **********************************************
C
40    X3 = 0.5*(X1+X2)
      Y3 = SSH_FUNCH6A (X3)
C 
C *** CALCULATE HSO4 SPECIATION AND RETURN *******************************
C
50    CONTINUE
      IF (MOLAL(1).GT.TINY .AND. MOLAL(5).GT.TINY) THEN
         CALL SSH_CALCHS4 (MOLAL(1), MOLAL(5), ZERO, DELTA)
c         MOLAL(1) = MOLAL(1) - DELTA                     ! H+   EFFECT
         MOLAL(1) = MAX(MOLAL(1) - DELTA,ZERO) ! YK   
c      write(*,*) "isofwd.f: L3585 MOLAL(1): ",molal(1)
c         MOLAL(5) = MOLAL(5) - DELTA                     ! SO4  EFFECT
         MOLAL(5) = MAX(MOLAL(5) - DELTA,ZERO) ! YK 
         MOLAL(6) = DELTA                                ! HSO4 EFFECT
      ENDIF
C
      RETURN
C
C *** END OF SUBROUTINE SSH_CALCH6 ******************************************
C
      END




C=======================================================================
C
C *** ISORROPIA CODE
C *** SUBROUTINE SSH_FUNCH6A
C *** CASE H6
C
C     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
C     1. SULFATE POOR (SULRAT > 2.0) ; SODIUM RICH (SODRAT >= 2.0)
C     2. THERE IS BOTH A LIQUID & SOLID PHASE
C     3. SOLIDS POSSIBLE : (NH4)2SO4, NH4CL, NA2SO4
C
C *** COPYRIGHT 1996-2006, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY,
C *** GEORGIA INSTITUTE OF TECHNOLOGY
C *** WRITTEN BY ATHANASIOS NENES
C *** UPDATED BY CHRISTOS FOUNTOUKIS
C
C=======================================================================
C
      DOUBLE PRECISION FUNCTION SSH_FUNCH6A (X)
      INCLUDE 'isrpia.inc'
C
      COMMON /SSHSOLUT/ CHI1, CHI2, CHI3, CHI4, CHI5, CHI6, CHI7, CHI8,
     &               PSI1, PSI2, PSI3, PSI4, PSI5, PSI6, PSI7, PSI8,
     &               A1,   A2,   A3,   A4,   A5,   A6,   A7,   A8


C
C *** SETUP PARAMETERS ************************************************
C
      PSI6   = X
      PSI1   = CHI1
      PSI2   = ZERO
      PSI3   = ZERO
      PSI7   = CHI7
      PSI8   = CHI8 
      FRST   = .TRUE.
      CALAIN = .TRUE. 
C
C *** SOLVE EQUATIONS ; WITH ITERATIONS FOR ACTIVITY COEF. ************
C
      DO 10 I=1,NSWEEP
C
      A1  = XK5 *(WATER/GAMA(2))**3.0
      A4  = (XK2/XKW)*R*TEMP*(GAMA(10)/GAMA(5))**2.0
      A5  = XK4 *R*TEMP*(WATER/GAMA(10))**2.0
      A6  = XK3 *R*TEMP*(WATER/GAMA(11))**2.0
      A7  = XK8 *(WATER/GAMA(1))**2.0
      A8  = XK9 *(WATER/GAMA(3))**2.0
      A9  = XK1*WATER/GAMA(7)*(GAMA(8)/GAMA(7))**2.
C
C  CALCULATE DISSOCIATION QUANTITIES
C
      PSI5 = CHI5*(PSI6+PSI7) - A6/A5*PSI8*(CHI6-PSI6-PSI3)
      PSI5 = PSI5/(A6/A5*(CHI6-PSI6-PSI3) + PSI6 + PSI7)
      PSI5 = MAX(PSI5, TINY)
C
      IF (W(3).GT.TINY .AND. WATER.GT.TINY) THEN  ! First try 3rd order soln
         BB   =-(CHI4 + PSI6 + PSI5 + 1.d0/A4)
         CC   = CHI4*(PSI5+PSI6)
         DD   = BB*BB-4.d0*CC
         PSI4 =0.5d0*(-BB - SQRT(DD))
         PSI4 = MIN(PSI4,CHI4)
      ELSE
         PSI4 = TINY
      ENDIF
C
C *** CALCULATE SPECIATION ********************************************
C
      MOLAL (2) = PSI8 + PSI7 + 2.D0*PSI1               ! NAI
      MOLAL (3) = PSI4                                  ! NH4I
      MOLAL (4) = PSI6 + PSI7                           ! CLI
      MOLAL (5) = PSI2 + PSI1                           ! SO4I
      MOLAL (6) = ZERO                                  ! HSO4I
      MOLAL (7) = PSI5 + PSI8                           ! NO3I
C
C *** FOR INTERACTION WITH HYDROPHILIC ORGANIC COMPOUNDS
      SMIN      = 2.d0*MOLAL(5)+MOLAL(7)+MOLAL(4)-MOLAL(2)-MOLAL(3)
     &     + ORGANION
      CALL SSH_CALCPH (SMIN, HI, OHI)
      MOLAL (1) = HI
C 
      GNH3      = MAX(CHI4 - PSI4, TINY)
      GHNO3     = MAX(CHI5 - PSI5, TINY)
      GHCL      = MAX(CHI6 - PSI6, TINY)
C
      CNH42S4   = ZERO
      CNH4NO3   = ZERO
      CNACL     = MAX(CHI7 - PSI7, ZERO)
      CNANO3    = MAX(CHI8 - PSI8, ZERO)
      CNA2SO4   = MAX(CHI1 - PSI1, ZERO) 
C
      CALL SSH_CALCMR                                    ! Water content
C
C *** CALCULATE ACTIVITIES OR TERMINATE INTERNAL LOOP *****************
C
      IF (FRST.AND.CALAOU .OR. .NOT.FRST.AND.CALAIN) THEN
         CALL SSH_CALCACT     
      ELSE
         GOTO 20
      ENDIF
10    CONTINUE
C
C *** CALCULATE FUNCTION VALUE FOR OUTER LOOP ***************************
C
20    SSH_FUNCH6A = MOLAL(3)*MOLAL(4)/GHCL/GNH3/A6/A4 - ONE
C
      RETURN
C
C *** END OF FUNCTION SSH_FUNCH6A *******************************************
C
      END

C=======================================================================
C
C *** ISORROPIA CODE
C *** SUBROUTINE SSH_CALCH5
C *** CASE H5
C
C     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
C     1. SULFATE POOR (SULRAT > 2.0) ; SODIUM RICH (SODRAT >= 2.0)
C     2. THERE IS BOTH A LIQUID & SOLID PHASE
C     3. SOLIDS POSSIBLE : (NH4)2SO4, NH4CL, NA2SO4
C
C *** COPYRIGHT 1996-2006, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY,
C *** GEORGIA INSTITUTE OF TECHNOLOGY
C *** WRITTEN BY ATHANASIOS NENES
C *** UPDATED BY CHRISTOS FOUNTOUKIS
C
C=======================================================================
C
      SUBROUTINE SSH_CALCH5
      INCLUDE 'isrpia.inc'
C
      COMMON /SSHSOLUT/ CHI1, CHI2, CHI3, CHI4, CHI5, CHI6, CHI7, CHI8,
     &               PSI1, PSI2, PSI3, PSI4, PSI5, PSI6, PSI7, PSI8,
     &               A1,   A2,   A3,   A4,   A5,   A6,   A7,   A8


C
C *** REGIME DEPENDS ON THE EXISTANCE OF NITRATES ***********************
C
      IF (W(4).LE.TINY .AND. W(5).LE.TINY) THEN  
         SCASE = 'H5'  
         CALL SSH_CALCH1A
         SCASE = 'H5'  
         RETURN
      ENDIF
C
C *** SETUP PARAMETERS ************************************************
C
      CALAOU = .TRUE.   
      CHI1   = W(2)                                ! CNA2SO4
      CHI2   = ZERO                                ! CNH42S4
      CHI3   = ZERO                                ! CNH4CL
      FRNA   = MAX (W(1)-2.D0*CHI1, ZERO)       
      CHI8   = MIN (FRNA, W(4))                    ! CNANO3
      CHI4   = W(3)                                ! NH3(g)
      CHI5   = MAX (W(4)-CHI8, ZERO)               ! HNO3(g)
      CHI7   = MIN (MAX(FRNA-CHI8, ZERO), W(5))    ! CNACL
      CHI6   = MAX (W(5)-CHI7, ZERO)               ! HCL(g)
C
      PSI6LO = TINY                  
      PSI6HI = CHI6-TINY    ! MIN(CHI6-TINY, CHI4)
C
C *** INITIAL VALUES FOR BISECTION ************************************
C
      X1 = PSI6LO
      Y1 = SSH_FUNCH5A (X1)
      IF (ABS(Y1).LE.EPS .OR. CHI6.LE.TINY) GOTO 50  
C
C *** ROOT TRACKING ; FOR THE RANGE OF HI AND LO **********************
C
      DX = (PSI6HI-PSI6LO)/FLOAT(NDIV)
      DO 10 I=1,NDIV
         X2 = X1+DX 
         Y2 = SSH_FUNCH5A (X2)
         IF (SIGN(1.d0,Y1)*SIGN(1.d0,Y2).LT.ZERO) GOTO 20  ! (Y1*Y2.LT.ZERO)
         X1 = X2
         Y1 = Y2
10    CONTINUE
C
C *** NO SUBDIVISION WITH SOLUTION; IF ABS(Y2)<EPS SOLUTION IS ASSUMED
C
      IF (ABS(Y2) .GT. EPS) Y2 = SSH_FUNCH5A (PSI6LO)
      GOTO 50
C
C *** PERFORM BISECTION ***********************************************
C
20    DO 30 I=1,MAXIT
         X3 = 0.5*(X1+X2)
         Y3 = SSH_FUNCH5A (X3)
         IF (SIGN(1.d0,Y1)*SIGN(1.d0,Y3) .LE. ZERO) THEN  ! (Y1*Y3 .LE. ZERO)
            Y2    = Y3
            X2    = X3
         ELSE
            Y1    = Y3
            X1    = X3
         ENDIF
         IF (ABS(X2-X1) .LE. EPS*X1) GOTO 40
30    CONTINUE
      CALL SSH_PUSHERR (0002, 'CALCH5')    ! WARNING ERROR: NO CONVERGENCE
C
C *** CONVERGED ; RETURN **********************************************
C
40    X3 = 0.5*(X1+X2)
      Y3 = SSH_FUNCH5A (X3)
C 
C *** CALCULATE HSO4 SPECIATION AND RETURN *******************************
C
50    CONTINUE
      IF (MOLAL(1).GT.TINY .AND. MOLAL(5).GT.TINY) THEN
         CALL SSH_CALCHS4 (MOLAL(1), MOLAL(5), ZERO, DELTA)
c         MOLAL(1) = MOLAL(1) - DELTA                     ! H+   EFECT
         MOLAL(1) = MAX(MOLAL(1) - DELTA,ZERO) ! YK 
c      write(*,*) "isofwd.f: L3808 MOLAL(1): ",molal(1)
c         MOLAL(5) = MOLAL(5) - DELTA                     ! SO4  EFFECT
         MOLAL(5) = MAX(MOLAL(5) - DELTA,ZERO) ! YK
         MOLAL(6) = DELTA                                ! HSO4 EFFECT
      ENDIF
C
      RETURN
C
C *** END OF SUBROUTINE SSH_CALCH5 ******************************************
C
      END




C=======================================================================
C
C *** ISORROPIA CODE
C *** SUBROUTINE SSH_FUNCH5A
C *** CASE H5
C
C     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
C     1. SULFATE POOR (SULRAT > 2.0) ; SODIUM RICH (SODRAT >= 2.0)
C     2. THERE IS BOTH A LIQUID & SOLID PHASE
C     3. SOLIDS POSSIBLE : NONE
C
C *** COPYRIGHT 1996-2006, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY,
C *** GEORGIA INSTITUTE OF TECHNOLOGY
C *** WRITTEN BY ATHANASIOS NENES
C *** UPDATED BY CHRISTOS FOUNTOUKIS
C
C=======================================================================
C
      DOUBLE PRECISION FUNCTION SSH_FUNCH5A (X)
      INCLUDE 'isrpia.inc'
C
      COMMON /SSHSOLUT/ CHI1, CHI2, CHI3, CHI4, CHI5, CHI6, CHI7, CHI8,
     &               PSI1, PSI2, PSI3, PSI4, PSI5, PSI6, PSI7, PSI8,
     &               A1,   A2,   A3,   A4,   A5,   A6,   A7,   A8


C
C *** SETUP PARAMETERS ************************************************
C
      PSI6   = X
      PSI1   = CHI1
      PSI2   = ZERO
      PSI3   = ZERO
      PSI7   = CHI7
      PSI8   = CHI8 
      FRST   = .TRUE.
      CALAIN = .TRUE. 
C
C *** SOLVE EQUATIONS ; WITH ITERATIONS FOR ACTIVITY COEF. ************
C
      DO 10 I=1,NSWEEP
C
      A1  = XK5 *(WATER/GAMA(2))**3.0
      A4  = (XK2/XKW)*R*TEMP*(GAMA(10)/GAMA(5))**2.0
      A5  = XK4 *R*TEMP*(WATER/GAMA(10))**2.0
      A6  = XK3 *R*TEMP*(WATER/GAMA(11))**2.0
      A7  = XK8 *(WATER/GAMA(1))**2.0
      A8  = XK9 *(WATER/GAMA(3))**2.0
      A9  = XK1*WATER/GAMA(7)*(GAMA(8)/GAMA(7))**2.
C
C  CALCULATE DISSOCIATION QUANTITIES
C
      PSI5 = CHI5*(PSI6+PSI7) - A6/A5*PSI8*(CHI6-PSI6-PSI3)
      PSI5 = PSI5/(A6/A5*(CHI6-PSI6-PSI3) + PSI6 + PSI7)
      PSI5 = MAX(PSI5, TINY)
C
      IF (W(3).GT.TINY .AND. WATER.GT.TINY) THEN  ! First try 3rd order soln
         BB   =-(CHI4 + PSI6 + PSI5 + 1.d0/A4)
         CC   = CHI4*(PSI5+PSI6)
         DD   = BB*BB-4.d0*CC
         PSI4 =0.5d0*(-BB - SQRT(DD))
         PSI4 = MIN(PSI4,CHI4)
      ELSE
         PSI4 = TINY
      ENDIF
C
      IF (CHI1.GT.TINY .AND. WATER.GT.TINY) THEN     ! NA2SO4 DISSOLUTION
         AA = PSI7+PSI8
         BB = AA*AA
         CC =-A1/4.D0
         CALL SSH_POLY3 (AA, BB, CC, PSI1, ISLV)
         IF (ISLV.EQ.0) THEN
             PSI1 = MIN (PSI1, CHI1)
         ELSE
             PSI1 = ZERO
         ENDIF
      ENDIF
C
C *** CALCULATE SPECIATION ********************************************
C
      MOLAL (2) = PSI8 + PSI7 + 2.D0*PSI1                ! NAI
      MOLAL (3) = PSI4                                   ! NH4I
      MOLAL (4) = PSI6 + PSI7                            ! CLI
      MOLAL (5) = PSI2 + PSI1                            ! SO4I
      MOLAL (6) = ZERO
      MOLAL (7) = PSI5 + PSI8                            ! NO3I
C
C *** FOR INTERACTION WITH HYDROPHILIC ORGANIC COMPOUNDS
      SMIN      = 2.d0*MOLAL(5)+MOLAL(7)+MOLAL(4)-MOLAL(2)-MOLAL(3)
     &     + ORGANION
      CALL SSH_CALCPH (SMIN, HI, OHI)
      MOLAL (1) = HI
C 
      GNH3      = MAX(CHI4 - PSI4, TINY)
      GHNO3     = MAX(CHI5 - PSI5, TINY)
      GHCL      = MAX(CHI6 - PSI6, TINY)
C
      CNH42S4   = ZERO
      CNH4NO3   = ZERO
      CNACL     = MAX(CHI7 - PSI7, ZERO)
      CNANO3    = MAX(CHI8 - PSI8, ZERO)
      CNA2SO4   = MAX(CHI1 - PSI1, ZERO) 
C
      CALL SSH_CALCMR                               ! Water content
C
C *** CALCULATE ACTIVITIES OR TERMINATE INTERNAL LOOP *****************
C
      IF (FRST.AND.CALAOU .OR. .NOT.FRST.AND.CALAIN) THEN
         CALL SSH_CALCACT     
      ELSE
         GOTO 20
      ENDIF
10    CONTINUE
C
C *** CALCULATE FUNCTION VALUE FOR OUTER LOOP ***************************
C
20    SSH_FUNCH5A = MOLAL(3)*MOLAL(4)/GHCL/GNH3/A6/A4 - ONE
C
      RETURN
C
C *** END OF FUNCTION SSH_FUNCH5A *******************************************
C
      END

C=======================================================================
C
C *** ISORROPIA CODE
C *** SUBROUTINE SSH_CALCH4
C *** CASE H4
C
C     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
C     1. SULFATE POOR (SULRAT > 2.0) ; SODIUM RICH (SODRAT >= 2.0)
C     2. THERE IS BOTH A LIQUID & SOLID PHASE
C     3. SOLIDS POSSIBLE : NA2SO4
C
C *** COPYRIGHT 1996-2006, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY,
C *** GEORGIA INSTITUTE OF TECHNOLOGY
C *** WRITTEN BY ATHANASIOS NENES
C *** UPDATED BY CHRISTOS FOUNTOUKIS
C
C=======================================================================
C
      SUBROUTINE SSH_CALCH4
      INCLUDE 'isrpia.inc'
C
      COMMON /SSHSOLUT/ CHI1, CHI2, CHI3, CHI4, CHI5, CHI6, CHI7, CHI8,
     &               PSI1, PSI2, PSI3, PSI4, PSI5, PSI6, PSI7, PSI8,
     &               A1,   A2,   A3,   A4,   A5,   A6,   A7,   A8


C
C *** REGIME DEPENDS ON THE EXISTANCE OF NITRATES ***********************
C
      IF (W(4).LE.TINY .AND. W(5).LE.TINY) THEN  
         SCASE = 'H4'  
         CALL SSH_CALCH1A
         SCASE = 'H4'  
         RETURN
      ENDIF
C
C *** SETUP PARAMETERS ************************************************
C
      CALAOU = .TRUE.   
      CHI1   = W(2)                                ! CNA2SO4
      CHI2   = ZERO                                ! CNH42S4
      CHI3   = ZERO                                ! CNH4CL
      FRNA   = MAX (W(1)-2.D0*CHI1, ZERO)       
      CHI8   = MIN (FRNA, W(4))                    ! CNANO3
      CHI4   = W(3)                                ! NH3(g)
      CHI5   = MAX (W(4)-CHI8, ZERO)               ! HNO3(g)
      CHI7   = MIN (MAX(FRNA-CHI8, ZERO), W(5))    ! CNACL
      CHI6   = MAX (W(5)-CHI7, ZERO)               ! HCL(g)
C
      PSI6LO = TINY                  
      PSI6HI = CHI6-TINY    ! MIN(CHI6-TINY, CHI4)
C
C *** INITIAL VALUES FOR BISECTION ************************************
C
      X1 = PSI6LO
      Y1 = SSH_FUNCH4A (X1)
      IF (ABS(Y1).LE.EPS .OR. CHI6.LE.TINY) GOTO 50  
C
C *** ROOT TRACKING ; FOR THE RANGE OF HI AND LO **********************
C
      DX = (PSI6HI-PSI6LO)/FLOAT(NDIV)
      DO 10 I=1,NDIV
         X2 = X1+DX 
         Y2 = SSH_FUNCH4A (X2)
         IF (SIGN(1.d0,Y1)*SIGN(1.d0,Y2).LT.ZERO) GOTO 20  ! (Y1*Y2.LT.ZERO)
         X1 = X2
         Y1 = Y2
10    CONTINUE
C
C *** NO SUBDIVISION WITH SOLUTION; IF ABS(Y2)<EPS SOLUTION IS ASSUMED
C
      IF (ABS(Y2) .GT. EPS) Y2 = SSH_FUNCH4A (PSI6LO)
      GOTO 50
C
C *** PERFORM BISECTION ***********************************************
C
20    DO 30 I=1,MAXIT
         X3 = 0.5*(X1+X2)
         Y3 = SSH_FUNCH4A (X3)
         IF (SIGN(1.d0,Y1)*SIGN(1.d0,Y3) .LE. ZERO) THEN  ! (Y1*Y3 .LE. ZERO)
            Y2    = Y3
            X2    = X3
         ELSE
            Y1    = Y3
            X1    = X3
         ENDIF
         IF (ABS(X2-X1) .LE. EPS*X1) GOTO 40
30    CONTINUE
      CALL SSH_PUSHERR (0002, 'CALCH4')    ! WARNING ERROR: NO CONVERGENCE
C
C *** CONVERGED ; RETURN **********************************************
C
40    X3 = 0.5*(X1+X2)
      Y3 = SSH_FUNCH4A (X3)
C 
C *** CALCULATE HSO4 SPECIATION AND RETURN *******************************
C
50    CONTINUE
      IF (MOLAL(1).GT.TINY .AND. MOLAL(5).GT.TINY) THEN
         CALL SSH_CALCHS4 (MOLAL(1), MOLAL(5), ZERO, DELTA)
c         MOLAL(1) = MOLAL(1) - DELTA                      ! H+   EFFECT
         MOLAL(1) = MAX(MOLAL(1) - DELTA,ZERO) ! YK  
c      write(*,*) "isofwd.f: L4043 MOLAL(1): ",molal(1)
c         MOLAL(5) = MOLAL(5) - DELTA                      ! SO4  EFFECT
         MOLAL(5) = MAX(MOLAL(5) - DELTA,ZERO) ! YK
         MOLAL(6) = DELTA                                 ! HSO4 EFFECT
      ENDIF
C
      RETURN
C
C *** END OF SUBROUTINE SSH_CALCH4 ******************************************
C
      END




C=======================================================================
C
C *** ISORROPIA CODE
C *** SUBROUTINE SSH_FUNCH4A
C *** CASE H4
C
C     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
C     1. SULFATE POOR (SULRAT > 2.0) ; SODIUM RICH (SODRAT >= 2.0)
C     2. THERE IS BOTH A LIQUID & SOLID PHASE
C     3. SOLIDS POSSIBLE : (NH4)2SO4, NH4CL, NA2SO4
C
C *** COPYRIGHT 1996-2006, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY,
C *** GEORGIA INSTITUTE OF TECHNOLOGY
C *** WRITTEN BY ATHANASIOS NENES
C *** UPDATED BY CHRISTOS FOUNTOUKIS
C
C=======================================================================
C
      DOUBLE PRECISION FUNCTION SSH_FUNCH4A (X)
      INCLUDE 'isrpia.inc'
C
      COMMON /SSHSOLUT/ CHI1, CHI2, CHI3, CHI4, CHI5, CHI6, CHI7, CHI8,
     &               PSI1, PSI2, PSI3, PSI4, PSI5, PSI6, PSI7, PSI8,
     &               A1,   A2,   A3,   A4,   A5,   A6,   A7,   A8


C
C *** SETUP PARAMETERS ************************************************
C
      PSI6   = X
      PSI1   = CHI1
      PSI2   = ZERO
      PSI3   = ZERO
      PSI7   = CHI7
      PSI8   = CHI8 
      FRST   = .TRUE.
      CALAIN = .TRUE. 
C
C *** SOLVE EQUATIONS ; WITH ITERATIONS FOR ACTIVITY COEF. ************
C
      DO 10 I=1,NSWEEP
C
      A1  = XK5 *(WATER/GAMA(2))**3.0
      A4  = (XK2/XKW)*R*TEMP*(GAMA(10)/GAMA(5))**2.0
      A5  = XK4 *R*TEMP*(WATER/GAMA(10))**2.0
      A6  = XK3 *R*TEMP*(WATER/GAMA(11))**2.0
      A7  = XK8 *(WATER/GAMA(1))**2.0
      A8  = XK9 *(WATER/GAMA(3))**2.0
      A9  = XK1*WATER/GAMA(7)*(GAMA(8)/GAMA(7))**2.
C
C  CALCULATE DISSOCIATION QUANTITIES
C
      PSI5 = CHI5*(PSI6+PSI7) - A6/A5*PSI8*(CHI6-PSI6-PSI3)
      PSI5 = PSI5/(A6/A5*(CHI6-PSI6-PSI3) + PSI6 + PSI7)
      PSI5 = MAX(PSI5, TINY)
C
      IF (W(3).GT.TINY .AND. WATER.GT.TINY) THEN  ! First try 3rd order soln
         BB   =-(CHI4 + PSI6 + PSI5 + 1.d0/A4)
         CC   = CHI4*(PSI5+PSI6)
         DD   = BB*BB-4.d0*CC
         PSI4 =0.5d0*(-BB - SQRT(DD))
         PSI4 = MIN(PSI4,CHI4)
      ELSE
         PSI4 = TINY
      ENDIF
C
      IF (CHI1.GT.TINY .AND. WATER.GT.TINY) THEN     ! NA2SO4 DISSOLUTION
         AA = PSI7+PSI8
         BB = AA*AA
         CC =-A1/4.D0
         CALL SSH_POLY3 (AA, BB, CC, PSI1, ISLV)
         IF (ISLV.EQ.0) THEN
             PSI1 = MIN (PSI1, CHI1)
         ELSE
             PSI1 = ZERO
         ENDIF
      ENDIF
C
C *** CALCULATE SPECIATION ********************************************
C
      MOLAL (2) = PSI8 + PSI7 + 2.D0*PSI1                ! NAI
      MOLAL (3) = PSI4                                   ! NH4I
      MOLAL (4) = PSI6 + PSI7                            ! CLI
      MOLAL (5) = PSI2 + PSI1                            ! SO4I
      MOLAL (6) = ZERO
      MOLAL (7) = PSI5 + PSI8                            ! NO3I
C
C *** FOR INTERACTION WITH HYDROPHILIC ORGANIC COMPOUNDS
      SMIN      = 2.d0*MOLAL(5)+MOLAL(7)+MOLAL(4)-MOLAL(2)-MOLAL(3)
     &     + ORGANION
      CALL SSH_CALCPH (SMIN, HI, OHI)
      MOLAL (1) = HI
C 
C *** CALCULATE GAS / SOLID SPECIES (LIQUID IN MOLAL ALREADY) *********
C
      GNH3      = MAX(CHI4 - PSI4, TINY)
      GHNO3     = MAX(CHI5 - PSI5, TINY)
      GHCL      = MAX(CHI6 - PSI6, TINY)
C
      CNH42S4   = ZERO
      CNH4NO3   = ZERO
      CNACL     = MAX(CHI7 - PSI7, ZERO)
      CNANO3    = MAX(CHI8 - PSI8, ZERO)
      CNA2SO4   = MAX(CHI1 - PSI1, ZERO) 
C      
C *** NH4Cl(s) calculations
C
      A3   = XK6 /(R*TEMP*R*TEMP)
      DELT = MIN(GNH3, GHCL)
      BB = -(GNH3+GHCL)
      CC = GNH3*GHCL-A3
      DD = BB*BB - 4.D0*CC
      PSI31 = 0.5D0*(-BB + SQRT(DD))
      PSI32 = 0.5D0*(-BB - SQRT(DD))
      IF (DELT-PSI31.GT.ZERO .AND. PSI31.GT.ZERO) THEN
         PSI3 = PSI31
      ELSEIF (DELT-PSI32.GT.ZERO .AND. PSI32.GT.ZERO) THEN
         PSI3 = PSI32
      ELSE
         PSI3 = ZERO
      ENDIF
C 
C *** CALCULATE GAS / SOLID SPECIES (LIQUID IN MOLAL ALREADY) *********
C
      GNH3    = MAX(GNH3 - PSI3, TINY)
      GHCL    = MAX(GHCL - PSI3, TINY)
      CNH4CL  = PSI3
C 
      CALL SSH_CALCMR                           ! Water content
C
C *** CALCULATE ACTIVITIES OR TERMINATE INTERNAL LOOP *****************
C
      IF (FRST.AND.CALAOU .OR. .NOT.FRST.AND.CALAIN) THEN
         CALL SSH_CALCACT     
      ELSE
         GOTO 20
      ENDIF
10    CONTINUE
C
C *** CALCULATE FUNCTION VALUE FOR OUTER LOOP ***************************
C
20    SSH_FUNCH4A = MOLAL(3)*MOLAL(4)/GHCL/GNH3/A6/A4 - ONE
C
      RETURN
C
C *** END OF FUNCTION SSH_FUNCH4A *******************************************
C
      END

C=======================================================================
C
C *** ISORROPIA CODE
C *** SUBROUTINE SSH_CALCH3
C *** CASE H3
C
C     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
C     1. SULFATE POOR (SULRAT > 2.0) ; SODIUM RICH (SODRAT >= 2.0)
C     2. THERE IS BOTH A LIQUID & SOLID PHASE
C     3. SOLIDS POSSIBLE : NH4CL, NA2SO4
C
C *** COPYRIGHT 1996-2006, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY,
C *** GEORGIA INSTITUTE OF TECHNOLOGY
C *** WRITTEN BY ATHANASIOS NENES
C *** UPDATED BY CHRISTOS FOUNTOUKIS
C
C=======================================================================
C
      SUBROUTINE SSH_CALCH3
      INCLUDE 'isrpia.inc'
C
      COMMON /SSHSOLUT/ CHI1, CHI2, CHI3, CHI4, CHI5, CHI6, CHI7, CHI8,
     &               PSI1, PSI2, PSI3, PSI4, PSI5, PSI6, PSI7, PSI8,
     &               A1,   A2,   A3,   A4,   A5,   A6,   A7,   A8


C
C *** REGIME DEPENDS ON THE EXISTANCE OF NITRATES ***********************
C
      IF (W(4).LE.TINY) THEN        ! NO3 NOT EXIST, WATER NOT POSSIBLE
         SCASE = 'H3'  
         CALL SSH_CALCH1A
         SCASE = 'H3'  
         RETURN
      ENDIF
C
C *** SETUP PARAMETERS ************************************************
C
      CALAOU = .TRUE.   
      CHI1   = W(2)                                ! CNA2SO4
      CHI2   = ZERO                                ! CNH42S4
      CHI3   = ZERO                                ! CNH4CL
      FRNA   = MAX (W(1)-2.D0*CHI1, ZERO)       
      CHI8   = MIN (FRNA, W(4))                    ! CNANO3
      CHI4   = W(3)                                ! NH3(g)
      CHI5   = MAX (W(4)-CHI8, ZERO)               ! HNO3(g)
      CHI7   = MIN (MAX(FRNA-CHI8, ZERO), W(5))    ! CNACL
      CHI6   = MAX (W(5)-CHI7, ZERO)               ! HCL(g)
C
      PSI6LO = TINY                  
      PSI6HI = CHI6-TINY    ! MIN(CHI6-TINY, CHI4)
C
C *** INITIAL VALUES FOR BISECTION ************************************
C
      X1 = PSI6LO
      Y1 = SSH_FUNCH3A (X1)
      IF (ABS(Y1).LE.EPS .OR. CHI6.LE.TINY) GOTO 50  
C
C *** ROOT TRACKING ; FOR THE RANGE OF HI AND LO **********************
C
      DX = (PSI6HI-PSI6LO)/FLOAT(NDIV)
      DO 10 I=1,NDIV
         X2 = X1+DX 
         Y2 = SSH_FUNCH3A (X2)
         IF (SIGN(1.d0,Y1)*SIGN(1.d0,Y2).LT.ZERO) GOTO 20  ! (Y1*Y2.LT.ZERO)
         X1 = X2
         Y1 = Y2
10    CONTINUE
C
C *** NO SUBDIVISION WITH SOLUTION; IF ABS(Y2)<EPS SOLUTION IS ASSUMED
C
      IF (ABS(Y2) .GT. EPS) Y2 = SSH_FUNCH3A (PSI6LO)
      GOTO 50
C
C *** PERFORM BISECTION ***********************************************
C
20    DO 30 I=1,MAXIT
         X3 = 0.5*(X1+X2)
         Y3 = SSH_FUNCH3A (X3)
         IF (SIGN(1.d0,Y1)*SIGN(1.d0,Y3) .LE. ZERO) THEN  ! (Y1*Y3 .LE. ZERO)
            Y2    = Y3
            X2    = X3
         ELSE
            Y1    = Y3
            X1    = X3
         ENDIF
         IF (ABS(X2-X1) .LE. EPS*X1) GOTO 40
30    CONTINUE
      CALL SSH_PUSHERR (0002, 'CALCH3')    ! WARNING ERROR: NO CONVERGENCE
C
C *** CONVERGED ; RETURN **********************************************
C
40    X3 = 0.5*(X1+X2)
      Y3 = SSH_FUNCH3A (X3)
C 
C *** CALCULATE HSO4 SPECIATION AND RETURN *******************************
C
50    CONTINUE
      IF (MOLAL(1).GT.TINY .AND. MOLAL(5).GT.TINY) THEN
         CALL SSH_CALCHS4 (MOLAL(1), MOLAL(5), ZERO, DELTA)
c         MOLAL(1) = MOLAL(1) - DELTA                     ! H+   EFFECT
         MOLAL(1) = MAX(MOLAL(1) - DELTA,ZERO) !YK 
c      write(*,*) "isofwd.f: L4303 MOLAL(1): ",molal(1)
c         MOLAL(5) = MOLAL(5) - DELTA                     ! SO4  EFFECT
         MOLAL(5) = MAX(MOLAL(5) - DELTA,ZERO) !YK   
         MOLAL(6) = DELTA                                ! HSO4 EFFECT
      ENDIF
C
      RETURN
C
C *** END OF SUBROUTINE SSH_CALCH3 ******************************************
C
      END




C=======================================================================
C
C *** ISORROPIA CODE
C *** SUBROUTINE SSH_FUNCH3A
C *** CASE H3
C
C     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
C     1. SULFATE POOR (SULRAT > 2.0) ; SODIUM RICH (SODRAT >= 2.0)
C     2. THERE IS BOTH A LIQUID & SOLID PHASE
C     3. SOLIDS POSSIBLE : (NH4)2SO4, NH4CL, NA2SO4
C
C *** COPYRIGHT 1996-2006, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY,
C *** GEORGIA INSTITUTE OF TECHNOLOGY
C *** WRITTEN BY ATHANASIOS NENES
C *** UPDATED BY CHRISTOS FOUNTOUKIS
C
C=======================================================================
C
      DOUBLE PRECISION FUNCTION SSH_FUNCH3A (X)
      INCLUDE 'isrpia.inc'
C
      COMMON /SSHSOLUT/ CHI1, CHI2, CHI3, CHI4, CHI5, CHI6, CHI7, CHI8,
     &               PSI1, PSI2, PSI3, PSI4, PSI5, PSI6, PSI7, PSI8,
     &               A1,   A2,   A3,   A4,   A5,   A6,   A7,   A8


C
C *** SETUP PARAMETERS ************************************************
C
      PSI6   = X
      PSI1   = CHI1
      PSI2   = ZERO
      PSI3   = ZERO
      PSI7   = CHI7
      PSI8   = CHI8 
      FRST   = .TRUE.
      CALAIN = .TRUE. 
C
C *** SOLVE EQUATIONS ; WITH ITERATIONS FOR ACTIVITY COEF. ************
C
      DO 10 I=1,NSWEEP
C
      A1  = XK5 *(WATER/GAMA(2))**3.0
      A4  = (XK2/XKW)*R*TEMP*(GAMA(10)/GAMA(5))**2.0
      A5  = XK4 *R*TEMP*(WATER/GAMA(10))**2.0
      A6  = XK3 *R*TEMP*(WATER/GAMA(11))**2.0
      A7  = XK8 *(WATER/GAMA(1))**2.0
      A8  = XK9 *(WATER/GAMA(3))**2.0
      A9  = XK1*WATER/GAMA(7)*(GAMA(8)/GAMA(7))**2.
C
C  CALCULATE DISSOCIATION QUANTITIES
C
      PSI5 = CHI5*(PSI6+PSI7) - A6/A5*PSI8*(CHI6-PSI6-PSI3)
      PSI5 = PSI5/(A6/A5*(CHI6-PSI6-PSI3) + PSI6 + PSI7)
      PSI5 = MAX(PSI5, TINY)
C
      IF (W(3).GT.TINY .AND. WATER.GT.TINY) THEN  ! First try 3rd order soln
         BB   =-(CHI4 + PSI6 + PSI5 + 1.d0/A4)
         CC   = CHI4*(PSI5+PSI6)
         DD   = BB*BB-4.d0*CC
         PSI4 =0.5d0*(-BB - SQRT(DD))
         PSI4 = MIN(PSI4,CHI4)
      ELSE
         PSI4 = TINY
      ENDIF
C
      IF (CHI7.GT.TINY .AND. WATER.GT.TINY) THEN     ! NACL DISSOLUTION
         DIAK = (PSI8-PSI6)**2.D0 + 4.D0*A7
         PSI7 = 0.5D0*( -(PSI8+PSI6) + SQRT(DIAK) )
         PSI7 = MAX(MIN(PSI7, CHI7), ZERO)
      ENDIF
C
      IF (CHI1.GT.TINY .AND. WATER.GT.TINY) THEN     ! NA2SO4 DISSOLUTION
         AA = PSI7+PSI8
         BB = AA*AA
         CC =-A1/4.D0
         CALL SSH_POLY3 (AA, BB, CC, PSI1, ISLV)
         IF (ISLV.EQ.0) THEN
             PSI1 = MIN (PSI1, CHI1)
         ELSE
             PSI1 = ZERO
         ENDIF
      ENDIF
C
C *** CALCULATE SPECIATION ********************************************
C
      MOLAL (2) = PSI8 + PSI7 + 2.D0*PSI1             ! NAI
      MOLAL (3) = PSI4                                ! NH4I
      MOLAL (4) = PSI6 + PSI7                         ! CLI
      MOLAL (5) = PSI2 + PSI1                         ! SO4I
      MOLAL (6) = ZERO
      MOLAL (7) = PSI5 + PSI8                         ! NO3I
C
C *** FOR INTERACTION WITH HYDROPHILIC ORGANIC COMPOUNDS
      SMIN      = 2.d0*MOLAL(5)+MOLAL(7)+MOLAL(4)-MOLAL(2)-MOLAL(3)
     &     + ORGANION
      CALL SSH_CALCPH (SMIN, HI, OHI)
      MOLAL (1) = HI
C 
C *** CALCULATE GAS / SOLID SPECIES (LIQUID IN MOLAL ALREADY) *********
C
      GNH3      = MAX(CHI4 - PSI4, TINY)
      GHNO3     = MAX(CHI5 - PSI5, TINY)
      GHCL      = MAX(CHI6 - PSI6, TINY)
C
      CNH42S4   = ZERO
      CNH4NO3   = ZERO
      CNACL     = MAX(CHI7 - PSI7, ZERO)
      CNANO3    = MAX(CHI8 - PSI8, ZERO)
      CNA2SO4   = MAX(CHI1 - PSI1, ZERO) 
C      
C *** NH4Cl(s) calculations
C
      A3   = XK6 /(R*TEMP*R*TEMP)
      DELT = MIN(GNH3, GHCL)
      BB = -(GNH3+GHCL)
      CC = GNH3*GHCL-A3
      DD = BB*BB - 4.D0*CC
      PSI31 = 0.5D0*(-BB + SQRT(DD))
      PSI32 = 0.5D0*(-BB - SQRT(DD))
      IF (DELT-PSI31.GT.ZERO .AND. PSI31.GT.ZERO) THEN
         PSI3 = PSI31
      ELSEIF (DELT-PSI32.GT.ZERO .AND. PSI32.GT.ZERO) THEN
         PSI3 = PSI32
      ELSE
         PSI3 = ZERO
      ENDIF
C 
C *** CALCULATE GAS / SOLID SPECIES (LIQUID IN MOLAL ALREADY) *********
C
      GNH3    = MAX(GNH3 - PSI3, TINY)
      GHCL    = MAX(GHCL - PSI3, TINY)
      CNH4CL  = PSI3
C
      CALL SSH_CALCMR                                 ! Water content
C
C *** CALCULATE ACTIVITIES OR TERMINATE INTERNAL LOOP *****************
C
      IF (FRST.AND.CALAOU .OR. .NOT.FRST.AND.CALAIN) THEN
         CALL SSH_CALCACT     
      ELSE
         GOTO 20
      ENDIF
10    CONTINUE
C
C *** CALCULATE FUNCTION VALUE FOR OUTER LOOP ***************************
C
20    SSH_FUNCH3A = MOLAL(3)*MOLAL(4)/GHCL/GNH3/A6/A4 - ONE
C
      RETURN
C
C *** END OF FUNCTION SSH_FUNCH3A *******************************************
C
      END

C=======================================================================
C
C *** ISORROPIA CODE
C *** SUBROUTINE SSH_CALCH2
C *** CASE H2
C
C     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
C     1. SULFATE POOR (SULRAT > 2.0) ; SODIUM RICH (SODRAT >= 2.0)
C     2. SOLID & LIQUID AEROSOL POSSIBLE
C     3. SOLIDS POSSIBLE : NH4Cl, NA2SO4, NANO3, NACL
C
C     THERE ARE THREE REGIMES IN THIS CASE:
C     1. NH4NO3(s) POSSIBLE. LIQUID & SOLID AEROSOL (SUBROUTINE SSH_CALCH2A)
C     2. NH4NO3(s) NOT POSSIBLE, AND RH < MDRH. SOLID AEROSOL ONLY 
C     3. NH4NO3(s) NOT POSSIBLE, AND RH >= MDRH. (MDRH REGION)
C
C     REGIMES 2. AND 3. ARE CONSIDERED TO BE THE SAME AS CASES H1A, H2B
C     RESPECTIVELY (BECAUSE MDRH POINTS COINCIDE).
C
C *** COPYRIGHT 1996-2006, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY,
C *** GEORGIA INSTITUTE OF TECHNOLOGY
C *** WRITTEN BY ATHANASIOS NENES
C *** UPDATED BY CHRISTOS FOUNTOUKIS
C
C=======================================================================
C
      SUBROUTINE SSH_CALCH2
      INCLUDE 'isrpia.inc'
      EXTERNAL SSH_CALCH1A, SSH_CALCH3
C
C *** REGIME DEPENDS ON THE EXISTANCE OF NITRATES ***********************
C
      IF (W(4).GT.TINY) THEN        ! NO3 EXISTS, WATER POSSIBLE
         SCASE = 'H2 ; SUBCASE 1'  
         CALL SSH_CALCH2A                                   
         SCASE = 'H2 ; SUBCASE 1'  
      ELSE                          ! NO3 NON EXISTANT, WATER NOT POSSIBLE
         SCASE = 'H2 ; SUBCASE 1'  
         CALL SSH_CALCH1A
         SCASE = 'H2 ; SUBCASE 1'  
      ENDIF
C
      IF (WATER.LE.TINY .AND. RH.LT.DRMH2) THEN      ! DRY AEROSOL
         SCASE = 'H2 ; SUBCASE 2'  
C
      ELSEIF (WATER.LE.TINY .AND. RH.GE.DRMH2) THEN  ! MDRH OF H2
         SCASE = 'H2 ; SUBCASE 3'
         CALL SSH_CALCMDRH (RH, DRMH2, DRNANO3, SSH_CALCH1A, SSH_CALCH3)
         SCASE = 'H2 ; SUBCASE 3'
      ENDIF
C 
      RETURN
C
C *** END OF SUBROUTINE SSH_CALCH2 ******************************************
C
      END




C=======================================================================
C
C *** ISORROPIA CODE
C *** SUBROUTINE SSH_CALCH2A
C *** CASE H2 ; SUBCASE 1
C
C     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
C     1. SULFATE POOR (SULRAT > 2.0) ; SODIUM RICH (SODRAT >= 2.0)
C     2. THERE IS BOTH A LIQUID & SOLID PHASE
C     3. SOLIDS POSSIBLE : NH4CL, NA2SO4
C
C *** COPYRIGHT 1996-2006, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY,
C *** GEORGIA INSTITUTE OF TECHNOLOGY
C *** WRITTEN BY ATHANASIOS NENES
C *** UPDATED BY CHRISTOS FOUNTOUKIS
C
C=======================================================================
C
      SUBROUTINE SSH_CALCH2A
      INCLUDE 'isrpia.inc'
C
      COMMON /SSHSOLUT/ CHI1, CHI2, CHI3, CHI4, CHI5, CHI6, CHI7, CHI8,
     &               PSI1, PSI2, PSI3, PSI4, PSI5, PSI6, PSI7, PSI8,
     &               A1,   A2,   A3,   A4,   A5,   A6,   A7,   A8


C
C *** SETUP PARAMETERS ************************************************
C
      CALAOU = .TRUE.   
      CHI1   = W(2)                                ! CNA2SO4
      CHI2   = ZERO                                ! CNH42S4
      CHI3   = ZERO                                ! CNH4CL
      FRNA   = MAX (W(1)-2.D0*CHI1, ZERO)       
      CHI8   = MIN (FRNA, W(4))                    ! CNANO3
      CHI4   = W(3)                                ! NH3(g)
      CHI5   = MAX (W(4)-CHI8, ZERO)               ! HNO3(g)
      CHI7   = MIN (MAX(FRNA-CHI8, ZERO), W(5))    ! CNACL
      CHI6   = MAX (W(5)-CHI7, ZERO)               ! HCL(g)
C
      PSI6LO = TINY                  
      PSI6HI = CHI6-TINY    ! MIN(CHI6-TINY, CHI4)
C
C *** INITIAL VALUES FOR BISECTION ************************************
C
      X1 = PSI6LO
      Y1 = SSH_FUNCH2A (X1)
      IF (ABS(Y1).LE.EPS .OR. CHI6.LE.TINY) GOTO 50  
C
C *** ROOT TRACKING ; FOR THE RANGE OF HI AND LO **********************
C
      DX = (PSI6HI-PSI6LO)/FLOAT(NDIV)
      DO 10 I=1,NDIV
         X2 = X1+DX 
         Y2 = SSH_FUNCH2A (X2)
         IF (SIGN(1.d0,Y1)*SIGN(1.d0,Y2).LT.ZERO) GOTO 20  ! (Y1*Y2.LT.ZERO)
         X1 = X2
         Y1 = Y2
10    CONTINUE
C
C *** NO SUBDIVISION WITH SOLUTION; IF ABS(Y2)<EPS SOLUTION IS ASSUMED
C
      IF (Y2 .GT. EPS) Y2 = SSH_FUNCH2A (PSI6LO)
      GOTO 50
C
C *** PERFORM BISECTION ***********************************************
C
20    DO 30 I=1,MAXIT
         X3 = 0.5*(X1+X2)
         Y3 = SSH_FUNCH2A (X3)
         IF (SIGN(1.d0,Y1)*SIGN(1.d0,Y3) .LE. ZERO) THEN  ! (Y1*Y3 .LE. ZERO)
            Y2    = Y3
            X2    = X3
         ELSE
            Y1    = Y3
            X1    = X3
         ENDIF
         IF (ABS(X2-X1) .LE. EPS*X1) GOTO 40
30    CONTINUE
      CALL SSH_PUSHERR (0002, 'CALCH2A')    ! WARNING ERROR: NO CONVERGENCE
C
C *** CONVERGED ; RETURN **********************************************
C
40    X3 = 0.5*(X1+X2)
      Y3 = SSH_FUNCH2A (X3)
C 
C *** CALCULATE HSO4 SPECIATION AND RETURN *******************************
C
50    CONTINUE
      IF (MOLAL(1).GT.TINY .AND. MOLAL(5).GT.TINY) THEN
         CALL SSH_CALCHS4 (MOLAL(1), MOLAL(5), ZERO, DELTA)
c         MOLAL(1) = MOLAL(1) - DELTA                    ! H+   EFFECT
         MOLAL(1) = MAX(MOLAL(1) - DELTA,ZERO) ! YK
c      write(*,*) "isofwd.f: L4625 MOLAL(1): ",molal(1)
c         MOLAL(5) = MOLAL(5) - DELTA                    ! SO4  EFFECT
         MOLAL(5) = MAX(MOLAL(5) - DELTA,ZERO) ! YK    
         MOLAL(6) = DELTA                               ! HSO4 EFFECT
      ENDIF
C
      RETURN
C
C *** END OF SUBROUTINE SSH_CALCH2A ******************************************
C
      END




C=======================================================================
C
C *** ISORROPIA CODE
C *** SUBROUTINE SSH_FUNCH2A
C *** CASE H2 ; SUBCASE 1
C
C     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
C     1. SULFATE POOR (SULRAT > 2.0) ; SODIUM RICH (SODRAT >= 2.0)
C     2. THERE IS BOTH A LIQUID & SOLID PHASE
C     3. SOLIDS POSSIBLE : (NH4)2SO4, NH4CL, NA2SO4
C
C *** COPYRIGHT 1996-2006, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY,
C *** GEORGIA INSTITUTE OF TECHNOLOGY
C *** WRITTEN BY ATHANASIOS NENES
C *** UPDATED BY CHRISTOS FOUNTOUKIS
C
C=======================================================================
C
      DOUBLE PRECISION FUNCTION SSH_FUNCH2A (X)
      INCLUDE 'isrpia.inc'
C
      COMMON /SSHSOLUT/ CHI1, CHI2, CHI3, CHI4, CHI5, CHI6, CHI7, CHI8,
     &               PSI1, PSI2, PSI3, PSI4, PSI5, PSI6, PSI7, PSI8,
     &               A1,   A2,   A3,   A4,   A5,   A6,   A7,   A8


C
C *** SETUP PARAMETERS ************************************************
C
      PSI6   = X
      PSI1   = CHI1
      PSI2   = ZERO
      PSI3   = ZERO
      PSI7   = CHI7
      PSI8   = CHI8 
      FRST   = .TRUE.
      CALAIN = .TRUE. 
C
C *** SOLVE EQUATIONS ; WITH ITERATIONS FOR ACTIVITY COEF. ************
C
      DO 10 I=1,NSWEEP
C
      A1  = XK5 *(WATER/GAMA(2))**3.0
      A4  = (XK2/XKW)*R*TEMP*(GAMA(10)/GAMA(5))**2.0
      A5  = XK4 *R*TEMP*(WATER/GAMA(10))**2.0
      A6  = XK3 *R*TEMP*(WATER/GAMA(11))**2.0
      A7  = XK8 *(WATER/GAMA(1))**2.0
      A8  = XK9 *(WATER/GAMA(3))**2.0
      A64 = (XK3*XK2/XKW)*(GAMA(10)/GAMA(5)/GAMA(11))**2.0
      A64 = A64*(R*TEMP*WATER)**2.0
      A9  = XK1*WATER/GAMA(7)*(GAMA(8)/GAMA(7))**2.
C
C  CALCULATE DISSOCIATION QUANTITIES
C
      PSI5 = CHI5*(PSI6+PSI7) - A6/A5*PSI8*(CHI6-PSI6-PSI3)
      PSI5 = PSI5/(A6/A5*(CHI6-PSI6-PSI3) + PSI6 + PSI7)
      PSI5 = MAX(PSI5, TINY)
C
      IF (W(3).GT.TINY .AND. WATER.GT.TINY) THEN  ! First try 3rd order soln
         BB   =-(CHI4 + PSI6 + PSI5 + 1.d0/A4)
         CC   = CHI4*(PSI5+PSI6)
         DD   = BB*BB-4.d0*CC
         PSI4 =0.5d0*(-BB - SQRT(DD))
         PSI4 = MIN(PSI4,CHI4)
      ELSE
         PSI4 = TINY
      ENDIF
C
      IF (CHI7.GT.TINY .AND. WATER.GT.TINY) THEN     ! NACL DISSOLUTION
         DIAK = (PSI8-PSI6)**2.D0 + 4.D0*A7
         PSI7 = 0.5D0*( -(PSI8+PSI6) + SQRT(DIAK) )
         PSI7 = MAX(MIN(PSI7, CHI7), ZERO)
      ENDIF
C
      IF (CHI8.GT.TINY .AND. WATER.GT.TINY) THEN     ! NANO3 DISSOLUTION
         DIAK = (PSI7-PSI5)**2.D0 + 4.D0*A8
         PSI8 = 0.5D0*( -(PSI7+PSI5) + SQRT(DIAK) )
         PSI8 = MAX(MIN(PSI8, CHI8), ZERO)
      ENDIF
C
      IF (CHI1.GT.TINY .AND. WATER.GT.TINY) THEN     ! NA2SO4 DISSOLUTION
         AA = PSI7+PSI8
         BB = AA*AA
         CC =-A1/4.D0
         CALL SSH_POLY3 (AA, BB, CC, PSI1, ISLV)
         IF (ISLV.EQ.0) THEN
             PSI1 = MIN (PSI1, CHI1)
         ELSE
             PSI1 = ZERO
         ENDIF
      ENDIF
C
C *** CALCULATE SPECIATION ********************************************
C
      MOLAL (2) = PSI8 + PSI7 + 2.D0*PSI1                 ! NAI
      MOLAL (3) = PSI4                                    ! NH4I
      MOLAL (4) = PSI6 + PSI7                             ! CLI
      MOLAL (5) = PSI2 + PSI1                             ! SO4I
      MOLAL (6) = ZERO                                    ! HSO4I
      MOLAL (7) = PSI5 + PSI8                             ! NO3I
C
C *** FOR INTERACTION WITH HYDROPHILIC ORGANIC COMPOUNDS
      SMIN = 2.d0*MOLAL(5)+MOLAL(7)+MOLAL(4)-MOLAL(2)-MOLAL(3)
     &     + ORGANION
      CALL SSH_CALCPH (SMIN, HI, OHI)
      MOLAL (1) = HI
C 
      GNH3      = MAX(CHI4 - PSI4, TINY)
      GHNO3     = MAX(CHI5 - PSI5, TINY)
      GHCL      = MAX(CHI6 - PSI6, TINY)
C
      CNH42S4   = ZERO
      CNH4NO3   = ZERO
      CNACL     = MAX(CHI7 - PSI7, ZERO)
      CNANO3    = MAX(CHI8 - PSI8, ZERO)
      CNA2SO4   = MAX(CHI1 - PSI1, ZERO) 
C      
C *** NH4Cl(s) calculations
C
      A3   = XK6 /(R*TEMP*R*TEMP)
      DELT = MIN(GNH3, GHCL)
      BB = -(GNH3+GHCL)
      CC = GNH3*GHCL-A3
      DD = BB*BB - 4.D0*CC
      PSI31 = 0.5D0*(-BB + SQRT(DD))
      PSI32 = 0.5D0*(-BB - SQRT(DD))
      IF (DELT-PSI31.GT.ZERO .AND. PSI31.GT.ZERO) THEN
         PSI3 = PSI31
      ELSEIF (DELT-PSI32.GT.ZERO .AND. PSI32.GT.ZERO) THEN
         PSI3 = PSI32
      ELSE
         PSI3 = ZERO
      ENDIF
C 
C *** CALCULATE GAS / SOLID SPECIES (LIQUID IN MOLAL ALREADY) *********
C
      GNH3    = MAX(GNH3 - PSI3, TINY)
      GHCL    = MAX(GHCL - PSI3, TINY)
      CNH4CL  = PSI3
C
      CALL SSH_CALCMR                        ! Water content
C
C *** CALCULATE ACTIVITIES OR TERMINATE INTERNAL LOOP *****************
C
      IF (FRST.AND.CALAOU .OR. .NOT.FRST.AND.CALAIN) THEN
         CALL SSH_CALCACT     
      ELSE
         GOTO 20
      ENDIF
10    CONTINUE
C
C *** CALCULATE FUNCTION VALUE FOR OUTER LOOP ***************************
C
20    SSH_FUNCH2A = MOLAL(3)*MOLAL(4)/GHCL/GNH3/A64 - ONE
C
      RETURN
C
C *** END OF FUNCTION SSH_FUNCH2A *******************************************
C
      END


C=======================================================================
C
C *** ISORROPIA CODE
C *** SUBROUTINE SSH_CALCH1
C *** CASE H1
C
C     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
C     1. SULFATE POOR (SULRAT > 2.0) ; SODIUM RICH (SODRAT >= 2.0)
C     2. SOLID AEROSOL ONLY
C     3. SOLIDS POSSIBLE : NH4NO3, NH4CL, NA2SO4
C
C     THERE ARE TWO POSSIBLE REGIMES HERE, DEPENDING ON RELATIVE HUMIDITY:
C     1. WHEN RH >= MDRH ; LIQUID PHASE POSSIBLE (MDRH REGION)
C     2. WHEN RH < MDRH  ; ONLY SOLID PHASE POSSIBLE (SUBROUTINE SSH_CALCH1A)
C
C *** COPYRIGHT 1996-2006, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY,
C *** GEORGIA INSTITUTE OF TECHNOLOGY
C *** WRITTEN BY ATHANASIOS NENES
C *** UPDATED BY CHRISTOS FOUNTOUKIS
C
C=======================================================================
C
      SUBROUTINE SSH_CALCH1
      INCLUDE 'isrpia.inc'
      EXTERNAL SSH_CALCH1A, SSH_CALCH2A
C
C *** REGIME DEPENDS UPON THE AMBIENT RELATIVE HUMIDITY *****************
C
      IF (RH.LT.DRMH1) THEN    
         SCASE = 'H1 ; SUBCASE 1'  
         CALL SSH_CALCH1A              ! SOLID PHASE ONLY POSSIBLE
         SCASE = 'H1 ; SUBCASE 1'
      ELSE
         SCASE = 'H1 ; SUBCASE 2'  ! LIQUID & SOLID PHASE POSSIBLE
         CALL SSH_CALCMDRH (RH, DRMH1, DRNH4NO3, SSH_CALCH1A,
     c        SSH_CALCH2A)
         SCASE = 'H1 ; SUBCASE 2'
      ENDIF
C 
      RETURN
C
C *** END OF SUBROUTINE SSH_CALCH1 ******************************************
C
      END


C=======================================================================
C
C *** ISORROPIA CODE
C *** SUBROUTINE SSH_CALCH1A
C *** CASE H1 ; SUBCASE 1
C
C     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
C     1. SULFATE POOR (SULRAT > 2.0) ; SODIUM RICH (SODRAT >= 2.0)
C     2. SOLID AEROSOL ONLY
C     3. SOLIDS POSSIBLE : NH4NO3, NH4CL, NANO3, NA2SO4
C
C *** COPYRIGHT 1996-2006, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY,
C *** GEORGIA INSTITUTE OF TECHNOLOGY
C *** WRITTEN BY ATHANASIOS NENES
C *** UPDATED BY CHRISTOS FOUNTOUKIS
C
C=======================================================================
C
      SUBROUTINE SSH_CALCH1A
      INCLUDE 'isrpia.inc'
      DOUBLE PRECISION LAMDA, LAMDA1, LAMDA2, KAPA, KAPA1, KAPA2, NAFR,
     &                 NO3FR
C
C *** CALCULATE NON VOLATILE SOLIDS ***********************************
C
      CNA2SO4 = W(2)
      CNH42S4 = ZERO
      NAFR    = MAX (W(1)-2*CNA2SO4, ZERO)
      CNANO3  = MIN (NAFR, W(4))
      NO3FR   = MAX (W(4)-CNANO3, ZERO)
      CNACL   = MIN (MAX(NAFR-CNANO3, ZERO), W(5))
      CLFR    = MAX (W(5)-CNACL, ZERO)
C
C *** CALCULATE VOLATILE SPECIES **************************************
C
      ALF     = W(3)                     ! FREE NH3
      BET     = CLFR                     ! FREE CL
      GAM     = NO3FR                    ! FREE NO3
C
      RTSQ    = R*TEMP*R*TEMP
      A1      = XK6/RTSQ
      A2      = XK10/RTSQ
C
      THETA1  = GAM - BET*(A2/A1)
      THETA2  = A2/A1
C
C QUADRATIC EQUATION SOLUTION
C
      BB      = (THETA1-ALF-BET*(ONE+THETA2))/(ONE+THETA2)
      CC      = (ALF*BET-A1-BET*THETA1)/(ONE+THETA2)
      DD      = BB*BB - 4.0D0*CC
      IF (DD.LT.ZERO) GOTO 100   ! Solve each reaction seperately
C
C TWO ROOTS FOR KAPA, CHECK AND SEE IF ANY VALID
C
      SQDD    = SQRT(DD)
      KAPA1   = 0.5D0*(-BB+SQDD)
      KAPA2   = 0.5D0*(-BB-SQDD)
      LAMDA1  = THETA1 + THETA2*KAPA1
      LAMDA2  = THETA1 + THETA2*KAPA2
C
      IF (KAPA1.GE.ZERO .AND. LAMDA1.GE.ZERO) THEN
         IF (ALF-KAPA1-LAMDA1.GE.ZERO .AND.
     &       BET-KAPA1.GE.ZERO .AND. GAM-LAMDA1.GE.ZERO) THEN
             KAPA = KAPA1
             LAMDA= LAMDA1
             GOTO 200
         ENDIF
      ENDIF
C
      IF (KAPA2.GE.ZERO .AND. LAMDA2.GE.ZERO) THEN
         IF (ALF-KAPA2-LAMDA2.GE.ZERO .AND. 
     &       BET-KAPA2.GE.ZERO .AND. GAM-LAMDA2.GE.ZERO) THEN
             KAPA = KAPA2
             LAMDA= LAMDA2
             GOTO 200
         ENDIF
      ENDIF
C
C SEPERATE SOLUTION OF NH4CL & NH4NO3 EQUILIBRIA 
C 
100   KAPA  = ZERO
      LAMDA = ZERO
      DD1   = (ALF+BET)*(ALF+BET) - 4.0D0*(ALF*BET-A1)
      DD2   = (ALF+GAM)*(ALF+GAM) - 4.0D0*(ALF*GAM-A2)
C
C NH4CL EQUILIBRIUM
C
      IF (DD1.GE.ZERO) THEN
         SQDD1 = SQRT(DD1)
         KAPA1 = 0.5D0*(ALF+BET + SQDD1)
         KAPA2 = 0.5D0*(ALF+BET - SQDD1)
C
         IF (KAPA1.GE.ZERO .AND. KAPA1.LE.MIN(ALF,BET)) THEN
            KAPA = KAPA1 
         ELSE IF (KAPA2.GE.ZERO .AND. KAPA2.LE.MIN(ALF,BET)) THEN
            KAPA = KAPA2
         ELSE
            KAPA = ZERO
         ENDIF
      ENDIF
C
C NH4NO3 EQUILIBRIUM
C
      IF (DD2.GE.ZERO) THEN
         SQDD2 = SQRT(DD2)
         LAMDA1= 0.5D0*(ALF+GAM + SQDD2)
         LAMDA2= 0.5D0*(ALF+GAM - SQDD2)
C
         IF (LAMDA1.GE.ZERO .AND. LAMDA1.LE.MIN(ALF,GAM)) THEN
            LAMDA = LAMDA1 
         ELSE IF (LAMDA2.GE.ZERO .AND. LAMDA2.LE.MIN(ALF,GAM)) THEN
            LAMDA = LAMDA2
         ELSE
            LAMDA = ZERO
         ENDIF
      ENDIF
C
C IF BOTH KAPA, LAMDA ARE > 0, THEN APPLY EXISTANCE CRITERION
C
      IF (KAPA.GT.ZERO .AND. LAMDA.GT.ZERO) THEN
         IF (BET .LT. LAMDA/THETA1) THEN
            KAPA = ZERO
         ELSE
            LAMDA= ZERO
         ENDIF
      ENDIF
C
C *** CALCULATE COMPOSITION OF VOLATILE SPECIES ***********************
C
200   CONTINUE
      CNH4NO3 = LAMDA
      CNH4CL  = KAPA
C
      GNH3    = ALF - KAPA - LAMDA
      GHNO3   = GAM - LAMDA
      GHCL    = BET - KAPA
C
      RETURN
C
C *** END OF SUBROUTINE SSH_CALCH1A *****************************************
C
      END
C=======================================================================
C
C *** ISORROPIA CODE
C *** SUBROUTINE SSH_CALCI6
C *** CASE I6
C
C     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
C     1. SULFATE RICH, NO FREE ACID (1.0 <= SULRAT < 2.0)
C     2. SOLID & LIQUID AEROSOL POSSIBLE
C     3. SOLIDS POSSIBLE : (NH4)2SO4, NA2SO4
C
C *** COPYRIGHT 1996-2006, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY,
C *** GEORGIA INSTITUTE OF TECHNOLOGY
C *** WRITTEN BY ATHANASIOS NENES
C *** UPDATED BY CHRISTOS FOUNTOUKIS
C
C=======================================================================
C
      SUBROUTINE SSH_CALCI6
      INCLUDE 'isrpia.inc'
      COMMON /SSHSOLUT/ CHI1, CHI2, CHI3, CHI4, CHI5, CHI6, CHI7, CHI8,
     &               PSI1, PSI2, PSI3, PSI4, PSI5, PSI6, PSI7, PSI8,
     &               A1,   A2,   A3,   A4,   A5,   A6,   A7,   A8


C
C *** FIND DRY COMPOSITION **********************************************
C
      CALL SSH_CALCI1A
C
C *** SETUP PARAMETERS ************************************************
C
      CHI1 = CNH4HS4               ! Save from CALCI1 run
      CHI2 = CLC    
      CHI3 = CNAHSO4
      CHI4 = CNA2SO4
      CHI5 = CNH42S4
C
      PSI1 = CNH4HS4               ! ASSIGN INITIAL PSI's
      PSI2 = CLC   
      PSI3 = CNAHSO4
      PSI4 = CNA2SO4
      PSI5 = CNH42S4
C
      CALAOU = .TRUE.              ! Outer loop activity calculation flag
      FRST   = .TRUE.
      CALAIN = .TRUE.
C
C *** SOLVE EQUATIONS ; WITH ITERATIONS FOR ACTIVITY COEF. ************
C
      DO 10 I=1,NSWEEP
C
      A6 = XK1 *WATER/GAMA(7)*(GAMA(8)/GAMA(7))**2.
C
C  CALCULATE DISSOCIATION QUANTITIES
C
C *** FOR INTERACTION WITH HYDROPHILIC ORGANIC COMPOUNDS
      BB   = PSI2 + PSI4 + PSI5 + A6 + ORGANION           ! PSI6
      CC   =-A6*(PSI2 + PSI3 + PSI1) + ORGANION*(PSI2+PSI4+PSI5)
      DD   = BB*BB - 4.D0*CC
      PSI6 = 0.5D0*(-BB + SQRT(DD))
C
C *** CALCULATE SPECIATION ********************************************
C
C *** FOR INTERACTION WITH HYDROPHILIC ORGANIC COMPOUNDS
      MOLAL (1) = PSI6 + ORGANION                         ! HI
      MOLAL (2) = 2.D0*PSI4 + PSI3                        ! NAI
      MOLAL (3) = 3.D0*PSI2 + 2.D0*PSI5 + PSI1            ! NH4I
      MOLAL (5) = PSI2 + PSI4 + PSI5 + PSI6               ! SO4I
      MOLAL (6) = PSI2 + PSI3 + PSI1 - PSI6               ! HSO4I
      CLC       = ZERO
      CNAHSO4   = ZERO
      CNA2SO4   = CHI4 - PSI4
      CNH42S4   = ZERO
      CNH4HS4   = ZERO
      CALL SSH_CALCMR                                         ! Water content
C
C *** CALCULATE ACTIVITIES OR TERMINATE INTERNAL LOOP *****************
C
      IF (FRST.AND.CALAOU .OR. .NOT.FRST.AND.CALAIN) THEN
         CALL SSH_CALCACT     
      ELSE
         GOTO 20
      ENDIF
10    CONTINUE
C 
20    RETURN
C
C *** END OF SUBROUTINE SSH_CALCI6 *****************************************
C
      END

C=======================================================================
C
C *** ISORROPIA CODE
C *** SUBROUTINE SSH_CALCI5
C *** CASE I5
C
C     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
C     1. SULFATE RICH, NO FREE ACID (1.0 <= SULRAT < 2.0)
C     2. SOLID & LIQUID AEROSOL POSSIBLE
C     3. SOLIDS POSSIBLE : (NH4)2SO4, NA2SO4
C
C *** COPYRIGHT 1996-2006, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY,
C *** GEORGIA INSTITUTE OF TECHNOLOGY
C *** WRITTEN BY ATHANASIOS NENES
C *** UPDATED BY CHRISTOS FOUNTOUKIS
C
C=======================================================================
C
      SUBROUTINE SSH_CALCI5
      INCLUDE 'isrpia.inc'
      COMMON /SSHSOLUT/ CHI1, CHI2, CHI3, CHI4, CHI5, CHI6, CHI7, CHI8,
     &               PSI1, PSI2, PSI3, PSI4, PSI5, PSI6, PSI7, PSI8,
     &               A1,   A2,   A3,   A4,   A5,   A6,   A7,   A8


C
C *** FIND DRY COMPOSITION **********************************************
C
      CALL SSH_CALCI1A
C
C *** SETUP PARAMETERS ************************************************
C
      CHI1 = CNH4HS4               ! Save from CALCI1 run
      CHI2 = CLC    
      CHI3 = CNAHSO4
      CHI4 = CNA2SO4
      CHI5 = CNH42S4
C
      PSI1 = CNH4HS4               ! ASSIGN INITIAL PSI's
      PSI2 = CLC   
      PSI3 = CNAHSO4
      PSI4 = ZERO
      PSI5 = CNH42S4
C
      CALAOU =.TRUE.               ! Outer loop activity calculation flag
      PSI4LO = ZERO                ! Low  limit
      PSI4HI = CHI4                ! High limit
C    
C *** IF NA2SO4(S) =0, CALL SSH_SSH_FUNCI5B FOR Y4=0 ***************************
C
      IF (CHI4.LE.TINY) THEN
         Y1 = SSH_FUNCI5A (ZERO)
         GOTO 50
      ENDIF
C
C *** INITIAL VALUES FOR BISECTION ************************************
C
      X1 = PSI4HI
      Y1 = SSH_FUNCI5A (X1)
      YHI= Y1                      ! Save Y-value at HI position
C
C *** YHI < 0.0 THE SOLUTION IS ALWAYS UNDERSATURATED WITH NA2SO4 **
C
      IF (ABS(Y1).LE.EPS .OR. YHI.LT.ZERO) GOTO 50
C
C *** ROOT TRACKING ; FOR THE RANGE OF HI AND LO **********************
C
      DX = (PSI4HI-PSI4LO)/FLOAT(NDIV)
      DO 10 I=1,NDIV
         X2 = X1-DX
         Y2 = SSH_FUNCI5A (X2)
         IF (SIGN(1.d0,Y1)*SIGN(1.d0,Y2).LT.ZERO) GOTO 20  ! (Y1*Y2.LT.ZERO)
         X1 = X2
         Y1 = Y2
10    CONTINUE
C
C *** { YLO, YHI } > 0.0 THE SOLUTION IS ALWAYS SUPERSATURATED WITH NH4CL  
C
      YLO= Y1                      ! Save Y-value at Hi position
      IF (YLO.GT.ZERO .AND. YHI.GT.ZERO) THEN
         Y3 = SSH_FUNCI5A (ZERO)
         GOTO 50
      ELSE IF (ABS(Y2) .LT. EPS) THEN   ! X2 IS A SOLUTION 
         GOTO 50
      ELSE
         CALL SSH_PUSHERR (0001, 'CALCI5')    ! WARNING ERROR: NO SOLUTION
         GOTO 50
      ENDIF
C
C *** PERFORM BISECTION ***********************************************
C
20    DO 30 I=1,MAXIT
         X3 = 0.5*(X1+X2)
         Y3 = SSH_FUNCI5A (X3)
         IF (SIGN(1.d0,Y1)*SIGN(1.d0,Y3) .LE. ZERO) THEN  ! (Y1*Y3 .LE. ZERO)
            Y2    = Y3
            X2    = X3
         ELSE
            Y1    = Y3
            X1    = X3
         ENDIF
         IF (ABS(X2-X1) .LE. EPS*X1) GOTO 40
30    CONTINUE
      CALL SSH_PUSHERR (0002, 'CALCI5')    ! WARNING ERROR: NO CONVERGENCE
C
C *** CONVERGED ; RETURN **********************************************
C
40    X3 = 0.5*(X1+X2)
      Y3 = SSH_FUNCI5A (X3)
C 
50    RETURN

C *** END OF SUBROUTINE SSH_CALCI5 *****************************************
C
      END




C=======================================================================
C
C *** ISORROPIA CODE
C *** SUBROUTINE SSH_FUNCI5A
C *** CASE I5
C
C     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
C     1. SULFATE RICH, NO FREE ACID (1.0 <= SULRAT < 2.0)
C     2. SOLID & LIQUID AEROSOL POSSIBLE
C     3. SOLIDS POSSIBLE : NA2SO4
C
C *** COPYRIGHT 1996-2006, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY,
C *** GEORGIA INSTITUTE OF TECHNOLOGY
C *** WRITTEN BY ATHANASIOS NENES
C *** UPDATED BY CHRISTOS FOUNTOUKIS
C
C=======================================================================
C
      DOUBLE PRECISION FUNCTION SSH_FUNCI5A (P4)
      INCLUDE 'isrpia.inc'
      COMMON /SSHSOLUT/ CHI1, CHI2, CHI3, CHI4, CHI5, CHI6, CHI7, CHI8,
     &               PSI1, PSI2, PSI3, PSI4, PSI5, PSI6, PSI7, PSI8,
     &               A1,   A2,   A3,   A4,   A5,   A6,   A7,   A8


C
C *** SETUP PARAMETERS ************************************************
C
      PSI4   = P4     ! PSI3 already assigned in FUNCI5A
      FRST   = .TRUE.
      CALAIN = .TRUE.
C
C *** SOLVE EQUATIONS ; WITH ITERATIONS FOR ACTIVITY COEF. ************
C
      DO 10 I=1,NSWEEP
C
      A4 = XK5 *(WATER/GAMA(2))**3.0
      A5 = XK7 *(WATER/GAMA(4))**3.0
      A6 = XK1 *WATER/GAMA(7)*(GAMA(8)/GAMA(7))**2.
C
C  CALCULATE DISSOCIATION QUANTITIES
C
C *** FOR INTERACTION WITH HYDROPHILIC ORGANIC COMPOUNDS
      BB   = PSI2 + PSI4 + PSI5 + A6 + ORGANION   ! PSI6
      CC   =-A6*(PSI2 + PSI3 + PSI1)+ORGANION*(PSI2+PSI4+PSI5)
      DD   = BB*BB - 4.D0*CC
      PSI6 = 0.5D0*(-BB + SQRT(DD))
C
C *** CALCULATE SPECIATION ********************************************
C
C *** FOR INTERACTION WITH HYDROPHILIC ORGANIC COMPOUNDS
      MOLAL (1) = PSI6 + ORGANION                 ! HI
      MOLAL (2) = 2.D0*PSI4 + PSI3                ! NAI
      MOLAL (3) = 3.D0*PSI2 + 2.D0*PSI5 + PSI1    ! NH4I
      MOLAL (5) = PSI2 + PSI4 + PSI5 + PSI6       ! SO4I
      MOLAL (6) = PSI2 + PSI3 + PSI1 - PSI6       ! HSO4I
      CLC       = ZERO
      CNAHSO4   = ZERO
      CNA2SO4   = CHI4 - PSI4
      CNH42S4   = ZERO
      CNH4HS4   = ZERO
      CALL SSH_CALCMR                                 ! Water content
C
C *** CALCULATE ACTIVITIES OR TERMINATE INTERNAL LOOP *****************
C
      IF (FRST.AND.CALAOU .OR. .NOT.FRST.AND.CALAIN) THEN
         CALL SSH_CALCACT     
      ELSE
         GOTO 20
      ENDIF
10    CONTINUE
C
C *** CALCULATE OBJECTIVE FUNCTION ************************************
C
20    A4     = XK5 *(WATER/GAMA(2))**3.0    
      SSH_FUNCI5A= MOLAL(5)*MOLAL(2)*MOLAL(2)/A4 - ONE
      RETURN
C
C *** END OF FUNCTION SSH_FUNCI5A ********************************************
C
      END
C=======================================================================
C
C *** ISORROPIA CODE
C *** SUBROUTINE SSH_CALCI4
C *** CASE I4
C
C     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
C     1. SULFATE RICH, NO FREE ACID (1.0 <= SULRAT < 2.0)
C     2. SOLID & LIQUID AEROSOL POSSIBLE
C     3. SOLIDS POSSIBLE : (NH4)2SO4, NA2SO4
C
C *** COPYRIGHT 1996-2006, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY,
C *** GEORGIA INSTITUTE OF TECHNOLOGY
C *** WRITTEN BY ATHANASIOS NENES
C *** UPDATED BY CHRISTOS FOUNTOUKIS
C
C=======================================================================
C
      SUBROUTINE SSH_CALCI4
      INCLUDE 'isrpia.inc'
      COMMON /SSHSOLUT/ CHI1, CHI2, CHI3, CHI4, CHI5, CHI6, CHI7, CHI8,
     &               PSI1, PSI2, PSI3, PSI4, PSI5, PSI6, PSI7, PSI8,
     &               A1,   A2,   A3,   A4,   A5,   A6,   A7,   A8


C
C *** FIND DRY COMPOSITION **********************************************
C
      CALL SSH_CALCI1A
C
C *** SETUP PARAMETERS ************************************************
C
      CHI1 = CNH4HS4               ! Save from CALCI1 run
      CHI2 = CLC    
      CHI3 = CNAHSO4
      CHI4 = CNA2SO4
      CHI5 = CNH42S4
C
      PSI1 = CNH4HS4               ! ASSIGN INITIAL PSI's
      PSI2 = CLC   
      PSI3 = CNAHSO4
      PSI4 = ZERO  
      PSI5 = ZERO
C
      CALAOU = .TRUE.              ! Outer loop activity calculation flag
      PSI4LO = ZERO                ! Low  limit
      PSI4HI = CHI4                ! High limit
C    
C *** IF NA2SO4(S) =0, CALL SSH_SSH_FUNCI4B FOR Y4=0 ***************************
C
      IF (CHI4.LE.TINY) THEN
         Y1 = SSH_FUNCI4A (ZERO)
         GOTO 50
      ENDIF
C
C *** INITIAL VALUES FOR BISECTION ************************************
C
      X1 = PSI4HI
      Y1 = SSH_FUNCI4A (X1)
      YHI= Y1                      ! Save Y-value at HI position
C
C *** YHI < 0.0 THE SOLUTION IS ALWAYS UNDERSATURATED WITH NA2SO4 **
C
      IF (ABS(Y1).LE.EPS .OR. YHI.LT.ZERO) GOTO 50
C
C *** ROOT TRACKING ; FOR THE RANGE OF HI AND LO **********************
C
      DX = (PSI4HI-PSI4LO)/FLOAT(NDIV)
      DO 10 I=1,NDIV
         X2 = X1-DX
         Y2 = SSH_FUNCI4A (X2)
         IF (SIGN(1.d0,Y1)*SIGN(1.d0,Y2).LT.ZERO) GOTO 20  ! (Y1*Y2.LT.ZERO)
         X1 = X2
         Y1 = Y2
10    CONTINUE
C
C *** { YLO, YHI } > 0.0 THE SOLUTION IS ALWAYS SUPERSATURATED WITH NH4CL  
C
      YLO= Y1                      ! Save Y-value at Hi position
      IF (YLO.GT.ZERO .AND. YHI.GT.ZERO) THEN
         Y3 = SSH_FUNCI4A (ZERO)
         GOTO 50
      ELSE IF (ABS(Y2) .LT. EPS) THEN   ! X2 IS A SOLUTION 
         GOTO 50
      ELSE
         CALL SSH_PUSHERR (0001, 'CALCI4')    ! WARNING ERROR: NO SOLUTION
         GOTO 50
      ENDIF
C
C *** PERFORM BISECTION ***********************************************
C
20    DO 30 I=1,MAXIT
         X3 = 0.5*(X1+X2)
         Y3 = SSH_FUNCI4A (X3)
         IF (SIGN(1.d0,Y1)*SIGN(1.d0,Y3) .LE. ZERO) THEN  ! (Y1*Y3 .LE. ZERO)
            Y2    = Y3
            X2    = X3
         ELSE
            Y1    = Y3
            X1    = X3
         ENDIF
         IF (ABS(X2-X1) .LE. EPS*X1) GOTO 40
30    CONTINUE
      CALL SSH_PUSHERR (0002, 'CALCI4')    ! WARNING ERROR: NO CONVERGENCE
C
C *** CONVERGED ; RETURN **********************************************
C
40    X3 = 0.5*(X1+X2)
      Y3 = SSH_FUNCI4A (X3)
C
50    RETURN

C *** END OF SUBROUTINE SSH_CALCI4 *****************************************
C
      END




C=======================================================================
C
C *** ISORROPIA CODE
C *** SUBROUTINE SSH_FUNCI4A
C *** CASE I4
C
C     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
C     1. SULFATE RICH, NO FREE ACID (1.0 <= SULRAT < 2.0)
C     2. SOLID & LIQUID AEROSOL POSSIBLE
C     3. SOLIDS POSSIBLE : (NH4)2SO4, NA2SO4
C
C *** COPYRIGHT 1996-2006, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY,
C *** GEORGIA INSTITUTE OF TECHNOLOGY
C *** WRITTEN BY ATHANASIOS NENES
C *** UPDATED BY CHRISTOS FOUNTOUKIS
C
C=======================================================================
C
      DOUBLE PRECISION FUNCTION SSH_FUNCI4A (P4)
      INCLUDE 'isrpia.inc'
      COMMON /SSHSOLUT/ CHI1, CHI2, CHI3, CHI4, CHI5, CHI6, CHI7, CHI8,
     &               PSI1, PSI2, PSI3, PSI4, PSI5, PSI6, PSI7, PSI8,
     &               A1,   A2,   A3,   A4,   A5,   A6,   A7,   A8


C
C *** SETUP PARAMETERS ************************************************
C
      PSI4   = P4     ! PSI3 already assigned in FUNCI4A
      FRST   = .TRUE.
      CALAIN = .TRUE.
C
C *** SOLVE EQUATIONS ; WITH ITERATIONS FOR ACTIVITY COEF. ************
C
      DO 10 I=1,NSWEEP
C
      A4 = XK5 *(WATER/GAMA(2))**3.0
      A5 = XK7 *(WATER/GAMA(4))**3.0
      A6 = XK1 *WATER/GAMA(7)*(GAMA(8)/GAMA(7))**2.
      A7 = SQRT(A4/A5)
C
C  CALCULATE DISSOCIATION QUANTITIES
C
C *** FOR INTERACTION WITH HYDROPHILIC ORGANIC COMPOUNDS
      BB   = PSI2 + PSI4 + PSI5 + A6 + ORGANION   ! PSI6
      CC   =-A6*(PSI2 + PSI3 + PSI1)+ORGANION*(PSI2+PSI4+PSI5)
      DD   = BB*BB - 4.D0*CC
      PSI6 = 0.5D0*(-BB + SQRT(DD))
C
      PSI5 = (PSI3 + 2.D0*PSI4 - A7*(3.D0*PSI2 + PSI1))/2.D0/A7 
      PSI5 = MAX (MIN (PSI5, CHI5), ZERO)
C
C *** CALCULATE SPECIATION ********************************************
C
C *** FOR INTERACTION WITH HYDROPHILIC ORGANIC COMPOUNDS
      MOLAL (1) = PSI6 + ORGANION                 ! HI
      MOLAL (2) = 2.D0*PSI4 + PSI3                ! NAI
      MOLAL (3) = 3.D0*PSI2 + 2.D0*PSI5 + PSI1    ! NH4I
      MOLAL (5) = PSI2 + PSI4 + PSI5 + PSI6       ! SO4I
      MOLAL (6) = PSI2 + PSI3 + PSI1 - PSI6       ! HSO4I
      CLC       = ZERO
      CNAHSO4   = ZERO
      CNA2SO4   = CHI4 - PSI4
      CNH42S4   = CHI5 - PSI5
      CNH4HS4   = ZERO
      CALL SSH_CALCMR                                 ! Water content
C
C *** CALCULATE ACTIVITIES OR TERMINATE INTERNAL LOOP *****************
C
      IF (FRST.AND.CALAOU .OR. .NOT.FRST.AND.CALAIN) THEN
         CALL SSH_CALCACT     
      ELSE
         GOTO 20
      ENDIF
10    CONTINUE
C
C *** CALCULATE OBJECTIVE FUNCTION ************************************
C
20    A4     = XK5 *(WATER/GAMA(2))**3.0    
      SSH_FUNCI4A= MOLAL(5)*MOLAL(2)*MOLAL(2)/A4 - ONE
      RETURN
C
C *** END OF FUNCTION SSH_FUNCI4A ********************************************
C
      END
C=======================================================================
C
C *** ISORROPIA CODE
C *** SUBROUTINE SSH_CALCI3
C *** CASE I3
C
C     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
C     1. SULFATE RICH, NO FREE ACID (1.0 <= SULRAT < 2.0)
C     2. SOLID & LIQUID AEROSOL POSSIBLE
C     3. SOLIDS POSSIBLE : (NH4)2SO4, NA2SO4, NAHSO4, LC
C
C     THERE ARE THREE REGIMES IN THIS CASE:
C     1.(NA,NH4)HSO4(s) POSSIBLE. LIQUID & SOLID AEROSOL (SUBROUTINE SSH_CALCI3A)
C     2.(NA,NH4)HSO4(s) NOT POSSIBLE, AND RH < MDRH. SOLID AEROSOL ONLY 
C     3.(NA,NH4)HSO4(s) NOT POSSIBLE, AND RH >= MDRH. SOLID & LIQUID AEROSOL 
C
C     REGIMES 2. AND 3. ARE CONSIDERED TO BE THE SAME AS CASES I1A, I2B
C     RESPECTIVELY
C
C *** COPYRIGHT 1996-2006, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY,
C *** GEORGIA INSTITUTE OF TECHNOLOGY
C *** WRITTEN BY ATHANASIOS NENES
C *** UPDATED BY CHRISTOS FOUNTOUKIS
C
C=======================================================================
C
      SUBROUTINE SSH_CALCI3
      INCLUDE 'isrpia.inc'
      EXTERNAL SSH_CALCI1A, SSH_CALCI4
C
C *** FIND DRY COMPOSITION **********************************************
C
      CALL SSH_CALCI1A
C
C *** REGIME DEPENDS UPON THE POSSIBLE SOLIDS & RH **********************
C
      IF (CNH4HS4.GT.TINY .OR. CNAHSO4.GT.TINY) THEN
         SCASE = 'I3 ; SUBCASE 1'  
         CALL SSH_CALCI3A                     ! FULL SOLUTION
         SCASE = 'I3 ; SUBCASE 1'  
      ENDIF
C
      IF (WATER.LE.TINY) THEN
         IF (RH.LT.DRMI3) THEN         ! SOLID SOLUTION
            WATER = TINY
            DO 10 I=1,NIONS
               MOLAL(I) = ZERO
10          CONTINUE
            CALL SSH_CALCI1A
            SCASE = 'I3 ; SUBCASE 2'  
C
         ELSEIF (RH.GE.DRMI3) THEN     ! MDRH OF I3
            SCASE = 'I3 ; SUBCASE 3'
            CALL SSH_CALCMDRH (RH, DRMI3, DRLC, SSH_CALCI1A, SSH_CALCI4)
            SCASE = 'I3 ; SUBCASE 3'
         ENDIF
      ENDIF
C 
      RETURN
C
C *** END OF SUBROUTINE SSH_CALCI3 ******************************************
C
      END



C=======================================================================
C
C *** ISORROPIA CODE
C *** SUBROUTINE SSH_CALCI3A
C *** CASE I3 ; SUBCASE 1
C
C     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
C     1. SULFATE RICH, NO FREE ACID (1.0 <= SULRAT < 2.0)
C     2. SOLID & LIQUID AEROSOL POSSIBLE
C     3. SOLIDS POSSIBLE : (NH4)2SO4, NA2SO4, LC
C
C *** COPYRIGHT 1996-2006, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY,
C *** GEORGIA INSTITUTE OF TECHNOLOGY
C *** WRITTEN BY ATHANASIOS NENES
C *** UPDATED BY CHRISTOS FOUNTOUKIS
C
C=======================================================================
C
      SUBROUTINE SSH_CALCI3A
      INCLUDE 'isrpia.inc'
      COMMON /SSHSOLUT/ CHI1, CHI2, CHI3, CHI4, CHI5, CHI6, CHI7, CHI8,
     &               PSI1, PSI2, PSI3, PSI4, PSI5, PSI6, PSI7, PSI8,
     &               A1,   A2,   A3,   A4,   A5,   A6,   A7,   A8


C
C *** FIND DRY COMPOSITION **********************************************
C
      CALL SSH_CALCI1A         ! Needed when called from CALCMDRH
C
C *** SETUP PARAMETERS ************************************************
C
      CHI1 = CNH4HS4               ! Save from CALCI1 run
      CHI2 = CLC    
      CHI3 = CNAHSO4
      CHI4 = CNA2SO4
      CHI5 = CNH42S4
C
      PSI1 = CNH4HS4               ! ASSIGN INITIAL PSI's
      PSI2 = ZERO   
      PSI3 = CNAHSO4
      PSI4 = ZERO  
      PSI5 = ZERO
C
      CALAOU = .TRUE.              ! Outer loop activity calculation flag
      PSI2LO = ZERO                ! Low  limit
      PSI2HI = CHI2                ! High limit
C
C *** INITIAL VALUES FOR BISECTION ************************************
C
      X1 = PSI2HI
      Y1 = SSH_FUNCI3A (X1)
      YHI= Y1                      ! Save Y-value at HI position
C
C *** YHI < 0.0 THE SOLUTION IS ALWAYS UNDERSATURATED WITH LC *********
C
      IF (YHI.LT.EPS) GOTO 50
C
C *** ROOT TRACKING ; FOR THE RANGE OF HI AND LO **********************
C
      DX = (PSI2HI-PSI2LO)/FLOAT(NDIV)
      DO 10 I=1,NDIV
         X2 = MAX(X1-DX, PSI2LO)
         Y2 = SSH_FUNCI3A (X2)
         IF (SIGN(1.d0,Y1)*SIGN(1.d0,Y2).LT.ZERO) GOTO 20  ! (Y1*Y2.LT.ZERO)
         X1 = X2
         Y1 = Y2
10    CONTINUE
C
C *** { YLO, YHI } > 0.0 THE SOLUTION IS ALWAYS SUPERSATURATED WITH LC  
C
      IF (Y2.GT.EPS) Y2 = SSH_FUNCI3A (ZERO)
      GOTO 50
C
C *** PERFORM BISECTION ***********************************************
C
20    DO 30 I=1,MAXIT
         X3 = 0.5*(X1+X2)
         Y3 = SSH_FUNCI3A (X3)
         IF (SIGN(1.d0,Y1)*SIGN(1.d0,Y3) .LE. ZERO) THEN  ! (Y1*Y3 .LE. ZERO)
            Y2    = Y3
            X2    = X3
         ELSE
            Y1    = Y3
            X1    = X3
         ENDIF
         IF (ABS(X2-X1) .LE. EPS*X1) GOTO 40
30    CONTINUE
      CALL SSH_PUSHERR (0002, 'CALCI3A')    ! WARNING ERROR: NO CONVERGENCE
C
C *** CONVERGED ; RETURN **********************************************
C
40    X3 = 0.5*(X1+X2)
      Y3 = SSH_FUNCI3A (X3)
C 
50    RETURN

C *** END OF SUBROUTINE SSH_CALCI3A *****************************************
C
      END
C
C=======================================================================
C
C *** ISORROPIA CODE
C *** SUBROUTINE SSH_FUNCI3A
C *** CASE I3 ; SUBCASE 1
C
C     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
C     1. SULFATE RICH, NO FREE ACID (1.0 <= SULRAT < 2.0)
C     2. SOLID & LIQUID AEROSOL POSSIBLE
C     3. SOLIDS POSSIBLE : (NH4)2SO4, NA2SO4, LC
C
C *** COPYRIGHT 1996-2006, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY,
C *** GEORGIA INSTITUTE OF TECHNOLOGY
C *** WRITTEN BY ATHANASIOS NENES
C *** UPDATED BY CHRISTOS FOUNTOUKIS
C
C=======================================================================
C
      DOUBLE PRECISION FUNCTION SSH_FUNCI3A (P2)
      INCLUDE 'isrpia.inc'
      COMMON /SSHSOLUT/ CHI1, CHI2, CHI3, CHI4, CHI5, CHI6, CHI7, CHI8,
     &               PSI1, PSI2, PSI3, PSI4, PSI5, PSI6, PSI7, PSI8,
     &               A1,   A2,   A3,   A4,   A5,   A6,   A7,   A8


C
C *** SETUP PARAMETERS ************************************************
C
      PSI2   = P2                  ! Save PSI2 in COMMON BLOCK
      PSI4LO = ZERO                ! Low  limit for PSI4
      PSI4HI = CHI4                ! High limit for PSI4

C    
C *** IF NH3 =0, CALL SSH_SSH_FUNCI3B FOR Y4=0 ********************************
C
      IF (CHI4.LE.TINY) THEN
         SSH_FUNCI3A = SSH_FUNCI3B (ZERO)
         GOTO 50
      ENDIF
C
C *** INITIAL VALUES FOR BISECTION ************************************
C
      X1 = PSI4HI
      Y1 = SSH_FUNCI3B (X1)
      IF (ABS(Y1).LE.EPS) GOTO 50
      YHI= Y1                      ! Save Y-value at HI position
C
C *** YHI < 0.0 THE SOLUTION IS ALWAYS UNDERSATURATED WITH NA2SO4 *****
C
      IF (YHI.LT.ZERO) GOTO 50
C
C *** ROOT TRACKING ; FOR THE RANGE OF HI AND LO **********************
C
      DX = (PSI4HI-PSI4LO)/FLOAT(NDIV)
      DO 10 I=1,NDIV
         X2 = MAX(X1-DX, PSI4LO)
         Y2 = SSH_FUNCI3B (X2)
         IF (SIGN(1.d0,Y1)*SIGN(1.d0,Y2).LT.ZERO) GOTO 20  ! (Y1*Y2.LT.ZERO)
         X1 = X2
         Y1 = Y2
10    CONTINUE
C
C *** { YLO, YHI } > 0.0 THE SOLUTION IS ALWAYS SUPERSATURATED WITH NA2SO4
C
      IF (Y2.GT.EPS) Y2 = SSH_FUNCI3B (PSI4LO)
      GOTO 50
C
C *** PERFORM BISECTION ***********************************************
C
20    DO 30 I=1,MAXIT
         X3 = 0.5*(X1+X2)
         Y3 = SSH_FUNCI3B (X3)
         IF (SIGN(1.d0,Y1)*SIGN(1.d0,Y3) .LE. ZERO) THEN  ! (Y1*Y3 .LE. ZERO)
            Y2    = Y3
            X2    = X3
         ELSE
            Y1    = Y3
            X1    = X3
         ENDIF
         IF (ABS(X2-X1) .LE. EPS*X1) GOTO 40
30    CONTINUE
      CALL SSH_PUSHERR (0004, 'SSH_FUNCI3A')    ! WARNING ERROR: NO CONVERGENCE
C
C *** INNER LOOP CONVERGED **********************************************
C
40    X3 = 0.5*(X1+X2)
      Y3 = SSH_FUNCI3B (X3)
C
C *** CALCULATE FUNCTION VALUE FOR OUTER LOOP ***************************
C
50    A2      = XK13*(WATER/GAMA(13))**5.0
      SSH_FUNCI3A = MOLAL(5)*MOLAL(6)*MOLAL(3)**3.D0/A2 - ONE
      RETURN
C
C *** END OF FUNCTION SSH_FUNCI3A *******************************************
C
      END



C=======================================================================
C
C *** ISORROPIA CODE
C *** FUNCTION SSH_FUNCI3B
C *** CASE I3 ; SUBCASE 2
C
C     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
C     1. SULFATE RICH, NO FREE ACID (1.0 <= SULRAT < 2.0)
C     2. SOLID & LIQUID AEROSOL POSSIBLE
C     3. SOLIDS POSSIBLE : (NH4)2SO4, NA2SO4, LC
C
C     SOLUTION IS SAVED IN COMMON BLOCK /CASE/
C
C=======================================================================

C
      DOUBLE PRECISION FUNCTION SSH_FUNCI3B (P4)
      INCLUDE 'isrpia.inc'
C
      COMMON /SSHSOLUT/ CHI1, CHI2, CHI3, CHI4, CHI5, CHI6, CHI7, CHI8,
     &               PSI1, PSI2, PSI3, PSI4, PSI5, PSI6, PSI7, PSI8,
     &               A1,   A2,   A3,   A4,   A5,   A6,   A7,   A8


C
C *** SETUP PARAMETERS ************************************************
C
      PSI4   = P4   
C
C *** SETUP PARAMETERS ************************************************
C
      FRST   = .TRUE.
      CALAIN = .TRUE.
C
C *** SOLVE EQUATIONS ; WITH ITERATIONS FOR ACTIVITY COEF. ************
C
      DO 10 I=1,NSWEEP
C
      A4 = XK5*(WATER/GAMA(2))**3.0
      A5 = XK7*(WATER/GAMA(4))**3.0
      A6 = XK1*WATER/GAMA(7)*(GAMA(8)/GAMA(7))**2.
      A7 = SQRT(A4/A5)
C
C  CALCULATE DISSOCIATION QUANTITIES
C
C *** FOR INTERACTION WITH HYDROPHILIC ORGANIC COMPOUNDS
      BB   = PSI2 + PSI4 + PSI5 + A6 + ORGANION        ! PSI6
      CC   =-A6*(PSI2 + PSI3 + PSI1)+ORGANION*(PSI2+PSI4+PSI5)
      DD   = BB*BB - 4.D0*CC
      PSI6 = 0.5D0*(-BB + SQRT(DD))
C
      PSI5 = (PSI3 + 2.D0*PSI4 - A7*(3.D0*PSI2 + PSI1))/2.D0/A7 
      PSI5 = MAX (MIN (PSI5, CHI5), ZERO)
C
C *** CALCULATE SPECIATION ********************************************
C
C *** FOR INTERACTION WITH HYDROPHILIC ORGANIC COMPOUNDS
      MOLAL(1) = PSI6 + ORGANION                       ! HI
      MOLAL(2) = 2.D0*PSI4 + PSI3                      ! NAI
      MOLAL(3) = 3.D0*PSI2 + 2.D0*PSI5 + PSI1          ! NH4I
      MOLAL(5) = PSI2 + PSI4 + PSI5 + PSI6             ! SO4I
      MOLAL(6) = MAX(PSI2 + PSI3 + PSI1 - PSI6, TINY)  ! HSO4I
      CLC      = MAX(CHI2 - PSI2, ZERO)
      CNAHSO4  = ZERO
      CNA2SO4  = MAX(CHI4 - PSI4, ZERO)
      CNH42S4  = MAX(CHI5 - PSI5, ZERO)
      CNH4HS4  = ZERO
      CALL SSH_CALCMR                                       ! Water content
C
C *** CALCULATE ACTIVITIES OR TERMINATE INTERNAL LOOP *****************
C
      IF (FRST.AND.CALAOU .OR. .NOT.FRST.AND.CALAIN) THEN
         CALL SSH_CALCACT     
      ELSE
         GOTO 20
      ENDIF
10    CONTINUE
C
C *** CALCULATE OBJECTIVE FUNCTION ************************************
C
20    A4     = XK5 *(WATER/GAMA(2))**3.0    
      SSH_FUNCI3B= MOLAL(5)*MOLAL(2)*MOLAL(2)/A4 - ONE
      RETURN
C
C *** END OF FUNCTION SSH_FUNCI3B ********************************************
C
      END
C=======================================================================
C
C *** ISORROPIA CODE
C *** SUBROUTINE SSH_CALCI2
C *** CASE I2
C
C     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
C     1. SULFATE RICH, NO FREE ACID (1.0 <= SULRAT < 2.0)
C     2. SOLID & LIQUID AEROSOL POSSIBLE
C     3. SOLIDS POSSIBLE : (NH4)2SO4, NA2SO4, NAHSO4, LC
C
C     THERE ARE THREE REGIMES IN THIS CASE:
C     1. NH4HSO4(s) POSSIBLE. LIQUID & SOLID AEROSOL (SUBROUTINE SSH_CALCI2A)
C     2. NH4HSO4(s) NOT POSSIBLE, AND RH < MDRH. SOLID AEROSOL ONLY 
C     3. NH4HSO4(s) NOT POSSIBLE, AND RH >= MDRH. SOLID & LIQUID AEROSOL 
C
C     REGIMES 2. AND 3. ARE CONSIDERED TO BE THE SAME AS CASES I1A, I2B
C     RESPECTIVELY
C
C *** COPYRIGHT 1996-2006, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY,
C *** GEORGIA INSTITUTE OF TECHNOLOGY
C *** WRITTEN BY ATHANASIOS NENES
C *** UPDATED BY CHRISTOS FOUNTOUKIS
C
C=======================================================================
C
      SUBROUTINE SSH_CALCI2
      INCLUDE 'isrpia.inc'
      EXTERNAL SSH_CALCI1A, SSH_CALCI3A
C
C *** FIND DRY COMPOSITION **********************************************
C
      CALL SSH_CALCI1A
C
C *** REGIME DEPENDS UPON THE POSSIBLE SOLIDS & RH **********************
C
      IF (CNH4HS4.GT.TINY) THEN
         SCASE = 'I2 ; SUBCASE 1'  
         CALL SSH_CALCI2A                       
         SCASE = 'I2 ; SUBCASE 1'  
      ENDIF
C
      IF (WATER.LE.TINY) THEN
         IF (RH.LT.DRMI2) THEN         ! SOLID SOLUTION ONLY
            WATER = TINY
            DO 10 I=1,NIONS
               MOLAL(I) = ZERO
10          CONTINUE
            CALL SSH_CALCI1A
            SCASE = 'I2 ; SUBCASE 2'  
C
         ELSEIF (RH.GE.DRMI2) THEN     ! MDRH OF I2
            SCASE = 'I2 ; SUBCASE 3'
            CALL SSH_CALCMDRH (RH, DRMI2, DRNAHSO4, SSH_CALCI1A,
     c           SSH_CALCI3A)
            SCASE = 'I2 ; SUBCASE 3'
         ENDIF
      ENDIF
C 
      RETURN
C
C *** END OF SUBROUTINE SSH_CALCI2 ******************************************
C
      END


C=======================================================================
C
C *** ISORROPIA CODE
C *** SUBROUTINE SSH_CALCI2A
C *** CASE I2 ; SUBCASE A
C
C     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
C     1. SULFATE RICH, NO FREE ACID (1.0 <= SULRAT < 2.0)
C     2. SOLID & LIQUID AEROSOL POSSIBLE
C     3. SOLIDS POSSIBLE : (NH4)2SO4, NA2SO4, NAHSO4, LC
C
C *** COPYRIGHT 1996-2006, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY,
C *** GEORGIA INSTITUTE OF TECHNOLOGY
C *** WRITTEN BY ATHANASIOS NENES
C *** UPDATED BY CHRISTOS FOUNTOUKIS
C
C=======================================================================
C
      SUBROUTINE SSH_CALCI2A
      INCLUDE 'isrpia.inc'
      COMMON /SSHSOLUT/ CHI1, CHI2, CHI3, CHI4, CHI5, CHI6, CHI7, CHI8,
     &               PSI1, PSI2, PSI3, PSI4, PSI5, PSI6, PSI7, PSI8,
     &               A1,   A2,   A3,   A4,   A5,   A6,   A7,   A8


C
C *** FIND DRY COMPOSITION **********************************************
C
      CALL SSH_CALCI1A    ! Needed when called from CALCMDRH
C
C *** SETUP PARAMETERS ************************************************
C
      CHI1 = CNH4HS4               ! Save from CALCI1 run
      CHI2 = CLC    
      CHI3 = CNAHSO4
      CHI4 = CNA2SO4
      CHI5 = CNH42S4
C
      PSI1 = CNH4HS4               ! ASSIGN INITIAL PSI's
      PSI2 = ZERO   
      PSI3 = ZERO   
      PSI4 = ZERO  
      PSI5 = ZERO
C
      CALAOU = .TRUE.              ! Outer loop activity calculation flag
      PSI2LO = ZERO                ! Low  limit
      PSI2HI = CHI2                ! High limit
C
C *** INITIAL VALUES FOR BISECTION ************************************
C
      X1 = PSI2HI
      Y1 = SSH_FUNCI2A (X1)
      YHI= Y1                      ! Save Y-value at HI position
C
C *** YHI < 0.0 THE SOLUTION IS ALWAYS UNDERSATURATED WITH LC *********
C
      IF (YHI.LT.EPS) GOTO 50
C
C *** ROOT TRACKING ; FOR THE RANGE OF HI AND LO **********************
C
      DX = (PSI2HI-PSI2LO)/FLOAT(NDIV)
      DO 10 I=1,NDIV
         X2 = MAX(X1-DX, PSI2LO)
         Y2 = SSH_FUNCI2A (X2)
         IF (SIGN(1.d0,Y1)*SIGN(1.d0,Y2).LT.ZERO) GOTO 20  ! (Y1*Y2.LT.ZERO)
         X1 = X2
         Y1 = Y2
10    CONTINUE
C
C *** { YLO, YHI } > 0.0 THE SOLUTION IS ALWAYS SUPERSATURATED WITH LC  
C
      IF (Y2.GT.EPS) Y2 = SSH_FUNCI2A (ZERO)
      GOTO 50
C
C *** PERFORM BISECTION ***********************************************
C
20    DO 30 I=1,MAXIT
         X3 = 0.5*(X1+X2)
         Y3 = SSH_FUNCI2A (X3)
         IF (SIGN(1.d0,Y1)*SIGN(1.d0,Y3) .LE. ZERO) THEN  ! (Y1*Y3 .LE. ZERO)
            Y2    = Y3
            X2    = X3
         ELSE
            Y1    = Y3
            X1    = X3
         ENDIF
         IF (ABS(X2-X1) .LE. EPS*X1) GOTO 40
30    CONTINUE
      CALL SSH_PUSHERR (0002, 'CALCI2A')    ! WARNING ERROR: NO CONVERGENCE
C
C *** CONVERGED ; RETURN **********************************************
C
40    X3 = 0.5*(X1+X2)
      Y3 = SSH_FUNCI2A (X3)
C
50    RETURN

C *** END OF SUBROUTINE SSH_CALCI2A *****************************************
C
      END




C=======================================================================
C
C *** ISORROPIA CODE
C *** SUBROUTINE SSH_FUNCI2A
C *** CASE I2 ; SUBCASE 1
C
C     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
C     1. SULFATE RICH, NO FREE ACID (1.0 <= SULRAT < 2.0)
C     2. SOLID & LIQUID AEROSOL POSSIBLE
C     3. SOLIDS POSSIBLE : (NH4)2SO4, NA2SO4, NAHSO4, LC
C
C *** COPYRIGHT 1996-2006, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY,
C *** GEORGIA INSTITUTE OF TECHNOLOGY
C *** WRITTEN BY ATHANASIOS NENES
C *** UPDATED BY CHRISTOS FOUNTOUKIS
C
C=======================================================================
C
      DOUBLE PRECISION FUNCTION SSH_FUNCI2A (P2)
      INCLUDE 'isrpia.inc'
      COMMON /SSHSOLUT/ CHI1, CHI2, CHI3, CHI4, CHI5, CHI6, CHI7, CHI8,
     &               PSI1, PSI2, PSI3, PSI4, PSI5, PSI6, PSI7, PSI8,
     &               A1,   A2,   A3,   A4,   A5,   A6,   A7,   A8


C
C *** SETUP PARAMETERS ************************************************
C
      FRST   = .TRUE.
      CALAIN = .TRUE.
      PSI2   = P2                  ! Save PSI2 in COMMON BLOCK
      PSI3   = CHI3
      PSI4   = CHI4

      PSI5   = CHI5
      PSI6   = ZERO
C
C *** SOLVE EQUATIONS ; WITH ITERATIONS FOR ACTIVITY COEF. ************
C
      DO 10 I=1,NSWEEP
C
      A3 = XK11*(WATER/GAMA(12))**2.0
      A4 = XK5 *(WATER/GAMA(2))**3.0
      A5 = XK7 *(WATER/GAMA(4))**3.0
      A6 = XK1 *WATER/GAMA(7)*(GAMA(8)/GAMA(7))**2.
      A7 = SQRT(A4/A5)
C
C  CALCULATE DISSOCIATION QUANTITIES
C
      IF (CHI5.GT.TINY .AND. WATER.GT.TINY) THEN     
         PSI5 = (PSI3 + 2.D0*PSI4 - A7*(3.D0*PSI2 + PSI1))/2.D0/A7 
         PSI5 = MAX(MIN (PSI5, CHI5), TINY)
      ENDIF
C
      IF (CHI4.GT.TINY .AND. WATER.GT.TINY) THEN     
         AA   = PSI2+PSI5+PSI6+PSI3
         BB   = PSI3*AA
         CC   = 0.25D0*(PSI3*PSI3*(PSI2+PSI5+PSI6)-A4)
         CALL SSH_POLY3 (AA, BB, CC, PSI4, ISLV)
         IF (ISLV.EQ.0) THEN
            PSI4 = MIN (PSI4, CHI4)
         ELSE
            PSI4 = ZERO
         ENDIF
      ENDIF
C
      IF (CHI3.GT.TINY .AND. WATER.GT.TINY) THEN     
         AA   = 2.D0*PSI4 + PSI2 + PSI1 - PSI6
         BB   = 2.D0*PSI4*(PSI2 + PSI1 - PSI6) - A3
         CC   = ZERO
         CALL SSH_POLY3 (AA, BB, CC, PSI3, ISLV)
         IF (ISLV.EQ.0) THEN
            PSI3 = MIN (PSI3, CHI3)
         ELSE
            PSI3 = ZERO
         ENDIF
      ENDIF
C
C *** FOR INTERACTION WITH HYDROPHILIC ORGANIC COMPOUNDS
      BB   = PSI2 + PSI4 + PSI5 + A6 + ORGANION  ! PSI6
      CC   =-A6*(PSI2 + PSI3 + PSI1)+ORGANION*(PSI2+PSI4+PSI5)
      DD   = BB*BB - 4.D0*CC
      PSI6 = 0.5D0*(-BB + SQRT(DD))
C
C *** CALCULATE SPECIATION ********************************************
C
C *** FOR INTERACTION WITH HYDROPHILIC ORGANIC COMPOUNDS
      MOLAL (1) = PSI6 + ORGANION                ! HI
      MOLAL (2) = 2.D0*PSI4 + PSI3               ! NAI
      MOLAL (3) = 3.D0*PSI2 + 2.D0*PSI5 + PSI1   ! NH4I
      MOLAL (5) = PSI2 + PSI4 + PSI5 + PSI6      ! SO4I
      MOLAL (6) = PSI2 + PSI3 + PSI1 - PSI6      ! HSO4I
      CLC       = CHI2 - PSI2
      CNAHSO4   = CHI3 - PSI3
      CNA2SO4   = CHI4 - PSI4
      CNH42S4   = CHI5 - PSI5
      CNH4HS4   = ZERO
      CALL SSH_CALCMR                                ! Water content
C
C *** CALCULATE ACTIVITIES OR TERMINATE INTERNAL LOOP *****************
C
      IF (FRST.AND.CALAOU .OR. .NOT.FRST.AND.CALAIN) THEN
         CALL SSH_CALCACT     
      ELSE
         GOTO 20
      ENDIF
10    CONTINUE
C
C *** CALCULATE FUNCTION VALUE FOR OUTER LOOP ***************************
C
20    A2      = XK13*(WATER/GAMA(13))**5.0
      SSH_FUNCI2A = MOLAL(5)*MOLAL(6)*MOLAL(3)**3.D0/A2 - ONE
      

      RETURN
C
C *** END OF FUNCTION SSH_FUNCI2A *******************************************
C
      END

C=======================================================================
C
C *** ISORROPIA CODE
C *** SUBROUTINE SSH_CALCI1
C *** CASE I1
C
C     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
C     1. SULFATE RICH, NO FREE ACID (1.0 <= SULRAT < 2.0)
C     2. SOLID AEROSOL ONLY
C     3. SOLIDS POSSIBLE : NH4NO3, NH4CL, NA2SO4
C
C     THERE ARE TWO POSSIBLE REGIMES HERE, DEPENDING ON RELATIVE HUMIDITY:
C     1. WHEN RH >= MDRH ; LIQUID PHASE POSSIBLE (MDRH REGION)
C     2. WHEN RH < MDRH  ; ONLY SOLID PHASE POSSIBLE (SUBROUTINE SSH_CALCI1A)
C
C *** COPYRIGHT 1996-2006, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY,
C *** GEORGIA INSTITUTE OF TECHNOLOGY
C *** WRITTEN BY ATHANASIOS NENES
C *** UPDATED BY CHRISTOS FOUNTOUKIS
C
C=======================================================================
C
      SUBROUTINE SSH_CALCI1
      INCLUDE 'isrpia.inc'
      EXTERNAL SSH_CALCI1A, SSH_CALCI2A
C
C *** REGIME DEPENDS UPON THE AMBIENT RELATIVE HUMIDITY *****************
C
      IF (RH.LT.DRMI1) THEN    
         SCASE = 'I1 ; SUBCASE 1'  
         CALL SSH_CALCI1A              ! SOLID PHASE ONLY POSSIBLE
         SCASE = 'I1 ; SUBCASE 1'
      ELSE
         SCASE = 'I1 ; SUBCASE 2'  ! LIQUID & SOLID PHASE POSSIBLE
         CALL SSH_CALCMDRH (RH, DRMI1, DRNH4HS4, SSH_CALCI1A,
     c        SSH_CALCI2A)
         SCASE = 'I1 ; SUBCASE 2'
      ENDIF
C 
C *** AMMONIA IN GAS PHASE **********************************************
C
C      CALL SSH_CALCNH3
C 
      RETURN
C
C *** END OF SUBROUTINE SSH_CALCI1 ******************************************
C
      END


C=======================================================================
C
C *** ISORROPIA CODE
C *** SUBROUTINE SSH_CALCI1A
C *** CASE I1 ; SUBCASE 1
C
C     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
C     1. SULFATE RICH, NO FREE ACID (1.0 <= SULRAT < 2.0)
C     2. SOLID AEROSOL ONLY
C     3. SOLIDS POSSIBLE : NH4HSO4, NAHSO4, (NH4)2SO4, NA2SO4, LC
C
C *** COPYRIGHT 1996-2006, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY,
C *** GEORGIA INSTITUTE OF TECHNOLOGY
C *** WRITTEN BY ATHANASIOS NENES
C *** UPDATED BY CHRISTOS FOUNTOUKIS
C
C=======================================================================
C
      SUBROUTINE SSH_CALCI1A
      INCLUDE 'isrpia.inc'
C
C *** CALCULATE NON VOLATILE SOLIDS ***********************************
C
      CNA2SO4 = 0.5D0*W(1)
      CNH4HS4 = ZERO
      CNAHSO4 = ZERO
      CNH42S4 = ZERO
      FRSO4   = MAX(W(2)-CNA2SO4, ZERO)
C
      CLC     = MIN(W(3)/3.D0, FRSO4/2.D0)
      FRSO4   = MAX(FRSO4-2.D0*CLC, ZERO)
      FRNH4   = MAX(W(3)-3.D0*CLC,  ZERO)
C
      IF (FRSO4.LE.TINY) THEN
         CLC     = MAX(CLC - FRNH4, ZERO)
         CNH42S4 = 2.D0*FRNH4

      ELSEIF (FRNH4.LE.TINY) THEN
         CNH4HS4 = 3.D0*MIN(FRSO4, CLC)
         CLC     = MAX(CLC-FRSO4, ZERO)
         IF (CNA2SO4.GT.TINY) THEN
            FRSO4   = MAX(FRSO4-CNH4HS4/3.D0, ZERO)
            CNAHSO4 = 2.D0*FRSO4
            CNA2SO4 = MAX(CNA2SO4-FRSO4, ZERO)
         ENDIF
      ENDIF
C
C *** CALCULATE GAS SPECIES *********************************************
C
      GHNO3 = W(4)
      GHCL  = W(5)
      GNH3  = ZERO
C
      RETURN
C
C *** END OF SUBROUTINE SSH_CALCI1A *****************************************
C
      END
C=======================================================================
C
C *** ISORROPIA CODE
C *** SUBROUTINE SSH_CALCJ3
C *** CASE J3
C
C     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
C     1. SULFATE RICH, FREE ACID (SULRAT < 1.0)
C     2. THERE IS ONLY A LIQUID PHASE
C
C *** COPYRIGHT 1996-2006, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY,
C *** GEORGIA INSTITUTE OF TECHNOLOGY
C *** WRITTEN BY ATHANASIOS NENES
C *** UPDATED BY CHRISTOS FOUNTOUKIS
C
C=======================================================================
C
      SUBROUTINE SSH_CALCJ3
      INCLUDE 'isrpia.inc'
C
      DOUBLE PRECISION LAMDA, KAPA
C
C *** SETUP PARAMETERS ************************************************
C
      CALAOU = .TRUE.              ! Outer loop activity calculation flag
      FRST   = .TRUE.
      CALAIN = .TRUE.
C
      LAMDA  = MAX(W(2) - W(3) - W(1), TINY)  ! FREE H2SO4
      CHI1   = W(1)                           ! NA TOTAL as NaHSO4
      CHI2   = W(3)                           ! NH4 TOTAL as NH4HSO4
      PSI1   = CHI1
      PSI2   = CHI2                           ! ALL NH4HSO4 DELIQUESCED
C
C *** SOLVE EQUATIONS ; WITH ITERATIONS FOR ACTIVITY COEF. ************
C
      DO 10 I=1,NSWEEP
C
      A3 = XK1  *WATER/GAMA(7)*(GAMA(8)/GAMA(7))**2.0
C
C  CALCULATE DISSOCIATION QUANTITIES
C
C *** FOR INTERACTION WITH HYDROPHILIC ORGANIC COMPOUNDS
      BB   = A3+LAMDA + ORGANION               ! KAPA
      CC   =-A3*(LAMDA + PSI1 + PSI2)
      DD   = BB*BB-4.D0*CC
      KAPA = 0.5D0*(-BB+SQRT(DD))
C
C *** CALCULATE SPECIATION ********************************************
C
C *** FOR INTERACTION WITH HYDROPHILIC ORGANIC COMPOUNDS
      MOLAL (1) = LAMDA + KAPA + ORGANION      ! HI
      MOLAL (2) = PSI1                         ! NAI
      MOLAL (3) = PSI2                         ! NH4I
      MOLAL (4) = ZERO                         ! CLI
      MOLAL (5) = KAPA                         ! SO4I
      MOLAL (6) = LAMDA + PSI1 + PSI2 - KAPA   ! HSO4I
      MOLAL (7) = ZERO                         ! NO3I
C
      CNAHSO4   = ZERO
      CNH4HS4   = ZERO
C
      CALL SSH_CALCMR                              ! Water content
C
C *** CALCULATE ACTIVITIES OR TERMINATE INTERNAL LOOP *****************
C
      IF (FRST.AND.CALAOU .OR. .NOT.FRST.AND.CALAIN) THEN
         CALL SSH_CALCACT     
      ELSE
         GOTO 50
      ENDIF
10    CONTINUE
C 
50    RETURN
C
C *** END OF SUBROUTINE SSH_CALCJ3 ******************************************
C
      END
C=======================================================================
C
C *** ISORROPIA CODE
C *** SUBROUTINE SSH_CALCJ2
C *** CASE J2
C
C     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
C     1. SULFATE RICH, FREE ACID (SULRAT < 1.0)
C     2. THERE IS BOTH A LIQUID & SOLID PHASE
C     3. SOLIDS POSSIBLE : NAHSO4
C
C *** COPYRIGHT 1996-2006, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY,
C *** GEORGIA INSTITUTE OF TECHNOLOGY
C *** WRITTEN BY ATHANASIOS NENES
C *** UPDATED BY CHRISTOS FOUNTOUKIS
C
C=======================================================================
C
      SUBROUTINE SSH_CALCJ2
      INCLUDE 'isrpia.inc'
C
      DOUBLE PRECISION LAMDA, KAPA
      COMMON /SSHCASEJ/ CHI1, CHI2, CHI3, LAMDA, KAPA, PSI1, PSI2, PSI3, 
     &               A1,   A2,   A3


C
C *** SETUP PARAMETERS ************************************************
C
      CALAOU = .TRUE.              ! Outer loop activity calculation flag
      CHI1   = W(1)                ! NA TOTAL
      CHI2   = W(3)                ! NH4 TOTAL
      PSI1LO = TINY                ! Low  limit
      PSI1HI = CHI1                ! High limit
C
C *** INITIAL VALUES FOR BISECTION ************************************
C
      X1 = PSI1HI
      Y1 = SSH_FUNCJ2 (X1)
      YHI= Y1                      ! Save Y-value at HI position
C
C *** YHI < 0.0 THE SOLUTION IS ALWAYS UNDERSATURATED WITH NH42SO4 ****
C
      IF (ABS(Y1).LE.EPS .OR. YHI.LT.ZERO) GOTO 50
C
C *** ROOT TRACKING ; FOR THE RANGE OF HI AND LO **********************
C
      DX = (PSI1HI-PSI1LO)/FLOAT(NDIV)
      DO 10 I=1,NDIV
         X2 = X1-DX
         Y2 = SSH_FUNCJ2 (X2)
         IF (SIGN(1.d0,Y1)*SIGN(1.d0,Y2).LT.ZERO) GOTO 20  ! (Y1*Y2.LT.ZERO)
         X1 = X2
         Y1 = Y2
10    CONTINUE
C
C *** { YLO, YHI } > 0.0 THE SOLUTION IS ALWAYS SUPERSATURATED WITH NH42SO4
C
      YLO= Y1                      ! Save Y-value at Hi position
      IF (YLO.GT.ZERO .AND. YHI.GT.ZERO) THEN
         Y3 = SSH_FUNCJ2 (ZERO)
         GOTO 50
      ELSE IF (ABS(Y2) .LT. EPS) THEN   ! X2 IS A SOLUTION 
         GOTO 50
      ELSE
         CALL SSH_PUSHERR (0001, 'CALCJ2')    ! WARNING ERROR: NO SOLUTION
         GOTO 50
      ENDIF
C
C *** PERFORM BISECTION ***********************************************
C
20    DO 30 I=1,MAXIT
         X3 = 0.5*(X1+X2)
         Y3 = SSH_FUNCJ2 (X3)
         IF (SIGN(1.d0,Y1)*SIGN(1.d0,Y3) .LE. ZERO) THEN  ! (Y1*Y3 .LE. ZERO)
            Y2    = Y3
            X2    = X3
         ELSE
            Y1    = Y3
            X1    = X3
         ENDIF
         IF (ABS(X2-X1) .LE. EPS*X1) GOTO 40
30    CONTINUE
      CALL SSH_PUSHERR (0002, 'CALCJ2')    ! WARNING ERROR: NO CONVERGENCE
C
C *** CONVERGED ; RETURN **********************************************
C
40    X3 = 0.5*(X1+X2)
      Y3 = SSH_FUNCJ2 (X3)
C 
50    RETURN
C
C *** END OF SUBROUTINE SSH_CALCJ2 ******************************************
C
      END




C=======================================================================
C
C *** ISORROPIA CODE
C *** SUBROUTINE SSH_FUNCJ2
C *** CASE J2
C
C     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
C     1. SULFATE RICH, FREE ACID (SULRAT < 1.0)
C     2. THERE IS BOTH A LIQUID & SOLID PHASE
C     3. SOLIDS POSSIBLE : (NH4)2SO4, NH4CL, NA2SO4
C
C *** COPYRIGHT 1996-2006, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY,
C *** GEORGIA INSTITUTE OF TECHNOLOGY
C *** WRITTEN BY ATHANASIOS NENES
C *** UPDATED BY CHRISTOS FOUNTOUKIS
C
C=======================================================================
C
      DOUBLE PRECISION FUNCTION SSH_FUNCJ2 (P1)
      INCLUDE 'isrpia.inc'
C
      DOUBLE PRECISION LAMDA, KAPA
      COMMON /SSHCASEJ/ CHI1, CHI2, CHI3, LAMDA, KAPA, PSI1, PSI2, PSI3, 
     &               A1,   A2,   A3


C
C *** SETUP PARAMETERS ************************************************
C
      FRST   = .TRUE.
      CALAIN = .TRUE.
C
      LAMDA  = MAX(W(2) - W(3) - W(1), TINY)  ! FREE H2SO4
      PSI1   = P1
      PSI2   = CHI2                           ! ALL NH4HSO4 DELIQUESCED
C
C *** SOLVE EQUATIONS ; WITH ITERATIONS FOR ACTIVITY COEF. ************
C
      DO 10 I=1,NSWEEP
C
      A1 = XK11 *(WATER/GAMA(12))**2.0
      A3 = XK1  *WATER/GAMA(7)*(GAMA(8)/GAMA(7))**2.0
C
C  CALCULATE DISSOCIATION QUANTITIES
C
C *** FOR INTERACTION WITH HYDROPHILIC ORGANIC COMPOUNDS
      BB   = A3+LAMDA + ORGANION                ! KAPA
      CC   =-A3*(LAMDA + PSI1 + PSI2)
      DD   = BB*BB-4.D0*CC
      KAPA = 0.5D0*(-BB+SQRT(DD))
C
C *** CALCULATE SPECIATION ********************************************
C
C *** FOR INTERACTION WITH HYDROPHILIC ORGANIC COMPOUNDS
      MOLAL (1) = LAMDA + KAPA + ORGANION       ! HI
      MOLAL (2) = PSI1                          ! NAI
      MOLAL (3) = PSI2                          ! NH4I
      MOLAL (4) = ZERO                          ! CLI
      MOLAL (5) = KAPA                          ! SO4I
      MOLAL (6) = LAMDA + PSI1 + PSI2 - KAPA    ! HSO4I
      MOLAL (7) = ZERO                          ! NO3I
C
      CNAHSO4   = MAX(CHI1-PSI1,ZERO)
      CNH4HS4   = ZERO
C
      CALL SSH_CALCMR                               ! Water content
C
C *** CALCULATE ACTIVITIES OR TERMINATE INTERNAL LOOP *****************
C
      IF (FRST.AND.CALAOU .OR. .NOT.FRST.AND.CALAIN) THEN
         CALL SSH_CALCACT     
      ELSE
         GOTO 20
      ENDIF
10    CONTINUE
C
C *** CALCULATE OBJECTIVE FUNCTION ************************************
C
20    SSH_FUNCJ2 = MOLAL(2)*MOLAL(6)/A1 - ONE
C
C *** END OF FUNCTION SSH_FUNCJ2 *******************************************
C
      END

C=======================================================================
C
C *** ISORROPIA CODE
C *** SUBROUTINE SSH_CALCJ1
C *** CASE J1
C
C     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
C     1. SULFATE RICH, FREE ACID (SULRAT < 1.0)
C     2. THERE IS BOTH A LIQUID & SOLID PHASE
C     3. SOLIDS POSSIBLE : NH4HSO4, NAHSO4
C
C *** COPYRIGHT 1996-2006, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY,
C *** GEORGIA INSTITUTE OF TECHNOLOGY
C *** WRITTEN BY ATHANASIOS NENES
C *** UPDATED BY CHRISTOS FOUNTOUKIS
C
C=======================================================================
C
      SUBROUTINE SSH_CALCJ1
      INCLUDE 'isrpia.inc'
C
      DOUBLE PRECISION LAMDA, KAPA
      COMMON /SSHCASEJ/ CHI1, CHI2, CHI3, LAMDA, KAPA, PSI1, PSI2, PSI3, 
     &               A1,   A2,   A3


C
C *** SETUP PARAMETERS ************************************************
C
      CALAOU =.TRUE.               ! Outer loop activity calculation flag
      CHI1   = W(1)                ! Total NA initially as NaHSO4
      CHI2   = W(3)                ! Total NH4 initially as NH4HSO4
C
      PSI1LO = TINY                ! Low  limit
      PSI1HI = CHI1                ! High limit
C
C *** INITIAL VALUES FOR BISECTION ************************************
C
      X1 = PSI1HI
      Y1 = SSH_FUNCJ1 (X1)
      YHI= Y1                      ! Save Y-value at HI position
C
C *** YHI < 0.0 THE SOLUTION IS ALWAYS UNDERSATURATED WITH NH42SO4 ****
C
      IF (ABS(Y1).LE.EPS .OR. YHI.LT.ZERO) GOTO 50
C
C *** ROOT TRACKING ; FOR THE RANGE OF HI AND LO **********************
C
      DX = (PSI1HI-PSI1LO)/FLOAT(NDIV)
      DO 10 I=1,NDIV
         X2 = X1-DX
         Y2 = SSH_FUNCJ1 (X2)
         IF (SIGN(1.d0,Y1)*SIGN(1.d0,Y2).LT.ZERO) GOTO 20  ! (Y1*Y2.LT.ZERO)
         X1 = X2
         Y1 = Y2
10    CONTINUE
C
C *** { YLO, YHI } > 0.0 THE SOLUTION IS ALWAYS SUPERSATURATED WITH NH42SO4
C
      YLO= Y1                      ! Save Y-value at Hi position
      IF (YLO.GT.ZERO .AND. YHI.GT.ZERO) THEN
         Y3 = SSH_FUNCJ1 (ZERO)
         GOTO 50
      ELSE IF (ABS(Y2) .LT. EPS) THEN   ! X2 IS A SOLUTION 
         GOTO 50
      ELSE
         CALL SSH_PUSHERR (0001, 'CALCJ1')    ! WARNING ERROR: NO SOLUTION
         GOTO 50
      ENDIF
C
C *** PERFORM BISECTION ***********************************************
C
20    DO 30 I=1,MAXIT
         X3 = 0.5*(X1+X2)
         Y3 = SSH_FUNCJ1 (X3)
         IF (SIGN(1.d0,Y1)*SIGN(1.d0,Y3) .LE. ZERO) THEN  ! (Y1*Y3 .LE. ZERO)
            Y2    = Y3
            X2    = X3
         ELSE
            Y1    = Y3
            X1    = X3
         ENDIF
         IF (ABS(X2-X1) .LE. EPS*X1) GOTO 40
30    CONTINUE
      CALL SSH_PUSHERR (0002, 'CALCJ1')    ! WARNING ERROR: NO CONVERGENCE
C
C *** CONVERGED ; RETURN **********************************************
C
40    X3 = 0.5*(X1+X2)
      Y3 = SSH_FUNCJ1 (X3)
C 
50    RETURN
C
C *** END OF SUBROUTINE SSH_CALCJ1 ******************************************
C
      END




C=======================================================================
C
C *** ISORROPIA CODE
C *** SUBROUTINE SSH_FUNCJ1
C *** CASE J1
C
C     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
C     1. SULFATE RICH, FREE ACID (SULRAT < 1.0)
C     2. THERE IS BOTH A LIQUID & SOLID PHASE
C     3. SOLIDS POSSIBLE : (NH4)2SO4, NH4CL, NA2SO4
C
C *** COPYRIGHT 1996-2006, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY,
C *** GEORGIA INSTITUTE OF TECHNOLOGY
C *** WRITTEN BY ATHANASIOS NENES
C *** UPDATED BY CHRISTOS FOUNTOUKIS
C
C=======================================================================
C
      DOUBLE PRECISION FUNCTION SSH_FUNCJ1 (P1)
      INCLUDE 'isrpia.inc'
      DOUBLE PRECISION LAMDA, KAPA
      COMMON /SSHCASEJ/ CHI1, CHI2, CHI3, LAMDA, KAPA, PSI1, PSI2, PSI3, 
     &               A1,   A2,   A3


C
C *** SETUP PARAMETERS ************************************************
C
      FRST   = .TRUE.
      CALAIN = .TRUE.
C
      LAMDA  = MAX(W(2) - W(3) - W(1), TINY)  ! FREE H2SO4
      PSI1   = P1
C
C *** SOLVE EQUATIONS ; WITH ITERATIONS FOR ACTIVITY COEF. ************
C
      DO 10 I=1,NSWEEP
C
      A1 = XK11 *(WATER/GAMA(12))**2.0
      A2 = XK12 *(WATER/GAMA(09))**2.0
      A3 = XK1  *WATER/GAMA(7)*(GAMA(8)/GAMA(7))**2.0
C
      PSI2 = 0.5*(-(LAMDA+PSI1) + SQRT((LAMDA+PSI1)**2.D0+4.D0*A2))  ! PSI2
      PSI2 = MIN (PSI2, CHI2)
C
C *** FOR INTERACTION WITH HYDROPHILIC ORGANIC COMPOUNDS
      BB   = A3+LAMDA + ORGANION                ! KAPA
      CC   =-A3*(LAMDA + PSI2 + PSI1)
      DD   = BB*BB-4.D0*CC
      KAPA = 0.5D0*(-BB+SQRT(DD))    
C
C *** SAVE CONCENTRATIONS IN MOLAL ARRAY ******************************
C
C *** FOR INTERACTION WITH HYDROPHILIC ORGANIC COMPOUNDS      
      MOLAL (1) = LAMDA + KAPA + ORGANION       ! HI
      MOLAL (2) = PSI1                          ! NAI
      MOLAL (3) = PSI2                          ! NH4I
      MOLAL (4) = ZERO
      MOLAL (5) = KAPA                          ! SO4I
      MOLAL (6) = LAMDA + PSI1 + PSI2 - KAPA    ! HSO4I
      MOLAL (7) = ZERO
C
      CNAHSO4   = MAX(CHI1-PSI1,ZERO)
      CNH4HS4   = MAX(CHI2-PSI2,ZERO)
C
      CALL SSH_CALCMR                               ! Water content
C
C *** CALCULATE ACTIVITIES OR TERMINATE INTERNAL LOOP *****************
C
      IF (FRST.AND.CALAOU .OR. .NOT.FRST.AND.CALAIN) THEN
         CALL SSH_CALCACT     
      ELSE
         GOTO 20
      ENDIF
10    CONTINUE
C
C *** CALCULATE OBJECTIVE FUNCTION ************************************
C
20    SSH_FUNCJ1 = MOLAL(2)*MOLAL(6)/A1 - ONE
C
C *** END OF FUNCTION SSH_FUNCJ1 *******************************************
C
      END

