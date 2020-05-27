!!-----------------------------------------------------------------------
!!     Copyright (C) 2019 CEREA (ENPC) - INERIS
!!     SSH-aerosol is distributed under the GNU General Public License v3
!!-----------------------------------------------------------------------

SUBROUTINE SSH_EULER_COUPLED(ns, nesp, eh2o,dbound, grand, d, diam, kloc, alpha, LMD, Q_esp, N)

!!$------------------------------------------------------------------------
!!$     
!!$     -- INPUT VARIABLES
!!$     
!!$     
!!$     ns             : number of sections
!!$     nesp           : number of species
!!$     dbound         : list of limit bound diameter [\mu m]
!!$     grand          : list of 0 or 1
!!$                      1 = cutting with the upper box
!!$                      0 = cutting with the lower box
!!$     alpha          : list of fraction of each species in Q
!!$     diam           : list of mean diameter after condensation/evaporation
!!$     kloc           : list of bin where is diam
!!$     d              : list of mean diameter before condensation/evaporation
!!$     j              : time integration
!!$     LMD            : list of liquid mass density of each species
!!$ 
!!$     -- VARIABLES
!!$     
!!$     Q              : Mass concentration
!!$     rho            : density per bin
!!$     Eps_machine    : tolerance due to the lack of precision of the machine     
!!$     Q_1_esp, N_1_esp : Concentration stay in the current bin
!!$     Q_2evap_esp, N_2evap_esp : Concentration give at the lower bin
!!$     Q_2cond_esp, N_2cond_esp : Concentration give at the upper bin
!!$
!!$     -- INPUT/OUTPUT VARIABLES
!!$      
!!$     N            : Number concentration by bin
!!$     Qesp         : Mass concentration by bin and species
!!$      
!!$     -- OUTPUT VARIABLES
!!$     
!!$     
!!$------------------------------------------------------------------------


  IMPLICIT NONE

  INCLUDE "../INC/parameuler.inc"

  ! ------ Input
  INTEGER, INTENT(in) :: ns,  nesp 
  DOUBLE PRECISION, DIMENSION(ns+1),  INTENT(in) :: dbound
  INTEGER, DIMENSION(ns), INTENT(in) :: grand, kloc
  DOUBLE PRECISION, DIMENSION(ns), INTENT(in) :: diam, d
  DOUBLE PRECISION, DIMENSION(ns,nesp), INTENT(in) :: alpha 
  DOUBLE PRECISION, DIMENSION(ns), INTENT(in) :: LMD
  integer eh2o

  ! ------ Input/Output
  DOUBLE PRECISION, DIMENSION(ns), INTENT(inout) :: N
  DOUBLE PRECISION, DIMENSION(ns, nesp), INTENT(inout)  ::Q_esp 

  ! ------ 
  INTEGER k, jesp
  DOUBLE PRECISION, DIMENSION(ns) :: rho
  DOUBLE PRECISION, DIMENSION(ns) :: Q

  DOUBLE PRECISION, DIMENSION(ns) :: Q_1, Q_2, N_1, N_2
  DOUBLE PRECISION :: AQ, BQ, CQ, DQ, AN, BN, CN, DN
  DOUBLE PRECISION :: N_ancien, N_nouveau, Q_ancien, Q_nouveau

  DOUBLE PRECISION, DIMENSION(ns, nesp) :: N_esp
  DOUBLE PRECISION, DIMENSION(ns, nesp) :: Q_1_esp, N_1_esp
  DOUBLE PRECISION, DIMENSION(ns, nesp) :: Q_2evap_esp, N_2evap_esp 
  DOUBLE PRECISION, DIMENSION(ns, nesp) :: Q_2cond_esp, N_2cond_esp
  DOUBLE PRECISION, DIMENSION(ns) :: RQ, RN

  !! ~~~~~~~~~~~ INITIALISATION

  Q_1_esp =0D0
  N_1_esp =0D0

  Q_2evap_esp = 0D0
  N_2evap_esp = 0D0

  Q_2cond_esp = 0D0
  N_2cond_esp = 0D0

  DO k = 1,ns
     Q(k) = 0D0
     DO jesp = 1, nesp
        N_esp(k, jesp) = alpha(k,jesp) * N(k)
        Q(k) = Q(k) + Q_esp(k, jesp)
     ENDDO
  ENDDO

  RQ = 0.d0
  RN = 0.d0

 !! ~~~~~~~~~~~~~ tests
  N_ancien = 0D0
  N_nouveau = 0D0
  Q_ancien = 0D0
  Q_nouveau = 0D0

  DO k = 1, ns
     Q_ancien = Q_ancien + Q(k)
     N_ancien = N_ancien + N(k)
  ENDDO

  !! ~~~~~~~~~~~~ ALGORITHM

  DO k = 1,ns
     IF (grand(k) == 0) THEN

        IF (kloc(k) .GT. 1) THEN

           AQ = 1D0 - ((d(kloc(k)-1)/diam(k))*(d(kloc(k)-1)/diam(k)) &
                *(d(kloc(k)-1)/diam(k)))
           BQ = 1D0 - ((d(kloc(k)-1)/d(kloc(k)))*(d(kloc(k)-1)/d(kloc(k))) &
                *(d(kloc(k)-1)/d(kloc(k))))
           CQ = 1D0 - ((d(kloc(k))/diam(k))*(d(kloc(k))/diam(k)) &
                *(d(kloc(k))/diam(k)))
           DQ = 1D0 - ((d(kloc(k))/d(kloc(k)-1))*(d(kloc(k))/d(kloc(k)-1)) &
                *(d(kloc(k))/d(kloc(k)-1)))

           AN = (diam(k)*diam(k)*diam(k)) - &
                (d(kloc(k)-1)*d(kloc(k)-1)*d(kloc(k)-1))
           BN = (d(kloc(k))*d(kloc(k))*d(kloc(k))) - &
                (d(kloc(k)-1)*d(kloc(k)-1)*d(kloc(k)-1))
           CN = (diam(k)*diam(k)*diam(k)) - &
                (d(kloc(k))*d(kloc(k))*d(kloc(k)))
           DN = (d(kloc(k)-1)*d(kloc(k)-1)*d(kloc(k)-1)) - &
                (d(kloc(k))*d(kloc(k))*d(kloc(k)))


           DO jesp = 1, nesp

              Q_1_esp(kloc(k), jesp) = Q_esp(k, jesp) * AQ/BQ &
                   + Q_1_esp(kloc(k), jesp)
              Q_2evap_esp(kloc(k), jesp) = Q_esp(k, jesp) * CQ/DQ &
                   + Q_2evap_esp(kloc(k), jesp)
              N_1_esp(kloc(k), jesp) = N_esp(k, jesp) * AN/BN &
                   + N_1_esp(kloc(k), jesp)
              N_2evap_esp(kloc(k), jesp) = N_esp(k, jesp) * CN/DN &
                   + N_2evap_esp(kloc(k), jesp)
           ENDDO

           !!~~~~~~~~~~~~~~~~~~~~ kloc(k) == 1
        ELSE

           AQ = 1D0 -((dbound(1)/diam(k))*(dbound(1)/diam(k)) &
                *(dbound(1)/diam(k)))
           BQ = 1D0 - ((dbound(1)/d(kloc(k)))*(dbound(1)/d(kloc(k))) &
                *(dbound(1)/d(kloc(k))))   

           AN = (diam(k)*diam(k)*diam(k)) &
                -(dbound(1)*dbound(1)*dbound(1))
           BN = (d(kloc(k))*d(kloc(k))*d(kloc(k))) &
                - (dbound(1)*dbound(1)*dbound(1))

           DO jesp = 1, nesp
              !Q_1_esp(kloc(k), jesp) = Q_esp(k, jesp) * AQ/BQ + Q_1_esp(kloc(k), jesp)
              !N_1_esp(kloc(k), jesp) = N_esp(k, jesp) * AN/BN + N_1_esp(kloc(k), jesp)
              Q_1_esp(kloc(k), jesp) = Q_esp(k, jesp) &
                   + Q_1_esp(kloc(k), jesp)
              N_1_esp(kloc(k), jesp) = N_esp(k, jesp) &
                   + N_1_esp(kloc(k), jesp)

           ENDDO

        ENDIF
        !!~~~~~~~~~~ grand(k) == 1

     ELSE

        IF ( kloc(k) .LT. ns) THEN

           AQ = 1D0 - ((d(kloc(k)+1)/diam(k))*(d(kloc(k)+1)/diam(k))* &
                (d(kloc(k)+1)/diam(k)))
           BQ = 1D0 - ((d(kloc(k)+1)/d(kloc(k)))*(d(kloc(k)+1)/d(kloc(k)))* &
                (d(kloc(k)+1)/d(kloc(k))))
           CQ = 1D0 -((d(kloc(k))/diam(k))*(d(kloc(k))/diam(k))* &
                (d(kloc(k))/diam(k)))
           DQ = 1D0 -((d(kloc(k))/d(kloc(k)+1))*(d(kloc(k))/d(kloc(k)+1)) &
                *(d(kloc(k))/d(kloc(k)+1)))

           AN = diam(k)*diam(k)*diam(k) &
                - (d(kloc(k)+1)*d(kloc(k)+1)*d(kloc(k)+1))
           BN = d(kloc(k))*d(kloc(k))*d(kloc(k)) &
                - (d(kloc(k)+1)*d(kloc(k)+1)*d(kloc(k)+1))
           CN = diam(k)*diam(k)*diam(k) &
                - (d(kloc(k))*d(kloc(k))*d(kloc(k)))
           DN = d(kloc(k)+1)*d(kloc(k)+1)*d(kloc(k)+1) &
                - (d(kloc(k))*d(kloc(k))*d(kloc(k)))

           DO jesp = 1, nesp
              Q_1_esp(kloc(k), jesp) = Q_esp(k, jesp) * AQ/BQ &
                   + Q_1_esp(kloc(k), jesp)
              Q_2cond_esp(kloc(k), jesp) = Q_esp(k, jesp) * CQ/DQ &
                   + Q_2cond_esp(kloc(k), jesp)
              N_1_esp(kloc(k), jesp) = N_esp(k, jesp) * AN/BN &
                   + N_1_esp(kloc(k), jesp)
              N_2cond_esp(kloc(k), jesp) = N_esp(k, jesp) * CN/DN &
                   + N_2cond_esp(kloc(k), jesp)

           ENDDO

           !!~~~~~~~~~~~~~~~~~~~~ kloc(k) == ns
        ELSE

           AQ = 1D0 - ((dbound(ns+1)/diam(k)) &
                *(dbound(ns+1)/diam(k))*(dbound(ns+1)/diam(k)))
           BQ = 1D0 - ((dbound(ns+1)/d(kloc(k))) &
                *(dbound(ns+1)/d(kloc(k)))*(dbound(ns+1)/d(kloc(k))))

           AN = diam(k)*diam(k)*diam(k) &
                - (dbound(ns+1)*dbound(ns+1)*dbound(ns+1))
           BN = d(kloc(k))*d(kloc(k))*d(kloc(k)) &
                - (dbound(ns+1)*dbound(ns+1)*dbound(ns+1))

           DO jesp = 1, nesp
              !Q_1_esp(kloc(k), jesp) = Q_esp(k, jesp) * AQ/BQ + Q_1_esp(kloc(k), jesp)
              !N_1_esp(kloc(k), jesp) = N_esp(k, jesp) * AN/BN + N_1_esp(kloc(k), jesp)
              Q_1_esp(kloc(k), jesp) = Q_esp(k, jesp) &
                   + Q_1_esp(kloc(k), jesp)
              N_1_esp(kloc(k), jesp) = N_esp(k, jesp) &
                   + N_1_esp(kloc(k), jesp)

           ENDDO

        ENDIF

     ENDIF

     !! ~~~~~~~~~~~~~~~~~~~~ TESTS
     IF(k .NE. 1 .AND. k .NE. ns)THEN

        DO jesp = 1, nesp
           RQ(k) = RQ(k) + Q_1_esp(k, jesp) &
                + Q_2evap_esp(k, jesp) + Q_2cond_esp(k, jesp)
           RN(k) = RN(k) + N_1_esp(k, jesp) &
                + N_2evap_esp(k, jesp) + N_2cond_esp(k, jesp)
        ENDDO

!!$        IF ( DABS(1D0 - (RQ(k)/Q(k))) .GE. Eps_machine &
!!$             .AND. DABS(1D0 - (RQ(k)/Q(k))).NE. 1d0 &
!!$             .AND. Q(k) .NE. 0D0) THEN
!!$
!!$           PRINT*,k,"non conservation de la masse ds algo!!"
!!$           PRINT*,"1 - RQ(k)/Qgros(k) = ",1D0 - RQ(k)/Q(k)
!!$           PRINT*, "RQ",RQ(k), "Q",Q(k)
!!$
!!$        ENDIF
!!$
!!$        IF ( DABS(1D0 - RN(k)/N(k)) .GE. Eps_machine &
!!$             .AND. DABS(1D0 - (RN(k)/N(k))).NE. 1d0&
!!$             .AND. N(k) .NE. 0D0) THEN
!!$
!!$           PRINT*,k,"non conservation du nombre ds algo!!"
!!$           PRINT*,"1 - RN(k)/N(k)",1D0 - RN(k)/N(k)
!!$
!!$        ENDIF
        DO jesp = 1, nesp

           IF ( N_1_esp(k, jesp) .LT. 0D0 ) THEN
              PRINT*, jesp,k,"->",kloc(k), grand(k),"N_1_esp(k,jesp) negatif"
              PRINT*, "N_1_esp(k,jesp)",N_1_esp(k,jesp),"AN/BN", AN/BN, BN, 'N', N_esp(k,jesp)
              PRINT*, "d(kloc(k)+1)",d(kloc(k)+1),"d(kloc(k))",d(kloc(k)), "diam(k)", diam(k)
              PRINT*, 'dbis',diam
              PRINT*, "kloc", kloc
!!              STOP
              N_1_esp(k, jesp) = 0.D0
              Q_1_esp(k, jesp) = 0.D0
           ENDIF

           IF ( Q_1_esp(k,jesp) .LT. 0D0 ) THEN
              PRINT*, jesp, k,"->",kloc(k), grand(k),"Q_1_esp(k,jesp) negatif"
              PRINT*, "Q_1_esp(k,jesp)",Q_1_esp(k, jesp), "AQ/BQ", AQ/BQ, BQ, 'Q', Q_esp(k,jesp)
              PRINT*, "d(kloc(k)+1)",d(kloc(k)+1),"d(kloc(k))",d(kloc(k)), "diam(k)", diam(k)
              PRINT*, "kloc", kloc
!!              STOP
              N_1_esp(k, jesp) = 0.D0
              Q_1_esp(k, jesp) = 0.D0
           ENDIF
           IF (kloc(k) .LT. ns .AND. kloc(k) .GT. 1) THEN
              IF ( N_2evap_esp(kloc(k),jesp) .LT. 0D0 ) THEN
                 PRINT*, k,"->",kloc(k),grand(k),"N_2evap(k) negatif"
                 PRINT*, "N_2evap_esp(k)",N_2evap_esp(kloc(k),jesp)
                 PRINT*, "CN",CN
                 PRINT*, "diam(k)",diam(k)
                 PRINT*, "d(kloc(k))",d(kloc(k))
                 PRINT*, "DN",DN, "d(k+1)",d(k+1)
!!              STOP
                 N_2evap_esp(k, jesp) = 0.D0
                 Q_2evap_esp(k, jesp) = 0.D0
              ENDIF
              

              IF ( Q_2evap_esp(kloc(k),jesp) .LT. 0D0 ) THEN
                 PRINT*, k,jesp,"Q_2evap_esp(k,jesp) negatif"
                 PRINT*, "Q_2evap_esp(k,jesp)",Q_2evap_esp(k,jesp)
                 print*, "Q_esp(k,jesp)",Q_esp(k, jesp), "CQ/DQ", CQ/DQ, CQ, DQ
!!              STOP
                 N_2evap_esp(k, jesp) = 0.D0
                 Q_2evap_esp(k, jesp) = 0.D0
              ENDIF

              IF ( N_2cond_esp(kloc(k),jesp) .LT. 0D0 ) THEN
                 PRINT*, k,"->",kloc(k),grand(k),"N_2cond(k) negatif"
                 PRINT*, "N_2cond_esp(k)",N_2cond_esp(kloc(k),jesp)
                 PRINT*, "CN",CN
                 PRINT*, "diam(k)",diam(k)
                 PRINT*, "d(kloc(k))",d(kloc(k))
                 PRINT*, "DN",DN, "d(k+1)",d(k+1)
 !!              STOP
                 N_2cond_esp(k, jesp) = 0.D0
                 Q_2cond_esp(k, jesp) = 0.D0
              ENDIF
              

              IF ( Q_2cond_esp(kloc(k),jesp) .LT. 0D0 ) THEN
                 PRINT*, k,jesp,"Q_2cond_esp(k,jesp) negatif"
                 PRINT*, "Q_2cond_esp(k,jesp)",Q_2cond_esp(k,jesp)
                 print*, "Q_esp(k,jesp)",Q_esp(k, jesp), "CQ/DQ", CQ/DQ, CQ, DQ
!!              STOP
                 N_2cond_esp(k, jesp) = 0.D0
                 Q_2cond_esp(k, jesp) = 0.D0
              ENDIF

           ENDIF
        ENDDO

     ELSE 
        DO jesp =1,nesp
           IF ( N_1_esp(k,jesp) .LT. 0D0 ) THEN
              PRINT*,  k,"N_1(k) negatif here"
              PRINT*, "N_1(k)",N_1_esp(k,jesp)
!!              STOP
                 N_1_esp(k, jesp) = 0.D0
                 Q_1_esp(k, jesp) = 0.D0
           ENDIF

           IF ( Q_1_esp(k,jesp) .LT. 0D0 ) THEN
              PRINT*, k,"Q_1(k) negatif here"
              PRINT*, "Q_1(k)",Q_1_esp(k,jesp)
 !!              STOP
                 N_1_esp(k, jesp) = 0.D0
                 Q_1_esp(k, jesp) = 0.D0
           ENDIF
        ENDDO

     ENDIF
     !! ~~~~~~~~~~~~~ Fin tests
  ENDDO
  !! ~~~~~~~~~~~~~~ REPARTITION
  DO k = 1, ns
     DO jesp = 1, nesp
        Q_esp(k, jesp) = Q_1_esp(k, jesp)
        N_esp(k, jesp) = N_1_esp(k, jesp)
     ENDDO
  ENDDO


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  N = 0.d0
  Q = 0.d0
  N_2 = 0.d0

  DO k = 1, ns
     DO jesp = 1, nesp
        N(k) = N(k) + N_esp(k, jesp)
        Q(k) = Q(k) + Q_esp(k, jesp)
        N_2(k) = N_2(k) + N_2evap_esp(k,jesp)+ N_2cond_esp(k,jesp)

     ENDDO
  ENDDO

  DO k = 1, ns
     N_nouveau = N_nouveau + N(k)
     Q_nouveau = Q_nouveau + Q(k)
  ENDDO


  !test de conservation du nombre et de la masse, se verifie si on ferme les bornes

!!$  IF ( DABS (1D0 - (N_nouveau/N_ancien)) .GE. Eps_machine &
!!$       .AND. N_ancien .GT. 0.d0) THEN            
!!$     PRINT*,"non conservation du nombre total !!"
!!$     PRINT*,"1 - N_nouveau/N_ancien ", 1D0 - N_nouveau/N_ancien   
!!$     PRINT*, "N_nouveau", N_nouveau 
!!$     PRINT*, "N_ancien", N_ancien
!!$     !STOP
!!$  ENDIF
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  DO k = 1, ns
     IF (grand(k) == 0) THEN
        IF (kloc(k) .NE. 1) THEN
           N = 0.d0
           Q = 0.d0
           DO jesp = 1, nesp
              Q_esp(kloc(k)-1, jesp) = Q_esp(kloc(k)-1, jesp) &
                   + Q_2evap_esp(kloc(k), jesp) 
              N_esp(kloc(k)-1, jesp) = N_esp(kloc(k)-1, jesp) &
                   + N_2evap_esp(kloc(k), jesp)             
           ENDDO
        ENDIF
     ELSE
        IF (kloc(k) .NE. ns) THEN
           DO jesp = 1, nesp
              Q_esp(kloc(k)+1, jesp) = Q_esp(kloc(k)+1, jesp) &
                   + Q_2cond_esp(kloc(k), jesp) 
              N_esp(kloc(k)+1, jesp) = N_esp(kloc(k)+1, jesp) &
                   + N_2cond_esp(kloc(k), jesp)
           ENDDO
        ENDIF
     ENDIF
  ENDDO

  N = 0.d0
  Q = 0.d0
  DO k = 1, ns
     DO jesp = 1, nesp
        N(k) = N(k) + N_esp(k, jesp)
        Q(k) = Q(k) + Q_esp(k, jesp)
     ENDDO
  ENDDO


  N_nouveau = 0.d0
  Q_nouveau = 0.d0
  DO k = 1, ns
     N_nouveau = N_nouveau + N(k)
     Q_nouveau = Q_nouveau + Q(k)
  ENDDO

END SUBROUTINE SSH_EULER_COUPLED


