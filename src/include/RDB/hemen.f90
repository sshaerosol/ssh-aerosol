!!-----------------------------------------------------------------------
!!     Copyright (C) 2019 CEREA (ENPC) - INERIS
!!     SSH-aerosol is distributed under the GNU General Public License v3
!!-----------------------------------------------------------------------

SUBROUTINE SSH_HEMEN(ns, nesp,eh2o, grand, section_pass, d, dbound, &
     X,diam , logdiam,kloc, alpha, LMD, rho, Qesp, N)

!!$------------------------------------------------------------------------
!!$     
!!$     -- DESCRIPTION 
!!$     
!!$     This subroutine redistribute the concentrations after the GDE.
!!$     Hybrid Euler Mass Euler Number.     
!!$     
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
!!$     logdiam        : log(d)
!!$     X              : log(diam)
!!$     d              : list of mean diameter before condensation/evaporation
!!$     j              : time integration
!!$     section_pass   : bin include 100nm  
!!$     LMD            : list of liquid mass density of each species
!!$ 
!!$     -- VARIABLES
!!$     
!!$     Q              : Mass concentration
!!$     N_esp          : Number concentration by bin and species
!!$     rho            : density per bin
!!$     Eps_machine    : tolerance due to the lack of precision of the machine     
!!$     Qdonne_esp     : Temporary Mass concentration 
!!$     Ndonne_esp     : Temporary number concentration 
!!$     frac           : fraction define by X and logdiam
!!$     Qd             : fraction of mass concentration give at the adjacent bin
!!$     Nd             : fraction of number concentration give at the adjacent bin
!!$
!!$     -- INPUT/OUTPUT VARIABLES
!!$      
!!$     N              : Number concentration by bin
!!$     Qesp           : Mass concentration by bin and species
!!$      
!!$     -- OUTPUT VARIABLES
!!$     
!!$     
!!$------------------------------------------------------------------------



  IMPLICIT NONE

  INCLUDE '../INC/parameuler.inc'


!!! ------ Input 
  INTEGER, INTENT(in) :: ns, nesp, section_pass
  INTEGER, DIMENSION(ns), INTENT(in) :: grand
  DOUBLE PRECISION, DIMENSION(ns), INTENT(in) :: X, logdiam
  DOUBLE PRECISION, DIMENSION(ns), INTENT(in) :: d, diam
  DOUBLE PRECISION, DIMENSION(ns+1), INTENT(in) :: dbound
  INTEGER, DIMENSION(ns), INTENT(in) :: kloc
  DOUBLE PRECISION, DIMENSION(ns, nesp), INTENT(in) :: alpha
  DOUBLE PRECISION, DIMENSION(nesp), INTENT(in) :: LMD
  integer eh2o

!!! ------ Input/Output 
  DOUBLE PRECISION, DIMENSION(ns), INTENT(inout) :: N
  DOUBLE PRECISION, DIMENSION(ns, nesp), INTENT(inout) :: Qesp

!!! ------  
  INTEGER k, jesp, ik
  DOUBLE PRECISION, DIMENSION(ns) :: rho
  DOUBLE PRECISION, DIMENSION(ns, nesp) :: N_esp, Qd_tmp, Nd_tmp
  DOUBLE PRECISION, DIMENSION(ns) :: Q
  DOUBLE PRECISION frac, Nd, Qd
  DOUBLE PRECISION, DIMENSION(ns, nesp) :: Qdonne_esp, Ndonne_esp

!!! ~~~~~~ Distribution of number concentration by species
!!! ~~~~~ alpha is a fraction of each species except water
  Q = 0
  DO k = 1,ns
     DO jesp = 1, nesp
        N_esp(k, jesp) = alpha(k,jesp) * N(k)
        if (jesp .ne. eh2o) then
        Q(k) = Q(k) + Qesp(k, jesp)
    endif   
     ENDDO
     CALL SSH_COMPUTE_DENSITY(ns, nesp,eh2o, TINYN, N_esp, LMD, k, rho(k))
  ENDDO

!!! ~~~~~~ initialization transfert vector
  Ndonne_esp = 0.d0
  Qdonne_esp = 0.d0

  DO k = 1, ns 

!!! ~~~~~~~~~~~~~~~~ redistribute according to Euler number algorithm ~~~~~~~~~~~~~~~~~~~
     IF(kloc(k) .LE. section_pass-1) THEN
        !! ~~~~~~ Redistribute between bins kloc(k)-1 and kloc(k) ~~~~~~~~~~~~~~~~~~~~~~~
        IF (grand(k) == 0) THEN 

           IF (kloc(k) .NE. 1) THEN
              frac = (logdiam(k) - X(k)) &
                   /(logdiam(k) - DLOG10(d(kloc(k)-1)))

           ELSE
              frac = (logdiam(k) - X(k)) &
                   /(logdiam(k) - DLOG10(dbound(1)))/2.d0
           ENDIF
           IF (frac .GT. 1) THEN
              PRINT * , "In SIREAM/hemen.f90: frac > 1."
              STOP
           ENDIF

           DO jesp = 1, nesp
              Nd = N_esp(k, jesp) * frac

              IF (kloc(k) .NE. 1) THEN 

                 Ndonne_esp(kloc(k)-1, jesp)  = &
                      Ndonne_esp(kloc(k)-1, jesp) + Nd 
                 Ndonne_esp(kloc(k), jesp) = &
                      Ndonne_esp(kloc(k), jesp) + N_esp(k, jesp) - Nd
              ELSE
                 Ndonne_esp(kloc(k), jesp) = &
                      Ndonne_esp(kloc(k), jesp) + N_esp(k, jesp) 
              ENDIF
              N_esp(k,jesp) = 0.d0
           ENDDO


           !! ~~~~~~ Redistribute between bins kloc(k)+1 and kloc(k) ~~~~~~~~~~~~~~~~~~~~
        ELSE  ! ~~~~~ grand(k) == 1
           IF(kloc(k) .NE. ns) THEN
              frac =  (X(k) - logdiam(k))/&
                   (DLOG10(d(kloc(k)+1))-logdiam(k))

           ELSE
              frac =  (X(k) - logdiam(k))/&
                   (DLOG10(dbound(ns+1)) - logdiam(k))/2.d0

           ENDIF
           IF (frac .GT. 1) THEN
              PRINT * , "In SIREAM/hemen.f90: frac > 1."
              STOP
           ENDIF

           DO jesp = 1, nesp
              Nd =  N_esp(k, jesp) * frac

              IF (kloc(k).NE.ns) THEN
                 Ndonne_esp(kloc(k), jesp) = &
                      Ndonne_esp(kloc(k), jesp) + N_esp(k, jesp) - Nd
                 Ndonne_esp(kloc(k)+1, jesp) = &
                      Ndonne_esp(kloc(k)+1, jesp) + Nd
              ELSE
                 Ndonne_esp(kloc(k), jesp) = &
                      Ndonne_esp(kloc(k), jesp) + N_esp(k, jesp)
              ENDIF


              N_esp(k,jesp) = 0.d0
           ENDDO
        ENDIF

!!! ~~~~~~~~~~~~~~~~ redistribute according to Euler mass algorithm ~~~~~~~~~~~~~~~~~~~
     ELSEIF (kloc(k) .GE. section_pass+1) THEN

        !! ~~~~~~ Redistribute between bins kloc(k)-1 and kloc(k) ~~~~~~~~~~~~~~~~~~~~~~~
        IF (grand(k) == 0) THEN
           frac = (logdiam(k) - X(k))/ &
                (logdiam(k)-DLOG10(d(kloc(k)-1)))
           IF (frac .GT. 1) THEN
              PRINT * , "In SIREAM/hemen.f90: frac > 1."
              STOP
           ENDIF


           DO jesp = 1, nesp
              Qd = Qesp(k,jesp) * frac


              IF (kloc(k).NE.1) THEN
                 Qdonne_esp(kloc(k)-1, jesp) = &
                      Qdonne_esp(kloc(k)-1, jesp) + Qd
                 Qdonne_esp(kloc(k), jesp) = &
                      Qdonne_esp(kloc(k), jesp) + Qesp(k, jesp) - Qd
              ELSE
                 Qdonne_esp(kloc(k), jesp) = &
                      Qdonne_esp(kloc(k), jesp) + Qesp(k, jesp)
              ENDIF

              Qesp(k, jesp) = 0.d0 
           ENDDO
           !! ~~~~~~ Redistribute between bins kloc(k)+1 and kloc(k) ~~~~~~~~~~~~~~~~~~~~~~~~~~
        ELSE  ! ~~~~~ grand(k) == 1

           IF (kloc(k).NE.ns) THEN
              frac =  (X(k) - logdiam(k))/&
                   (DLOG10(d(kloc(k)+1))-logdiam(k))

           ELSE
              frac =  (X(k) - logdiam(k))/&
                   (DLOG10(dbound(ns+1)) - logdiam(k))/2.d0

           ENDIF
           IF (frac .GT. 1) THEN
              PRINT * , "In SIREAM/hemen.f90: frac > 1."
              STOP
           ENDIF

           DO jesp = 1, nesp  
              Qd = Qesp(k, jesp) * frac
              IF (kloc(k).NE.ns) THEN 
                 Qdonne_esp(kloc(k)+1, jesp) = &
                      Qdonne_esp(kloc(k)+1, jesp) + Qd
                 Qdonne_esp(kloc(k), jesp) = &
                      Qdonne_esp(kloc(k), jesp) + Qesp(k, jesp) - Qd
              ELSE
                 Qdonne_esp(kloc(k), jesp) = &
                      Qdonne_esp(kloc(k), jesp) + Qesp(k, jesp) 

              ENDIF
              Qesp(k, jesp) = 0.d0           
           ENDDO

        ENDIF

     ELSE !! kloc(k) in section_pass: 
          !! choose between Euler number and Euler mass to redistribute the section k 
          !! in section_pass and redistribute mass in the section just above 
          !! or number in the section just below



        IF (grand(k) == 0) THEN 
           frac = (logdiam(k) - X(k))/ &
                (logdiam(k)-DLOG10(d(kloc(k)-1)))
           IF (frac .GT. 1) THEN
              PRINT * , "In SIREAM/hemen.f90: frac > 1."
              STOP
           ENDIF

           DO jesp = 1, nesp
              Nd = N_esp(k, jesp) * frac
              Ndonne_esp(kloc(k)-1, jesp) = &
                   Ndonne_esp(kloc(k)-1, jesp) + Nd

              IF (d(section_pass).LE.diam_pass) THEN ! ~~~~ Euler number
                 Ndonne_esp(kloc(k), jesp) = &
                      Ndonne_esp(kloc(k), jesp) + N_esp(k, jesp) - Nd
              ELSE
!!! rho(k) is the density after condensation/evaporation and before redistribution
                 Qdonne_esp(kloc(k), jesp) = Qdonne_esp(kloc(k), jesp) &
                      + (PI/6D0) * rho(k) * diam(k) * diam(k) * diam(k) &
                      * (N_esp(k,jesp)-Nd)

              ENDIF

              N_esp(k, jesp) = 0d0
              Qesp(k, jesp) = 0d0
           ENDDO

        ELSE
           frac = (X(k) - logdiam(k))/ &
                (DLOG10(d(kloc(k)+1)) - logdiam(k))
           DO jesp = 1, nesp
              Qd = Qesp(k,jesp) * frac
              Qdonne_esp(kloc(k)+1, jesp) = &
                   Qdonne_esp(kloc(k)+1, jesp) + Qd

              IF (d(section_pass).LE.diam_pass) THEN ! ~~~~ Euler number
!!! rho(k) is the density after condensation/evaporation and before redistribution
                 Ndonne_esp(kloc(k), jesp) = Ndonne_esp(kloc(k), jesp) &
                      + (6.D0/PI) / rho(k) / diam(k) / diam(k) / diam(k) &
                      * (Qesp(k,jesp)-Qd)

              ELSE
                 Qdonne_esp(kloc(k), jesp) = &
                      Qdonne_esp(kloc(k), jesp) + Qesp(k, jesp) - Qd
              ENDIF

              Qesp(k, jesp) = 0d0
              N_esp(k, jesp) = 0d0
           ENDDO

        ENDIF
     ENDIF
  ENDDO



  !! ~~~~~~~~~~~~~~~~ REDIAGNOTIC ~~~~~~~~~~~~~~~~~~~~~~~~
  ! recompute total mass for each species and total number
  N = 0.d0
  Q = 0.d0
  DO k = 1,section_pass-1
     DO jesp = 1, nesp
        N_esp(k, jesp) = N_esp(k, jesp) + Ndonne_esp(k, jesp)
        if (jesp .ne. eh2o) then
        N(k) = N(k) + N_esp(k, jesp)
        endif
     ENDDO
  ENDDO
  DO k = section_pass + 1, ns
     DO jesp = 1, nesp
        Qesp(k, jesp) = Qesp(k, jesp) + Qdonne_esp(k, jesp)
        if (jesp .ne. eh2o) then
        Q(k) = Q(k) + Qesp(k,jesp)
        endif
     ENDDO
  ENDDO

  k = section_pass
  IF (d(section_pass).LE.diam_pass) then 
     ! ~~~~ Euler number
     DO jesp = 1, nesp
        N_esp(k, jesp) = N_esp(k, jesp) + Ndonne_esp(k, jesp)
        if (jesp .ne. eh2o) then
        N(k) = N(k) + N_esp(k, jesp)
        endif
     ENDDO
  ELSE 
     ! ~~~~ Euler mass
     DO jesp = 1, nesp
        Qesp(k, jesp) = Qesp(k, jesp) + Qdonne_esp(k, jesp)
        if (jesp .ne. eh2o) then
        Q(k) = Q(k) + Qesp(k,jesp)
        endif
     ENDDO
  ENDIF

  DO k = 1,section_pass-1
     CALL SSH_COMPUTE_DENSITY(ns, nesp,eh2o,TINYN, N_esp, LMD, k, rho(k))
     DO jesp = 1, nesp
        Qesp(k, jesp) = (PI/6D0) * rho(k) * d(k) * d(k) * d(k) * N_esp(k, jesp)
        if (jesp .ne. eh2o) then
        Q(k) = Q(k) + Qesp(k,jesp)
        endif
     ENDDO
  ENDDO

  DO k = section_pass+1 , ns
     CALL SSH_COMPUTE_DENSITY(ns, nesp,eh2o,TINYM, Qesp, LMD, k, rho(k))
     N(k) = Q(k) * 6.D0 / (PI * rho(k) * d(k) * d(k) * d(k))
  ENDDO

  k = section_pass
  IF (d(section_pass).LE.diam_pass) then 
     ! ~~~~ Euler number
     CALL SSH_COMPUTE_DENSITY(ns, nesp, eh2o,TINYN,N_esp, LMD, k, rho(k))
     DO jesp = 1, nesp
        Qesp(k, jesp) = (PI/6D0) * rho(k) * d(k) * d(k) * d(k) * N_esp(k, jesp)
        if (jesp .ne. eh2o) then
        Q(k) = Q(k) + Qesp(k,jesp)
        endif
     ENDDO
  ELSE 
     ! ~~~~ Euler mass
     CALL SSH_COMPUTE_DENSITY(ns, nesp,eh2o,TINYM, Qesp, LMD, k, rho(k))
     N(k) = Q(k) * 6.D0 / (PI * rho(k) * d(k) * d(k) * d(k))
  ENDIF

  CALL SSH_TEST_MASS_NB(ns,nesp,rho,dbound,Q,N,Qesp)

END SUBROUTINE SSH_HEMEN
