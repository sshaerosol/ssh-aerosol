!!-----------------------------------------------------------------------
!!     Copyright (C) 2019 CEREA (ENPC) - INERIS
!!     SSH-aerosol is distributed under the GNU General Public License v3
!!-----------------------------------------------------------------------

SUBROUTINE SSH_EULER_NUMBER(ns, nesp, eh2o, dbound, grand, alpha, &
     fixed_diameter, diameter_before_redist, X, log_fixed_diameter, kloc, LMD, rho, Qesp, N)

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
!!$     diameter_before_redist           : list of mean diameter after condensation/evaporation
!!$     kloc           : list of bin where is diameter_before_redist
!!$     log_fixed_diameter        : log(fixed_diameter)
!!$     X              : log(diameter_before_redist)
!!$     fixed_diameter              : list of mean diameter before condensation/evaporation
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
!!$     Ndonne_esp     : Temporary number concentration 
!!$     frac           : fraction define by X and log_fixed_diameter
!!$     Nd             : fraction of number concentration give at the adjacent bin
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
  INCLUDE '../INC/parameuler.inc'

  ! ------ Input 
  INTEGER, INTENT(in) :: ns, nesp
  INTEGER, DIMENSION(ns), INTENT(in) :: grand
  DOUBLE PRECISION, DIMENSION(ns), INTENT(in) :: X, log_fixed_diameter 
  DOUBLE PRECISION, DIMENSION(ns), INTENT(in) :: fixed_diameter , diameter_before_redist
  DOUBLE PRECISION, DIMENSION(ns+1), INTENT(in) ::dbound
  INTEGER, DIMENSION(ns), INTENT(in) :: kloc
  DOUBLE PRECISION, DIMENSION(ns, nesp), INTENT(in) :: alpha 
  DOUBLE PRECISION, DIMENSION(nesp), INTENT(in) :: LMD
 integer eh2o

  ! ------ Input/Output
  DOUBLE PRECISION, DIMENSION(ns), INTENT(inout) :: N
  DOUBLE PRECISION, DIMENSION(ns, nesp), INTENT(inout) :: Qesp

  ! ------
  INTEGER k, jesp
  DOUBLE PRECISION, DIMENSION(ns) :: rho
  DOUBLE PRECISION, DIMENSION(ns) :: Q, QT
  DOUBLE PRECISION, DIMENSION(ns, nesp) :: N_esp, Ndonne_esp
  DOUBLE PRECISION Nd, frac


  Q = 0.d0
  DO k = 1,ns
     DO jesp = 1, nesp
        if (jesp .ne. eh2o) then
        Q(k) = Q(k) + Qesp(k, jesp) 
        endif 
        N_esp(k, jesp) = alpha(k,jesp) * N(k)
     ENDDO
  ENDDO

  !***** Calcul total mass per species

  DO k = 1,ns
     DO jesp = 1, nesp
        Ndonne_esp(k, jesp) = 0d0
     ENDDO
  ENDDO

  DO k = 1,ns 

     IF (grand(k) == 0)THEN

        IF (kloc(k) .NE. 1) THEN
           frac = (log_fixed_diameter(k) - X(k))/ &
                (log_fixed_diameter(k) - DLOG10(fixed_diameter(kloc(k)-1))) 
        ELSE
           frac = (log_fixed_diameter(k) - X(k))/ &
                (log_fixed_diameter(k) - DLOG10(dbound(1)))/2.d0
        ENDIF
           IF (frac .GT. 1) THEN
           PRINT * , "In SIREAM/euler_number.f90: frac > 1."
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
           N_esp(k, jesp) = 0.d0
        ENDDO

     ELSE
        IF (kloc(k) .NE. ns) THEN
           frac = (X(k) - log_fixed_diameter(k))/ &
                (DLOG10(fixed_diameter(kloc(k)+1)) - log_fixed_diameter(k)) 
        ELSE
           frac = (X(k) - log_fixed_diameter(k))/ &
                (DLOG10(dbound(ns+1)) - log_fixed_diameter(k))/2.d0
        ENDIF
           IF (frac .GT. 1) THEN
           PRINT * , "In SIREAM/euler_number.f90: frac > 1."
              STOP
           ENDIF

        DO jesp = 1, nesp

           Nd =  N_esp(k, jesp) * frac

           IF (kloc(k) .NE. ns) THEN
              Ndonne_esp(kloc(k)+1,jesp) = &
                   Ndonne_esp(kloc(k)+1, jesp) + Nd
              Ndonne_esp(kloc(k), jesp) = &
                   Ndonne_esp(kloc(k), jesp) + N_esp(k, jesp) - Nd
           ELSE
              Ndonne_esp(kloc(k), jesp) = &
                   Ndonne_esp(kloc(k),jesp) + N_esp(k, jesp)
           ENDIF
           N_esp(k, jesp)=0.d0
        ENDDO
     ENDIF

  ENDDO



  N = 0.d0
  Q = 0.d0
  DO k = 1,ns
     DO jesp = 1, nesp
        N_esp(k, jesp) = N_esp(k, jesp) + Ndonne_esp(k, jesp)
        if (jesp .ne. eh2o) then
        N(k) = N(k) + N_esp(k, jesp)
        endif
     ENDDO

     !***** Recalculation of mass concentration from number concentration 
     IF (N(k) .LT. TINYN) THEN
        rho(k) = 1.d0
     ELSE
        CALL SSH_COMPUTE_DENSITY(ns,nesp,eh2o, TINYN,N_esp,LMD,k,rho(k))
     ENDIF
     DO jesp = 1, nesp
        Qesp(k, jesp) = rho(k) * (PI/6D0) * N_esp(k,jesp) &
             * (fixed_diameter(k)*fixed_diameter(k)*fixed_diameter(k))
        if (jesp .ne. eh2o) then
        Q(k) = Q(k) + Qesp(k,jesp)
        endif
     ENDDO
  ENDDO



  CALL SSH_TEST_MASS_NB(ns,nesp,rho,dbound,Q,N,Qesp)

END SUBROUTINE SSH_EULER_NUMBER


