!!-----------------------------------------------------------------------
!!     Copyright (C) 2019 CEREA (ENPC) - INERIS
!!     SSH-aerosol is distributed under the GNU General Public License v3
!!-----------------------------------------------------------------------

SUBROUTINE SSH_EULER_MASS(ns, nesp, EH2O,dbound, grand, X, diameter_before_redist, fixed_diameter, log_fixed_diameter, &
     kloc, LMD, Qesp, N)

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
!!$     diameter_before_redist : list of mean diameter after 
!!$                              condensation/evaporation
!!$     kloc           : list of bin where is diameter_before_redist
!!$     log_fixed_diameter : log(fixed_diameter)
!!$     X              : log(diameter_before_redist)
!!$     fixed_diameter : list of mean diameter before 
!!$                      condensation/evaporation
!!$     j              : time integration
!!$     section_pass   : bin include 100nm  
!!$     LMD            : list of liquid mass density of each species
!!$ 
!!$     -- VARIABLES
!!$     
!!$     Q              : Mass concentration
!!$     rho            : density per bin
!!$     Eps_machine    : tolerance due to the lack of 
!!$                      precision of the machine     
!!$     Qdonne_esp     : Temporary Mass concentration 
!!$     frac           : fraction define by X and log_fixed_diameter
!!$     Qd             : fraction of mass concentration give at 
!!$                      the adjacent bin
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


  ! ------ Input 
  INTEGER, INTENT(in)::ns, nesp
  DOUBLE PRECISION, DIMENSION(ns), INTENT(in) :: X, log_fixed_diameter
  DOUBLE PRECISION, DIMENSION(ns), INTENT(in) :: fixed_diameter , diameter_before_redist
  DOUBLE PRECISION, DIMENSION(ns+1), INTENT(in) ::dbound
  INTEGER, DIMENSION(ns), INTENT(in) :: grand
  INTEGER, DIMENSION(ns), INTENT(in) :: kloc
  DOUBLE PRECISION, DIMENSION(nesp), INTENT(in) :: LMD
  integer eh2o
  ! ------ Input/Output
  DOUBLE PRECISION, DIMENSION(ns, nesp), INTENT(inout) :: Qesp
  DOUBLE PRECISION, DIMENSION(ns), INTENT(inout) :: N
  ! ------ 
  INTEGER k, iesp
  DOUBLE PRECISION, DIMENSION(ns) :: rho
  DOUBLE PRECISION, DIMENSION(ns,nesp) :: Qdonne_esp
  DOUBLE PRECISION, DIMENSION(ns) :: Q
  DOUBLE PRECISION Qd, frac

  DO k = 1, ns
     DO iesp=1,nesp
        Qdonne_esp(k,iesp) = 0D0
     ENDDO
  ENDDO


  DO k = 1, ns
     Q(k) =  0d0
     DO iesp = 1, nesp
        if (iesp .ne. eh2o) then
        Q(k) = Q(k) + Qesp(k,iesp)
        endif    
     ENDDO
  ENDDO

  DO k = 1,ns 

     IF (grand(k) == 0) THEN

        IF (kloc(k) .NE. 1) THEN
           frac = (log_fixed_diameter(k) - X(k))/&
                (log_fixed_diameter(k) - DLOG10(fixed_diameter(kloc(k)-1))) 
        ELSE
           frac = (log_fixed_diameter(k) - X(k))/&
                (log_fixed_diameter(k) - DLOG10(dbound(1)))/2.d0
        ENDIF
        IF (frac .GT. 1) THEN
           PRINT * , "In SIREAM/euler_mass.f90: frac > 1."
           STOP
        ENDIF
        DO iesp = 1, nesp 
           Qd = Qesp(k,iesp) * frac   
           IF (kloc(k).NE.1) THEN 
              Qdonne_esp(kloc(k)-1, iesp) = &
                   Qdonne_esp(kloc(k)-1, iesp) + Qd
              Qdonne_esp(kloc(k), iesp) = &
                   Qdonne_esp(kloc(k), iesp) + Qesp(k, iesp) - Qd 
           ELSE              
              Qdonne_esp(kloc(k), iesp) = &
                   Qdonne_esp(kloc(k), iesp) + Qesp(k, iesp)
           ENDIF
           Qesp(k, iesp) = 0.d0 
        Enddo

     ELSE 
        IF (kloc(k) .NE. ns) THEN
           frac = (X(k) - log_fixed_diameter(k))/&
                (DLOG10(fixed_diameter(kloc(k)+1)) - log_fixed_diameter(k)) 
        ELSE
           frac = (X(k) - log_fixed_diameter(k))/&
                (DLOG10(dbound(ns+1)) - log_fixed_diameter(k))/2.d0
        ENDIF
        IF (frac .GT. 1) THEN
           PRINT * , "In SIREAM/euler_mass.f90: frac > 1."
           STOP
        ENDIF

        DO iesp = 1, nesp  
           Qd = Qesp(k, iesp) * frac
           IF (kloc(k).NE.ns) THEN 
              Qdonne_esp(kloc(k)+1, iesp) = &
                   Qdonne_esp(kloc(k)+1, iesp) + Qd
              Qdonne_esp(kloc(k), iesp) = &
                   Qdonne_esp(kloc(k), iesp) + Qesp(k, iesp) - Qd
           ELSE              
              Qdonne_esp(kloc(k), iesp) = &
                   Qdonne_esp(kloc(k), iesp) + Qesp(k, iesp)
           ENDIF
           Qesp(k, iesp) = 0.d0           
        ENDDO
     ENDIF
  ENDDO

  !***** Recalculation of number concentration from mass concentration
  Q = 0.D0
  N = 0.d0
  DO k = 1,ns
     DO iesp = 1,nesp
        Qesp(k, iesp) = Qesp(k, iesp) + Qdonne_esp(k, iesp) 
        if (iesp .ne. eh2o) then
            Q(k) = Q(k) + Qesp(k, iesp)
        endif
     ENDDO
     
     CALL SSH_COMPUTE_DENSITY(ns,nesp,EH2O,TINYM,Qesp,LMD,k,rho(k))

     N(k) = (Q(k) * 6.D0)/(PI * rho(k) * &
          (fixed_diameter(k)*fixed_diameter(k)*fixed_diameter(k)))
     
  ENDDO  
 
  CALL SSH_TEST_MASS_NB(ns,nesp,rho,dbound,Q,N,Qesp)

END SUBROUTINE SSH_EULER_MASS



