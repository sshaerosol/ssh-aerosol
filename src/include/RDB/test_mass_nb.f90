!!-----------------------------------------------------------------------
!!     Copyright (C) 2019 CEREA (ENPC) - INERIS
!!     SSH-aerosol is distributed under the GNU General Public License v3
!!-----------------------------------------------------------------------

SUBROUTINE SSH_TEST_MASS_NB(ns,nesp,rho,dbound,Q,N,Qesp)

  IMPLICIT NONE

  INCLUDE '../INC/parameuler.inc'

  INTEGER k, jesp

  INTEGER, INTENT(in) :: ns, nesp
  DOUBLE PRECISION, DIMENSION(ns), INTENT(in) :: rho
  DOUBLE PRECISION, DIMENSION(ns), INTENT(in) :: dbound
  DOUBLE PRECISION, DIMENSION(ns), INTENT(inout) :: Q, N 
!  DOUBLE PRECISION :: Qmin
  DOUBLE PRECISION, DIMENSION(ns, nesp), INTENT(inout) :: Qesp
  
 
  DO k = 1, ns

!     Qmin = rho(k) * (PI/6.D0) * (dbound(1) ** 3D0)
     IF (N(k) .LT. 0D0)THEN
        PRINT*, k,"stop, N negatif after redistribution"
        PRINT*, "N",N
        STOP
     ENDIF

     IF (Q(k) .LT. 0D0)THEN
        PRINT*, k,"stop, Q negatif after redistribution"
        PRINT*, "Q",Q
        STOP
     ENDIF
     
! !    IF(Q(k) .LT. Qmin * N(k) .OR. N(k) .LT. 1D0 .OR. Q(k) .LT. nesp*TINYM) THEN
!     IF (N(k).LE.TINYN) THEN
!        !Q(k) = 0.D0
!        N(k) = TINYN
!        DO jesp = 1, nesp
!           Qesp(k, jesp) = TINYM
!        ENDDO
!     ENDIF
  ENDDO

END SUBROUTINE SSH_TEST_MASS_NB
