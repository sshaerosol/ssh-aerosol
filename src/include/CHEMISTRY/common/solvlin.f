
      SUBROUTINE solvlin (NS,Kindlu,DLa,DLalu,DLx,DLb)

C------------------------------------------------------------------------
C
C     -- DESCRIPTION
C
C     This routine solves DLA * DLX = DLB where DLA is an input matrix,
C     and DLB is an input vector.
C
C------------------------------------------------------------------------
C
C     -- INPUT VARIABLES
C
C     KINDLU: 0 if DLALU is not a LU factorization of DLA. If KINDLU is
C     # not zero, DLALU is assumed to be a LU factorization of DLA.
C     DLA: matrix (NESP x NESP).
C     DLB: right-hand-side vector (NEPS) of the equation to be solved.
C
C     -- INPUT/OUTPUT VARIABLES
C
C     DLALU: if KINDLU is not zero, DLALU is an LU factorisation of DLA.
C     # Otherwise, on exit, DLALU is an LU factorization of DLA.
C     IPVT: pivot indices; for 1 <= i <= NESP, row i of the
C     # matrix was interchanged with row IPVT(i).
C
C     -- OUTPUT VARIABLES
C
C     DLX: solution of DLA * DLX = DLB.
C
C------------------------------------------------------------------------
C
C     -- REMARKS
C
C------------------------------------------------------------------------
C
C     -- MODIFICATIONS
C
C     2004/01/28: NGAS and NSRCGAS rather than NESP and NSRC
C     for the split of chemistry and aerosol (Hadjira Foudhil).
C------------------------------------------------------------------------
C
C     -- AUTHOR(S)
C
C     Denis Quélo, CEREA, June 2001.
C
C------------------------------------------------------------------------

      IMPLICIT NONE

C     -- INCLUDE FILES
C     PARADOM: parameters for all forced fields.
C     PARACHEM: parameters for sizes of 'chemical' arrays.

      INTEGER NS


      INTEGER Kindlu
      INTEGER Ji, Jj
      DOUBLE PRECISION DLa(NS,NS)
      DOUBLE PRECISION DLalu(NS,NS)
      DOUBLE PRECISION DLx(NS), DLb(NS)

      DO Ji=1,NS
         DLx(Ji)=DLb(Ji)
      ENDDO

C------------------------------------------------------------------------
C     1 - Solve DLa * Dlx = Dlb

      IF (Kindlu .EQ. 0) THEN   ! DLalu is not
                                ! an LU factorization of DLa.
         DO Jj=1,NS
            DO Ji=1,NS
               DLalu(Ji,Jj)=DLa(Ji,Jj)
            ENDDO
         ENDDO

         CALL LU_decompose(NS,DLalu)
         CALL LU_solve(NS,DLalu,DLx)

      ELSE                      ! DLalu is an LU factorization of DLa.

         CALL LU_solve(NS,DLalu,DLx)

      ENDIF

      END
