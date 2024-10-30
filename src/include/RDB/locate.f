C!!-----------------------------------------------------------------------
C!!     Copyright (C) 2019-2024 CEREA (ENPC) - INERIS
C!!     SSH-aerosol is distributed under the GNU General Public License v3
C!!-----------------------------------------------------------------------
      
      SUBROUTINE SSH_LOCATE(n,x,xx,j)

C------------------------------------------------------------------------
C     
C     -- DESCRIPTION 
C     
C     This subroutine locates the bin j where xx is located in the 
C     discretization x(*).     
C     
C------------------------------------------------------------------------
C     
C     -- INPUT VARIABLES
C     
C     N : length of array x(*).
C     X : array of discretization where xx has to be located.
C     XX: point to be located.
C     
C     -- INPUT/OUTPUT VARIABLES
C     
C     
C     -- OUTPUT VARIABLES
C     
C     J : interval/bin where xx is located.
C     
C------------------------------------------------------------------------
C     
C     -- REMARKS
C     
C------------------------------------------------------------------------
C     
C     -- MODIFICATIONS
C     
C     2005/3/23: cleaning (Bruno Sportisse, CEREA).
C     
C------------------------------------------------------------------------
C     
C     -- AUTHOR(S)
C     
C     2004: Edouard Debry, CEREA.
C     
C------------------------------------------------------------------------

      IMPLICIT NONE

      INTEGER j,n
      DOUBLE PRECISION xx,x(n)

      INTEGER jl,jm,ju

      IF (xx.LT.x(1)) THEN
         j=0
      ELSEIF (xx.GT.x(n)) THEN
         j=n
      ELSE
         jl=1
         ju=n

         DO WHILE (ju-jl.GT.1)
            jm=(ju+jl)/2
            
            IF (xx.GT.x(jm)) THEN
               jl=jm
            ELSE 
               ju=jm
            ENDIF
         END DO
         
         j=jl
      ENDIF

      END

      
      SUBROUTINE SSH_LOCATE2(n,dbound,DSF,xx,j)

C------------------------------------------------------------------------
C     
C     -- DESCRIPTION 
C     
C     This subroutine locates the bin j where xx is located in the 
C     discretization x(*).     
C     
C------------------------------------------------------------------------
C     
C     -- INPUT VARIABLES
C     
C     N : length of array x(*).
C     X : array of discretization where xx has to be located.
C     XX: point to be located.
C     
C     -- INPUT/OUTPUT VARIABLES
C     
C     
C     -- OUTPUT VARIABLES
C     
C     J : interval/bin where xx is located.
C     
C------------------------------------------------------------------------
C     
C     -- REMARKS
C     
C------------------------------------------------------------------------
C     

      IMPLICIT NONE

      INTEGER j,n,i
      DOUBLE PRECISION xx,DSF(n),dbound(n+1)

      INTEGER jl,jm,ju

      j = 0
      IF (xx.LT.dbound(1)) THEN
         j=0
      ELSEIF (xx.GT.dbound(n+1)) THEN
         j=n+1
      ELSE
         do i=1,n
            if((xx.GE.dbound(i)).AND.(xx.LT.dbound(i+1))) then
               j = i
            endif
         enddo
      ENDIF

      If(xx.LE.DSF(j)) j = j-1

      END

