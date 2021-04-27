C!-----------------------------------------------------------------------
C!     Copyright (C) 2019 CEREA (ENPC) - INERIS
C!     SSH-aerosol is distributed under the GNU General Public License v3
C!-----------------------------------------------------------------------

      SUBROUTINE SSH_COMPUTE_DENSITY(nbin_aer, nesp_aer, EH2O, tinyc,
     &     Cesp, rho_esp, k, rho)

C------------------------------------------------------------------------
C     
C     -- DESCRIPTION 
C     
C     This is the main subroutine compute the density for a bin k
C     
C------------------------------------------------------------------------
C     
C     -- INPUT VARIABLES
C
C     nbin_aer : number of bin
C     nesp_aer : number of aerosol species
C     Cesp : concentration
C     rho_esp : list of density of each species
C     conctot : total concentration
C     k : bin where is the calcul of rho
C
C     -- INPUT/OUTPUT VARIABLES
C
C     rho : density of bin k
C
C------------------------------------------------------------------------
C
C     -- REMARKS
C
C     rho = sum_j(Q(k,j)) / sum_j(Q(k,j)/rho_esp(j))
C
C------------------------------------------------------------------------      

      

      IMPLICIT NONE
      !INCLUDE '../INC/dynaero.inc'
      INTEGER i,j,k
      
      INTEGER nbin_aer, nesp_aer, EH2O 
      DOUBLE PRECISION Cesp(nbin_aer,nesp_aer) 
      DOUBLE PRECISION rho_esp(nesp_aer)
      
      DOUBLE PRECISION rho
      DOUBLE PRECISION subrho, conctot
      DOUBLE PRECISION tinyc

      rho = 0.d0
      subrho = 0.d0
      conctot = 0.d0
      
      do j = 1, nesp_aer
        if (j .NE. EH2O) then
            subrho = subrho + Cesp(k,j)/rho_esp(j)
            conctot = conctot + Cesp(k,j)
       endif
      enddo

      if (subrho.LE.tinyc) then
         subrho = 0.d0
      endif

      if (conctot.GE.tinyc .AND. subrho.GT.0.d0) then
         rho = conctot/subrho
      else
         rho = 1.d-6 
      endif
      
      if (rho .EQ. 0.d0) then
         print *, "rho = ", rho
         print *, "subrho", subrho
         print *, "conctot", conctot
         stop
      endif
      END
