C-----------------------------------------------------------------------
C     Copyright (C) 2003-2007, ENPC - INRIA - EDF R&D
C     Author(s): Edouard Debry
C     
C     This file is part of the Size Resolved Aerosol Model (SIREAM), a
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

      SUBROUTINE LOCATE(n,x,xx,j)

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

