!C-----------------------------------------------------------------------
!C     Copyright (C) 2019 CEREA (ENPC) - INERIS
!C     SSH-aerosol is distributed under the GNU General Public License v3
!C-----------------------------------------------------------------------

!!$ ---------------------------------------------------------
!!$     Eps_machine	: Tolerance due to the lack of 
!!$	  		  precision of the machine
!!$
!!$	Eps_dbl_prec 	: Lower tolerance if the variation 
!!$	  		  of diameter after cond/evap is too small, 
!!$ 			  i.e. has no effect on the redistribution
!!$  
!!$     Eps_diam	: Tolerance for the diameter 
!!$  
!!$     PI		: PI number
!!$
!!$	diam_pass       : bin in which there is 100nm 
!!$			  (also defined in paraero.inc)
!!$
!!$     TINYM           : threshold mass conc (µg.m-3)     
!!$
!!$     TINYN           : threshold num conc (#part.m-3)   
!!$
!!$	DQLIMIT         : threshold (mass final - mass init) (µg.m-3)
!!$ ---------------------------------------------------------  

	DOUBLE PRECISION Eps_machine, Eps_dbl_prec, Eps_diam
	         PARAMETER(Eps_machine = 1.d-11)
	  	 PARAMETER(Eps_dbl_prec = 1.D-15)
	  	 PARAMETER(Eps_diam = 1.d-11)

 	DOUBLE PRECISION PI
      		 PARAMETER(PI = 3.14159265358979323846D0)

	DOUBLE PRECISION diam_pass   
	         PARAMETER(diam_pass = 1.D-01)

	DOUBLE PRECISION TINYM, TINYN
      	         PARAMETER(TINYM = 1.D-18) !1.D-19) !!1.D-15)
      		 PARAMETER(TINYN = 1.D-10) !1.D-06) !!1.D-03) 

	DOUBLE PRECISION DQLIMIT
!	       	PARAMETER(DQLIMIT = 0.D0)
