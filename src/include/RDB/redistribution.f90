!!-----------------------------------------------------------------------
!!     Copyright (C) 2019 CEREA (ENPC) - INERIS
!!     SSH-aerosol is distributed under the GNU General Public License v3
!!-----------------------------------------------------------------------

SUBROUTINE REDISTRIBUTION(ns, naer, EH2O, dbound, fixed_diameter, scheme, &
     section_pass, LMD, &
     DQLIMIT, Qesp, N, Q, with_fixed_density,fixed_density,d_after_cond) 

!!$------------------------------------------------------------------------
!!$     
!!$     -- INPUT VARIABLES
!!$     
!!$     
!!$     ns             : number of sections
!!$     naer           : number of species
!!$     dbound         : list of limit bound diameter [\mu m]
!!$     fixed_diameter : list of mean diameter 
!!$     scheme         : redistribution scheme
!!$                      3 = euler_mass
!!$                      4 = euler_number
!!$                      5 = hemen 
!!$                      6 = euler_coupled
!!$     section_pass   : bin include 100nm  
!!$     LMD            : list of liquid mass density of each species
!!$ 
!!$     -- VARIABLES
!!$     
!!$     rho            : density per bin
!!$     
!!$     Eps_diam       : tolerance to the diameter which may be slightly above 
!!$                      the edge of the section
!!$     Eps_dbl_prec   : tolerance lower in the case where the diameter so little 
!!$                      increases or decreases as one considers that it is 
!!$                      not happening 
!!$                    : it can not work with a time not too restrictive 
!!$     grand          : list of 0 or 1
!!$                      1 = cutting with the upper box
!!$                      0 = cutting with the lower box
!!$     d_after_cond   : list of mean diameter after condensation/evaporation
!!$     kloc           : list of bin where is diam
!!$     logdiam        : log(d)
!!$     X              : log(diam)
!!$     alpha          : list of fraction of each species in Q
!!$
!!$     -- INPUT/OUTPUT VARIABLES
!!$      
!!$     N            : Number concentration per bin
!!$     Qesp         : Mass concentration per bin and species
!!$     Q            : Total mass concentration per bin
!!$      
!!$     -- OUTPUT VARIABLES
!!$     
!!$     
!!$------------------------------------------------------------------------

  IMPLICIT NONE
  INCLUDE '../INC/parameuler.inc'

  ! ------ Input
  INTEGER, INTENT(in) :: ns, naer, section_pass, scheme
  DOUBLE PRECISION, DIMENSION(ns+1), INTENT(in) :: dbound
  DOUBLE PRECISION, DIMENSION(ns), INTENT(in) :: fixed_diameter
  DOUBLE PRECISION, DIMENSION(naer), INTENT(in) :: LMD
  INTEGER, INTENT(in) :: EH2O

  ! ------ Input/Output
  DOUBLE PRECISION, DIMENSION(ns), INTENT(inout) :: N, Q
  DOUBLE PRECISION, DIMENSION(ns, naer), INTENT(inout) :: Qesp

  ! ------ 
  INTEGER k, js, jaer, loc
  DOUBLE PRECISION, DIMENSION(ns) :: rho
  INTEGER, DIMENSION(ns) :: grand
  INTEGER, DIMENSION(ns) :: kloc
  DOUBLE PRECISION, DIMENSION(ns) :: dbis_tmp
  DOUBLE PRECISION :: dbis
  DOUBLE PRECISION , DIMENSION(ns) :: X, d_after_cond, logfixed_diameter
  DOUBLE PRECISION, DIMENSION(ns,naer) :: alpha 
  
  DOUBLE PRECISION :: frac
  integer with_fixed_density ! YK
  double precision :: fixed_density

  !! Calcul of d(k) = dsqrt(dbound(k)*dbound(k+1))
  !! Reestimate d(k) from number and mass concentration BEFORE condensation/evaporation
  if(scheme.EQ.11) then
    call REDIST(ns,naer,N,Qesp,fixed_density,dbound,d_after_cond)
  else 
  if(scheme.EQ.12) then
    DO k = 1, ns
       if (with_fixed_density .eq. 0) then
          CALL compute_density(ns,naer, EH2O,TINYM,Qesp,LMD,k,rho(k))
       else
          rho(k) = fixed_density  !! YK
       endif
    ENDDO
    call REDIST_EULERCOUPLE(ns,naer,N,Qesp,fixed_density,rho,dbound,d_after_cond)
    !call REDIST_MOVINGDIAM(ns,naer,N,Qesp,fixed_density,rho,dbound,d_after_cond)
  else
  if(scheme.EQ.13) then
    DO k = 1, ns
       if (with_fixed_density .eq. 0) then
          CALL compute_density(ns,naer, EH2O,TINYM,Qesp,LMD,k,rho(k))
       else
          rho(k) = fixed_density  !! YK
       endif
    ENDDO
    call REDIST_MOVINGDIAM(ns,naer,N,Qesp,fixed_density,rho,dbound,d_after_cond)
  else
  !***** Calcul dry total mass per bin
  Q=0.D0
  DO js = 1, ns
     DO jaer = 1, naer-1 !KS Water is last
        if (jaer .NE. EH2O) then
        Q(js)= Q(js) + Qesp(js, jaer) 
        endif
     ENDDO
  ENDDO


  !***** Calcul fraction of each composant of mass
  DO js = 1, ns
     kloc(js) = js
     DO jaer = 1, naer
        IF (Q(js) .GT. TINYM) THEN
           alpha(js,jaer) = Qesp(js, jaer)/Q(js)
        ELSE
           alpha(js, jaer) = 0.D0
        ENDIF
     ENDDO
  ENDDO

  !****** Calcul of new mean diameter after c/e
  DO k = 1, ns
     if (isNaN(N(k))) then
        print *, k, "N(k) NaN dans redist"
     endif
     
     if (with_fixed_density .eq. 0) then
        CALL compute_density(ns,naer, EH2O,TINYM,Qesp,LMD,k,rho(k))
     else
	rho(k) = fixed_density  !! YK
     endif
     IF ( N(k) .GT. TINYN .and. Q(k) .GT. TINYM) THEN
        dbis = ((Q(k) * 6D0)/(PI * rho(k) * N(k))) **(1D0/3D0)        
        !call locate(ns,dbound,dbis,kloc(k))
        IF (dbis .LT. dbound(1)) THEN
           Q(k) = 0.D0
           N(k) = 0.D0
           dbis = dbound(1) !fixed_diameter(1)
           kloc(k) = 1 
        ELSEIF (dbis .GT. dbound(ns+1)) THEN
           dbis = dbound(ns+1)
           kloc(k)= ns 
        ELSE
           loc = 1
           DO WHILE(dbis .GT. dbound(loc+1))
              loc=loc+1
           ENDDO
           kloc(k) = loc
        ENDIF
     ELSE
        dbis = fixed_diameter(k)
        kloc(k)=k
     ENDIF
     dbis_tmp(kloc(k)) = dbis
  Enddo

  DO k = 1, ns
     d_after_cond(k) = dbis_tmp(kloc(k))
     X(k) = DLOG10(dbis_tmp(kloc(k)))
     logfixed_diameter(k) = DLOG10(fixed_diameter(kloc(k)))
     IF(dbis_tmp(kloc(k)) .GT. fixed_diameter(kloc(k))) THEN
        grand(k) = 1
     ELSE
        grand(k) = 0 
     ENDIF


  ENDDO
 
  !****** Select the redistribution method
  SELECT CASE (scheme)
  CASE (3)
     CALL EULER_MASS(ns, naer, eh2o,dbound, grand, &
          X, d_after_cond, fixed_diameter, logfixed_diameter, kloc, LMD, Qesp, N)

  CASE (4)
     CALL EULER_NUMBER(ns, naer,eh2o, dbound, grand, alpha, &
          fixed_diameter, d_after_cond, X, logfixed_diameter, kloc, LMD, rho, Qesp, N)
 
  CASE (5)
     CALL HEMEN(ns, naer, eh2o,grand, section_pass, fixed_diameter, &
          dbound, X, d_after_cond, logfixed_diameter, kloc, alpha, LMD, &
          rho, Qesp, N)

  CASE (6)
     CALL EULER_COUPLED(ns, naer,eh2o,dbound, grand, fixed_diameter, d_after_cond, &
           kloc, alpha, LMD, Qesp, N)

  CASE (7)
     CALL EULER_NUMBER_NORM(ns, naer, eh2o,dbound, grand, alpha, &
          fixed_diameter, d_after_cond, X, logfixed_diameter, kloc, LMD,DQLIMIT, rho, Qesp, N)
 
  CASE (8)
     CALL HEMEN_NORM(ns, naer,  EH2O, grand, section_pass, &
          dbound, X, d_after_cond, fixed_diameter, kloc, alpha, LMD, &
          DQLIMIT, rho, Qesp, N)

  CASE (9)
     CALL EULER_COUPLED_NORM(ns, naer, eh2o,dbound, grand, fixed_diameter, d_after_cond, &
           kloc, alpha, LMD, DQLIMIT, Qesp, N)

  CASE (10)
     CALL MOVING_DIAM(ns, naer, eh2o, dbound,kloc,Qesp, N)  
  
  CASE DEFAULT
     PRINT*, "Please choose from the following redistribution methods : ", &
          "number-conserving, interpolation, euler-mass, euler-number, ",&
          "hemen, euler-coupled."
  END SELECT

 endif
 endif
 endif
! Total dry mass
!  Q=0.D0
!  DO js = 1, ns
!     DO jaer = 1, naer-1
!      !  if (jaer .NE. EH2O) then
!        Q(js)= Q(js) + Qesp(js, jaer) 
!       ! endif
!     ENDDO
!  ENDDO

END SUBROUTINE REDISTRIBUTION

!**************************************************
!*     *
!*     ROUTINE REDISTRIBUTING AEROSOL NUMBER AND  *
!*     MASS CONCENTRATIONS ON THE FIXED GRID SIZE *
!*     *
!**************************************************
      SUBROUTINE REDIST(ns,naer,n,q,fixed_density,dbound,diam)
!**************************************************
!*     *
!*     ns   (IN)    running number of eqs.       *
!*     q(*)  (INOUT) gas/aero conc (µg.m-3)       *
!*     *
!**************************************************
      IMPLICIT NONE

      include '../Module/CONST_B.inc'
      INTEGER ns,naer,NB
      DOUBLE PRECISION q(ns,naer),n(ns),qt(ns)
      DOUBLE PRECISION Nnew(ns)
      INTEGER ji,jj,jesp,j1lo,j1hi,js1,js2,icpt
      DOUBLE PRECISION qnew(ns,naer),qtnew(ns)
      DOUBLE PRECISION x2lo,x2hi,xmin,xmax
      DOUBLE PRECISION xx,tmp1,tmp2,frac
      DOUBLE PRECISION CST,fixed_density
      DOUBLE PRECISION dbound(ns+1),XSD(ns),XBF(ns+1),XBD(ns+1)
      DOUBLE PRECISION HSD(ns),XSF(ns),diam(ns)

      CST  =0.523598775598D0
      NB = NS+1

!******zero init
      DO js1=1,NS
         qt(js1) = 0.d0
         Do jj=1,naer-1
             qt(js1) = qt(js1) + q(js1,jj)
             qnew(js1,jj)=0.D0
         Enddo
         qnew(js1,naer)=0.D0
         qtnew(js1)=0.D0
         Nnew(js1) = 0.D0
      END DO

      DO js1=1,NS+1
         XBF(js1)=DLOG(fixed_density*CST*dbound(js1)**3)
      ENDDO
      DO js1=1,NS
         XSF(js1) = (XBF(js1)+XBF(js1+1))*5.D-01
         if(N(js1).GT.TINYN) then
            XSD(js1) = DLOG(qt(js1)/N(js1))
         else
            XSD(js1) = XSF(js1)
         endif
      ENDDO

      call SIZEBND(ns,XSD,XSF,XBF,HSD,XBD)

!******compute redistribution
      x2lo=XBD(1)
      
      CALL LOCATE(NB,XBF,x2lo,j1lo)
      
      DO js2=1,NS
         x2hi=XBD(js2+1)
         CALL LOCATE(NB,XBF,x2hi,j1hi)
         IF (N(js2).GT.TINYN) THEN
                                ! redistribute over fixed sections
            DO js1=MAX0(j1lo,1),MIN0(j1hi,NS)
               
               xmin=DMAX1(x2lo,XBF(js1))
               xmax=DMIN1(x2hi,XBF(js1+1))

               frac=(xmax-xmin)/HSD(js2)

               Nnew(js1)=Nnew(js1)+N(js2)*frac

               DO jesp=1,naer
                  qnew(js1,jesp)=qnew(js1,jesp)+q(js2,jesp)*frac
               END DO
            END DO
         ENDIF
         
         x2lo=x2hi
         j1lo=j1hi
      END DO
     
!!$!******check if quantities are not too small
!!$      icpt=1                    ! to enter in the loop
!!$
!!$      DO WHILE (icpt.GT.0)
!!$         icpt = 0
!!$         
!!$         DO js1=1,NS
!!$            DO jesp=1,naer-1
!!$               qtnew(js1)=qtnew(js1)+qnew(js1,jesp)
!!$            END DO
!!$            
!!$            tmp1= Nnew(js1)*( Nnew(js1)-TINYN)
!!$            tmp2=qtnew(js1)*(qtnew(js1)-TINYM)
!!$            
!!$            IF ( tmp1.LT.0.D0.OR.tmp2.LT.0.D0 ) THEN
!!$               if (Nnew(js1).eq.0.) then
!!$                  write(*,*)'pb1 in redist',Nnew(js1)
!!$               endif
!!$               if ((qtnew(js1)/Nnew(js1)).le.0.) then
!!$                  write(*,*)'pb2 in redist',qtnew(js1)/Nnew(js1)
!!$               endif
!!$               xx=DLOG(qtnew(js1)/Nnew(js1))
!!$               frac=(xx-XBD(js1))/HSD(js1)
!!$
!!$               js2=js1+1
!!$               IF (frac.LT.5.D-01) js2=js1-1
!!$               
!!$               IF (js1.EQ.1)  js2=js1+1
!!$               IF (js1.EQ.NS) js2=js1-1
!!$               
!!$               Nnew(js2)=Nnew(js2)+Nnew(js1)
!!$               qtnew(js2)=qtnew(js2)+qtnew(js1)
!!$               Nnew(js1)=0.d0
!!$               qtnew(js1)=0.d0
!!$
!!$               DO jesp=1,naer
!!$                  qnew(js2,jesp)=qnew(js2,jesp)+qnew(js1,jesp)
!!$                  qnew(js1,jesp)=0.d0
!!$               END DO
!!$
!!$               icpt=icpt+1
!!$
!!$            ENDIF
!!$         END DO
!!$      END DO
      
!******turn back to conc vector
      
      DO js1=1,NS
         Do jj=1,naer
             q(js1,jj)=qnew(js1,jj)
         Enddo
         N(js1) = Nnew(js1)
      END DO

!**************************************************
      END SUBROUTINE REDIST
!**************************************************

!**************************************************
!*                                                *
!*     ROUTINE COMPUTING SIZE BOUNDS VARIABLES    *
!*                                                *
!**************************************************
      SUBROUTINE SIZEBND(ns2,XSD,XSF,XBF,HSD,XBD)
!**************************************************
      IMPLICIT NONE

!******
      INTEGER ns2,jb,js,itest,ixbd(NS2+1),NB2
      DOUBLE PRECISION dx(NS2),d2x(NS2)
      DOUBLE PRECISION dxb,dvar(2)
      DOUBLE PRECISION XSD(ns2),XSF(ns2),XBF(ns2+1)
      DOUBLE PRECISION HSD(ns2),XBD(ns2+1)

      NB2 = NS2 + 1
!****** zero init
      DO js=1,NS2
         dx(js) =0.D0
         d2x(js)=0.D0
      END DO

!****** bound mass calculation 
      DO js=1,NS2
         dx(js)=XSD(js)-XSF(js)
      END DO

      CALL SPLINE(NS2,XSF,dx,d2x)

      DO jb=1,NB2
         CALL SPLINT(NS2,XSF,dx,d2x,XBF(jb),dxb)
         XBD(jb)=XBF(jb)+dxb
      END DO

!!****** nucleation bound
      XBD(1)=XBF(1)
      ! the first size is kept equal to that of
      ! nucleation : an improvement would be 
      ! to use the computed size given by
      ! nucleation routines, but coagulation
      ! coefficients should be re-computed

!****** check bound order
      CALL ASCORDER(NB2,XBD,ixbd)

      itest=0
      DO jb=1,NB2
         IF (ixbd(jb).NE.jb) THEN
            dvar(1)=DBLE(ixbd(jb))
            dvar(2)=DBLE(jb)

            itest=itest+1
         ENDIF
      END DO

      ! if there had cross bounds then
      ! recompute them with a less accurate but
      ! safer algorithm ( i.e. which cannot cross )
      IF (itest.GT.0) THEN
         dvar(1)=DBLE(itest)

         DO jb=2,NS2
            XBD(jb)=(XSD(jb-1)+XSD(jb))*5.D-01
         END DO

         XBD(1)=2.D0*XSD(1)-XBD(2)
         XBD(NB2)=2.D0*XSD(NS2)-XBD(NS2)
      ENDIF

!****** sctn width
      ! this is a last test in case
      ! 0 width sctn occur, or something else
      ! that would have gone through later tests
      itest=0
      DO js=1,NS2
         HSD(js)=XBD(js+1)-XBD(js)

         IF (HSD(js).LE.0.D0) THEN
            dvar(1)=DBLE(js)
            !write(*,*) 'w:sizebnd:bin n',dvar,HSD(js)
            !HSD(js) = abs(HSD(js))
            itest=itest+1
         ENDIF
      END DO

      IF (itest.GT.0) THEN
         dvar(1)=DBLE(itest)
         ! write(*,*)'ws:sizebnd:sctn width <=0',dvar
      ENDIF

!****** compute mass and diameter
      ! once we are sure of bounds order
      ! we compute associated mass and diam
  !    DO jb=1,NB2
  !       MBD(jb)=DEXP(XBD(jb))
  !       DBD(jb)=(MBD(jb)/RHOA/CST)**FRAC3
  !    END DO
!**************************************************
      END SUBROUTINE SIZEBND
!**************************************************
!*************************************************
!*                                               *
!*     ROUTINE COMPUTING 2nd DERIVATIVES OF y(*) *
!*     AT x(*) POINTS FOR NATURAL CUBIC SPLINES  *
!*                                               *
!*************************************************
      SUBROUTINE SPLINE(n,x,y,y2)
!*************************************************
!*                                               *
!*     n     (IN)  number of points              *
!*     x(*)  (IN)  points of collocation         *
!*     y(*)  (IN)  value at points               *
!*     y2(*) (OUT) 2nd derivative of y(*)        *
!*                                               *
!*************************************************
      INTEGER n
      DOUBLE PRECISION x(n),y(n),y2(n)
!******
      INTEGER i,k
      DOUBLE PRECISION p,qn,sig,un,u(n)
     
      y2(1)=0.D0
      u(1)=0.D0

      DO i=2,n-1
         sig=(x(i)-x(i-1))/(x(i+1)-x(i-1))
         p=sig*Y2(i-1)+2.D0
         y2(i)=(sig- 1.D0)/p
         u(i)=(y(i+1)-y(i))/(x(i+1)-x(i)) &
             -(y(i)-y(i-1))/(x(i)-x(i-1))
         u(i)=u(i)*6.D0/(x(i+1)-x(i-1))
         u(i)=(u(i)-sig*u(i-1))/p
      END DO
      
      qn=0.D0
      un=0.D0
      
      y2(N)=(un-qn*u(n-1))/(qn*y2(n-1)+1.D0)
      
      DO k=n-1,1,-1
         y2(k)=y2(k)*y2(k+1)+u(k)
      END DO
!*************************************************
      END SUBROUTINE SPLINE
!*************************************************
!**************************************************
!*                                                *
!*     ROUTINE PERFORMING CUBIC SPLINE EVALUATION *
!*                                                *
!**************************************************
      SUBROUTINE SPLINT(n,xa,ya,y2a,x,y)
!**************************************************
!*                                                *
!*     n       (IN)   integer array size          *
!*     xa(*)   (IN)   abscisse function           *
!*     ya(*)   (IN)   function value at abscisse  *
!*     y2a(*)  (IN)   2nd deriv of function at xa *   
!*     x       (IN)   abscisse where to evaluate  *
!*     y       (OUT)  function evaluation at x    *
!*                                                *
!**************************************************
      IMPLICIT NONE
      
      INTEGER n
      DOUBLE PRECISION x,y,xa(n),y2a(n),ya(n)
!******
      INTEGER k,khi,klo
      DOUBLE PRECISION a,b,h
      
      klo=1
      khi=n

      DO WHILE (khi-klo.GT.1)
         k=(khi+klo)/2
         IF (xa(k).GT.x) THEN
            khi=k
         ELSE
            klo=k
         ENDIF
      END DO

      h=xa(khi)-xa(klo)
      IF (h.EQ.0.D0) THEN
         write(*,*) '< bad xa input in splint >'
         stop
      ENDIF

      a=(xa(khi)-x)/h
      b=(x-xa(klo))/h
      y= a*ya(klo)+b*ya(khi)  +( a*(a*a-1.D0)*y2a(klo)&
         +b*(b*b-1.D0)*y2a(khi) ) *(h*h)/6.D0
!**************************************************
      END SUBROUTINE SPLINT
!**************************************************
!**************************************************
!*                                                *
!*     ROUTINE COMPUTING THE ASCENDANT            *
!*     ORDER OF A GIVEN ARRAY                     *
!*                                                *
!**************************************************
      SUBROUTINE ASCORDER(n,x,ix)
!**************************************************
!*                                                *
!*     n      (IN)    array size                  *
!*     x(*)   (IN)    double array to order       *
!*     ix(*)  (OUT)   order in x(*) array         *
!*                                                *
!*     .e.g. ix(3)=4 means that x(4) is the       *
!*     3rd value in ascendant order, then         *            
!*     ix(1)=j => x(j) min, ix(n)=j => x(j) max   *
!*                                                *
!**************************************************
      IMPLICIT NONE

      INTEGER n,ix(n)
      DOUBLE PRECISION x(n)
!******
      INTEGER j,k,kmin,tix
      DOUBLE PRECISION xx(n),xmin,txx

!****** init
      DO j=1,n
         xx(j)=x(j)
         ix(j)=j
      END DO

!****** compute ascendant order on xx(*)
      DO j=1,n-1
         kmin=n
         xmin=xx(n)

         ! search min between j to n-1 values
         DO k=j,n-1
            IF (xx(k).LT.xmin) THEN
               kmin=k
               xmin=xx(k)
            ENDIF
         END DO

         ! save value for next loop
         tix=ix(j)
         txx=xx(j)

         ! store order
         ix(j)=ix(kmin)
         xx(j)=xmin

         ! exchange j and kmin values
         ! in ix(*) and xx(*)
         ix(kmin)=tix
         xx(kmin)=txx
      END DO
!**************************************************
      END SUBROUTINE ASCORDER
!**************************************************
!**************************************************
      SUBROUTINE REDIST_EULERCOUPLE(ns,naer,n,q,fixed_density,rho,dbound,diam)
!**************************************************
!*     *
!*     ns   (IN)    running number of eqs.       *
!*     q(*)  (INOUT) gas/aero conc (µg.m-3)       *
!*     *
!**************************************************
      IMPLICIT NONE

      include '../Module/CONST_B.inc'
      INTEGER ns,naer,NB
      DOUBLE PRECISION q(ns,naer),n(ns),qt(ns),diam(ns)
      DOUBLE PRECISION CST,Nnew(ns)
      INTEGER ji,jj,jesp,js1,js2,icpt,j1hi
      DOUBLE PRECISION qnew(ns,naer),qtnew(ns)
      DOUBLE PRECISION xx,tmp1,tmp2,frac
      DOUBLE PRECISION fixed_density,rho(ns),diam1,diam2
      DOUBLE PRECISION dbound(ns+1),XSD(ns),XBF(ns+1),XBD(ns+1)
      DOUBLE PRECISION x2hi,HSD(ns),XSF(ns)
      DOUBLE PRECISION DSF(ns),MSF(ns+1),MSD(ns+1),Qtotal

      NB = NS+1
      CST  =0.523598775598D0

!******zero init
      DO js1=1,NS
         qt(js1) = 0.d0
         Do jj=1,naer-1
             qt(js1) = qt(js1) + q(js1,jj)
             qnew(js1,jj)=0.D0
         Enddo
         qnew(js1,naer)=0.D0
         qtnew(js1)=0.D0
         Nnew(js1) = 0.D0
      END DO
    
      DO js1=1,NS+1
         XBF(js1)=DLOG(fixed_density*CST*dbound(js1)**3)
      ENDDO
      DO js1=1,NS
         XSF(js1) = (XBF(js1)+XBF(js1+1))*5.D-01
         if(N(js1).GT.TINYN) then
            XSD(js1) = DLOG(qt(js1)/N(js1))
            MSD(js1) = DEXP(XSD(js1))
         else
            XSD(js1) = XSF(js1)
            MSD(js1) = DEXP(XSD(js1))
         endif
         MSF(js1) = DEXP(XSF(js1))
         DSF(js1) = (dbound(js1)*dbound(js1+1))**0.5d0
      ENDDO

      call SIZEBND(ns,XSD,XSF,XBF,HSD,XBD)

      DO js2=1,NS

         x2hi=(MSD(js2)/(rho(js2)*CST))**(1.d0/3.D0) ! lagrangian diameter
         CALL LOCATE(NS,DSF,x2hi,j1hi)

         IF (N(js2).GT.TINYN) THEN
                                ! redistribute over fixed sections
            IF (j1hi.EQ.0) THEN ! number change
               Qtotal = 0.D0
               DO jesp=1,naer-1
                  !!qnew(1,jesp)=q(js2,jesp) !+ qnew(1,jesp)+q(js2,jesp)
                  Qtotal = Qtotal+qnew(1,jesp)
               END DO
               jesp = naer
               !qnew(1,jesp)=q(js2,jesp) !qnew(1,jesp)+q(js2,jesp)
               Nnew(1) = Qtotal/MSF(1)
            ELSEIF (j1hi.EQ.NS) THEN ! number change
               Qtotal = 0.D0
               DO jesp=1,naer-1
                  qnew(NS,jesp)=qnew(NS,jesp)+q(js2,jesp)
                  Qtotal = Qtotal+qnew(NS,jesp)
               END DO
               jesp = naer
               qnew(NS,jesp)=qnew(NS,jesp)+q(js2,jesp)
               Nnew(NS) = Qtotal/MSF(NS)
            ELSE
               diam1=DSF(j1hi)
               diam2=DSF(j1hi+1)
               DO jesp=1,naer
                  qnew(j1hi,jesp)=qnew(j1hi,jesp)+ q(js2,jesp)*(1-(diam2/x2hi)**3.D0)/ &
                      (1-(diam2/diam1)**3.D0)
                  qnew(j1hi+1,jesp)=qnew(j1hi+1,jesp)+ q(js2,jesp)*(1-(diam1/x2hi)**3.D0)/ &
                      (1-(diam1/diam2)**3.D0)
               END DO
               Nnew(j1hi)=Nnew(j1hi)+N(js2)*(x2hi**3-diam2**3.D0)/ &
                   (diam1**3-diam2**3)
               Nnew(j1hi+1)=Nnew(j1hi+1)+N(js2)*(x2hi**3-diam1**3.D0)/ &
                   (diam2**3-diam1**3)

            ENDIF
         ENDIF
      END DO

!     ******check if quantities are not too small
      icpt=1                    ! to enter in the loop
      icpt = 0
      DO WHILE (icpt.GT.0)
         icpt = 0

         DO js1=1,NS
            DO jesp=1,naer-1
               qtnew(js1)=qtnew(js1)+qnew(js1,jesp)
            END DO

            tmp1= Nnew(js1)*( Nnew(js1)-TINYN)
            tmp2=qtnew(js1)*(qtnew(js1)-TINYM)

            IF ( tmp1.LT.0.D0.OR.tmp2.LT.0.D0 ) THEN

               IF (Nnew(js1).EQ.0.d0) THEN
                  WRITE(6,*)'SIREAM (redist.f): (1)Q<0 ',Nnew(js1)
                  STOP
               ENDIF
               IF ((qtnew(js1)/Nnew(js1)).LE.0.d0) THEN
                  WRITE(6,*)'SIREAM (redist.f): (2)Q<0 ', qtnew(js1)/Nnew(js1)
                  STOP
               ENDIF

               xx=DLOG(qtnew(js1)/Nnew(js1))
               frac=(xx-XBD(js1))/HSD(js1)

               js2=js1+1
               IF (frac.LT.5.D-01) js2=js1-1

               IF (js1.EQ.1)   js2=js1+1
               IF (js1.EQ.NS) js2=js1-1

               Nnew(js2)=Nnew(js2)+Nnew(js1)
               qtnew(js2)=qtnew(js2)+qtnew(js1)
               Nnew(js1)=0.d0
               qtnew(js1)=0.d0

               DO jesp=1,naer
                  qnew(js2,jesp)=qnew(js2,jesp)+qnew(js1,jesp)
                  qnew(js1,jesp)=0.d0
               END DO
               icpt=icpt+1
            ENDIF
         END DO

      END DO

!******turn back to conc vector
      
      DO js1=1,NS
         Do jj=1,naer
             q(js1,jj)=qnew(js1,jj)
         Enddo
         N(js1) = Nnew(js1)
      END DO

  END SUBROUTINE REDIST_EULERCOUPLE
!**************************************************
  SUBROUTINE REDIST_MOVINGDIAM(ns,naer,n,q,fixed_density,rho,dbound,diam)
!**************************************************
!*     *
!*     ns   (IN)    running number of eqs.       *
!*     q(*)  (INOUT) gas/aero conc (µg.m-3)       *
!*     *
!**************************************************
      IMPLICIT NONE

      include '../Module/CONST_B.inc'
      INTEGER ns,naer,NB
      DOUBLE PRECISION q(ns,naer),n(ns),qt(ns),diam(ns)
      DOUBLE PRECISION CST,Nnew(ns)
      INTEGER ji,jj,jesp,js1,js2,icpt,j1hi
      DOUBLE PRECISION qnew(ns,naer),qtnew(ns)
      DOUBLE PRECISION xx,tmp1,tmp2,frac
      DOUBLE PRECISION fixed_density,rho(ns),diam1,diam2
      DOUBLE PRECISION dbound(ns+1),XSD(ns),XBF(ns+1),XBD(ns+1)
      DOUBLE PRECISION x2hi,HSD(ns),XSF(ns)
      DOUBLE PRECISION DSF(ns),MSF(ns+1),MSD(ns+1),Qtotal

      NB = NS+1
      CST  =0.523598775598D0

!******zero init
      DO js1=1,NS
         qt(js1) = 0.d0
         Do jj=1,naer-1
             qt(js1) = qt(js1) + q(js1,jj)
             qnew(js1,jj)=0.D0
         Enddo
         qnew(js1,naer)=0.D0
         qtnew(js1)=0.D0
         Nnew(js1) = 0.D0
      END DO
    
      DO js1=1,NS+1
         XBF(js1)=DLOG(fixed_density*CST*dbound(js1)**3)
      ENDDO
      DO js1=1,NS
         XSF(js1) = (XBF(js1)+XBF(js1+1))*5.D-01
         if(N(js1).GT.TINYN) then
            XSD(js1) = DLOG(qt(js1)/N(js1))
            MSD(js1) = DEXP(XSD(js1))
         else
            XSD(js1) = XSF(js1)
         endif
         MSF(js1) = DEXP(XSF(js1))
         DSF(js1) = (dbound(js1)*dbound(js1+1))**0.5d0
      ENDDO

      call SIZEBND(ns,XSD,XSF,XBF,HSD,XBD)

      DO js2=1,NS

         !x2hi=(MSD(js2)/(rho(js2)*CST))**(1.d0/3.D0) ! lagrangian diameter
         x2hi=XSD(js2) ! lagrangian diameter
         CALL LOCATE(NS,XBF,x2hi,j1hi)

         IF (N(js2).GT.TINYN) THEN
                                ! redistribute over fixed sections
            IF (j1hi.EQ.0) THEN ! number change
               Qtotal = 0.D0
               DO jesp=1,naer-1
                  !!qnew(1,jesp)=qnew(1,jesp)+q(js2,jesp)
                  Qtotal = Qtotal+qnew(1,jesp)
               END DO
               jesp = naer
               !!qnew(1,jesp)=qnew(1,jesp)+q(js2,jesp)
               Nnew(1) = Qtotal/MSF(1)
            ELSEIF (j1hi.EQ.NS) THEN ! number change
               Qtotal = 0.D0
               DO jesp=1,naer-1
                  qnew(NS,jesp)=qnew(NS,jesp)+q(js2,jesp)
                  Qtotal = Qtotal+qnew(NS,jesp)
               END DO
               jesp = naer
               qnew(NS,jesp)=qnew(NS,jesp)+q(js2,jesp)
               Nnew(NS) = Qtotal/MSF(NS)
            ELSE
               DO jesp=1,naer
                  qnew(j1hi,jesp)=qnew(j1hi,jesp)+ q(js2,jesp)
               END DO
               Nnew(j1hi)=Nnew(j1hi)+N(js2)

            ENDIF
         ENDIF
      END DO

!     ******check if quantities are not too small
!      icpt=1                    ! to enter in the loop
!      icpt = 0
!      DO WHILE (icpt.GT.0)
!         icpt = 0
!
!         DO js1=1,NS
!            DO jesp=1,naer-1
!               qtnew(js1)=qtnew(js1)+qnew(js1,jesp)
!            END DO
!
!            tmp1= Nnew(js1)*( Nnew(js1)-TINYN)
!            tmp2=qtnew(js1)*(qtnew(js1)-TINYM)
!
!            IF ( tmp1.LT.0.D0.OR.tmp2.LT.0.D0 ) THEN
!
!               IF (Nnew(js1).EQ.0.d0) THEN
!                  WRITE(6,*)'SIREAM (redist.f): (1)Q<0 ',Nnew(js1)
!                  STOP
!               ENDIF
!               IF ((qtnew(js1)/Nnew(js1)).LE.0.d0) THEN
!                  WRITE(6,*)'SIREAM (redist.f): (2)Q<0 ', qtnew(js1)/Nnew(js1)
!                  STOP
!               ENDIF
!
!               xx=DLOG(qtnew(js1)/Nnew(js1))
!               frac=(xx-XBD(js1))/HSD(js1)
!
!               js2=js1+1
!               IF (frac.LT.5.D-01) js2=js1-1
!
!               IF (js1.EQ.1)   js2=js1+1
!               IF (js1.EQ.NS) js2=js1-1
!
!               Nnew(js2)=Nnew(js2)+Nnew(js1)
!               qtnew(js2)=qtnew(js2)+qtnew(js1)
!               Nnew(js1)=0.d0
!               qtnew(js1)=0.d0
!
!               DO jesp=1,naer
!                  qnew(js2,jesp)=qnew(js2,jesp)+qnew(js1,jesp)
!                  qnew(js1,jesp)=0.d0
!               END DO
!               icpt=icpt+1
!            ENDIF
!         END DO
!
!      END DO

!******turn back to conc vector
      
      DO js1=1,NS
         Do jj=1,naer
             q(js1,jj)=qnew(js1,jj)
         Enddo
         N(js1) = Nnew(js1)
      END DO

  END SUBROUTINE REDIST_MOVINGDIAM
!**************************************************
