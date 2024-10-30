!-----------------------------------------------------------------------
!     Copyright (C) 2024 CEREA (ENPC) - INERIS
!     SSH-aerosol is distributed under the GNU General Public License v3
!-----------------------------------------------------------------------


MODULE mod_cubicspline

  IMPLICIT NONE

CONTAINS

!-----------------------------------------------------------------------
! Over the interval xp(i) to xp(i+1), the interpolating polynomial               
!  y = yp(i) + a(i)*z + b(i)*z**2 + c(i)*z**3     
! where z=(x-xp(i))/(xp(i+1)-xp(i)) is used. 
! The coefficients a(i),b(i) and c(i) are computed by cspline and stored
! in locations coe(3*i-2),coe(3*i-1) and coe(3*i) respectively. 
! While working in the ith interval,the variable xdel will represent 
! xdel = xp(i+1) - xp(i), and yp(i) will represent yp(i+1) - yp(i).
!-----------------------------------------------------------------------
SUBROUTINE ssh_gck_cspline(ndat,xp,yp,yout)
 
  IMPLICIT NONE

  INTEGER, INTENT(IN) :: ndat ! nsza
  DOUBLE PRECISION, INTENT(IN)  :: xp(ndat)  ! szas(nsza)
  DOUBLE PRECISION, INTENT(IN)  :: yp(ndat)  ! a_tmp(nsza) read from file
  DOUBLE PRECISION, INTENT(OUT) :: yout(4,ndat-1) ! tmp_inter(4, nsza-1)

! local 
  DOUBLE PRECISION     :: w(ndat*3), coe(ndat*3)
  DOUBLE PRECISION     :: deriv(2)      

  DOUBLE PRECISION    :: ydel, xdel, ratio, xdel1, zp
  INTEGER :: m, ii, i, j
   
! Initialize 
  j=1          ! derivative at end point provided
  deriv(:)=0.  ! value of derivative at the end points (set to 0)
  coe(:)=0.

  xdel=xp(2)-xp(1)                
  ydel=yp(2)-yp(1)              

  IF (j==2) THEN 
!-----------------------------------------------------------------------
! If the second derivative at the end points is given b(1) is known, 
! the second equation becomes a(1)+c(1)=ydel-0.5*xdel*xdel*deriv(1). 
! During the solution of the 3ndat-4 equations, a(1) will be kept in 
! cell coe(2) instead of coe(1) to retain the tridiagonal form of the 
! coefficient matrix.    
!-----------------------------------------------------------------------
    coe(2)=0.0                     
    w(2)=0.5*xdel*xdel*deriv(1)       

  ELSE     
! ----------------------------------------------------------------------
! If the first derivative at the end points is given, a(1) is known, 
! and the second equation becomes merely b(1)+c(1)=ydel-xdel*deriv(1).     
! ----------------------------------------------------------------------
    coe(1)=xdel*deriv(1)      
    coe(2)=1.0     
    w(2)=ydel-coe(1)
  ENDIF      

  m=ndat-2     
  IF (m > 0) THEN     
!-----------------------------------------------------------------------
! Upper triangularization of the tridiagonal system of equations for  
! the coefficient matrix follows.    
!-----------------------------------------------------------------------
    DO i=1,m                     
      xdel1=xdel      
      xdel=xp(i+2)-xp(i+1)      
      ratio=xdel1/xdel 
      coe(3*i)=-ratio/(2.0-coe(3*i-1))     
      w(3*i)=(-ydel-w(3*i-1))/(2.0-coe(3*i-1))   
      coe(3*i+1)=-ratio*ratio/(ratio-coe(3*i))     
      w(3*i+1)=(ydel-w(3*i))/(ratio-coe(3*i))    
      ydel=yp(i+2)-yp(i+1)       
      coe(3*i+2)=1.0/(1.0-coe(3*i+1))   
      w(3*i+2)=(ydel-w(3*i+1))/(1.0-coe(3*i+1))
    ENDDO
  ENDIF    

!-----------------------------------------------------------------------
! c(ndat-1) is determined directly from the last equation obtained  
! above, and the first or second derivative value given at the end point.     
!-----------------------------------------------------------------------
  IF (j==1) THEN     
    coe(3*ndat-3)=(xdel*deriv(2)-ydel-w(3*ndat-4))/(2.0-coe(3*ndat-4))    
  ELSE
    coe(3*ndat-3)=(xdel*xdel*deriv(2)/2.0 - w(3*ndat-4)) / &
                  (3.0-coe(3*ndat-4))  
  ENDIF

  m=3*ndat-6       
  IF (m > 0) THEN    
! ----------------------------------------------------------------------
! Back solution for all coefficents except a(1) and b(1) follows.    
! ----------------------------------------------------------------------
    DO ii=1,m      
      i=m-ii+3     
      coe(i)=w(i)-coe(i)*coe(i+1)
    ENDDO
  ENDIF     

  IF (j==1) THEN      
    coe(2)=w(2)-coe(3)      
!-----------------------------------------------------------------------
! If the second derivative is given at the end points, a(1) can now be 
! computed from the known values of b(1) and c(1). Then a(1) and b(1) 
! are put into their proper places in the coe array.   
!-----------------------------------------------------------------------
  ELSE
    coe(1)=yp(2)-yp(1)-w(2)-coe(3)     
    coe(2)=w(2)     
  ENDIF

! from csplint(ndat,xp,yp,coe,xin,yout,ierr)
! ---------------------------------------------------------------------
! Given xin, return yout according to the 3 order polynomial : 
!   y = yp(i) + A(i)*zp + B(i)*zp**2 + C(i)*zp**3 
! where zp=(x-xp(i))/(xp(i+1)-xp(i)) 
!
! The polynomial coefficients are compute by cspline
! ---------------------------------------------------------------------

  ! compute yout
  DO i=2,ndat

      xdel = xp(i)-xp(i-1)                   
      !zp=(xin-xp(i-1))/xdel                   
      !yout(i-1)=((zp*coe(3*i-3) + coe(3*i-4))*zp + coe(3*i-5))*zp + yp(i-1)
      do j = 1, 3
          yout(j,i-1) = coe(3*i-2-j)
          !if (i.eq.9) print*,j,yout(i-1,j)
      enddo
      yout(4,i-1) = yp(i-1)
  ENDDO

END SUBROUTINE ssh_gck_cspline                      

subroutine ssh_SPL3(n1,a,b,c)

    implicit none
    
    INTEGER, INTENT(IN) :: n1
    integer i,n0
     
    DOUBLE PRECISION, INTENT(IN) :: a(:), b(:)
    DOUBLE PRECISION, INTENT(OUT) :: c(:,:)

    double precision :: c20(n1),c21(n1)
    double precision :: c30(n1),d(n1)
    double precision :: c31(n1),c40(n1)
    double precision :: c41(n1),cc31

      n0=n1-1
      do i = 1,n0
         c(1,i) = b(i)
         d(i) = a(i+1)-a(i)
      enddo
      c20(1) = 0.D0
      c21(1) = 0.D0
      c30(1) = 0.D0
      c31(1) = 1.D0
      c40(1) = (b(2  )-c(1,1)-c20(1)*d(1)-c30(1)*d(1)**2)/d(1)**3
      c41(1) = (             -c21(1)*d(1)-c31(1)*d(1)**2)/d(1)**3
      do i = 2,n0
      
         c20(i) = c20(i-1)+2.D0*c30(i-1)*d(i-1)+3.D0*c40(i-1)*d(i-1)**2
         c21(i) = c21(i-1)+2.D0*c31(i-1)*d(i-1)+3.D0*c41(i-1)*d(i-1)**2

         c30(i) = c30(i-1) +3.D0*c40(i-1)*d(i-1)
         c31(i) = c31(i-1) +3.D0*c41(i-1)*d(i-1)


         c40(i) = (b(i+1)-c(1,i)-c20(i)*d(i)-c30(i)*d(i)**2)/d(i)**3
         c41(i) = (             -c21(i)*d(i)-c31(i)*d(i)**2)/d(i)**3
      enddo

      cc31 = c20(n0)+2.D0*c30(n0)*d(n0)+3.D0*c40(n0)*d(n0)**2
      cc31 = -cc31/(c21(n0)+2.D0*c31(n0)*d(n0)+3.D0*c41(n0)*d(n0)**2)

      do i = 1,n0
         c(2,i) = c20(i)+cc31*c21(i)
         c(3,i) = c30(i)+cc31*c31(i)
         c(4,i) = c40(i)+cc31*c41(i)
      enddo

END subroutine ssh_SPL3

END MODULE mod_cubicspline
