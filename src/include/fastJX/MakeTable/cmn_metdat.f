c-----------------------------------------------------------------------
c    include 'cmn_metdat.f'  for fast-JX code v5.3 (prather 6/05)
c
c        needs 'parm_ctm.f' for dimensions
c        delivers p, T, Surf Albedo, and Optical Depth from CTM to fastJX
c        >>>>this is for standalone fast-JX ver 5.3   (6.05)
c-----------------------------------------------------------------------
c
      real*8  P(I_,J_)         !  Surface pressure
      real*8  T(I_,J_,L_)      !  Temperature profile
      real*8  OD(I_,J_,LWE_)   !  Optical Depth profile
c
      real*8  XGRD(I_)         !  Longitude (midpoint, radians)
      real*8  XDGRD(I_)
      real*8  YGRD(J_)         !  Latitude  (midpoint, radians)
      real*8  YDGRD(J_) 
      real*8  ETAA(L_+1)       !  Eta(a) value for level boundaries
      real*8  ETAB(L_+1)       !  Eta(b) value for level boundaries
      real*8  AREAXY(I_,J_)    !  area (m^2)
      real*8  ZCTM(L_+1)         !  Altitudes
      integer  MONTH
      integer  NSLAT         ! Latitude(J) index of current column
      integer  NSLON         ! Longitude(I) index of current column

      real*8, dimension(I_,J_,L1_) :: TJ, DM, DO3, ZH
      real*8, dimension(I_,J_,L1_) :: DAER1, DAER2, DAER3, ODCLD, ICECLD
      integer,dimension(I_,J_,L1_) :: NAER1, NAER2, NAER3, NCLDX, NICEX
      real*8, dimension(I_,J_)     :: PMEAN, SA
c
      real*8   STT(I_,J_,L_,NTR_)
      real*8   TREF(51,18,12),OREF(51,18,12)
      character*10  TCNAME(NTR_)
c
      common/metdat/P,T,OD, STT,
     &  XGRD,YGRD,XDGRD,YDGRD,ETAA,ETAB,AREAXY,ZCTM,
     &  TREF,OREF,PMEAN,SA, MONTH,NSLAT,NSLON, TCNAME
c
      common /jvdatIJ/TJ,DM,DO3,ZH,
     &                  DAER1,DAER2,DAER3,ODCLD,
     &                  NAER1,NAER2,NAER3,NCLDX
c-----------------------------------------------------------------------

