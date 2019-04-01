c-----------------------------------------------------------------------
c    include 'cmn_metdat.f'  for fast-JX code v5.3 (prather 6/05)
c
c        needs 'parm_ctm.f' for dimensions
c        delivers p, T, Surf Albedo, and Optical Depth from CTM to fastJX
c        >>>>this is for standalone fast-JX ver 5.3   (6.05)
c-----------------------------------------------------------------------
c
c      integer  IMAX
      real  P(IMAX,JMAX)         !  Surface pressure
      real  T(IMAX,JMAX,LMAX)      !  Temperature profile
      real  OD(IMAX,JMAX,LWEMAX)   !  Optical Depth profile
c
      real  XGRD(IMAX)         !  Longitude (midpoint, radians)
      real  XDGRD(IMAX)
      real  YGRD(JMAX)         !  Latitude  (midpoint, radians)
      real  YDGRD(JMAX)
      real  ETAA(LMAX+1)       !  Eta(a) value for level boundaries
      real  ETAB(LMAX+1)       !  Eta(b) value for level boundaries
      real  AREAXY(IMAX,JMAX)    !  area (m^2)
      real  ZCTM(LMAX+1)         !  Altitudes
      integer  MONTH
      integer  NSLAT         ! Latitude(J) index of current column
      integer  NSLON         ! Longitude(I) index of current column

      real, dimension(IMAX,JMAX,LMAX+1) :: TJ, DM, DO3, ZH
      real, dimension(IMAX,JMAX,LMAX+1) :: DAER1, DAER2, DAER3
      real, dimension(IMAX,JMAX,LMAX+1) :: ODCLD, ICECLD
      integer,dimension(IMAX,JMAX,LMAX+1) :: NAER1, NAER2, NAER3
      integer,dimension(IMAX,JMAX,LMAX+1) :: NCLDX, NICEX
      real, dimension(IMAX,JMAX)     :: PMEAN, SA
c
      real   STT(IMAX,JMAX,LMAX,NTR_)
      real   TREF(51,18,12),OREF(51,18,12)
      character*10  TCNAME(NTR_)
c
      common/metdat/P,T,OD, STT,
     &  XGRD,YGRD,XDGRD,YDGRD,ETAA,ETAB,AREAXY,ZCTM,
     &  TREF,OREF,PMEAN,SA, MONTH,NSLAT,NSLON, TCNAME
c
      common /jvdatIJ/TJ,DM,DO3,ZH,
     &                  DAER1,DAER2,DAER3,ODCLD,ICECLD,
     &                  NAER1,NAER2,NAER3,NCLDX,NICEX
c-----------------------------------------------------------------------

