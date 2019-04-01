c    include 'parm_CTM.f'  for fast-JX code v5.3 (prather 6/05)
c
c     IMAX, JMAX, LMAX= longitude, latitude, altitude maximums of the grid
c     LWEMAX= altitude maximum for trop processes (clouds, rain)
c     I_ = longitude dim of CTM grid
c     J_ = latitude  dim of CTM grid
c     L_ = altitude(levels) dim of CTM grid
c     LWE_ = altitude(level) dim for trop processes (clouds, rain)
c     JVL_ = vertical(levels) dim for J-values
c     L2_  = 2*L1_ = 2*L_ + 2 = no. levels in the basic Fast-JX grid (mid-level_
c     W_   = dim = no. of Wavelength bins
c     X_   = dim = no. of X-section data sets (input data)
c     A_   = dim = no. of Aerosol/cloud Mie sets (input data)
c     MX   = no. of aerosol/cloud types supplied from CTM
c     NTR_ = no. of CTM tracers
c     SZAMAX    Solar zenith angle cut-off, above which to skip calculation
c
c-----------------------------------------------------------------------
      integer, parameter ::  IMAX=210, JMAX=115, LMAX=15  !for 13 CTM grid
      integer, parameter ::  LWEMAX= LMAX
      integer, parameter ::  Nt=1
      integer, parameter ::  W_=18, X_=51, A_=40
      integer, parameter ::  MX=4, NTR_=1
      real,  parameter ::  SZAMAX=98.0d0
      integer I_, J_, L_, LWE_, L1_, L2_, JVL_
      common /ArraySize/I_, J_, L_, LWE_, L1_, L2_, JVL_

c-----------------------------------------------------------------------

