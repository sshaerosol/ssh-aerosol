c-----------------------------------------------------------------------
c    include 'cmn_JVdat.f'  for fast-JX code v5.3 (prather 6/05)
c
c NB - ALL of these common variables are set paramters,
c    They are NOT to be used as variables for a local solution
c    Thus this entire set is 'in' only after it is initialized
c-----------------------------------------------------------------------
      real  RAD,ZZHT,ATAU,ATAU0
      real  WBIN(W_+1),WL(W_),FL(W_),QO2(W_,3),QO3(W_,3),Q1D(W_,3)
      real  QQQ(W_,2,X_),QRAYL(W_+1),TQQ(3,X_)
      real  WAA(5,A_),QAA(5,A_),PAA(8,5,A_),RAA(5,A_),SSA(5,A_)
c      integer JIND(JVN_),NRATJ,NJVAL,NW1,NW2,NAA,JTAUMX
      integer NJVAL,NW1,NW2,NAA,JTAUMX
      character*20 TITLEA(A_)
      character*78 TITLE0
      character*10  TITLEJ(X_),TITLEJ2,TITLEJ3
c
      common /jvchem/TITLEA
c
      common /jvdat/WBIN,WL,FL,QO2,QO3,Q1D,QQQ,QRAYL,TQQ,
     &     WAA,QAA,PAA,RAA,SSA,RAD,ZZHT,ATAU,ATAU0,
     &     JTAUMX, NJVAL,NW1,NW2,NAA ,TITLE0,TITLEJ,
     &     TITLEJ2,TITLEJ3
C-----------------------------------------------------------------------









