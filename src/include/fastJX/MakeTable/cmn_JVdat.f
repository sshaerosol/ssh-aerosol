c-----------------------------------------------------------------------
c    include 'cmn_JVdat.f'  for fast-JX code v5.3 (prather 6/05)
c 
c NB - ALL of these common variables are set paramters, 
c    They are NOT to be used as variables for a local solution
c    Thus this entire set is 'in' only after it is initialized
c-----------------------------------------------------------------------
      real*8  RAD,ZZHT,ATAU,ATAU0
      real*8  WBIN(W_+1),WL(W_),FL(W_),QO2(W_,3),QO3(W_,3),Q1D(W_,3)
      real*8  QQQ(W_,2,X_),QRAYL(W_+1),TQQ(3,X_)
      real*8  WAA(5,A_),QAA(5,A_),PAA(8,5,A_),RAA(5,A_),SSA(5,A_)
      real*8  JFACTA(JVN_)
      integer JIND(JVN_),NRATJ,NJVAL,NW1,NW2,NAA,JTAUMX
      character*20 TITLEA(A_)
      character*78 TITLE0
      character*7  TITLEJ(X_),TITLEJ2,TITLEJ3
      character*10 JLABEL(JVN_)
c
      common /jvchem/JFACTA,JIND,NRATJ,JLABEL,TITLEA 
c
      common /jvdat/WBIN,WL,FL,QO2,QO3,Q1D,QQQ,QRAYL,TQQ,
     &             WAA,QAA,PAA,RAA,SSA,RAD,ZZHT,ATAU,ATAU0,
     &             JTAUMX, NJVAL,NW1,NW2,NAA ,TITLE0,TITLEJ
c-----------------------------------------------------------------------









