      implicit none
      integer     NF,I,J,JQ,K,L,N,M,N1,N2,IW,NWWW, NT
      integer     NW1,NW2,NWSRB,NSR,NODF,ISR,NBIN,NSPEC,NQQQ,NJVAL
      real*8      WBIN,WL,FL,QO2,O2X,ODF,FNO,QNO,QO3,Q1D,QQQ,TQQ,QRL
      real*8      CNO

      character*20 TITLEJ(3,99)
      character*78 TITLE0,TITLE
      common / TITLES/ TITLE0,TITLE,TITLEJ
      common /CCWVL/ WBIN(77+1),WL(77),FL(77),QO2(77,3),O2X(6,15,3),
     & ODF(6,15),FNO(15),QNO(6,15),  QO3(77,3),Q1D(77,3),QQQ(77,2,99),
     & TQQ(3,99),QRL(77),
     & NW1,NW2,NWSRB,NSR,NODF,ISR(15),NT,NBIN,NSPEC, NQQQ,NJVAL

c  Re-distribution common block (expand SR ODFs to get 145 W bins)
      real           FFL(145),QQNO(145),QQRL(145),WWL(145)
      real           QQO3(145,3),QQ1D(145,3),QQO2(145,3),
     &               QQQQ(145,2,99)
      common /REDIS/  FFL,QQNO,QQRL,WWL,QQO3,QQ1D,QQO2,QQQQ
c
c  Average into Fast-J2/X 18 bins  
      real           SFL(18),SQNO(18),SQRL(18),SWL(18)
      real           SQO3(18,3),SQ1D(18,3),SQO2(18,3),
     &               SQQQ(18,2,95)
      common /AVGJ2/  SFL,SQNO,SQRL,SWL,SQO3,SQ1D,SQO2,SQQQ

