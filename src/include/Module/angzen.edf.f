      
C     **********************************************************************
C
      DOUBLE PRECISION function ssh_tsolv(tu,long)
C
C     Routine DIFEUL    DE D.Wendum /EDF
C     **********************************************************************
C
C     calcul de l'heure solaire "vraie"
C     fonction de l'heure tu (gmt)
C     exprimee en secondes dans l'annee
C     et de la longitude en degres
C
C     retour = valeur en heures
C
      implicit none
C
      DOUBLE PRECISION tu,long
      DOUBLE PRECISION t
      DOUBLE PRECISION ssh_aent,DLaent,ssh_amodp
      DOUBLE PRECISION ssh_corheu,DLcorheu
      DOUBLE PRECISION pi
      parameter(pi=3.141592653589793238462D0)
      external ssh_aent,ssh_amodp,ssh_corheu
C
      t=(tu/86400.D0)
      DLaent=ssh_aent(t)
      t=t-DLaent
      t=t*24.D0
      DLcorheu=ssh_corheu(tu)
      ssh_tsolv=t+long/15.D0+DLcorheu
      ssh_tsolv=ssh_amodp(ssh_tsolv,24.D0)
      ssh_tsolv=ssh_tsolv-12.D0
      ssh_tsolv=ssh_tsolv*pi/12.D0
C
      end

C
C     **********************************************************************
C
      DOUBLE PRECISION function ssh_declin(tsec)

C     Routine DIFEUL    DE D.Wendum /EDF
C     **********************************************************************
C
C     calcul de la declinaison du soleil fonction
C     du temps tu en secondes dans l'annee
C
      implicit none
C
      DOUBLE PRECISION tsec,t
      DOUBLE PRECISION pi
      parameter(pi=3.141592653589793238462D0)
C
      t=1.D0+(tsec/86400.D0)
      t=t+0.1D0
      t=2*pi*t/365.D0
      ssh_declin=0.006918D0
      ssh_declin=ssh_declin-0.399912D0*dcos(t)+0.070257D0*dsin(t)
      ssh_declin=ssh_declin-0.006758D0*dcos(2.D0*t)+
     &     0.000907D0*dsin(2.D0*t)
      ssh_declin=ssh_declin-0.002697D0*dcos(3.D0*t)+
     &     0.001480D0*dsin(3.D0*t)
C
      end

C     **********************************************************************
C
      DOUBLE PRECISION function ssh_muzero(tu,long,lat)
C
C     Routine DIFEUL    DE D.Wendum /EDF
C     **********************************************************************
C
C     calcul du cosinus de l'angle zenithal
C     fonction de l'heure tu (gmt)
C     exprimee en secondes dans l'annee
C     et de la longitude,latitude en degres
C
      implicit none
C
      DOUBLE PRECISION tu,tul,long,hr,ssh_tsolv
      DOUBLE PRECISION decl,ssh_declin
      DOUBLE PRECISION lat,flat
      DOUBLE PRECISION pi
      parameter(pi=3.141592653589793238462D0)
      external ssh_declin,ssh_tsolv
C
      tul=tu
      if (tu.gt.31536000.d0) then
         tul=mod(tu,31536000.d0)
      endif

      hr=ssh_tsolv(tul,long)
      decl=ssh_declin(tul)
      flat=lat*pi/180.D0
      ssh_muzero=dsin(decl)*dsin(flat)+dcos(decl)*dcos(flat)*dcos(hr)
C
      end
C
C     **********************************************************************
C
      DOUBLE PRECISION function ssh_corheu(tsec)
C
C     Routine DIFEUL    DE D.Wendum /EDF
C     **********************************************************************
C     calcul de la difference temps "vrai" - temps "moyen"
C     fonction du temps tu en secondes dans l'annee
C
C     retour = valeur en heures
C
      implicit none
C
      DOUBLE PRECISION tsec,t
      DOUBLE PRECISION pi
      parameter(pi=3.141592653589793238462D0)
C
      t=1.D0+(tsec/86400.D0)
      t=t+.1D0
      t=2.D0*pi*t/365.D0
C
      ssh_corheu=0.000075D0
      ssh_corheu=ssh_corheu+0.001868D0*dcos(t)-0.032077D0*dsin(t)
      ssh_corheu=ssh_corheu-0.014615D0*dcos(2.D0*t)-
     &     0.040849D0*dsin(2.D0*t)
      ssh_corheu=ssh_corheu*12.D0/pi
C
      end

C     **********************************************************************
C
      DOUBLE PRECISION FUNCTION SSH_AENT(X)
C
C     PARTIE ENTIERE "USUELLE" : UNIQUE ENTIER K TEL QUE K =<X < K+1
C
      IMPLICIT NONE
      DOUBLE PRECISION X
C
      IF(X.GE.0.D0)THEN
         SSH_AENT=DINT(X)
      ELSE
         SSH_AENT=DINT(X)-1.D0
      ENDIF
C
      END
C
C     **********************************************************************
C
      DOUBLE PRECISION FUNCTION SSH_AMODP(X,Y)
C
C     L'UNIQUE REEL F TEL QUE X= N*|Y| +F, 0<= F < |Y|
C
      IMPLICIT NONE
      DOUBLE PRECISION X,Y,Z
      DOUBLE PRECISION SSH_AENT,DLAENT
      EXTERNAL SSH_AENT
C
      Z=DMAX1(Y,-Y)
      IF(Z.GT.0)THEN
         DLAENT=SSH_AENT(X/Z)
         SSH_AMODP=X-Z*DLAENT
      ELSE
         STOP 'APPEL DE SSH_AMODP(X,Y) AVEC Y=0'
      ENDIF
C
      END
C
C     **********************************************************************

