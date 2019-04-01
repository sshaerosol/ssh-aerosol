C-----------------------------------------------------------------------
C     Copyright (C) 2003-2007, ENPC - INRIA - EDF R&D
C     Author(s): Denis Wendum
C
C     This file is part of the Size Resolved Aerosol Model (SIREAM_AEC),
C     a component of the air quality modeling system Polyphemus.
C
C     Polyphemus is developed in the INRIA - ENPC joint project-team
C     CLIME and in the ENPC - EDF R&D joint laboratory CEREA.
C
C     Polyphemus is free software; you can redistribute it and/or modify
C     it under the terms of the GNU General Public License as published
C     by the Free Software Foundation; either version 2 of the License,
C     or (at your option) any later version.
C
C     Polyphemus is distributed in the hope that it will be useful, but
C     WITHOUT ANY WARRANTY; without even the implied warranty of
C     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
C     General Public License for more details.
C
C     For more information, visit the Polyphemus web site:
C     http://cerea.enpc.fr/polyphemus/
C-----------------------------------------------------------------------

C     **********************************************************************
C
      DOUBLE PRECISION function tsolv(tu,long)
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
      DOUBLE PRECISION aent,DLaent,amodp
      DOUBLE PRECISION corheu,DLcorheu
      DOUBLE PRECISION pi
      parameter(pi=3.141592653589793238462D0)
      external aent,amodp,corheu
C
      t=(tu/86400.D0)
      DLaent=aent(t)
      t=t-DLaent
      t=t*24.D0
      DLcorheu=corheu(tu)
      tsolv=t+long/15.D0+DLcorheu
      tsolv=amodp(tsolv,24.D0)
      tsolv=tsolv-12.D0
      tsolv=tsolv*pi/12.D0
C
      end

C
C     **********************************************************************
C
      DOUBLE PRECISION function declin(tsec)

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
      declin=0.006918D0
      declin=declin-0.399912D0*dcos(t)+0.070257D0*dsin(t)
      declin=declin-0.006758D0*dcos(2.D0*t)+0.000907D0*dsin(2.D0*t)
      declin=declin-0.002697D0*dcos(3.D0*t)+0.001480D0*dsin(3.D0*t)
C
      end

C     **********************************************************************
C
      DOUBLE PRECISION function muzero(tu,long,lat)
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
      DOUBLE PRECISION tu,tul,long,hr,tsolv
      DOUBLE PRECISION decl,declin
      DOUBLE PRECISION lat,flat
      DOUBLE PRECISION pi
      parameter(pi=3.141592653589793238462D0)
      external declin,tsolv
C
      tul=tu
      if (tu.gt.31536000.d0) then
         tul=mod(tu,31536000.d0)
      endif

      hr=tsolv(tul,long)
      decl=declin(tul)
      flat=lat*pi/180.D0
      muzero=dsin(decl)*dsin(flat)+dcos(decl)*dcos(flat)*dcos(hr)
C
      end
C
C     **********************************************************************
C
      DOUBLE PRECISION function corheu(tsec)
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
      corheu=0.000075D0
      corheu=corheu+0.001868D0*dcos(t)-0.032077D0*dsin(t)
      corheu=corheu-0.014615D0*dcos(2.D0*t)-0.040849D0*dsin(2.D0*t)
      corheu=corheu*12.D0/pi
C
      end

C     **********************************************************************
C
      DOUBLE PRECISION FUNCTION AENT(X)
C
C     PARTIE ENTIERE "USUELLE" : UNIQUE ENTIER K TEL QUE K =<X < K+1
C
      IMPLICIT NONE
      DOUBLE PRECISION X
C
      IF(X.GE.0.D0)THEN
         AENT=DINT(X)
      ELSE
         AENT=DINT(X)-1.D0
      ENDIF
C
      END
C
C     **********************************************************************
C
      DOUBLE PRECISION FUNCTION AMODP(X,Y)
C
C     L'UNIQUE REEL F TEL QUE X= N*|Y| +F, 0<= F < |Y|
C
      IMPLICIT NONE
      DOUBLE PRECISION X,Y,Z
      DOUBLE PRECISION AENT,DLAENT
      EXTERNAL AENT
C
      Z=DMAX1(Y,-Y)
      IF(Z.GT.0)THEN
         DLAENT=AENT(X/Z)
         AMODP=X-Z*DLAENT
      ELSE
         STOP 'APPEL DE AMODP(X,Y) AVEC Y=0'
      ENDIF
C
      END
C
C     **********************************************************************

