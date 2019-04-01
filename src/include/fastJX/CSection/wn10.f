      PROGRAM WNTEN
C-----------------------------------------------------------------------
C
C  Pre-processing routine to read in supplied cross-section data and 
C  generate standardized table of cross-sections at 10-waveneumber
C  intervals. The wide range of formats used in the literature when 
C  giving cross-sections prevents this routine from being entirely
C  general - the reading section will need to be tailored to the
C  cross-section of interest. In addition, the interpolation routine
C  may need to be altered depending on the behaviour of the x-sections,
C  and care may be needed in choosing where to allow the x-sections
C  to tend to zero. This is NOT a "black-box" routine...
C
C                                      Author:        Xin Zhu
C                                      Distribution:  Oliver Wild
C                                      Date:          Sept 98
C
C Adapted by Elsa Real to read the cross-section as used in JProc (pre-
C processor used to generate photolysis rates in CMAQ and Polyphemus) 
C
C-----------------------------------------------------------------------
C     NWN     Maximum number of standard wavenumbers to use
C     NXS     Maximum number of supplied wavelengths (increase if needed)
C     WLEN    Supplied wavelengths (nm)
C     SEGM    Supplied cross-sections (cm2)
C     TIN     Temperature of supplied cross-sections (max 2 for Fast-J)
C     WN      Standard wavenumbers
C     W_REF   Standard wavelengths
C     ACS1    Interpolated cross-sections (temp 1)
C     ACS2    Interpolated cross-sections (temp 2)
C-----------------------------------------------------------------------
      implicit none
      integer        NWN, NXS, NSPEC
      parameter      (NWN=4600, NXS=240, NSPEC=1) !!!!! to modify
      common /ORI/   WN,W_REF,NMAX
      character*20   LABEL,LABEL_TAB(NSPEC) 
      character*25   FNAMIN, FNAMIN_TAB(NSPEC)
      character*34   FNAMOU,FNAMOU_TAB(NSPEC)
      double precision    WLEN(NXS,2),SEGM(NXS,2),WN(NWN),W_REF(NWN)
C modif YS
      double precision    QY(NXS,2), FACT
      double precision    ACS1(NWN),ACS2(NWN),TIN(2),TIN_TAB(2,NSPEC)
      integer        I,K,L,N,N1,N2,NE,NMAX,ICOLS, ICOLS_TAB(NSPEC)
C-----------------------------------------------------------------------
      CALL INPUT
      N1           = 1
      N2           = NMAX
      write(6,1100) NMAX

c
c  Read data from control file

c  modif YS
      open(2,file='wn10.in_CB05')   
        read(2,*)
        read(2,*)
        read(2,*)
        read(2,*)
        read(2,*)
        read(2,*)
      do L=1,NSPEC
        read(2,1500) FNAMIN_TAB(L), FNAMOU_TAB(L), LABEL_TAB(L), 
     $   (TIN_TAB(i,L),i=1,2), ICOLS_TAB(L)
      enddo 
      close(2)
      write(*,*) FNAMIN_TAB

      do L=1,NSPEC

      FNAMIN=FNAMIN_TAB(L)
      FNAMOU=FNAMOU_TAB(L)
      LABEL=LABEL_TAB(L)
      TIN(1)=TIN_TAB(1,L)
      TIN(2)=TIN_TAB(2,L)
      ICOLS=ICOLS_TAB(L)
      write(6,1000) FNAMIN, FNAMOU
c
c  Read supplied cross-sections
      open (9,file=FNAMIN)
      read(9,*) FACT
      write(*,*) FACT
      do K = 1,NXS
C modif YS
        read(9,*,end=99) WLEN(K,1),(SEGM(K,I),I=1,ICOLS-1),
     $   (QY(K,I),I=1,ICOLS-1)
        WLEN(K,1) = 10.d0 * WLEN(K,1)
        do I=1,ICOLS-1
           SEGM(K,I)=SEGM(K,I)*QY(K,I)*FACT
        enddo
      enddo
 99   NE = K-1
      write(6,1200) NE
      close(9)
c
c  Interpolate
      call INTPLA(N1,N2,NE,W_REF,ACS1,WLEN(1,1),SEGM(1,1))
      if(ICOLS.le.2.and.TIN(1).eq.TIN(2)) then
        do K = N1,N2
          ACS2(K)  = ACS1(K)
        enddo
      else
        call INTPLA(N1,N2,NE,W_REF,ACS2,WLEN(1,1),SEGM(1,2))
      endif
c
c  Output to new file
      NMAX = N2 - N1 + 1
      open (9,file=FNAMOU)
      write(9,202) LABEL,NMAX,(TIN(N),N=1,2)
      do N = N2,N1,-1
        write(9,'(I5,1P2E12.5)') int(WN(N)),ACS1(N),ACS2(N)
      enddo

      close(9)

      enddo
      
  202 format(A20,I10,2F7.1,'  /total #/T1/T2'/'WN      ACS(T1)',
     &       3X,'ACS(T2)')
c
 1000 format(' Reading file ',a,' and writing to ',a)
 1100 format(' Number of wavenumber bins: ',i5)
 1200 format(' Number of x-sections input:',i5)
 1500 format(a25,a34,a20,2(f5.1,2x),i1)
 2000 format(' Need to alter format')
      STOP
      END
c
c

      SUBROUTINE INPUT
C-----------------------------------------------------------------------
C
C  Routine to input standard wavelengths. Reads from the solar flux
C  file, but the fluxes themselves aren't used at this stage.
C
C-----------------------------------------------------------------------
      implicit none
      integer        NWN
      parameter      (NWN=4600)
      integer        NMAX, N, ICNST
      double precision           WN(NWN),W_REF(NWN), FNMAX
      COMMON /ORI/   WN, W_REF, NMAX
C
      OPEN (1,FILE='avcs_in/sf.d')
      READ (1,501) FNMAX
      NMAX    = int(FNMAX)
      DO N = NMAX, 1,-1
        READ (1,502) ICNST
        WN(N) = real(ICNST)
      ENDDO
      CLOSE (1)
C
      DO N = 1, NMAX
        W_REF(N)  = 1.E+8 /WN(N)
      ENDDO
C
      RETURN
 501  FORMAT(8F10.1)
 502  FORMAT(I6,f10.1)
      END
c
c
      SUBROUTINE INTPLA(N1,N2,JAR,X,Y,XAR,YAR) 
C-----------------------------------------------------------------------
C
C  Power interpolation of the cross-sections onto the 10-wavenumber
C  grid; linear interpolation is used at either end of the data.
C  Replace with your favourite interpolation method, if desired.
C
C-----------------------------------------------------------------------
      implicit none
      integer     N1,N2,JAR
      double precision        X(*),Y(*),XAR(*),YAR(*) 
      integer     N,NB,NE,J,NZ
C
      DO N = N1, N2
        NZ         = N1 + N2 - N
        if (X(N) .lt. XAR(2))  then
          if (X(N) .lt. XAR(1))  then
            Y(N)   = 0.d0
           else
             Y(N)  = ((X(N) - XAR(1))*YAR(2) - (X(N) - XAR(2))*YAR(1))
     &               / (XAR(2) - XAR(1))
             NB    = N
          endif
        endif
        if (X(NZ) .gt. XAR(JAR-1))  then
          if (X(NZ) .gt. XAR(JAR))  then
            Y(NZ)  = 0.d0
           else
             Y(NZ) = ((X(NZ) - XAR(JAR-1))*YAR(JAR) -
     &        (X(NZ) - XAR(JAR))*YAR(JAR-1)) / (XAR(JAR) - XAR(JAR-1))
             NE    = NZ
          endif
        endif
      ENDDO
      J        = 2
      do N = NB+1,NE-1
        if (X(N) .lt. XAR(J))  then
c         A     = (X(N) - XAR(J-1)) / (XAR(J) - XAR(J-1))
c         Y(N)  = YAR(J-1)**(1.d0 - A) * YAR(J)**A
          Y(N)  = ((X(N) -XAR(J-1)) *YAR(J) -(X(N) -XAR(J)) *YAR(J-1))
     &             /(XAR(J) -XAR(J-1))
         else
c          A    = (X(N) - XAR(J)) / (XAR(J+1) - XAR(J))
c          Y(N) = YAR(J)**(1.d0 - A) * YAR(J+1)**A
           Y(N) = ((X(N) -XAR(J)) *YAR(J+1) -(X(N) -XAR(J+1)) *YAR(J))
     &             /(XAR(J+1) -XAR(J))
           J    = J + 1
         endif
       enddo
c
      RETURN
      END
