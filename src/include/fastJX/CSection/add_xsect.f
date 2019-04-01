      program ADD_XSECT
C
C     ******************************************************************
C     *                                                                *
C     *   GENERATING SOLAR FLUX AVERAGED CROSS SECTION IN GIVEN BINS   *
C     *                                                                *
C     ******************************************************************
C
      implicit none
      integer, parameter ::  NIB=100, NWN=4600
      real*8   WKSP(NWN),SF(NWN),ACS1(NWN),ACS2(NWN),WTBIN(NWN)
     +        ,BIN(NIB),FBAR(NIB),WTBINB(2,NIB),WTBINE(2,NIB)
     +        ,TEM(2)
      integer  NTRA(NIB),NMAX,NBIN,I,K,N,M,N1,N2,NE
      character*20 XNAME
      character*15 FNAME
      integer, parameter ::  NSMX=40
      character*20  FNAMIN(NSMX), DIRNAM, DIRNAM_OUT
      integer  NSPEC, IDIR, IDIR_OUT, NJVAL , NP
c
C-------- input X-section ----------------------------------------------
c  program ADD_XSECT.f should be compiled and run using std i/o as below
c  program needs the data files sflx_10cm.dat (solar flux at 10 wave numbers)
c                               wbins_77.dat  (edges of the 77 wavelength bins)

c      open(5,file='NO2-IUPAC.in')
c      open(7,file='NO2-IUPAC.out')

C-------- list of species ----------------------------------------------
      OPEN (1, FILE='list_spec_10wnto77.d_cb05')

      DO N = 1,4
        READ (1,*)
      ENDDO
      READ (1,*) DIRNAM
      READ (1,*) DIRNAM_OUT
      READ (1,*)
      READ (1,*)
      DO N = 1, NSMX
        READ (1,1000,err=99,end=99) FNAMIN(N)
      ENDDO
 99   NSPEC=N-1
      NJVAL=NSPEC+4
      IDIR = 20

      DO N = 20,1,-1
        IF(DIRNAM(N:N).NE.' ') THEN
           IDIR = N
           GO TO 30
        ENDIF
      ENDDO
30    IDIR_OUT = 20
      DO N = 20,1,-1
        IF(DIRNAM_OUT(N:N).NE.' ') THEN
           IDIR_OUT = N
           GO TO 40
        ENDIF
      ENDDO
 40   CLOSE(1)

       write(6,1100) nspec, dirnam(1:idir)

      do  NP = 1,NSPEC
         open(7,file=DIRNAM_OUT(1:IDIR_OUT)//'/CB05/'//FNAMIN(NP))
         open (5,file=DIRNAM(1:IDIR)//'/CB05/'//FNAMIN(NP))

         read (5,'(A20,I10,2F7.1)') XNAME,NE,TEM
         write(7,'(A20,A,2F7.1)') XNAME,' at temp. T1/T2 ',TEM

         read (5,*)
         do N = NE,1,-1
            read (5,'(5X,2F12.1)')  ACS1(N),ACS2(N)
         enddo

         close(5)
c
         call INPUT(SF,BIN,WTBIN,WTBINB,WTBINE,NTRA,NMAX,NBIN)
c
         do K = 1, NBIN
            N1      = NTRA(K) + 1
            N2      = NTRA(K+1) - 2
            FBAR(K) = SF(N1-2)*WTBINB(1,K) + SF(N1-1)*WTBINB(2,K)
     &           +SF(N2+1)*WTBINE(1,K) + SF(N2+2)*WTBINE(2,K)
            do N = N1, N2
               FBAR(K) = FBAR(K) + SF(N)*WTBIN(N)
            enddo
         enddo
c
         call WEIGHT(ACS1,SF,FBAR,WKSP,WTBIN,WTBINB,WTBINE,NTRA,NBIN)
         call WEIGHT(ACS2,SF,FBAR,WKSP,WTBIN,WTBINB,WTBINE,NTRA,NBIN)
c
         write(7,*)
     +   'Bin  Wavelength (nm)   Solar Flux    X-sect(T1) X-sect(T2)'
         do N = 1,NBIN
            write(7,103) N,.1d0*BIN(N),.1d0*BIN(N+1),
     +       FBAR(N),ACS1(N),ACS2(N)
         enddo
c
         WRITE (7,'(/10HSOLAR FLUX)')
         WRITE (7,102)  (FBAR(N), N=1,NBIN)
         WRITE (7,'(A20,F5.0)') XNAME,TEM(1)
         WRITE (7,102)  (ACS1(N), N=1, NBIN)
         WRITE (7,'(A20,F5.0)') XNAME,TEM(2)
         WRITE (7,102)  (ACS2(N), N=1, NBIN)
         close(7)
C
      enddo
      stop
  102 format(1P,8E10.3)
  103 format(I3,2X,F7.2,2H -,F7.2,1P,E12.3,3X,2E11.3)
 1000 format(a20,a20)
 1100 format(2x,i3,' x-section files read from directory ''',a,'''')
      end
C
C     ******************************************************************
C
      subroutine INPUT(SF,BIN,WTBIN,WTBINB,WTBINE,NTRA,NMAX,NBIN)
C
C     ******************************************************************
C     *                                                                *
C     *   DATA INPUT                                                   *
C     *                                                                *
C     ******************************************************************
C
C-----------------------------------------------------------------------
      implicit none
      integer, parameter ::   NIB=100, NWN=4600
      real*8, intent(out) ::  SF(NWN),WTBIN(NWN),BIN(NIB)
     +                       ,WTBINB(2,NIB),WTBINE(2,NIB)
      integer, intent(out) :: NTRA(NIB),NMAX,NBIN
C
      real*8   WN(NWN),W_REF(NWN),CON,FNMAX,FNBIN,WNMID,DWEDG
      integer  I,K,N,M,N1,N2
c
C-------- solar flux in 10 wave number intervals -----------------------
  501 format(8F10.1)
      CON    = 6.6256d0 *2.9979d0
      open (1,file='sflx_10cm.dat')
      read (1,501) FNMAX
      NMAX    = FNMAX
      do N = NMAX, 1,-1
        read (1,'(F6.0,F10.1)') WN(N), SF(N)
        SF(N) = SF(N) *1.d+20 /CON /WN(N)
      enddo
      close (1)
C---------------------------------------- bins -------------------------
      open (1,file='wbins_77.dat')
      read (1,501) FNBIN
      NBIN    = FNBIN
      read (1,501) (BIN(N), N=1,NBIN+1)
      close (1)
C
      do N = 1, NMAX
        W_REF(N)  = 1.E+8 /WN(N)
      enddo
      K           = 0
      do 10 N = 1, NBIN+1
   20   K           = K +1
        if (W_REF(K) .ge. BIN(N))  THEN
          NTRA(N)   = K
          GO TO 10
         else
           GO TO 20
        endif
   10 continue
      do N = 2,NMAX-1
        WTBIN(N)  = 0.5d0*(W_REF(N+1)-W_REF(N-1))
      enddo
      WTBIN(1)    = 0.5d0*(W_REF(2)-W_REF(1))
      WTBIN(NMAX) = 0.5d0*(W_REF(NMAX)-W_REF(NMAX-1))
      do K = 1,NBIN
        N2    = NTRA(K)
        N1    = max(1,N2-1)
        WNMID = 0.5d0*(W_REF(N1)+W_REF(N2))
        DWEDG = BIN(K) - WNMID
        if (DWEDG .lt. 0.d0)  then
          WTBINB(1,K) = -DWEDG
          WTBINB(2,K) = WTBIN(N2)
        else
          WTBINB(1,K) = 0.d0
          WTBINB(2,K) = WTBIN(N2) - DWEDG
        endif
        N2    = NTRA(K+1)
        N1    = N2 - 1
        WNMID = 0.5d0*(W_REF(N1)+W_REF(N2))
        DWEDG = WNMID - BIN(K+1)
        if (DWEDG .lt. 0.d0)  then
          WTBINE(1,K) = WTBIN(N1)
          WTBINE(2,K) = -DWEDG
        else
          WTBINE(1,K) = WTBIN(N1) - DWEDG
          WTBINE(2,K) = 0.d0
        endif
      enddo
c
      return
      end
C
C     ******************************************************************
C
      subroutine WEIGHT(ACS,SF,FBAR,WKSP,WTBIN,WTBINB,WTBINE,NTRA,NBIN)
C
C     ******************************************************************
C     *                                                                *
C     *   CALCULATE SOLAR FLUX WEIGHTED AVERAGE VALUE                  *
C     *                                                                *
C     ******************************************************************
c
      implicit none
      integer, parameter ::   NIB=100, NWN=4600
      real*8, intent(in) ::   WTBIN(NWN),WTBINB(2,NIB),WTBINE(2,NIB)
     +                       ,SF(NWN),FBAR(NIB)
      integer, intent(in) ::  NTRA(NIB),NBIN
      real*8, intent(inout) ::  ACS(NWN),WKSP(NWN)
c
      integer  K,N,N1,N2
C-----------------------------------------------------------------------
      do K = 1, NBIN
        N1      = NTRA(K) + 1
        N2      = NTRA(K+1) - 2
        WKSP(K) = ACS(N1-2)*SF(N1-2)*WTBINB(1,K)
     &          + ACS(N1-1)*SF(N1-1)*WTBINB(2,K)
     &          + ACS(N2+1)*SF(N2+1)*WTBINE(1,K)
     &          + ACS(N2+2)*SF(N2+2)*WTBINE(2,K)
        do N = N1, N2
          WKSP(K) = WKSP(K) + ACS(N)*SF(N)*WTBIN(N)
        enddo
      enddo
      do K = 1, NBIN
        ACS(K)     = WKSP(K)/FBAR(K)
      enddo
c
      return
      end
