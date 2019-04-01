      program FJX1_bins
C
c
c********** this shortened version reads in the single pair of Xsects
c********** from the output of add_xsect.f - just use XXX.out as stdin
c********** and the fast-JX tables are printed to stdout.
c
C   ******************************************************************
C   *                                                                *
C   *   GENERATING SOLAR FLUX AVERAGED CROSS SECTION IN GIVEN BINS   *
C   *   >>>>>>>from ref code 77-bin to 18-bin fast-JX v5.3 <<<<<<<   *
C   *                                                                *
C   ******************************************************************
C
C-----------------------------------------------------------------------
C
C   Combines the 77-bin std Xsects (including XX SR ODF bins) = 145 total
C      into the fast-J(2) bins (7 at trop wavels, 18 for strat+trop)
C      according to a supplied "FJX_bins.dat" file for the Fast-JX code
C
C   history          Standard model:        Xin Zhu
C                    Fast-J & Distribution: Oliver Wild   (9/98)
C                    Fast-J2:               Huisheng Bian (8/02)
C                    new Fast-J(2)          Prather (7/04)
C                    new Fast-JX v-5.3      Prather (6/05)
C
c      open ( 2, file='FJX_bins.dat',status='OLD')

C-----------------------------------------------------------------------
      include 'cmn_bins.f'

      integer, parameter ::  NSMX=40
      character*20  FNAMIN(NSMX), DIRNAM, DIRNAM_OUT
      integer IDIR, IDIR_OUT, NP

C---------77-bin spectral data----------UNIT=10----------------------------

      do J=1,40
       do K=1,3
        TQQ(K,J) = 0.0
       enddo
      enddo

      open(1,file='list_spec_77to18.d_cb05')

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
 30   IDIR_OUT = 20
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

c      open(7,file='no2_fjx.out')
c      open(5,file='NO2-IUPAC.out')

         read( 5,'(A)') TITLE0
         write(7,'(a)') TITLE0
         read( 5,*)

         NWWW = 77
         NW1 = 1
         NW2 = NWWW
         NQQQ = 1

         do IW=1,NWWW
            read(5,'(5x,f7.2,2x,f7.2,e12.3)') WBIN(IW),WBIN(IW+1),FL(IW)
         enddo

         do IW=1,12
            read( 5,*)
         enddo

        do IW=1,NWWW
          WL(IW) = 0.5*(WBIN(IW)+WBIN(IW+1))
        enddo

C-------rest of X-sections
        do JQ=1,NQQQ
           read( 5,'(a20,f5.0)') TITLEJ(1,JQ),TQQ(1,JQ)
           read( 5,'(8e10.3)') (QQQ(IW,1,JQ),IW=1,NWWW)
           read( 5,'(a20,f5.0)') TITLEJ(2,JQ),TQQ(2,JQ)
           read( 5,'(8e10.3)') (QQQ(IW,2,JQ),IW=1,NWWW)
        enddo

        NSPEC = NQQQ
        NBIN = NWWW

        close(5)

C****************now combine 77-bins+ODFs (145 bins) to get fast-J*********

        call COMB1(18)

        close(7)

      enddo
 
      stop     
 1000 format(a20,a20)
 1100 format(2x,i3,' x-section files read from directory ''',a,'''')
       end


      subroutine COMB1(NN)
C-----------------------------------------------------------------------
C  Combine x-sections from 'NBIN' to 'NN' bins  (77 to 18 for Fast-J2)
C-----------------------------------------------------------------------
      include 'cmn_bins.f'
      integer   NN,KT,IG,JG,NG(25,25)
      real*8    RAYLAY
c
C-----Read in wavelength bins to combine:
      open (2, file='FJX_bins.dat',status='OLD')
      do IG = 1,NN
        read(2,'(25I4)') (NG(JG,IG),JG=1,25)
      enddo

      read(2,'(30x,4i5)') NWWW,NWSRB,NSR,NODF
      read(2,'(a)') 
      read(2,'(8e10.3)') (WBIN(IW),IW=1,NWWW+1)
        do IW=1,NWWW
          WL(IW) = 0.5*(WBIN(IW)+WBIN(IW+1))
        enddo
      read(2,'(a)') 
      read(2,'(8e10.3)') (FL(IW),IW=1,NWWW)
      read(2,'(a)') 
      do L = 1,NSR
        read(2,'(6F10.1,I3)')  (ODF(I,L),I=1,NODF), ISR(L)
      enddo

      close(2)
c
C----re-distribution of S-R sub-bins into one dimension
      N=0
      do I=1,NSR
        M=N
        do J=1,ISR(I)
          N=M+J
          FFL(N)=FL(I)*ODF(J,I)
          do JQ=1,NQQQ
            QQQQ(N,1,JQ)=QQQ(I,1,JQ)
            QQQQ(N,2,JQ)=QQQ(I,2,JQ)
          enddo
        enddo
      enddo
c
C----distribute the rest of 77 bins into one dimension above
      NT = N + NBIN - NSR
c      WRITE(6,*) ' TOTAL BINS',NT,NSR
      do I=N+1,NT
        FFL(I)=FL(I-N+NSR)
        do JQ=1,NQQQ
          QQQQ(I,1,JQ)=QQQ((I-N+NSR),1,JQ)
          QQQQ(I,2,JQ)=QQQ((I-N+NSR),2,JQ)
        enddo
      enddo
C
      do I=1,NN
        SFL(I)=0.D0
        do JQ=1,NQQQ
          SQQQ(I,1,JQ)=0.D0
          SQQQ(I,2,JQ)=0.D0
        enddo
      enddo
C
C----replace all X-sect_s by X-sect * Flux
      do N=1,NT
        do IG=1,NN
          do JG=1,25
            if(N.EQ.NG(JG,IG)) then
              call SUMM(IG,N)
            endif
          enddo
        enddo
      enddo
C
      do I=1,NN
        do JQ=1,NQQQ
          QQQ(I,1,JQ)=SQQQ(I,1,JQ)/SFL(I)
          QQQ(I,2,JQ)=SQQQ(I,2,JQ)/SFL(I)
        enddo
      enddo

C-----output in Fast-JX 18-wavelength format

c  (JX_spec.dat) UCI fastJX-5.3: JPL02+irHNO4+IUPAC/NO2/VOC+Bltz (J-8.6, 6/05)

      do JQ=1,NQQQ
        do J=1,2
         write(7,1100)  TITLEJ(J,JQ)(1:7), nint(TQQ(J,JQ)),
     &                                    (QQQ(K,J,JQ),K=NW1,NN)
 1100 format(a7,i3,6(1pE10.3)/(10x,6(1pE10.3))/(10x,6(1pE10.3)))
        enddo
      enddo

c
   99 continue
      return
      end

      
      subroutine SUMM(I,N)
      include 'cmn_bins.f'
         SFL(I)=SFL(I)+FFL(N)
        do J=1,NQQQ
          SQQQ(I,1,J)=SQQQ(I,1,J)+QQQQ(N,1,J)*FFL(N)
          SQQQ(I,2,J)=SQQQ(I,2,J)+QQQQ(N,2,J)*FFL(N)
        enddo
      return
      end


