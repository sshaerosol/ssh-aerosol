      subroutine lectcinet(ifdth,indicaq)
C------------------------------------------------------------------------
C
C     -- DESCRIPTION
C
C     Read chemical mechanism.
C------------------------------------------------------------------------
C
C     -- INPUT VARIABLES
C
C     IFDTH: input file
C     FPAR: physical parameters file
C
C     -- INPUT/OUTPUT VARIABLES
C
C     -- OUTPUT VARIABLES
C
C     NESP: number of species.
C     NR: number of reactions.
C     XLPHY: physical lumping.
C     INDPUR: ii=indpur(i,j) true label of j-th species in i-th lumping.
C     NALG: number of algebraic constraints.
C
C------------------------------------------------------------------------
C
C     -- REMARKS
C
C------------------------------------------------------------------------
C
C     -- MODIFICATIONS
C
C------------------------------------------------------------------------
C
C     -- AUTHOR(S)
C
C     Bruno Sportisse, CEREA, 2003.
C
C------------------------------------------------------------------------

      implicit double precision (a-h,o-z)
      include 'parametre'
      parameter (nbmot=100)
c
      include 'ficcom'
      common/nblanc/nblanc
c
      character *500 chdon
      character *500 mot(nbmot)
      dimension imot(nbmot)
c
      imp=6
      iattente=0
      iarret=0
c
c     Loop for reading the file.
c     Write the first lines of the output routines

      call write_header
C
 100  read(ifdth,'(a)')chdon
c      write(imp,'(a)')chdon
c     YK
c
c     Build the sequence and decompose in words.

c
      nblanc=nbmot
      call part(chdon,mot,imot,nmot)
c
c     Case of long sequences (two lines).
c
 102  continue
      if(mot(nmot)(1:2).eq.'//')then
         nchold=nmot
         read(ifdth,'(a)')chdon
c         write(imp,'(a)')chdon
c     YK
         nblanc=nbmot-nchold+1
c
         call part(chdon,mot(nchold),imot(nchold),nchplus)
         nmot=nchold-1+nchplus
      endif
c
      if (mot(nmot)(1:2).eq.'//') go to 102
c
c     Call the different routines according to the keyword.
c
      if(mot(1)(1:3).eq.'KIN')  then
         if (iattente.eq.0) then
            write(*,*)'ERROR: kinetics before reaction ',nr
            stop 1
         elseif (iattente.eq.1) then
            call kinreac(mot,imot,nmot)
         elseif (iattente.eq.2) then
            call kindis(mot,imot,nmot)
         elseif (iattente.eq.3) then
            call kinhenry(mot,imot,nmot)
         endif
         iattente=0
c
      elseif (mot(1)(1:3).eq.'SET') then
c
         if (mot(2)(1:4).eq.'UNIT')then
            call lectunit(mot,imot,nmot)
c
         elseif (mot(2)(1:10).eq.'TABULATION') then
            call lect_tabulation(mot,imot,nmot)
c
         else
            write(*,*)'ERROR: UNKNOWN SET FUNCTIONS'
            stop 1
         endif
c     Symbols for commented lines: %, !,
      elseif (mot(1)(1:1).eq.'%')then
      elseif (mot(1)(1:1).eq.'!')then

c     END.
      elseif (mot(1)(1:3).eq.'END')then
         iarret=1
      elseif (iattente.ne.0) then
         write(*,*)'ERROR: I wait for kinetics'
         stop 1
      else
         call reaction(mot,imot,nmot,iattente)
      endif
c
      if (iarret.eq.0) go to 100
c
      CALL write_end
c
c     Write file non_zero.dat
c
      call WNONZERO(s,nr,jer)
C
      write(*,*)'########################################'
      write(*,*)'########################################'
      write(6,*)'Summary for the kinetic scheme'
      write(*,*)'Total number of reactions =',nr
      write(*,*)'Number of photolytic reactions =',nrphot
      write(*,*)'Number of dissociation equilibria =',nequil
c
      if (indicaq.eq.0) then
         write(*,*)'Gas-phase chemistry'
      endif
c
      if (indicaq.eq.1) then
         write(*,*)'Multiphase (gas-phase and aqueous-phase) chemistry'
      endif
c
      call initphase
c
      return
      end
c
c---------------------------------------------------
      subroutine lectunit(mot,imot,nmot)
C------------------------------------------------------------------------
C
C     -- DESCRIPTION
C
C     Read units.
C------------------------------------------------------------------------
C
C     -- INPUT VARIABLES
C
C     -- INPUT/OUTPUT VARIABLES
C
C     -- OUTPUT VARIABLES
C
C------------------------------------------------------------------------
C
C     -- REMARKS
C
C------------------------------------------------------------------------
C
C     -- MODIFICATIONS
C
C------------------------------------------------------------------------
C
C     -- AUTHOR(S)
C
C     Bruno Sportisse, CEREA, 2003.
C
C------------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      include 'parametre'
      parameter (nbmot=100)
      include 'ficcom'
c
      character *500 mot(nbmot)
      dimension imot(nbmot)
c

      iunitgas = 0
      if (mot(2)(1:imot(2)).eq.'GAS') then
         if (mot(3)(1:imot(3)).eq.'MOLCM3') then
            iunitgas=0
         elseif (mot(3)(1:imot(3)).eq.'PPB') then
            iunitgas=1
         else
            write(*,*)'ERROR: unknown units for gas-phase kinetics'
            stop 1
         endif
      elseif (mot(2)(1:imot(2)).eq.'AQ') then
         if (mot(3)(1:imot(3)).eq.'MOLL') then
            iunitaq=2
         else
            write(*,*)'ERROR: unknown units for aqueous-phase kinetics'
            stop 1
         endif
      endif
      return
      end
c---------------------------------------------------
      subroutine lect_tabulation(mot,imot,nmot)
C------------------------------------------------------------------------
C
C     -- DESCRIPTION
C
C     Read tabulation for photolysis.
C     The format is: TABULATION N DEGREES D1 D2 ... DN
C     N is the number of tabulated angles of values Di (in degrees).
C     The sequence may be increasing or decreasing.
C------------------------------------------------------------------------
C
C     -- INPUT VARIABLES
C
C     -- INPUT/OUTPUT VARIABLES
C
C     -- OUTPUT VARIABLES
C
C------------------------------------------------------------------------
C
C     -- REMARKS
C
C------------------------------------------------------------------------
C
C     -- MODIFICATIONS
C
C------------------------------------------------------------------------
C
C     -- AUTHOR(S)
C
C     Bruno Sportisse, CEREA, 2003.
C
C------------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      include 'parametre'
      parameter (nbmot=100)
      include 'ficcom'
c
      character *500 mot(nbmot)
      dimension imot(nbmot)
      dimension t(ntabphotmax)
c
      ireversetab=0
c
      call entier(ntabphot,mot(3),imot(3))
      if (ntabphot.eq.0) then
         write(*,*)'ERROR: number of tabulated angles to be checked in'
         write(*,*)'subroutine lect-tabulation'
         stop 1
      endif
      if (ntabphot.gt.ntabphotmax) then
         write(*,*)'ERROR: ntabphot>ntabphotmax'
         stop 1
      endif

      do i=1,ntabphot
         call reel(t(i),mot(i+4),imot(i+4))
      enddo
c     Check increasing order
      if (t(1).gt.t(2)) then
         ireversetab=1
         do j=1,ntabphot
            tabphot(j)=t(ntabphot+1-j)
         enddo
      else
         do j=1,ntabphot
            tabphot(j)=t(j)
         enddo
      endif
      do j=1,ntabphot-1
         if (tabphot(j).ge.tabphot(j+1)) then
            write(*,*)
     &           'ERROR: the tabulation has to be strictly monotonic'
            stop 1
         endif
      enddo
c
      return
      end
c---------------------------------------------------
      subroutine reaction(mot,imot,nmot,iattente)
C
C     -- DESCRIPTION
C
C     Read reactions.
C------------------------------------------------------------------------
C
C     -- INPUT VARIABLES
C
C     -- INPUT/OUTPUT VARIABLES
C
C     -- OUTPUT VARIABLES
C
C------------------------------------------------------------------------
C
C     -- REMARKS
C
C------------------------------------------------------------------------
C
C     -- MODIFICATIONS
C
C------------------------------------------------------------------------
C
C     -- AUTHOR(S)
C
C     Bruno Sportisse, CEREA, 2003.
C
C------------------------------------------------------------------------
c     implicit double precision (a-h,o-z)
      include 'parametre'
      parameter (nbmot=100)
c
      character *500 mot(nbmot)
      dimension imot(nbmot)
c
      do i=1,nmot
         if (mot(i)(1:imot(i)).eq.'>') then
            if (iattente.ne.0) then
               write(*,*)'ERROR: too many symbols >'
               stop 1
            else
               iattente=1
            endif
         elseif (mot(i)(1:imot(i)).eq.'->') then
            if (iattente.ne.0) then
               write(*,*)'ERROR: too many symbols ->'
               stop 1
            else
               iattente=1
            endif
         elseif (mot(i)(1:imot(i)).eq.'=') then
            if (iattente.ne.0) then
               write(*,*)'ERROR: too many symbols ='
               stop 1
            else
               iattente=2
            endif
         elseif (mot(i)(1:imot(i)).eq.'<H>') then
            if (iattente.ne.0) then
               write(*,*)'ERROR: too many symbols <H>'
               stop 1
            else
               iattente=3
            endif
         elseif (mot(i)(1:imot(i)).eq.'=H=') then
            if (iattente.ne.0) then
               write(*,*)'ERROR: too many symbols =H='
               stop 1
            else
               iattente=4
            endif
         endif
      enddo
c
      if (iattente.eq.0) then
         write(*,*)'ERROR: lack of symbol for reaction'
         stop 1
      elseif (iattente.eq.1) then
         call creac(mot,imot,nmot)
      elseif (iattente.eq.2) then
         call cdis(mot,imot,nmot)
      elseif (iattente.eq.3) then
         call chenry(mot,imot,nmot,0)
      elseif (iattente.eq.4) then
         call chenry(mot,imot,nmot,1)
         iattente=3
      endif
c
      return
      end














