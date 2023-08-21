      subroutine ssh_initconst
C------------------------------------------------------------------------
C
C     -- DESCRIPTION
C
C     Initialization of constants.
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
      include 'ficcom'
c
      av = 6.022d23
c
      return
      end
c
C------------------------------------------------------------------------
      subroutine ssh_lectdata(y0,neq,indicaq)
C------------------------------------------------------------------------
C
C     -- DESCRIPTION
C
C     Read data.
C------------------------------------------------------------------------
C
C     -- INPUT VARIABLES
C
C     -- INPUT/OUTPUT VARIABLES
C
C     -- OUTPUT VARIABLES
C
C     SLUMP: lumped stoichiometric matrix.
C     XLPHY: physical lumping.
C     INDPUR: ii=indpur(i,j) true label of J-th  species in lumping I.
C     Y0: initial conditions.
C     S: stoichiometrix matrix.
C     NALG: physical algebraic onstraints.
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
      include 'ficcom'
      include 'nficfort'
      character *20 fd
c
      dimension y0(nespmax)
c
      ipiste=15
      nmaster=ipiste
      ipiste=ipiste+1
C     Initialization of constants and main arrays.
      call ssh_initconst
      call ssh_initcinet

      nalg=0
      iunit = 0

c
      call ssh_openfic(ipiste,ifdth,filemeca,0)
      call ssh_openfic(ipiste,ifdin,filespecies,0)
c
c     Modif BS multiphase: units molec. cm-3 (ippb=0) ou ppb (1)

c
      read(Ifdin,*)
      read(Ifdin,*)
C     read(Ifdin,*)ippb
C     read(Ifdin,*)
C     End modif BS

      ippb = 0
      read(Ifdin,*)nesp(2),nesp(3)
c
      nesp(1)=nesp(2)+nesp(3)
      indicaq=0
      if (nesp(3).gt.0) indicaq=1
      indicaqcom=indicaq
C!genoa      write(6,*)'Number of multiphase species (nesp): ',nesp(1)
C!genoa      write(6,*)'Number of gas-phase species: ',nesp(2)
C!genoa      write(6,*)'Number of aqueous-phase species: ',nesp(3)
C!genoa      write(6,*)'Max number of species (nespmax): ',nespmax
C!genoa      write(6,*)'Max number of reactions (nrmax): ',nrmax
      if (nesp(1).gt.nespmax) then
         write(*,*)'ERROR: dimension, nesp>nespmax'
         call ssh_halte
      endif

      call ssh_lectci(ifdin)
c
C     ..v.7..x....v....x....v....x....v....x....v....x....v....x....v....x.I
C     PP 12 02 2002
C     Check species name
c      do ie = 1,nesp(1)
c         write(6,777)ie,NOM(ie)
c      enddo
c 777  format(2x,i4,2x,a10)
C     ..v.7..x....v....x....v....x....v....x....v....x....v....x....v....x.I
c     Read chemical mechanism

c      write(*,*)
c      write(*,*)'-----------Chemical mechanism---------'
c      write(*,*)
      call ssh_lectcinet(ifdth,indicaq)
C     Modif BS for gas-phase only
C     call ssh_convcinet(iunit,indicaq)

C
C     Modif BS for gas-phase only
c     call inition
c     call initlphy
c     call initphot
c     End modif gas-phase only.

c     Dimension
c
      neq=ndiff(1)*nbrem
c
 100  format(A10)
c
      return
      end

      subroutine ssh_openfic(ipiste,ifd,fd,inew)
c
      character *256 fd
      IFD =ipiste
      ipiste=ipiste+1
      if (inew.eq.0) open(IFD,file = FD,status = 'old')
      if (inew.eq.1) open(IFD,file = FD,status = 'new')
      return
      end

      subroutine ssh_convcinet(iunit,indicaq)
C------------------------------------------------------------------------
C
C     -- DESCRIPTION
C
C     Conversion of units: define conversion factors.
C     Change units for concentrations (kinetic rates): debug(1) ((2),(3))
C     Initial conditions in molec.cm-3
C
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
c
c
      implicit double precision (a-h,o-z)
      include 'parametre'
      include 'ficcom'

c     Conversion factors ppm <-> molec/cm3
c
      debug(1) = 2.46d10
      debug(2) = 2.46d13
      debug(3) = 6.05d26
c
      if (iunitgas.eq.-1) then
         write(*,*)'ERROR:  gas units undefined.'
         stop 1
      endif
c
      if ((indicaq.eq.1).and.(iunitaq.eq.-1)) then
         write(*,*)'ERROR: liquid units undefined'
         stop 1
      endif
c
      if (iunitgas.ne.iunit) then
         if ((iunitgas.eq.0).and.(iunit.eq.1)) then
            debug2=1.D0/debug(1)
            debug3=1.D0/debug(1)**2
         elseif ((iunitgas.eq.1).and.(iunit.eq.0)) then
            debug2=debug(1)
            debug3=debug(1)**2
         else
            write(*,*)'ERROR: gas units conversion'
            stop 1
         endif
      endif
c
      if ((iunit.eq.1).and.(iunitaq.ne.-1)) then
         write(*,*)'ERROR: liquid units conversion'
         stop 1
      endif
c
      return
      end
C------------------------------------------------------------------------
      subroutine ssh_convci(ippb,iunit,nn,y0)
C------------------------------------------------------------------------
C
C     -- DESCRIPTION
C
C     Conversion of units for initial conditions.
C     Change units for concentrations (kinetic rates): debug(1) ((2),(3))
C
C------------------------------------------------------------------------
C
C     -- INPUT VARIABLES
C
C     IPPB: 1 (0) if initial conditions in ppb (mplec/cm3)
C     IUNIT: 1(2) if computation in molec/cm3 (ppm)
C     NN: number of multiphase species
C
C     -- INPUT/OUTPUT VARIABLES
C
C     Y0: initial conditions
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
      include 'ficcom'
c
      dimension y0(nespmax)
c
c     Conversion factor for concentrations: debug(1)
c     Conversion factor for kinetic rates: debug(2) and debug(3)
c     Initial condtions in molec/cm3
c
c     Conversion ppm <-> molec/cm3
c
      debug(1) = 2.46d10
      debug(2) = 2.46d13
      debug(3) = 6.05d26
c
      if (iunit.ne.ippb) then
         write(*,*)'Initial concentrations without unit conversion'
         if (iunit.eq.0) then
            if (ippb.eq.1) then
               write(*,*)
     $              'concentrations initiales en ppb converties mol/cm3'
c
               do i=1,nn
                  y0(i) = y0(i) * debug(1)
               enddo
            else
               write(*,*)'ERROR: initial concentrations w/o conversion'
               stop 1
            endif
         elseif (iunit.eq.1) then
            if (ippb.eq.0) then
               write(*,*)'Initial concentrations in mol/cm3 -> ppb'
               do i=1,nn
                  y0(i) = y0(i) / debug(1)
               enddo
            else
               write(*,*)
     $              'ERROR: no conversion for initial concentrations'
               stop 1
            endif
         else
            write(*,*)'ERROR: no conversion for initial concentrations'
            stop 1
         endif
      endif
c
      return
      end
c--------------------------------------
      subroutine ssh_halte
      implicit double precision (a-h,o-z)
      stop 1
      return
      end
c--------------------------------------
C     subroutine ssh_precalcul
ccccccccccccccccccccccCCCcccccccccccccc
c     routine de precalcul des parametres
c     a partir des donnees d'entree
ccccccccccccccccccccccccccccccccccccccc
C     implicit double precision (a-h,o-z)
C     include 'parametre'
C     include 'ficcom'
C     common/comprec/bpsave(4,nrmax)
c
ccccc
c     preprocess des ctes cinetiques
ccccc
C     do i=1,nr
C     if (iprecalc(i).eq.2) then
C     bp(1,i)=bpsave(1,i)*dexp(-bp(2,i)/bp(3,i))
C     elseif (iprecalc(i).eq.3) then
C     bp(1,i)=bpsave(1,i)*dexp(-bp(3,i)/bp(4,i))
C     endif
C     enddo
c
C     return
C     end
