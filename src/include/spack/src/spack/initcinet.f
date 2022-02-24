      subroutine  ssh_initcinet
C------------------------------------------------------------------------
C
C     -- DESCRIPTION
C
C     Initialization.
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
      iunitaq=-1
      iunitgas=-1
      nr=0
      nrphot = 0
      nequil=0
      ntabphot=0

      do i=1,nphase
         nrp(i)=0
         nrmol1(i)=0
         nrmol2(i)=0
         nrmol3(i)=0
         naq2(i)=0
         narr1(i)=0
         narr2(i)=0
         narr3(i)=0
         narr4(i)=0
         narr5(i)=0
         narr6(i)=0
         narr7(i)=0
         narr8(i)=0
         nthird(i)=0
      enddo
c
      do i=1,nespmax
         rmol(i)=0.D0
         ihreac(i)=0
         do j=1,nionx
            seqion(i,j)=0.0d0
         enddo
      enddo
c
      do ir=1,nrmax
         iphotinv(ir)=0
         molec(ir)=0
         bp(1,ir)=0.D0
         bp(2,ir)=0.D0
         bp(3,ir)=0.D0
         bp(4,ir)=0.D0
         bp(5,ir)=0.D0
         bp(6,ir)=0.D0
         ispebp(ir)=0
         ihetero(ir)=0
         jer(1,ir)=0
         jer(2,ir)=0
         jer(3,ir)=0
         ittb(ir)=0
         iprecalc(ir)=0
c
         do i=1,nphase
            irmonodi(ir,i)=0
            imolec1(ir,i)=0
            imolec2(ir,i)=0
            imolec3(ir,i)=0
            iarr1(ir,i)=0
            iarr2(ir,i)=0
            iarr3(ir,i)=0
            iarr4(ir,i)=0
            iarr5(ir,i)=0
            iarr6(ir,i)=0
            iarr7(ir,i)=0
            iarr8(ir,i)=0
            inaq2(ir,i)=0
            indthird(ir,i)=0
         enddo
c
         irmonodi(ir,1)=ir
         iheq(ir)=0
         nb(ir)=0
         nneq(ir)=0
c
         do i=1,nespmax
            s(i,ir) = 0.0d0
         enddo
c
      enddo
c
      do i=1,nespmax
         iemonodi(i,1)=i
         iedimono(i,1)=i
         do j=2,nphase
            iemonodi(i,j)=0
            iedimono(i,j)=0
         enddo
      enddo
c
      do i=1,nphase
         ndiff(i)=0
      enddo
c
      nequil11=0
      nequil12=0
      nequil13=0
      nequil21=0
      nequil22=0
      nequil23=0
      nequil3=0
      nequil41=0
      nequil42=0
      nequil51=0
      nequil52=0
c
      jhplus = 0
      johmoin = 0

      do i=1,nequilx
c
         jhpoh(i)=0
         idifford(i)=0
         jion1(i) = 0
         jion2(i) = 0
         jaq(i) = 0
c
      enddo

      RETURN
      END




