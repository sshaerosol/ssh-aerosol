      subroutine ssh_entier(i,ch,ich)
      implicit double precision (a-h,o-z)
c
c     Read an integer in a string

      character *500 ch
      character *1 chtemp
      character *20 chtemp2
      write(chtemp,'(i1)')ich
      chtemp2='(i'//chtemp//')'
      read(ch,chtemp2)i
      return
      end
c

c
      subroutine ssh_part(chdon,mot,imot,nmot)
c
      implicit double precision (a-h,o-z)
c     Decomposition of a sequence (CHDON) in NMOT words (MOT)
c
      parameter (nbmot=100)
      common/nblanc/nblanc
      character *500 chdon
      character *500 mot(nbmot)
      dimension imot(nbmot)
c

      do 10 i0=1,nblanc
         mot(i0)='                                                    '
 10   continue
      nmot=0
      iblanc=1
      ilettre=0
      do 100 nl=1,500
         if(chdon(nl:nl).ne.' ')then
            if(iblanc.eq.1)then
               iblanc=0
               nmot=nmot+1
               ilettre=1
               mot(nmot)(ilettre:ilettre)=chdon(nl:nl)
            else
               ilettre=ilettre+1
               mot(nmot)(ilettre:ilettre)=chdon(nl:nl)
            endif
         else
            if(iblanc.eq.0)imot(nmot)=ilettre
            iblanc=1
         endif
 100  continue
c
      return
      end
c
c
c     *****************************
      subroutine ssh_reel(r,ch,ich)
      implicit double precision (a-h,o-z)
c
c     Read a real in a string

      character *500 ch
      character *2 chtemp
      character *20 chtemp2
      write(chtemp,'(i2)')ich
      chtemp2='(e'//chtemp//'.0)'
      read(ch,chtemp2)r
      return
      end





