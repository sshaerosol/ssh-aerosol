      subroutine ssh_lectci(ifdin)
C------------------------------------------------------------------------
C
C     -- DESCRIPTION
C
C     Read chemical species.
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

c     CHDON:  string variable corresponding to one line of the file.
C     MOT: array composed of the words in CHDON.
C     IMOT: array of sizes for words of MOT.
c     NMOT: number of words.

      implicit double precision (a-h,o-z)
      parameter (nbmot=100)
      include 'parametre'
      include 'ficcom'
c
      common/nblanc/nblanc
c
      character *500 chdon
      character *500 mot(nbmot)
      dimension imot(nbmot)
c

      imp=6
c
c     Loop for reading the input file.
c
      read(ifdin,*)
c
      do ie=1,nesp(2)
         read(ifdin,'(a)')chdon
         nblanc=nbmot
         call ssh_part(chdon,mot,imot,nmot)
         indaq(ie)=0
         nom(ie)(1:imot(1))=mot(1)(1:imot(1))
         inom(ie)=imot(1)
         do k=1,ie-1
            if (nom(k)(1:inom(k)).eq.nom(ie)(1:inom(ie))) then
               write(*,*)'ERROR: species already initialized: ',nom(ie)
               stop 1
            endif
         enddo
      enddo
c
c      read(ifdin,*)
c YK
c
      do ie=nesp(2)+1,nesp(2)+nesp(3)
         read(ifdin,'(a)')chdon
c     write(imp,'(a)')chdon
         nblanc=nbmot
         call ssh_part(chdon,mot,imot,nmot)
         indaq(ie)=1
         nom(ie)(1:imot(1))=mot(1)(1:imot(1))
         inom(ie)=imot(1)
         do k=1,ie-1
            if (nom(k)(1:inom(k)).eq.nom(ie)(1:inom(ie))) then
               write(*,*)'ERROR: species already initialized: ',nom(ie)
               stop 1
            endif
         enddo
      enddo

      return
      end




































