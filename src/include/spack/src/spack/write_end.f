      subroutine write_end
C------------------------------------------------------------------------
C
C     -- DESCRIPTION
C
C     Write end of subroutines fexchem.f, jacdchemdc.f and kinetic.f.
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
C     Label 300: blank line.
C     Label 150: comments.
C     Label 200: C-------...
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
      include 'nficfort'
      include 'parametre'
      include 'ficcom'

      integer HR1,HR2,HR3,HR4
      common /HR/ HR1,HR2,HR3,HR4

C     Routine dimensions.f
      nwrite=nficdim
      write(nwrite,301) function_suffix
 301  format(6x,'subroutine dimensions',a,'(')
      write(nwrite,3011)
 3011 format(5x,'$',4x,'Ns, Nr, Nr_photolysis)')
      write(nwrite,300)
      write(nwrite,200)
      write(nwrite,150)
      write(nwrite,302)
 302  format(6x,'integer Ns, Nr, Nr_photolysis')
      write(nwrite,300)
      write(nwrite,303) nesp(1)
 303  format(6x,'Ns = ',i5)
      write(nwrite,304) nr
 304  format(6x,'Nr = ',i5)
      write(nwrite,305) nrphot
 305  format(6x,'Nr_photolysis = ',i5)

      call write_end2(nwrite)

C     Routine fexchem.f
      nwrite=nficf

      write(nwrite,300)
      write(nwrite,401)
 401  format('C    Conversion molecules/cm3 to mug/m3.')
      write(nwrite,300)
      write(nwrite,402)
 402  format('      do i = 1, ns')
      write(nwrite,403)
 403  format('         chem(i) = chem(i) / convers_factor(i)')
      write(nwrite,404)
 404  format('      enddo')
      write(nwrite,300)
      write(nwrite,101)
 101  format('C     Volumic source terms.')
      write(nwrite,300)
C =============================
      if (aerosol_formation) then
      write(nwrite,102)
 102  format('      IF(nemis.gt.0)THEN')
      endif
C ===================
      write(nwrite,103)
 103  format('      do i=1,ns')
      write(nwrite,104)
 104  format('      chem(i)=chem(i)+ZCsourc(i)')
      write(nwrite,105)
 105  format('      enddo')
C =============================
      if (aerosol_formation) then
      write(nwrite,106)
 106  format('      endif')
      endif
C =============================
      call write_end2(nwrite)

C     Routine kinetic.f
      nwrite=nfick
C ======================
      if (aerosol_formation) then
      write(nwrite,601)
 601  format('C    Heteroheneous Reactions on aerosol surface:')
      write(nwrite,150)
      write(nwrite,602) HR1
 602  format('C    Reaction ',I3,': HO2  -->  0.5 H2O2')
      write(nwrite,603) HR2
 603  format('C    Reaction ',I3,': NO2  -->  0.5 HONO + 0.5 HNO3')
      write(nwrite,604) HR3
 604  format('C    Reaction ',I3,': NO3  -->  HNO3')
      write(nwrite,605) HR4
 605  format('C    Reaction ',I3,': N2O5 -->  2 HNO3')
      write(nwrite,300)
      write(nwrite,606) HR1
 606  format(6x,'rk',I3,'=0.D0')
      write(nwrite,607) HR2
 607  format(6x,'rk',I3,'=0.D0')
      write(nwrite,608) HR3
 608  format(6x,'rk',I3,'=0.D0')
      write(nwrite,609) HR4
 609  format(6x,'rk',I3,'=0.D0')
      write(nwrite,300)
      write(nwrite,610)
 610  format(6x,'if (IHETER.eq.1) then')
      write(nwrite,611)
 611  format(6x,'   call HETRXN(Ns,Nbin_aer,temp,press,ICLD,lwctmp,')
      write(nwrite,612) HR1,HR2,HR3,HR4
 612  format(5x,'&        WetDiam,granulo,rk',I3,',rk',I3,',rk',I3,
     2      ',rk',I3,',')
      write(nwrite,613)
 613  format(5x,'&        dsf_aero,ispeclost,Wmol,LWCmin)')
      write(nwrite,614)
 614  format(6x,'endif')
      write(nwrite,300)
      write(nwrite,615) HR1,HR1
 615  format(6x,'rk(',I3,') = rk',I3)
      write(nwrite,616) HR2,HR2
 616  format(6x,'rk(',I3,') = rk',I3)
      write(nwrite,617) HR3,HR3
 617  format(6x,'rk(',I3,') = rk',I3)
      write(nwrite,618) HR4,HR4
 618  format(6x,'rk(',I3,') = rk',I3)
      endif
C =====================================


      call write_end2(nwrite)

C     Routine jacdchemdc.f
      nwrite=nficj
      write(nwrite,300)
      write(nwrite,502)
 502  format('      do j = 1, ns')
      write(nwrite,503)
 503  format('         do i = 1, ns')
      write(nwrite,504)
 504  format('            JacC(i,j) = JacC(i,j) *',
     2      ' convers_factor_jac(i,j)')
      write(nwrite,505)
 505  format('         enddo')
      write(nwrite,506)
 506  format('      enddo')
      write(nwrite,300)

      call write_end2(nwrite)

C     Routine fexloss.f
      nwrite=nficloss
      call write_end2(nwrite)

C     Routine fexprod.f
      nwrite=nficprod
      call write_end2(nwrite)

C     Routine rates.f
      nwrite=nficw
      call write_end2(nwrite)

C     Routine dratedc.f
      nwrite=nficdw
      call write_end2(nwrite)



 150  format('C',6x,a65)
 200  format('C-----------------------------------',
     2     '-------------------------------------')
 300  format(' ')
 110  format(6x,'RETURN',/,6x,'END')

      RETURN
      END

C-----------------------------------------------------

      subroutine write_end2(nwrite)
      integer nwrite

      write(nwrite,300)
      write(nwrite,110)
      write(nwrite,300)

 150  format('C',6x,a65)
 200  format('C-----------------------------------',
     2     '-------------------------------------')
 300  format(' ')
 110  format(6x,'RETURN',/,6x,'END')

      RETURN
      END
