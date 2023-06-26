C------------------------------------------------------------------------
C     Subroutines for automatic generation of fortran files
C     needed for integration of gas-phase chemistry.
C     Authors: Bruno Sportisse and Pierre Plion, CEREA/ENPC
C     Date: February 2003.
C------------------------------------------------------------------------
C     genoa add for Complex rate coefficients in MCM v3.3.1
      subroutine ssh_WKMCM1(nr,b1,b2,b3,b4,b5,b6,b7,b8)
C------------------------------------------------------------------------
C
C     -- DESCRIPTION
C
C     http://mcm.leeds.ac.uk/MCMv3.3.1/parameters/complex.htt
C
C     ka = a1 * M * (T / a2) ** a3
C
C     if b2 != 0.0:
C       kb = b1 * (T / b2) ** b3
C     else:
C       kb = b1  ###(KMT07, 09, 13)
C
C     kc = 0.75 - 1.27 * (LOG10(c1))
C     kd = 10 ** (LOG10(c1)/(1+(LOG10(ka/kb)/kc)**2))
C     k  = (ka * kb) * kd/(ka + kb) * b8 (b8 in case of ratio)
C
C     Use for: KFPAN, KMT01, 02, 03, 07, 08, 09, 13, 15, 16
C              (12 and 17 are (ka * kb * kd)/(ka + kb))
C
C     specific case: KMT04, KMT05, KMT06, KMT08, KMT11, KMT12, KMT17, KMT18
C------------------------------------------------------------------------
      include 'nficfort'
      integer nr
      double precision b1,b2,b3,b4,b5,b6,b7,b8
C
      write(nficK90,10)b1,b2,b3 !ka Effko

      if (b5 .ne. 0.0d0) then
        write(nficK90,11)b4,b5,b6 !kb Rapk
      else
        write(nficK90,12)b4 !kb Rapk
      endif

      write(nficK90,13)b7,b7    !kd facteur
      write(nficK90,14)nr,b8    !k

 10   format('  Effko = ',D23.16,'* SumM * &',/
     &     '    (temp / ',D23.16,') &',/
     &     '    **(',D23.16,')')
 11   format('  Rapk = ',D23.16,'* (temp / ',D23.16,') &',/
     &     '    **(',D23.16,')')
 12   format('  Rapk = ',D23.16)
 13   format('  facteur = 1.0d1 ** (LOG10(',D23.16,') &',/
     &     '    /(1.0d0 +(LOG10( Effko / Rapk ) &',/
     &     '    / (0.75d0 - 1.27d0 * &',/
     &     '    (LOG10(',D23.16,'))) )**2.0d0))')
 14   format('  rk(',i4,') = (Effko*Rapk)*facteur &',/
     &     '    /(Effko+Rapk)*',D23.16)
      return
      end
C------------------------------------------------------------------------
      subroutine ssh_WKMCM2(nr,b1,b2,b3,b4,b5,b6,b7,b8)
C------------------------------------------------------------------------
C
C     -- DESCRIPTION
C
C     ka = a1 * M * dexp(a2 / T) ** a3
C     kb = b1 * dexp(b2 / T) ** b3
C     kc = 0.75 - 1.27 * (LOG10(c1))
C     kd = 10 ** (LOG10(c1)/(1+(LOG10(ka/kb)/kc)**2))
C     k  = (ka * kb) * kd/(ka + kb) * b8 (b8 in case of ratio)
C
C     Use for: KBPAN, KBPPN, KMT10, 14
C------------------------------------------------------------------------
      include 'nficfort'
      integer nr
      double precision b1,b2,b3,b4,b5,b6,b7,b8
C
      write(nficK90,10)b1,b2,b3 !ka Effko
      write(nficK90,11)b4,b5,b6 !kb Rapk
      write(nficK90,12)b7,b7    !kd facteur
      write(nficK90,13)nr,b8    !k

 10   format('  Effko = ',D23.16,'* SumM * &',/
     &     '    dexp(',D23.16,'/temp) &',/
     &     '    **(',D23.16,')')
 11   format('  Rapk = ',D23.16,'* dexp(',D23.16,'/temp) &',/
     &     '    **(',D23.16,')')
 12   format('  facteur = 1.0d1 ** (LOG10(',D23.16,') &',/
     &     '    /(1.0d0 +(LOG10( Effko / Rapk ) &',/
     &     '    / (0.75d0 - 1.27d0 * &',/
     &     '    (LOG10(',D23.16,'))) )**2.0d0))')
 13   format('  rk(',i4,') = (Effko*Rapk)*facteur &',/
     &     '    /(Effko+Rapk)*',D23.16)
      return
      end
C------------------------------------------------------------------------
C     genoa add for photolysis in MCM
      subroutine ssh_WKMCM3(nr,b1,b2,b3)
C------------------------------------------------------------------------
C
C     -- DESCRIPTION
C
C     IN THE FORMAT J = l* cosX ** m *dexp(-n*secX)
C
C------------------------------------------------------------------------
      include 'nficfort'
      integer nr
      double precision b1,b2,b3
C
      write(nficK90,10) nr,b1,b2,b3 ! l,m,n
 10   format('  rk(',i4,') = ',D23.16,' * cosX **(',D23.16,') &',/
     &     '* dexp( -',D23.16,' * secX) ')
      return
      end

C------------------------------------------------------------------------
      subroutine ssh_WK190(nr,k)
C------------------------------------------------------------------------
C
C     -- DESCRIPTION
C
C     Write kinetic rates case ARR1:
C     k = constant
C
C------------------------------------------------------------------------
C
C     -- AUTHOR(S)
C
C     Pierre Plion, 2002.
C
C------------------------------------------------------------------------
      implicit none
      double precision k
      integer nr
      include 'nficfort'

      write(nficK90,10)nr,k

 10   format('  rk(',i4,') = ',D23.16)

      return
      end
C------------------------------------------------------------------------
      subroutine ssh_WK290(nr,k,Tact)
C------------------------------------------------------------------------
C
C     -- DESCRIPTION
C
C     Write kinetic rates case ARR2:
C     k(T) = b1 * dexp(-b2/T)
C
C------------------------------------------------------------------------
C
C     -- AUTHOR(S)
C
C     Pierre Plion, 2002.
C
C------------------------------------------------------------------------
      implicit none
      double precision logk,k,Tact
      integer nr
      include 'nficfort'

      if(abs(Tact).lt.petit)then
         call ssh_WK190(nr,k)
      else
         logk = dlog(k)
         write(nficK90,10)nr,logk,Tact
      endif

 10   format('  rk(',i4,') =  dexp( ',D23.16,' &',/
     &     '   - ( ',D23.16,' )/temp)')

      return
      end
C------------------------------------------------------------------------
      subroutine ssh_WK390(nr,k,expT,Tact)
C------------------------------------------------------------------------
C
C     -- DESCRIPTION
C
C     Write kinetic rates case ARR3:
C     k(T) = b1 * (T**b2) * dexp(-b3/T)
C
C------------------------------------------------------------------------
C
C     -- AUTHOR(S)
C
C     Pierre Plion, 2002.
C
C------------------------------------------------------------------------
      implicit none
      double precision logk,k,Tact,expT
      integer nr
      include 'nficfort'

      if(abs(expT).lt.petit) then
         call ssh_WK290(nr,k,Tact)
      else
         logk = dlog(k)
         if(Tact.gt.petit)then
            write(nficK90,10)nr,logk,expT,Tact
         elseif(Tact .lt. -petit) then
            write(nficK90,11)nr,logk,expT,(-Tact)
         else
            write(nficK90,12)nr,logk,expT
         endif
      endif

 10   format('  rk(',i4,') =  dexp(',D23.16,' &',/
     &     '   + (',D23.16,' * dlog(temp)) &',/
     &     '   - ',D23.16,'/temp)')
 11   format('  rk(',i4,') =  dexp(',D23.16,' &',/
     &     '   + (',D23.16,' * dlog(temp)) &',/
     &     '   + ',D23.16,'/temp)')
 12   format('  rk(',i4,') =  dexp(',D23.16,' &',/
     &     '   + (',D23.16,' * dlog(temp)) ) ')

      return
      end
C------------------------------------------------------------------------
      subroutine ssh_WTROE90(nr,b1,b2,b3,b4,b5)
C------------------------------------------------------------------------
C
C     -- DESCRIPTION
C
C     Write kinetic rates case for Troe reactions:
C
C     k(T)=k0*M/(1+r)  * b5**(1/(1+[log10 r]**2))
C     with r= k0*M/kinf, k0=b1*(T/300)**(-b2), kinf=b3*(T/300)**(-b4)
C     b5=0.6 if not specified (case of KINETIC TROE4 ...).
C------------------------------------------------------------------------
C
C     -- AUTHOR(S)
C
C     Pierre Plion and Bruno Sportisse, 2002.
C
C------------------------------------------------------------------------
      include 'nficfort'
      integer nr
      double precision b1,b2,b3,b4,b5

      write(nficK90,10)b1,b2
      write(nficK90,11)b3,b4
      write(nficK90,12)nr,b5

c     Modification BS/KS 06/2002
c     Replace log10 by log and Effko*Rapk by Effko/Rapk

 10   format('  Effko = ',D23.16,'* (temp / 3.d2) &',/
     &     '    **(- (',D23.16,'))')
 11   format('  Rapk = ',D23.16,'* (temp / 3.d2) &',/
     &     '    **(- (',D23.16,'))')
 12   format('  rk(',i4,') = (Effko * SumM / &',/
     &     '    ( 1.0d0 + Effko * SumM / Rapk)) * &',/
     &     '    ',D10.4,'** (1.0d0 / (1.0d0 + &',/
     &     '    (log10(Effko * SumM / Rapk))**2))')
c
      return
      end
C------------------------------------------------------------------------
      subroutine ssh_WTROE790(nr,b1,b2,b3,b4,b5,b6,b7)
C------------------------------------------------------------------------
C
C     -- DESCRIPTION
C
C     Write kinetics for reactions computed from equilibria.
C     k(T)=k0*M/(1+r)  * b7**(1/(1+[log10 r]**2))
C     with r= k0*M/kinf
C     k0=b1*dexp(-b2/T)*(T/300)**(-b3)
C     kinf=b4*dexp(-b5/T)*(T/300)**(-b6)
C
C------------------------------------------------------------------------
C
C     -- AUTHOR(S)
C
C     Pierre Plion and Bruno Sportisse 2002.
C
C------------------------------------------------------------------------
      include 'nficfort'
      integer nr
      double precision b1,b2,b3,b4,b5,b6,b7
C
      write(nficK90,10)b1,b3,b2
      write(nficK90,11)b4,b6,b5
      write(nficK90,12)nr,b7

 10   format('  Effko = ',D23.16,'* (temp / 3.d2) &',/
     &     '    **(- (',D23.16,'))* &',/
     &     '    dexp(-',D23.16,'/temp)')
 11   format('  Rapk = ',D23.16,'* (temp / 3.d2) &',/
     &     '    **(- (',D23.16,'))* &',/
     &     '    dexp(-',D23.16,'/temp)')
 12   format('  rk(',i4,') = (Effko*SumM / ( 1.0d0 + Effko * SumM / &',/
     &     '    Rapk)) * ',D10.4,'** (1.0d0 / (1.0d0 + &',/
     &     '    (log10(Effko * SumM / Rapk))**2))')
C     13    format(6x,'rk(',i4,') = facteur * (',D23.16,') * '/
C     &       5x,'&      dexp((',D23.16,')/temp)')
c
      return
      end
C------------------------------------------------------------------------
      subroutine ssh_WTROE1090(nr,b1,b2,b3,b4,b5,b6,b7,b8,b9,b10)
C------------------------------------------------------------------------
C
C     -- DESCRIPTION
C
C     Write kinetic rates case for Troe reactions/MOCA.
C
C------------------------------------------------------------------------
C
C     -- AUTHOR(S)
C
C     Pierre Plion and Bruno Sportisse, 2002.
C
C------------------------------------------------------------------------
      include 'nficfort'
      integer nr
      double precision b1,b2,b3,b4,b5,b6,b7,b8,b9,b10

C     ! ! ! MOCA computes QLN = dlog(k)
c     Modification BS/KS: 22/03/2002
c     Error for b5,b6 -> b6,b5
c
      if(b6.gt.petit)then
         write(nficK90,10)b4,b6,b5
      elseif(b6.lt.-petit)then
         write(nficK90,20)b4,(-b6),(b5)
      else
         write(nficK90,30)b4,b5
      endif
c
c     End modification.
c
c     Modification BS/KS: 20/03/2002
c     Error b1 --> dlog(b1)
c
      if(b3.gt.petit)then
         write(nficK90,11)dlog(b1),b2,b3
      elseif(b3.lt.-petit)then
         write(nficK90,21)dlog(b1),b2,(-b3)
      else
         write(nficK90,31)dlog(b1),b2
      endif
c
c     End modification.
      write(nficK90,12)
      IF(abs(b8).lt.petit)then
         write(nficK90,14)nr,b7
      else
         write(nficK90,15)(1.d0-b7),b8,b7,b9,b10
         if(abs(b10).lt.petit) write(nficK90,16)
         write(nficK90,17)nr
      endif

 10   format('  Effko = ',D23.16,'* SumM * &',/
     &     '    dexp(-(',D23.16,')/temp + &',/
     &     '    (', D23.16,') * dlog(temp/3.d2))')

 20   format('  Effko = ',D23.16,'*dexp(',D23.16,'/temp + &',/
     &     '    (', D23.16,') * dlog(temp/3.d2)) * SumM')

 30   format('  Effko = ',D23.16,'*dexp( &',/
     &     '    ',D23.16,' * dlog(temp/3.d2)) * SumM')

 11   format('  Rapk = Effko / dexp(',D23.16,' + &',/
     &     '    (',D23.16,')* dlog(temp/3.d2) - &',/
     &     '    (',D23.16,')/temp)')

 21   format('  Rapk = Effko / dexp(',D23.16,' + &',/
     &     '    (',D23.16,')* dlog(temp/3.d2) + &',/
     &     '    (',D23.16,')/temp)')

 31   format('  Rapk = Effko / dexp(',D23.16,' + &',/
     &     '    (',D23.16,')* dlog(temp/3.d2) ) ')

 12   format('  facteur = 1.d0/(1.d0+log10(Rapk)**2)')
C     12    format(6x,'facteur = 0.d0',/
C     &   6x,'if(Rapk.gt.0.d0) facteur = 1.d0/(1.d0+log10(Rapk)**2)')
 14   format('  rk(',i4,') = dexp(dlog(Effko/(1.d0+Rapk)) + &',/
     &     '    dlog(',D23.16,') * facteur)')

 15   format('  Fcent = ',D23.16,' * dexp(-temp/ &',/
     &     '    (', D23.16,')) &',/
     &     '    + (',D23.16,') * dexp(-temp/ &', /
     &     '    (',D23.16,')) &',/
     &     '    +  exp (-(',D23.16,')/temp)')
 16   format('  Fcent = Fcent -1.d0')
 17   format('  rk(',i4,') = dexp(dlog(Effko/(1.d0+Rapk)) &',/
     &     '    + facteur*dlog(Fcent))')

      return
      end
C------------------------------------------------------------------------
      subroutine ssh_WSPEC_RACM90 (nr,ispebp)
C------------------------------------------------------------------------
C
C     -- DESCRIPTION
C
C     Write kinetic rates case for specific reactions in RACM mechanism.
C
C------------------------------------------------------------------------
C
C     -- AUTHOR(S)
C
C     Pierre Plion and Bruno Sportisse, 2002.
C
C------------------------------------------------------------------------
      include 'nficfort'
      integer nr,ispebp

      If(ispebp.EQ.-1) then
         write(nficK90,11)nr
      elseif(ispebp.EQ.-2) then
         write(nficK90,12)nr
      elseif(ispebp.EQ.-3) then
         write(nficK90,13)nr
      elseif(ispebp.EQ.-4) then
         write(nficK90,14)nr
      elseif(ispebp.EQ.-5) then
         write(nficK90,15)nr
      elseif(ispebp.EQ.-6) then
         write(nficK90,16)nr
      elseif(ispebp.EQ.-7) then
         write(nficK90,17)nr
      else
         write(*,*) 'ERROR: unknown specific reaction ',ispebp
         stop 1
      endif

 11   format('  rk(',i4,') = SumM * 6.0d-34 * (temp/3.d2) ** (-2.3d0)')
C     BS 05/02/2003 values given by RACM
 12   format('  rk(',i4,') = 2.3d-13 * dexp(600.0d0 / temp) &',/
     &     '    + 1.73d-33* SumM * dexp(1000.0d0 / temp)')
 13   format('  rk(',i4,') = 3.22d-34 * dexp(2800.0d0 / temp) &',/
     &     '    + 2.38d-54 * SumM * dexp(3200.0d0 / temp)')
C     MODIF BS 06/06/2003 on the basis of CMAQ
 14   format('  Effko = 7.2d-15 * dexp(785.0d0 / temp)',/
     &     '  Rapk = 4.1d-16 * dexp(1440.0d0 / temp)',/
     &     '  facteur =1.9d-33 * dexp(725.0d0 / temp) * SumM',/
     &     '  rk(',i4,') = Effko + facteur/(1.0d0 + facteur / Rapk)')
 15   format('  rk(',i4,') = 1.5d-13 * (1.0d0 + 2.439d-20 * SumM)')
 16   format('  Rapk = 3.4d-30 * (300./temp)**(3.2)*SumM',/
     &     '  Effko = Rapk/(4.77D-11*(300./temp)**1.4)',/
     &     '  rk(',i4,')=(Rapk/(1.+Effko))*0.3** &',/
     &     '    (1.0d0/(1.0d0+ ((log10(Effko)-0.12)/1.2)**2))')
 17   format('  rk(',i4,') = 2.0d-39 * YlH2O * YlH2O')
      return
      end

C------------------------------------------------------------------------
      subroutine ssh_WSPEC_CB0590 (nr,ispebp)
C------------------------------------------------------------------------
C
C     -- DESCRIPTION
C
C     Write kinetic rates case for specific reactions in CB05 mechanism.
C
C------------------------------------------------------------------------
C
C     -- AUTHOR(S)
C
C     Youngseob KIM, 2010.
C
C------------------------------------------------------------------------
      include 'nficfort'
      integer nr,ispebp

      If(ispebp.EQ.-1) then
         write(nficK90,11)nr
      elseif(ispebp.EQ.-2) then
         write(nficK90,12)nr
      elseif(ispebp.EQ.-3) then
         write(nficK90,13)nr
      elseif(ispebp.EQ.-4) then
         write(nficK90,14)nr
      elseif(ispebp.EQ.-5) then
         write(nficK90,15)nr
      elseif(ispebp.EQ.-6) then
         write(nficK90,16)nr
      elseif(ispebp.EQ.-7) then
         write(nficK90,17)nr
      else
         write(*,*) 'ERROR: unknown specific reaction ',ispebp
         stop 1
      endif

 11   format('  rk(',i4,') = SumM * 6.0d-34 * (temp/3.d2) ** (-2.4d0)')
C     YS 17/11/2008 values given by NASA/JPL 2003
 12   format('  rk(',i4,') = 2.3d-13 * dexp(600.0d0 / temp) &',/
     &     '    + 1.7d-33* SumM * dexp(1000.0d0 / temp)')
 13   format('  rk(',i4,') = 3.22d-34 * dexp(2800.0d0 / temp) &',/
     &     '    + 2.38d-54 * SumM * dexp(3200.0d0 / temp)')
 14   format('  Effko = 2.4d-14 * dexp(460.0d0 / temp)',/
     &     '  Rapk = 2.7d-17 * dexp(2199.0d0 / temp)',/
     &     '  facteur =6.5d-34 * dexp(1335.0d0 / temp) * SumM',/
     &     '  rk(',i4,') = Effko + facteur/(1.0d0 + facteur / Rapk)')
C     YS 17/11/2008 values given by NASA/JPL 2003
 15   format('  rk(',i4,') = 1.44d-13 * (1.0d0 + 2.381d-20 * &', /
     &     '    8.0d-1 * SumM)')
C     YS 19/11/2008 value given by IUPAC 2005
 16   format('  Rapk = 3.4d-30 * (300./temp)**(3.2)*SumM',/
     &     '  Effko = Rapk/(4.77D-11*(300./temp)**1.4)',/
     &     '  rk(',i4,')=(Rapk/(1.+Effko))*0.3** &',/
     &     '    (1.0d0/(1.0d0+ ((log10(Effko)-0.12)/1.2)**2))')
 17   format('  rk(',i4,') = 1.8d-39 * YlH2O * YlH2O')
C     YS 26/11/2008 value given by IUPAC 2005
      return
      end

C------------------------------------------------------------------------
      subroutine ssh_WSPEC_RACM290 (nr,ispebp)
C------------------------------------------------------------------------
C
C     -- DESCRIPTION
C
C     Write kinetic rates case for specific reactions in RACM2 mechanism.
C
C------------------------------------------------------------------------
C
C     -- AUTHOR(S)
C
C     Youngseob KIM, 2010.
C     Victor Lannuque, 2020.
C
C------------------------------------------------------------------------
      include 'nficfort'
      integer nr,ispebp

      If(ispebp.EQ.-1) then
         write(nficK90,11)nr
      elseif(ispebp.EQ.-2) then
         write(nficK90,12)nr
      elseif(ispebp.EQ.-3) then
         write(nficK90,13)nr
      elseif(ispebp.EQ.-4) then
         write(nficK90,14)nr
      elseif(ispebp.EQ.-5) then
         write(nficK90,15)nr
      elseif(ispebp.EQ.-6) then
         write(nficK90,16)nr
      elseif(ispebp.EQ.-7) then
         write(nficK90,17)nr
      elseif(ispebp.EQ.-8) then
         write(nficK90,18)nr
      else
         write(*,*) 'ERROR: unknown specific reaction ',ispebp
         stop 1
      endif

 11   format('  rk(',i4,') = SumM * 5.74d-34 * (temp/3.d2) &',/
     &     '     ** (-2.6d0)')
C     YS(16/02/2009): values given by IUPAC 2006
 12   format('  rk(',i4,') = 2.2d-13 * dexp(600.0d0 / temp) &',/
     &     '    + 1.9d-33* SumM * dexp(980.0d0 / temp)')
 13   format('  rk(',i4,') = 3.08d-34 * dexp(2800.0d0 / temp) &',/
     &     '    + 2.59d-54 * SumM * dexp(3180.0d0 / temp)')
 14   format('  Effko = 2.4d-14 * dexp(460.0d0 / temp)',/
     &     '  Rapk = 2.7d-17 * dexp(2199.0d0 / temp)',/
     &     '  facteur =6.5d-34 * dexp(1335.0d0 / temp) * SumM',/
     &     '  rk(',i4,') = Effko + facteur/(1.0d0 + facteur / Rapk)')
C     Modif (YK:201108/09)): 2.4d-17 to 2.7d-17 based on RACM2_5L version
 15   format('  rk(',i4,') = 1.44d-13 * (1.0d0 + 8.0d-1 * SumM &', /
     &     '    / 4.0d19)')
C     VICTOR2020: modified according to RACM2 v2013 :
 16   format('  Rapk = 3.43d-12 * dexp(270./temp)',/
     &     '  Effko = (530./temp) + 4.8d-6 * Press - 1.73',/
     &     '  rk(',i4,')= Rapk * Effko / 100.')
 17   format('  rk(',i4,') = 1.8d-39 * YlH2O * YlH2O')
C     VICTOR2020: modified according to RACM2 v2013 :
 18   format('  rk(',i4,') = 1.39d-13 + 3.72d-11 * dexp(-2044./temp)')
C     YK 30/08/2010
      return
      end

C------------------------------------------------------------------------
      subroutine ssh_WSPEC_MCM (nr,ispebp)
C------------------------------------------------------------------------
C
C     -- DESCRIPTION
C
C     Write kinetic rates case for specific reactions in MCM mechanism.
C
C------------------------------------------------------------------------
C
C     -- AUTHOR(S)
C
C     Zhizhao WANG, 2022
C
C------------------------------------------------------------------------
      include 'nficfort'
      integer nr,ispebp

      If(ispebp.EQ.-8) then !genoa KMT13
         write(nficK90,18)nr
      elseif(ispebp.EQ.-9) then !genoa KMT14
         write(nficK90,19)nr
      elseif(ispebp.EQ.-10) then !genoa add SPEC '3.8D-13*dexp(780/TEMP)*(1-1/(1+498*dexp(-1160/TEMP)))'
         write(nficK90,20)nr
      elseif(ispebp.EQ.-11) then !genoa add SPEC '3.8D-13*dexp(780/TEMP)*(1/(1+498*dexp(-1160/TEMP)))'
         write(nficK90,21)nr
      elseif(ispebp.EQ.-12) then !genoa add SPEC '1.03D-13*math.dexp(365/TEMP)*(1-7.18*dexp(-885/TEMP))'
         write(nficK90,22)nr
      elseif(ispebp.EQ.-13) then !genoa add SPEC '8.8D-12*dexp(-1320/TEMP) + 1.7D-14*dexp(423/TEMP)'
         write(nficK90,23)nr
      elseif(ispebp.EQ.-14) then !genoa add SPEC '5.00E-12*O2*3.2*(1-dexp(-550/TEMP))'
         write(nficK90,24)nr

      elseif(ispebp.EQ.-15) then !genoa add SPEC '1.03D-13*0.5*math.dexp(365/TEMP)*(1-7.18*dexp(-885/TEMP))' ! for seperate
         write(nficK90,25)nr

      ! for MT - MCM
      elseif(ispebp.EQ.-16) then !genoa add SPEC '2.20E+10*EXP(-8174/TEMP)*EXP(1.00E+8/TEMP@3)'
         write(nficK90,26)nr
      elseif(ispebp.EQ.-17) then !genoa add SPEC '8.14E+9*EXP(-8591/TEMP)*EXP(1.00E+8/TEMP@3)'
         write(nficK90,27)nr
      endif

C genoa add formula for KMT13
 18   format('  rk(',i4,') = (4.45D-41*SumM*(temp/300) &',/
     &     '    **(-5.5d0))*(0.36d0/(1d0+(log10(1.39D-19* &',/
     &     '    SumM*(temp/300d0)**(-5.5d0))/1.3135d0)**2d0))/ &',/
     &     '    (2.5D-30*SumM*(temp/300d0)**(-5.5d0)+1.8D-11)')
C genoa add formula for KMT14
 19   format('  rk(',i4,') = (9.0D-5*dexp(-9690d0/temp)*SumM)* &',/
     &     '    (1.1D+16d*dexp(-1056d1/temp))*10d**(dlog10(4d-1)/&',/
     &     '    (1d0+(dlog10((9.0D-5*dexp(-9690d0/temp)*SumM)/ &',/
     &     '    (1.1D+16*dexp(-10560d0/temp)))/1.2554d0)**2d0))/ &',/
     &     '    (9.0D-5*dexp(-9690d0/temp) &',/
     &     '    *SumM+1.1D+16*dexp(-10560d0/temp))')
C genoa add formula for '3.8D-13*dexp(780/TEMP)*(1-1/(1+498*dexp(-1160/TEMP)))'
 20   format('  rk(',i4,') = 3.8D-13*dexp(780d0/temp)* &',/
     &     '    (1d0-1d0/(1d0+498d0*dexp(-1160/temp)))')
C genoa add formula for '3.8D-13*dexp(780/TEMP)*(1/(1+498*dexp(-1160/TEMP)))'
 21   format('  rk(',i4,') = 3.8D-13*dexp(780d0/TEMP)* &',/
     &     '    (1d0/(1d0+498d0*dexp(-1160d0/TEMP)))')
C genoa add formula for '2*1.03D-13*dexp(365/TEMP)*(1-7.18*dexp(-885/TEMP))'
 22   format('  rk(',i4,') = 2d0*1.03D-13*dexp(365d0/TEMP)* &',/
     &     '    (1d0-7.18d0*dexp(-885d0/TEMP))')
C genoa add formula for '8.8D-12*dexp(-1320/TEMP) + 1.7D-14*dexp(423/TEMP)'
 23   format('  rk(',i4,') = 8.8D-12*dexp(-1320d0/TEMP) +  &',/
     &     '    1.7D-14*dexp(423d0/TEMP)')
C genoa add formula for '5.00E-12*O2*3.2*(1-dexp(-550/TEMP))'
 24   format('  rk(',i4,') = 5.00D-12*3.2d0  &',/
     &     '    *(1d0-dexp(-550d0/TEMP))')
C for seperate
C genoa add formula for '2*1.03D-13*0.5*dexp(365/TEMP)*(1-7.18*dexp(-885/TEMP))'
 25   format('  rk(',i4,') = 1.03D-13*dexp(365d0/TEMP)* &',/
     &     '    (1d0-7.18d0*dexp(-885d0/TEMP))')
C for Limonene species INDO
C genoa add formula for '2.20E+10*EXP(-8174/TEMP)*EXP(1.00E+8/TEMP@3)'
 26   format('  rk(',i4,') =  2.2d10*dexp(-8174d0/TEMP)* &',/
     &     '    dexp(1d8/TEMP**3)')
C genoa add formula for '8.14E+9*EXP(-8591/TEMP)*EXP(1.00E+8/TEMP@3)'
 27   format('  rk(',i4,') = 8.14d9*dexp(-8591/TEMP)* &',/
     &     '    dexp(1d8/TEMP**3)')

      return
      end
C      
C   

C------------------------------------------------------------------------
      subroutine ssh_WWALL90 (nr,b1,b2,b3,b4,b5)
C------------------------------------------------------------------------
C
C     -- DESCRIPTION
C
C     Write kinetic rates case for wall transfer.
C     possible to consider only tranfer from gas toward wall or also from wall toward gas
C     See Huang et al. (2018) for more informations (doi: 10.1021/acs.est.7b05575)
C     
C     WARNING: Should be implemented in SOAP and not here 
C
C     PARAMETERS:
C
C      options:
C       - flag on(0)/off(1) wall  = b1
C
C      species depending:
C       - Psat298K(i)             = b2
C       - DHvap(i)                = b3
C       - MW(i)                   = b5
C
C      chamber depending:
C       - surface A [m2]
C       - volume V [m3]
C       - Mw equivalent of wall Mwall [g/mol]
C       - activity coef in wall gammawall
C       - concentration equivalent of wall Cwall [g/m3]
C       - wall accomodation param1 avar
C       - wall accomodation param2 bvar
C       - eddy diffusivity ke (V) [s-1] = b4
C       - AonV = A/V [m-1]
C
C      depending of both:
C       - wall accomodation coef aw (Psat,dHvap,Mw,T,avar,bvar)
C
C      other:
C        - diffusivity in gas Dgas [m2/s] 

C
C------------------------------------------------------------------------
C
C     -- AUTHOR(S)
C
C     VICTOR LANNUQUE, 2021.
C
C------------------------------------------------------------------------

      include 'nficfort'
      integer nr,b1
      double precision b2,b3,b4,b5,Mwall,gammawall,Cwall
      double precision AonV,cbar,Dgas,avar,bvar,aw
      
      Dgas = 5.0d-6

c For barbara D'Anna AFT (A=0.5730m2  V=0.01719m3) :
c   -Gas to wall
      AonV = 33.33
      avar = -2.744
      bvar = -1.4066
c   -Wall to gas (Ne pas utiliser tel quel, A REVOIR!!!)
      Mwall = 200.
      gammawall = 1.0
      Cwall = 0.360


C= gas vers wall kw,on      
      If(b1.EQ.0) then
         write(nficK90,11)b2,b3,b5,avar,bvar,b4,Dgas,AonV
         write(nficK90,12)nr
         
C= wall vers gas kw,off
      elseif(b1.EQ.1) then
         write(nficK90,11)b2,b3,b5,avar,bvar,b4,Dgas,AonV
         write(nficK90,13)Mwall,gammawall,Cwall
         write(nficK90,14)nr         
         
      else
         write(*,*) 'ERROR: unknown wall transfer reaction ',nr
         stop 1
      endif         
               
 11   format('  Psat = ',D12.5,'/760.*exp(',D12.5,'*1000. &',/
     & '    / 8.314*(1./298.-1./temp))',/
     & '  Masmol = ',D12.5,' ',/
     & '  cstar = Masmol*1d+6*Psat/(8.205d-5*temp)',/
     & '  aw = 10**(',D12.5,')*cstar**(',D12.5,')',/
     & '  cbar = (8000*8.314*temp/(3.14*Masmol))**0.5',/  
     & '  awc = (aw*cbar/4.)',/
     & '  denom = awc / ((',D12.5,'*',D12.5,')**0.5)',/
     & '  kwon = ',D12.5,' * awc / (1.+1.5708*denom)')
     
 12   format('  rk(',i4,') = kwon')     

 13   format('  facteur=',D12.5,'*',D12.5,'/(8.205d-5 * temp &',/
     & '    * ',D12.5,')')
 
 14   format('  rk(',i4,') = kwon * Psat * facteur')      

      return
      end
      



C------------------------------------------------------------------------
      subroutine ssh_IRDICARB (nr,b1)
C------------------------------------------------------------------------
C
C     -- DESCRIPTION
C
C     Write kinetic rate corresponding to irreversible pathway for MGLY condensation.
C     depending on RH
C     
C     See SI of Lannuque et al. (2023) for more informations
C
C
C------------------------------------------------------------------------
C
C     -- AUTHOR(S)
C
C     VICTOR LANNUQUE, 2022.
C
C------------------------------------------------------------------------
      include 'nficfort'
      integer nr,b1

      If(b1.EQ.-1) then
         write(nficK90,11)nr
      endif  

 11   format('  Psat = 611.2d0 * dexp(17.67d0 * (temp - 273.15d0) &',/
     & '    /(temp - 29.65d0))',/
     & '  facteur = xlw*Press / ((0.62197d0 * (1.d0-xlw)+xlw) &',/
     & '  * Psat ) * 100.',/
     & '  rk(',i4,') = 1d-9*facteur**3 - 1d-7*facteur**2 &',/
     & '  + 3.0d-7*facteur + 0.0003 ') 
     
      return
      end
      

   
C------------------------------------------------------------------------
      subroutine ssh_WHETERO90 (nr,ihetero)
C------------------------------------------------------------------------
C
C     -- DESCRIPTION
C
C     Associate kinetic rates of heterogeneous reactions with right variable
C     for hetrxn.f.
C
C     HETERO -1 : HO2  -->  0.5 H2O2
C     HETERO -2 : NO2  -->  0.5 HONO + 0.5 HNO3
C     HETERO -3 : NO3  -->  HNO3
C     HETERO -4 : N2O5 -->  2 HNO3
C
C------------------------------------------------------------------------
C
C     -- AUTHOR(S)
C
C     VICTOR LANNUQUE, 2020.
C
C------------------------------------------------------------------------
      include 'nficfort'
      integer nr,ihetero

      If(ihetero.EQ.-1) then
         write(nficK90,11)nr
      elseif(ihetero.EQ.-2) then
         write(nficK90,12)nr
      elseif(ihetero.EQ.-3) then
         write(nficK90,13)nr
      elseif(ihetero.EQ.-4) then
         write(nficK90,14)nr

      else
         write(*,*) 'ERROR: unknown heterogeneous reaction ',ihetero
         stop 1
      endif
      
      
!     HETERO -1 : HO2  -->  0.5 H2O2
 11   format('  rk(',i4,') = rkhHO2')
!     HETERO -2 : NO2  -->  0.5 HONO + 0.5 HNO3
 12   format('  rk(',i4,') = rkhNO2')
!     HETERO -3 : NO3  -->  HNO3
 13   format('  rk(',i4,') = rkhNO3')
!     HETERO -4 : N2O5 -->  2 HNO3
 14   format('  rk(',i4,') = rkhN2O5')

      return
      end
      
      
      subroutine ssh_WSPEC_MELCHIOR2 (nr,ispebp)
C------------------------------------------------------------------------
C
C     -- DESCRIPTION
C
C     Write kinetic rates case for specific reactions in
C     MELCHIOR2 mechanism.
C     
C------------------------------------------------------------------------
C
C     -- AUTHOR(S)
C
C     Yelva Roustan, 2020.
C
C------------------------------------------------------------------------
      include 'nficfort'
      integer nr,ispebp

      If(ispebp.EQ.-6) then
         write(nficK90, 16) nr
      elseif(ispebp.EQ.-7) then
         write(nficK90, 17) nr
      elseif(ispebp.EQ.-9) then
         write(nficK90, 19) nr
      else
         write(*,*) 'ERROR: unknown specific reaction ',ispebp
         stop 1
      endif

C     YS 19/11/2008 value given by IUPAC 2005
C     NO2 + OH + M -> HNO3
C     k(T, M) = mtroe(3.4e-30, 0, 3.2, 4.77e-11, 0, 1.4, 0.30)
 16   format(6x,'Rapk = 3.4d-30 * (300./temp)**(3.2)*SumM',/
     &     6x,'Effko = Rapk/(4.77D-11*(300./temp)**1.4)',/
     &     6x,'rk(',i3,')=(Rapk/(1.+Effko))*0.3** &',/
     &     5x,'(1.0d0/(1.0d0+ ((dlog10(Effko)-0.12)/1.2)**2))')
C     YR 28/04/2020 value given by REACTIONS.univ.melchior2
C     N2O5 + H2O + H2O -> 2*HNO3      k=2e-39
 17   format(6x,'rk(',i3,') = 2.0d-39 * YlH2O * YlH2O')
C     YR 28/04/2020 value given by REACTIONS.univ.melchior2
C     Surface production by heterogeneous reaction (Aumont et al., 2003)
C     NO2 -> HONO + NO2      ks = 0.5 * vd_NO2 / h
C     with vd_NO2 = 0.033 cm / s in urban area.
C     with h_mixing = 15 m in urban area.
 19   format(6x,'h_mixing  = 15.0d0',/
     &     6x,'vd_NO2 = 0.00033d0',/
     &     6x,'rk(',i3,') = 0.5d0 * vd_NO2 / h_mixing')
      return
      end







C------------------------------------------------------------------------
      subroutine ssh_WTB90(nr,ittb)
C------------------------------------------------------------------------
C
C     -- DESCRIPTION
C
C     Write third body kinetics.
C
C------------------------------------------------------------------------
C
C     -- AUTHOR(S)
C
C     Pierre Plion and Bruno Sportisse, 2002.
C
C------------------------------------------------------------------------
      include 'nficfort'
      integer nr,ittb

      If(ittb.EQ.1) then
         write(nficK90,11)nr,nr
      elseif (ittb.EQ.2) then
         write(nficK90,12)nr,nr
      elseif (ittb.EQ.3) then
         write(nficK90,13)nr,nr
      elseif (ittb.EQ.4) then
         write(nficK90,14)nr,nr
      elseif (ittb.EQ.5) then
         write(nficK90,15)nr,nr
      elseif (ittb.EQ.6) then !genoa add TB RO2 
         write(nficK90,16)nr,nr
      endif

 11   format('  rk(',i4,') = rk(',i4,') * SumM')
C     Seinfeld pp 22: N2 0.78; O2 0.21
 12   format('  rk(',i4,') = rk(',i4,') * SumM * 0.2d0')
 13   format('  rk(',i4,') = rk(',i4,') * SumM * 0.8d0')
 14   format('  rk(',i4,') = rk(',i4,') * YlH2O')
 15   format('  rk(',i4,') = rk(',i4,') * SumM * 5.8d-7')
 16   format('  rk(',i4,') = rk(',i4,') * RO2') !genoa TB RO2
      return
      end
C------------------------------------------------------------------------
      subroutine ssh_WCV90(nr,b1,b2,b3,b4,b5,b6,b7)
C------------------------------------------------------------------------
C
C     -- DESCRIPTION
C
C     Write kinetics for variable stoichiometry for MOCA.
C
C------------------------------------------------------------------------
C
C     -- AUTHOR(S)
C
C     Pierre Plion, 2002.
C
C------------------------------------------------------------------------
      include 'nficfort'
      integer nr
      double precision b1,b2,b3,b4,b5,b6,b7
      double precision d54,d65,d76

c      if(abs(b2).lt.petit)then
c         call ssh_WK290(nr,b1,b3)
c      else
c         call ssh_WK390(nr,b1,b2,b3)
c      endif
c     write(nficK90,10)nr,b4,b5,b6,b7
      d54 = (b5-b4)*.05D0
      d65 = (b6-b5)*.05D0
      d76 = (b7-b6)*.05D0
      write(nficK90,12)
      if(abs(b4).gt.petit)then
         write(nficK90,13)nr,nr,b4
      else
         write(nficK90,18)nr
      endif
      write(nficK90,14)260.,280.
      if(abs(b4).gt.petit .or. abs(d54).gt. petit)then
         write(nficK90,15)nr,nr,b4,260.,d54
      else
         write(nficK90,18)nr
      endif
      write(nficK90,14)280.,300.
      if(abs(b5).gt.petit .or. abs(d65).gt. petit)then
         write(nficK90,15)nr,nr,b5,280.,d65
      else
         write(nficK90,18)nr
      endif
      write(nficK90,14)300.,320.
      if(abs(b6).gt.petit .or. abs(d76).gt. petit)then
         write(nficK90,15)nr,nr,b6,300.,d76
      else
         write(nficK90,18)nr
      endif
      write(nficK90,16)
      if(abs(b7).gt.petit)then
         write(nficK90,13)nr,nr,b7
      else
         write(nficK90,18)nr
      endif
      write(nficK90,17)

 10   format('  !',6x,i4,4(2x,D23.16))
 12   format('  if (temp.le.260.d0) then ')
 13   format('    rk(',i4,') = rk(',i4,') *  ',D23.16)
 14   format('  elseif(temp.gt.',D10.3,'.and.temp.le.',
     &     D10.3,') then')
 15   format('    rk(',i4,') = rk(',i4,') * (',
     &     D23.16,' &',/' + (temp-',D10.3,') * (',D23.16,'))')
 16   format('  else')
 17   format('  endif',/'!')
 18   format('  rk(',i4,') = 0.d0')

      return
      end
C------------------------------------------------------------------------
      subroutine ssh_WRCFE90(nr,b1,b2,b3,b4,b5,b6)
C------------------------------------------------------------------------
C
C     -- DESCRIPTION
C
C     Write kinetics for reactions computed from equilibria.
C     k(T)=k0*M/(1+r)  * 0.6**(1/(1+[log10 r]**2)) * b5* dexp(-b6/T)
C     with r= k0*M/kinf, k0=b1*(T/300)**(-b2), kinf=b3*(T/300)**(-b4)

C------------------------------------------------------------------------
C
C     -- AUTHOR(S)
C
C     Pierre Plion and Bruno Sportisse, 2002.
C
C------------------------------------------------------------------------
C
C     -- Comment
C     2009/03/13(YS): Make sure that b5 and b6 are the parameters for
C                     the inverse reaction which means the thermal
C                     decomposition: C -> A + B
C
C------------------------------------------------------------------------
      include 'nficfort'
      integer nr
      double precision b1,b2,b3,b4,b5,b6

C
      call ssh_WK290(nr,b5,b6)
C
      write(nficK90,10)b1,b2
      write(nficK90,11)b3,b4
      write(nficK90,12)
      write(nficK90,13)nr,nr

 10   format('  Effko = ',D23.16,'* (temp / 3.d2) &',/
     &     '    **(- (',D23.16,'))')
 11   format('  Rapk = ',D23.16,'* (temp / 3.d2) &',/
     &     '    **(- (',D23.16,'))')
 12   format('  facteur = (Effko * SumM / ( 1.0d0 + Effko * SumM / &',/
     &     '    Rapk)) * 0.6d0 ** (1.0d0 / (1.0d0 + &',/
     &     '    (log10(Effko * SumM / Rapk))**2))')
 13   format('  rk(',i4,') = facteur * rk(',i4,')')
C     13    format(6x,'rk(',i4,') = facteur * (',D23.16,') * '/
C     &       5x,'&      dexp((',D23.16,')/temp)')


      return
      end
C------------------------------------------------------------------------
      subroutine ssh_WPHOT90(nr,ntabphot,b,tabphot,it)
C------------------------------------------------------------------------
C
C     -- DESCRIPTION
C
C     Write kinetics for tabulated photolysis at angles TABPHOT(I)
C     for 1<=I<=NTABPHOT.
C     The corresponding values for photolysis rates are BP(I,NR)
C
C     at 90: no photolysis.
C     at 0 and 90: first derivative = 0.
C
C     Interpolation with standard cubic spline (second derivative =0
C     at boundaries).
C
C     Correction factor according to IT (IT =1 : O3 >2 OH).
C
C------------------------------------------------------------------------
C
C     -- AUTHOR(S)
C
C     Bruno Sportisse, 2003.
C
C------------------------------------------------------------------------
      include 'parametre'
      include 'nficfort'
      integer nr,i,j,it
      double precision a(ntabphotmax),b(ntabphotmax),c(4,nintphotmax)
      double precision tabphot(ntabphotmax)

      do i=1,ntabphot
         a(i)=tabphot(i)
      enddo
C
      call ssh_SPL3(ntabphot,a,b,c)
C
c      do i = 1,ntabphot-1
c         write(6,100)b(i),(c(j,i),j=1,4),b(i+1)
c      enddo

c loop for photolysis options
      write(nficK90,131)
      write(nficK90,110)
      write(nficK90,111) nr
      write(nficK90,112)

c

      write(nficK90,12)
      do i = 1,ntabphot-1
         write(nficK90,11)a(i),a(i+1)
         write(nficK90,14)nr,c(4,i)
         write(nficK90,15)nr,c(3,i),a(i),nr
         write(nficK90,15)nr,c(2,i),a(i),nr
         write(nficK90,15)nr,c(1,i),a(i),nr
      enddo
      write(nficK90,16)nr,b(ntabphot)
      write(nficK90,13)
      write(nficK90,10)nr,nr
      write(nficK90,13)
      write(nficK90,131)
C
      if (it.eq.1) then
         write(nficK90,103)nr,nr
      endif

 110  format('  if(option_photolysis.eq.2) then')
 111  format('    rk(',i4,')= 0.D0')
 112  format('  elseif(option_photolysis.eq.1) then')

 10   format('  if(att.lt.0.99999) rk(',i4,') = rk(',i4,') * att')
 11   format('  elseif(azi.ge.',D9.2,
     &     ' .and. azi.lt.',D9.2,') then')
 14   format('    rk(',i4,')=',D23.16)
 15   format('    rk(',i4,')=',D23.16,'+(azi-',D9.2,') * rk(',i4,')')
 12   format('  if(azi.lt.0.D0)then',/,'    STOP')
 16   format('  elseif(azi.ge.90.D0)then',/,'    rk(',i4,')=',D23.16)
 13   format('  endif')
 131  format('  !')
 100  format(6(2x,D23.16))
C
 103  format('  VO2  = 3.2D-11 * dexp(70.D0 /temp) * (0.21D0*SumM)',/
     &     '  VN2  = 1.8D-11 * dexp(110.D0/temp) * (0.79D0*SumM)',/
     &     '  VH2O = 2.2D-10 * YlH2O'/
     &     '  rk(',i4,') = rk(',i4,') * VH2O / (VH2O + VN2 + VO2) ')

      return
      end
C------------------------------------------------------------------------
      subroutine ssh_WW90(nr,ne,ie1,ie2,ie3)
C------------------------------------------------------------------------
C
C     -- DESCRIPTION
C
C     Write reaction rates.
C
C------------------------------------------------------------------------
C
C     -- AUTHOR(S)
C
C     Pierre Plion, 2002.
C
C------------------------------------------------------------------------
      implicit none
      integer nr,ne,ie1,ie2,ie3
      include 'nficfort'

      if    (ne.eq.1)then
         write(nficw90,10)nr,nr,ie1
C     write(nficFF,100)nr,ie1
      elseif(ne.eq.2)then
         write(nficw90,11)nr,nr,ie1,ie2
C     write(nficFF,110)nr,ie1,ie2
      elseif(ne.eq.3)then
         write(nficw90,12)nr,nr,ie1,ie2,ie3
C     write(nficFF,120)nr,ie1,ie2,ie3
      endif

 10   format('  w(',i4,') =  rk(',i4,') * Y(',i4,')')
 11   format('  w(',i4,') =  rk(',i4,') * Y(',i4,')',
     $     ' * Y(',i4,')' )
 12   format('  w(',i4,') =  rk(',i4,') * Y(',i4,')',
     $     ' * Y(',i4,')',
     $     ' * Y(',i4,')' )

      return
      end
C------------------------------------------------------------------------
      subroutine ssh_DW90(nr,ne,ie1,ie2,ie3)
C------------------------------------------------------------------------
C
C     -- DESCRIPTION
C
C     Write derivative of reaction rates.
C
C------------------------------------------------------------------------
C
C     -- AUTHOR(S)
C
C     Pierre Plion, 2002.
C
C------------------------------------------------------------------------
      implicit none
      integer nr,ne,ie1,ie2,ie3
      include 'nficfort'

      if    (ne.eq.1)then
         write(nficdw90,10)nr,ie1,nr
C     write(nficJJ,100)nr
      elseif(ne.eq.2)then
         write(nficdw90,11)nr,ie1,nr,ie2
C     write(nficJJ,110)nr,ie2
         write(nficdw90,11)nr,ie2,nr,ie1
C     write(nficJJ,111)nr,ie1
      elseif(ne.eq.3)then
         write(nficdw90,12)nr,ie1,nr,ie2,ie3
C     write(nficJJ,120)nr,ie2,ie3
         write(nficdw90,12)nr,ie2,nr,ie1,ie3
C     write(nficJJ,121)nr,ie1,ie3
         write(nficdw90,12)nr,ie3,nr,ie1,ie2
C     write(nficJJ,122)nr,ie1,ie2
      endif

 10   format('  dw(',i4,',',i4,') =  rk(',i4,')')
 11   format('  dw(',i4,',',i4,') =  rk(',i4,') * Y(',i4,')' )
 12   format('  dw(',i4,',',i4,') =  rk(',i4,') * Y(',i4,')',
     $     ' * Y(',i4,')' )

      return
      end
C------------------------------------------------------------------------
      subroutine ssh_WFJ90(s,nr,jer)
C------------------------------------------------------------------------
C
C     -- DESCRIPTION
C
C     Update chemical production term and Jacobian Matrix.
C
C------------------------------------------------------------------------
C
C     -- AUTHOR(S)
C
C     Pierre Plion, 2002.
C
C------------------------------------------------------------------------
      implicit none
      include 'nficfort'
      include 'parametre'
      double precision s(nespmax,nrmax)
      integer nr,jer(3,nrmax)
      integer ie,je

      do ie = 1,nespmax
         if(abs(s(ie,nr)).gt.petit)then
C     Products
            if(s(ie,nr).gt.0.d0)then
C     Products with  stoichiometry = 1
               if(abs(s(ie,nr)-1.d0).lt.petit)then
                  write(nficF90,30)ie,ie,nr
                  do je = 1,3
                     if(jer(je,nr).gt.0) then
                        write(nficJ90,31)
     $                       ie,jer(je,nr),ie,jer(je,nr),nr,jer(je,nr)
                     endif
                  enddo
               else
                  write(nficF90,10)ie,ie,s(ie,nr),nr
                  do je = 1,3
                     if(jer(je,nr).gt.0)
     &                    write(nficJ90,11)ie,jer(je,nr),ie,jer(je,nr),
     &                    s(ie,nr),nr,jer(je,nr)
                  enddo
               endif
            else
C     Reactants
               if(abs(s(ie,nr)+1.d0).lt.petit)then
C     Reactants with stochiometry = 1
                  write(nficF90,40)ie,ie,nr
                  do je = 1,3
                     if(jer(je,nr).gt.0)
     &                    write(nficJ90,41)
     $                    ie,jer(je,nr),ie,jer(je,nr),nr,jer(je,nr)
                  enddo
               else
                  write(nficF90,20)ie,ie,(-s(ie,nr)),nr
                  do je = 1,3
                     if(jer(je,nr).gt.0)
     &                    write(nficJ90,21)ie,jer(je,nr),ie,jer(je,nr),
     &                    (-s(ie,nr)),nr,jer(je,nr)
                  enddo
               endif
            endif
         endif
      enddo

 10   format('  chem(',i4,') = chem(',i4,') + ',D23.16,' * w(',i4,')')
 11   format('  JacC(',i4,',',i4,') = JacC(',i4,',',i4,')+',
     &     D23.16,'*dw(',i4,',',i4,')')

 20   format('  chem(',i4,') = chem(',i4,') - ',D23.16,' * w(',i4,')')
 21   format('  JacC(',i4,',',i4,') = JacC(',i4,',',i4,')-',
     &     D23.16,'*dw(',i4,',',i4,')')

 30   format('  chem(',i4,') = chem(',i4,') + w(',i4,')')
 40   format('  chem(',i4,') = chem(',i4,') - w(',i4,')')
 31   format('  JacC(',i4,',',i4,') = JacC(',i4,',',i4,
     &     ') + dw(',i4,',',i4,')')
 41   format('  JacC(',i4,',',i4,') = JacC(',i4,',',i4,
     &     ') - dw(',i4,',',i4,')')


      RETURN
      END
C------------------------------------------------------------------------
      subroutine ssh_WPL90(s,nr,jer)
C------------------------------------------------------------------------
C
C     -- DESCRIPTION
C
C     Update production and loss terms (P-Lc formulation).
C
C------------------------------------------------------------------------
C
C     -- AUTHOR(S)
C
C     Bruno Sportisse, 2003.
C
C------------------------------------------------------------------------
      implicit none
      include 'nficfort'
      include 'parametre'
      double precision s(nespmax,nrmax)
      integer nr,jer(3,nrmax)
      integer ie

      do ie = 1,nespmax
         if(abs(s(ie,nr)).gt.petit)then
C     Products
            if(s(ie,nr).gt.0.d0)then
C     Products with  stoichiometry = 1
               if(abs(s(ie,nr)-1.d0).lt.petit)then
                  write(nficprod90,30)ie,ie,nr
               else
                  write(nficprod90,10)ie,ie,s(ie,nr),nr
               endif
            else
C     Reactants
               if(abs(s(ie,nr)+1.d0).lt.petit)then
C     Reactants with stochiometry = 1
                  write(nficloss90,40)ie,ie,nr,ie
               else
                  write(nficloss90,20)ie,ie,(-s(ie,nr)),nr,ie
               endif
C     Correction if bimolecular reactant: none
C     Because d(c*c)/dc=c and not 2*c as computed !
C     nm=0
C     do je=1,3
C     if (jer(je,nr).eq.ie) nm=nm+1
C     enddo
C     if (nm.gt.1) write(nficloss90,50)ie,ie,nm*1.
            endif
         endif
      enddo

 10   format('  prod(',i4,') = prod(',i4,') + ',D23.16,' * w(',i4,')')
 30   format('  prod(',i4,') = prod(',i4,') + w(',i4,')')

 20   format('  loss(',i4,') = loss(',i4,') + ',D23.16,' * dw(',i4,',',
     &     i4,')')
 40   format('  loss(',i4,') = loss(',i4,') + dw(',i4,',',i4,')')

 50   format('  loss(',i4,') = loss(',i4,')/',D23.16)

      RETURN
      END
C------------------------------------------------------------------------
      subroutine ssh_WNONZERO90(s,nrtot,jer)
C------------------------------------------------------------------------
C
C     -- DESCRIPTION
C
C     Create file for nonzero entries of Jacobian matrix:
C     non_zero.dat
C
C------------------------------------------------------------------------
C
C     -- AUTHOR(S)
C
C     Bruno Sportisse, 2003.
C
C------------------------------------------------------------------------
      implicit none
      include 'nficfort'
      include 'parametre'
      double precision s(nespmax,nrmax),jnz(nespmax,nespmax)
      integer nr,nrtot,jer(3,nrmax)
      integer ie,je,nztot

      do ie=1,nespmax
         do je=1,nespmax
            jnz(ie,je)=0
         enddo
      enddo
c
      do nr=1,nrtot
         do ie = 1,nespmax
            if(abs(s(ie,nr)).gt.petit)then
C     Products
               if(s(ie,nr).gt.0.d0)then
                  do je = 1,3
                     if(jer(je,nr).gt.0) then
                        jnz(ie,jer(je,nr))=1
                     endif
                  enddo
               else
C     Reactants
                  do je = 1,3
                     if(jer(je,nr).gt.0) then
                        jnz(ie,jer(je,nr))=1
                     endif
                  enddo
               endif
            endif
         enddo
      enddo
C
      nztot=0
      do ie=1,nespmax
         do je=1,nespmax
            if (jnz(ie,je).eq.1) nztot=nztot+1
         enddo
      enddo
C
      write(nficnz90,*)nztot
      do ie=1,nespmax
         do je=1,nespmax
            if (jnz(ie,je).eq.1) write(nficnz90,*)ie,je
         enddo
      enddo

      RETURN
      END
C------------------------------------------------------------------------
      subroutine ssh_SPL3(n1,a,b,c)
C------------------------------------------------------------------------
C     Determination des coefficients des arcs de cubiques permettant
C     d'interpoler b(a) en respectant la condition de derivee premiere
C     nulle au bord
C
C     n1 nombre de points, n0 nombre de segments
C
C     Pour a(i) < aa < a(i+1)
C     On ecrira
C     bb = c(1,i) + c(2,i)*aa + c(3,i)*aa**2 + c(4,i)*aa**3
C     ou aa = angle -a(i)
C
C     c(1,i) = b(i)
C
C     bp = c(2,i) + 2*c(3,i)*aa + 3*c(4,i)*aa**2
C     bs = 2*c(3,i) + 6*c(4,i)*aa
C
C     l'annulation de la derivee premiere au premier point fournit
C     c(2,1) = 0
C
C     1<i<n0
C     la continuite de la fonction s'ecrit
C     bb(i,a(i+1) = bb(i+1,a(i+1)) = b(i+1)
C     c(1,i)+c(2,i)*d(i)+c(3,i)*d(i)**2+c(4,i)*d(i)**3 = b(i+1)
C     ou d(i) = a(i+1)-a(i)
C
C     la continuite des derivees premieres s'ecrit
C
C     1<i<(n0-1)
C     bp(i,a(i+1)) = bp(i+1,a(i+1))
C     c(2,i) +2*c(3,i)*d(i)+3*c(4,i)*d(i)**2 = c2(i+1)
C
C     la continuite des derivees secondes s'ecrit
C
C     1<i<(n0-1)
C     bs(i,a(i+1)) = bs(i+1,a(i+1))
C     c(3,i)+3.*c(4,i)*d(i) = c(3,i+1)
C
C     L'annulation de la derivee premiere au dernier point s'ecrit
C
C     bp(n0,a(n1)) = 0
C     c(2,n0)+2.*c(3,n0)*d(n0)+3.*c(4,n0)*d(n0)**2 = 0
C
C     i/ La continuite de la fonction sert de relation de recurrence
C     pour calculer c(4,i) en fonction de c(3,i)
C     ii/ La continuite de la derivee premiere sert de relation de recurrence
C     pour calculer c(2,i+1) en fonction de c(2,i),c(3,i),c(4,i)
C     iii/ La continuite de la derivee seconde sert de relation de recurrence
C     pour calculer c(3,i+1) en fonction de c(3,i),c(4,i)
C
C     On utilise donc comme inconnue auxilliaire c(3,1) qui sera determine
C     par la condition sur la derivee premiere au dernier point
C------------------------------------------------------------------------
      implicit none
      include 'parametre'
      integer i,n1,n0
      double precision a(ntabphotmax),b(ntabphotmax)
      double precision c(4,nintphotmax),d(ntabphotmax)
      double precision c20(nintphotmax),c21(nintphotmax)
      double precision c30(nintphotmax)
      double precision c31(nintphotmax),c40(nintphotmax)
      double precision c41(nintphotmax),cc31

      n0=n1-1
      do i = 1,n0
         c(1,i) = b(i)
         d(i) = a(i+1)-a(i)
      enddo
      c20(1) = 0.D0
      c21(1) = 0.D0
      c30(1) = 0.D0
      c31(1) = 1.D0
      c40(1) = (b(2  )-c(1,1)-c20(1)*d(1)-c30(1)*d(1)**2)/d(1)**3
      c41(1) = (             -c21(1)*d(1)-c31(1)*d(1)**2)/d(1)**3
      do i = 2,n0
C     c(2,i) +2*c(3,i)*d(i)+3*c(4,i)*d(i)**2 = c2(i+1)
         c20(i) = c20(i-1)+2.D0*c30(i-1)*d(i-1)+3.D0*c40(i-1)*d(i-1)**2
         c21(i) = c21(i-1)+2.D0*c31(i-1)*d(i-1)+3.D0*c41(i-1)*d(i-1)**2
C     c(3,i)+3.*c(4,i)*d(i) = c(3,i+1)
         c30(i) = c30(i-1) +3.D0*c40(i-1)*d(i-1)
         c31(i) = c31(i-1) +3.D0*c41(i-1)*d(i-1)
C     c(1,i)+c(2,i)*d(i)+c(3,i)*d(i)**2+c(4,i)*d(i)**3 = b(i+1)
C     c4(i) = (b(i+1)-c(1,i)-c(2,i)*d(i)-c(3,i)*d(i)**2)/d(i)**3
         c40(i) = (b(i+1)-c(1,i)-c20(i)*d(i)-c30(i)*d(i)**2)/d(i)**3
         c41(i) = (             -c21(i)*d(i)-c31(i)*d(i)**2)/d(i)**3
      enddo
C     c(2,n0)+2.*c(3,n0)*d(n0)+3.*c(4,n0)*d(n0)**2 = 0
C     (c20(n0)+c21(n0)*cc31) + 2.*(c30(n0)+c31(n0)*cc31)....
      cc31 = c20(n0)+2.D0*c30(n0)*d(n0)+3.D0*c40(n0)*d(n0)**2
      cc31 = -cc31
     &     /(c21(n0)+2.D0*c31(n0)*d(n0)+3.D0*c41(n0)*d(n0)**2)
C
      do i = 1,n0
         c(2,i) = c20(i)+cc31*c21(i)
         c(3,i) = c30(i)+cc31*c31(i)
         c(4,i) = c40(i)+cc31*c41(i)
      enddo

      RETURN
      END
