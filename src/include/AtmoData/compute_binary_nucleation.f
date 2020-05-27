C-----------------------------------------------------------------------
C     Copyright (C) 2019 CEREA (ENPC) - INERIS
C     SSH-aerosol is distributed under the GNU General Public License v3
C-----------------------------------------------------------------------

c     Function: compute_binary_nucleation_kernel
c
c     Computes binary nucleation rate.
c     Based on the parameterization of
c     Veahkamaki et al, JGR Atmosphere, vol 107, D22, p4622, 2002.
c     Modifications:
c     2005/3/23: cleaning Bruno Sportisse.
c
c     Parameters:
c     rh - Relative humidity 0< <1 ().
c     temperature - Temperature (K).
c     natmp - Gas-phase H2SO4 concentration (#molec/cm^3).
c
c     Returns:
c     jnucl - Nucleation rate (#part/cm^3/s).
c     ntot - Number of molecules in the critical cluster.
c     xstar - mol fraction of H2SO4 in the nucleated aerosol (mol).
c     dpnucl -  Nucleation diameter (nm).
      subroutine ssh_compute_binary_nucleation_kernel(rh,
     &     temp, natmp, jnucl,
     &     ntot, xstar, dpnucl)

      double precision rh, temp, na, natmp, xstar
      double precision jnucl, ntot, dpnucl

      double precision lnrh, lnrh2, lnrh3, fa(10)
      double precision t2, t3, lnna, lnna2, lnna3

      double precision tmin, tmax,
     &     rhmin, rhmax,
     &     namin, namax
      parameter(tmin = 230.15d0, tmax = 300.15d0,
     &     rhmin = 0.01d0, rhmax = 1.0d0,
     &     namin = 1.d04, namax = 1.d11)

C     test if parameterization is valid

      na = natmp
      if (na.gt.namax) then
         na = namax
      endif
      if (temp.lt.tmin.or.temp.gt.tmax) then
         dpnucl = 0.d0
         jnucl = 0.d0
         xstar = 0.d0
         ntot = 0.d0
      elseif(rh.gt.rhmax) then
         dpnucl = 0.d0
         jnucl = 0.d0
         xstar = 0.d0
         ntot = 0.d0
      elseif(na.lt.namin.or.rh.lt.rhmin) then
         dpnucl = 0.d0
         jnucl = 0.d0
         xstar = 0.d0
         ntot = 0.d0
      else
         lnna = dlog(na)
         lnna2 = lnna * lnna
         lnna3 = lnna2 * lnna
         lnrh = dlog(rh)
         lnrh2 = lnrh * lnrh
         lnrh3 = lnrh2 * lnrh
         t2 = temp * temp
         t3 = t2 * temp

c     Compute mole fraction of h2so4
c     in the critical cluster XSTAR.

         xstar = 0.740997d0
     &        - 2.66379d-3 * temp
     &        - 3.49998d-3 * lnna
     &        + 5.04022d-5 * temp * lnna
     &        + 2.01048d-3 * lnrh
     &        - 1.83289d-4 * temp * lnrh
     &        + 1.57407d-3 * lnrh2
     &        - 1.79059d-5 * temp * lnrh2
     &        + 1.84403d-4 * lnrh3
     &        - 1.50345d-6 * temp * lnrh3

c     Compute nucleation rate.

         call ssh_veahkamaki_coefficients(temp, t2, t3, xstar, fa)

         jnucl = dexp( fa(1)
     &        + fa(2) * lnrh
     &        + fa(3) * lnrh2
     &        + fa(4) * lnrh3
     &        + fa(5) * lnna
     &        + fa(6) * lnrh * lnna
     &        + fa(7) * lnrh2 * lnna
     &        + fa(8) * lnna2
     &        + fa(9) * lnrh * lnna2
     &        + fa(10) * lnna3 )

c     Compute total number of molecules
c     in the critical cluster ntot

         call ssh_veahkamaki_coefficients_number(temp, t2, t3, xstar,
     &        fa)

         ntot = dexp( fa(1)
     &        + fa(2) * lnrh
     &        + fa(3) * lnrh2
     &        + fa(4) * lnrh3
     &        + fa(5) * lnna
     &        + fa(6) * lnrh * lnna
     &        + fa(7) * lnrh2 * lnna
     &        + fa(8) * lnna2
     &        + fa(9) * lnrh * lnna2
     &        + fa(10) * lnna3 )

c     Compute cluster diameter in nm

         dpnucl = 2.d0 * dexp( - 1.6524245d0
     &        + 0.42316402d0 * xstar
     &        + 0.3346648d0 * dlog(ntot) )
      endif

      end


c     Function: veahkamaki_coefficients
c
c     Computes nucleation kernel coefficients.
c     Based on the parameterization of
c     Veahkamaki et al, JGR Atmosphere, vol 107, D22, p4622, 2002.
c
c     Parameters:
c     temperature - Temperature (K).
c     t2 - square of the temperature.
c     t3 - temperature to power three.
c     xstar - mole fraction  (mol).
c
c     Returns:
c     fa : nucleation kernel coefficients.
      subroutine ssh_veahkamaki_coefficients(temp, t2, t3, xstar, fa)

      double precision xstar, temp, t2, t3, fa(10)

      integer jj, nn
      double precision wa(50)

      data wa / 0.14309d0 , 2.21956d0 , - 2.73911d-2,
     &     7.22811d-5, 5.91822d0 , 0.117489d0,
     &     0.462532d0, - 1.18059d-2, 4.04196d-5,
     &     1.57963d01, - 0.215554d0, - 8.10269d-2,
     &     1.43581d-3, - 4.7758d-6 , - 2.91297d0,
     &     - 3.58856d0 , 4.9508d-2 , - 2.1382d-4,
     &     3.10801d-7, - 2.93333d-2, 1.14598d0,
     &     - 0.600796d0, 8.64245d-3, - 2.28947d-5,
     &     - 8.44985d0 , 2.15855d0 , 8.08121d-2,
     &     - 4.07382d-4, - 4.01957d-7, 0.721326d0,
     &     1.6241d0 , - 1.60106d-2, 3.77124d-5,
     &     3.21794d-8, - 1.13255d-2, 9.71682d0,
     &     - 0.115048d0, 1.57098d-4, 4.00914d-7,
     &     0.71186d0 , - 1.05611d0 , 9.03378d-3,
     &     - 1.98417d-5, 2.46048d-8, - 5.79087d-2,
     &     - 0.148712d0, 2.83508d-3, - 9.24619d-6,
     &     5.00427d-9, - 1.27081d-2 /

      do jj = 1, 10
         nn = (jj - 1) * 5

         fa(jj) = wa(nn + 1)
     &        + wa(nn + 2) * temp
     &        + wa(nn + 3) * t2
     &        + wa(nn + 4) * t3
     &        + wa(nn + 5) / xstar
      end do

      end


c     Function: veahkamaki_coefficients_number
c
c     Computes coefficients for the computation
c     of the number of nucleated molecules.
c     Based on the parameterization of
c     Veahkamaki et al, JGR Atmosphere, vol 107, D22, p4622, 2002.
c
c     Parameters:
c     temperature - Temperature (K).
c     t2 - square of the temperature.
c     t3 - temperature to power three.
c     xstar - mole fraction  (mol).
c
c     Returns:
c     fa : nucleation kernel coefficients for number.
      subroutine ssh_veahkamaki_coefficients_number(temp, t2, t3, xstar,
     &     fa)

      double precision xstar, temp, t2, t3, fa(10)

      integer jj, nn
      double precision wa(50)

      data wa / - 2.95413d-3 , - 9.76834d-2 , 1.02485d-3,
     &     - 2.18646d-6 , - 0.101717d0 , - 2.05064d-3,
     &     - 7.58504d-3 , 1.92654d-4 , - 6.7043d-7,
     &     - 0.255774d0 , 3.22308d-3 , 8.52637d-4,
     &     - 1.54757d-5 , 5.66661d-8 , 3.38444d-2,
     &     4.74323d-2 , - 6.25104d-4 , 2.65066d-6,
     &     - 3.67471d-9 , - 2.67251d-4 , - 1.25211d-2,
     &     5.80655d-3 , - 1.01674d-4 , 2.88195d-7,
     &     9.42243d-2 , - 3.8546d-2 , - 6.72316d-4,
     &     2.60288d-6 , 1.19416d-8 , - 8.51515d-3,
     &     - 1.83749d-2 , 1.722072d-4, - 3.71766d-7,
     &     - 5.14875d-10, 2.6866d-4 , - 6.19974d-2,
     &     9.06958d-4 , - 9.11728d-7 , - 5.36796d-9,
     &     - 7.74234d-3 , 1.21827d-2 , - 1.0665d-4,
     &     2.5346d-7 , - 3.63519d-10, 6.10065d-4,
     &     3.20184d-4 , - 1.74762d-5 , 6.06504d-8,
     &     - 1.42177d-11, 1.35751d-4 /

      do jj = 1, 10
         nn = (jj - 1) * 5

         fa(jj) = wa(nn + 1)
     &        + wa(nn + 2) * temp
     &        + wa(nn + 3) * t2
     &        + wa(nn + 4) * t3
     &        + wa(nn + 5) / xstar
      end do

      end


c     Function: na_threshold_veahkamaki.
c
c     Computes the critical H2SO4 concentration for
c     nucleation rate.
c     Based on the parameterization of
c     Veahkamaki et al, JGR Atmosphere, vol 107, D22, p4622, 2002.
c     Modified:
c     2004 : optimization (Kathleen Fahey).
c     2005/3/23: cleaning (Bruno Sportisse).
c
c     Parameters:
c     rh - Relative humidity 0< <1 ().
c     temperature - Temperature (K).
c
c     Returns:
c     nanucl : gas threshold H2SO4 concentration (#molec/cm^3).
      subroutine ssh_na_threshold_veahkamaki(rh, temp, nanucl)

      double precision rh, temp, nanucl

      double precision lnrh, invtemp

      lnrh = dlog(rh)
      invtemp = 1.d0 / temp

c     Compute threshold h2so4 concentration that produces the nucleation
c     rate 1.#part.cm-3.

      nanucl = dexp( - 2.79243d02
     &     + 1.17300d+01 * rh
     &     + ( 2.27009d+04
     &     - 1.08864d+03 * rh ) * invtemp
     &     + ( 1.14436d+00
     &     - 3.02331d-02 * rh
     &     - 1.30254d-03 * temp ) * temp
     &     + ( - 6.38697d+00
     &     + 8.54980d+02 * invtemp
     &     + 8.79662d-03 * temp ) * lnrh )

      end
