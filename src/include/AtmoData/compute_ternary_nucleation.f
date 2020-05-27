C-----------------------------------------------------------------------
C     Copyright (C) 2019 CEREA (ENPC) - INERIS
C     SSH-aerosol is distributed under the GNU General Public License v3
C-----------------------------------------------------------------------

c     Function: compute_ternary_nucleation
c
c     Computes ternary nucleation rate.
c     Parameterization of the nucleation rate fit
c     in Napari et al. JGR 2002, vol 107, D19, 4381
c
c     Parameters:
c     rh - Relative humidity 0< <1 ().
c     temperature - Temperature (K).
c     natmp - Gas-phase H2SO4 concentration (#molec/cm^3).
c     mrtmp - Mixing ratio of NH3 (ppt).
c
c     Returns:
c     jnucl - Nucleation rate (#part/cm^3/s).
c     ntoth2so4 - Number of molecules of H2SO4 in the critical cluster.
c     ntotnh3 - Number of molecules of NH3 in the critical cluster.
c     dpnucl -  Nucleation diameter (nm).
      subroutine ssh_compute_ternary_nucleation(rh, temp, natmp, mrtmp,
     $     jnucl, ntoth2so4, ntotnh3, dpnucl)

      double precision a0(20), a1(20), a2(20), a3(20)

      data a0 / - 0.355297d0, 3.13735d0, 19.0359d0,
     $     1.07605d0, 6.0916d0, 0.31176d0,
     $     - 0.0200738d0 , 0.165536d0 , 6.52645d0,
     $     3.68024d0, - 0.066514d0, 0.65874d0,
     $     0.0599321d0, - 0.732731d0, 0.728429d0 ,
     $     41.3016d0, - 0.160336d0, 8.57868d0,
     $     0.0530167d0, - 2.32736d0 /
      data a1 / - 33.8449d0, - 0.772861d0, - 0.170957d0,
     $     1.48932d0, - 1.25378d0, 1.64009d0,
     $     - 0.752115d0, 3.26623d0, - 0.258002d0,
     $     - 0.204098d0, - 7.82382d0, 0.190542d0,
     $     5.96475d0, - 0.0184179d0, 3.64736d0,
     $     - 0.35752d0, 0.00889881d0, - 0.112358d0,
     $     - 1.98815d0, 0.0234646d0 /
      data a2 / 0.34536d0, 0.00561204d0, 0.000479808d0,
     $     - 0.00796052d0, 0.00939836d0, - 0.00343852d0,
     $     0.00525813d0, - 0.0489703d0, 0.00143456d0,
     $     0.00106259d0, 0.0122938d0, - 0.00165718d0,
     $     - 0.0362432d0, 0.000147186d0, - 0.027422d0,
     $     0.000904383d0, - 0.0000539514d0, 0.000472626d0,
     $     0.0157827d0, - 0.000076519d0 /
      data a3 / - 0.000824007d0, - 9.74576d-06, - 4.14699d-07,
     $     7.61229d-06, - 0.0000174927d0, - 0.0000109753d0,
     $     - 8.98038d-06, 0.000146967d0, - 2.02036d-06,
     $     - 1.2656d-06, 0.0000618554d0, 3.41744d-06,
     $     0.0000493337d0, - 2.37711d-07, 0.0000493478d0,
     $     - 5.73788d-07, 8.39522d-08, - 6.48365d-07,
     $     - 0.0000293564d0, 8.0459d-08 /

      double precision rh, temp, natmp, na, mrtmp, mr
      double precision jnucl, ntoth2so4, ntotnh3, dpnucl

      double precision ssh_fa
      external ssh_fa
      double precision lnrh, lnmr2, lnmr, lnj, lnj2
      double precision t2, t3, lnna, lnna2

      double precision tmin, tmax,
     $     rhmin, rhmax,
     $     namin, namax, mrmin, mrmax
      parameter(tmin = 240.d0, tmax = 300.d0,
     $     rhmin = 0.05d0, rhmax = 0.95d0,
     $     namin = 1.d04, namax = 1.d09,
     s     mrmin = 0.1d0, mrmax = 100.d0)

c     test if parameterization is valid.

      na = natmp
      mr = mrtmp
      if (temp.lt.tmin.or.temp.gt.tmax
     $     .or.rh.lt.rhmin.or.rh.gt.rhmax
     $     .or.na.lt.namin.or.mr.lt.mrmin) then
         dpnucl = 0.d0
         jnucl = 0.d0
         ntoth2so4 = 0.d0
         ntotnh3 = 0.d0

      else
         na = dmin1(na, namax)
         mr = dmin1(mr, mrmax)

         lnna = dlog(na)
         lnna2 = lnna * lnna
         lnrh = dlog(rh)
         lnmr = dlog(mr)
         lnmr2 = lnmr * lnmr
         t2 = temp * temp
         t3 = t2 * temp

         lnj = ( - 84.7551d0
     $        + ssh_fa(1, temp, t2, t3, a0, a1, a2, a3) / lnna +
     $        ssh_fa(2, temp, t2, t3, a0, a1, a2, a3) * lnna +
     $        ssh_fa(3, temp, t2, t3, a0, a1, a2, a3) * lnna2 +
     $        ssh_fa(4, temp, t2, t3, a0, a1, a2, a3) * lnmr +
     $        ssh_fa(5, temp, t2, t3, a0, a1, a2, a3) * lnmr2 +
     $        ssh_fa(6, temp, t2, t3, a0, a1, a2, a3) * rh +
     $        ssh_fa(7, temp, t2, t3, a0, a1, a2, a3) * lnrh +
     $        ssh_fa(8, temp, t2, t3, a0, a1, a2, a3) * lnmr / lnna +
     $        ssh_fa(9, temp, t2, t3, a0, a1, a2, a3) * lnna * lnmr +
     $        ssh_fa(10, temp, t2, t3, a0, a1, a2, a3) * rh * lnna +
     $        ssh_fa(11, temp, t2, t3, a0, a1, a2, a3) * rh / lnna +
     $        ssh_fa(12, temp, t2, t3, a0, a1, a2, a3) * rh * lnmr +
     $        ssh_fa(13, temp, t2, t3, a0, a1, a2, a3) * lnrh / lnna +
     $        ssh_fa(14, temp, t2, t3, a0, a1, a2, a3) * lnrh * lnmr +
     $        ssh_fa(15, temp, t2, t3, a0, a1, a2, a3) * lnmr2 / lnna +
     $        ssh_fa(16, temp, t2, t3, a0, a1, a2, a3) * lnna * lnmr2 +
     $        ssh_fa(17, temp, t2, t3, a0, a1, a2, a3) * lnna2 * lnmr +
     $        ssh_fa(18, temp, t2, t3, a0, a1, a2, a3) * rh * lnmr2 +
     $        ssh_fa(19, temp, t2, t3, a0, a1, a2, a3) * rh * lnmr/lnna+
     $        ssh_fa(20, temp, t2, t3, a0, a1, a2, a3) * lnna2 * lnmr2)

         jnucl = dexp(lnj)

         if(jnucl.lt.1.d-5) then
             jnucl = 0.d0
             ntoth2so4 = 0.d0
             ntotnh3 = 0.d0
             dpnucl = 1.d0
         else 
           if(jnucl.gt.1.d6) then
              jnucl = 1.d6
              lnj = dlog(jnucl)
            endif
c     compute total number of molecules in the critical cluster: this is
c     ntot.

           lnj2 = lnj * lnj
         ntoth2so4 = 38.1645d0 + 0.774106d0 * lnj + 0.00298879d0 * lnj2
     s        - 0.357605d0 * temp - 0.00366358d0 * temp * lnj
     s        + 0.0008553d0 * temp ** 2

           ntotnh3 = 26.8982d0 + 0.682905d0 * lnj + 0.00357521d0 * lnj2
     s        - 0.265748d0 * temp - 0.00341895d0 * temp * lnj
     s        + 0.000673454d0 * temp ** 2

c     compute cluster diameter in nm.

           dpnucl = 2.d0 * (0.141027d0 - 0.00122625d0 * lnj
     s        - 7.82211d-06 * lnj2 - 0.00156727d0 * temp
     s        - 0.00003076d0 * temp * lnj
     s        + 0.0000108375d0 * temp ** 2)

        endif
      endif

      end

c-------------------------------------------------

      function ssh_fa(i, temp, t2, t3, a0, a1, a2, a3)

      double precision ssh_fa, temp, t2, t3
      double precision a0(20), a1(20), a2(20), a3(20)

      ssh_fa = a0(i) + a1(i) * temp + a2(i) * t2 + a3(i) * t3

      end
