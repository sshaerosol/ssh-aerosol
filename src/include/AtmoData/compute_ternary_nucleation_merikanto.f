C-----------------------------------------------------------------------
C     Copyright (C) 2019 CEREA (ENPC) - INERIS
C     SSH-aerosol is distributed under the GNU General Public License v3
C-----------------------------------------------------------------------

c     Function: compute_ternary_nucleation
c     
c     Computes ternary nucleation rate.
c     Parameterization of the nucleation rate fit
c     in Merikanto et al. 2007, corrigendum 2009
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
      subroutine ssh_compute_ternary_nucleation_merikanto(rhtmp,
     $     temptmp, natmp, mrtmp, jnucl, ntoth2so4, ntotnh3, dpnucl)

      double precision a0(20), a1(20), a2(20), a3(20)

      data a0 / -358.2337705052991,-980.923146020468,1200.472096232311,
     $     -14.833042158178936,-4.39129415725234D6, 4.905527742256349,
     $     -231375.56676032578,75061.15281456841, -3180.5610833308,
     $     -100.21645273730675, 5599.912337254629, 2.360931724951942D6,
     $     16597.75554295064,-89.38961120336789, -629.7882041830943,
     $     -732006.8180571689, 40751.075322248245, -1911.0303773001353,
     $     2.792313345723013,3.1712136610383244/
      data a1 /4.8630382337426985,10.054155220444462,-17.37107890065621,
     $     0.2932631303555295, 56383.93843154586, -0.05463019231872484,
     $     2919.2852552424706,-931.8802278173565, 39.08268568672095,
     $     0.977886555834732,-70.70896612937771,-29752.130254319443,
     $     -175.2365504237746,1.153344219304926, 7.772806552631709,
     $     9100.06398573816,-501.66977622013934,23.6903969622286,
     $     -0.03422552111802899,-0.037822330602328806/
      data a2 /-0.02175548069741675,-0.03306644502023841,
     $     0.08170681335921742,-0.0016497524241142845,
     $     -239.835990963361,0.00020258394697064567,
     $     -12.286497122264588,3.863266220840964,-0.16048521066690752,
     $     -0.0030511783284506377,0.2978801613269466,125.04965118142027,
     $     0.6033215603167458,-0.004954549700267233,
     $     -0.031974053936299256,
     $     -37.771091915932004,2.063469732254135,-0.09807872005428583,
     $     0.00014019195277521142,0.0001500555743561457/
      data a3 /0.00003212869941055865,0.000034274041225891804,
     $     -0.00012534476159729881,2.844074805239367D-6,
     $     0.33765136625580167,-2.502406532869512D-7,
     $     0.017249301826661612,-0.005349472062284983,
     $     0.00022031380023793877,2.967320346100855D-6,
     $     -0.00041866525019504,-0.1752996881934318,
     $     -0.0006731787599587544,7.096309866238719D-6,
     $     0.00004383764128775082,0.05235455395566905,
     $     -0.002836873785758324, 0.00013564560238552576,
     $     -1.9201227328396297D-7,-1.9828365865570703D-7/

      double precision rhtmp, temptmp
      double precision rh, temp, natmp, na, mrtmp, mr
      double precision jnucl, ntoth2so4, ntotnh3, dpnucl

      double precision ssh_fa2
      external ssh_fa2
      double precision lnrh, lnmr2, lnmr, lnj, lnj2
      double precision t2, t3, lnna, lnna2, lnmr3

      double precision c2,c3,t,Tonset

      double precision tmin, tmax,
     $     rhmin, rhmax,
     $     namin, namax, mrmin, mrmax
      parameter(tmin = 235.d0, tmax = 295.d0,
     $     rhmin = 0.05d0, rhmax = 0.95d0,
     $     namin = 5.d04, namax = 1.d09,
     s     mrmin = 0.1d0, mrmax = 1000.d0)

c     test if parameterization is valid.

      temp = dmin1(temptmp, tmax)
      rh = rhtmp
      na = dmin1(natmp, namax)
      mr = dmin1(mrtmp, mrmax)
      if (temp.lt.tmin.or.temp.gt.tmax
     $     .or.rh.lt.rhmin.or.rh.gt.rhmax
     $     .or.na.lt.namin.or.mr.lt.mrmin) then
         dpnucl = 0.d0
         jnucl = 0.d0
         ntoth2so4 = 0.d0
         ntotnh3 = 0.d0

      else
         lnna = dlog(na)
         lnna2 = lnna * lnna
         lnmr = dlog(mr)
         lnmr2 = lnmr * lnmr
         Tonset = 143.6 + 1.01789 * rh + 10.1964 * lnna - 0.184988*lnna2
     $        - 17.1618 * lnmr + 109.9247 * lnmr/lnna + 
     $        0.773412 * lnna * lnmr - 0.155764 * lnmr2
         if(temp.LT.Tonset) then
            dpnucl = 0.d0
            jnucl = 0.d0
            ntoth2so4 = 0.d0
            ntotnh3 = 0.d0
         else
            if(na.GT.namax) write(*,*) 'nucleation - sulfate too high',
     $           na,namax
!     na = dmin1(na, namax)
            mr = dmin1(mr, mrmax)
            lnrh = dlog(rh)
            lnmr3 = lnmr2 * lnmr
            t2 = temp * temp
            t3 = t2 * temp

            lnj = ( -12.86185d0
     $           + ssh_fa2(1, temp, t2, t3, a0, a1, a2, a3) *rh +
     $           ssh_fa2(2, temp, t2, t3, a0, a1, a2, a3) * lnrh +
     $           ssh_fa2(3, temp, t2, t3, a0, a1, a2, a3) * lnna +
     $           ssh_fa2(4, temp, t2, t3, a0, a1, a2, a3) * lnna2 +
     $           ssh_fa2(5, temp, t2, t3, a0, a1, a2, a3) / lnna2 +
     $           ssh_fa2(6, temp, t2, t3, a0, a1, a2, a3) * mr +
     $           ssh_fa2(7, temp, t2, t3, a0, a1, a2, a3) * lnmr +
     $           ssh_fa2(8, temp, t2, t3, a0, a1, a2, a3) * lnmr2 +
     $           ssh_fa2(9, temp, t2, t3, a0, a1, a2, a3) * lnmr3 +
     $           ssh_fa2(10, temp, t2, t3, a0, a1, a2, a3) * rh * lnmr +
     $           ssh_fa2(11, temp, t2, t3, a0, a1, a2, a3) * lnmr* lnna+
     $           ssh_fa2(12, temp, t2, t3, a0, a1, a2, a3) * lnmr/ lnna+
     $           ssh_fa2(13, temp, t2, t3, a0, a1, a2, a3) * lnrh /lnna+
     $           ssh_fa2(14, temp, t2, t3, a0, a1, a2, a3) * lnrh*lnmr +
     $           ssh_fa2(15, temp,t2,t3,a0,a1, a2, a3)*rh/mr/mr/mr/lnna+
     $           ssh_fa2(16, temp, t2, t3, a0, a1, a2, a3)/lnna * lnmr2+
     $           ssh_fa2(17, temp, t2, t3, a0, a1, a2, a3)/lnna * lnmr3+
     $           ssh_fa2(18, temp, t2, t3, a0, a1, a2, a3) * lnna*lnmr2+
     $           ssh_fa2(19, temp, t2, t3, a0, a1, a2, a3)*lnmr3* lnna2+
     $           ssh_fa2(20, temp, t2, t3, a0, a1, a2, a3) * lnrh*lnmr3)

            jnucl = dexp(lnj)

c$$$  C  Implementation as in the paper of Merikanto CCCCCCCCCCCCCCCCCCC
c$$$  c3 = mr
c$$$  c2 = na
c$$$  t = temp
c$$$  lnj = -12.861848898625231 + 4.905527742256349*c3 -  
c$$$  $  358.2337705052991*rh - 
c$$$  $  0.05463019231872484*c3*t + 4.8630382337426985*rh*t + 
c$$$  $  0.00020258394697064567*c3*t**2 - 0.02175548069741675*rh*t**2 -
c$$$  $  2.502406532869512e-7*c3*t**3 + 0.00003212869941055865*rh*t**3 - 
c$$$  $  4.39129415725234e6/Log(c2)**2+(56383.93843154586*t)/Log(c2)**2-
c$$$  $  (239.835990963361*t**2)/Log(c2)**2 + 
c$$$  $  (0.33765136625580167*t**3)/Log(c2)**2 - 
c$$$  $  (629.7882041830943*rh)/(c3**3*Log(c2)) + 
c$$$  $(7.772806552631709*rh*t)/(c3**3*Log(c2)) - 
c$$$  $  (0.031974053936299256*rh*t**2)/(c3**3*Log(c2)) + 
c$$$  $ (0.00004383764128775082*rh*t**3)/(c3**3*Log(c2)) + 
c$$$  $ 1200.472096232311*Log(c2) - 17.37107890065621*t*Log(c2) + 
c$$$  $  0.08170681335921742*t**2*Log(c2) - 
c$$$  $ 0.00012534476159729881*t**3*Log(c2) - 
c$$$  $14.833042158178936*Log(c2)**2 + 0.2932631303555295*t*Log(c2)**2 - 
c$$$  $ 0.0016497524241142845*t**2*Log(c2)**2 + 
c$$$  $2.844074805239367e-6*t**3*Log(c2)**2-231375.56676032578*Log(c3) - 
c$$$  $100.21645273730675*rh*Log(c3)+2919.2852552424706*t*Log(c3) + 
c$$$  $0.977886555834732*rh*t*Log(c3)-12.286497122264588*t**2*Log(c3) - 
c$$$  $0.0030511783284506377*rh*t**2*Log(c3) + 
c$$$  $0.017249301826661612*t**3*Log(c3) + 
c$$$  $2.967320346100855e-6*rh*t**3*Log(c3) + 
c$$$  $(2.360931724951942e6*Log(c3))/Log(c2) - 
c$$$  $(29752.130254319443*t*Log(c3))/Log(c2) + 
c$$$  $(125.04965118142027*t**2*Log(c3))/Log(c2) - 
c$$$  $(0.1752996881934318*t**3*Log(c3))/Log(c2) +
c$$$  $5599.912337254629*Log(c2)*Log(c3) - 
c$$$  $70.70896612937771*t*Log(c2)*Log(c3) + 
c$$$  $0.2978801613269466*t**2*Log(c2)*Log(c3) - 
c$$$  $0.00041866525019504*t**3*Log(c2)*Log(c3) + 
c$$$  $75061.15281456841*Log(c3)**2 - 
c$$$  $931.8802278173565*t*Log(c3)**2+3.863266220840964*t**2*Log(c3)**2 -
c$$$  $0.005349472062284983*t**3*Log(c3)**2 -
c$$$  $(732006.8180571689*Log(c3)**2)/Log(c2) + 
c$$$  $(9100.06398573816*t*Log(c3)**2)/Log(c2) - 
c$$$  $(37.771091915932004*t**2*Log(c3)**2)/Log(c2) + 
c$$$  $(0.05235455395566905*t**3*Log(c3)**2)/Log(c2) - 
c$$$  $1911.0303773001353*Log(c2)*Log(c3)**2 + 
c$$$  $23.6903969622286*t*Log(c2)*Log(c3)**2 - 
c$$$  $0.09807872005428583*t**2*Log(c2)*Log(c3)**2 + 
c$$$  $0.00013564560238552576*t**3*Log(c2)*Log(c3)**2 - 
c$$$  $3180.5610833308*Log(c3)**3 + 39.08268568672095*t*Log(c3)**3 -
c$$$  $0.16048521066690752*t**2*Log(c3)**3 + 
c$$$  $0.00022031380023793877*t**3*Log(c3)**3 + 
c$$$  $(40751.075322248245*Log(c3)**3)/Log(c2) - 
c$$$  $(501.66977622013934*t*Log(c3)**3)/Log(c2) + 
c$$$  $(2.063469732254135*t**2*Log(c3)**3)/Log(c2) - 
c$$$  $(0.002836873785758324*t**3*Log(c3)**3)/Log(c2) + 
c$$$  $2.792313345723013*Log(c2)**2*Log(c3)**3 - 
c$$$  $0.03422552111802899*t*Log(c2)**2*Log(c3)**3 + 
c$$$  $0.00014019195277521142*t**2*Log(c2)**2*Log(c3)**3 - 
c$$$  $1.9201227328396297e-7*t**3*Log(c2)**2*Log(c3)**3 - 
c$$$  $980.923146020468*Log(rh) + 10.054155220444462*t*Log(rh) - 
c$$$  $0.03306644502023841*t**2*Log(rh) + 
c$$$  $0.000034274041225891804*t**3*Log(rh) + 
c$$$  $(16597.75554295064*Log(rh))/Log(c2) - 
c$$$  $(175.2365504237746*t*Log(rh))/Log(c2) + 
c$$$  $(0.6033215603167458*t**2*Log(rh))/Log(c2) - 
c$$$  $(0.0006731787599587544*t**3*Log(rh))/Log(c2) - 
c$$$  $89.38961120336789*Log(c3)*Log(rh) + 
c$$$  $1.153344219304926*t*Log(c3)*Log(rh) - 
c$$$  $0.004954549700267233*t**2*Log(c3)*Log(rh) + 
c$$$  $7.096309866238719e-6*t**3*Log(c3)*Log(rh) + 
c$$$  $3.1712136610383244*Log(c3)**3*Log(rh) - 
c$$$  $0.037822330602328806*t*Log(c3)**3*Log(rh) + 
c$$$  $0.0001500555743561457*t**2*Log(c3)**3*Log(rh) - 
c$$$  $1.9828365865570703e-7*t**3*Log(c3)**3*Log(rh)
c$$$  $         jnucl = dexp(lnj)

            if(jnucl.lt.1.d-5) then
               jnucl = 0.d0
               ntoth2so4 = 0.0
               ntotnh3 = 0.0
            else
c     compute total number of molecules in the critical cluster: this is
c     ntot.
               if(jnucl.gt.1.d6) then
                  jnucl = 1.d6
                  lnj = log(jnucl)
               endif
               lnj2 = lnj * lnj
               ntoth2so4 = -4.71542 + 0.134364*temp - 0.000471847 * t2
     $      -2.56401 * lnna + 0.0113533 * temp* lnna + 0.00108019*lnna2 
     $      + 0.517137* lnmr - 0.00278825 * temp * lnmr + 0.806697*lnmr2
     $              - 0.00318491 * temp * lnmr2 - 0.0995118 * lnmr3
     $              + 0.000400728 * temp * lnmr3 + 1.32765 * lnj 
     $              - 0.00616765*temp*lnj - 0.110614 * lnmr * lnj 
     $              + 0.000436758 * temp* lnmr * lnj
     $              + 0.000916366 * lnj2

               ntotnh3 = 71.2007 - 0.840960 * temp + 0.00248030 * t2 
     $      + 2.77986 * lnna - 0.0147502 * temp * lnna + 0.0122645*lnna2
     $              - 2.00993 * lnmr + 0.00868912 * temp * lnmr 
     $              - 0.00914118 * lnna * lnmr + 0.137412 * lnmr2 
     $              - 0.000625323 * temp * lnmr2 + 0.0000937733 * lnmr3 
     $              + 0.520297 * lnj - 0.00241987 * temp * lnj
     $         + 0.0791639 * lnmr * lnj - 0.000302159 * temp * lnmr* lnj
     $              + 0.00469770 * lnj2

               if((ntoth2so4.le.0.).OR.(ntotnh3.le.0.)) then
                  jnucl = 0.d0
                  ntoth2so4 = 0.0
                  ntotnh3 = 0.0
               endif
            endif
c     compute cluster diameter in nm.

            dpnucl = 1.1

         endif
      endif

      end

c-------------------------------------------------

      function ssh_fa2(i, temp, t2, t3, a0, a1, a2, a3)

      double precision ssh_fa2, temp, t2, t3
      double precision a0(20), a1(20), a2(20), a3(20)

      ssh_fa2 = a0(i) + a1(i) * temp + a2(i) * t2 + a3(i) * t3

      end
