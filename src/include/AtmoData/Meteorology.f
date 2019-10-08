C-----------------------------------------------------------------------
C     Copyright (C) 2019 CEREA (ENPC) - INERIS
C     SSH-aerosol is distributed under the GNU General Public License v3
C-----------------------------------------------------------------------
      
c     File: Meteorology.f

c     Function: compute_relative_humidity
c     Computes relative humidity from specific humidity.
c
c     Parameters:
c     specific_humidity - Specific humidity (kg/kg).
c     temperature - Temperature (K).
c     pressure - Pressure (Pa).
c
c     Returns:
c     relative_humidity - relative humidity.
      subroutine compute_relative_humidity(sh,
     $     temp, pres, rh)

      double precision sh
      double precision temp
      double precision pres

      double precision rh

      double precision pressure_sat

      pressure_sat = 611.2d0 * dexp(17.67d0 * (temp - 273.15d0) /
     $     (temp - 29.65d0))

      rh = sh * pres / ( (0.62197d0
     $     * (1.d0 - sh) + sh)
     $     * pressure_sat)

      end
