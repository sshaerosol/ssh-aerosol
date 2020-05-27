C-----------------------------------------------------------------------
C     Copyright (C) 2020 CEREA (ENPC) - INERIS
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
      subroutine ssh_compute_relative_humidity(sh,
     $     temp, pres, rh)

      double precision sh
      double precision temp
      double precision pres

      double precision rh

      double precision pressure_sat

      call ssh_compute_psat_rh(rh, temp, pres, pressure_sat, rh)

      end

c     Function: compute_psat
c     Computes vapor saturation pressure
c          from temperature.
c
c     Parameters:
c     temp - Temperature (K).
c
c     Returns:
c     psat - saturation vapor pressure (Pa)

      subroutine ssh_compute_psat(temp, psat)

      double precision temp

      double precision psat

      psat = 611.2d0 * dexp(17.67d0 * (temp - 273.15d0) /
     $     (temp - 29.65d0))

      end

c     Function: compute_psat_rh
c     Computes vapor saturation pressure
c          and relative humidity from specific humidity.
c
c     Parameters:
c     sh - Specific humidity (kg/kg).
c     temp - Temperature (K).
c     pres - Pressure (Pa).
c
c     Returns:
c     psat - saturation vapor pressure (Pa)
c     rh - relative humidity.

      subroutine ssh_compute_psat_rh(sh, temp, pres, psat, rh)

      double precision sh
      double precision temp
      double precision pres

      double precision psat
      double precision rh

      call ssh_compute_psat(temp, psat)

      rh = sh * pres
     $   / ( (0.62197d0 * (1.d0 - sh) + sh) * psat )

      end

c     Function: compute_psat_sh
c     Computes vapor saturation pressure
c          and specific humidity from relative humidity.
c
c     Parameters:
c     sh - Specific humidity (kg/kg).
c     temp - Temperature (K).
c     pres - Pressure (Pa).
c
c     Returns:
c     psat - saturation vapor pressure (Pa)
c     rh - relative humidity.

      subroutine ssh_compute_psat_sh(rh, temp, pres, psat, sh)

      double precision rh
      double precision temp
      double precision pres

      double precision psat
      double precision sh

      call ssh_compute_psat(temp, psat)

      sh = 1.d0 / (   1.d0
     $              + pres / (psat * rh * 0.62197d0)
     $              - 1.d0 / 0.62197d0)

      end
