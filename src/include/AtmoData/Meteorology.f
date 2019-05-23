C-----------------------------------------------------------------------
C     Copyright (C) 2006-2007, ENPC - INRIA - EDF R&D
C     Author(s): Vivien Mallet
C
C     This file is part of AtmoData library, a tool for data processing
C     in atmospheric sciences.
C
C     AtmoData is developed in the INRIA - ENPC joint project-team CLIME
C     and in the ENPC - EDF R&D joint laboratory CEREA.
C
C     AtmoData is free software; you can redistribute it and/or modify
C     it under the terms of the GNU General Public License as published
C     by the Free Software Foundation; either version 2 of the License,
C     or (at your option) any later version.
C
C     AtmoData is distributed in the hope that it will be useful, but
C     WITHOUT ANY WARRANTY; without even the implied warranty of
C     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
C     General Public License for more details.
C
C     For more information, visit the AtmoData home page:
C          http://cerea.enpc.fr/polyphemus/atmodata.html
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
