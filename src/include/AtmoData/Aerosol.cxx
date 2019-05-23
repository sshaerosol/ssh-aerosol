// Copyright (C) 2007, ENPC - INRIA - EDF R&D
// Author(s): Vivien Mallet, Elsa Real
//
// This file is part of AtmoData library, a tool for data processing in
// atmospheric sciences.
//
// AtmoData is developed in the INRIA - ENPC joint project-team CLIME and in
// the ENPC - EDF R&D joint laboratory CEREA.
//
// AtmoData is free software; you can redistribute it and/or modify it under
// the terms of the GNU General Public License as published by the Free
// Software Foundation; either version 2 of the License, or (at your option)
// any later version.
//
// AtmoData is distributed in the hope that it will be useful, but WITHOUT ANY
// WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
// FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
// details.
//
// For more information, visit the AtmoData home page:
//      http://cerea.enpc.fr/polyphemus/atmodata.html


#ifndef ATMODATA_FILE_AEROSOL_CXX

#include "Aerosol.hxx"

#include "Talos.hxx"
using namespace Talos;

#include <cmath>

namespace AtmoData
{

  template<class T>
  T wpower(T a, T b)
  {
    return exp(b * log(a));
  }

  //! Computes the aerosol wet diameter following the Hanel formula.
  /*!
    \param dry_radius aerosol dry radius (in m).
    \param relative_humidity relative humidity.
    \return The aerosol wet radius (in m).
  */
  template<class T>
  T compute_Hanel_diameter(T dry_radius, T relative_humidity)
  {
    // epsilon = 0.25 for organics and 0.285 for sulfate.
    const T epsilon = 0.25;
    return dry_radius * exp(-epsilon * log(1. - relative_humidity));
  }


  //! Computes the aerosol wet diameter following the Gerber formula.
  /*!
    \param dry_diameter aerosol dry diameter (in m).
    \param relative_humidity relative humidity.
    \param temperature temperature (in K).
    \return The aerosol wet radius (in m).
  */
  template<class T>
  T compute_Gerber_diameter(T dry_diameter, T relative_humidity,
                            T temperature)
  {
    // 'dry_diameter' and 'wet_diameter' unit is \mu m.
    const T C1(0.4989352162271429),
      C2(0.3026183900844475e1),
      C3(0.5372215625062934e-12),
      C4(-0.1371059101078550e1),
      C5(0.3942463621284677e-02);
    const T temperature_ref(298.0);

    T dry_radius = (dry_diameter / 2.) * 1.e-4;
    T aa = C3 * (1. + C5 * (temperature_ref - temperature));
    T wet_radius = wpower(C1 * wpower(dry_radius, C2) /
                          abs(aa * wpower(dry_radius, C4)
                              - log(relative_humidity))
                          + wpower<double>(T(dry_radius), T(3.)), T(1.) / 3.);
    return 2. * wet_radius * 1.e4;
  }


  //! Computes the aerosol wet diameter from the aerosol water content.
  /*!
    \param dry_diameter aerosol dry diameter (in m).
    \param dry_aerosol_concentration total dry aerosol concentration (in microg/m^3).
    \param water_concentration water concentration in aerosols (in microg/m^3).
    \return The aerosol wet radius (in m).
  */
  template<class T>
  T compute_wet_diameter_from_water_content(T dry_diameter,
                                            T dry_aerosol_concentration,
                                            T water_concentration)
  {
    // Unit of densities is g.m^{-3}.
    const T dry_aerosol_density = 1.4e6;
    const T water_density = 1.0e6;
    const T pi = 3.14159265358979323846264;

    T global_aerosol_concentration =
      dry_aerosol_concentration + water_concentration;

    T aerosol_number = 6. * dry_aerosol_concentration /
      (pi * wpower(dry_diameter, T(3.)) * dry_aerosol_density);

    T global_aerosol_density = (dry_aerosol_density * dry_aerosol_concentration
                                + water_density * water_concentration)
      / global_aerosol_concentration;

    return wpower(6. * global_aerosol_concentration
                  / (pi * global_aerosol_density * aerosol_number), T(1.) / 3.);
  }

}  // namespace AtmoData.

#define ATMODATA_FILE_AEROSOL_CXX
#endif
