// Copyright (C) 2003-2007, ENPC - INRIA - EDF R&D
// Author(s): Vivien Mallet
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


#ifndef ATMODATA_FILE_DEPOSITION_CXX

#include "Deposition.hxx"

namespace AtmoData
{


  //! Computes the bulk surface resistance according to Wesely (1989) and
  //! its revision in Walmsley and Wesely (1996).
  /*!
    \param surface_temperature surface temperature (K).
    \param solar_radiation solar radiation (W.m^{-2}).
    \param species species name, useful only if species=="O3"
    or species=="SO2".
    \param reactivity reactivity factor for the species.
    \param diffusivity molecular diffusivity for the species in air.
    \param Henry Henry constant for the species.
    \param Ri minimum bulk canopy stomatal resistance for water vapour.
    \param Rlu resistance of the outer surfaces in the upper canopy.
    \param Rac in-canopy aerodynamic resistance.
    \param RgsS soil resistance for SO2.
    \param RgsO soil resistance for O3.
    \param RclS resistance due to exposed surface in the lower canopy,
    for SO2.
    \param RclO resistance due to exposed surface in the lower canopy, for O3.
    \param limit (optional) upper limit for resistances. Default: 1.e10.
    \param D_H2O (optional) molecular diffusivity for water vapor.
    Default: 0.25.
    \return the bulk surface resistance.
  */
  template <class T>
  inline T ComputeWesely(T surface_temperature, T solar_radiation,
                         string species, T reactivity, T diffusivity, T Henry,
                         T Ri, T Rlu, T Rac, T RgsS, T RgsO, T RclS, T RclO,
                         T limit, T D_H2O)
  {

    T Rsx, Rmx, Rlux, Rdc, Rclx, Rgx;

    surface_temperature -= 273.15;

    // Stomatal resistance.
    if ((surface_temperature < 0.01) || (surface_temperature > 39.99))
      Rsx = 1.e30;
    else
      Rsx = Ri * (1. + 400. / ((solar_radiation + 0.1)
                               * (solar_radiation + 0.1)))
        * (400. / (surface_temperature * (40. - surface_temperature)));
    Rsx *= D_H2O / diffusivity;

    // Mesophyll resistance.
    Rmx = 1. / (Henry / 3000. + 100. * reactivity);

    // Cuticle resistance.
    Rlux = Rlu / (1.e-5 * Henry + reactivity);

    // Buoyant-convection resistance.
    Rdc = 100. * (1. + 1000. / (solar_radiation + 10.));

    // Rclx.
    if (species == "O3")
      Rclx = RclO;
    else if (species == "SO2")
      Rclx = RclS;
    else
      Rclx = 1. / (1.e-5 * Henry / RclS + reactivity / RclO);

    // Rgx.
    if (species == "O3")
      Rgx = RgsO;
    else if (species == "SO2")
      Rgx = RgsS;
    else
      Rgx = 1. / (1.e-5 * Henry / RgsS + reactivity / RgsO);

    // Limit.
    if (Rsx > limit) Rsx = limit;
    if (Rmx > limit) Rmx = limit;
    if (Rlux > limit) Rlux = limit;
    if (Rdc > limit) Rdc = limit;
    if (Rclx > limit) Rclx = limit;
    if (Rgx > limit) Rgx = limit;

    // Rc.
    return  1. / (1. / (Rsx + Rmx) + 1. / Rlux
                  + 1. / (Rdc + Rclx) + 1. / (Rac + Rgx));

  }


}  // namespace AtmoData.

#define ATMODATA_FILE_DEPOSITION_CXX
#endif
