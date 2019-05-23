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


#ifndef ATMODATA_FILE_POLAIR_CXX

#include "Polair.hxx"

namespace AtmoData
{


  //! Transforms the zonal wind for Polair.
  /*!
    Formula: ZonalWind = ZonalWind / cos(latitude).
    \param ZonalWind zonal wind.
    \note Coordinates associated with ZonalWind must be in degrees.
  */
  template<class T, class TG>
  void TransformZonalWind(Data<T, 4, TG>& ZonalWind)
  {

    int h, i, j, k;

    const T pi(3.14159265358979323846264);
    const T ratio = pi / 180.;

    int Nx = ZonalWind.GetLength(3);
    int Ny = ZonalWind.GetLength(2);
    int Nz = ZonalWind.GetLength(1);
    int Nt = ZonalWind.GetLength(0);

    for (h = 0; h < Nt; h++)
      for (k = 0; k < Nz; k++)
        for (j = 0; j < Ny; j++)
          for (i = 0; i < Nx; i++)
            ZonalWind(h, k, j, i) /= cos(ZonalWind[2].Value(h, k, j, i)
                                         * ratio);

  }


  //! Transforms the meridional wind for Polair.
  /*!
    Formula: MeridionalWind = MeridionalWind * cos(latitude).
    \param MeridionalWind meridional wind.
    \note Coordinates associated with MeridionalWind must be in degrees.
  */
  template<class T, class TG>
  void TransformMeridionalWind(Data<T, 4, TG>& MeridionalWind)
  {

    int h, i, j, k;

    const T pi(3.14159265358979323846264);
    const T ratio = pi / 180.;

    int Nx = MeridionalWind.GetLength(3);
    int Ny = MeridionalWind.GetLength(2);
    int Nz = MeridionalWind.GetLength(1);
    int Nt = MeridionalWind.GetLength(0);

    for (h = 0; h < Nt; h++)
      for (k = 0; k < Nz; k++)
        for (j = 0; j < Ny; j++)
          for (i = 0; i < Nx; i++)
            MeridionalWind(h, k, j, i)
              *= cos(MeridionalWind[2].Value(h, k, j, i) * ratio);

  }


}  // namespace AtmoData.

#define ATMODATA_FILE_POLAIR_CXX
#endif
