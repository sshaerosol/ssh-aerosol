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


#ifndef ATMODATA_FILE_KZ_HXX

namespace AtmoData
{

  //! Returns the minimum of two elements.
  /*!
    \param x first number.
    \param y second number.
    \return The minimum of {x, y}.
  */
  template <class T, class T0>
  inline T min(T x, T0 y)
  {

    return x < y ? x : y;

  }

  //! Returns the maximum of two elements.
  /*!
    \param x first number.
    \param y second number.
    \return The maximum of {x, y}.
  */
  template <class T, class T0>
  inline T max(T x, T0 y)
  {

    return x > y ? x : y;

  }

  template<class TU, class TV, class TTp, class T, class TG>
  void ComputeLouisKz(Data<TU, 4, TG>& U, Data<TV, 4, TG>& V,
                      Data<TTp, 4, TG>& Tp, Data<T, 4, TG>& Kz,
                      T L0 = T(100), T B = T(5),
                      T C = T(5), T D = T(5), T z0 = T(1),
                      T a = T(0.115), T b = T(0.175), T delta_z0 = T(0.01),
                      T Ka = T(0.4));

  template<class TU, class TTp, class T, class TG>
  void ComputeLMO(const Data<TU, 3, TG>& FrictionModule,
                  const Data<TTp, 3, TG>& SurfacePotentialTemperature,
                  const Data<TTp, 4, TG>& PotentialTemperature,
                  const Data<T, 3, TG>& SensibleHeat,
                  const Data<T, 3, TG>& Evaporation,
                  Data<T, 3, TG>& LMO, T Ka = T(0.4));

  template<class TT, class TTp, class TU, class T, class TG>
  void ComputePBLH_TM(const Data<TT, 3, TG>& SurfaceTemperature,
                      const Data<TTp, 3, TG>& SurfacePotentialTemperature,
                      const Data<TTp, 4, TG>& PotentialTemperature,
                      const Data<TU, 3, TG>& FrictionModule,
                      const Data<TU, 4, TG>& WindModule,
                      const Data<T, 3, TG>& SensibleHeat,
                      const Data<T, 3, TG>& LMO,
                      const Grid<TG>& GridZ_interf,
                      Data<T, 3, TG>& BoundaryHeight,
                      T SBL = T(0.1), T Ric = T(0.21), T C = T(6.5),
                      T Ka = T(0.4));

  template<class T, class TG>
  void ComputePBLH_Richardson(const Data<T, 4, TG>& Richardson,
                              const Grid<TG>& GridZ_interf,
                              Data<T, 3, TG>& BoundaryHeight,
                              T Ric = T(0.21));

  template<class TU, class TT, class T, class TG>
  void ComputeTM_Kz(const Data<TT, 3, TG>& SurfaceTemperature,
                    const Data<TU, 3, TG>& FrictionModule,
                    const Data<T, 3, TG>& SensibleHeat,
                    const Data<T, 3, TG>& LMO,
                    const Data<T, 3, TG>& BoundaryHeight,
                    Data<T, 4, TG>& Kz,
                    bool TM_stable = true, T SBL = T(0.1), T p = T(2.0),
                    T Ka = T(0.4));

}  // namespace AtmoData.

#define ATMODATA_FILE_KZ_HXX
#endif
