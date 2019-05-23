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


#ifndef ATMODATA_FILE_METEOROLOGY_CXX

#include "Meteorology.hxx"

namespace AtmoData
{


  //! Computes the Richardson number.
  /*!
    Winds may be provided in two ways. The first option is to provide winds on
    interfaces (along x for the zonal wind, along y for the meridional wind).
    The second option is simply to provide winds at nodes (i.e. where the
    potential temperature and the Richardson number are defined).
    \param ZonalWind zonal wind.
    \param MeridionalWind meridional wind.
    \param PotentialTemperature potential temperature.
    \param Richardson (output) Richardson number.
    \param wind_threshold (optional) minimum of the wind shear. Default:0.001.
  */
  template<class TU, class TV, class TTp, class T, class TG>
  void ComputeRichardson(Data<TU, 4, TG>& ZonalWind,
                         Data<TV, 4, TG>& MeridionalWind,
                         Data<TTp, 4, TG>& PotentialTemperature,
                         Data<T, 4, TG>& Richardson, T wind_threshold)
  {

    int h, i, j, k;

    int Nx = Richardson.GetLength(3);
    int Ny = Richardson.GetLength(2);
    int Nz = Richardson.GetLength(1);
    int Nt = Richardson.GetLength(0);

    Grid<TG>& Levels = Richardson[1];
    Grid<TG>& MeridionalWindLevels = MeridionalWind[1];
    Grid<TG>& ZonalWindLevels = ZonalWind[1];

    const T g(9.81);
    T dudz, dvdz, dwinddz;
    int level;

    if (ZonalWind.GetLength(3) == Nx + 1
        && MeridionalWind.GetLength(2) == Ny + 1)
      for (h = 0; h < Nt; h++)
        for (k = 0; k < Nz; k++)
          for (j = 0; j < Ny; j++)
            for (i = 0; i < Nx; i++)
              {

                level = k == 0 ? 1 : k;

                dudz = 0.5
                  * ((ZonalWind(h, level, j, i + 1)
                      - ZonalWind(h, level - 1, j, i + 1))
                     / (ZonalWindLevels.Value(h, level, j, i + 1)
                        - ZonalWindLevels.Value(h, level - 1, j, i + 1))
                     + (ZonalWind(h, level, j, i)
                        - ZonalWind(h, level - 1, j, i))
                     / (ZonalWindLevels.Value(h, level, j, i)
                        - ZonalWindLevels.Value(h, level - 1, j, i)));
                dvdz = 0.5
                  * ((MeridionalWind(h, level, j + 1, i)
                      - MeridionalWind(h, level - 1, j + 1, i))
                     / (MeridionalWindLevels.Value(h, level, j + 1, i)
                        - MeridionalWindLevels.Value(h, level - 1, j + 1, i))
                     + (MeridionalWind(h, level, j, i)
                        - MeridionalWind(h, level - 1, j, i))
                     / (MeridionalWindLevels.Value(h, level, j, i)
                        - MeridionalWindLevels.Value(h, level - 1, j, i)));
                dwinddz = max(sqrt(dudz * dudz + dvdz * dvdz), wind_threshold);
                Richardson(h, k, j, i) =
                  g * (PotentialTemperature(h, level, j, i)
                       - PotentialTemperature(h, level - 1, j, i))
                  / (dwinddz * dwinddz
                     * PotentialTemperature(h, level - 1, j, i)
                     * (Levels.Value(h, level, j, i)
                        - Levels.Value(h, level - 1, j, i)));

              }
    else
      for (h = 0; h < Nt; h++)
        for (k = 0; k < Nz; k++)
          for (j = 0; j < Ny; j++)
            for (i = 0; i < Nx; i++)
              {

                level = k == 0 ? 1 : k;

                dudz = (ZonalWind(h, level, j, i)
                        - ZonalWind(h, level - 1, j, i))
                  / (ZonalWindLevels.Value(h, level, j, i)
                     - ZonalWindLevels.Value(h, level - 1, j, i));
                dvdz = (MeridionalWind(h, level, j, i)
                        - MeridionalWind(h, level - 1, j, i))
                  / (MeridionalWindLevels.Value(h, level, j, i)
                     - MeridionalWindLevels.Value(h, level - 1, j, i));
                dwinddz = max(sqrt(dudz * dudz + dvdz * dvdz), wind_threshold);
                Richardson(h, k, j, i) =
                  g * (PotentialTemperature(h, level, j, i)
                       - PotentialTemperature(h, level - 1, j, i))
                  / (dwinddz * dwinddz
                     * PotentialTemperature(h, level - 1, j, i)
                     * (Levels.Value(h, level, j, i)
                        - Levels.Value(h, level - 1, j, i)));

              }

  }


  //! Computes the Richardson number in the first layer.
  /*!
    Winds may be provided in two ways. The first option is to provide winds on
    interfaces (along x for the zonal wind, along y for the meridional wind).
    The second option is simply to provide winds at nodes (i.e. where the
    potential temperature and the Richardson number are defined).
    \param ZonalWind zonal wind.
    \param MeridionalWind meridional wind.
    \param PotentialTemperature potential temperature.
    \param Richardson (output) Richardson number in the first layer.
    \param wind_threshold (optional) minimum of the wind shear. Default:0.001.
  */
  template<class TU, class TV, class TTp, class T, class TG>
  void ComputeRichardson(Data<TU, 4, TG>& ZonalWind,
                         Data<TV, 4, TG>& MeridionalWind,
                         Data<TTp, 4, TG>& PotentialTemperature,
                         Data<T, 3, TG>& Richardson, T wind_threshold)
  {

    int h, i, j;

    int Nx = Richardson.GetLength(2);
    int Ny = Richardson.GetLength(1);
    int Nt = Richardson.GetLength(0);

    Grid<TG>& Levels = PotentialTemperature[1];
    Grid<TG>& MeridionalWindLevels = MeridionalWind[1];
    Grid<TG>& ZonalWindLevels = ZonalWind[1];

    const T g(9.81);
    T dudz, dvdz, dwinddz;

    if (ZonalWind.GetLength(3) == Nx + 1
        && MeridionalWind.GetLength(2) == Ny + 1)
      for (h = 0; h < Nt; h++)
        for (j = 0; j < Ny; j++)
          for (i = 0; i < Nx; i++)
            {

              dudz = 0.5
                * (ZonalWind(h, 0, j, i + 1)
                   / ZonalWindLevels.Value(h, 0, j, i + 1)
                   + ZonalWind(h, 0, j, i)
                   / ZonalWindLevels.Value(h, 0, j, i));
              dvdz = 0.5 *
                (MeridionalWind(h, 0, j + 1, i)
                 / MeridionalWindLevels.Value(h, 0, j + 1, i)
                 + MeridionalWind(h, 0, j, i)
                 / MeridionalWindLevels.Value(h, 0, j, i));
              dwinddz = max(sqrt(dudz * dudz + dvdz * dvdz), wind_threshold);
              Richardson(h, j, i) =
                g * (PotentialTemperature(h, 1, j, i)
                     - PotentialTemperature(h, 0, j, i))
                / (dwinddz * dwinddz * PotentialTemperature(h, 0, j, i)
                   * (Levels.Value(h, 1, j, i) - Levels.Value(h, 0, j, i)));

            }
    else
      for (h = 0; h < Nt; h++)
        for (j = 0; j < Ny; j++)
          for (i = 0; i < Nx; i++)
            {

              dudz = ZonalWind(h, 0, j, i)
                / ZonalWindLevels.Value(h, 0, j, i);
              dvdz = MeridionalWind(h, 0, j, i)
                / MeridionalWindLevels.Value(h, 0, j, i);
              dwinddz = max(sqrt(dudz * dudz + dvdz * dvdz), wind_threshold);
              Richardson(h, j, i) =
                g * (PotentialTemperature(h, 1, j, i)
                     - PotentialTemperature(h, 0, j, i))
                / (dwinddz * dwinddz * PotentialTemperature(h, 0, j, i)
                   * (Levels.Value(h, 1, j, i) - Levels.Value(h, 0, j, i)));

            }

  }


  //! Computes the surface Richardson number.
  /*!
    Computes the surface Richardson number.
    \param WindModule wind module in the first layer.
    \param SurfacePotentialTemperature surface potential temperature.
    \param PotentialTemperature potential temperature.
    \param SurfaceRichardson (output) surface Richardson number.
    \param wind_threshold (optional) minimum of the wind shear.
    Default: 0.001.
  */
  template<class TU, class TTp, class T, class TG>
  void ComputeRichardson(Data<TU, 3, TG>& WindModule,
                         Data<TTp, 3, TG>& SurfacePotentialTemperature,
                         Data<TTp, 4, TG>& PotentialTemperature,
                         Data<T, 3, TG>& SurfaceRichardson, T wind_threshold)
  {

    int h, i, j;

    int Nx = SurfaceRichardson.GetLength(2);
    int Ny = SurfaceRichardson.GetLength(1);
    int Nt = SurfaceRichardson.GetLength(0);

    Grid<TG>& Levels = PotentialTemperature[1];

    const T g(9.81);
    T wind;

    for (h = 0; h < Nt; h++)
      for (j = 0; j < Ny; j++)
        for (i = 0; i < Nx; i++)
          {
            wind = max(WindModule(h, j, i), wind_threshold);
            SurfaceRichardson(h, j, i) =
              2. * g * (PotentialTemperature(h, 0, j, i)
                        - SurfacePotentialTemperature(h, j, i))
              * Levels.Value(h, 0, j, i)
              / (wind * wind * (PotentialTemperature(h, 0, j, i)
                                + SurfacePotentialTemperature(h, j, i)));
          }

  }


  //! Computes the surface Richardson number.
  /*!
    Computes the surface Richardson number.
    \param Roughness roughness height.
    \param WindModule wind module in the first layer.
    \param SurfacePotentialTemperature surface potential temperature.
    \param PotentialTemperature potential temperature.
    \param SurfaceRichardson (output) surface Richardson number.
    \param wind_threshold (optional) minimum of the wind module.
    Default: 0.001.
  */
  template<class TR, class TU, class TTp, class T, class TG>
  void ComputeRichardson(Data<TR, 2, TG>& Roughness,
                         Data<TU, 3, TG>& WindModule,
                         Data<TTp, 3, TG>& SurfacePotentialTemperature,
                         Data<TTp, 4, TG>& PotentialTemperature,
                         Data<T, 3, TG>& SurfaceRichardson, T wind_threshold)
  {

    int h, i, j;

    int Nx = SurfaceRichardson.GetLength(2);
    int Ny = SurfaceRichardson.GetLength(1);
    int Nt = SurfaceRichardson.GetLength(0);

    Grid<TG>& Levels = PotentialTemperature[1];

    const T g(9.81);
    T wind;

    for (h = 0; h < Nt; h++)
      for (j = 0; j < Ny; j++)
        for (i = 0; i < Nx; i++)
          {
            wind = max(WindModule(h, j, i), wind_threshold);
            SurfaceRichardson(h, j, i) =
              2. * g * (PotentialTemperature(h, 0, j, i)
                        - SurfacePotentialTemperature(h, j, i))
              * (Levels.Value(h, 0, j, i) - Roughness(j, i))
              / (wind * wind * (PotentialTemperature(h, 0, j, i)
                                + SurfacePotentialTemperature(h, j, i))
                 * Levels.Value(h, 0, j, i) * Levels.Value(h, 0, j, i));
          }

  }


  //! Computes the potential temperature.
  /*!
    Formula: PotentialTemperature = Temperature * (Pressure / P0)^(-r/cp).
    \param Temperature temperature (or virtual temperature).
    \param Pressure pressure.
    \param PotentialTemperature (output) potential temperature.
    \param P0 (optional) standard pressure. Default: 101325 Pa.
    \param cp (optional) specific heat of dry air at constant pressure.
    Default: 1005 J.kg^{-1}.K^{-1}.
    \param r (optional) molar gas constant for air.
    Default: 287.0 J.kg^{-1}.K^{-1}.
  */
  template<class TT, class TP, class T, class TG>
  void ComputePotentialTemperature(Data<TT, 4, TG>& Temperature,
                                   Data<TP, 4, TG>& Pressure,
                                   Data<T, 4, TG>& PotentialTemperature,
                                   T P0, T cp, T r)
  {

    int h, i, j, k;

    int Nx = PotentialTemperature.GetLength(3);
    int Ny = PotentialTemperature.GetLength(2);
    int Nz = PotentialTemperature.GetLength(1);
    int Nt = PotentialTemperature.GetLength(0);

    T ratio = -r / cp;

    for (h = 0; h < Nt; h++)
      for (k = 0; k < Nz; k++)
        for (j = 0; j < Ny; j++)
          for (i = 0; i < Nx; i++)
            PotentialTemperature(h, k, j, i) =
              Temperature(h, k, j, i)
              * pow(Pressure(h, k, j, i) / P0, ratio);

  }


  //! Computes the potential temperature.
  /*!
    Formula: PotentialTemperature = Temperature * (Pressure / P0)^(-r/cp).
    \param Temperature temperature (or virtual temperature).
    \param Pressure pressure.
    \param PotentialTemperature (output) potential temperature.
    \param P0 (optional) standard pressure. Default: 101325 Pa.
    \param cp (optional) specific heat of dry air at constant pressure.
    Default: 1005 J.kg^{-1}.K^{-1}.
    \param r (optional) molar gas constant for air.
    Default: 287.0 J.kg^{-1}.K^{-1}.
  */
  template<class TT, class TP, class T, class TG>
  void ComputePotentialTemperature(Data<TT, 3, TG>& Temperature,
                                   Data<TP, 3, TG>& Pressure,
                                   Data<T, 3, TG>& PotentialTemperature,
                                   T P0, T cp, T r)
  {

    int h, i, j;

    int Nx = PotentialTemperature.GetLength(2);
    int Ny = PotentialTemperature.GetLength(1);
    int Nt = PotentialTemperature.GetLength(0);

    T ratio = -r / cp;

    for (h = 0; h < Nt; h++)
      for (j = 0; j < Ny; j++)
        for (i = 0; i < Nx; i++)
          PotentialTemperature(h, j, i) =
            Temperature(h, j, i) * pow(Pressure(h, j, i) / P0, ratio);

  }


  //! Computes the temperature.
  /*!
    Formula: Temperature = PotentialTemperature * (Pressure / P0)^(r/cp).
    \param PotentialTemperature potential temperature.
    \param Pressure pressure.
    \param Temperature (output) temperature (or virtual temperature).
    \param P0 (optional) standard pressure. Default: 101325 Pa.
    \param cp (optional) specific heat of dry air at constant pressure.
    Default: 1005 J.kg^{-1}.K^{-1}.
    \param r (optional) molar gas constant for air.
    Default: 287.0 J.kg^{-1}.K^{-1}.
  */
  template<class TT, class TP, class T, class TG>
  void ComputeTemperature(Data<TT, 4, TG>& PotentialTemperature,
                          Data<TP, 4, TG>& Pressure,
                          Data<T, 4, TG>& Temperature,
                          T P0, T cp, T r)
  {

    int h, i, j, k;

    int Nx = Temperature.GetLength(3);
    int Ny = Temperature.GetLength(2);
    int Nz = Temperature.GetLength(1);
    int Nt = Temperature.GetLength(0);

    T ratio = r / cp;

    for (h = 0; h < Nt; h++)
      for (k = 0; k < Nz; k++)
        for (j = 0; j < Ny; j++)
          for (i = 0; i < Nx; i++)
            Temperature(h, k, j, i) =
              PotentialTemperature(h, k, j, i)
              * pow(Pressure(h, k, j, i) / P0, ratio);

  }


  //! Computes the temperature.
  /*!
    Formula: Temperature = PotentialTemperature * (Pressure / P0)^(r/cp).
    \param PotentialTemperature potential temperature.
    \param Pressure pressure.
    \param Temperature (output) temperature (or virtual temperature).
    \param P0 (optional) standard pressure. Default: 101325 Pa.
    \param cp (optional) specific heat of dry air at constant pressure.
    Default: 1005 J.kg^{-1}.K^{-1}.
    \param r (optional) molar gas constant for air.
    Default: 287.0 J.kg^{-1}.K^{-1}.
  */
  template<class TT, class TP, class T, class TG>
  void ComputeTemperature(Data<TT, 3, TG>& PotentialTemperature,
                          Data<TP, 3, TG>& Pressure,
                          Data<T, 3, TG>& Temperature,
                          T P0, T cp, T r)
  {

    int h, i, j;

    int Nx = Temperature.GetLength(2);
    int Ny = Temperature.GetLength(1);
    int Nt = Temperature.GetLength(0);

    T ratio = r / cp;

    for (h = 0; h < Nt; h++)
      for (j = 0; j < Ny; j++)
        for (i = 0; i < Nx; i++)
          Temperature(h, j, i) = PotentialTemperature(h, j, i)
            * pow(Pressure(h, j, i) / P0, ratio);

  }


  //! Computes the saturation humidity.
  /*!
    \param Temperature temperature (or virtual temperature) (K).
    \param Pressure pressure (Pa).
    \param SaturationHumidity (output) saturation humidity (kg/kg).
  */
  template<class TT, class TP, class T, class TG>
  void ComputeSaturationHumidity(Data<TT, 3, TG>& Temperature,
                                 Data<TP, 3, TG>& Pressure,
                                 Data<T, 3, TG>& SaturationHumidity)
  {
    int h, j, i;

    int Nt(SaturationHumidity.GetLength(0));
    int Ny(SaturationHumidity.GetLength(1));
    int Nx(SaturationHumidity.GetLength(2));

    T P_sat;
    for (h = 0; h < Nt; h++)
      for (j = 0; j < Ny; j++)
        for (i = 0; i < Nx; i++)
          {
            P_sat = 611.2 * exp(17.67 * (Temperature(h, j, i) - 273.15)
                                / (Temperature(h, j, i) - 29.65));
            SaturationHumidity(h, j, i) = 0.622 * P_sat
              / (Pressure(h, j, i) - 0.378 * P_sat);
          }
  }


  //! Computes the saturation humidity.
  /*!
    \param Temperature temperature (or virtual temperature) (K).
    \param Pressure pressure (Pa).
    \param SaturationHumidity (output) saturation humidity (kg/kg).
  */
  template<class TT, class TP, class T, class TG>
  void ComputeSaturationHumidity(Data<TT, 4, TG>& Temperature,
                                 Data<TP, 4, TG>& Pressure,
                                 Data<T, 4, TG>& SaturationHumidity)
  {
    int h, k, j, i;

    int Nt(SaturationHumidity.GetLength(0));
    int Nz(SaturationHumidity.GetLength(1));
    int Ny(SaturationHumidity.GetLength(2));
    int Nx(SaturationHumidity.GetLength(3));

    T P_sat;
    for (h = 0; h < Nt; h++)
      for (k = 0; k < Nz; k++)
        for (j = 0; j < Ny; j++)
          for (i = 0; i < Nx; i++)
            {
              P_sat = 611.2 * exp(17.67 * (Temperature(h, k, j, i) - 273.15)
                                  / (Temperature(h, k, j, i) - 29.65));
              SaturationHumidity(h, k, j, i) = 0.622 * P_sat
                / (Pressure(h, k, j, i) - 0.378 * P_sat);
            }
  }


  //! Computes the relative humidity from the specific humidity.
  /*!
    \param SpecificHumidity specific humidity (kg/kg).
    \param Temperature temperature (or virtual temperature) (K).
    \param Pressure pressure (Pa).
    \param RelativeHumidity (output) relative humidity.
  */
  template<class TS, class TT, class TP, class T, class TG>
  void ComputeRelativeHumidity(Data<TS, 4, TG>& SpecificHumidity,
                               Data<TT, 4, TG>& Temperature,
                               Data<TP, 4, TG>& Pressure,
                               Data<T, 4, TG>& RelativeHumidity)
  {
    int h, k, j, i;
    int Nt(RelativeHumidity.GetLength(0));
    int Nz(RelativeHumidity.GetLength(1));
    int Ny(RelativeHumidity.GetLength(2));
    int Nx(RelativeHumidity.GetLength(3));

    T P_sat;
    for (h = 0; h < Nt; h++)
      for (k = 0; k < Nz; k++)
        for (j = 0; j < Ny; j++)
          for (i = 0; i < Nx; i++)
            {
              P_sat = 611.2 * exp(17.67 * (Temperature(h, k, j, i) - 273.15)
                                  / (Temperature(h, k, j, i) - 29.65));
              RelativeHumidity(h, k, j, i) = SpecificHumidity(h, k, j, i)
                * Pressure(h, k, j, i)
                / ((0.62197 * (1.0 - SpecificHumidity(h, k, j, i))
                    + SpecificHumidity(h, k, j, i)) * P_sat);
            }
  }


  //! Computes the relative humidity from the specific humidity.
  /*!
    \param SpecificHumidity specific humidity (kg/kg).
    \param Temperature temperature (or virtual temperature) (K).
    \param Pressure pressure (Pa).
    \param RelativeHumidity (output) relative humidity.
  */
  template<class TS, class TT, class TP, class T, class TG>
  void ComputeRelativeHumidity(Data<TS, 3, TG>& SpecificHumidity,
                               Data<TT, 3, TG>& Temperature,
                               Data<TP, 3, TG>& Pressure,
                               Data<T, 3, TG>& RelativeHumidity)
  {
    int k, j, i;
    int Nz(RelativeHumidity.GetLength(0));
    int Ny(RelativeHumidity.GetLength(1));
    int Nx(RelativeHumidity.GetLength(2));

    T P_sat;
    for (k = 0; k < Nz; k++)
      for (j = 0; j < Ny; j++)
        for (i = 0; i < Nx; i++)
          {
            P_sat = 611.2 * exp(17.67 * (Temperature(k, j, i) - 273.15)
                                / (Temperature(k, j, i) - 29.65));
            RelativeHumidity(k, j, i) = SpecificHumidity(k, j, i)
              * Pressure(k, j, i)
              / ((0.62197 * (1.0 - SpecificHumidity(k, j, i))
                  + SpecificHumidity(k, j, i)) * P_sat);
          }
  }


  //! Computes the surface specific humidity.
  /*!
    \param SpecificHumidity specific humidity (kg/kg).
    \param SaturationHumidity saturation humidity (kg/kg).
    \param SoilWater volumetric soil water content (m^3/m^3).
    \param LUC land use coverage in format Nc x Ny x Nx where Nc is the
    number of land use categories. LUC(c, j, i) is the relative surface
    (in [0, 1]) of the category c in the cell (j, i).
    \param sea_index index of the sea category in LUC.
    \param SurfaceHumidity (output) surface specific humidity (kg/kg).
    \param veg (optional) the vegetation proportion (in [0, 1]) on the ground.
    Default: 1.0.
    \param theta_cap (optional) soil moisture at field capacity (m^3/m^3).
    Default: 0.323 m^3/m^3.
  */
  template<class TH, class TS, class TW, class TL, class T, class TG>
  void ComputeSurfaceHumidity_diag(Data<TH, 4, TG>& SpecificHumidity,
                                   Data<TS, 3, TG>& SaturationHumidity,
                                   Data<TW, 3, TG>& SoilWater,
                                   Data<TL, 3, TG>& LUC, int sea_index,
                                   Data<T, 3, TG>& SurfaceHumidity,
                                   T veg, T theta_cap)
  {
    int h, j, i;

    int Nt(SurfaceHumidity.GetLength(0));
    int Ny(SurfaceHumidity.GetLength(1));
    int Nx(SurfaceHumidity.GetLength(2));

    T q_sat, alpha;
    for (h = 0; h < Nt; h++)
      for (j = 0; j < Ny; j++)
        for (i = 0; i < Nx; i++)
          {
            q_sat = SaturationHumidity(h, j, i);

            if (SoilWater(h, j, i) < theta_cap)
              // 1.963495 = pi / 1.6.
              alpha = 0.5 * (1. - cos(1.963495 * SoilWater(h, j, i)
                                      / theta_cap));
            else
              alpha = 1.0;

            SurfaceHumidity(h, j, i) = LUC(sea_index, j, i) * q_sat
              + (1.0 - LUC(sea_index, j, i))
              * (alpha * q_sat + veg * (1. - alpha)
                 * min(q_sat, SpecificHumidity(h, 0, j, i)));
          }
  }


  //! Computes the critical relative humidity.
  /*!
    Formula: CriticalRelativeHumidity = 1.0 - coeff0 * sig
    * (1.0 - sig) * (1.0 + (sig - 0.5) * coeff1) where
    sig = Pressure / SurfacePressure.
    \param SurfacePressure surface pressure (Pa).
    \param Pressure pressure (Pa).
    \param CriticalRelativeHumidity (output) critical relative humidity.
    \param coeff0 coefficient (see the formula). Default: 2.0.
    \param coeff1 coefficient (see the formula). Default: sqrt(3.).
  */
  template<class TS, class TP, class T, class TG>
  void ComputeCriticalRelativeHumidity(Data<TS, 3, TG>& SurfacePressure,
                                       Data<TP, 4, TG>& Pressure,
                                       Data<T, 4, TG>&
                                       CriticalRelativeHumidity,
                                       T coeff0, T coeff1)
  {
    int h, k, j, i;
    int Nt(CriticalRelativeHumidity.GetLength(0));
    int Nz(CriticalRelativeHumidity.GetLength(1));
    int Ny(CriticalRelativeHumidity.GetLength(2));
    int Nx(CriticalRelativeHumidity.GetLength(3));

    T sig;
    for (h = 0; h < Nt; h++)
      for (k = 0; k < Nz; k++)
        for (j = 0; j < Ny; j++)
          for (i = 0; i < Nx; i++)
            {
              sig = Pressure(h, k, j, i) / SurfacePressure(h, j, i);
              CriticalRelativeHumidity(h, k, j, i) = 1.0 - coeff0 * sig
                * (1.0 - sig) * (1.0 + (sig - 0.5) * coeff1);
            }
  }


  //! Computes the critival relative humidity.
  /*!
    Extended formula: CriticalRelativeHumidity = 1.0 - coeff0 * sig^a0
    * (1.0 - sig)^a1 * (1.0 + (sig - 0.5) * coeff1) where
    sig = Pressure / SurfacePressure.
    \param SurfacePressure surface pressure (Pa).
    \param Pressure pressure (Pa).
    \param CriticalRelativeHumidity (output) critical relative humidity.
    \param coeff0 coefficient (see the formula). Default: 1.1.
    \param coeff1 coefficient (see the formula). Default: sqrt(1.3).
    \param a0 exponent (see the formula). Default: 0.0.
    \param a1 exponent (see the formula). Default: 1.1.
  */
  template<class TS, class TP, class T, class TG>
  void ComputeCriticalRelativeHumidity_extended(Data<TS, 3, TG>&
                                                SurfacePressure,
                                                Data<TP, 4, TG>& Pressure,
                                                Data<T, 4, TG>&
                                                CriticalRelativeHumidity,
                                                T coeff0, T coeff1,
                                                T a0, T a1)
  {
    int h, k, j, i;
    int Nt(CriticalRelativeHumidity.GetLength(0));
    int Nz(CriticalRelativeHumidity.GetLength(1));
    int Ny(CriticalRelativeHumidity.GetLength(2));
    int Nx(CriticalRelativeHumidity.GetLength(3));

    T sig;
    for (h = 0; h < Nt; h++)
      for (k = 0; k < Nz; k++)
        for (j = 0; j < Ny; j++)
          for (i = 0; i < Nx; i++)
            {
              sig = Pressure(h, k, j, i) / SurfacePressure(h, j, i);
              CriticalRelativeHumidity(h, k, j, i) =
                1.0 - coeff0 * pow(sig, a0)
                * pow(T(1.) - sig, a1) * (1.0 + (sig - 0.5) * coeff1);
            }
  }


  //! Computes the critical relative humidity.
  /*!
    Formula: inside the boundary layer,
    CriticalRelativeHumidity = BL_CRH and, above the boundary layer,
    CriticalRelativeHumidity = 1.0 - coeff0 * sig
    * (1.0 - sig) * (1.0 + (sig - 0.5) * coeff1) where
    sig = Pressure / SurfacePressure.
    \param BoundaryLayerHeight boundary layer height (m).
    \param SurfacePressure surface pressure (Pa).
    \param Pressure pressure (Pa).
    \param CriticalRelativeHumidity (output) critical relative humidity.
    \param coeff0 coefficient (see the formula). Default: 2.0.
    \param coeff1 coefficient (see the formula). Default: sqrt(3.).
    \param BL_CRH critical relative humidity within the boundary layer.
    Default: 0.98.
  */
  template<class TB, class TS, class TP, class T, class TG>
  void ComputeCriticalRelativeHumidity(Data<TB, 3, TG>& BoundaryLayerHeight,
                                       Data<TS, 3, TG>& SurfacePressure,
                                       Data<TP, 4, TG>& Pressure,
                                       Data<T, 4, TG>&
                                       CriticalRelativeHumidity,
                                       T coeff0, T coeff1, T BL_CRH)
  {
    int h, k, j, i;
    int Nt(CriticalRelativeHumidity.GetLength(0));
    int Nz(CriticalRelativeHumidity.GetLength(1));
    int Ny(CriticalRelativeHumidity.GetLength(2));
    int Nx(CriticalRelativeHumidity.GetLength(3));

    T sig;
    for (h = 0; h < Nt; h++)
      for (k = 0; k < Nz; k++)
        for (j = 0; j < Ny; j++)
          for (i = 0; i < Nx; i++)
            if (CriticalRelativeHumidity[1].Value(h, k, j, i)
                < BoundaryLayerHeight(h, j, i))
              CriticalRelativeHumidity(h, k, j, i) = BL_CRH;
            else
              {
                sig = Pressure(h, k, j, i) / SurfacePressure(h, j, i);
                CriticalRelativeHumidity(h, k, j, i) = 1.0 - coeff0 * sig
                  * (1.0 - sig) * (1.0 + (sig - 0.5) * coeff1);
              }
  }


  //! Computes the critical relative humidity.
  /*!
    The relative humidity is set to CRH_0 if Pressure > P_0,
    to CRH_1 if P_0 >= Pressure > P_1 and to CRH_2 otherwise.
    \param Pressure pressure (Pa).
    \param CriticalRelativeHumidity (output) critical relative humidity.
    \param CRH_0 critical relative humidity in the first layer. Default: 0.75.
    \param CRH_1 critical relative humidity in the second layer.
    Default: 0.95.
    \param CRH_2 critical relative humidity in the third layer. Default: 0.95.
    \param P_0 first pressure limit. Default: 70 000 Pa.
    \param P_1 second pressure limit. Default: 40 000 Pa.
  */
  template<class TP, class T, class TG>
  void ComputeCriticalRelativeHumidity(Data<TP, 4, TG>& Pressure,
                                       Data<T, 4, TG>&
                                       CriticalRelativeHumidity,
                                       T CRH_0, T CRH_1, T CRH_2,
                                       T P_0, T P_1)
  {
    int h, k, j, i;
    int Nt(CriticalRelativeHumidity.GetLength(0));
    int Nz(CriticalRelativeHumidity.GetLength(1));
    int Ny(CriticalRelativeHumidity.GetLength(2));
    int Nx(CriticalRelativeHumidity.GetLength(3));

    for (h = 0; h < Nt; h++)
      for (k = 0; k < Nz; k++)
        for (j = 0; j < Ny; j++)
          for (i = 0; i < Nx; i++)
            if (Pressure(h, k, j, i) > P_0)
              CriticalRelativeHumidity(h, k, j, i) = CRH_0;
            else if (Pressure(h, k, j, i) > P_1)
              CriticalRelativeHumidity(h, k, j, i) = CRH_1;
            else
              CriticalRelativeHumidity(h, k, j, i) = CRH_2;
  }


  //! Computes the cloud fraction.
  /*!
    Formula: CloudFraction = [ (RelativeHumidity - CriticalRelativeHumidity)
    / (1.0 - CriticalRelativeHumidity) ]^2
    \param RelativeHumidity relative humidity.
    \param CriticalRelativeHumidity critical relative humidity.
    \param CloudFraction (output) cloud fraction.
  */
  template<class TR, class TC, class T, class TG>
  void ComputeCloudFraction(Data<TR, 4, TG>& RelativeHumidity,
                            Data<TC, 4, TG>& CriticalRelativeHumidity,
                            Data<T, 4, TG>& CloudFraction)
  {
    int h, k, j, i;
    int Nt(CloudFraction.GetLength(0));
    int Nz(CloudFraction.GetLength(1));
    int Ny(CloudFraction.GetLength(2));
    int Nx(CloudFraction.GetLength(3));

    T tmp, crh;
    for (h = 0; h < Nt; h++)
      for (k = 0; k < Nz; k++)
        for (j = 0; j < Ny; j++)
          for (i = 0; i < Nx; i++)
            {
              crh = CriticalRelativeHumidity(h, k, j, i);
              tmp = RelativeHumidity(h, k, j, i) - crh;
              if (tmp < 0. || crh == 1.)
                CloudFraction(h, k, j, i) = 0.;
              else
                {
                  tmp = tmp / (1. - crh);
                  CloudFraction(h, k, j, i) = tmp * tmp;
                }
            }
  }


  //! Computes the cloud fraction.
  /*!
    Formula: inside the boundary layer,
    CloudFraction = 0.34 * [ (RelativeHumidity - CriticalRelativeHumidity)
    / (1.0 - CriticalRelativeHumidity) ] and, above the boundary layer,
    CloudFraction = [ (RelativeHumidity - CriticalRelativeHumidity)
    / (1.0 - CriticalRelativeHumidity) ]^2
    \param BoundaryLayerHeight boundary layer height (m).
    \param RelativeHumidity relative humidity.
    \param CriticalRelativeHumidity critical relative humidity.
    \param CloudFraction (output) cloud fraction.
  */
  template<class TP, class TR, class TC, class T, class TG>
  void ComputeCloudFraction(Data<TP, 3, TG>& BoundaryLayerHeight,
                            Data<TR, 4, TG>& RelativeHumidity,
                            Data<TC, 4, TG>& CriticalRelativeHumidity,
                            Data<T, 4, TG>& CloudFraction)
  {
    int h, k, j, i;
    int Nt(CloudFraction.GetLength(0));
    int Nz(CloudFraction.GetLength(1));
    int Ny(CloudFraction.GetLength(2));
    int Nx(CloudFraction.GetLength(3));

    T tmp, crh;
    for (h = 0; h < Nt; h++)
      for (k = 0; k < Nz; k++)
        for (j = 0; j < Ny; j++)
          for (i = 0; i < Nx; i++)
            {
              crh = CriticalRelativeHumidity(h, k, j, i);
              tmp = RelativeHumidity(h, k, j, i) - crh;
              if (tmp < 0. || crh == 1.)
                CloudFraction(h, k, j, i) = 0.;
              else if (CloudFraction[1].Value(h, k, j, i)
                       < BoundaryLayerHeight(h, j, i))
                CloudFraction(h, k, j, i) = 0.34 * tmp / (1. - crh);
              else
                {
                  tmp = tmp / (1. - crh);
                  CloudFraction(h, k, j, i) = tmp * tmp;
                }
            }
  }


  //! Computes the module of a 2D-vectors field.
  /*!
    This function was initially dedicated to winds. In this case, zonal winds
    and meridional winds are provided and the module of the wind is computed
    (assuming that the vertical wind is zero). Winds may be provided in two
    ways. The first option is to provide winds on interfaces (along x for
    the zonal wind, along y for the meridional wind). The second option is
    simply to provide winds at nodes (i.e. where the module is defined).
    \param U first component of vectors.
    \param V second component of vectors.
    \param Module (output) module.
  */
  template<class TU, class TV, class T, class TG>
  void ComputeModule(Data<TU, 4, TG>& U, Data<TV, 4, TG>& V,
                     Data<T, 4, TG>& Module)
  {

    int h, i, j, k;

    int Nx = Module.GetLength(3);
    int Ny = Module.GetLength(2);
    int Nz = Module.GetLength(1);
    int Nt = Module.GetLength(0);

    T u, v;

    if ((U.GetLength(3) == Nx + 1) && (V.GetLength(2) == Ny + 1))
      for (h = 0; h < Nt; h++)
        for (k = 0; k < Nz; k++)
          for (j = 0; j < Ny; j++)
            for (i = 0; i < Nx; i++)
              {
                u = 0.5 * (U(h, k, j, i + 1) + U(h, k, j, i));
                v = 0.5 * (V(h, k, j + 1, i) + V(h, k, j, i));
                Module(h, k, j, i) = sqrt(u * u + v * v);
              }
    else
      for (h = 0; h < Nt; h++)
        for (k = 0; k < Nz; k++)
          for (j = 0; j < Ny; j++)
            for (i = 0; i < Nx; i++)
              {
                u = U(h, k, j, i);
                v = V(h, k, j, i);
                Module(h, k, j, i) = sqrt(u * u + v * v);
              }

  }


  //! Computes the module of a 2D-vectors field on the surface.
  /*!
    This function was initially dedicated to winds. In this case, zonal winds
    and meridional winds are provided and the module of the wind is computed
    (assuming that the vertical wind is zero). Winds may be provided in two
    ways. The first option is to provide winds on interfaces (along x
    for the zonal wind, along y for the meridional wind). The second option is
    simply to provide winds at nodes (i.e. where the module is defined).
    \param U first component of vectors.
    \param V second component of vectors.
    \param Module (output) surface module.
  */
  template<class TU, class TV, class T, class TG>
  void ComputeModule(Data<TU, 4, TG>& U, Data<TV, 4, TG>& V,
                     Data<T, 3, TG>& Module)
  {

    int h, i, j;

    int Nx = Module.GetLength(2);
    int Ny = Module.GetLength(1);
    int Nt = Module.GetLength(0);

    T u, v;

    if ((U.GetLength(3) == Nx + 1) && (V.GetLength(2) == Ny + 1))
      for (h = 0; h < Nt; h++)
        for (j = 0; j < Ny; j++)
          for (i = 0; i < Nx; i++)
            {
              u = 0.5 * (U(h, 0, j, i + 1) + U(h, 0, j, i));
              v = 0.5 * (V(h, 0, j + 1, i) + V(h, 0, j, i));
              Module(h, j, i) = sqrt(u * u + v * v);
            }
    else
      for (h = 0; h < Nt; h++)
        for (j = 0; j < Ny; j++)
          for (i = 0; i < Nx; i++)
            {
              u = U(h, 0, j, i);
              v = V(h, 0, j, i);
              Module(h, j, i) = sqrt(u * u + v * v);
            }

  }


  //! Computes the module of a 2D-vectors field on the surface.
  /*!
    This function was initially dedicated to winds. In this case, zonal winds
    and meridional winds are provided and the module of the wind is computed
    (assuming that the vertical wind is zero). Winds may be provided in two
    ways. The first option is to provide winds on interfaces (along x for
    the zonal wind, along y for the meridional wind). The second option is
    simply to provide winds at nodes (i.e. where the module is defined).
    \param U first component of vectors.
    \param V second component of vectors.
    \param Module (output) surface module.
  */
  template<class TU, class TV, class T, class TG>
  void ComputeModule(Data<TU, 3, TG>& U, Data<TV, 3, TG>& V,
                     Data<T, 3, TG>& Module)
  {

    int h, i, j;

    int Nx = Module.GetLength(2);
    int Ny = Module.GetLength(1);
    int Nt = Module.GetLength(0);

    T u, v;

    if ((U.GetLength(2) == Nx + 1) && (V.GetLength(1) == Ny + 1))
      for (h = 0; h < Nt; h++)
        for (j = 0; j < Ny; j++)
          for (i = 0; i < Nx; i++)
            {
              u = 0.5 * (U(h, j, i + 1) + U(h, j, i));
              v = 0.5 * (V(h, j + 1, i) + V(h, j, i));
              Module(h, j, i) = sqrt(u * u + v * v);
            }
    else
      for (h = 0; h < Nt; h++)
        for (j = 0; j < Ny; j++)
          for (i = 0; i < Nx; i++)
            {
              u = U(h, j, i);
              v = V(h, j, i);
              Module(h, j, i) = sqrt(u * u + v * v);
            }

  }


  //! Computes low, medium and high cloudiness, cloud base and top.
  /*!
    \param CloudFraction cloud fraction.
    \param Pressure pressure (Pa).
    \param GridZ_interf altitudes of interfaces (m).
    \param LowIndices (output) vertical indices of base and top of low clouds.
    \param MediumIndices (output) vertical indices of base and
    top of medium clouds.
    \param HighIndices (output) vertical indices of base and
    top of high clouds.
    \param LowCloudiness (output) low cloudiness.
    \param MediumCloudiness (output) medium cloudiness.
    \param HighCloudiness (output) high cloudiness.
    \param P_0 first pressure limit. Default: 80 000 Pa.
    \param P_1 second pressure limit. Default: 45 000 Pa.
    \note Dimensions of LowIndices, MediumIndices and HighIndices are
    Nt x Ny x Nx x 2. Along the last dimension, those arrays store the index
    of the cloud base and the index of the cloud top (in this order). Those
    indices are indices of interfaces. E.g., if LowIndices(t, y, x, 0) equals
    2 and  LowIndices(t, y, x, 1) equals 4, then a cloud lies in layers 2
    and 3.
  */
  template<class TC, class TP, class T, class TG>
  void ComputeCloudiness(Data<TC, 4, TG>& CloudFraction,
                         Data<TP, 4, TG>& Pressure,
                         Grid<TG>& GridZ_interf,
                         Data<int, 4>& LowIndices,
                         Data<int, 4>& MediumIndices,
                         Data<int, 4>& HighIndices,
                         Data<T, 3, TG>& LowCloudiness,
                         Data<T, 3, TG>& MediumCloudiness,
                         Data<T, 3, TG>& HighCloudiness,
                         T P_0, T P_1)
  {
    int h, k, j, i;
    int Nt(CloudFraction.GetLength(0));
    int Nz(CloudFraction.GetLength(1));
    int Ny(CloudFraction.GetLength(2));
    int Nx(CloudFraction.GetLength(3));

    LowCloudiness.SetZero();
    MediumCloudiness.SetZero();
    HighCloudiness.SetZero();

    TC cloud_max;
    bool above, below;
    int k_base, k_top, k_max;
    for (h = 0; h < Nt; h++)
      for (j = 0; j < Ny; j++)
        for (i = 0; i < Nx; i++)
          {

            /*** Low clouds ***/

            cloud_max = 0;
            // The first level is excluded.
            for (k = 1; k < Nz && Pressure(h, k, j, i) > P_0; k++)
              cloud_max = max(cloud_max, CloudFraction(h, k, j, i));
            below = true;
            above = false;
            k_base = 0;
            k_top = 0;
            for (k = 1; k < Nz && Pressure(h, k, j, i) > P_0 && !above; k++)
              {
                below = below && (CloudFraction(h, k, j, i) < 0.5 * cloud_max
                                  || CloudFraction(h, k, j, i) == 0);
                above = !below && CloudFraction(h, k, j, i) < 0.5 * cloud_max;
                if (!below && k_base == 0)
                  k_base = k;
                if (above)
                  k_top = k;
              }
            if (above)
              while (k < Nz && Pressure(h, k, j, i) > P_0)
                k++;
            k_max = k - 1;
            // Goes up to P_0.
            if (k_base > k_top)
              k_top = k_max + 1;
            LowIndices(h, j, i, 0) = k_base;
            LowIndices(h, j, i, 1) = k_top;
            k = k_base;
            // k_top == 0 means no cloud.
            while (k < k_top && k_top != 0)
              {
                LowCloudiness(h, j, i) += CloudFraction(h, k, j, i)
                  * (GridZ_interf.Value(h, k + 1, j, i)
                     - GridZ_interf.Value(h, k, j, i));
                k++;
              }
            if (k_top != 0)
              LowCloudiness(h, j, i) /=
                GridZ_interf.Value(h, k_top, j, i)
                - GridZ_interf.Value(h, k_base, j, i);

            /*** Medium clouds ***/

            cloud_max = 0;
            // Starts above low clouds.
            for (k = k_max + 1; k < Nz && Pressure(h, k, j, i) > P_1; k++)
              cloud_max = max(cloud_max, CloudFraction(h, k, j, i));
            below = true;
            above = false;
            k_base = 0;
            k_top = 0;
            for (k = k_max + 1; k < Nz && Pressure(h, k, j, i) > P_1
                   && !above; k++)
              {
                below = below && (CloudFraction(h, k, j, i) < 0.5 * cloud_max
                                  || CloudFraction(h, k, j, i) == 0);
                above = !below && CloudFraction(h, k, j, i) < 0.5 * cloud_max;
                if (!below && k_base == 0)
                  k_base = k;
                if (above)
                  k_top = k;
              }
            if (above)
              while (k < Nz && Pressure(h, k, j, i) > P_1)
                k++;
            k_max = k - 1;
            // Goes up to P_1.
            if (k_base > k_top)
              k_top = k_max + 1;
            MediumIndices(h, j, i, 0) = k_base;
            MediumIndices(h, j, i, 1) = k_top;
            k = k_base;
            // k_top == 0 means no cloud.
            while (k < k_top && k_top != 0)
              {
                MediumCloudiness(h, j, i) += CloudFraction(h, k, j, i)
                  * (GridZ_interf.Value(h, k + 1, j, i)
                     - GridZ_interf.Value(h, k, j, i));
                k++;
              }
            if (k_top != 0)
              MediumCloudiness(h, j, i) /=
                GridZ_interf.Value(h, k_top, j, i)
                - GridZ_interf.Value(h, k_base, j, i);

            /*** High clouds ***/

            cloud_max = 0;
            // Starts above low clouds.
            for (k = k_max + 1; k < Nz; k++)
              cloud_max = max(cloud_max, CloudFraction(h, k, j, i));
            below = true;
            above = false;
            k_base = 0;
            k_top = 0;
            for (k = k_max + 1; k < Nz && !above; k++)
              {
                below = below && (CloudFraction(h, k, j, i) < 0.5 * cloud_max
                                  || CloudFraction(h, k, j, i) == 0);
                above = !below && CloudFraction(h, k, j, i) < 0.5 * cloud_max;
                if (!below && k_base == 0)
                  k_base = k;
                if (above)
                  k_top = k;
              }
            k_max = k - 1;
            // Goes up to the top.
            if (k_base > k_top)
              k_top = k_max + 1;
            HighIndices(h, j, i, 0) = k_base;
            HighIndices(h, j, i, 1) = k_top;
            k = k_base;
            // k_top == 0 means no cloud.
            while (k < k_top && k_top != 0)
              {
                HighCloudiness(h, j, i) += CloudFraction(h, k, j, i)
                  * (GridZ_interf.Value(h, k + 1, j, i)
                     - GridZ_interf.Value(h, k, j, i));
                k++;
              }
            if (k_top != 0)
              HighCloudiness(h, j, i) /=
                GridZ_interf.Value(h, k_top, j, i)
                - GridZ_interf.Value(h, k_base, j, i);
          }
  }


  //! Computes the height of cloud basis.
  /*!
    \param Pressure pressure (Pa).
    \param RelativeHumidity relative humidity (kg/kg).
    \param CriticalRelativeHumidity function that returns the critical
    relative humidity as function of the altitude, the pressure
    and reference pressure.
    \param CloudBaseHeight (output) altitudes of cloud basis.
  */
  template < class TP, class TH,
             class T, class TG >
  void ComputeCloudBaseHeight(Data<TP, 4, TG>& Pressure,
                              Data<TH, 4, TG>& RelativeHumidity,
                              T(CriticalRelativeHumidity)(const T&, const T&,
                                                          const T&),
                              Data<T, 3, TG>& CloudBaseHeight)
  {

    int h, k, j, i;
    int Nt(CloudBaseHeight.GetLength(0));
    int Nz(Pressure.GetLength(1));
    int Ny(CloudBaseHeight.GetLength(1));
    int Nx(CloudBaseHeight.GetLength(2));

    // Index "0" and "1" refer to two contiguous levels.
    T rh0, rh1, rhc;

    T max_height = 2. * Pressure[1].Value(0, Nz - 1, 0, 0);
    CloudBaseHeight.Fill(max_height);

    for (h = 0; h < Nt; h++)
      for (j = 0; j < Ny; j++)
        for (i = 0; i < Nx; i++)
          {

            rh0 = RelativeHumidity(h, Nz - 1, j, i);

            k = 0;
            while ((k < Nz) && (CloudBaseHeight(h, j, i) == max_height))
              {

                rh1 = RelativeHumidity(h, k, j, i);

                // Critical relative humidity.
                rhc = CriticalRelativeHumidity(Pressure[1].Value(h, k, j, i),
                                               Pressure(h, k, j, i),
                                               Pressure(h, 0, j, i));

                if (rh1 >= rhc)  // Above a cloud.
                  CloudBaseHeight(h, j, i) = Pressure[1].Value(h, k, j, i);

                // For the next level.
                rh0 = rh1;

                k++;

              }

          }
  }


  //! Computes the height of cloud basis.
  /*!
    \param RelativeHumidity relative humidity (kg/kg).
    \param CriticalRelativeHumidity critical relative humidity.
    \param CloudBaseHeight (output) altitudes of cloud basis.
  */
  template <class TH, class TCRH, class T, class TG>
  void ComputeCloudBaseHeight(Data<TH, 4, TG>& RelativeHumidity,
                              Data<TCRH, 4, TG>& CriticalRelativeHumidity,
                              Data<T, 3, TG>& CloudBaseHeight)
  {

    int h, k, j, i;
    int Nt(CloudBaseHeight.GetLength(0));
    int Nz(RelativeHumidity.GetLength(1));
    int Ny(CloudBaseHeight.GetLength(1));
    int Nx(CloudBaseHeight.GetLength(2));

    // Index "0" and "1" refer to two contiguous levels.
    T rh0, rh1, rhc;

    T max_height = 2. * RelativeHumidity[1].Value(0, Nz - 1, 0, 0);
    CloudBaseHeight.Fill(max_height);

    for (h = 0; h < Nt; h++)
      for (j = 0; j < Ny; j++)
        for (i = 0; i < Nx; i++)
          {

            rh0 = RelativeHumidity(h, Nz - 1, j, i);

            k = 0;
            while ((k < Nz) && (CloudBaseHeight(h, j, i) == max_height))
              {

                rh1 = RelativeHumidity(h, k, j, i);

                // Critical relative humidity.
                rhc = CriticalRelativeHumidity(h, k, j, i);

                if (rh1 >= rhc)  // Above a cloud.
                  CloudBaseHeight(h, j, i) =
                    RelativeHumidity[1].Value(h, k, j, i);

                // For the next level.
                rh0 = rh1;

                k++;

              }

          }
  }


  //! Computes the height of cloud basis.
  /*!
    \param LowIndices vertical indices of base and top of low clouds.
    \param MediumIndices vertical indices of base and top of medium clouds.
    \param HighIndices vertical indices of base and top of high clouds.
    \param GridZ_interf altitudes of interfaces (m).
    \param CloudBaseHeight (output) altitudes of cloud basis.
    \note Dimensions of LowIndices, MediumIndices and HighIndices are
    Nt x Ny x Nx x 2. Along the last dimension, those arrays store the index
    of the cloud base and the index of the cloud top (in this order). Those
    indices are indices of interfaces. E.g., if LowIndices(t, y, x, 0) equals
    2 and  LowIndices(t, y, x, 1) equals 4, then a cloud lies in layers 2
    and 3.
  */
  template <class T, class TG>
  void ComputeCloudBaseHeight(Data<int, 4>& LowIndices,
                              Data<int, 4>& MediumIndices,
                              Data<int, 4>& HighIndices,
                              Grid<TG>& GridZ_interf,
                              Data<T, 3, TG>& CloudBaseHeight)
  {
    int h, j, i;
    int Nt(CloudBaseHeight.GetLength(0));
    int Ny(CloudBaseHeight.GetLength(1));
    int Nx(CloudBaseHeight.GetLength(2));

    CloudBaseHeight.SetZero();

    for (h = 0; h < Nt; h++)
      for (j = 0; j < Ny; j++)
        for (i = 0; i < Nx; i++)
          if (LowIndices(h, j, i, 0) != 0)
            CloudBaseHeight(h, j, i) =
              GridZ_interf.Value(h, LowIndices(h, j, i, 0), j, i);
          else if (MediumIndices(h, j, i, 0) != 0)
            CloudBaseHeight(h, j, i) =
              GridZ_interf.Value(h, MediumIndices(h, j, i, 0), j, i);
          else if (HighIndices(h, j, i, 0) != 0)
            CloudBaseHeight(h, j, i) =
              GridZ_interf.Value(h, HighIndices(h, j, i, 0), j, i);
  }


  //! Computes the height of cloud top (only for in-cloud scavenging purposes)
  /*!
    \param LowIndices vertical indices of base and top of low clouds.
    \param MediumIndices vertical indices of base and top of medium clouds.
    \param HighIndices vertical indices of base and top of high clouds.
    \param GridZ_interf altitudes of interfaces (m).
    \param CloudTopHeight (output) altitudes of cloud tops.
    \note Dimensions of LowIndices, MediumIndices and HighIndices are
    Nt x Ny x Nx x 2. Along the last dimension, those arrays store the index
    of the cloud base and the index of the cloud top (in this order). Those
    indices are indices of interfaces. E.g., if LowIndices(t, y, x, 0) equals
    2 and  LowIndices(t, y, x, 1) equals 4, then a cloud lies in layers 2
    and 3.
  */
  template <class T, class TG>
  void ComputeCloudTopHeight(Data<int, 4>& LowIndices,
                             Data<int, 4>& MediumIndices,
                             Data<int, 4>& HighIndices,
                             Grid<TG>& GridZ_interf,
                             Data<T, 3, TG>& CloudTopHeight)
  {
    int h, j, i;
    int Nt(CloudTopHeight.GetLength(0));
    int Ny(CloudTopHeight.GetLength(1));
    int Nx(CloudTopHeight.GetLength(2));

    CloudTopHeight.SetZero();

    for (h = 0; h < Nt; h++)
      for (j = 0; j < Ny; j++)
        for (i = 0; i < Nx; i++)
          if (LowIndices(h, j, i, 1) != 0)
            {
              if (LowIndices(h, j, i, 1) != MediumIndices(h, j, i, 0))
                CloudTopHeight(h, j, i) =
                  GridZ_interf.Value(h, LowIndices(h, j, i, 1), j, i);
              else if (MediumIndices(h, j, i, 1) != HighIndices(h, j, i, 0))
                CloudTopHeight(h, j, i) =
                  GridZ_interf.Value(h, MediumIndices(h, j, i, 1), j, i);
              else
                CloudTopHeight(h, j, i) =
                  GridZ_interf.Value(h, HighIndices(h, j, i, 1), j, i);
            }
          else if (MediumIndices(h, j, i, 1) != 0)
            {
              if (MediumIndices(h, j, i, 1) != HighIndices(h, j, i, 0))
                CloudTopHeight(h, j, i) =
                  GridZ_interf.Value(h, MediumIndices(h, j, i, 1), j, i);
              else
                CloudTopHeight(h, j, i) =
                  GridZ_interf.Value(h, HighIndices(h, j, i, 1), j, i);
            }
          else if (HighIndices(h, j, i, 1) != 0)
            CloudTopHeight(h, j, i) =
              GridZ_interf.Value(h, HighIndices(h, j, i, 1), j, i);
  }


  //! Computes the total cloudiness.
  /*!
    \param LowCloudiness low cloudiness in [0, 1].
    \param MediumCloudiness medium cloudiness in [0, 1].
    \param HighCloudiness high cloudiness in [0, 1].
    \param Cloudiness (output) total cloudiness in [0, 1].
    \return the total cloudiness in [0, 1].
  */
  template<class T, class TLC, class TC>
  void ComputeTotalCloudiness(Data<TLC, 3, T>& LowCloudiness,
                              Data<TLC, 3, T>& MediumCloudiness,
                              Data<TLC, 3, T>& HighCloudiness,
                              Data<TC, 3, T>& Cloudiness)
  {
    int h, j, i;
    int Nt(Cloudiness.GetLength(0));
    int Ny(Cloudiness.GetLength(1));
    int Nx(Cloudiness.GetLength(2));

    for (h = 0; h < Nt; h++)
      for (j = 0; j < Ny; j++)
        for (i = 0; i < Nx; i++)
          Cloudiness(h, j, i) = 1. - (1. - LowCloudiness(h, j, i))
            * (1. - MediumCloudiness(h, j, i))
            * (1. - HighCloudiness(h, j, i));
  }


  //! Computes the total cloudiness at a given point.
  /*!
    \param low_cloudiness low cloudiness in [0, 1].
    \param medium_cloudiness medium cloudiness in [0, 1].
    \param high_cloudiness high cloudiness in [0, 1].
    \return the total cloudiness in [0, 1].
  */
  template<class TLC, class TC>
  inline TC ComputeTotalCloudiness(TLC low_cloudiness, TLC medium_cloudiness,
                                   TLC high_cloudiness)
  {
    return 1. - (1. - low_cloudiness) * (1. - medium_cloudiness)
      * (1. - high_cloudiness);
  }


  //! Computes the Pasquill stability class (source: Turner 1969).
  /*!
    \param surface_wind wind speed at 10 meters (m/s).
    \param solar_radiation solar radiation (W.m^{-2}).
    \param cloudiness cloudiness in [0, 1].
    \param isday Boolean equal to true if it is daytime, false otherwise.
    \return the Pasquill stability class (A, B, C, D, E or F).
    \note the conditions of dual stability classes like A-B, B-C or C-D are
    considered as B, C and D respectively.
  */
  template<class T>
  string ComputePasquillStabilityClass(T surface_wind, T solar_radiation,
                                       T cloudiness, bool isday)
  {
    if (isday)
      if (surface_wind < 2.)
        if (solar_radiation > 700.)
          return "A";
        else
          return "B";
      else if (surface_wind < 3.)
        if (solar_radiation >= 350.)
          return "B";
        else
          return "C";
      else if (surface_wind < 5.)
        if (solar_radiation > 700.)
          return "B";
        else
          return "C";
      else if (surface_wind < 6.)
        if (solar_radiation > 700.)
          return "C";
        else
          return "D";
      else
        return "D";
    else if (surface_wind < 2.)
      return "F";
    else if (surface_wind < 3.)
      if (cloudiness > 0.5)
        return "E";
      else
        return "F";
    else if (surface_wind < 5.)
      if (cloudiness > 0.5)
        return "D";
      else
        return "E";
    else
      return "D";
  }


  //! Computes the pressure from the surface pressure.
  /*!

    For MACC files, computes with 4D surface pressure !

    Formula: Pressure_k = alpha_k * P0 + beta_k * SurfacePressure,
    where k is the level index.
    \param alpha coefficients.
    \param beta coefficients.
    \param SurfacePressure surface pressure.
    \param Pressure (output) pressure.
    \param P0 (optional) standard pressure. Default: 101325 Pa.
  */
  template < class Ta, class Tb, class TSP,
             class T, class TG >
  void Compute4DPressure(Data<Ta, 1, TG>& alpha, Data<Tb, 1, TG>& beta,
                         Data<TSP, 4, TG>& SurfacePressure,
                         Data<T, 4, TG>& Pressure, T P0)
  {

    int h, i, j, k;

    int Nx = Pressure.GetLength(3);
    int Ny = Pressure.GetLength(2);
    int Nz = Pressure.GetLength(1);
    int Nt = Pressure.GetLength(0);

    for (h = 0; h < Nt; h++)
      for (k = 0; k < Nz; k++)
        for (j = 0; j < Ny; j++)
          for (i = 0; i < Nx; i++)
            Pressure(h, k, j, i) = alpha(k)
              + beta(k) * SurfacePressure(h, 0, j, i);

  }


  //! Computes the pressure from the surface pressure.
  /*!
    Formula: Pressure_k = alpha_k * P0 + beta_k * SurfacePressure,
    where k is the level index.
    \param alpha coefficients.
    \param beta coefficients.
    \param SurfacePressure surface pressure.
    \param Pressure (output) pressure.
    \param P0 (optional) standard pressure. Default: 101325 Pa.
  */
  template < class Ta, class Tb, class TSP,
             class T, class TG >
  void ComputePressure(Data<Ta, 1, TG>& alpha, Data<Tb, 1, TG>& beta,
                       Data<TSP, 3, TG>& SurfacePressure,
                       Data<T, 4, TG>& Pressure, T P0)
  {

    int h, i, j, k;

    int Nx = Pressure.GetLength(3);
    int Ny = Pressure.GetLength(2);
    int Nz = Pressure.GetLength(1);
    int Nt = Pressure.GetLength(0);

    for (h = 0; h < Nt; h++)
      for (k = 0; k < Nz; k++)
        for (j = 0; j < Ny; j++)
          for (i = 0; i < Nx; i++)
            Pressure(h, k, j, i) = alpha(k) * P0
              + beta(k) * SurfacePressure(h, j, i);

  }


  //! Computes the altitudes from pressure fields.
  /*!
    Level heights are computed according to:
    Z_k = (-r/a)*[-1 + (P_k/P0)^(1/5.255)]
    where Z is the altitude, P the pressure,
    k the level index, r the normalized temperature,
    and a the vertical temperature gradient.
    \param Pressure pressure (Pa).
    \param Height (output) altitudes (m).
    \param a (optional) Vertical temperature gradient. Default: 0.65K for 100m.
    \param r (optional) Normalized temperature. Default: 288.15.
    \note Pressure and Height must be defined on the same grid.
  */
  template<class TPS, class TP, class TT, class T, class TG>
  void Compute4DHeight(Data<TPS, 4, TG>& SurfacePressure,
                       Data<TP, 4, TG>& Pressure,
                       Data<TT, 4, TG>& Temperature,
                       Grid<T>& Height, T g, T r)
  {

    int h, i, j, k;

    int Nx = Height.GetLength(3);
    int Ny = Height.GetLength(2);
    int Nz = Height.GetLength(1);
    int Nt = Height.GetLength(0);

    for (h = 0; h < Nt; h++)
      for (k = 0; k < Nz; k++)
        for (j = 0; j < Ny; j++)
          for (i = 0; i < Nx; i++)
            if (k == 0)
              Height.Value(h, k, j, i) = r / g * Temperature(h, k, j, i)
                * log(SurfacePressure(h, 0, j, i) / Pressure(h, k, j, i));
            else
              Height.Value(h, k, j, i) = Height.Value(h, k - 1, j, i)
                - r / g * Temperature(h, k, j, i)
                * log(Pressure(h, k, j, i) / Pressure(h, k - 1, j, i));

  }


  //! Computes the altitudes from pressure and temperature fields.
  /*!
    Level heights are computed according to:
    Z_{k+1} = Z_k + (r * T_k / g) * log(P_k/P_{k+1})
    where Z is the altitude, T the temperature, P the pressure,
    k the level index, r the molar gas constant for dry air
    and g the standard gravity.
    \par For the first level, Z_0 = r * T_0 / g * log(PS / P_0)
    where PS is the surface pressure.
    \param SurfacePressure surface pressure (Pa).
    \param Pressure pressure (Pa).
    \param Temperature temperature (K).
    \param Height (output) altitudes (m).
    \param g (optional) standard gravity. Default: 9.80665.
    \param r (optional) molar gas constant for dry air. Default: 287.0.
    \note Temperature, Pressure and Height must be defined on the same grid.
  */
  template<class TPS, class TP, class TT, class T, class TG>
  void ComputeHeight(Data<TPS, 3, TG>& SurfacePressure,
                     Data<TP, 4, TG>& Pressure,
                     Data<TT, 4, TG>& Temperature,
                     Grid<T>& Height, T g, T r)
  {

    int h, i, j, k;

    int Nx = Height.GetLength(3);
    int Ny = Height.GetLength(2);
    int Nz = Height.GetLength(1);
    int Nt = Height.GetLength(0);

    for (h = 0; h < Nt; h++)
      for (k = 0; k < Nz; k++)
        for (j = 0; j < Ny; j++)
          for (i = 0; i < Nx; i++)
            if (k == 0)
              Height.Value(h, k, j, i) = r / g * Temperature(h, k, j, i)
                * log(SurfacePressure(h, j, i) / Pressure(h, k, j, i));
            else
              Height.Value(h, k, j, i) = Height.Value(h, k - 1, j, i)
                - r / g * Temperature(h, k, j, i)
                * log(Pressure(h, k, j, i) / Pressure(h, k - 1, j, i));

  }


  //! Computes the altitudes at interfaces from pressure and
  //! temperature fields.
  /*!
    Level heights are computed according to:
    Z_{k+1/2} = Z_{k-1/2} - (r * T_k / g) * log(P_{k+1/2}/P_{k-1/2})
    where Z is the altitude, T the temperature, P the pressure,
    k the level index, r the molar gas constant for dry air
    and g the standard gravity.
    \param Pressure pressure (Pa).
    \param Temperature temperature (K).
    \param Height (output) altitudes (m).
    \param g (optional) standard gravity. Default: 9.80665.
    \param r (optional) molar gas constant for dry air. Default: 287.0.
    \param ground_set (optional) true if ground-level altitudes are set,
    false if they have to be set (they are set to zero in this case).
    Default: false.
    \note Temperature is provided at middle points (not interfaces).
    Pressure and Height are defined at interfaces (including ground-level).
  */
  template<class TP, class TT, class T, class TG>
  void ComputeInterfHeight(Data<TP, 4, TG>& Pressure,
                           Data<TT, 4, TG>& Temperature,
                           Grid<T>& Height, bool ground_set, T g, T r)
  {

    int h, i, j, k;

    int Nx = Height.GetLength(3);
    int Ny = Height.GetLength(2);
    int Nz = Height.GetLength(1) - 1;
    int Nt = Height.GetLength(0);

    if (!ground_set)
      for (h = 0; h < Nt; h++)
        for (j = 0; j < Ny; j++)
          for (i = 0; i < Nx; i++)
            Height.Value(h, 0, j, i) = T(0);

    for (h = 0; h < Nt; h++)
      for (k = 0; k < Nz; k++)
        for (j = 0; j < Ny; j++)
          for (i = 0; i < Nx; i++)
            Height.Value(h, k + 1, j, i) = Height.Value(h, k, j, i)
              - r / g * Temperature(h, k, j, i)
              * log(Pressure(h, k + 1, j, i) / Pressure(h, k, j, i));

  }


  //! Computes the altitudes at middle levels from altitudes of interfaces and
  //! pressure and temperature fields.
  /*!
    Level heights are computed according to:
    Z_k = Z_{k-1/2} + a_k * r / g * T_k
    where Z is the altitude, T the temperature, k the level index,
    r the molar gas constant for dry air, g the standard gravity and:
    a_k = 1 - P_{k+1/2} / (P_{k-1/2} - P_{k+1/2}) * log(P_{k-1/2} / P_{k+1/2})
    \param Pressure pressure (Pa).
    \param Temperature temperature (Pa).
    \param InterfHeight (output) altitudes of interfaces (m).
    \param MiddleHeight (output) altitudes of middle points (m).
    \param g (optional) standard gravity. Default: 9.80665.
    \param r (optional) molar gas constant for dry air. Default: 287.0.
    \note Temperature is provided at middle points (not interfaces).
    Pressure is defined at interfaces (including ground-level).
  */
  template<class TP, class TT, class T, class TG>
  void ComputeMiddleHeight(Data<TP, 4, TG>& Pressure,
                           Data<TT, 4, TG>& Temperature,
                           Grid<T>& InterfHeight, Grid<T>& MiddleHeight,
                           T g, T r)
  {

    int h, i, j, k;

    int Nx = MiddleHeight.GetLength(3);
    int Ny = MiddleHeight.GetLength(2);
    int Nz = MiddleHeight.GetLength(1);
    int Nt = MiddleHeight.GetLength(0);

    for (h = 0; h < Nt; h++)
      for (k = 0; k < Nz; k++)
        for (j = 0; j < Ny; j++)
          for (i = 0; i < Nx; i++)
            MiddleHeight.Value(h, k, j, i) = InterfHeight.Value(h, k, j, i)
              + r / g * Temperature(h, k, j, i)
              * (1. - Pressure(h, k + 1, j, i)
                 / (Pressure(h, k, j, i) - Pressure(h, k + 1, j, i))
                 * log(Pressure(h, k, j, i) / Pressure(h, k + 1, j, i)));

  }


  //! Computes the virtual temperature.
  /*!
    The virtual temperature is computed according to: T_v = (1 + c * q) * T
    where T_s is the virtual temperature, q the specific humidity,
    T the temperature and c a coefficient (0.608, usually).
    \param Temperature temperature (K).
    \param SpecificHumidity specific humidity (kg/kg).
    \param VirtualTemperature (output) virtual temperature (K).
    \param c (optional) coefficient. Default: 0.608.
    \note Temperature and VirtualTemperature may be the same object.
  */
  template <class TT, class TH, class T, class TG>
  void ComputeVirtualTemperature(Data<TT, 4, TG>& Temperature,
                                 Data<TH, 4, TG>& SpecificHumidity,
                                 Data<T, 4, TG>& VirtualTemperature, T c)
  {

    int h, i, j, k;

    int Nx = VirtualTemperature.GetLength(3);
    int Ny = VirtualTemperature.GetLength(2);
    int Nz = VirtualTemperature.GetLength(1);
    int Nt = VirtualTemperature.GetLength(0);

    for (h = 0; h < Nt; h++)
      for (k = 0; k < Nz; k++)
        for (j = 0; j < Ny; j++)
          for (i = 0; i < Nx; i++)
            VirtualTemperature(h, k, j, i) = Temperature(h, k, j, i)
              * (1 + c * SpecificHumidity(h, k, j, i)
                 / (1. - SpecificHumidity(h, k, j, i)));

  }


}  // namespace AtmoData.

#define ATMODATA_FILE_METEOROLOGY_CXX
#endif
