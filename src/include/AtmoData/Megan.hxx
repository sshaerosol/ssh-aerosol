// Copyright (C) 2013, ENPC
//    Author(s): Florian Couvidat
//
// This file is part of the air quality modeling system Polyphemus.
//
// Polyphemus is developed in the INRIA - ENPC joint project-team CLIME and in
// the ENPC - EDF R&D joint laboratory CEREA.
//
// Polyphemus is free software; you can redistribute it and/or modify it under
// the terms of the GNU General Public License as published by the Free
// Software Foundation; either version 2 of the License, or (at your option)
// any later version.
//
// Polyphemus is distributed in the hope that it will be useful, but WITHOUT
// ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
// FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for
// more details.
//
// For more information, visit the Polyphemus web site:
//      http://cerea.enpc.fr/polyphemus/


#ifndef ATMODATA_FILE_MEGAN_HXX


#include <string>

#include "Data.hxx"


namespace AtmoData
{

  using namespace std;
  using namespace SeldonData;

  inline void assert_defined_species(const string& species);

  template <class T>
  T light_dependent_factor(const string& species);

  template <class T>
  T GammaTemp(const string& species, T temperature);

  template <class T>
  T GammaLAI(T LAI);

  template <class T>
  T SolarAngle(int iday, T localhour, T lat);

  template <class T>
  T GammaLight(int iday, T sinangle, T PAR, T daily_PAR);

  template <class T>
  T GammaAGE(const string& species, T LAIp, T LAIc,
             T daily_temperature, T t);

  template <class T>
  T GammaSM(T WP, T SoilMoisture);

  template <class T>
  void MeganAggregation(const Data<T, 4>& emis_out,
                        const string& file_aggregation,
                        const vector<string>& biogenic_names,
                        const vector<string>& output_names,
                        Data<T, 4>& emissions);

}  // namespace AtmoData.

#define ATMODATA_FILE_MEGAN_HXX
#endif
