// Copyright (C) 2007, ENPC - INRIA - EDF R&D
// Author(s): Meryem Ahmed de Biasi
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

#ifndef COMMON_FILE_COMMON_HXX

#include <cmath>
#include <cstring>
#include <iostream>
#include <vector>
using namespace std;

#include "TalosHeader.hxx"
using namespace Talos;

#include "AtmoDataHeader.hxx"
using namespace AtmoData;

namespace Polyphemus
{
  void parse_argument(int argc, char** argv, string& config_file,
                      string default_name = "");
  void parse_argument(int argc, char** argv, string& first_file,
                      string& second_file, string default_name = "");
  void parse_argument(int argc, char** argv, string& first_file,
                      string& second_file, Date& date,
                      string default_name = "");
  void parse_argument(int argc, char** argv, string& first_file,
                      string& second_file, Date& date_beg, Date& date_end,
                      string default_name = "");

  Date read_date_MM5(string FileName);

#ifdef COMMON_WITH_NETCDF

  Date read_date_WRF(string FileName);
  float read_delta_t_WRF(string FileName);

#endif

  Date convert_delta(string delta, Date date_beg);

  template <class T>
  int compute_Nt(Date date_beg, Date date_end, T Delta_t);

  template<class T>
  void abs_(T& x);

  template<class T>
  void read_refractiveindex_tabulation(string file_species,
                                       string file_water_refractive_index,
                                       string directory_OPAC,
                                       int& black_carbon_index,
                                       int Nwater_wavelength,
                                       int N_OPAC_wavelength,
                                       RegularGrid<T> GridWavelength,
                                       Data<T, 2>& PureSpeciesIndexReal,
                                       Data<T, 2>& PureSpeciesIndexImaginary,
                                       Data<T, 1>& WaterIndexReal,
                                       Data<T, 1>& WaterIndexImaginary,
                                       Data<string, 1, T>& SpeciesNames);

  template<class T>
  void read_Mie_tabulation(string directory_efficiency_factor,
                           RegularGrid<T>& GridIndexReal,
                           RegularGrid<T>& GridIndexImaginary,
                           RegularGrid<T>& GridDiameter,
                           RegularGrid<T> GridWavelength,
                           Data<T, 4>& AbsorptionEfficiencyFactor,
                           Data<T, 4>& ExtinctionEfficiencyFactor,
                           Data<T, 5>& PhaseFunctionTable);

  template<class T>
  void compute_optical_properties(Data<T, 4>& OpticalThickness,
                                  Data<T, 4>& SingleScatteringAlbedo,
                                  Data<T, 4>& MeanExtinctionEfficiencyFactor,
                                  Data<T, 4>& MeanAbsorbtionEfficiencyFactor,
                                  Data<T, 5>& PhaseFunction,
                                  RegularGrid<T> GridZ_interf_out,
                                  RegularGrid<T> GridIndexReal,
                                  RegularGrid<T> GridIndexImaginary,
                                  RegularGrid<T> GridDiameter,
                                  Data<T, 4> AbsorptionEfficiencyFactor,
                                  Data<T, 4> ExtinctionEfficiencyFactor,
                                  Data<T, 5> PhaseFunctionTable,
                                  Data<T, 5> ConcentrationAerosol,
                                  Data<T, 4> ConcentrationWater,
                                  Data<T, 2> PureSpeciesIndexReal,
                                  Data<T, 2> PureSpeciesIndexImaginary,
                                  Data<T, 1> WaterIndexReal,
                                  Data<T, 1> WaterIndexImaginary,
                                  Data<T, 1> dry_diameter_computed,
                                  Data<T, 4> Wet_diameter,
                                  T dry_aerosol_density,
                                  int black_carbon_index,
                                  int option_black_carbon_treatment,
                                  int option_wet_index,
                                  int option_well_mixed_index);


} //namespace Polyphemus

#define COMMON_FILE_COMMON_HXX
#endif
