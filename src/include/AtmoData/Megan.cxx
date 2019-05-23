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


#ifndef ATMODATA_FILE_MEGAN_CXX


#include "Megan.hxx"

#include "Files.hxx"
using namespace Talos;


namespace AtmoData
{

  inline void assert_defined_species(const string& species)
  {
    if (species != "ISOP" && species != "MBO" && species != "FORM"
        && species != "CO" && species != "MYRC" && species != "SABI"
        && species != "LIMO" && species != "3CAR" && species != "OCIM"
        && species != "BPIN" && species != "APIN" && species != "OMTP"
        && species != "FARN" && species != "BCAR" && species != "OSQT"
        && species != "MEOH" && species != "ACTO" && species != "NO"
        && species != "CH4" && species != "ACTA")
      throw "[ERROR] Species " + species + " undefined.";
  }


  template <class T>
  T light_dependent_factor(const string& species)
  {
    T ldf;
    if (species == "ISOP" || species == "MBO")
      ldf = 0.9999;
    else if (species == "MYRC" || species == "LIMO" || species == "3CAR")
      ldf = 0.05;
    else if (species == "SABI" || species == "BPIN" || species == "APIN"
             || species == "OMTP")
      ldf = 0.1;
    else if (species == "OCIM")
      ldf = 0.8;
    else if (species == "MEOH" || species == "CH4")
      ldf = 0.75;
    else if (species == "ACTO")
      ldf = 0.25;
    else if (species == "NO")
      ldf = 0.0;
    else if (species == "ACTA" || species == "FARN" || species == "BCAR"
             || species == "OSQT" || species == "FORM" || species == "CO")
      ldf = 0.5;
    else
      throw "[ERROR] Species " + species
        + " light dependent factor undefined.";
    return ldf;
  }


  template <class T>
  T GammaTemp(const string& species, T temperature)
  {
    T temperature_factor;
    const T Ts = 303.0;

    if (species == "ISOP" || species == "MBO" || species == "FORM"
        || species == "CO")
      temperature_factor = 0.09;
    else if (species == "MYRC" || species == "SABI" || species == "LIMO"
             || species == "3CAR" || species == "OCIM" || species == "BPIN"
             || species == "APIN" || species == "OMTP")
      temperature_factor = 0.1;
    else if (species == "FARN" || species == "BCAR" || species == "OSQT")
      temperature_factor = 0.17;
    else if (species == "MEOH")
      temperature_factor = 0.08;
    else if (species == "ACTO" || species == "NO")
      temperature_factor = 0.11;
    else if (species == "CH4")
      temperature_factor = 0.05;
    else if (species == "ACTA")
      temperature_factor = 0.13;
    else
      throw "[ERROR] Species " + species + " temperature factor undefined.";

    return exp(temperature_factor * (temperature - Ts));
  }


  template <class T>
  T GammaTempISOP(T temperature, T daily_temperature)
  {
    const T CT1 = 80.;
    const T CT2 = 200.;
    T Eopt = 1.75 * exp(0.08 * (daily_temperature - 297));
    T Topt = 313. + 0.6 * (daily_temperature - 297.);
    T x = ((1 / Topt) - (1 / temperature)) / 0.00831;
    return Eopt * CT2 * exp(CT1 * x)
      / (CT2 - CT1 * (1. - exp(CT2 * x)));
  }


  template <class T>
  T GammaLAI(T LAI)
  {
    return 0.49 * LAI / (sqrt((1. + 0.2 * LAI * LAI)));
  }


  template <class T>
  void SolarAngle(int iday, T localhour, T lat)
  {
    const T pi = 3.14159265358979323846264;
    T sindelta = -sin(0.40907) * cos(6.28 * float((iday + 10)) / 365.0);
    T cosdelta = pow((1 - pow(sindelta, 2)), 0.5);
    T A = sin(lat * pi / 180.0) * sindelta;
    T B = cos(lat * pi / 180.0) * cosdelta;
    return A + B * cos(2.0 * pi * (localhour - 12.0) / 24.0);
  }


  template <class T>
  T GammaLight(int iday, T sinangle, T PAR, T daily_PAR)
  {
    T gamma = 0.0;
    if (sinangle >= 0.0)
      {
        T top_PAR = 3000.0
          + 99.0 * cos(2.0 * 3.14 - float((iday - 10)) / 365.0);
        T phi = PAR / (sinangle * top_PAR);
        T B = 1.0 + 0.0005 * (daily_PAR - 400.0);
        T A = 2.46 * B * phi - 0.9 * pow(phi, 2);
        gamma = sinangle * A;
        if (gamma < 0.0)
          gamma = 0.0;
      }
    return gamma;
  }


  template <class T>
  T GammaAGE(const string& species, T LAIp, T LAIc,
             T daily_temperature, T t)
  {
    T Fnew, Fgro, Fmat, Fold;
    T ti, tm;
    T Anew, Agro, Amat, Aold;

    if (species == "ISOP" || species == "MBO")
      {
        Anew = 0.05;
        Agro = 0.6;
        Amat = 1.125;
        Aold = 1.0;
      }
    else if (species == "ACTO" || species == "ACTA" || species == "FORM"
             || species == "CH4" || species == "NO" || species == "CO")
      {
        Anew = 1.0;
        Agro = 1.0;
        Amat = 1.0;
        Aold = 1.0;
      }
    else if (species == "MYRC" || species == "SABI" || species == "LIMO"
             || species == "3CAR" || species == "OCIM" || species == "BPIN"
             || species == "APIN" || species == "OMTP")
      {
        Anew = 2.0;
        Agro = 1.8;
        Amat = 0.95;
        Aold = 1.0;
      }
    else if (species == "FARN" || species == "BCAR" || species == "OSQT")
      {
        Anew = 0.4;
        Agro = 0.6;
        Amat = 1.075;
        Aold = 1.0;
      }
    else if (species == "MEOH")
      {
        Anew = 3.0;
        Agro = 2.6;
        Amat = 0.85;
        Aold = 1.0;
      }
    else
      throw "[ERROR] Species " + species + " canopy factors undefined.";

    if (daily_temperature > 303.0)
      ti = 2.9;
    else
      ti = 5.0 + (0.7 * (300.0 - daily_temperature));
    tm = 2.3 * ti;

    if (LAIp < LAIc)
      {
        if (t <= ti)
          Fnew = 1 - LAIp / LAIc;
        else
          Fnew = (ti / t) * (1 - LAIp / LAIc);
        if (t <= tm)
          Fmat = LAIp / LAIc;
        else
          Fmat = LAIp / LAIc + (t - tm) / t * (1 - LAIp / LAIc);
        Fgro = 1 - Fnew - Fmat;
        Fold = 0.0;
      }
    else if (LAIp == LAIc)
      {
        Fmat = 0.8;
        Fnew = 0.0;
        Fgro = 0.1;
        Fold = 0.1;
      }
    else
      {
        Fnew = 0.0;
        Fgro = 0.0;
        Fold = (LAIp - LAIc) / LAIp;
        Fmat = 1 - Fold;
      }
    return Fold * Aold + Fnew * Anew + Fgro * Agro + Fmat * Amat;
  }


  template <class T>
  T GammaSM(T WP, T SoilMoisture)
  {
    T deltaSM = 0.06;
    T thetal = WP + deltaSM;
    T gamma;
    if (SoilMoisture > thetal)
      gamma = 1.0;
    else if (SoilMoisture < WP)
      gamma = 0.0;
    else
      gamma = (SoilMoisture - WP) / deltaSM;
    return gamma;
  }


  template <class T>
  void MeganAggregation(const Data<T, 4>& emis_out,
                        const string& file_aggregation,
                        const vector<string>& biogenic_names,
                        const vector<string>& output_names,
                        Data<T, 4>& emissions)
  {
    int Nsp_in = emis_out.GetLength(0);
    int Ny = emis_out.GetLength(1);
    int Nx = emis_out.GetLength(2);
    int Nt = emis_out.GetLength(3);
    int Nsp_out = emissions.GetLength(0);
    int Nmegan = 20;

    string name_in, name_out;
    RegularGrid<int> GridInSpecies(Nsp_in);
    Data<int, 1> indices(GridInSpecies);

    int i, j, k, x, y, t;
    T coefficient;

    emissions.Fill(0.0);

    ExtStream aggregation_stream(file_aggregation);
    if (!aggregation_stream.is_open())
      throw "[ERROR] File \"" + file_aggregation + "\" does not exist.";

    while (!aggregation_stream.IsEmpty())
      {
        aggregation_stream >> name_out;
        if (name_out == "SPECIES")
          for (i = 0; i < Nmegan; i++)
            {
              aggregation_stream >> name_in;
              for (j = 0; j < Nsp_in; j++)
                if (name_in == biogenic_names[j])
                  indices(j) = i;
            }
        else
          {
            for (i = 0; i < Nsp_out; i++)
              if (name_out == output_names[i])
                for (j = 0; j < Nmegan; j++)
                  {
                    aggregation_stream >> coefficient;
                    if (coefficient != 0.0)
                      for (k = 0; k < Nsp_in; k++)
                        if (indices(k) == j)
                          for (x = 0; x < Nx; x++)
                            for (y = 0; y < Ny; y++)
                              for (t = 0; t < Nt; t++)
                                emissions(i, t, y, x) += coefficient
                                  * emis_out(k, y, x, t);
                  }
          }
      }
  }


}  // namespace AtmoData.

#define ATMODATA_FILE_MEGAN_CXX
#endif
