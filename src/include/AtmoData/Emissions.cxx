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


#ifndef ATMODATA_FILE_EMISSIONS_CXX

#include "Emissions.hxx"

namespace AtmoData
{


  //////////////
  // BIOGENIC //
  //////////////


  //! Computes biogenic emissions rates.
  /*!
    Computes biogenic emission rates according to Simpson et al. (1999).
    \param LUC Land use coverage ([0, 1]).
    \param Density densities of the land categories.
    \param EF_isoprene isoprene emission factors (\mu g . m^{-2} . s^{-1}).
    \param EF_terpenes terpenes emission factors (\mu g . m^{-2} . s^{-1}).
    \param EF_NO NO emission factors (\mu g . m^{-2} . s^{-1}).
    \param Isoprene (output) isoprene emission rates (\mu g . m^{-2} . s^{-1}).
    \param Terpenes (output) terpenes emission rates (\mu g . m^{-2} . s^{-1}).
    \param NO (output) NO emission rates (\mu g . m^{-2} . s^{-1}).
  */
  template < class TL, class TD, class TEFI, class TEFT,
             class TEFN, class TI, class TT, class TN, class TG >
  void ComputeBiogenicRates(Data<TL, 3, TG>& LUC, Data<TD, 1, TG>& Density,
                            Data<TEFI, 1, TG>& EF_isoprene,
                            Data<TEFT, 1, TG>& EF_terpenes,
                            Data<TEFN, 1, TG>& EF_NO,
                            Data<TI, 2, TG>& Isoprene,
                            Data<TT, 2, TG>& Terpenes,
                            Data<TN, 2, TG>& NO)
  {

    int i, j, k;

    int Nx = Isoprene.GetLength(1);
    int Ny = Isoprene.GetLength(0);
    int Nc = LUC.GetLength(0);

    Isoprene.SetZero();
    Terpenes.SetZero();
    NO.SetZero();

    for (k = 0; k < Nc; k++)
      for (j = 0; j < Ny; j++)
        for (i = 0; i < Nx; i++)
          {
            Isoprene(j, i) += Density(k) * EF_isoprene(k) * LUC(k, j, i);
            Terpenes(j, i) += Density(k) * EF_terpenes(k) * LUC(k, j, i);
            NO(j, i) += Density(k) * EF_NO(k) * LUC(k, j, i);
          }

  }


  //! Computes biogenic emissions.
  /*!
    Computes biogenic emissions according to Simpson et al. (1999).
    \param Temperature soil or leaf temperature.
    \param PAR Photosynthetically active radiation.
    \param LUC Land use coverage.
    \param Density densities of the land categories.
    \param EF_isoprene isoprene emission factors (\mu g . m^{-2} . s^{-1}).
    \param EF_terpenes terpenes emission factors (\mu g . m^{-2} . s^{-1}).
    \param EF_NO NO emission factors (\mu g . m^{-2} . s^{-1}).
    \param Isoprene (output) isoprene emissions (\mu g . m^{-2} . s^{-1}).
    \param Terpenes (output) terpenes emissions (\mu g . m^{-2} . s^{-1}).
    \param NO (output) NO emissions (\mu g . m^{-2} . s^{-1}).
  */
  template < class TTemp, class TP, class TL, class TD, class TEFI, class TEFT,
             class TEFN, class TI, class TT, class TN, class TG >
  void ComputeBiogenicEmissions(Data<TTemp, 3, TG>& Temperature,
                                Data<TP, 3, TG>& PAR,
                                Data<TL, 3, TG>& LUC,
                                Data<TD, 1, TG>& Density,
                                Data<TEFI, 1, TG>& EF_isoprene,
                                Data<TEFT, 1, TG>& EF_terpenes,
                                Data<TEFN, 1, TG>& EF_NO,
                                Data<TI, 3, TG>& Isoprene,
                                Data<TT, 3, TG>& Terpenes,
                                Data<TN, 3, TG>& NO)
  {

    int h, i, j, k;

    int Nx = Isoprene.GetLength(2);
    int Ny = Isoprene.GetLength(1);
    int Nt = Isoprene.GetLength(0);
    int Nc = LUC.GetLength(0);

    TTemp T, Ts_NO;
    TP par;
    double c_l, c_t;

    // For environmental correction factors.
    const double alpha(0.0027), c_l1(1.066);
    const double Ts(303), Tm(314), R(8.314), c_t1(95000), c_t2(230000);
    const double alpha2(alpha * alpha), ratio1(c_t1 / (R * Ts)),
      ratio2(c_t2 / (R * Ts));

    Isoprene.SetZero();
    Terpenes.SetZero();
    NO.SetZero();

    for (h = 0; h < Nt; h++)
      for (i = 0; i < Nx; i++)
        for (j = 0; j < Ny; j++)
          for (k = 0; k < Nc; k++)
            {
              // Temperature.
              T = Temperature(h, j, i);
              // Photosynthetically active radiation.
              par = PAR(h, j, i);
              // Light dependence.
              c_l = alpha * c_l1 * par
                / sqrt(1 + alpha2 * par * par);
              // Temperature dependence.
              c_t = exp(ratio1 * (T - Ts) / T)
                / (1. + exp(ratio2 * (T - Tm) / T));
              // Emission.
              Isoprene(h, j, i) += Density(k) * EF_isoprene(k)
                * LUC(k, j, i) * c_l * c_t;

              // Emission.
              Terpenes(h, j, i) += Density(k) * EF_terpenes(k)
                * LUC(k, j, i) * exp(0.09 * (T - Ts));

              // Emission.
              if (EF_NO(k) == 0.9)
                Ts_NO = .67 * (T - 273.15) + 8.8;
              else
                Ts_NO = .84 * (T - 273.15) + 3.6;
              if (Ts_NO > 35)
                Ts_NO = 35;
              if (Ts_NO < 15)
                Ts_NO = 15;
              NO(h, j, i) += EF_NO(k) * LUC(k, j, i) * exp(0.071 * (Ts_NO));
            }

  }


  ///////////////////
  // ANTHROPOGENIC //
  ///////////////////


  //! Constructor.
  /*!
    \param emission emission rate.
    \param country EMEP country code number.
  */
  template <class T>
  EmepCountryEmission<T>::EmepCountryEmission(T emission, int country):
    emission_(emission), country_(country)
  {
  }


  //! Constructor.
  /*!
    \param N number of countries.
  */
  TimeZone::TimeZone(int N): countries_(N), local_times_(N)
  {
    for (vector<int>::iterator iter = countries_.begin();
         iter != countries_.end(); ++iter)
      {
        *iter = iter - countries_.begin();
        local_times_.at(iter - countries_.begin()) = 1;
      }
  }


  //! Initializes from an input file.
  /*!
    \param file_name input file name.
  */
  void TimeZone::Init(string file_name)
  {
    int country, time;
    ExtStream file_stream(file_name);
    while (!file_stream.IsEmpty())
      {
        file_stream.GetElement();
        file_stream.GetElement(country);
        file_stream.GetElement(time);
        if (country < int(local_times_.size()))
          local_times_.at(country) = time;
        file_stream.GetLine();
      }
  }


  //! Returns the time zone offset for a given country number.
  /*!
    \param i country number.
    \return The time zone offset from GMT.
  */
  int TimeZone::operator()(int i) const
  {
    return local_times_.at(i);
  }


  //! Reads temporal factors from a file.
  /*!
    Reads monthly profiles for a given month or daily profiles for a week. In
    the file, each line first stores the EMEP country code, then the SNAP
    sector and finally the factors for the 12 months or the 7 weekdays.  If no
    factor is available for a given country and a given SNAP sector, it is set
    to 1.
    \param input_file input file.
    \param Factors (output) temporal factors.
  */
  template <class real>
  void GetTemporalFactors(string input_file, Data<real, 3>& Factors)
  {

    string line;
    int s, Ncountry;

    Factors.Fill(1.0);

    ConfigStream input_stream(input_file);
    if (!input_stream.is_open())
      throw string("Error in GetTemporalFactors: \"") + input_file
        + "\" cannot be opened.";

    while (input_stream.GetLine(line))
      {
        vector<real> factors;
        split(line, factors);

        if (int(factors.size()) != 2 + Factors.GetLength(2))
          throw string("Error in GetTemporalFactors: \"") + input_file
            + "\" badly formatted.";

        Ncountry = int(factors[0]);
        s = int(factors[1]);

        if (s - 1 < Factors.GetLength(1))
          for (int i = 0; i < Factors.GetLength(2); i++)
            Factors(Ncountry, s - 1, i) = factors[i + 2];
      }

  }


  //! Computes speciation/aggregation factors.
  /*!
    Computes the factor for each sector to be applied to inventory species
    emissions in order to get model species emissions. For NOX, the factor
    takes into account the fact that NOX emissions are given in NO2 equivalent
    units. First, One file per inventory species is read for speciation in the
    speciation directory (format "<species>.dat") and provides percentages to
    split inventory species emissions into real-species emissions. Then, an
    aggregation file is read. It contains a matrix that gives molecular
    weights, reactivities and aggregation coefficients for real (lines) and
    model (columns) species. VOCs aggregation coefficients are computed
    according to Middleton et al. (1990) with the integral of [OH] set to 1.0.
    \param Sp_emis_names inventory species names.
    \param aggregation_matrix aggregation-matrix file.
    \param speciation_directory speciation directory.
    \param Sp_model_names (output) model species names.
    \param Species_factor (output) species factors indexed by the model
    species, the inventory species and the sector.
  */
  template <class real>
  void SpeciationAggregation(const vector<string>& Sp_emis_names,
                             string aggregation_matrix,
                             string speciation_directory,
                             vector<string>& Sp_model_names,
                             Data<real, 3>& Species_factor)
  {

    int Nsectors = Species_factor[2].GetLength();
    int Nsp_emis = Species_factor[1].GetLength();

    string Sp_real_names;

    string input_file, line, stmp;
    string species_names;

    vector<real> vtmp;
    vector<string> Sp_model_names_tmp;
    vector<real> molecular_weights_model;
    vector<real> React_model;

    /** Reads model species names, molecular weight and reactivity ***/

    ExtStream aggregation_stream(aggregation_matrix);
    if (!aggregation_stream.is_open())
      throw string(" File ") + input_file + " doesn't exist.";

    /*** Reads speciation ***/

    RegularGrid<real> GridSectors(Species_factor[2]);
    Data<real, 1> Speciation_coeff(GridSectors);

    real Aggregation_coeff;
    real React_real;
    real molecular_weights_real;

    const real M_NO2 = 46.0;
    const real M_SO2 = 64.0;
    const real int_OH(1.0);

    Species_factor.Fill(0.0);

    int mm = 0;
    for (int l = 0; l < Nsp_emis; l++)
      {
        input_file = speciation_directory + Sp_emis_names[l] + ".dat";
        ExtStream speciation_stream(input_file);
        if (!speciation_stream.is_open())
          throw string("File \"") + input_file + "\" doesn't exist.";

        while (!speciation_stream.IsEmpty())
          {
            Sp_real_names = speciation_stream.GetElement();
            speciation_stream >> molecular_weights_real;
            for (int k = 0; k < Nsectors; k++)
              speciation_stream >> Speciation_coeff(k);

            // NOX is given in NO2 equivalent units.
            if (Sp_emis_names[l] == "NOX")
              for (int k = 0; k < Nsectors; k++)
                Speciation_coeff(k) *= molecular_weights_real / M_NO2;

            // SOX is given in SO2 equivalent units.
            if (Sp_emis_names[l] == "SOX")
              for (int k = 0; k < Nsectors; k++)
                Speciation_coeff(k) *= molecular_weights_real / M_SO2;

            if (Sp_emis_names[l] != "NMVOC")
              {
                Sp_model_names[mm] = Sp_real_names;
                for (int k = 0; k < Nsectors; k++)
                  Species_factor(mm, l, k) = Speciation_coeff(k) * 0.01;
                mm++;
              }
            else
              {
                ExtStream aggregation_stream(aggregation_matrix);

                // Species names.
                aggregation_stream.SkipElements(3);
                aggregation_stream.GetLine(line);
                split(line, Sp_model_names_tmp);
                for (int m = 0; m < int(Sp_model_names_tmp.size()); m++)
                  Sp_model_names[mm + m] = Sp_model_names_tmp[m];
                // Molecular weights.
                aggregation_stream.SkipElements(3);
                aggregation_stream.GetLine(line);
                split(line, molecular_weights_model);
                // Reactivities.
                aggregation_stream.SkipElements(3);
                aggregation_stream.GetLine(line);
                split(line, React_model);

                stmp = "";
                while (!aggregation_stream.IsEmpty() && Sp_real_names != stmp)
                  {
                    aggregation_stream >> stmp;
                    line = aggregation_stream.GetLine();
                  }
                if (aggregation_stream.IsEmpty() && Sp_real_names != stmp)
                  throw string("Species ") + Sp_real_names + " not found.";

                split(line, vtmp);
                molecular_weights_real = vtmp[0];
                React_real = vtmp[1];

                for (int m = 0; m < int(Sp_model_names_tmp.size()); m++)
                  {
                    Aggregation_coeff = vtmp[m + 2];
                    for (int k = 0; k < Nsectors; k++)
                      {
                        Species_factor(m + mm, l, k) += Speciation_coeff(k) * 0.01
                          * Aggregation_coeff
                          * molecular_weights_model[m] / molecular_weights_real
                          * (1. - exp(-int_OH * React_real))
                          / (1. - exp(-int_OH * React_model[m]));
                      }
                  }
                if (aggregation_stream.bad())
                  throw "Aggregation file is badly formatted";
              }
          }
      }

  }


  //! Computes vertical distribution factors.
  /*!
    Maps the vertical distribution 'vertical_distribution_in' given on the
    vertical mesh 'GridZ_interf_in' to the output mesh 'GridZ_interf_out'.
    \param vertical_distribution_in vertical distribution on input grid.
    \param GridZ_interf_in altitudes of interfaces of input grid (m).
    \param GridZ_interf_out altitudes of interfaces of output grid (m).
    \param vertical_distribution_out (output) vertical distribution on output
    grid.
    \note The vertical distributions are indexed by the activity sector and
    the height.
  */
  template <class real>
  void ComputeVerticalDistribution(const Data<real, 2>& vertical_distribution_in,
                                   const RegularGrid<real>& GridZ_interf_in,
                                   const RegularGrid<real>& GridZ_interf_out,
                                   Data<real, 2>& vertical_distribution_out)
  {

    int k, s;

    int Nsectors = vertical_distribution_in[0].GetLength();

    // Vertical Grids.
    int Nz_in = vertical_distribution_in[1].GetLength();
    int Nz_out = vertical_distribution_out[1].GetLength();

    RegularGrid<real> DeltaZ_in(Nz_in);
    RegularGrid<real> DeltaZ_out(Nz_out);

    string input_file, line;

    /*** Reads vertical distribution heights ***/

    // Sets values at nodes.
    for (k = 0; k < Nz_in; k++)
      DeltaZ_in(k) = GridZ_interf_in(k + 1) - GridZ_interf_in(k);
    for (k = 0; k < Nz_out; k++)
      DeltaZ_out(k) = GridZ_interf_out(k + 1) - GridZ_interf_out(k);

    vertical_distribution_out.Fill(0);

    /*** Linear interpolation of the vertical distribution ***/

    int k_in = 0;
    for (k = 0; k < Nz_out && k_in < Nz_in; k++)
      {
        if (GridZ_interf_out(k + 1) <= GridZ_interf_in(k_in + 1))
          for (s = 0; s < Nsectors; s++)
            vertical_distribution_out(s, k) =
              vertical_distribution_in(s, k_in) / DeltaZ_in(k_in)
              * DeltaZ_out(k);
        else
          {
            for (s = 0; s < Nsectors; s++)
              vertical_distribution_out(s, k) =
                vertical_distribution_in(s, k_in) / DeltaZ_in(k_in)
                * (GridZ_interf_in(k_in + 1) - GridZ_interf_out(k));
            k_in++;
            while (k_in < Nz_in
                   && GridZ_interf_out(k + 1) > GridZ_interf_in(k_in + 1))
              {
                for (s = 0; s < Nsectors; s++)
                  vertical_distribution_out(s, k) +=
                    vertical_distribution_in(s, k_in);
                k_in++;
              }
            if (k_in < Nz_in)
              for (s = 0; s < Nsectors; s++)
                vertical_distribution_out(s, k) +=
                  vertical_distribution_in(s, k_in) / DeltaZ_in(k_in)
                  * (GridZ_interf_out(k + 1) - GridZ_interf_in(k_in));
          }
      }
  }


  //! Divides the vertical distribution by the heights.
  /*!
    \param GridZ_interf_out heights of the output grid.
    \param vertical_distribution_out (output) vertical distribution on the
    output grid indexed by the activity sector and the height.
  */
  template <class real>
  void DividesByHeights(const RegularGrid<real>& GridZ_interf_out,
                        Data<real, 2>& vertical_distribution_out)
  {

    int k, s;

    int Nsectors = vertical_distribution_out[0].GetLength();

    // Vertical Grids.
    int Nz_out = vertical_distribution_out[1].GetLength();

    RegularGrid<real> DeltaZ_out(Nz_out);

    // Sets values at nodes.
    for (k = 0; k < Nz_out; k++)
      DeltaZ_out(k) = GridZ_interf_out(k + 1) - GridZ_interf_out(k);

    for (s = 0; s < Nsectors; s++)
      for (k = 0; k < Nz_out; k++)
        vertical_distribution_out(s, k) /= DeltaZ_out(k);
  }


  //! Correspondence between Emep and Polair3D grids.
  /*!
    Computes (1) the number of urban, water, forest and other types of LUC that
    are included in each EMEP cell, and (2) the total number of LUC cells
    included in each Polair3D cell.
    \param LUC land use coverage on GLCF categories.
    \param Nurb_emep (output) number of urban LUC cells in EMEP cells.
    \param Nwat_emep (output) number of water LUC cells in EMEP cells.
    \param Nfor_emep (output) number of forest LUC cells in EMEP cells.
    \param Noth_emep (output) number of other LUC cells in EMEP cells.
    \param Ntot_polair (output) total number of LUC cells in each Polair3D
    cell.
    \note The arrays are stored in (Y, X) format.
  */
  template <class real>
  void GridCorrespondences(const Data<int, 2, real>& LUC,
                           Data<int, 2, real>& Nurb_emep,
                           Data<int, 2, real>& Nwat_emep,
                           Data<int, 2, real>& Nfor_emep,
                           Data<int, 2, real>& Noth_emep,
                           Data<int, 2, real>& Ntot_polair)
  {

    const real pi = 3.14159265358979323846264;

    const real xpol(8.);
    const real ypol(110.);
    const real M = 6370. / 50. * (1. + sin(pi / 3.));

    /*** Coordinates ***/

    int Nx_luc = LUC[1].GetLength();
    real x_min_center_luc = LUC[1](0);
    real delta_x_luc = (LUC[1](Nx_luc - 1) - LUC[1](0)) / real(Nx_luc - 1);
    real x_min_luc = x_min_center_luc - delta_x_luc / 2.0;

    int Ny_luc = LUC[0].GetLength();
    real y_min_center_luc = LUC[0](0);
    real delta_y_luc = (LUC[0](Ny_luc - 1) - LUC[0](0)) / real(Ny_luc - 1);
    real y_min_luc = y_min_center_luc - delta_y_luc / 2.0;

    int Nx_emep = Nurb_emep[1].GetLength();
    int Ny_emep = Nurb_emep[0].GetLength();

    int Nx = Ntot_polair[1].GetLength();
    real x_min_center = Ntot_polair[1](0);
    real Delta_x = (Ntot_polair[1](Nx - 1) - Ntot_polair[1](0)) / real(Nx - 1);
    real x_min = x_min_center - Delta_x / 2.0;

    int Ny = Ntot_polair[0].GetLength();
    real y_min_center = Ntot_polair[0](0);
    real Delta_y = (Ntot_polair[0](Ny - 1) - Ntot_polair[0](0)) / real(Ny - 1);
    real y_min = y_min_center - Delta_y / 2.0;

    /*** Computes output data ***/

    Ntot_polair.Fill(0);

    Nurb_emep.Fill(0);
    Nwat_emep.Fill(0);
    Nfor_emep.Fill(0);
    Noth_emep.Fill(0);

    for (int i = 0; i < Nx_luc; i++)
      for (int j = 0; j < Ny_luc; j++)
        {
          real lon = (x_min_center_luc + i * delta_x_luc) * pi / 180.; // rad
          real lat = (y_min_center_luc + j * delta_y_luc) * pi / 180.; // rad

          int i_emep = int(xpol + M * tan(pi / 4. - lat / 2.)
                           * sin(lon + 32. * pi / 180.) - 0.5);
          int j_emep = int(ypol - M * tan(pi / 4. - lat / 2.)
                           * cos(lon + 32. * pi / 180.) - 0.5);

          if (i_emep >= 0 && i_emep < Nx_emep && j_emep >= 0
              && j_emep < Ny_emep)
            if (LUC(j, i) == 13)
              Nurb_emep(j_emep, i_emep) += 1;
            else if (LUC(j, i) == 0)
              Nwat_emep(j_emep, i_emep) += 1;
            else if (LUC(j, i) >= 1 && LUC(j, i) <= 6)
              Nfor_emep(j_emep, i_emep) += 1;
            else if (LUC(j, i) >= 7 && LUC(j, i) <= 12)
              Noth_emep(j_emep, i_emep) += 1;
            else
              throw "Error in GridCorrespondences: LUC index out of range.";

          lon = x_min_luc + i * delta_x_luc; // deg
          lat = y_min_luc + j * delta_y_luc; // deg

          int i_polair = -999;
          int j_polair = -999;
          if (lon >= x_min)
            i_polair = int((lon - x_min) / Delta_x);
          if (lat >= y_min)
            j_polair = int((lat - y_min) / Delta_y);

          if (i_polair >= 0 && i_polair < Nx && j_polair >= 0
              && j_polair < Ny)
            Ntot_polair(j_polair, i_polair) += 1;
        }
  }


  //! Reads Emep standard files and applies monthly and daily factors.
  /*!
    First reads Emep standard files located in the input directory (format
    <species>.dat). In these files, the Emep emissions (Tons) correspond to a
    given year, and are provided for European countries and for each
    SNAP. Then applies monthly and daily coefficients in order to provide
    emissions over land and over water for a given day.
    \param date date.
    \param Sp_emis_names inventory emission species names.
    \param input_directory input directory.
    \param input_file input file.
    \param MonthlyFactors monthly factors indexed by the country, the sector
    and the month.
    \param DailyFactors daily factors indexed by the country, the sector and
    the day.
    \param deposition_factor_nh3 local deposition factor for NH3 (range 0-1).
    \param Emis_land (output) emissions over land (Tons) for a given day.
    \param Emis_water (output) emissions over water (Tons) for a given day.
    \note Output data is indexed by the species, Y, X and the sector.
  */
  template <class real>
  void ReadEmep(Date date, const vector<string>& Sp_emis_names,
                string input_directory, string input_file,
                const Data<real, 3>& MonthlyFactors,
                const Data<real, 3>& DailyFactors,
                const real& deposition_factor_nh3,
                Data<list<EmepCountryEmission<real> >, 4, real>& Emis_land,
                Data<list<EmepCountryEmission<real> >, 4, real>& Emis_water)
  {

    int Ncountry_max = MonthlyFactors.GetLength(0);
    int Ncountry, i, j, s;
    string country, line, sector;
    vector<string> v;
    real quantity, monthly, daily;

    int Nsp_emis = Emis_land[0].GetLength();

    RegularGrid<real> GridSectors(Emis_land[3]);
    Data<real, 1> Emis_emep(GridSectors);

    // Reads the country codes and country numbers.
    ExtStream file_stream(input_file);
    vector<string> CountryCode;
    vector<int> CountryNumber;
    while (!file_stream.IsEmpty())
      {
        file_stream.GetElement(country);
        CountryCode.push_back(country);
        file_stream.GetElement(Ncountry);
        CountryNumber.push_back(Ncountry);
        file_stream.GetLine();
      }

    // Finds the day of the week.
    int day = date.GetWeekDay();

    for (int l = 0; l < Nsp_emis; l++)
      {
        string emis_file = input_directory + Sp_emis_names[l] + ".dat";
        ExtStream EmepEmisStream(emis_file);
        if (!EmepEmisStream.is_open())
          throw string("File ") + emis_file + " doesn't exist.";

        while (has_element(EmepEmisStream))
          {
            EmepEmisStream.GetLine(line);
            v = split(line, ";");
            country = v[0];
            sector = v[2];
            s = to_num<int>(sector.erase(0, 1)) - 1;
            i = to_num<int>(v[4]);
            j = to_num<int>(v[5]);

            quantity = to_num<real>(v[7]);

            // NH3 local deposition.
            // 9 is agricultural sector.
            if (Sp_emis_names[l] == "NH3" && s == 9)
              quantity *= (1. - deposition_factor_nh3);

            if (i < 1 || j < 1)
              continue;

            int n = 0;
            while (CountryCode[n] != country)
              {
                n++;
                if (n >= int(CountryCode.size()))
                  throw string("Error in ReadEmep: country code \"") + country
                    + string("\" not found in \"") + input_file + "\".";
              }
            Ncountry = CountryNumber[n];

            if (Ncountry >= Ncountry_max)
              throw string("Country code number ") + to_str(Ncountry)
                + string(" (country \"")
                + country + string("\") is greater than ")
                + string("or equal to the maximum number of countries (")
                + to_str(Ncountry_max) + ").";

            if (s < 10)
              {
                monthly = MonthlyFactors(Ncountry, s, date.GetMonth() - 1);
                daily = DailyFactors(Ncountry, s, day);

                // In water.
                if ((Ncountry >= 30 && Ncountry <= 35) || Ncountry == 70)
                  Emis_water(l, j - 1, i - 1, s).
                    push_back(EmepCountryEmission<real>(quantity *
                                                        monthly * daily /
                                                        365, Ncountry));

                // In land.
                else
                  Emis_land(l, j - 1, i - 1, s).
                    push_back(EmepCountryEmission<real>(quantity *
                                                        monthly * daily /
                                                        365, Ncountry));
              }
          }

        if (EmepEmisStream.bad())
          throw string("EMEP emission file \"") + emis_file
            + "\" is badly formatted.";
      }
  }


  //! Reads Emep modified files and applies monthly and daily factors.
  /*!
    First reads Emep modified files located in the input directory (format
    <species>.dat). In these files, the Emep emissions (Tons) correspond to a
    given year, and are provided for European countries and for each
    SNAP. Then applies monthly and daily coefficients in order to provide
    emissions over land and over water for a given day.
    \param date date.
    \param Sp_emis_names inventory emission species names.
    \param input_directory input directory.
    \param MonthlyFactors monthly factors indexed by the country, the sector
    and the month.
    \param DailyFactors daily factors indexed by the country, the sector and
    the day.
    \param deposition_factor_nh3 local deposition factor for NH3 (range 0-1).
    \param Emis_land (output) emissions over land (Tons) for a given day.
    \param Emis_water (output) emissions over water (Tons) for a given day.
    \note Output data is indexed by the species, Y, X and the sector.
    \warning The input files read by this function are not the standard Emep
    emission files. This function is therefore deprecated. Please use
    ReadEmep.
  */
  template <class real>
  void ReadModifiedEmep(Date date, const vector<string>& Sp_emis_names,
                        string input_directory,
                        const Data<real, 3>& MonthlyFactors,
                        const Data<real, 3>& DailyFactors,
                        const real& deposition_factor_nh3,
                        Data<list<EmepCountryEmission<real> >, 4, real>& Emis_land,
                        Data<list<EmepCountryEmission<real> >, 4, real>& Emis_water)
  {

    int Ncountry, i, j;
    string line;

    int Nsp_emis = Emis_land[0].GetLength();
    int Nsectors = Emis_land[3].GetLength();
    RegularGrid<real> GridSectors(Emis_land[3]);
    Data<real, 1> Emis_emep(GridSectors);

    int AgriculturalSector = 10;

    // Finds the day of the week.
    int day = date.GetWeekDay();

    for (int l = 0; l < Nsp_emis; l++)
      {
        string input_file = input_directory + Sp_emis_names[l] + ".dat";
        ExtStream EmepEmisStream(input_file);
        if (!EmepEmisStream.is_open())
          throw string("File ") + input_file + " doesn't exist.";

        while (has_element(EmepEmisStream))
          {
            EmepEmisStream >> Ncountry >> i >> j;
            for (int s = 0; s < Nsectors; s++)
              EmepEmisStream >> Emis_emep(s);
            getline(EmepEmisStream, line);

            // NH3 local deposition
            if (Sp_emis_names[l] == "NH3")
              Emis_emep(AgriculturalSector - 1) *= (1. - deposition_factor_nh3);

            // In water.
            if (Ncountry >= 30 && Ncountry <= 35)
              for (int s = 0; s < Nsectors; s++)
                Emis_water(l, j - 1, i - 1, s).push_back(EmepCountryEmission<real>(Emis_emep(s)
                                                                                   * MonthlyFactors(Ncountry, s, date.GetMonth() - 1)
                                                                                   * DailyFactors(Ncountry, s, day) / 365.,
                                                                                   Ncountry));
            // In land.
            else
              for (int s = 0; s < Nsectors; s++)
                Emis_land(l, j - 1, i - 1, s).push_back(EmepCountryEmission<real>(Emis_emep(s)
                                                                                  * MonthlyFactors(Ncountry, s, date.GetMonth() - 1)
                                                                                  * DailyFactors(Ncountry, s, day) / 365.,
                                                                                  Ncountry));
          }
        if (EmepEmisStream.bad())
          throw string("EMEP emission file \"") + input_file
            + "\" is badly formatted.";
      }
  }


  //! Computes emissions on a latitude-longitude grid.
  /*!
    Computes emissions on a latitude-longitude grid (like Polair3D grid) on
    the basis of EMEP inventory and weights associated with aggregated (in
    urban, forest and miscellaneous) LUC categories.
    \param LUC land use coverage.
    \param Ratio_urb weight of urban cells.
    \param Ratio_for weight of forest cells.
    \param Ratio_oth weight of other cells.
    \param Nurb_emep number of urban LUC cells in EMEP cells.
    \param Nwat_emep number of water LUC cells in EMEP cells.
    \param Nfor_emep number of forest LUC cells in EMEP cells.
    \param Noth_emep number of other LUC cells in EMEP cells.
    \param Ntot_polair total number of LUC cells in each Polair3D cell.
    \param Emis_land emissions over land (Tons).
    \param Emis_water emissions over water (Tons).
    \param Emis_out (output) output emissions (Tons . m^{-2}) for a given day.
  */
  template <class real>
  void EmepToLatLon(const Data<int, 2, real>& LUC,
                    const real Ratio_urb,
                    const real Ratio_for,
                    const real Ratio_oth,
                    const Data<int, 2, real>& Nurb_emep,
                    const Data<int, 2, real>& Nwat_emep,
                    const Data<int, 2, real>& Nfor_emep,
                    const Data<int, 2, real>& Noth_emep,
                    const Data<int, 2, real>& Ntot_polair,
                    Data<list<EmepCountryEmission<real> >, 4, real>& Emis_land,
                    Data<list<EmepCountryEmission<real> >, 4, real>& Emis_water,
                    Data<list<EmepCountryEmission<real> >, 4, real>& Emis_out)
  {
    int l, i, j, s;
    typename list<EmepCountryEmission<real> >::const_iterator iter;
    typename list<EmepCountryEmission<real> >::iterator iter_out;

    const real pi = 3.14159265358979323846264;
    const real Earth_radius = 6370997.;

    const real xpol(8.);
    const real ypol(110.);
    const real M = 6370. / 50. * (1. + sin(pi / 3.));

    int Nsectors = Emis_land.GetLength(3);
    int Nsp_emis = Emis_land.GetLength(0);

    int Nx_luc = LUC[1].GetLength();
    real x_min_center_luc = LUC[1](0);
    real delta_x_luc = (LUC[1](Nx_luc - 1) - LUC[1](0)) / real(Nx_luc - 1);
    real x_min_luc = x_min_center_luc - delta_x_luc / 2.0;

    int Ny_luc = LUC[0].GetLength();
    real y_min_center_luc = LUC[0](0);
    real delta_y_luc = (LUC[0](Ny_luc - 1) - LUC[0](0)) / real(Ny_luc - 1);
    real y_min_luc = y_min_center_luc - delta_y_luc / 2.0;

    int Nx_emep = Nurb_emep[1].GetLength();
    int Ny_emep = Nurb_emep[0].GetLength();

    int Nx = Ntot_polair[1].GetLength();
    real x_min_center = Ntot_polair[1](0);
    real Delta_x = (Ntot_polair[1](Nx - 1) - Ntot_polair[1](0)) / real(Nx - 1);
    real x_min = x_min_center - Delta_x / 2.0;

    int Ny = Ntot_polair[0].GetLength();
    real y_min_center = Ntot_polair[0](0);
    real Delta_y = (Ntot_polair[0](Ny - 1) - Ntot_polair[0](0)) / real(Ny - 1);
    real y_min = y_min_center - Delta_y / 2.0;

    real Ratio(0.0);

    for (i = 0; i < Nx_luc; i++)
      for (j = 0; j < Ny_luc; j++)
        {
          real lon = (x_min_center_luc + i * delta_x_luc) * pi / 180.; // rad
          real lat = (y_min_center_luc + j * delta_y_luc) * pi / 180.; // rad

          int i_emep = int(xpol + M * tan(pi / 4. - lat / 2.)
                           * sin(lon + 32. * pi / 180.) - 0.5);
          int j_emep = int(ypol - M * tan(pi / 4. - lat / 2.)
                           * cos(lon + 32. * pi / 180.) - 0.5);

          lon = x_min_luc + i * delta_x_luc; // deg
          lat = y_min_luc + j * delta_y_luc; // deg

          int i_polair = int ((lon - x_min) / Delta_x);
          int j_polair = int ((lat - y_min) / Delta_y);

          if (i_emep >= 0 && i_emep < Nx_emep && j_emep >= 0 && j_emep < Ny_emep
              && i_polair >= 0 && i_polair < Nx && j_polair >= 0 && j_polair < Ny)
            {
              if (LUC(j, i) == 0)
                {
                  Ratio = 1. / real(Nwat_emep(j_emep, i_emep));
                  for (l = 0; l < Nsp_emis; l++)
                    for (s = 0; s < Nsectors; s++)
                      for (iter = Emis_water(l, j_emep, i_emep, s).begin();
                           iter != Emis_water(l, j_emep, i_emep, s).end(); ++iter)
                        {
                          iter_out = Emis_out(l, s, j_polair, i_polair).begin();
                          while (iter_out != Emis_out(l, s, j_polair, i_polair).end()
                                 && iter_out->country_ != iter->country_)
                            ++iter_out;
                          if (iter_out != Emis_out(l, s, j_polair, i_polair).end())
                            iter_out->emission_ += iter->emission_ * Ratio;
                          else
                            Emis_out(l, s, j_polair, i_polair).push_back(EmepCountryEmission<real>(iter->emission_ * Ratio, iter->country_));
                        }
                }
              else
                {
                  if (LUC(j, i) == 13)
                    Ratio = Ratio_urb
                      / real(Ratio_urb * Nurb_emep(j_emep, i_emep)
                             + Ratio_for * Nfor_emep(j_emep, i_emep)
                             + Ratio_oth * Noth_emep(j_emep, i_emep));

                  if (LUC(j, i) >= 1 && LUC(j, i) <= 6)
                    Ratio = Ratio_for
                      / real(Ratio_urb * Nurb_emep(j_emep, i_emep)
                             + Ratio_for * Nfor_emep(j_emep, i_emep)
                             + Ratio_oth * Noth_emep(j_emep, i_emep));

                  if (LUC(j, i) >= 7 && LUC(j, i) <= 12)
                    Ratio = Ratio_oth
                      / real(Ratio_urb * Nurb_emep(j_emep, i_emep)
                             + Ratio_for * Nfor_emep(j_emep, i_emep)
                             + Ratio_oth * Noth_emep(j_emep, i_emep));

                  for (l = 0; l < Nsp_emis; l++)
                    for (s = 0; s < Nsectors; s++)
                      for (iter = Emis_land(l, j_emep, i_emep, s).begin();
                           iter != Emis_land(l, j_emep, i_emep, s).end(); ++iter)
                        {
                          iter_out = Emis_out(l, s, j_polair, i_polair).begin();
                          while (iter_out != Emis_out(l, s, j_polair, i_polair).end()
                                 && iter_out->country_ != iter->country_)
                            ++iter_out;
                          if (iter_out != Emis_out(l, s, j_polair, i_polair).end())
                            iter_out->emission_ += iter->emission_ * Ratio;
                          else
                            Emis_out(l, s, j_polair, i_polair).push_back(EmepCountryEmission<real>(iter->emission_ * Ratio, iter->country_));
                        }
                }
            }
        }

    const real ratio_pi = pi / 180.;
    const real factor = ratio_pi * ratio_pi * Earth_radius * Earth_radius
      * delta_x_luc * delta_y_luc;
    real surface;

    // Divides by the surface. Emissions are then given in Tons/m^2.
    for (i = 0; i < Nx; i++)
      for (j = 0; j < Ny; j++)
        {
          surface = Ntot_polair(j, i) * factor
            * cos(Ntot_polair[0](j) * ratio_pi);
          for (l = 0; l < Nsp_emis; l++)
            for (s = 0 ; s < Nsectors; s++)
              for (iter_out = Emis_out(l, s, j, i).begin();
                   iter_out != Emis_out(l, s, j, i).end(); ++iter_out)
                iter_out->emission_ /= surface;
        }

  }


}  // namespace AtmoData.

#define ATMODATA_FILE_EMISSIONS_CXX
#endif
