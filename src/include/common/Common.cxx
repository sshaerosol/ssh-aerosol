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


#ifndef COMMON_FILE_COMMON_CXX

#include "Common.hxx"

#include "Talos.hxx"

#include "AtmoData.hxx"

namespace Polyphemus
{

  //! Parse the arguments of a preprocessing program.
  /*!
    \param argc number of system arguments.
    \param argv table of system arguments.
    \param config_file (output) configuration file for the program.
    \param default_name (optional) default name of the configuration file.
    Default: "".
  */
  void parse_argument(int argc, char** argv, string& config_file,
                      string default_name)
  {
    string msg;
    if (default_name != "")
      {
        msg =  "Usage:\n";
        msg += string("  ") + argv[0] + " [configuration file]\n";
        msg += string("  ") + argv[0] + "\n\n";
        msg += "Arguments:\n";
        msg += "  [configuration file] (optional): configuration";
        msg += " file. Default: " + default_name + "\n";
      }
    else
      {
        msg =  "Usage:\n";
        msg += string("  ") + argv[0] + " [configuration file]\n\n";
        msg += "Arguments:\n";
        msg += "  [main configuration file] : main configuration file.\n";
      }

    if (argc > 2 || (default_name == "" && argc < 2))
      throw msg;

    if (argc == 1)
      config_file = default_name;
    else
      config_file = argv[1];

  }


  //! Parse the arguments of a preprocessing program.
  /*!
    \param argc number of system arguments.
    \param argv table of system arguments.
    \param first_file (output) main configuration file for the program.
    \param second_file (output) secondary configuration file for the program.
    \param default_name (optional) default name of the main configuration
    file. Default: "".
  */
  void parse_argument(int argc, char** argv, string& first_file,
                      string& second_file, string default_name)
  {
    string msg;
    if (default_name != "")
      {
        msg =  "Usage:\n";
        msg += string("  ") + argv[0] + " [main configuration file]";
        msg += "[secondary config file] \n";
        msg += string("  ") + argv[0] + " [main configuration file] \n";
        msg += string("  ") + argv[0] + "\n\n";
        msg += "Arguments:\n";
        msg += "  [main configuration file] (optional): main configuration";
        msg += " file. Default: " + default_name + "\n";
        msg += "  [secondary configuration file] (optional): secondary";
        msg += " configuration file.\n";
      }
    else
      {
        msg =  "Usage:\n";
        msg += string("  ") + argv[0] + " [main configuration file]";
        msg += "[secondary config file] \n";
        msg += string("  ") + argv[0] + " [main configuration file] \n\n";
        msg += "Arguments:\n";
        msg += "  [main configuration file] : main configuration file.\n";
        msg += "  [secondary configuration file] (optional): secondary";
        msg += " configuration file. \n";
      }

    if (argc > 3 || (default_name == "" && argc < 2))
      throw msg;

    if (argc == 1)
      first_file = default_name;
    else
      first_file = argv[1];

    if (argc != 3)
      second_file = "";
    else
      second_file = argv[2];
  }


  //! Parse the arguments of a preprocessing program.
  /*!
    \param argc number of system arguments.
    \param argv table of system arguments.
    \param first_file (output) main configuration file for the program.
    \param second_file (output) secondary configuration file for the program.
    \param date (output) beginning date for the program.
    \param default_name (optional) default name of the main configuration
    file. Default: "".
  */
  void parse_argument(int argc, char** argv, string& first_file,
                      string& second_file, Date& date,
                      string default_name)
  {
    string msg;
    if (default_name != "")
      {
        msg =  "Usage:\n";
        msg += string("  ") + argv[0] + " [main configuration file] ";
        msg += "[secondary config file] [date]\n";
        msg += string("  ") + argv[0] + " [main configuration file] [date]\n";
        msg += string("  ") + argv[0] + " [date] \n\n";
        msg += "Arguments:\n";
        msg += "  [main configuration file] (optional): main configuration";
        msg += " file. Default: " + default_name + "\n";
        msg += "  [secondary configuration file] (optional): secondary";
        msg += " configuration file.\n";
        msg += "  [date]: date in any valid format. \n";
      }
    else
      {
        msg =  "Usage:\n";
        msg += string("  ") + argv[0] + " [main configuration file] ";
        msg += "[secondary config file] [date]\n";
        msg += string("  ") + argv[0] + " [main configuration file]"
          + " [date]\n\n";
        msg += "Arguments:\n";
        msg += "  [main configuration file] : main configuration file.\n";
        msg += "  [secondary configuration file] (optional): secondary";
        msg += " configuration file. \n";
        msg += "  [date]: date in any valid format. \n";
      }

    if (argc > 4 || !is_date(argv[argc - 1]) ||
        (default_name == "" && argc < 3) ||
        (default_name != "" && argc < 2))
      throw msg;

    if (argc == 2)
      first_file = default_name;
    else
      first_file = string(argv[1]);

    if (argc != 4)
      second_file = "";
    else
      second_file = string(argv[2]);

    if (!is_date(argv[argc - 1]))
      throw "Date " + string(argv[argc - 1]) + "is not in a valid format.";
    else
      date = string(argv[argc - 1]);
  }

  //! Parse the arguments of a preprocessing program.
  /*!
    \param argc number of system arguments.
    \param argv table of system arguments.
    \param first_file (output) main configuration file for the program.
    \param second_file (output) secondary configuration file for the program.
    \param date_beg (output) beginning date for the program.
    \param date_end (output) end date for the program.
    \param default_name (optional) default name of the main configuration
    file. Default: "".
  */
  void parse_argument(int argc, char** argv, string& first_file,
                      string& second_file, Date& date_beg, Date& date_end,
                      string default_name)
  {
    string msg;
    if (default_name != "")
      {
        msg =  "Usage:\n";
        msg += string("  ") + argv[0] + " [main configuration file]";
        msg += " [secondary config file] [first date]";
        msg += " [second date/interval]\n";
        msg += string("  ") + argv[0] + " [main configuration file]";
        msg += " [first date] [second date/interval] \n";
        msg += string("  ") + argv[0] + " [main configuration file]";
        msg += " [secondary config file] [first date] \n";
        msg += string("  ") + argv[0] + " [first date]";
        msg += " [second date/interval] \n";
        msg += string("  ") + argv[0] + " [first date]  \n\n";
        msg += "Arguments:\n";
        msg += "  [main configuration file] (optional): main configuration";
        msg += " file. Default: " + default_name + "\n";
        msg += "  [secondary configuration file] (optional): secondary";
        msg += " configuration file. Default: \"\".\n";
        msg += "  [first date]: beginning date in any valid format. \n";
        msg += "  [second date]: end date in any valid format. \n";
        msg += "  [interval] (optional): Interval in format NdMh or Nd-Mh";
        msg += " or Nd or Mh where N is the number of\n";
        msg += "  days and M the number of hours. Default: 1d. \n\n";
        msg += "Note:\n";
        msg += "  The end date, whether it is given directly or computed by";
        msg += " adding the time interval to the beginning\n";
        msg += "  date, is always considered as excluded. \n";

      }
    else
      {
        msg =  "Usage:\n";
        msg += string("  ") + argv[0] + " [main configuration file]";
        msg += " [secondary config file] [first date]";
        msg += " [second date/interval]\n";
        msg += string("  ") + argv[0] + " [main configuration file]";
        msg += " [first date] [second date/interval] \n";
        msg += string("  ") + argv[0] + " [main configuration file]";
        msg += " [secondary config file] [first date] \n";
        msg += string("  ") + argv[0] + " [main configuration file]";
        msg += " [first date] \n\n";
        msg += "Arguments:\n";
        msg += "  [main configuration file] : main configuration file. \n";
        msg += "  [secondary configuration file] (optional): secondary";
        msg += " configuration file. Default: \"\".\n";
        msg += "  [first date]: beginning date in any valid format. \n";
        msg += "  [second date]: end date in any valid format. \n";
        msg += "  [interval] (optional): Interval in format NdMh or Nd-Mh";
        msg += " or Nd or Mh where N is the number of days and M the number";
        msg += " of hours. Default: 1d. \n\n";
        msg += "Note:\n";
        msg += "The end date, whether it is given directly or computed by";
        msg += " adding the time interval to the beginning date, is always";
        msg += " considered as excluded. \n";
      }

    if (argc > 5 || (default_name == "" && argc < 3) ||
        (default_name != "" && argc < 2) || (!is_date(argv[argc - 1]) &&
                                             !is_delta(argv[argc - 1])))
      throw msg;

    if (argc == 2)
      {
        first_file = default_name;
        second_file = "";
        date_beg = string(argv[1]);
        date_end = date_beg;
        date_end.AddDays(1);
      }
    else if (argc == 3)
      {
        if (is_date(argv[1]))
          {
            first_file = default_name;
            second_file = "";
            date_beg = string(argv[1]);
            if (is_date(argv[2]))
              date_end = string(argv[2]);
            else
              date_end = convert_delta(string(argv[2]), date_beg);
          }
        else
          {
            first_file = string(argv[1]);
            second_file = "";
            date_beg = string(argv[2]);
            date_end = date_beg;
            date_end.AddDays(1);
          }
      }
    else if (argc == 4)
      {
        first_file = string(argv[1]);
        if (is_date(argv[2]))
          {
            second_file = "";
            date_beg = string(argv[2]);
            if (is_date(argv[3]))
              date_end =  string(argv[3]);
            else
              date_end = convert_delta(string(argv[3]), date_beg);
          }
        else
          {
            second_file = string(argv[2]);
            date_beg = string(argv[3]);
            date_end = date_beg;
            date_end.AddDays(1);
          }
      }
    else // argc == 5
      {
        first_file = string(argv[1]);
        second_file = string(argv[2]);
        date_beg = string(argv[3]);
        if (is_date(argv[4]))
          date_end = string(argv[4]);
        else
          date_end = convert_delta(string(argv[4]), date_beg);
      }
  }


  //! Reads the beginning date of a MM5 file.
  /*! Note that the beginning date in the BigHeader is not necessarily the
    date we need.
    Indeed, data may have been generated for a long period of time but stored
    in several files. In that case, the date given in the BigHeader is the
    beginning date for the first file.
    \param FileName name of the MM5 file.
    \return The beginning  date.
  */
  Date read_date_MM5(string FileName)
  {
    ifstream file(FileName.c_str());
    FormatMM5 InputMeteo;
    MM5SubHeader SubHeader;

    while (InputMeteo.ReadFlag(file) != 1)
      {
        Array<int, 2> BHI;
        Array<float, 2> BHR;
        Array<string, 2> BHIC;
        Array<string, 2> BHRC;
        InputMeteo.ReadBigHeader(file, BHI, BHR, BHIC, BHRC);
      }
    InputMeteo.ReadSubHeader(file, SubHeader);
    string current_date = SubHeader.GetCurrentDate();
    current_date = current_date.substr(0, 19);
    current_date = find_replace(current_date, ":", "-");
    Date date(current_date);
    return date;
  }


#ifdef COMMON_WITH_NETCDF

  //! Reads the beginning date of a WRF file.
  /*!
    \param FileName name of the WRF file.
    \return The beginning  date.
  */
  Date read_date_WRF(string FileName)
  {
    string current_date = "";
    int i, Nt_in;
    FormatNetCDF<float> WRF;
    WRF.ReadDimension(FileName, "Times", 0, Nt_in);
    RegularGrid<char> GridTimeStep_in(Nt_in);
    RegularGrid<char> GridDateFormat_in(19);
    Data<char, 2> Times(GridTimeStep_in, GridDateFormat_in);
    WRF.Read(FileName, "Times", Times);
    for (i = 0; i < 19; i++)
      current_date += Times(0, i);
    current_date = find_replace(current_date, ":", "-");
    Date date(current_date);
    return date;
  }

  //! Reads the output timestep (hour) in a WRF file.
  /*!
    \param FileName name of the WRF file.
    \return The timestep delta_t
  */
  float read_delta_t_WRF(string FileName)
  {
    string first_date(""), second_date("");
    int i, Nt_in;
    float delta_t_in;
    FormatNetCDF<float> WRF;
    WRF.ReadDimension(FileName, "Times", 0, Nt_in);
    RegularGrid<char> GridTimeStep_in(Nt_in);
    RegularGrid<char> GridDateFormat_in(19);
    Data<char, 2> Times(GridTimeStep_in, GridDateFormat_in);
    WRF.Read(FileName, "Times", Times);
    for (i = 0; i < 19; i++)
      {
        first_date += Times(0, i);
        second_date += Times(1, i);
      }
    first_date = find_replace(first_date, ":", "-");
    second_date = find_replace(second_date, ":", "-");
    Date date1(first_date), date2(second_date);
    delta_t_in = float(date2.GetSecondsFrom(date1)) / 3600;
    return delta_t_in;
  }


#endif


  //! Computes the end date from the beginning date and the time interval.
  /*!
    \param delta time interval, in the format "NdMh" or "Nd-Mh" where N is the
    number of days (integer) and M the number of hours (integer). If N is
    zero, the time interval can be expressed as "Mh".
    \param date_beg beginning date.
    \return The end date.
  */
  Date convert_delta(string delta, Date date)
  {
    Date date_end = date;

    if (!is_delta(delta))
      throw "Wrong time period " + delta + string(".");
    else
      {
        vector<string> period = split(delta, "dh-_");
        int days = 0, hours = 0;
        if (period.size() == 2)
          {
            to_num(period[0], days);
            to_num(period[1], hours);
          }
        else if (delta.find('h', 0) != string::npos)
          to_num(period[0], hours);
        else if (delta.find('d', 0) != string::npos)
          to_num(period[0], days);
        date_end.AddDays(days);
        date_end.AddHours(hours);
        return date_end;
      }
  }


  //! Compute the number of iterations from two dates and a time step.
  /*! The number of iterations is computed to go from the beginning date
    included to the end date excluded.
    If the time between the beginning and end date of the preprocessing is not
    a multiple of the time step, an exception is raised.
    \param date_beg beginning date of the preprocessing.
    \param date_end end date of the preprocessing.
    \param Delta_t time step in hours.
    \return The number of iteration for the preprocessing (integer).
  */
  template<class T>
  int compute_Nt(Date date_beg, Date date_end, T Delta_t)
  {
    int distance = int(date_end.GetSecondsFrom(date_beg));
    int Delta_t_sec =  int(Delta_t * 3600);
    if (distance % Delta_t_sec != 0)
      throw string("The time between the beginning and end date must be a ")
        + "multiple of the time step (" + to_str(Delta_t) + " hours).";
    int Nt = distance / Delta_t_sec;
    return Nt;
  }

  template<class T>
  void abs_(T& x)
  {
    if (x < T(0))
      x = -x;
  }

  //! read the OPAC tabulation for pure species refractive index .
  /*!
    \param file_species file describing model - OPAC specie correspondance.
    \param file_water_refractive_index file where water refractive
    \                    index are stored for several wavelenghts.
    \param directory_OPAC OPAC directory.
    \param black_carbon_index index of the black carbon specie in
    \                               the Data ConcentrationAerosol.
    \param Nwater_wavelength number of wavelenght bins in the OPAC
    \                       tabulation for water refractive index.
    \param N_OPAC_wavelength number of wavelenght bins in the OPAC
    \                                                  tabulation.
    \param GridWavelength output grid wavelenght.
    \param PureSpeciesIndexReal  real part of the pure specie
    \             refractive indexes at the output wavelenght.
    \param PureSpeciesIndexImaginary imaginary part of the pure
    \        specie refractive indexes at the output wavelenght.
    \param WaterIndexReal real part of the water refractive index
    \                                    at the output wavelenght.
    \param WaterIndexImaginary imaginary part of the water refractive
    \                                 index at the output wavelenght.
    \param SpeciesNames model pure specie names.
  */
  template<class T>
  void read_refractiveindex_tabulation(int Nspecies,
                                       string file_species,
                                       string file_water_refractive_index,
                                       string directory_OPAC,
                                       int& black_carbon_index,
                                       int Nwater_wavelength,
                                       int N_OPAC_wavelength,
                                       RegularGrid<T> GridWavelength,
                                       Data<T, 2>& PureSpeciesIndexReal,
                                       Data<T, 2>& PureSpeciesIndexImaginary,
                                       Data<T, 1>& WaterIndexReal,
                                       Data<T, 1>& WaterIndexImaginary)
  {
    string description;
    FormatText input_tab(",");
    int i, j, k, Nwavelength;
    string file_efficiency_factor, file_grid_Mie;

    RegularGrid<T> GridSpecies(Nspecies);

    Data<string, 1, T> SpeciesNames(GridSpecies);
    Data<string, 1, T> OPACNames(GridSpecies);
    RegularGrid<T> GridWavelengthOPAC(N_OPAC_wavelength);
    Data<T, 2> IndexReal(GridSpecies, GridWavelengthOPAC);
    Data<T, 2> IndexImaginary(GridSpecies, GridWavelengthOPAC);

    Nwavelength = WaterIndexReal.GetLength(0);

    // READ PURE SPECIES REFRACTIVE INDEX (OPAC) + Define black_carbon_index

    FormatFormattedText match_OPAC(string("<e><e>"));
    match_OPAC.SetDelimiters("\t");
    match_OPAC.Read(file_species, "0", SpeciesNames);
    match_OPAC.Read(file_species, "1", OPACNames);

    for (i = 0; i < Nspecies; i++)
      {
        if (SpeciesNames(i) == "PBC")
          black_carbon_index = i;
        ExtStream  file_OPAC(directory_OPAC + string("/") + OPACNames(i));
        cout << "OPACNames" << directory_OPAC + string("/")
          + OPACNames(i) << endl;
        file_OPAC.SkipLines(17);
        for (j = 0; j < N_OPAC_wavelength; j++)
          {
            file_OPAC.SkipElements(1);
            file_OPAC.GetElement(GridWavelengthOPAC(j));
            file_OPAC.SkipElements(6);
            file_OPAC.GetElement(IndexReal(i, j));
            file_OPAC.GetElement(IndexImaginary(i, j));
          }
        file_OPAC.Close();
      }

    // Interpolates refractive index of pure species on desired wavelengths.
    for (k = 0; k < Nwavelength; k++)
      {
        int index_in_OPAC_table(-999);
        T coeff1(-999.), coeff2(-999.);
        for (i = 0; i < N_OPAC_wavelength - 1; i++)
          {
            if (GridWavelengthOPAC(i) <= GridWavelength(k)
                && GridWavelengthOPAC(i + 1) > GridWavelength(k))
              {
                index_in_OPAC_table = i;
                coeff2 = (GridWavelength(k) - GridWavelengthOPAC(i))
                  / (GridWavelengthOPAC(i + 1) - GridWavelengthOPAC(i));
                coeff1 = (GridWavelengthOPAC(i + 1) - GridWavelength(k))
                  / (GridWavelengthOPAC(i + 1) - GridWavelengthOPAC(i));

                for (j = 0; j < Nspecies; j++)
                  {
                    PureSpeciesIndexReal(j, k) =
                      coeff1 * IndexReal(j, index_in_OPAC_table)
                      + coeff2 * IndexReal(j, index_in_OPAC_table + 1);

                    PureSpeciesIndexImaginary(j, k) =
                      coeff1 * IndexImaginary(j, index_in_OPAC_table)
                      + coeff2 * IndexImaginary(j, index_in_OPAC_table + 1);
                  }
              }
          }

        if (index_in_OPAC_table == -999)
          throw string("Unable to find required wavelength in OPAC.");
      }
    // For Mie calculations, imaginary indices are positive.
    PureSpeciesIndexImaginary.Apply(abs_);

    // READ REFRACTIVE INDEX FOR PURE WATER
    RegularGrid<T> GridWaterWavelength(Nwater_wavelength);
    Data<T, 1> WaterIndexInputReal(GridWaterWavelength);
    Data<T, 1> WaterIndexInputImaginary(GridWaterWavelength);
    ifstream water_refractive_index_data
      (file_water_refractive_index.c_str());
    getline(water_refractive_index_data, description);
    input_tab.Read(water_refractive_index_data, GridWaterWavelength);
    input_tab.Read(water_refractive_index_data, WaterIndexInputReal);
    input_tab.Read(water_refractive_index_data, WaterIndexInputImaginary);
    water_refractive_index_data.close();

    LinearInterpolationRegular(WaterIndexInputReal, WaterIndexReal);
    LinearInterpolationRegular(WaterIndexInputImaginary,
                               WaterIndexImaginary);
  }

  //! Read tabulation for optical properties (extinction and
  //  absorbtion efficiency factors, single scattering albedo and
  //  optical thickness and the 8 first term of the phase function),
  //  giving the knowledge of the aerosol refractive index
  //  and its diameter. The tabulation is obtained from a Mie code
  //  (Mishchenko et al., 1999).
  /*!
    \param directory_efficiency_factor directory of the tabulation.
    \param GridIndexReal tabulation input: real part of the aerosol
    \                                              refractive index.
    \param GridIndexImaginary tabulation input: imaginary part
    \                         of the aerosol refractive index.
    \param GridDiameter tabulation input: aerosol diameter.
    \param GridWavelength output grid wavelenght.
    \param AbsorptionEfficiencyFactor absorbtion efficiency factor.
    \param ExtinctionEfficiencyFactor extinction efficiency factor.
    \param PhaseFunctionTable tabulation output: 8 first terms of
    \                                         the phase function.
  */
  template<class T>
  void read_Mie_tabulation(string directory_efficiency_factor,
                           RegularGrid<T>& GridIndexReal,
                           RegularGrid<T>& GridIndexImaginary,
                           RegularGrid<T>& GridDiameter,
                           RegularGrid<T> GridWavelength,
                           Data<T, 4>& AbsorptionEfficiencyFactor,
                           Data<T, 4>& ExtinctionEfficiencyFactor,
                           Data<T, 5>& PhaseFunctionTable)
  {

    string description;
    FormatText input_tab(",");
    int i, j, l, k, ld, Nwavelength, Ndiameter;
    int tabulation_index_real, tabulation_index_imaginary;
    string file_efficiency_factor, file_grid_Mie;

    RegularGrid<T> GridLegendre(8);

    Ndiameter = ExtinctionEfficiencyFactor.GetLength(3);
    Nwavelength = ExtinctionEfficiencyFactor.GetLength(0);
    tabulation_index_real = ExtinctionEfficiencyFactor.GetLength(1);
    tabulation_index_imaginary = ExtinctionEfficiencyFactor.GetLength(2);

    // READ EFFICICENCY AND ABSORPTION FACTOR IN TABULATION

    file_grid_Mie = directory_efficiency_factor + "/Mish_Grid_Mie.dat";

    ifstream stream_efficiency_grid(file_grid_Mie.c_str());
    cout << "Read Grid file" << endl;
    stream_efficiency_grid >> description;
    input_tab.Read(stream_efficiency_grid, GridIndexReal);
    input_tab.Read(stream_efficiency_grid, GridIndexImaginary);
    input_tab.Read(stream_efficiency_grid, GridDiameter);
    stream_efficiency_grid.close();


    for (k = 0; k < Nwavelength; k++)
      {
        cout << "Computing optical thickness at " <<
          GridWavelength(k) * 1.e3  << " nm..." << endl;
        cout << "    Reading efficiency factors file...";
        cout.flush();

        Data<T, 3> AbsorptionEfficiencyFactor_tmp(GridIndexReal,
                                                  GridIndexImaginary,
                                                  GridDiameter);
        Data<T, 3> ExtinctionEfficiencyFactor_tmp(GridIndexReal,
                                                  GridIndexImaginary,
                                                  GridDiameter);
        Data<T, 4> PhaseFunctionTable_tmp(GridIndexReal,
                                          GridIndexImaginary,
                                          GridDiameter, GridLegendre);

        file_efficiency_factor = directory_efficiency_factor
          + string("/Mish_efficiency_factors_tab_")
          + to_str(GridWavelength(k) * 1.e3) + ".dat";

        ifstream stream_efficiency_data(file_efficiency_factor.c_str());
        cout << "Read Data file" << endl;
        stream_efficiency_data >> description;
        stream_efficiency_data >> description;
        input_tab.Read(stream_efficiency_data,
                       AbsorptionEfficiencyFactor_tmp);
        input_tab.Read(stream_efficiency_data,
                       ExtinctionEfficiencyFactor_tmp);
        input_tab.Read(stream_efficiency_data,
                       PhaseFunctionTable_tmp);
        for (i = 0; i < tabulation_index_real; i++)
          for (j = 0; j < tabulation_index_imaginary; j++)
            for (l = 0; l < Ndiameter; l++)
              {
                AbsorptionEfficiencyFactor(k, i, j, l) =
                  AbsorptionEfficiencyFactor_tmp(i, j, l);
                ExtinctionEfficiencyFactor(k, i, j, l) =
                  ExtinctionEfficiencyFactor_tmp(i, j, l);
                for (ld = 0; ld < 8; ld++)
                  PhaseFunctionTable(k, i, j, l, ld)
                    = PhaseFunctionTable_tmp(i, j, l, ld);
              }
        stream_efficiency_data.close();
        cout << "done." << endl;
      }
  }


  //! Compute optical properties (extinction and absorbtion efficiency
  //  factors, single scattering albedo and optical thickness and
  //  the 8 first term of the phase function). This fonction make
  //  use of tabulation for extinction and absorbtion efficiency factors
  //  as well as phase function, giving the
  //  knowledge of the aerosol refractive index and its diameter.
  /*!
    \param SingleScatteringAlbedo albedo (in [0,1]).
    \param MeanExtinctionEfficiencyFactor extinction efficiency factor.
    \param MeanAbsorbtionEfficiencyFactor absorbtion efficiency factor.
    \param OpticalThickness opticalthickness (vertical integration of
    \                                     the extinction coefficient).
    \param PhaseFunction 8 first term of the phase function.
    \param GridZ_interf_out values of the vertical grid nodes (in m).
    \param GridIndexReal, GridIndexImaginary, GridDiameter values of
    \             input parameters for optical properties tabulation
    \             (aerosol refractive indexes and diameter).
    \param AbsorptionEfficiencyFactor, ExtinctionEfficiencyFactor,
    \    PhaseFunctionTable values of tabulated optical properties
    \    corresponding to previous input parameters.
    \param ConcentrationAerosol aerosol specie concentrations
    \                                        (in microg/m^3).
    \param ConcentrationWater aerosol water concentration in
    \                                            microg/m^3).
    \param PureSpeciesIndexReal, PureSpeciesIndexImaginary
    \                     refractive index of pure specie.
    \param WaterIndexReal, WaterIndexImaginary refractive index
    \                                             of pure water.
    \param dry_diameter_computed, Wet_diameter aerosol wet
    \                                    and dry diameters.
    \param dry_aerosol_density aerosol density.
    \param black_carbon_index index of the black carbon specie in
    \                              the Data ConcentrationAerosol.
    \param Nbins_in_simulation .
    \param option_black_carbon_treatment, option_wet_index,
    \   option_well_mixed_index option for optic treatment.
  */
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
                                  int option_well_mixed_index)
  {
    T wet_diameter;
    T number_part;
    int i, j, k, h, is, l, ld;
    int Nwavelength,  Ny_out, Nz_out, Nx_out;

    const T pi = 3.14159265358979323846264;

    Nwavelength = SingleScatteringAlbedo.GetLength(0);
    Nz_out = SingleScatteringAlbedo.GetLength(1);
    Ny_out = SingleScatteringAlbedo.GetLength(2);
    Nx_out = SingleScatteringAlbedo.GetLength(3);

    int Nspecies(ConcentrationAerosol.GetLength(0));
    int Nbins_in_simulation(ConcentrationAerosol.GetLength(1));

    RegularGrid<T> GridSpeciesNoBlackCarbon(Nspecies - 1);
    RegularGrid<T> GridSpecies(Nspecies);
    RegularGrid<T> GridLegendre(8);

    for (k = 0; k < Nwavelength; k++)
      {
        Data<T, 3> AbsorptionEfficiencyFactor_tmp(GridIndexReal,
                                                  GridIndexImaginary,
                                                  GridDiameter);
        Data<T, 3> ExtinctionEfficiencyFactor_tmp(GridIndexReal,
                                                  GridIndexImaginary,
                                                  GridDiameter);
        Data<T, 3> PhaseFunctionTable_tmp(GridIndexReal,
                                          GridIndexImaginary,
                                          GridDiameter);
        Data<T, 3> PhaseFunctionTable_tmp2(GridIndexReal,
                                           GridIndexImaginary,
                                           GridDiameter);
        Data<T, 3> PhaseFunctionTable_tmp3(GridIndexReal,
                                           GridIndexImaginary,
                                           GridDiameter);
        Data<T, 3> PhaseFunctionTable_tmp4(GridIndexReal,
                                           GridIndexImaginary,
                                           GridDiameter);
        Data<T, 3> PhaseFunctionTable_tmp5(GridIndexReal,
                                           GridIndexImaginary,
                                           GridDiameter);
        Data<T, 3> PhaseFunctionTable_tmp6(GridIndexReal,
                                           GridIndexImaginary,
                                           GridDiameter);
        Data<T, 3> PhaseFunctionTable_tmp7(GridIndexReal,
                                           GridIndexImaginary,
                                           GridDiameter);
        Data<T, 3> PhaseFunctionTable_tmp8(GridIndexReal,
                                           GridIndexImaginary,
                                           GridDiameter);
        AbsorptionEfficiencyFactor_tmp.SubData(AbsorptionEfficiencyFactor,
                                               k, Range::all(),
                                               Range::all(), Range::all());
        ExtinctionEfficiencyFactor_tmp.SubData(ExtinctionEfficiencyFactor,
                                               k, Range::all(),
                                               Range::all(), Range::all());
        PhaseFunctionTable_tmp.SubData(PhaseFunctionTable,
                                       k, Range::all(), Range::all(),
                                       Range::all(), 0);
        PhaseFunctionTable_tmp2.SubData(PhaseFunctionTable,
                                        k, Range::all(), Range::all(),
                                        Range::all(), 1);
        PhaseFunctionTable_tmp3.SubData(PhaseFunctionTable,
                                        k, Range::all(), Range::all(),
                                        Range::all(), 2);
        PhaseFunctionTable_tmp4.SubData(PhaseFunctionTable,
                                        k, Range::all(), Range::all(),
                                        Range::all(), 3);
        PhaseFunctionTable_tmp5.SubData(PhaseFunctionTable,
                                        k, Range::all(), Range::all(),
                                        Range::all(), 4);
        PhaseFunctionTable_tmp6.SubData(PhaseFunctionTable,
                                        k, Range::all(), Range::all(),
                                        Range::all(), 5);
        PhaseFunctionTable_tmp7.SubData(PhaseFunctionTable,
                                        k, Range::all(), Range::all(),
                                        Range::all(), 6);
        PhaseFunctionTable_tmp8.SubData(PhaseFunctionTable,
                                        k, Range::all(), Range::all(),
                                        Range::all(), 7);

        for (j = 0; j < Ny_out; j++)
          {
            for (i = 0; i < Nx_out; i++)
              {
                for (h = 0; h < Nz_out; h++)
                  {
                    T coeff_extinction(0.0);
                    T coeff_extinction_tmp(0.0);
                    T coeff_absorption(0.0);
                    T tot_number_part(0.0);
                    T surface_mean(0.0);
                    Data<T, 1> phase_function_tmp(8);
                    Data<T, 1> phase_function(8);
                    for (ld = 0; ld < 8; ld++)
                      phase_function(ld) = 0.0;

                    for (is = 0; is < Nbins_in_simulation; is++)
                      {
                        Data<T, 1> index_real(GridSpecies);
                        Data<T, 1> index_imaginary(GridSpecies);
                        T index_dry_real, index_dry_imaginary;
                        T index_wet_real, index_wet_imaginary;
                        Data<T, 1> concentration(GridSpecies);

                        wet_diameter = Wet_diameter(is, h, j, i);

                        // Black carbon treated as the other components.
                        if (option_black_carbon_treatment == 1)
                          {
                            for (l = 0; l < Nspecies; l++)
                              {
                                concentration(l) =
                                  ConcentrationAerosol(l, is, h, j, i);
                                index_real(l) = PureSpeciesIndexReal(l, k);
                                index_imaginary(l) =
                                  PureSpeciesIndexImaginary(l, k);
                              }


                            if (option_wet_index == 1) // species treats
                              // apart from water -
                              // use Hanel's relation to calcul Wet ACRI
                              {
                                if (option_well_mixed_index == 1)
                                  // Chemical formula.
                                  compute_refractive_index
                                    (concentration,
                                     index_real,
                                     index_imaginary,
                                     index_dry_real,
                                     index_dry_imaginary);
                                else if (option_well_mixed_index == 2)
                                  // Lorentz-Lorenz.
                                  compute_refractive_index_Lorentz_Lorenz
                                    (concentration,
                                     index_real,
                                     index_imaginary,
                                     index_dry_real,
                                     index_dry_imaginary);

                                compute_Hanel_index
                                  (index_dry_real,
                                   index_dry_imaginary,
                                   WaterIndexReal(k),
                                   WaterIndexImaginary(k),
                                   dry_diameter_computed(is),
                                   wet_diameter,
                                   index_wet_real,
                                   index_wet_imaginary);
                              }
                            else // option_wet_index == 2 - treats water
                              // as othersspecies for ACRI calculation
                              {
                                Data<T, 1> concentration_tmp(Nspecies + 1);
                                Data<T, 1> index_tmp_real(Nspecies + 1);
                                Data<T, 1>
                                  index_tmp_imaginary(Nspecies + 1);
                                for (l = 0; l < Nspecies; l++)
                                  {
                                    concentration_tmp(l) =
                                      concentration(l);
                                    index_tmp_real(l) =
                                      index_real(l);
                                    index_tmp_imaginary(l) =
                                      index_imaginary(l);
                                  }
                                concentration_tmp(Nspecies) =
                                  ConcentrationWater(is, h, j, i);
                                index_tmp_real(Nspecies)
                                  = WaterIndexReal(k);
                                index_tmp_imaginary(Nspecies) =
                                  WaterIndexImaginary(k);

                                if (option_well_mixed_index == 1)
                                  // Chemical formula.
                                  compute_refractive_index
                                    (concentration_tmp,
                                     index_tmp_real,
                                     index_tmp_imaginary,
                                     index_wet_real,
                                     index_wet_imaginary);
                                else if (option_well_mixed_index == 2)
                                  // Lorentz-Lorenz.
                                  compute_refractive_index_Lorentz_Lorenz
                                    (concentration_tmp,
                                     index_tmp_real,
                                     index_tmp_imaginary,
                                     index_wet_real,
                                     index_wet_imaginary);
                              }
                          }
                        else if (option_black_carbon_treatment == 2)
                          // Black carbon core treatement.
                          // Water is treat as other species for
                          // ACRI calculation and
                          // Maxwell-Garnett dipole formula is used.
                          // Wet diameter is computed for extinction
                          // coefficient calculation using model water
                          // if water < 1.e4)
                          {
                            T index_wet_no_black_carbon_real;
                            T index_wet_no_black_carbon_imaginary;
                            Data<T, 1> concentration_no_black_carbon
                              (GridSpeciesNoBlackCarbon);
                            Data<T, 1> concentration_water_no_black_carbon
                              (GridSpecies);
                            T concentration_black_carbon =
                              ConcentrationAerosol(black_carbon_index, is,
                                                   h, j, i);
                            int ll = 0;
                            for (l = 0; l < Nspecies; l++)
                              {
                                concentration(l) =
                                  ConcentrationAerosol(l, is, h, j, i);
                                if (l != black_carbon_index)
                                  {
                                    concentration_no_black_carbon(ll) =
                                      ConcentrationAerosol(l, is, h, j, i);
                                    concentration_water_no_black_carbon(ll)
                                      = ConcentrationAerosol(l, is,
                                                             h, j, i);
                                    index_real(ll) =
                                      PureSpeciesIndexReal(l, k);
                                    index_imaginary(ll) =
                                      PureSpeciesIndexImaginary(l, k);
                                    ll += 1;
                                  }
                              }

                            concentration_water_no_black_carbon
                              (Nspecies - 1) =
                              ConcentrationWater(is, h, j, i);

                            index_real(Nspecies - 1) = WaterIndexReal(k);
                            index_imaginary(Nspecies - 1) =
                              WaterIndexImaginary(k);

                            // Solution without black carbon.
                            if (option_well_mixed_index == 1)
                              // Chemical formula.
                              compute_refractive_index
                                (concentration_water_no_black_carbon,
                                 index_real,
                                 index_imaginary,
                                 index_wet_no_black_carbon_real,
                                 index_wet_no_black_carbon_imaginary);
                            else if (option_well_mixed_index == 2)
                              // Lorentz-Lorenz.
                              compute_refractive_index_Lorentz_Lorenz
                                (concentration_water_no_black_carbon,
                                 index_real,
                                 index_imaginary,
                                 index_wet_no_black_carbon_real,
                                 index_wet_no_black_carbon_imaginary);

                            compute_refractive_index_Maxwell_Garnet
                              (concentration_black_carbon,
                               concentration_water_no_black_carbon.Sum(),
                               PureSpeciesIndexReal(black_carbon_index, k),
                               PureSpeciesIndexImaginary(black_carbon_index
                                                         , k),
                               index_wet_no_black_carbon_real,
                               index_wet_no_black_carbon_imaginary,
                               index_wet_real, index_wet_imaginary);
                          }

                        // CALCUL OF AOT AND SSA

                        // 1)  efficiency from tabulation
                        RegularGrid<T> grid_tmp_index_real(1);
                        RegularGrid<T> grid_tmp_index_imaginary(1);
                        RegularGrid<T> grid_tmp_diameter(1);
                        grid_tmp_index_real(0) = index_wet_real;
                        grid_tmp_index_imaginary(0) = index_wet_imaginary;
                        grid_tmp_diameter(0) = wet_diameter;
                        Data<T, 3> data_absorption_efficiency_factor
                          (grid_tmp_index_real, grid_tmp_index_imaginary,
                           grid_tmp_diameter);
                        Data<T, 3> data_extinction_efficiency_factor
                          (grid_tmp_index_real, grid_tmp_index_imaginary,
                           grid_tmp_diameter);
                        Data<T, 3> data_phase_function
                          (grid_tmp_index_real,
                           grid_tmp_index_imaginary,
                           grid_tmp_diameter);
                        T efficiency_absorption_factor;

                        LinearInterpolationRegular
                          (AbsorptionEfficiencyFactor_tmp,
                           data_absorption_efficiency_factor);

                        LinearInterpolationRegular
                          (ExtinctionEfficiencyFactor_tmp,
                           data_extinction_efficiency_factor);

                        efficiency_absorption_factor =
                          (data_extinction_efficiency_factor(0, 0, 0) -
                           data_absorption_efficiency_factor(0, 0, 0));

                        LinearInterpolationRegular
                          (PhaseFunctionTable_tmp,
                           data_phase_function);
                        phase_function_tmp(0) = data_phase_function(0, 0, 0);
                        LinearInterpolationRegular
                          (PhaseFunctionTable_tmp2,
                           data_phase_function);
                        phase_function_tmp(1) = data_phase_function(0, 0, 0);
                        LinearInterpolationRegular
                          (PhaseFunctionTable_tmp3,
                           data_phase_function);
                        phase_function_tmp(2) = data_phase_function(0, 0, 0);
                        LinearInterpolationRegular
                          (PhaseFunctionTable_tmp4,
                           data_phase_function);
                        phase_function_tmp(3) = data_phase_function(0, 0, 0);
                        LinearInterpolationRegular
                          (PhaseFunctionTable_tmp5,
                           data_phase_function);
                        phase_function_tmp(4) = data_phase_function(0, 0, 0);
                        LinearInterpolationRegular
                          (PhaseFunctionTable_tmp6,
                           data_phase_function);
                        phase_function_tmp(5) = data_phase_function(0, 0, 0);
                        LinearInterpolationRegular
                          (PhaseFunctionTable_tmp7,
                           data_phase_function);
                        phase_function_tmp(6) = data_phase_function(0, 0, 0);
                        LinearInterpolationRegular
                          (PhaseFunctionTable_tmp8,
                           data_phase_function);
                        phase_function_tmp(7) = data_phase_function(0, 0, 0);

                        // 2) Calcul OF AOT, SSA
                        if (concentration.Sum() > 0.0)
                          {
                            number_part = 6 * concentration.Sum()
                              / (pi * dry_aerosol_density
                                 * dry_diameter_computed(is)
                                 * dry_diameter_computed(is)
                                 * dry_diameter_computed(is));

                            tot_number_part += number_part;

                            surface_mean +=  pi / 4
                              * wet_diameter * wet_diameter
                              * number_part;

                            coeff_extinction_tmp = pi / 4
                              * wet_diameter * wet_diameter
                              * data_extinction_efficiency_factor(0, 0, 0)
                              * number_part;
                            coeff_extinction += coeff_extinction_tmp;

                            coeff_absorption += pi / 4
                              * wet_diameter * wet_diameter
                              * efficiency_absorption_factor
                              * number_part;

                            for (ld = 0; ld < 8; ld++)
                              phase_function(ld) += coeff_extinction_tmp *
                                phase_function_tmp(ld);
                          }
                      }

                    surface_mean  = surface_mean / tot_number_part;
                    OpticalThickness(k, h, j, i) = coeff_extinction
                      * (GridZ_interf_out(h + 1) - GridZ_interf_out(h));
                    for (ld = 0; ld < 8; ld++)
                      PhaseFunction(k,  h, j, i, ld) = phase_function(ld) /
                        coeff_extinction;
                    SingleScatteringAlbedo(k, h, j, i) =
                      coeff_extinction
                      / (coeff_extinction
                         + coeff_absorption);
                    MeanExtinctionEfficiencyFactor(k, h, j, i) =
                      coeff_extinction / surface_mean;
                    MeanAbsorbtionEfficiencyFactor(k, h, j, i) =
                      coeff_absorption / surface_mean;
                  }
              }
          }
      }
  }


} //namespace Polyphemus

#define COMMON_FILE_COMMON_CXX
#endif
