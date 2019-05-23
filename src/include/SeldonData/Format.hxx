// Copyright (C) 2003-2007, INRIA
// Author(s): Vivien Mallet
//
// This file is part of SeldonData library, used for data processing.
//
// SeldonData is free software; you can redistribute it and/or modify it under
// the terms of the GNU General Public License as published by the Free
// Software Foundation; either version 2 of the License, or (at your option)
// any later version.
//
// SeldonData is distributed in the hope that it will be useful, but WITHOUT
// ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
// FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for
// more details.
//
// For more information, visit the SeldonData home page:
//      http://vivienmallet.net/lib/seldondata/

#ifndef FILE_SELDONDATA_FORMAT_HXX

#include <stdio.h>
#include <iostream>
using std::cout;
using std::endl;
#include <string>
#include <fstream>
using namespace std;

#ifdef SELDONDATA_WITH_NETCDF
#include "netcdfcpp.h"
#endif

#ifdef SELDONDATA_WITH_GRIB
#include "decode_grib.cpp"
#include "grib_api.h"
extern "C"
{
  void grib_multi_support_on(grib_context* c);
  grib_handle* grib_handle_new_from_file(grib_context* c, FILE* f, int* error);
  int grib_get_long(grib_handle* h, const char* key, long* value);
  int grib_get_double_array(grib_handle* h, const char* key, double* vals,
			    size_t *length);
  int grib_handle_delete(grib_handle* h);
}
#endif

namespace SeldonData
{

  //! Base class for input/output classes.
  class Format
  {

  protected:

  public:
    Format()  throw();
    ~Format()  throw();

  };


  //! Input/ouput class to read binary files.
  template<class T>
  class FormatBinary: public Format
  {

  protected:

  public:
    FormatBinary()  throw();
    ~FormatBinary()  throw();

    // Grid.

    template<class TG>
    void Read(string FileName, RegularGrid<TG>& G) const;
    template<class TG>
    void Read(ExtStream& FileStream, RegularGrid<TG>& G) const;
    template<class TG>
    void Read(ifstream& FileStream, RegularGrid<TG>& G) const;
    template<class TG, int N>
    void Read(string FileName, GeneralGrid<TG, N>& G) const;
    template<class TG, int N>
    void Read(ExtStream& FileStream, GeneralGrid<TG, N>& G) const;
    template<class TG, int N>
    void Read(ifstream& FileStream, GeneralGrid<TG, N>& G) const;

    template<class TG>
    void Write(RegularGrid<TG>& G, string FileName) const;
    template<class TG>
    void Write(RegularGrid<TG>& G, ofstream& FileStream) const;
    template<class TG, int N>
    void Write(GeneralGrid<TG, N>& G, string FileName) const;
    template<class TG, int N>
    void Write(GeneralGrid<TG, N>& G, ofstream& FileStream) const;

    template<class TG>
    void Append(RegularGrid<TG>& G, string FileName) const;
    template<class TG, int N>
    void Append(GeneralGrid<TG, N>& G, string FileName) const;

    // Data.

    template<class TD, int N, class TG>
    void Read(string FileName, Data<TD, N, TG>& D) const;
    template<class TD, int N, class TG>
    void Read(ExtStream& FileStream, Data<TD, N, TG>& D) const;
    template<class TD, int N, class TG>
    void Read(ifstream& FileStream, Data<TD, N, TG>& D) const;

    template<class TD, int N, class TG>
    void ReadSteps(string FileName, int steps, Data<TD, N, TG>& D) const;
    template<class TD, int N, class TG>
    void ReadRecord(string FileName, int steps, Data<TD, N, TG>& D) const;

    template<class TD, int N, class TG>
    void Write(Data<TD, N, TG>& D, string FileName) const;
    template<class TD, int N, class TG>
    void Write(Data<TD, N, TG>& D, ofstream& FileStream) const;

    template<class TD, int N, class TG>
    void Append(Data<TD, N, TG>& D, string FileName) const;

    // Array.

    template<class TA, int N>
    void Read(string FileName, Array<TA, N>& A) const;
    template<int N>
    void Read(ExtStream& FileStream, Array<T, N>& A) const;
    template<int N>
    void Read(ifstream& FileStream, Array<T, N>& A) const;
    template<class TA, int N>
    void Read(ExtStream& FileStream, Array<TA, N>& A) const;
    template<class TA, int N>
    void Read(ifstream& FileStream, Array<TA, N>& A) const;

    template<class TA, int N>
    void ReadSteps(string FileName, int steps, Array<TA, N>& A) const;
    template<class TA, int N>
    void ReadRecord(string FileName, int steps, Array<TA, N>& A) const;

    template<class TA, int N>
    void Write(Array<TA, N>& A, string FileName) const;
    template<int N>
    void Write(Array<T, N>& A, ofstream& FileStream) const;
    template<class TA, int N>
    void Write(Array<TA, N>& A, ofstream& FileStream) const;

    template<class TA, int N>
    void Append(Array<TA, N>& A, string FileName) const;

  };


  //! Input/ouput class to read text files.
  class FormatText: public Format
  {

  protected:
    string separator_;
    fstream::fmtflags flags_;
    streamsize precision_;
    streamsize width_;

  public:
    FormatText()  throw();
    FormatText(string separator)  throw();
    FormatText(fstream::fmtflags flags, string separator = "\t\t")  throw();
    FormatText(fstream::fmtflags flags, streamsize precision,
               streamsize width = -1, string separator = "\t\t")  throw();
    ~FormatText()  throw();

    void SetSeparator(string separator);
    void SetFlags(ofstream::fmtflags flags);
    void SetPrecision(streamsize precision);
    void SetWidth(streamsize width);

    // Grid.

    template<class TG>
    void Read(string FileName, RegularGrid<TG>& G) const;
    template<class TG>
    void Read(ifstream& FileStream, RegularGrid<TG>& G) const;
    template<class TG, int N>
    void Read(string FileName, GeneralGrid<TG, N>& G) const;
    template<class TG, int N>
    void Read(ifstream& FileStream, GeneralGrid<TG, N>& G) const;

    template<class TG>
    void Write(RegularGrid<TG>& G, string FileName) const;
    template<class TG>
    void Write(RegularGrid<TG>& G, ofstream& FileStream) const;
    template<class TG, int N>
    void Write(GeneralGrid<TG, N>& G, string FileName) const;
    template<class TG, int N>
    void Write(GeneralGrid<TG, N>& G, ofstream& FileStream) const;

    // Data.

    template<class TD, int N, class TG>
    void Read(string FileName, Data<TD, N, TG>& D) const;
    template<class TD, int N, class TG>
    void Read(ifstream& FileStream, Data<TD, N, TG>& D) const;

    template<class TD, int N, class TG>
    void Write(Data<TD, N, TG>& D, string FileName) const;
    template<class TD, int N, class TG>
    void Write(Data<TD, N, TG>& D, ofstream& FileStream) const;

    // Array.

    template<class TA, int N>
    void Read(string FileName, Array<TA, N>& A) const;
    template<class TA, int N>
    void Read(ifstream& FileStream, Array<TA, N>& A) const;
    template<class TA>
    void Read(ifstream& FileStream, Array<TA, 1>& A) const;

    template<class TA, int N>
    void Write(Array<TA, N>& A, string FileName) const;
    template<class TA, int N>
    void Write(Array<TA, N>& A, ofstream& FileStream) const;

  };


  //! Input/ouput class to read formatted text files.
  class FormatFormattedText: public Format
  {

  protected:
    //! Description of the file format.
    string format_;
    //! Characters that denote a comment line.
    string comments_;
    //! Characters considered as delimiters.
    string delimiters_;
    //! First vector describing the format.
    vector<string> info_str;
    //! Second vector describing the format.
    vector<int> info_nb0;
    //! Third vector describing the format.
    vector<int> info_nb1;

  private:
    void SetVectors();
    void SkipMarkup(ExtStream&, streampos pos, int) const;
    template <class T>
    int ReadMarkup(ExtStream&, streampos pos, int, T*, int) const;

  public:
    FormatFormattedText(string format,
                        string comments = "#%",
                        string delimiters = " \t:;,|\n");
    ~FormatFormattedText();

    string GetFormat() const;
    string GetDelimiters() const;
    string GetComments() const;

    void SetFormat(string format);
    void SetDelimiters(string delimiters);
    void SetComments(string comments);

    // Grid.

    template<class TG>
    void Read(string FileName, string extract, RegularGrid<TG>& G) const;
    template<class TG>
    void Read(ExtStream& FileStream, string extract, RegularGrid<TG>& G)
      const;
    template<class TG, int N>
    void Read(string FileName, string extract, GeneralGrid<TG, N>& G) const;
    template<class TG, int N>
    void Read(ExtStream& FileStream, string extract, GeneralGrid<TG, N>& G)
      const;

    // Data.

    template<class TD, int N, class TG>
    void Read(string FileName, string extract, Data<TD, N, TG>& D) const;
    template<class TD, int N, class TG>
    void Read(ExtStream& FileStream, string extract, Data<TD, N, TG>& D)
      const;

    // Array.

    template<class TA, int N>
    void Read(string FileName, string extract, Array<TA, N>& A) const;
    template<class TA, int N>
    void Read(ExtStream& FileStream, string extract, Array<TA, N>& A) const;

  };


#ifdef SELDONDATA_WITH_NETCDF
  //! Input/ouput class to read netCDF files.
  template<class T>
  class FormatNetCDF: public Format
  {

  protected:

  public:
    FormatNetCDF()  throw();
    ~FormatNetCDF()  throw();
    // For variables.
    // Grid.

    template<class TG>
    void Read(string FileName, string variable, RegularGrid<TG>& G) const;
    template<class TG, int N>
    void Read(string FileName, string variable, GeneralGrid<TG, N>& G) const;

    // Data.

    template<class TD, int N, class TG>
    void Read(string FileName, string variable, Data<TD, N, TG>& D) const;

    // Array.

    template<class TA, int N>
    void Read(string FileName, string variable, Array<TA, N>& A) const;

    // For dimensions.
    void ReadDimension(string FileName, string variable, int dim_num,
                       int& dim_value) const;

    // For attibutes.
    void ReadAttribute(string FileName, string attribute, float& value) const;
    void ReadAttribute(string FileName, string attribute, int& value) const;

  };
#endif


#ifdef SELDONDATA_WITH_GRIB
  //! Input/ouput class to read Grib files.
  class FormatGrib: public Format
  {

  protected:

  public:
    FormatGrib()  throw();
    ~FormatGrib()  throw();

    // Grid.

    template<class TG>
    void Read(string FileName, int variable, RegularGrid<TG>& G) const;
    template<class TG, int N>
    void Read(string FileName, int variable, GeneralGrid<TG, N>& G) const;

    // Data.

    template<class TD, int N, class TG>
    void Read(string FileName, int variable, Data<TD, N, TG>& D) const;

    // Array.

    template<class TA, int N>
    void Read(string FileName, int variable, Array<TA, N>& A) const;

  };

  //! Input/ouput class to read Grib2 files.
  class FormatGrib2: public Format
  {

  protected:

  public:
    FormatGrib2()  throw();
    ~FormatGrib2()  throw();

    // Grid.

    template<class TG>
    void Read(string FileName, int discipline, int parameterCategory,
	      int parameterNumber, RegularGrid<TG>& G) const;
    template<class TG, int N>
    void Read(string FileName, int discipline, int parameterCategory,
	      int parameterNumber, GeneralGrid<TG, N>& G) const;

    // Data.

    template<class TD, int N, class TG>
    void Read(string FileName, int discipline, int parameterCategory,
	      int parameterNumber, Data<TD, N, TG>& D) const;

    // Array.

    template<class TA, int N>
    void Read(string FileName, int discipline, int parameterCategory,
	      int parameterNumber, Array<TA, N>& A) const;
  };
#endif


  //! Input/ouput class to read files in Chimere format.
  class FormatChimere: public Format
  {

  protected:
    int date_;

  public:
    FormatChimere()  throw();
    FormatChimere(int date)  throw();
    ~FormatChimere()  throw();

    void SetDate(int date);
    int GetDate() const;

    // Data.

    template<class TD, int N, class TG>
    void Read(string FileName, Data<TD, N, TG>& D,
              int nb_lines = -1) const;
    template<class TD, int N, class TG>
    void Read(ifstream& FileStream, Data<TD, N, TG>& D,
              int nb_lines = -1) const;

    // Array.

    template<class TA, int N>
    void Read(string FileName, Array<TA, N>& A,
              int nb_lines = -1) const;
    template<class TA, int N>
    void Read(ifstream& FileStream, Array<TA, N>& A,
              int nb_lines = -1) const;

  };

}  // namespace SeldonData.

#define FILE_SELDONDATA_FORMAT_HXX
#endif
