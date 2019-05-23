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


#ifndef ATMODATA_FILE_FORMAT_HXX

namespace AtmoData
{

  //! Input/ouput class to read files in CSV format (Airparif).
  class FormatCSV: public Format
  {

  protected:

  public:
    FormatCSV()  throw();
    ~FormatCSV()  throw();

    // Data.

    template < class TD, int N, class TG,
               class TS, class TGS >
    void Read(string FileName, Data<TD, N, TG>& D,
              Data<TS, 1, TGS>& S) const;
    template < class TD, int N, class TG,
               class TS, class TGS >
    void Read(ifstream& FileStream, Data<TD, N, TG>& D,
              Data<TS, 1, TGS>& S) const;

    // Array.

    template < class TA, int N,
               class TS, class TGS >
    void Read(string FileName, Array<TA, N>& A,
              Data<TS, 1, TGS>& S) const;
    template < class TA, int N,
               class TS, class TGS >
    void Read(ifstream& FileStream, Array<TA, N>& A,
              Data<TS, 1, TGS>& S) const;

  };

  //! Input/ouput class to read files in binary format at ECMWF.
  template<class T>
  class FormatECMWF: public Format
  {

  protected:
    //! Date.
    int date_;

  public:
    FormatECMWF()  throw();
    FormatECMWF(int date)  throw();
    ~FormatECMWF()  throw();

    void SetDate(int date);
    int GetDate() const;

    // Data.

    template<class TD, int N, class TG>
    void Read(string FileName, Data<TD, N, TG>& D) const;
    template<class TD, int N, class TG>
    void Read(ifstream& FileStream, Data<TD, N, TG>& D) const;

    // Array.

    template<class TA, int N>
    void Read(string FileName, Array<TA, N>& A) const;
    template<int N>
    void Read(ifstream& FileStream, Array<T, N>& A) const;
    template<class TA, int N>
    void Read(ifstream& FileStream, Array<TA, N>& A) const;

  };

  // For MM5 sub-headers.
  class MM5SubHeader;

  //! Input/ouput class to read MM5 files (version 3).
  class FormatMM5: public Format
  {

  protected:

  public:
    FormatMM5()  throw();
    ~FormatMM5()  throw();

    // Flag.
    int ReadFlag(ifstream& FileStream) const;

    // Big header.
    void ReadBigHeader(string FileName,
                       Array<int, 2>& BHI, Array<float, 2>& BHR,
                       Array<string, 2>& BHIC, Array<string, 2>& BHRC) const;
    void ReadBigHeader(ifstream& FileStream,
                       Array<int, 2>& BHI, Array<float, 2>& BHR,
                       Array<string, 2>& BHIC, Array<string, 2>& BHRC) const;
    void ReadBigHeader(ifstream& FileStream) const;

    // Sub-header.
    void ReadSubHeader(ifstream& FileStream, MM5SubHeader& SH) const;
    void ReadSubHeader(ifstream& FileStream) const;

    // Field.
    template <int N, class TG>
    void ReadWholeField(string FileName, string FieldName,
                        Data<float, N, TG>& A) const;
    template <int N>
    void ReadWholeField(string FileName, string FieldName,
                        Array<float, N>& A) const;
    template <int N, class TG>
    void ReadWholeField(ifstream& FileStream, string FieldName,
                        Data<float, N, TG>& A) const;
    template <int N>
    void ReadWholeField(ifstream& FileStream, string FieldName,
                        Array<float, N>& A) const;
    template <int N>
    void ReadField(ifstream& FileStream, bool cross,
                   Array<float, N>& A) const;
    template <int N, class TG>
    void ReadField(ifstream& FileStream, Data<float, N, TG>& A) const;
    template <int N>
    void ReadField(ifstream& FileStream, MM5SubHeader& SH,
                   Array<float, N>& A) const;
    template <int N>
    void ReadField(ifstream& FileStream, Array<float, N>& A) const;
    void ReadField(ifstream& FileStream) const;

  };

  //! For MM5 sub-headers.
  class MM5SubHeader
  {

  public:
    //! Dimension of the field.
    int ndim;
    //! Starting indices of the field array.
    Array<int, 1> start_index;
    //! Endiing indices of the field array.
    Array<int, 1> end_index;
    //! Integration or forecast time for the field.
    float xtime;
    //! Field at dot or cross point (character C or D).
    string staggering;
    //! Order of the field array dimension.
    string ordering;
    //! Current date.
    string current_date;
    //! Name of the field.
    string name;
    //! Unit of the field.
    string unit;
    //! Field description.
    string description;

  public:
    MM5SubHeader()  throw();
    MM5SubHeader(const MM5SubHeader&)  throw();
    ~MM5SubHeader()  throw();

    void Init();
    MM5SubHeader& operator=(MM5SubHeader&);
    string GetCurrentDate();

  };

}  // namespace AtmoData.


#define ATMODATA_FILE_FORMAT_HXX
#endif
