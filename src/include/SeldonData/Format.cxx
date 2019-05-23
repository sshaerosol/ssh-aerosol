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

#ifndef FILE_SELDONDATA_FORMAT_CXX

#include "Format.hxx"

#ifndef CONST_SELDONDATA_BUFFER_SIZE
#define CONST_SELDONDATA_BUFFER_SIZE 16384
#endif

namespace SeldonData
{


  ////////////
  // FORMAT //
  ////////////

  //! Default constructor.
  Format::Format()  throw()
  {
  }

  //! Destructor.
  Format::~Format()  throw()
  {
  }


  //////////////////
  // FORMATBINARY //
  //////////////////

  //! Default constructor.
  template<class T>
  FormatBinary<T>::FormatBinary()  throw()
  {
  }

  //! Destructor.
  template<class T>
  FormatBinary<T>::~FormatBinary()  throw()
  {
  }

  /********/
  /* Grid */
  /********/

  //! Reads a binary file.
  template<class T>
  template<class TG>
  void FormatBinary<T>::Read(string FileName, RegularGrid<TG>& G) const
  {

    this->Read(FileName, G.GetArray());

  }

  //! Reads a binary file.
  template<class T>
  template<class TG>
  void FormatBinary<T>::Read(ExtStream& FileStream, RegularGrid<TG>& G) const
  {

    this->Read(FileStream, G.GetArray());

  }

  //! Reads a binary file.
  template<class T>
  template<class TG>
  void FormatBinary<T>::Read(ifstream& FileStream, RegularGrid<TG>& G) const
  {

    this->Read(FileStream, G.GetArray());

  }

  //! Reads a binary file.
  template<class T>
  template<class TG, int n>
  void FormatBinary<T>::Read(string FileName, GeneralGrid<TG, n>& G) const
  {

    this->Read(FileName, G.GetArray());

  }

  //! Reads a binary file.
  template<class T>
  template<class TG, int n>
  void FormatBinary<T>
  ::Read(ExtStream& FileStream, GeneralGrid<TG, n>& G) const
  {

    this->Read(FileStream, G.GetArray());

  }

  //! Reads a binary file.
  template<class T>
  template<class TG, int n>
  void FormatBinary<T>::Read(ifstream& FileStream,
                             GeneralGrid<TG, n>& G) const
  {

    this->Read(FileStream, G.GetArray());

  }

  //! Writes a binary file.
  template<class T>
  template<class TG>
  void FormatBinary<T>::Write(RegularGrid<TG>& G, string FileName) const
  {

    this->Write(G.GetArray(), FileName);

  }

  //! Writes a binary file.
  template<class T>
  template<class TG>
  void FormatBinary<T>::Write(RegularGrid<TG>& G, ofstream& FileStream) const
  {

    this->Write(G.GetArray(), FileStream);

  }

  //! Writes a binary file.
  template<class T>
  template<class TG, int n>
  void FormatBinary<T>::Write(GeneralGrid<TG, n>& G, string FileName) const
  {

    this->Write(G.GetArray(), FileName);

  }

  //! Writes a binary file.
  template<class T>
  template<class TG, int n>
  void FormatBinary<T>::Write(GeneralGrid<TG, n>& G,
                              ofstream& FileStream) const
  {

    this->Write(G.GetArray(), FileStream);

  }

  //! Appends data to a binary file.
  template<class T>
  template<class TG>
  void FormatBinary<T>::Append(RegularGrid<TG>& G, string FileName) const
  {

    this->Append(G.GetArray(), FileName);

  }

  //! Appends data to a binary file.
  template<class T>
  template<class TG, int n>
  void FormatBinary<T>::Append(GeneralGrid<TG, n>& G, string FileName) const
  {

    this->Append(G.GetArray(), FileName);

  }

  /********/
  /* Data */
  /********/

  //! Reads a binary file.
  template<class T>
  template<class TD, int N, class TG>
  void FormatBinary<T>::Read(string FileName, Data<TD, N, TG>& D) const
  {

    this->Read(FileName, D.GetArray());

  }

  //! Reads a binary file.
  template<class T>
  template<class TD, int N, class TG>
  void FormatBinary<T>::Read(ExtStream& FileStream, Data<TD, N, TG>& D) const
  {

    this->Read(FileStream, D.GetArray());

  }

  //! Reads a binary file.
  template<class T>
  template<class TD, int N, class TG>
  void FormatBinary<T>::Read(ifstream& FileStream, Data<TD, N, TG>& D) const
  {

    this->Read(FileStream, D.GetArray());

  }

  //! Reads given steps in a binary file.
  template<class T>
  template<class TD, int N, class TG>
  void FormatBinary<T>::ReadSteps(string FileName, int steps,
                                  Data<TD, N, TG>& D) const
  {

    this->ReadSteps(FileName, steps, D.GetArray());

  }

  //! Reads a given record in a binary file.
  template<class T>
  template<class TD, int N, class TG>
  void FormatBinary<T>::ReadRecord(string FileName, int record,
                                   Data<TD, N, TG>& D) const
  {

    this->ReadRecord(FileName, record, D.GetArray());

  }

  //! Writes a binary file.
  template<class T>
  template<class TD, int N, class TG>
  void FormatBinary<T>::Write(Data<TD, N, TG>& D, string FileName) const
  {

    this->Write(D.GetArray(), FileName);

  }

  //! Writes a binary file.
  template<class T>
  template<class TD, int N, class TG>
  void FormatBinary<T>::Write(Data<TD, N, TG>& D, ofstream& FileStream) const
  {

    this->Write(D.GetArray(), FileStream);

  }

  //! Appends data to a binary file.
  template<class T>
  template<class TD, int N, class TG>
  void FormatBinary<T>::Append(Data<TD, N, TG>& D, string FileName) const
  {

    this->Append(D.GetArray(), FileName);

  }

  /*********/
  /* Array */
  /*********/

  //! Reads a binary file.
  template<class T>
  template<class TA, int N>
  void FormatBinary<T>::Read(string FileName, Array<TA, N>& A) const
  {

    ExtStream FileStream(FileName);

#ifdef SELDONDATA_DEBUG_CHECK_IO
    // Checks if the file was opened.
    if (!FileStream.is_open())
      throw IOError("FormatBinary<T>::Read(string FileName, Array<TA, N>& A)",
                    "Unable to open file \"" + FileName + "\".");
#endif

    this->Read(FileStream, A);

    FileStream.close();

  }

  //! Reads a binary file.
  template<class T>
  template<int N>
  void FormatBinary<T>::Read(ExtStream& FileStream, Array<T, N>& A) const
  {

    unsigned long data_size = A.numElements() * sizeof(T);
    T* data = A.data();

#ifdef SELDONDATA_DEBUG_CHECK_IO
    // Checks if the file ready.
    if (!FileStream.good())
      throw IOError(string("FormatBinary<T>::Read(ExtStream& FileStream, ")
                    + "Array<T, N>& A)", string("Unable to read data in \"")
                    + FileStream.GetFileName() + "\".");

    // Checks file length.
    streampos position;
    position = FileStream.tellg();
    FileStream.seekg(0, ios::end);
    unsigned long file_size = FileStream.tellg() - position;

    if (position > FileStream.tellg())
      throw IOError(string("FormatBinary<T>::Read(ExtStream& FileStream,")
                    + " Array<T, N>& A)", string("Unable to read ")
                    + to_str(data_size) + string(" byte(s) in \"")
                    + FileStream.GetFileName()
                    + "\". The input stream is empty.");

    if (data_size > file_size)
      throw IOError(string("FormatBinary<T>::Read(ExtStream& FileStream,")
                    + " Array<T, N>& A)", string("Unable to read ")
                    + to_str(data_size) + string(" byte(s) in \"")
                    + FileStream.GetFileName()
                    + string("\". The input stream is only ")
                    + to_str(file_size) + " byte(s) long.");

    FileStream.seekg(position);
#endif

    FileStream.read(reinterpret_cast<char*>(data), data_size);

  }

  //! Reads a binary file.
  template<class T>
  template<int N>
  void FormatBinary<T>::Read(ifstream& FileStream, Array<T, N>& A) const
  {

    unsigned long data_size = A.numElements() * sizeof(T);
    T* data = A.data();

#ifdef SELDONDATA_DEBUG_CHECK_IO
    // Checks if the file ready.
    if (!FileStream.good())
      throw IOError("FormatBinary<T>::Read"
                    "(ifstream& FileStream, Array<T, N>& A)",
                    "File is not ready.");

    // Checks file length.
    streampos position;
    position = FileStream.tellg();
    FileStream.seekg(0, ios::end);
    unsigned long file_size = FileStream.tellg() - position;

    if (position > FileStream.tellg())
      throw IOError("FormatBinary<T>::Read"
                    "(ifstream& FileStream, Array<T, N>& A)",
                    "Unable to read " + to_str(data_size) + " byte(s)." +
                    " The input stream is empty.");

    if (data_size > file_size)
      throw IOError("FormatBinary<T>::Read"
                    "(ifstream& FileStream, Array<T, N>& A)",
                    "Unable to read " + to_str(data_size) + " byte(s)." +
                    " The input stream is only " + to_str(file_size) +
                    " byte(s) long.");

    FileStream.seekg(position);
#endif

    FileStream.read(reinterpret_cast<char*>(data), data_size);

  }

  //! Reads a binary file.
  template<class T>
  template<class TA, int N>
  void FormatBinary<T>::Read(ExtStream& FileStream, Array<TA, N>& A) const
  {

    unsigned long data_size = A.numElements() * sizeof(T);

    unsigned long length = CONST_SELDONDATA_BUFFER_SIZE / sizeof(T);
    T* data = new T[length];

    TA* data_output = A.data();

#ifdef SELDONDATA_DEBUG_CHECK_IO
    // Checks if the file ready.
    if (!FileStream.good())
      throw IOError(string("FormatBinary<T>::Read(ExtStream& FileStream,")
                    + " Array<TA, N>& A)", string("Unable to read data in \"")
                    + FileStream.GetFileName() + "\".");

    // Checks file length.
    streampos position;
    position = FileStream.tellg();
    FileStream.seekg(0, ios::end);
    unsigned long file_size = FileStream.tellg() - position;

    if (position > FileStream.tellg())
      throw IOError(string("FormatBinary<T>::Read(ExtStream& FileStream,")
                    + " Array<T, N>& A)", string("Unable to read ")
                    + to_str(data_size) + string(" byte(s) in \"")
                    + FileStream.GetFileName()
                    + "\". The input stream is empty.");

    if (data_size > file_size)
      throw IOError(string("FormatBinary<T>::Read(ExtStream& FileStream,")
                    + " Array<TA, N>& A)", string("Unable to read ")
                    + to_str(data_size) + string(" byte(s) in \"")
                    + FileStream.GetFileName()
                    + string("\". The input stream is only ")
                    + to_str(file_size) + " byte(s) long.");

    FileStream.seekg(position);
#endif

    int i = 0;
    int j = 0;
    for (i = 0; i < int(data_size / sizeof(T) / length); i++)
      {
        FileStream.read(reinterpret_cast<char*>(data), length * sizeof(T));
        for (j = 0; j < int(length); j++)
          data_output[j + i * length] = data[j];
      }

    if (data_size % (length * sizeof(T)) != 0)
      {
        FileStream.read(reinterpret_cast<char*>(data),
                        data_size - i * length * sizeof(T));
        for (j = 0; j < int((data_size % (length * sizeof(T))) / sizeof(T));
             j++)
          data_output[j + i * length] = data[j];
      }

    delete [] data;

  }

  //! Reads a binary file.
  template<class T>
  template<class TA, int N>
  void FormatBinary<T>::Read(ifstream& FileStream, Array<TA, N>& A) const
  {

    unsigned long data_size = A.numElements() * sizeof(T);

    unsigned long length = CONST_SELDONDATA_BUFFER_SIZE / sizeof(T);
    T* data = new T[length];

    TA* data_output = A.data();

#ifdef SELDONDATA_DEBUG_CHECK_IO
    // Checks if the file ready.
    if (!FileStream.good())
      throw IOError("FormatBinary<T>::Read"
                    "(ifstream& FileStream, Array<TA, N>& A)",
                    "File is not ready.");

    // Checks file length.
    streampos position;
    position = FileStream.tellg();
    FileStream.seekg(0, ios::end);
    unsigned long file_size = FileStream.tellg() - position;

    if (position > FileStream.tellg())
      throw IOError("FormatBinary<T>::Read"
                    "(ifstream& FileStream, Array<TA, N>& A)",
                    "Unable to read " + to_str(data_size) + " byte(s)." +
                    " The input stream is empty.");

    if (data_size > file_size)
      throw IOError("FormatBinary<T>::Read"
                    "(ifstream& FileStream, Array<TA, N>& A)",
                    "Unable to read " + to_str(data_size) + " byte(s)." +
                    " The input stream is only " + to_str(file_size)
                    + " byte(s) long.");

    FileStream.seekg(position);
#endif

    int i = 0;
    int j = 0;
    for (i = 0; i < int(data_size / sizeof(T) / length); i++)
      {
        FileStream.read(reinterpret_cast<char*>(data), length * sizeof(T));
        for (j = 0; j < int(length); j++)
          data_output[j + i * length] = data[j];
      }

    if (data_size % (length * sizeof(T)) != 0)
      {
        FileStream.read(reinterpret_cast<char*>(data),
                        data_size - i * length * sizeof(T));
        for (j = 0; j < int((data_size % (length * sizeof(T))) / sizeof(T));
             j++)
          data_output[j + i * length] = data[j];
      }

    delete [] data;

  }

  //! Reads given steps in a binary file.
  template<class T>
  template<class TA, int N>
  void FormatBinary<T>::ReadSteps(string FileName, int steps,
                                  Array<TA, N>& A) const
  {

    ExtStream FileStream(FileName);

#ifdef SELDONDATA_DEBUG_CHECK_IO
    // Checks if the file was opened.
    if (!FileStream.is_open())
      throw IOError("FormatBinary<T>::ReadSteps"
                    "(string FileName, int steps, Array<TA, N>& A)",
                    "Unable to open file \"" + FileName + "\".");
#endif

    size_t pos = A.numElements() / A.extent(0);
    pos *= steps * sizeof(T);
    FileStream.seekg(pos);

#ifdef SELDONDATA_DEBUG_CHECK_IO
    // Checks whether all steps were skipped.
    if (!FileStream.good())
      throw IOError("FormatBinary<T>::ReadSteps"
                    "(string FileName, int steps, Array<TA, N>& A)",
                    string("Unable to skip ") + to_str(steps)
                    + " steps in file \"" + FileName + "\".");
#endif

    this->Read(FileStream, A);

    FileStream.close();

  }

  //! Reads given steps in a binary file.
  template<class T>
  template<class TA, int N>
  void FormatBinary<T>::ReadRecord(string FileName, int record,
                                   Array<TA, N>& A) const
  {

    ExtStream FileStream(FileName);

#ifdef SELDONDATA_DEBUG_CHECK_IO
    // Checks if the file was opened.
    if (!FileStream.is_open())
      throw IOError("FormatBinary<T>::ReadRecord"
                    "(string FileName, int record, Array<TA, N>& A)",
                    "Unable to open file \"" + FileName + "\".");
#endif

    size_t pos = record * sizeof(T) * A.numElements();
    FileStream.seekg(pos);

#ifdef SELDONDATA_DEBUG_CHECK_IO
    // Checks whether all steps were skipped.
    if (!FileStream.good())
      throw IOError("FormatBinary<T>::ReadRecord"
                    "(string FileName, int record, Array<TA, N>& A)",
                    string("Unable to skip ") + to_str(record)
                    + " records in file \"" + FileName + "\".");
#endif

    this->Read(FileStream, A);

    FileStream.close();

  }

  //! Writes a binary file.
  template<class T>
  template<class TA, int N>
  void FormatBinary<T>::Write(Array<TA, N>& A, string FileName) const
  {

    ofstream FileStream;
    FileStream.open(FileName.c_str(), ofstream::binary);

#ifdef SELDONDATA_DEBUG_CHECK_IO
    // Checks if the file was opened.
    if (!FileStream.is_open())
      throw IOError("FormatBinary<T>::"
                    "Write(Array<TA, N>& A, string FileName)",
                    "Unable to open file \"" + FileName + "\".");
#endif

    this->Write(A, FileStream);

    FileStream.close();

  }

  //! Writes a binary file.
  template<class T>
  template<int N>
  void FormatBinary<T>::Write(Array<T, N>& A, ofstream& FileStream) const
  {

    unsigned long data_size = A.numElements() * sizeof(T);

    T* data = A.data();

#ifdef SELDONDATA_DEBUG_CHECK_IO
    // Checks if the file ready.
    if (!FileStream.good())
      throw IOError("FormatBinary<T>::"
                    "Write(Array<T, N>& A, ofstream& FileStream)",
                    "File is not ready.");
#endif

    FileStream.write(reinterpret_cast<char*>(data), data_size);

  }

  //! Writes a binary file.
  template<class T>
  template<class TA, int N>
  void FormatBinary<T>::Write(Array<TA, N>& A, ofstream& FileStream) const
  {

    unsigned long data_size = A.numElements() * sizeof(T);

    unsigned long length = CONST_SELDONDATA_BUFFER_SIZE / sizeof(T);
    T* data = new T[length];

    TA* data_input = A.data();

#ifdef SELDONDATA_DEBUG_CHECK_IO
    // Checks if the file ready.
    if (!FileStream.good())
      throw IOError("FormatBinary<T>::"
                    "Write(Array<TA, N>& A, ofstream& FileStream)",
                    "File is not ready.");
#endif

    int i = 0;
    int j = 0;
    for (i = 0; i < int(data_size / sizeof(T) / length); i++)
      {
        for (j = 0; j < int(length); j++)
          data[j] = data_input[j + i * length];
        FileStream.write(reinterpret_cast<char*>(data), length * sizeof(T));
      }

    if (data_size % (length * sizeof(T)) != 0)
      {
        for (j = 0; j < int((data_size % (length * sizeof(T))) / sizeof(T));
             j++)
          data[j] = data_input[j + i * length];
        FileStream.write(reinterpret_cast<char*>(data),
                         data_size - i * length * sizeof(T));
      }

    delete[] data;

  }

  //! Appends data to a binary file.
  template<class T>
  template<class TA, int N>
  void FormatBinary<T>::Append(Array<TA, N>& A, string FileName) const
  {

    ofstream FileStream;
    FileStream.open(FileName.c_str(), ofstream::binary | ios::app);

#ifdef SELDONDATA_DEBUG_CHECK_IO
    // Checks if the file was opened.
    if (!FileStream.is_open())
      throw IOError("FormatBinary<T>::"
                    "Append(Array<TA, N>& A, string FileName)",
                    "Unable to open file \"" + FileName + "\".");
#endif

    this->Write(A, FileStream);

    FileStream.close();

  }


  ////////////////
  // FORMATTEXT //
  ////////////////

  //! Default constructor.
  FormatText::FormatText()  throw()
  {
    separator_ = "\t\t";
    flags_ = fstream::scientific | fstream::skipws;
    precision_ = -1;
    width_ = -1;
  }

  //! Constructor.
  FormatText::FormatText(string separator)  throw()
  {
    separator_ = separator;
    flags_ = fstream::scientific | fstream::skipws;
    precision_ = -1;
    width_ = -1;
  }

  //! Constructor.
  FormatText::FormatText(fstream::fmtflags flags, string separator)  throw()
  {
    separator_ = separator;
    flags_ = flags | fstream::skipws;
    precision_ = -1;
    width_ = -1;
  }

  //! Constructor.
  FormatText::FormatText(fstream::fmtflags flags, streamsize precision,
                         streamsize width, string separator)  throw()
  {
    separator_ = separator;
    flags_ = flags | fstream::skipws;
    precision_ = precision;
    width_ = width;
  }

  //! Destructor.
  FormatText::~FormatText()  throw()
  {

  }

  //! Sets the separator.
  /*!
    \param separator the new separator.
  */
  void FormatText::SetSeparator(string separator)
  {
    separator_ = separator;
  }

  //! Sets format flags.
  /*!
    \param flags format flags.
  */
  void FormatText::SetFlags(fstream::fmtflags flags)
  {
    flags_ = flags;
  }

  //! Sets floating-point decimal presision.
  /*!
    \param precision floating-point decimal presision.
  */
  void FormatText::SetPrecision(streamsize precision)
  {
    precision_ = precision;
  }

  //! Sets field width.
  /*!
    \param width field width.
  */
  void FormatText::SetWidth(streamsize width)
  {
    width_ = width;
  }

  /********/
  /* Grid */
  /********/

  //! Reads a text file.
  template<class TG>
  void FormatText::Read(string FileName, RegularGrid<TG>& G) const
  {

    this->Read(FileName, G.GetArray());

  }

  //! Reads a text file.
  template<class TG>
  void FormatText::Read(ifstream& FileStream, RegularGrid<TG>& G) const
  {

    this->Read(FileStream, G.GetArray());

  }

  //! Reads a text file.
  template<class TG, int n>
  void FormatText::Read(string FileName, GeneralGrid<TG, n>& G) const
  {

    this->Read(FileName, G.GetArray());

  }

  //! Reads a text file.
  template<class TG, int n>
  void FormatText::Read(ifstream& FileStream, GeneralGrid<TG, n>& G) const
  {

    this->Read(FileStream, G.GetArray());

  }

  //! Writes a text file.
  template<class TG>
  void FormatText::Write(RegularGrid<TG>& G, string FileName) const
  {

    this->Write(G.GetArray(), FileName);

  }

  //! Writes a text file.
  template<class TG>
  void FormatText::Write(RegularGrid<TG>& G, ofstream& FileStream) const
  {

    this->Write(G.GetArray(), FileStream);

  }

  //! Writes a text file.
  template<class TG, int n>
  void FormatText::Write(GeneralGrid<TG, n>& G, string FileName) const
  {

    this->Write(G.GetArray(), FileName);

  }

  //! Writes a text file.
  template<class TG, int n>
  void FormatText::Write(GeneralGrid<TG, n>& G, ofstream& FileStream) const
  {

    this->Write(G.GetArray(), FileStream);

  }

  /********/
  /* Data */
  /********/

  //! Reads a text file.
  template<class TD, int N, class TG>
  void FormatText::Read(string FileName, Data<TD, N, TG>& D) const
  {

    this->Read(FileName, D.GetArray());

  }

  //! Reads a text file.
  template<class TD, int N, class TG>
  void FormatText::Read(ifstream& FileStream, Data<TD, N, TG>& D) const
  {

    this->Read(FileStream, D.GetArray());

  }

  //! Writes a text file.
  template<class TD, int N, class TG>
  void FormatText::Write(Data<TD, N, TG>& D, string FileName) const
  {

    this->Write(D.GetArray(), FileName);

  }

  //! Writes a text file.
  template<class TD, int N, class TG>
  void FormatText::Write(Data<TD, N, TG>& D, ofstream& FileStream) const
  {

    this->Write(D.GetArray(), FileStream);

  }

  /*********/
  /* Array */
  /*********/

  //! Reads a text file.
  template<class TA, int N>
  void FormatText::Read(string FileName, Array<TA, N>& A) const
  {

    ifstream FileStream;
    FileStream.open(FileName.c_str(), ifstream::in);

    FileStream.flags(flags_);

    if (precision_ != -1)
      FileStream.precision(precision_);

    if (width_ != -1)
      FileStream.width(width_);

#ifdef SELDONDATA_DEBUG_CHECK_IO
    // Checks if the file was opened.
    if (!FileStream.is_open())
      throw IOError("FormatText::Read(string FileName, Array<TA, N>& A)",
                    "Unable to open file \"" + FileName + "\".");
#endif

    this->Read(FileStream, A);

    FileStream.close();

  }

  //! Reads a text file.
  template<class TA, int N>
  void FormatText::Read(ifstream& FileStream, Array<TA, N>& A) const
  {

    Array<TA, 1> B(A.data(), shape(A.numElements()), neverDeleteData);
    this->Read(FileStream, B);

  }

  //! Reads a text file.
  template<class TA>
  void FormatText::Read(ifstream& FileStream, Array<TA, 1>& A) const
  {

#ifdef SELDONDATA_DEBUG_CHECK_IO
    // Checks if the file ready.
    if (!FileStream.good())
      throw IOError("FormatText::Read(ifstream& FileStream, Array<TA, N>& A)",
                    "File is not ready.");
#endif

    int nb_elements = A.numElements();
    char c;
    int i = 0;

    while ((i < nb_elements) && (FileStream.good()))
      {

        FileStream >> A(i);

        c = FileStream.peek();
        while ((FileStream.good())
               && ((c < '0') || (c > '9'))
               && (c != '.') && (c != '-')
               && (c != '+'))
          {
            FileStream.ignore(1);
            c = FileStream.peek();
          }

        i++;
      }

#ifdef SELDONDATA_DEBUG_CHECK_IO
    // Checks if all was read.
    if (i != nb_elements)
      throw IOError("FormatText::Read(ifstream& FileStream, Array<TA, N>& A)",
                    to_str(i) + " elements were read instead of "
                    + to_str(nb_elements) + ".");
#endif

  }

  //! Writes a text file.
  template<class TA, int N>
  void FormatText::Write(Array<TA, N>& A, string FileName) const
  {

    ofstream FileStream;
    FileStream.open(FileName.c_str(), ofstream::out);

    FileStream.flags(flags_);

    if (precision_ != -1)
      FileStream.precision(precision_);

    if (width_ != -1)
      FileStream.width(width_);

#ifdef SELDONDATA_DEBUG_CHECK_IO
    // Checks if the file was opened.
    if (!FileStream.is_open())
      throw IOError("FormatText::Write(Array<TA, N>& A, string FileName)",
                    "Unable to open file \"" + FileName + "\".");
#endif

    this->Write(A, FileStream);

    FileStream.close();

  }

  //! Writes a text file.
  template<class TA, int N>
  void FormatText::Write(Array<TA, N>& A, ofstream& FileStream) const
  {

#ifdef SELDONDATA_DEBUG_CHECK_IO
    // Checks if the file ready.
    if (!FileStream.good())
      throw IOError("FormatText::Write"
                    "(Array<TA, N>& A, ofstream& FileStream)",
                    "File is not ready.");
#endif

    int nb_elements = A.numElements();
    Data<TA, N> DA(A.data(), A.shape(), neverDeleteData);
    int i = 0;
    int j;
    Array<int, 1> Index(10), Length(10);

    for (j = 0; j < 10; j++)
      {
        Index(j) = 0;
        Length(j) = A.extent(j);
      }

    j = N - 1;
    while ((i < nb_elements) && (FileStream.good()))
      {

        FileStream << DA.Value(Index(0), Index(1), Index(2),
                               Index(3), Index(4), Index(5),
                               Index(6), Index(7), Index(8),
                               Index(9));

        j = N - 1;
        while ((j >= 0) && (Index(j) == Length(j) - 1))
          {
            Index(j) = 0;
            j--;
          }

        if (j != -1)
          Index(j)++;

        if ((j != N - 1) || (N == 1))
          FileStream << '\n';
        else
          FileStream << separator_;

        i++;
      }

  }


  /////////////////////////
  // FORMATFORMATTEDTEXT //
  /////////////////////////

  //! Main constructor.
  /*!
    \param format format of the file.
    \param commments characters that denote a comment line.
    \param delimiters characters used to delimit elements in the file.
  */
  FormatFormattedText::FormatFormattedText(string format,
                                           string comments,
                                           string delimiters):
    format_(format), comments_(comments), delimiters_(delimiters)
  {
    this->SetVectors();
  }

  //! Destructor.
  FormatFormattedText::~FormatFormattedText()
  {

  }

  //! Sets vectors associated with the format.
  void FormatFormattedText::SetVectors()
  {
    info_str.clear();
    info_nb0.clear();
    info_nb1.clear();

    string stmp;
    vector<string> markup, desc;

    split(format_, markup, ">");
    for (unsigned int i = 0; i < markup.size(); i++)
      {
        split(markup[i], desc, "<");
#ifdef SELDONDATA_DEBUG_CHECK_IO
        if (desc.size() == 0)
          throw IOError("FormatFormattedText::SetVectors()",
                        "Empty markup detected found.");
#endif
        split(desc[desc.size() - 1], desc, " \n\t-,|/");
#ifdef SELDONDATA_DEBUG_CHECK_IO
        if (desc.size() == 0)
          throw IOError("FormatFormattedText::SetVectors()",
                        "Empty markup found.");
#endif
        info_str.push_back(desc[0]);

        if (info_str[i] == "c")
          {

#ifdef SELDONDATA_DEBUG_CHECK_IO
            if (desc.size() != 3)
              throw IOError("FormatFormattedText::SetVectors()",
                            string("Column descriptor must be followed by"
                                   " two numbers. ")
                            + to_str(desc.size() - 1) + " elements were "
                            "provided.");
#endif

#ifdef SELDONDATA_DEBUG_CHECK_IO
            if (!is_unsigned_integer(desc[1]))
              throw IOError("FormatFormattedText::SetVectors()",
                            string("Column descriptor must be followed by two"
                                   " numbers. ")
                            + "The first element that follows 'c' is not an"
                            " unsigned integer.");
            if (!is_unsigned_integer(desc[2]))
              throw IOError("FormatFormattedText::SetVectors()",
                            string("Column descriptor must be followed by two"
                                   " numbers. ")
                            + "The second number that follows 'c' is not an"
                            " unsigned integer.");
#endif

            info_nb0.push_back(to_num<int>(desc[1]));
            info_nb1.push_back(to_num<int>(desc[2]));

          }
        else if (info_str[i] == "e")
          {

#ifdef SELDONDATA_DEBUG_CHECK_IO
            if ((desc.size() != 1) && (desc.size() != 2))
              throw IOError("FormatFormattedText::SetVectors()",
                            "Element descriptor must be followed by at most "
                            "one number.");
            if ((desc.size() == 2) && !is_unsigned_integer(desc[1]))
              throw IOError("FormatFormattedText::SetVectors()",
                            "If followed by something, element descriptor"
                            " must be followed an unsigned integer.");
#endif

            if (desc.size() == 1)
              info_nb0.push_back(1);
            else
              info_nb0.push_back(to_num<int>(desc[1]));

            info_nb1.push_back(-1);

          }
        else if (info_str[i] == "s")
          {

#ifdef SELDONDATA_DEBUG_CHECK_IO
            if ((desc.size() != 1) && (desc.size() != 2))
              throw IOError("FormatFormattedText::SetVectors()",
                            "Skip descriptor must be followed by at most one"
                            " number.");
            if ((desc.size() == 2) && !is_unsigned_integer(desc[1]))
              throw IOError("FormatFormattedText::SetVectors()",
                            "If followed by something, skip descriptor must"
                            " be followed an unsigned integer.");
#endif

            if (desc.size() == 1)
              info_nb0.push_back(1);
            else
              info_nb0.push_back(to_num<int>(desc[1]));

            info_nb1.push_back(-1);

          }
        else if (info_str[i] == "a" || info_str[i] == "A")
          {

#ifdef SELDONDATA_DEBUG_CHECK_IO
            if (desc.size() != 1)
              throw IOError("FormatFormattedText::SetVectors()",
                            "All descriptor must not be followed by"
                            " anything.");
#endif

            info_nb0.push_back(-1);
            info_nb1.push_back(-1);

          }
#ifdef SELDONDATA_DEBUG_CHECK_IO
        else
          throw IOError("FormatFormattedText::SetVectors()",
                        "Unknown delimiter.");
#endif

      }
  }

  //! Skips data.
  void FormatFormattedText::SkipMarkup(ExtStream& FileStream,
                                       streampos pos, int i) const
  {
    if (info_str[i] == "c")
      {
#ifdef SELDONDATA_DEBUG_CHECK_IO
        if (info_nb0[i] < pos)
          throw IOError("FormatFormattedText::SkipMarkup"
                        "(ExtStream& FileStream, streampos pos, int i)",
                        string("Unable to move forward to column #")
                        + to_str(info_nb0[i] + "."));
#endif
        char* buf = new char[info_nb0[i] - pos];
        FileStream.read(buf, info_nb0[i] - pos);
        delete [] buf;
        buf = new char[info_nb1[i] + 1];
        FileStream.read(buf, info_nb1[i]);
        buf[info_nb1[i]] = '\0';
        delete [] buf;
      }
    else if (info_str[i] == "e")
      {
        for (int j = 0; j < info_nb0[i]; j++)
          FileStream.GetElement();
      }
    else if (info_str[i] == "s")
      {
        char* buf = new char[info_nb0[i]];
        FileStream.read(buf, info_nb0[i]);
        delete [] buf;
      }
    else if (info_str[i] == "a" || info_str[i] == "A")
      {
        FileStream.GetFullLine();
        FileStream.seekg(-1, ExtStream::cur);
      }
  }

  //! Reads data.
  template <class T>
  int FormatFormattedText::ReadMarkup(ExtStream& FileStream,
                                      streampos pos, int i,
                                      T* value, int max_length) const
  {
    if (info_str[i] == "c")
      {
#ifdef SELDONDATA_DEBUG_CHECK_IO
        if (info_nb0[i] < pos)
          throw IOError("FormatFormattedText::ReadMarkup<T>"
                        "(ExtStream& FileStream, streampos pos, int i,"
                        " T* value, int max_length)",
                        string("Unable to move forward to column #")
                        + to_str(info_nb0[i] + "."));
#endif
        char* buf = new char[info_nb0[i] - pos];
        FileStream.read(buf, info_nb0[i] - pos);
        delete [] buf;
        buf = new char[info_nb1[i] + 1];
        FileStream.read(buf, info_nb1[i]);
        buf[info_nb1[i]] = '\0';
        string sbuf(buf);
        convert(sbuf, *value);
        delete [] buf;
        return 1;
      }
    else if (info_str[i] == "e")
      {
#ifdef SELDONDATA_DEBUG_CHECK_IO
        if (max_length < info_nb0[i])
          throw IOError("FormatFormattedText::ReadMarkup<T>"
                        "(ExtStream& FileStream, streampos pos, int i,"
                        " T* value, int max_length)",
                        string("Unable to read ") + to_str(info_nb0[i])
                        + " elements: data array is full.");
#endif
        for (int j = 0; j < info_nb0[i]; j++)
          FileStream.GetElement(value[j]);
        return info_nb0[i];
      }
    else if (info_str[i] == "s")
      {
#ifdef SELDONDATA_DEBUG_CHECK_IO
        throw IOError("FormatFormattedText::ReadMarkup<T>"
                      "(ExtStream& FileStream, streampos pos, int i,"
                      " T* value, int max_length)",
                      "Attempted to read data supposed to be skipped.");
#endif
      }
    else if (info_str[i] == "a")
      {
        convert(trim(FileStream.GetLine(), delimiters_), *value);
        FileStream.seekg(-1, ExtStream::cur);
        FileStream.GetFullLine();
        FileStream.seekg(-1, ExtStream::cur);
        pos = FileStream.tellg();
        return 1;
      }
    else if (info_str[i] == "A")
      {
        convert(FileStream.GetFullLine(), *value);
        FileStream.seekg(-1, ExtStream::cur);
        pos = FileStream.tellg();
        return 1;
      }
    return 0;
  }

  //! Returns the current format description.
  /*!
    \return The current format description.
  */
  string FormatFormattedText::GetFormat() const
  {
    return format_;
  }

  //! Returns the current delimiters.
  /*!
    \return The current delimiters.
  */
  string FormatFormattedText::GetDelimiters() const
  {
    return delimiters_;
  }

  //! Returns the current characters that denote a comment line.
  /*!
    \return The current characters that denote a comment line.
  */
  string FormatFormattedText::GetComments() const
  {
    return comments_;
  }

  //! Sets the format description.
  /*!
    \param format the new format description.
  */
  void FormatFormattedText::SetFormat(string format)
  {
    format_ = format;
    this->SetVectors();
  }

  //! Sets the delimiters.
  /*!
    \param delimiters the new delimiters.
  */
  void FormatFormattedText::SetDelimiters(string delimiters)
  {
    delimiters_ = delimiters;
  }

  //! Sets the characters that denote a comment line.
  /*!
    \param the new characters that denote a comment line.
  */
  void FormatFormattedText::SetComments(string comments)
  {
    comments_ = comments;
  }

  /********/
  /* Grid */
  /********/

  //! Reads a text file.
  template<class TG>
  void FormatFormattedText::Read(string FileName, string extract,
                                 RegularGrid<TG>& G) const
  {

    this->Read(FileName, extract, G.GetArray());

  }

  //! Reads a text file.
  template<class TG>
  void FormatFormattedText::Read(ExtStream& FileStream, string extract,
                                 RegularGrid<TG>& G) const
  {

    this->Read(FileStream, extract, G.GetArray());

  }

  //! Reads a text file.
  template<class TG, int n>
  void FormatFormattedText::Read(string FileName, string extract,
                                 GeneralGrid<TG, n>& G) const
  {

    this->Read(FileName, extract, G.GetArray());

  }

  //! Reads a text file.
  template<class TG, int n>
  void FormatFormattedText::Read(ExtStream& FileStream, string extract,
                                 GeneralGrid<TG, n>& G) const
  {

    this->Read(FileStream, extract, G.GetArray());

  }

  /********/
  /* Data */
  /********/

  //! Reads a text file.
  template<class TD, int N, class TG>
  void FormatFormattedText::Read(string FileName, string extract,
                                 Data<TD, N, TG>& D) const
  {

    this->Read(FileName, extract, D.GetArray());

  }

  //! Reads a text file.
  template<class TD, int N, class TG>
  void FormatFormattedText::Read(ExtStream& FileStream, string extract,
                                 Data<TD, N, TG>& D) const
  {

    this->Read(FileStream, extract, D.GetArray());

  }

  /*********/
  /* Array */
  /*********/

  //! Reads a text file.
  template<class TA, int N>
  void FormatFormattedText::Read(string FileName, string extract,
                                 Array<TA, N>& A) const
  {

    ExtStream FileStream(FileName, comments_, delimiters_);

#ifdef SELDONDATA_DEBUG_CHECK_IO
    // Checks if the file was opened.
    if (!FileStream.is_open())
      throw IOError("FormatFormattedText::Read"
                    "(string FileName, string extract, Array<TA, N>& A)",
                    "Unable to open file \"" + FileName + "\".");
#endif

    this->Read(FileStream, extract, A);

    FileStream.close();

  }

  //! Reads a text file.
  template<class TA, int N>
  void FormatFormattedText::Read(ExtStream& FileStream, string extract,
                                 Array<TA, N>& A) const
  {

#ifdef SELDONDATA_DEBUG_CHECK_IO
    // Checks if the file ready.
    if (!FileStream.good())
      throw IOError("FormatFormattedText::Read(ExtStream& FileStream,"
                    " string extract, Array<TA, N>& A)",
                    "File is not ready.");
#endif

    string delimiters(FileStream.GetDelimiters()),
      comments(FileStream.GetComments());

    FileStream.SetDelimiters(delimiters_);
    FileStream.SetComments(comments_);

    FileStream.Skip();

    streampos pos_beg = FileStream.tellg();
    streampos pos_cur(pos_beg);
    streampos position;

    vector<int> markups;
    split(extract, markups, " \t\n,;:/|");

    int nb_elements = A.numElements();
    TA* data = A.data();
    int i = 0;
    int k(0), l(0);
    bool newline;

    while ((i < nb_elements) && !is_emptystream(FileStream))
      {

        if (k < int(markups.size()))
          {
            while (l != markups[k])
              {
                this->SkipMarkup(FileStream, pos_cur - pos_beg, l);
                newline = false;
                while ((!is_emptystream(FileStream))
                       && FileStream.Discard(FileStream.
                                             PeekFullLine(position)))
                  {
                    newline = true;
                    FileStream.seekg(position);
                  }
                if (newline)
                  {
                    FileStream.seekg(-1, ExtStream::cur);
                    pos_beg = FileStream.tellg();
                  }
                pos_cur = FileStream.tellg();
                l++;
              }

            i += this->ReadMarkup(FileStream, pos_cur - pos_beg, l,
                                  &data[i], nb_elements - i);

#ifdef SELDONDATA_DEBUG_CHECK_IO
            // Checks if data was read.
            if (!FileStream.good())
              throw IOError("FormatFormattedText::Read(ExtStream& FileStream,"
                            "string extract, Array<TA, N>& A)",
                            string("Unable to read data associated with <")
                            + info_str[l] + " " + to_str(info_nb0[l]) + " "
                            + to_str(info_nb1[l]) + ">.");
#endif

            newline = false;
            while (!is_emptystream(FileStream)
                   && FileStream.Discard(FileStream.PeekFullLine(position)))
              {
                newline = true;
                FileStream.seekg(position);
              }
            if (newline)
              {
                FileStream.seekg(-1, ExtStream::cur);
                pos_beg = FileStream.tellg();
              }
            pos_cur = FileStream.tellg();

            k++;
            if (!is_emptystream(FileStream))
              l++;
          }
        else
          {
            k = 0;
            l = 0;
            FileStream.GetFullLine();
            while ((!is_emptystream(FileStream))
                   && FileStream.Discard(FileStream.PeekFullLine(position)))
              FileStream.seekg(position);
            pos_beg = pos_cur = FileStream.tellg();
          }
      }

#ifdef SELDONDATA_DEBUG_CHECK_IO
    // Checks if data was read.
    if (i < nb_elements)
      throw IOError("FormatFormattedText::Read"
                    "(ExtStream& FileStream, string extract,"
                    " Array<TA, N>& A)",
                    string("Only ") + to_str(i) +
                    " elements were read instead of "
                    + to_str(nb_elements) + ".");
#endif

    FileStream.SetDelimiters(delimiters);
    FileStream.SetComments(comments);

  }


  //////////////////
  // FORMATNETCDF //
  //////////////////

#ifdef SELDONDATA_WITH_NETCDF

  //! Default constructor.
  template<class T>
  FormatNetCDF<T>::FormatNetCDF()  throw()
  {
  }

  //! Destructor.
  template<class T>
  FormatNetCDF<T>::~FormatNetCDF()  throw()
  {
  }

  ////////////////////////////////////////////////
  // Reads one variable in a netcdf format file //
  ////////////////////////////////////////////////

  /********/
  /* Grid */
  /********/

  //! Reads a netCDF file.
  template<class T>
  template<class TG>
  void FormatNetCDF<T>::Read(string FileName, string variable,
                             RegularGrid<TG>& G) const
  {

    this->Read(FileName, variable, G.GetArray());

  }

  //! Reads a netCDF file.
  template<class T>
  template<class TG, int n>
  void FormatNetCDF<T>::Read(string FileName, string variable,
                             GeneralGrid<TG, n>& G) const
  {

    this->Read(FileName, variable, G.GetArray());

  }

  /********/
  /* Data */
  /********/

  //! Reads a netCDF file.
  template<class T>
  template<class TD, int N, class TG>
  void FormatNetCDF<T>::Read(string FileName, string variable,
                             Data<TD, N, TG>& D) const
  {

    this->Read(FileName, variable, D.GetArray());

  }

  /*********/
  /* Array */
  /*********/

  //! Reads a netCDF file.
  template<class T>
  template<class TA, int N>
  void FormatNetCDF<T>::Read(string FileName, string variable,
                             Array<TA, N>& A) const
  {

    NcFile File(FileName.c_str());

#ifdef SELDONDATA_DEBUG_CHECK_IO
    // Checks if the file is valid.
    if (!File.is_valid())
      throw IOError("FormatNetCDF<T>::Read(string FileName, Array<TA, N>& A)",
                    "\"" + FileName + "\" is not a valid netCDF file.");
#endif

    int Nb_vars = File.num_vars();

    int i(0);
    while ((i < Nb_vars) && (string(File.get_var(i)->name()) != variable))
      i++;

#ifdef SELDONDATA_DEBUG_CHECK_IO
    // Checks whether the variable was found.
    if (i == Nb_vars)
      throw IOError("FormatNetCDF<T>::Read(string FileName, Array<TA, N>& A)",
                    "Unable to find variable \"" + variable
                    + "\" in \"" + FileName + "\".");
#endif

    NcVar* var = File.get_var(i);

#ifdef SELDONDATA_DEBUG_CHECK_DIMENSIONS
    // Checks the dimension.
    if (var->num_dims() != N)
      throw WrongDim("FormatNetCDF<T>::"
                     "Read(string FileName, Array<TA, N>& A)",
                     "Data has " + to_str(N) +
                     "dimensions, but stored data has "
                     + to_str(var->num_dims()) + "dimensions.");
#endif

#ifdef SELDONDATA_DEBUG_CHECK_INDICES
    long* input_dimensions = var->edges();
    for (i = 0; i < var->num_dims(); i++)
      if (A.extent(i) > input_dimensions[i])
        throw WrongIndex("FormatNetCDF<T>::"
                         "Read(string FileName, Array<TA, N>& A)",
                         "Array extent is " + to_str(A.extent(i))
                         + " along dimension #" + to_str(i)
                         + " , but it should not be strictly more than "
                         + to_str(input_dimensions[i]) + ".");
    delete[] input_dimensions;
#endif

    long* extents = new long[N];
    for (i = 0; i < N; i++)
      extents[i] = A.extent(i);

    bool op = var->get(A.data(), extents);

    delete [] extents;

#ifdef SELDONDATA_DEBUG_CHECK_IO
    // Checks whether input operation succeeded.
    if (!op)
      throw IOError("FormatNetCDF<T>::Read(string FileName, Array<TA, N>& A)",
                    "Data type doesn't match type of stored values.");
#endif

  }

  /////////////////////////////////////////////////
  // Reads one dimension in a netcdf format file //
  /////////////////////////////////////////////////


  //! Reads the dim_num th dimension of the variable in a netCDF file.
  template<class T>
  void FormatNetCDF<T>::ReadDimension(string FileName, string variable,
                                      int dim_num, int& dim_value) const
  {

    NcFile File(FileName.c_str());

#ifdef SELDONDATA_DEBUG_CHECK_IO
    // Checks if the file is valid.
    if (!File.is_valid())
      throw IOError("FormatNetCDF<T>::ReadDimension(string, string, int, int&)",
                    "\"" + FileName + "\" is not a valid netCDF file.");
#endif

    int Nb_vars = File.num_vars();

    int i(0);
    while ((i < Nb_vars) && (string(File.get_var(i)->name()) != variable))
      i++;

#ifdef SELDONDATA_DEBUG_CHECK_IO
    // Checks whether the variable was found.
    if (i == Nb_vars)
      throw IOError("FormatNetCDF<T>::ReadDimension(string, string, int, int&)",
                    "Unable to find variable \"" + variable
                    + "\" in \"" + FileName + "\".");
#endif

    NcVar* var = File.get_var(i);

#ifdef SELDONDATA_DEBUG_CHECK_DIMENSIONS
    // Checks the dimension.
    if (var->num_dims() < dim_num)
      throw WrongDim("FormatNetCDF<T>::Read(string, string, int, int&)",
                     "Variable has " + to_str(var->num_dims()) + " dimensions,"
                     + "but requested dimension is " + to_str(dim_num));
#endif

    long* shape = var->edges();
    dim_value = shape[dim_num];
    delete[] shape;
  }

  /////////////////////////////////////////////////
  // Reads an attribute in a netcdf format file //
  /////////////////////////////////////////////////


  //! Reads the global attribute in a netCDF file
  //! and put the value of it in a float format.
  template<class T>
  void FormatNetCDF<T>::ReadAttribute(string FileName, string attribute,
                                      float& value) const
  {

    NcFile File(FileName.c_str());

#ifdef SELDONDATA_DEBUG_CHECK_IO
    // Checks if the file is valid.
    if (!File.is_valid())
      throw IOError("FormatNetCDF<T>::ReadAttribute(string, string, float&)",
                    "\"" + FileName + "\" is not a valid netCDF file.");
#endif

    NcAtt* att = 0;
    int Nb_atts = File.num_atts();
    for (int i = 0; i < Nb_atts; ++i)
      {
        att = File.get_att(i);
        if (attribute == att->name())
          break;
        delete att;
        att = 0;
      }

#ifdef SELDONDATA_DEBUG_CHECK_IO
    // Checks whether the attribute was found.
    if (att == 0)
      throw IOError("FormatNetCDF<T>::ReadAttribute(string, string, float&)",
                    "Unable to find global attribute \"" + attribute
                    + "\" in \"" + FileName + "\".");
#endif

#ifdef SELDONDATA_DEBUG_CHECK_DIMENSIONS
    // Checks if the attribute value is a float.
    if (att->type() != ncFloat)
      throw IOError("FormatNetCDF<T>::ReadAttribute(string, string, float&)",
                    "Attribute \"" + attribute + "\" is not a float"
                    + " in \"" + FileName + "\".");
#endif
    value = att->as_float(0);
    delete att;
  }

  //! Reads the global attribute in a netCDF file
  //! and put the value of it in a int format.
  template<class T>
  void FormatNetCDF<T>::ReadAttribute(string FileName, string attribute,
                                      int& value) const
  {

    NcFile File(FileName.c_str());

#ifdef SELDONDATA_DEBUG_CHECK_IO
    // Checks if the file is valid.
    if (!File.is_valid())
      throw IOError("FormatNetCDF<T>::ReadAttribute(string, string, int&)",
                    "\"" + FileName + "\" is not a valid netCDF file.");
#endif

    NcAtt* att = 0;
    int Nb_atts = File.num_atts();
    for (int i = 0; i < Nb_atts; ++i)
      {
        att = File.get_att(i);
        if (attribute == att->name())
          break;
        delete att;
        att = 0;
      }

#ifdef SELDONDATA_DEBUG_CHECK_IO
    // Checks whether the attribute was found.
    if (att == 0)
      throw IOError("FormatNetCDF<T>::ReadAttribute(string, string, int&)",
                    "Unable to find global attribute \"" + attribute
                    + "\" in \"" + FileName + "\".");
#endif

#ifdef SELDONDATA_DEBUG_CHECK_DIMENSIONS
    // Checks if the attribute value is a float.
    if (att->type() != ncInt)
      throw IOError("FormatNetCDF<T>::ReadAttribute(string, string, int&)",
                    "Attribute \"" + attribute + "\" is not a float"
                    + " in \"" + FileName + "\".");
#endif
    value = att->as_int(0);
    delete att;
  }

#endif


  ////////////////
  // FORMATGRIB //
  ////////////////

#ifdef SELDONDATA_WITH_GRIB

  //! Default constructor.
  FormatGrib::FormatGrib()  throw()
  {
  }

  //! Destructor.
  FormatGrib::~FormatGrib()  throw()
  {
  }

  /********/
  /* Grid */
  /********/

  //! Reads a Grib file.
  template<class TG>
  void FormatGrib::Read(string FileName, int variable,
                        RegularGrid<TG>& G) const
  {

    this->Read(FileName, variable, G.GetArray());

  }

  //! Reads a Grib file.
  template<class TG, int n>
  void FormatGrib::Read(string FileName, int variable,
                        GeneralGrid<TG, n>& G) const
  {

    this->Read(FileName, variable, G.GetArray());

  }

  /********/
  /* Data */
  /********/

  //! Reads a Grib file.
  template<class TD, int N, class TG>
  void FormatGrib::Read(string FileName, int variable,
                        Data<TD, N, TG>& D) const
  {

    this->Read(FileName, variable, D.GetArray());

  }

  /*********/
  /* Array */
  /*********/


  //! Reads a Grib file.
  template <class TA, int N>
  void FormatGrib::Read(string FileName, int variable, Array<TA, N>& A) const
  {
    int i, j, nx, ny;
    long int pos = 0;
    FILE *grib_file;

    int status(0), param(0);

    grib_file = fopen(FileName.c_str(), "r");

#ifdef SELDONDATA_DEBUG_CHECK_IO
    // Checks if the file was opened.
    if (grib_file == NULL)
      throw IOError("FormatGrib::Read"
                    "(string FileName, int variable, Array<TA, N>& A)",
                    "Unable to open file \"" + FileName + "\".");
#endif

    int nb_elements = A.numElements();
    int max_length(nb_elements);
    size_t nrec(0);

    float *data;

    while (max_length != 0 &&
           (status = decode_grib(grib_file, &pos, &param, &data, &nx, &ny))
           == 0)
      {
        j = nb_elements - max_length;
        if (param == variable)
          {
            for (i = 0; i < nx * ny; i++)
              A.data()[j + i] = data[i];
            max_length -= nx * ny;
          }
        nrec ++;
        free(data);
      }

#ifdef SELDONDATA_DEBUG_CHECK_IO
    if (status == 1)
      throw IOError("FormatGrib::Read(string FileName, Array<TA, N>& A)",
                    "Read error after " + to_str(nrec) + " records.");

    if (max_length != 0)
      {
        if (status == -2)
          throw IOError("FormatGrib::Read(string FileName, Array<TA, N>& A)",
                        "File \"" + FileName + "\" contains at least "
                        + to_str(nb_elements - max_length + nx * ny)
                        + " elements for field #" + to_str(variable)
                        + ", but data has only " + to_str(nb_elements)
                        + " elements. The current record (whose length is "
                        + to_str(nx * ny)
                        + " elements) must be completely read.");

        if (status == -1)
          throw IOError("FormatGrib::Read(string FileName, Array<TA, N>& A)",
                        "End of file found. "
                        + to_str(nb_elements - max_length)
                        + " elements were read for field #"
                        + to_str(variable)
                        + " in file \"" + FileName + "\", but data has "
                        + to_str(nb_elements) + " elements.");

        throw IOError("FormatGrib::Read(string FileName, Array<TA, N>& A)",
                      "Cannot find all values. File \"" + FileName
                      + "\" contains " + to_str(nb_elements - max_length)
                      + " elements for field #" + to_str(variable)
                      + ", but data has " + to_str(nb_elements)
                      + " elements.");
      }
#endif

    fclose(grib_file);
  }

  // Format Grib2.
  //! Default constructor.
  FormatGrib2::FormatGrib2()  throw()
  {
  }

  //! Destructor.
  FormatGrib2::~FormatGrib2()  throw()
  {
  }

  /********/
  /* Grid */
  /********/

  //! Reads a Grib2 file.
  template<class TG>
  void FormatGrib2::Read(string FileName, int discipline, int parameterCategory,
			 int parameterNumber, RegularGrid<TG>& G) const
  {

    this->Read(FileName, discipline, parameterCategory, parameterNumber,
	       G.GetArray());

  }

  //! Reads a Grib2 file.
  template<class TG, int n>
  void FormatGrib2::Read(string FileName, int discipline, int parameterCategory,
			 int parameterNumber, GeneralGrid<TG, n>& G) const
  {

    this->Read(FileName, discipline, parameterCategory, parameterNumber,
	       G.GetArray());

  }

  /********/
  /* Data */
  /********/
  
  //! Reads a Grib2 file.
  template<class TD, int N, class TG>
  void FormatGrib2::Read(string FileName, int discipline, int parameterCategory,
			 int parameterNumber, Data<TD, N, TG>& D) const
  {

    this->Read(FileName, discipline, parameterCategory, parameterNumber,
	       D.GetArray());

  }

  /*********/
  /* Array */
  /*********/

  //! Reads a Grib file.
  template <class TA, int N>
  void FormatGrib2::Read(string FileName, int discipline, int parameterCategory,
			 int parameterNumber, Array<TA, N>& A) const
  {
    FILE* grib_file;
 
    grib_file = fopen(FileName.c_str(), "r");
    grib_multi_support_on(0);

#ifdef SELDONDATA_DEBUG_CHECK_IO
    // Checks if the file was opened.
    if (grib_file == NULL)
      throw IOError("FormatGrib2::Read(string FileName, int discipline, int parameterCategory, int parameterNumber, Array<TA, N>& A)",
		    "Unable to open file \"" + FileName + "\".");
#endif

    int nb_elements = A.numElements();
    int max_length(nb_elements);

    int Nt = 1, Nz = 1;
    if (A.dimensions() == 4)
      {
	Nt = A.shape()[0];
	Nz = A.shape()[1];
      }
    else if (A.dimensions() == 3)
      Nz = A.shape()[0];

    RegularGrid<double> GridX(nb_elements / Nt / Nz);
    Data<double, 1> data(GridX);
    int data_size;

    int err = 0;
    long disc, paramCat, paramNum;
    grib_handle *h = NULL;
    size_t values_len = nb_elements / Nt / Nz;

    while ((h = grib_handle_new_from_file(0, grib_file, &err)) != NULL)
      {
	GRIB_CHECK(err, 0);
	GRIB_CHECK(grib_get_long(h, "discipline", &disc), 0);
	GRIB_CHECK(grib_get_long(h, "parameterCategory", &paramCat), 0);
	GRIB_CHECK(grib_get_long(h, "parameterNumber", &paramNum), 0);

	if (disc == discipline && paramCat == parameterCategory \
	    && paramNum == parameterNumber)
	  {
	    GRIB_CHECK(grib_get_double_array(h, "values", data.GetData(),
					     &values_len), 0);

	    int j = nb_elements - max_length;
	    data_size = data.GetNbElements();
	    for(int i = 0; i < data_size; i++)
	      A.data()[j + i] = data(i);
	    max_length -= data_size;
	  }
      }
#ifdef SELDONDATA_DEBUG_CHECK_IO
    if (max_length < 0)
      throw IOError("FormatGrib2::Read(string FileName, Array<TA, N>& A)",
		    "Cannot find all values. File \"" + FileName
		    + "\" contains " + to_str(nb_elements - max_length)
		    + " elements for field #" + to_str(discipline)
		    + ", " + to_str(parameterCategory) + ", "
		    + to_str(parameterCategory)
		    + ", but data has " + to_str(nb_elements)
		    + " elements.");

#endif
    grib_handle_delete(h);
    data.Resize();
    fclose(grib_file);
  }

#endif


  ///////////////////
  // FORMATCHIMERE //
  ///////////////////

  //! Default constructor.
  FormatChimere::FormatChimere()  throw()
  {
    date_ = -1;
  }

  //! Constructor.
  FormatChimere::FormatChimere(int date)  throw()
  {
    date_ = date;
  }

  //! Destructor.
  FormatChimere::~FormatChimere()  throw()
  {
  }

  //! Sets the date.
  /*!
    \param date date.
  */
  void FormatChimere::SetDate(int date)
  {
    date_ = date;
  }

  //! Get the date.
  /*!
    \return The date.
  */
  int FormatChimere::GetDate() const
  {
    return date_;
  }

  /********/
  /* Data */
  /********/

  //! Reads a file in "Chimere" format.
  template<class TD, int N, class TG>
  void FormatChimere::Read(string FileName, Data<TD, N, TG>& D,
                           int nb_lines) const
  {

    this->Read(FileName, D.GetArray(), nb_lines);

  }

  //! Reads a file in "Chimere" format.
  template<class TD, int N, class TG>
  void FormatChimere::Read(ifstream& FileStream, Data<TD, N, TG>& D,
                           int nb_lines) const
  {

    this->Read(FileStream, D.GetArray(), nb_lines);

  }

  /*********/
  /* Array */
  /*********/

  //! Reads a file in "Chimere" format.
  template<class TA, int N>
  void FormatChimere::Read(string FileName, Array<TA, N>& A,
                           int nb_lines) const
  {

    ifstream FileStream;
    FileStream.open(FileName.c_str(), ifstream::in);
    FileStream.flags(fstream::scientific | fstream::skipws);

#ifdef SELDONDATA_DEBUG_CHECK_IO
    // Checks if the file was opened.
    if (!FileStream.is_open())
      throw IOError("FormatChimere::Read"
                    "(string FileName, Array<TA, N>& A, int nb_lines)",
                    "Unable to open file \"" + FileName + "\".");
#endif

    this->Read(FileStream, A);

    FileStream.close();

  }

  //! Reads a file in "Chimere" format.
  template<class TA, int N>
  void FormatChimere::Read(ifstream& FileStream, Array<TA, N>& A,
                           int nb_lines) const
  {

#ifdef SELDONDATA_DEBUG_CHECK_IO
    // Checks if the file ready.
    if (!FileStream.good())
      throw IOError("FormatChimere::Read"
                    "(ifstream& FileStream, Array<TA, N>& A, int nb_lines)",
                    "File is not ready.");
#endif

    int nb_elements = A.numElements();
    int i, j, k;
    char c;
    bool reading = false;
    TA temp;
    int date, level;
    Array<int, 1> Index(10), Length(10);

    for (j = 0; j < 10; j++)
      {
        Index(j) = 0;
        Length(j) = A.extent(j);
      }

    int nb_elements_per_line = Length(N - 1);

    if (nb_lines == -1)
      if (N == 2)
        nb_lines = 1;
      else if (N == 3)
        nb_lines = Length(1);
      else
        nb_lines = nb_elements
          / nb_elements_per_line
          / Length(0)
          / Length(1);

    fstream::fmtflags flags = FileStream.flags();

    c = FileStream.peek();
    while ((FileStream.good())
           && ((c < '0') || (c > '9'))
           && (c != '.') && (c != '-')
           && (c != '+'))
      {
        FileStream.ignore(1);
        c = FileStream.peek();
      }

    FileStream.flags(ifstream::dec | fstream::skipws);
    FileStream >> date;
    FileStream.flags(flags);

    reading = ((date_ == -1) || (date == date_));

    string temp_str;
    while ((!reading) && (FileStream.good()))
      {
        for (i = 0; i < nb_lines + 1; i++)
          getline(FileStream, temp_str);

        FileStream.flags(ifstream::dec | fstream::skipws);
        FileStream >> date;
        FileStream.flags(flags);

        reading = ((date_ == -1) || (date == date_));
      }

    i = 0;
    while ((i < nb_elements) && (FileStream.good()))
      {

        // Reads level.
        FileStream.flags(ifstream::dec | fstream::skipws);
        FileStream >> level;
        FileStream.flags(flags);

        k = 0;
        j = N - 1;

        while ((k < nb_lines * nb_elements_per_line)
               && (FileStream.good()))
          {

            c = FileStream.peek();
            while ((FileStream.good())
                   && ((c < '0') || (c > '9'))
                   && (c != '.') && (c != '-')
                   && (c != '+'))
              {
                FileStream.ignore(1);
                c = FileStream.peek();
              }

            // Reads level.
            if (reading)
              FileStream >> A(Index(0), Index(1), Index(2),
                              Index(3), Index(4), Index(5),
                              Index(6), Index(7), Index(8),
                              Index(9));
            else
              FileStream >> temp;

            if (reading)
              {
                j = N - 1;
                while ((j >= 0) && (Index(j) == Length(j) - 1))
                  {
                    Index(j) = 0;
                    j--;
                  }

                if (j != -1)
                  Index(j)++;
              }

            k++;
          }

        if (reading)
          i += k;

        c = FileStream.peek();
        while ((FileStream.good())
               && ((c < '0') || (c > '9'))
               && (c != '.') && (c != '-')
               && (c != '+'))
          {
            FileStream.ignore(1);
            c = FileStream.peek();
          }
        FileStream.flags(ifstream::dec | fstream::skipws);
        FileStream >> date;
        FileStream.flags(flags);

        reading = (reading || (date_ == -1) || (date == date_));

      }

#ifdef SELDONDATA_DEBUG_CHECK_IO
    // Checks if all was read.
    if (!reading)
      throw IOError("FormatChimere::"
                    "Read(ifstream& FileStream, Array<TA, N>& A)",
                    "The date was not found.");
    // Checks if all was read.
    if (i != nb_elements)
      throw IOError("FormatChimere::"
                    "Read(ifstream& FileStream, Array<TA, N>& A)",
                    to_str(i) + " elements were read instead of "
                    + to_str(nb_elements) + ".");
#endif

  }


}  // namespace SeldonData.

#define FILE_SELDONDATA_FORMAT_CXX
#endif
