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

#ifndef FILE_SELDONDATA_FUNCTIONS_CXX

#include "Functions.hxx"
#include "stdio.h"

namespace SeldonData
{

  /////////////
  // REGULAR //

  //! Linear interpolation for data defined on regular grids.
  /*!
    Linear interpolation from data defined on regular grids to
    data defined on regular grids.
    \param dataIn reference data.
    \param dataOut interpolated data (on exit).
  */
  template < int N, class TIn, class TGIn,
             class TOut, class TGOut >
  void LinearInterpolationRegular(Data<TIn, N, TGIn>& dataIn,
                                  Data<TOut, N, TGOut>& dataOut)
  {

    int i, j, k, l, m;
    Array<int, 1> IndexIn(10), IndexIn0(10), IndexOut(10);
    Array<int, 1> LengthIn(10), LengthOut(10);
    Array<TIn, 1> Coeff(10), Coeff0(10);
    Array<bool, 1> Pos(10);
    TIn coeff;

    for (i = 0; i < 10; i++)
      {
        Pos(i) = 0;
        LengthIn(i) = dataIn.GetLength(i);
        LengthOut(i) = dataOut.GetLength(i);
        IndexOut(i) = 0;
        IndexIn(i) = 0;
        while ((IndexIn(i) < LengthIn(i))
               && (dataIn[i](IndexIn(i)) < dataOut[i](0)))
          IndexIn(i)++;
        if (IndexIn(i) == LengthIn(i))
          IndexIn(i) = LengthIn(i) - 1;
        else if (IndexIn(i) == 0)
          IndexIn(i) = 1;
        IndexIn0(i) = IndexIn(i);
        if (LengthIn(i) != 0)
          Coeff0(i) = (dataOut[i](0) - dataIn[i](IndexIn(i) - 1)) /
            (dataIn[i](IndexIn(i)) - dataIn[i](IndexIn(i) - 1));
        else
          Coeff0(i) = TIn(0);
        Coeff(i) = Coeff0(i);
      }

    j = N - 1;
    for (i = 0; i < dataOut.GetNbElements(); i++)
      {

        while ((IndexIn(j) < LengthIn(j))
               && (dataIn[j](IndexIn(j)) < dataOut[j](IndexOut(j))))
          IndexIn(j)++;

        if (IndexIn(j) == LengthIn(j))
          IndexIn(j) = LengthIn(j) - 1;
        else if (IndexIn(j) == 0)
          IndexIn(j) = 1;

        Coeff(j) = (dataOut[j](IndexOut(j)) - dataIn[j](IndexIn(j) - 1)) /
          (dataIn[j](IndexIn(j)) - dataIn[j](IndexIn(j) - 1));

        dataOut.Value(IndexOut(0), IndexOut(1), IndexOut(2),
                      IndexOut(3), IndexOut(4), IndexOut(5),
                      IndexOut(6), IndexOut(7), IndexOut(8),
                      IndexOut(9)) = TOut(0);


        for (k = 0; k < int(pow(2.0, double(N))); k++)
          {
            l = k;
            coeff = TIn(1);
            for (m = 0; m < N; m++)
              {
                Pos(m) = l % 2;
                if (l % 2 == 1)
                  coeff *= TIn(1) - Coeff(m);
                else
                  coeff *= Coeff(m);
                l = l / 2;
              }

            dataOut.Value(IndexOut(0), IndexOut(1), IndexOut(2),
                          IndexOut(3), IndexOut(4), IndexOut(5),
                          IndexOut(6), IndexOut(7), IndexOut(8),
                          IndexOut(9)) +=
              TOut(coeff *
                   dataIn.Value(IndexIn(0) - Pos(0),
                                IndexIn(1) - Pos(1),
                                IndexIn(2) - Pos(2),
                                IndexIn(3) - Pos(3),
                                IndexIn(4) - Pos(4),
                                IndexIn(5) - Pos(5),
                                IndexIn(6) - Pos(6),
                                IndexIn(7) - Pos(7),
                                IndexIn(8) - Pos(8),
                                IndexIn(9) - Pos(9)));

          }

        j = N - 1;
        while ((j >= 0) && (IndexOut(j) == LengthOut(j) - 1))
          {
            IndexOut(j) = 0;
            IndexIn(j) = IndexIn0(j);
            Coeff(j) = Coeff0(j);
            j--;
          }
        if (j != -1)
          IndexOut(j)++;

      }

  }

  // REGULAR //
  /////////////


  ////////////////////////
  // UNIFORM TO GENERAL //

  //! Linear interpolation from data defined on uniform grids
  //! to data defined on any grid.
  /*! Linear interpolation from an input data defined on uniform grids only
    to an output data defined on any type of grids. Uniform grids are
    regular grids with fixed steps.
    \param dataIn reference data defined on uniform grids.
    \param dataOut interpolated data (on exit) defined on any type of grids.
  */
  template < int N, class TIn, class TGIn,
             class TOut, class TGOut >
  void LinearInterpolationUniformToGeneral(Data<TIn, N, TGIn>& dataIn,
                                           Data<TOut, N, TGOut>& dataOut)
  {

    int i, j, k, l, m;
    Array<int, 1> IndexIn(10), IndexOut(10);
    Array<int, 1> LengthIn(10), LengthOut(10);
    Array<TIn, 1> MinIn(10), DeltaIn(10), Coeff(10);
    Array<bool, 1> Pos(10);
    TIn coord_out, coord_in, coeff;

    for (i = 0; i < 10; i++)
      {
        Pos(i) = 0;
        LengthIn(i) = dataIn.GetLength(i);
        LengthOut(i) = dataOut.GetLength(i);
        Coeff(i) = 0;
        IndexOut(i) = 0;
        IndexIn(i) = 0;
      }

    for (i = 0; i < N; i++)
      {
        MinIn(i) = dataIn[i](0);
        DeltaIn(i) = dataIn[i](1) - dataIn[i](0);
      }

    j = N - 1;
    for (i = 0; i < dataOut.GetNbElements(); i++)
      {

        for (k = 0; k < N; k++)
          {
            coord_out = dataOut[k].Value(IndexOut(0), IndexOut(1), IndexOut(2),
                                         IndexOut(3), IndexOut(4), IndexOut(5),
                                         IndexOut(6), IndexOut(7), IndexOut(8),
                                         IndexOut(9));
            l = int((coord_out - MinIn(k)) / DeltaIn(k));
            l = l < 1 ? 1 : l;
            l = l < LengthIn(k) ? l : (LengthIn(k) - 1);
            IndexIn(k) = l;

            coord_in = MinIn(k) + TIn(l - 1) * DeltaIn(k);
            Coeff(k) = (coord_out - coord_in) / DeltaIn(k);
          }

        dataOut.Value(IndexOut(0), IndexOut(1), IndexOut(2),
                      IndexOut(3), IndexOut(4), IndexOut(5),
                      IndexOut(6), IndexOut(7), IndexOut(8),
                      IndexOut(9)) = TOut(0);

        for (k = 0; k < int(pow(2.0, double(N))); k++)
          {
            l = k;
            coeff = TIn(1);
            for (m = 0; m < N; m++)
              {
                Pos(m) = l % 2;
                if (l % 2 == 1)
                  coeff *= TIn(1) - Coeff(m);
                else
                  coeff *= Coeff(m);
                l = l / 2;
              }

            dataOut.Value(IndexOut(0), IndexOut(1), IndexOut(2),
                          IndexOut(3), IndexOut(4), IndexOut(5),
                          IndexOut(6), IndexOut(7), IndexOut(8),
                          IndexOut(9)) +=
              TOut(coeff *
                   dataIn.Value(IndexIn(0) - Pos(0),
                                IndexIn(1) - Pos(1),
                                IndexIn(2) - Pos(2),
                                IndexIn(3) - Pos(3),
                                IndexIn(4) - Pos(4),
                                IndexIn(5) - Pos(5),
                                IndexIn(6) - Pos(6),
                                IndexIn(7) - Pos(7),
                                IndexIn(8) - Pos(8),
                                IndexIn(9) - Pos(9)));

          }

        j = N - 1;
        while ((j >= 0) && (IndexOut(j) == LengthOut(j) - 1))
          {
            IndexOut(j) = 0;
            j--;
          }
        if (j != -1)
          IndexOut(j)++;

      }

  }

  // UNIFORM TO GENERAL //
  ////////////////////////


  ////////////////////////
  // REGULAR TO GENERAL //

  //! Linear interpolation for data defined on regular grids to data
  //! defined on a general grids.
  /*!
    \param dataIn reference data defined on regular grids.
    \param dataOut interpolated data (on exit) defined on any kind of grids.
  */
  template < int N, class TIn, class TGIn,
             class TOut, class TGOut >
  void LinearInterpolationRegularToGeneral(Data<TIn, N, TGIn>& dataIn,
                                           Data<TOut, N, TGOut>& dataOut)
  {

    int i, j, k, l, m;
    Array<int, 1> IndexIn(10), IndexOut(10);
    Array<int, 1> LengthIn(10), LengthOut(10);
    Array<TIn, 1> Coeff(10);
    Array<bool, 1> Pos(10);
    TIn coeff;

    dataOut.SetZero();

    for (i = 0; i < 10; i++)
      {
        Pos(i) = 0;
        LengthIn(i) = dataIn.GetLength(i);
        LengthOut(i) = dataOut.GetLength(i);
        IndexOut(i) = 0;
      }

    for (i = 0; i < dataOut.GetNbElements(); i++)
      {

        for (k = 0; k < N; k++)
          {
            IndexIn(k) = 1;
            while ((IndexIn(k) < LengthIn(k) - 1)
                   && (dataIn[k](IndexIn(k))
                       < dataOut[k].Value(IndexOut(0), IndexOut(1),
                                          IndexOut(2), IndexOut(3),
                                          IndexOut(4), IndexOut(5),
                                          IndexOut(6), IndexOut(7),
                                          IndexOut(8), IndexOut(9))))
              IndexIn(k)++;

            Coeff(k) = (dataOut[k].Value(IndexOut(0), IndexOut(1),
                                         IndexOut(2), IndexOut(3),
                                         IndexOut(4), IndexOut(5),
                                         IndexOut(6), IndexOut(7),
                                         IndexOut(8), IndexOut(9))
                        - dataIn[k](IndexIn(k) - 1)) /
              (dataIn[k](IndexIn(k)) - dataIn[k](IndexIn(k) - 1));
          }

        for (k = 0; k < int(pow(2.0, double(N)) + 0.5); k++)
          {

            l = k;
            coeff = TIn(1);
            for (m = 0; m < N; m++)
              {
                Pos(m) = l % 2;
                if (l % 2 == 1)
                  coeff *= TIn(1) - Coeff(m);
                else
                  coeff *= Coeff(m);
                l = l / 2;
              }

            dataOut.Value(IndexOut(0), IndexOut(1), IndexOut(2),
                          IndexOut(3), IndexOut(4), IndexOut(5),
                          IndexOut(6), IndexOut(7), IndexOut(8),
                          IndexOut(9)) +=
              TOut(coeff *
                   dataIn.Value(IndexIn(0) - Pos(0),
                                IndexIn(1) - Pos(1),
                                IndexIn(2) - Pos(2),
                                IndexIn(3) - Pos(3),
                                IndexIn(4) - Pos(4),
                                IndexIn(5) - Pos(5),
                                IndexIn(6) - Pos(6),
                                IndexIn(7) - Pos(7),
                                IndexIn(8) - Pos(8),
                                IndexIn(9) - Pos(9)));

          }

        j = N - 1;
        while ((j >= 0) && (IndexOut(j) == LengthOut(j) - 1))
          {
            IndexOut(j) = 0;
            j--;
          }
        if (j != -1)
          IndexOut(j)++;

      }

  }

  // REGULAR TO GENERAL //
  ////////////////////////


  ////////////////
  // ONEGENERAL //

  //! Linear interpolation for data defined on regular grids,
  //! but one grid.
  /*!
    Performs linear interpolation on data defined on regular grids,
    except one grid which may be a general grid (i.e. depending on other coordinates).
    Both input and output data may be defined on a general grid along
    dimension 'dim', but only along this dimension. Moreover, input data or
    output data can still be defined on regular grids along dimension 'dim'.
    \param dataIn reference data.
    \param dataOut interpolated data (on exit).
    \param dim dimension related to the general grid.
  */
  template < int N, class TIn, class TGIn,
             class TOut, class TGOut >
  void LinearInterpolationOneGeneral(Data<TIn, N, TGIn>& dataIn,
                                     Data<TOut, N, TGOut>& dataOut,
                                     int dim)
  {

    int i, j, k, l, m;
    Array<int, 1> IndexIn(10), IndexIn0(10), IndexOut(10);
    Array<int, 1> LengthIn(10), LengthOut(10);
    Array<TIn, 1> Coeff(10), Coeff0(10);
    TIn temp;
    Array<bool, 1> Pos(10), Pos_dim(10);
    TIn coeff;

    dataOut.SetZero();

    for (i = 0; i < 10; i++)
      {
        Pos(i) = 0;
        Pos_dim(i) = 0;
      }

    for (i = 0; i < 10; i++)
      if ((i < N) && (i != dim))
        {
          LengthIn(i) = dataIn.GetLength(i);
          LengthOut(i) = dataOut.GetLength(i);
          IndexOut(i) = 0;
          IndexIn(i) = 1;
          while ((IndexIn(i) < LengthIn(i) - 1)
                 && (dataIn[i](IndexIn(i)) < dataOut[i](0)))
            IndexIn(i)++;
          IndexIn0(i) = IndexIn(i);
          Coeff0(i) = (dataOut[i](0) - dataIn[i](IndexIn(i) - 1)) /
            (dataIn[i](IndexIn(i)) - dataIn[i](IndexIn(i) - 1));
          Coeff(i) = Coeff0(i);
        }
      else
        {
          LengthIn(i) = dataIn.GetLength(i);
          LengthOut(i) = dataOut.GetLength(i);
          IndexOut(i) = 0;
          IndexIn(i) = 0;
          Coeff0(i) = 0;
          Coeff(i) = 0;
        }

    j = N - 1;
    for (i = 0; i < dataOut.GetNbElements(); i++)
      {

        if (j != dim)
          {
            while ((IndexIn(j) < LengthIn(j) - 1)
                   && (dataIn[j](IndexIn(j))
                       < dataOut[j](IndexOut(j))))
              IndexIn(j)++;

            Coeff(j) = (dataOut[j](IndexOut(j)) - dataIn[j](IndexIn(j) - 1)) /
              (dataIn[j](IndexIn(j)) - dataIn[j](IndexIn(j) - 1));
          }

        for (k = 0; k < int(pow(2.0, double(N)) + 0.5); k++)
          {

            l = k;
            for (m = 0; m < N; m++)
              {
                Pos_dim(m) = l % 2;
                l = l / 2;
              }
            Pos_dim(dim) = 0;

            IndexIn(dim) = 1;
            while ((IndexIn(dim) < LengthIn(dim) - 1)
                   && (dataIn[dim].Value(IndexIn(0) - Pos_dim(0), IndexIn(1) - Pos_dim(1),
                                         IndexIn(2) - Pos_dim(2), IndexIn(3) - Pos_dim(3),
                                         IndexIn(4) - Pos_dim(4), IndexIn(5) - Pos_dim(5),
                                         IndexIn(6) - Pos_dim(6), IndexIn(7) - Pos_dim(7),
                                         IndexIn(8) - Pos_dim(8), IndexIn(9) - Pos_dim(9))
                       < dataOut[dim].Value(IndexOut(0), IndexOut(1),
                                            IndexOut(2), IndexOut(3),
                                            IndexOut(4), IndexOut(5),
                                            IndexOut(6), IndexOut(7),
                                            IndexOut(8), IndexOut(9))))
              IndexIn(dim)++;

            IndexIn(dim)--;
            temp = dataIn[dim].Value(IndexIn(0) - Pos_dim(0), IndexIn(1) - Pos_dim(1),
                                     IndexIn(2) - Pos_dim(2), IndexIn(3) - Pos_dim(3),
                                     IndexIn(4) - Pos_dim(4), IndexIn(5) - Pos_dim(5),
                                     IndexIn(6) - Pos_dim(6), IndexIn(7) - Pos_dim(7),
                                     IndexIn(8) - Pos_dim(8), IndexIn(9) - Pos_dim(9));
            IndexIn(dim)++;

            Coeff(dim) = (dataOut[dim].Value(IndexOut(0), IndexOut(1),
                                             IndexOut(2), IndexOut(3),
                                             IndexOut(4), IndexOut(5),
                                             IndexOut(6), IndexOut(7),
                                             IndexOut(8), IndexOut(9))
                          - temp) /
              (dataIn[dim].Value(IndexIn(0) - Pos_dim(0), IndexIn(1) - Pos_dim(1),
                                 IndexIn(2) - Pos_dim(2), IndexIn(3) - Pos_dim(3),
                                 IndexIn(4) - Pos_dim(4), IndexIn(5) - Pos_dim(5),
                                 IndexIn(6) - Pos_dim(6), IndexIn(7) - Pos_dim(7),
                                 IndexIn(8) - Pos_dim(8), IndexIn(9) - Pos_dim(9)) - temp);

            l = k;
            coeff = TIn(1);
            for (m = 0; m < N; m++)
              {
                Pos(m) = l % 2;
                if (l % 2 == 1)
                  coeff *= TIn(1) - Coeff(m);
                else
                  coeff *= Coeff(m);
                l = l / 2;
              }

            dataOut.Value(IndexOut(0), IndexOut(1), IndexOut(2),
                          IndexOut(3), IndexOut(4), IndexOut(5),
                          IndexOut(6), IndexOut(7), IndexOut(8),
                          IndexOut(9)) +=
              TOut(coeff *
                   dataIn.Value(IndexIn(0) - Pos(0),
                                IndexIn(1) - Pos(1),
                                IndexIn(2) - Pos(2),
                                IndexIn(3) - Pos(3),
                                IndexIn(4) - Pos(4),
                                IndexIn(5) - Pos(5),
                                IndexIn(6) - Pos(6),
                                IndexIn(7) - Pos(7),
                                IndexIn(8) - Pos(8),
                                IndexIn(9) - Pos(9)));

          }

        j = N - 1;
        while ((j >= 0) && (IndexOut(j) == LengthOut(j) - 1))
          {
            IndexOut(j) = 0;
            IndexIn(j) = IndexIn0(j);
            Coeff(j) = Coeff0(j);
            j--;
          }
        if (j != -1)
          IndexOut(j)++;

      }

  }

  // ONEGENERAL //
  ////////////////


  /////////////////////////
  // ONEGENERALGETCOEFFS //

  //! Saves coefficients and indices of linear interpolation from data defined
  //! on regular grids, but one grid, to data defined on general grids.
  /*!
    Saves indices and coefficients of linear interpolation from data defined on regular
    grids, except one grid which may be a general grid (i.e. depending on other
    coordinates).
    Input data may be defined on a general grid along dimension 'dim', but only along
    this dimension. Output data may be defined on general or regular grids.
    Moreover, input data can still be defined on regular grids along dimension 'dim'.
    \param dataIn reference data.
    \param dataOut data to be interpolated.
    \param dim dimension related to the general grid.
    \param FileName file into which interpolation coefficients and indices are to be stored (on exit).
    \note dataOut is not modified. Interpolation oefficients and indices are just computed,
    not used for interpolation.
  */
  template < int N, class TIn, class TGIn,
             class TOut, class TGOut >
  void LinearInterpolationOneGeneralGetCoeffs(Data<TIn, N, TGIn>& dataIn,
                                              Data<TOut, N, TGOut>& dataOut,
                                              int dim, string FileName)
  {
    ofstream FileStream;
    FileStream.open(FileName.c_str(), ofstream::binary | ofstream::app);
#ifdef SELDONDATA_DEBUG_CHECK_IO
    // Checks if the file was opened.
    if (!FileStream.is_open())
      throw IOError("SeldonData::LinearInterpolationOneGeneralGetCoeffs(Data<TIn, N, TGIn>& dataIn, Data<TOut, N, TGOut>& dataOut, int dim, string FileName)",
                    "Unable to open file \"" + FileName + "\".");
#endif

    LinearInterpolationOneGeneralGetCoeffs(dataIn, dataOut, dim, FileStream);
    FileStream.close();
  }


  //! Saves coefficients and indices of linear interpolation from data defined
  //! on regular grids, but one grid, to data defined on general grids.
  /*!
    Saves indices and coefficients of linear interpolation from data defined on regular
    grids, except one grid which may be a general grid (i.e. depending on other
    coordinates).
    Input data may be defined on a general grid along dimension 'dim', but only along
    this dimension. Output data may be defined on general or regular grids.
    Moreover, input data can still be defined on regular grids along dimension 'dim'.
    \param dataIn reference data.
    \param dataOut data to be interpolated.
    \param dim dimension related to the general grid.
    \param FileStream stream into which interpolation coefficients and indices are to be stored (on exit).
    \note dataOut is not modified. Interpolation oefficients and indices are just computed,
    not used for interpolation.
  */
  template < int N, class TIn, class TGIn,
             class TOut, class TGOut >
  void LinearInterpolationOneGeneralGetCoeffs(Data<TIn, N, TGIn>& dataIn,
                                              Data<TOut, N, TGOut>& dataOut,
                                              int dim, ofstream& FileStream)
  {

#ifdef SELDONDATA_DEBUG_CHECK_IO
    // Checks if the file is ready.
    if (!FileStream.good())
      throw IOError("SeldonData::LinearInterpolationOneGeneralGetCoeffs(Data<TIn, N, TGIn>& dataIn, Data<TOut, N, TGOut>& dataOut, int dim, ofstream& FileStream",
                    "File is not ready.");
#endif

    Array<TIn, 2> RegularCoeffs;
    Array<TIn, 2> GeneralCoeffs;
    Array<int, 2> RegularIndices;
    Array<int, 2> GeneralIndices;
    FormatBinary<int> IntOutput;
    FormatBinary<TIn> TInOutput;
    LinearInterpolationOneGeneralGetCoeffs(dataIn, dataOut, dim,
                                           RegularCoeffs, GeneralCoeffs,
                                           RegularIndices, GeneralIndices);
    TInOutput.Write(RegularCoeffs, FileStream);
    TInOutput.Write(GeneralCoeffs, FileStream);
    IntOutput.Write(RegularIndices, FileStream);
    IntOutput.Write(GeneralIndices, FileStream);
  }


  //! Saves coefficients and indices of linear interpolation from data defined
  //! on regular grids, but one grid, to data defined on general grids.
  /*!
    Saves indices and coefficients of linear interpolation from data defined on regular
    grids, except one grid which may be a general grid (i.e. depending on other
    coordinates).
    Input data may be defined on a general grid along dimension 'dim', but only along
    this dimension. Output data may be defined on general or regular grids.
    Moreover, input data can still be defined on regular grids along dimension 'dim'.
    \param dataIn reference data.
    \param dataOut interpolated data (on exit).
    \param dim dimension related to the general grid.
    \param RegularCoeffs interpolation coefficients associated with regular grids (on exit).
    \param GeneralCoeffs interpolation coefficients associated with the input general grid (on exit).
    \param RegularIndices interpolation indices associated with regular grids (on exit).
    \param GeneralIndices interpolation indices associated with the input general grid (on exit).
    \note dataOut is not modified. Interpolation oefficients and indices are just computed,
    not used for interpolation.
  */
  template < int N, class TIn, class TGIn,
             class TOut, class TGOut >
  void LinearInterpolationOneGeneralGetCoeffs(Data<TIn, N, TGIn>& dataIn,
                                              Data<TOut, N, TGOut>& dataOut,
                                              int dim,
                                              Array<TIn, 2>& RegularCoeffs,
                                              Array<TIn, 2>& GeneralCoeffs,
                                              Array<int, 2>& RegularIndices,
                                              Array<int, 2>& GeneralIndices)
  {

    int i, j, k, l, m;
    Array<int, 1> IndexIn(10), IndexOut(10);
    Array<int, 1> LengthIn(10), LengthOut(10);
    Array<TIn, 1> Coeff(10);
    TIn temp;
    Array<bool, 1> Pos(10), Pos_dim(10);
    int NbElementsOut;

    dataOut.SetZero();
    NbElementsOut = dataOut.GetNbElements();

    RegularIndices.resize(NbElementsOut, N - 1);
    RegularIndices = 0;
    RegularCoeffs.resize(NbElementsOut, N - 1);
    RegularCoeffs = 0;
    GeneralIndices.resize(NbElementsOut, int(pow(2., N - 1) + 0.5));
    GeneralIndices = 0;
    GeneralCoeffs.resize(NbElementsOut, int(pow(2., N - 1) + 0.5));
    GeneralCoeffs = 0;

    for (i = 0; i < 10; i++)
      {
        Pos(i) = 0;
        Pos_dim(i) = 0;
      }

    for (i = 0; i < 10; i++)
      if (i < N)
        {
          LengthIn(i) = dataIn.GetLength(i);
          LengthOut(i) = dataOut.GetLength(i);
          IndexOut(i) = 0;
        }

    for (i = 0; i < NbElementsOut; i++)
      {
        // Loops over every regular dimensions to compute coefficients.
        for (j = 0; j < N; j++)
          if (j != dim)
            {
              // Searches for nearest coordinates beginning from 0.
              IndexIn(j) = 0;
              while ((IndexIn(j) < LengthIn(j) - 1)
                     && (dataIn[j](IndexIn(j))
                         < dataOut[j](IndexOut(j))))
                IndexIn(j)++;
              IndexIn(j) = max(1, IndexIn(j));

              Coeff(j) = (dataOut[j](IndexOut(j)) - dataIn[j](IndexIn(j) - 1)) /
                (dataIn[j](IndexIn(j)) - dataIn[j](IndexIn(j) - 1));

              // Saves coeffs and indices (values for dim are not saved).
              if (j < dim)
                {
                  RegularCoeffs(i, j) = Coeff(j);
                  RegularIndices(i, j) = IndexIn(j);
                }
              else
                {
                  RegularCoeffs(i, j - 1) = Coeff(j);
                  RegularIndices(i, j - 1) = IndexIn(j);
                }
            }

        for (k = 0; k < int(pow(2.0, double(N - 1)) + 0.5); k++)
          {
            l = k;
            for (m = 0; m < N; m++)
              if (m != dim)
                {
                  Pos_dim(m) = l % 2;
                  l = l / 2;
                }
            Pos_dim(dim) = 0;

            IndexIn(dim) = 1;
            while ((IndexIn(dim) < LengthIn(dim) - 1)
                   && (dataIn[dim].Value(IndexIn(0) - Pos_dim(0), IndexIn(1) - Pos_dim(1),
                                         IndexIn(2) - Pos_dim(2), IndexIn(3) - Pos_dim(3),
                                         IndexIn(4) - Pos_dim(4), IndexIn(5) - Pos_dim(5),
                                         IndexIn(6) - Pos_dim(6), IndexIn(7) - Pos_dim(7),
                                         IndexIn(8) - Pos_dim(8), IndexIn(9) - Pos_dim(9))
                       < dataOut[dim].Value(IndexOut(0), IndexOut(1),
                                            IndexOut(2), IndexOut(3),
                                            IndexOut(4), IndexOut(5),
                                            IndexOut(6), IndexOut(7),
                                            IndexOut(8), IndexOut(9))))
              IndexIn(dim)++;

            IndexIn(dim)--;
            temp = dataIn[dim].Value(IndexIn(0) - Pos_dim(0), IndexIn(1) - Pos_dim(1),
                                     IndexIn(2) - Pos_dim(2), IndexIn(3) - Pos_dim(3),
                                     IndexIn(4) - Pos_dim(4), IndexIn(5) - Pos_dim(5),
                                     IndexIn(6) - Pos_dim(6), IndexIn(7) - Pos_dim(7),
                                     IndexIn(8) - Pos_dim(8), IndexIn(9) - Pos_dim(9));
            IndexIn(dim)++;

            Coeff(dim) = (dataOut[dim].Value(IndexOut(0), IndexOut(1),
                                             IndexOut(2), IndexOut(3),
                                             IndexOut(4), IndexOut(5),
                                             IndexOut(6), IndexOut(7),
                                             IndexOut(8), IndexOut(9))
                          - temp) /
              (dataIn[dim].Value(IndexIn(0) - Pos_dim(0), IndexIn(1) - Pos_dim(1),
                                 IndexIn(2) - Pos_dim(2), IndexIn(3) - Pos_dim(3),
                                 IndexIn(4) - Pos_dim(4), IndexIn(5) - Pos_dim(5),
                                 IndexIn(6) - Pos_dim(6), IndexIn(7) - Pos_dim(7),
                                 IndexIn(8) - Pos_dim(8), IndexIn(9) - Pos_dim(9)) - temp);
            // Saves coeffs and indices (the index after).
            GeneralCoeffs(i, k) = Coeff(dim);
            GeneralIndices(i, k) = IndexIn(dim);
          }

        j = N - 1;
        while ((j >= 0) && (IndexOut(j) == LengthOut(j) - 1))
          {
            IndexOut(j) = 0;
            j--;
          }
        if (j != -1)
          IndexOut(j)++;

      }

  }

  // ONEGENERALGETCOEFFS //
  /////////////////////////


  ///////////////////////
  // ONEGENERALCOMPUTE //

  //! Compute linear interpolation from data defined
  //! on regular grids, but one grid, to data defined on general grids,
  //! using indices and coefficients previously computed.
  /*!
    Performs linear interpolation from data defined on regular grids,
    except one grid which may be a general grid.
    Input data may be defined on a general grid along dimension 'dim', but only along
    this dimension. Output data may be defined on general or regular grids.
    Moreover, input data can still be defined on regular grids along dimension 'dim'.
    \param dataIn reference data.
    \param dim dimension related to the general grid.
    \param FileName file into which interpolation coefficients and indices are stored.
    \param dataOut interpolated data (on exit).
  */
  template < int N, class TIn, class TGIn,
             class TOut, class TGOut >
  void LinearInterpolationOneGeneralCompute(Data<TIn, N, TGIn>& dataIn,
                                            int dim, string FileName,
                                            Data<TOut, N, TGOut>& dataOut)
  {
    ifstream FileStream;
    FileStream.open(FileName.c_str(), ifstream::binary);
#ifdef SELDONDATA_DEBUG_CHECK_IO
    // Checks if the file was opened.
    if (!FileStream.is_open())
      throw IOError("SeldonData::LinearInterpolationOneGeneralCompute(Data<TIn, N, TGIn>& dataIn, int dim, string FileName, Data<TOut, N, TGOut>& dataOut)",
                    "Unable to open file \"" + FileName + "\".");
#endif

    LinearInterpolationOneGeneralCompute(dataIn, dim, FileStream, dataOut);
    FileStream.close();
  }


  //! Compute linear interpolation from data defined
  //! on regular grids, but one grid, to data defined on general grids,
  //! using indices and coefficients previously computed.
  /*!
    Performs linear interpolation from data defined on regular grids,
    except one grid which may be a general grid.
    Input data may be defined on a general grid along dimension 'dim', but only along
    this dimension. Output data may be defined on general or regular grids.
    Moreover, input data can still be defined on regular grids along dimension 'dim'.
    \param dataIn reference data.
    \param dim dimension related to the general grid.
    \param FileStream stream into which interpolation coefficients and indices are stored.
    \param dataOut interpolated data (on exit).
  */
  template < int N, class TIn, class TGIn,
             class TOut, class TGOut >
  void LinearInterpolationOneGeneralCompute(Data<TIn, N, TGIn>& dataIn,
                                            int dim, ifstream& FileStream,
                                            Data<TOut, N, TGOut>& dataOut)
  {
#ifdef SELDONDATA_DEBUG_CHECK_IO
    // Checks if the file is ready.
    if (!FileStream.good())
      throw IOError("SeldonData::LinearInterpolationOneGeneralCompute(Data<TIn, N, TGIn>& dataIn, int dim, ifstream& FileStream, Data<TOut, N, TGOut>& dataOut)",
                    "File is not ready.");
#endif

    int NbElementsOut;
    FormatBinary<int> IntInput;
    FormatBinary<TIn> TInInput;

    Array<TIn, 2> RegularCoeffs;
    Array<TIn, 2> GeneralCoeffs;
    Array<int, 2> RegularIndices;
    Array<int, 2> GeneralIndices;

    NbElementsOut = dataOut.GetNbElements();

    // Resizes arrays to the right size.
    RegularIndices.resize(NbElementsOut, N - 1);
    RegularIndices = 0;
    RegularCoeffs.resize(NbElementsOut, N - 1);
    RegularCoeffs = 0;
    GeneralIndices.resize(NbElementsOut, int(pow(2., N - 1) + 0.5));
    GeneralIndices = 0;
    GeneralCoeffs.resize(NbElementsOut, int(pow(2., N - 1) + 0.5));
    GeneralCoeffs = 0;

    TInInput.Read(FileStream, RegularCoeffs);
    TInInput.Read(FileStream, GeneralCoeffs);
    IntInput.Read(FileStream, RegularIndices);
    IntInput.Read(FileStream, GeneralIndices);

    LinearInterpolationOneGeneralCompute(dataIn, dim,
                                         RegularCoeffs, GeneralCoeffs,
                                         RegularIndices, GeneralIndices,
                                         dataOut);
  }


  //! Compute linear interpolation from data defined
  //! on regular grids, but one grid, to data defined on general grids,
  //! using indices and coefficients previously computed.
  /*!
    Performs linear interpolation from data defined on regular grids,
    except one grid which may be a general grid.
    Input data may be defined on a general grid along dimension 'dim', but only along
    this dimension. Output data may be defined on general or regular grids.
    Moreover, input data can still be defined on regular grids along dimension 'dim'.
    \param dataIn reference data.
    \param dim dimension related to the general grid.
    \param RegularCoeffs interpolation coefficients associated with regular grids.
    \param GeneralCoeffs interpolation coefficients associated with the input general grid.
    \param RegularIndices interpolation indices associated with regular grids.
    \param GeneralIndices interpolation indices associated with the input general grid.
    \param dataOut interpolated data (on exit).
  */
  template < int N, class TIn, class TGIn,
             class TOut, class TGOut >
  void LinearInterpolationOneGeneralCompute(Data<TIn, N, TGIn>& dataIn,
                                            int dim,
                                            Array<TIn, 2>& RegularCoeffs,
                                            Array<TIn, 2>& GeneralCoeffs,
                                            Array<int, 2>& RegularIndices,
                                            Array<int, 2>& GeneralIndices,
                                            Data<TOut, N, TGOut>& dataOut)
  {
    int i, j, k, l, ldim, m;
    Array<int, 1> IndexIn(10), IndexOut(10);
    Array<int, 1> LengthOut(10);
    Array<TIn, 1> Coeff(10);
    TIn coeff;

    int limit = int(pow(2.0, double(dim)) + 0.5);

    IndexOut = 0;

    for (i = 0; i < 10; i++)
      if (i < N)
        LengthOut(i) = dataOut.GetLength(i);

    dataOut.SetZero();

    for (i = 0; i < dataOut.GetNbElements(); i++)
      {
        for (k = 0; k < int(pow(2.0, double(N)) + 0.5); k++)
          {
            l = k;
            if (k < limit)
              ldim = k;
            else
              ldim = (k / (limit * 2)) * limit + k % limit;

            coeff = TIn(1);
            for (m = 0; m < N; m++)
              {
                if (m < dim)
                  {
                    if (l % 2 == 0)
                      {
                        IndexIn(m) = RegularIndices(i, m);
                        coeff *= RegularCoeffs(i, m);
                      }
                    else
                      {
                        IndexIn(m) = RegularIndices(i, m) - 1;
                        coeff *= TIn(1) - RegularCoeffs(i, m);
                      }
                    l = l / 2;
                  }
                else if (m > dim)
                  {
                    if (l % 2 == 0)
                      {
                        IndexIn(m) = RegularIndices(i, m - 1);
                        coeff *= RegularCoeffs(i, m - 1);
                      }
                    else
                      {
                        IndexIn(m) = RegularIndices(i, m - 1) - 1;
                        coeff *= TIn(1) - RegularCoeffs(i, m - 1);
                      }
                    l = l / 2;
                  }
                else
                  {
                    if (l % 2 == 1)
                      {
                        IndexIn(m) = GeneralIndices(i, ldim);
                        coeff *= GeneralCoeffs(i, ldim);
                      }
                    else
                      {
                        IndexIn(m) = GeneralIndices(i, ldim) - 1;
                        coeff *= TIn(1) - GeneralCoeffs(i, ldim);
                      }
                    l = l / 2;
                  }
              }

            dataOut.Value(IndexOut(0), IndexOut(1), IndexOut(2),
                          IndexOut(3), IndexOut(4), IndexOut(5),
                          IndexOut(6), IndexOut(7), IndexOut(8),
                          IndexOut(9)) +=
              TOut(coeff *
                   dataIn.Value(IndexIn(0), IndexIn(1),
                                IndexIn(2), IndexIn(3),
                                IndexIn(4), IndexIn(5),
                                IndexIn(6), IndexIn(7),
                                IndexIn(8), IndexIn(9)));

          }

        j = N - 1;
        while ((j >= 0) && (IndexOut(j) == LengthOut(j) - 1))
          {
            IndexOut(j) = 0;
            j--;
          }
        if (j != -1)
          IndexOut(j)++;
      }
  }

  // ONEGENERALCOMPUTE //
  ///////////////////////


  ///////////////
  // DIMENSION //

  //! Linear interpolation along a given dimension.
  /*!
    Linear interpolation only along a given dimension.
    \param dataIn reference data.
    \param dataOut interpolated data (on exit).
    \param dim dimension along which data should be interpolated.
    \note 'dataIn' and 'dataOut' are assumed to be defined on
    the same grids except on the grid related to dimension 'dim'.
  */
  template < int N, class TIn, class TGIn,
             class TOut, class TGOut >
  void LinearInterpolationDimension(Data<TIn, N, TGIn>& dataIn,
                                    Data<TOut, N, TGOut>& dataOut, int dim)
  {

    int i, j, k;
    Array<int, 1> IndexIn(10), IndexOut(10);
    Array<int, 1> Length(10);
    TIn coeff, coord;
    int Ndim_in, Ndim_out;

    for (i = 0; i < 10; i++)
      {
        Length(i) = dataOut.GetLength(i);
        IndexOut(i) = 0;
        IndexIn(i) = 0;
      }

#ifdef SELDONDATA_DEBUG_CHECK_DIMENSIONS
    // Checks whether dimensions match.
    for (i = 0; i < 10; i++)
      if ((i != dim) && (Length(i) != dataIn.GetLength(i)))
        throw WrongDim("LinearInterpolationDimension(Data&, Data&, " + to_str(dim) + ")",
                       "Along dimension #" + to_str(i) + ", input data has "
                       + to_str(dataIn.GetLength(i)) + " elements and ouput data"
                       + " has " + to_str(Length(i)) + " elements. There must be the same"
                       + " number of elements.");
#endif

    Ndim_in = dataIn.GetLength(dim);
    Ndim_out = dataOut.GetLength(dim);

    for (i = 0; i < dataOut.GetNbElements() / Ndim_out; i++)
      {

        IndexIn(dim) = 1;
        for (k = 0; k < Ndim_out; k++)
          {
            IndexOut(dim) = k;
            coord = dataOut[dim].Value(IndexOut(0), IndexOut(1), IndexOut(2),
                                       IndexOut(3), IndexOut(4), IndexOut(5),
                                       IndexOut(6), IndexOut(7), IndexOut(8),
                                       IndexOut(9));
            while ((IndexIn(dim) < Ndim_in)
                   && (dataIn[dim].Value(IndexIn(0), IndexIn(1), IndexIn(2),
                                         IndexIn(3), IndexIn(4), IndexIn(5),
                                         IndexIn(6), IndexIn(7), IndexIn(8),
                                         IndexIn(9)) < coord))
              IndexIn(dim)++;

            if (IndexIn(dim) == Ndim_in)
              IndexIn(dim) = Ndim_in - 1;

            IndexIn(dim)--;
            coord = dataIn[dim].Value(IndexIn(0), IndexIn(1), IndexIn(2),
                                      IndexIn(3), IndexIn(4), IndexIn(5),
                                      IndexIn(6), IndexIn(7), IndexIn(8),
                                      IndexIn(9));
            IndexIn(dim)++;

            coeff = (dataOut[dim].Value(IndexOut(0), IndexOut(1), IndexOut(2),
                                        IndexOut(3), IndexOut(4), IndexOut(5),
                                        IndexOut(6), IndexOut(7), IndexOut(8),
                                        IndexOut(9)) - coord)
              / (dataIn[dim].Value(IndexIn(0), IndexIn(1), IndexIn(2),
                                   IndexIn(3), IndexIn(4), IndexIn(5),
                                   IndexIn(6), IndexIn(7), IndexIn(8),
                                   IndexIn(9)) - coord);

            dataOut.Value(IndexOut(0), IndexOut(1), IndexOut(2),
                          IndexOut(3), IndexOut(4), IndexOut(5),
                          IndexOut(6), IndexOut(7), IndexOut(8),
                          IndexOut(9)) =
              coeff * dataIn.Value(IndexIn(0), IndexIn(1), IndexIn(2),
                                   IndexIn(3), IndexIn(4), IndexIn(5),
                                   IndexIn(6), IndexIn(7), IndexIn(8),
                                   IndexIn(9));

            IndexIn(dim)--;
            dataOut.Value(IndexOut(0), IndexOut(1), IndexOut(2),
                          IndexOut(3), IndexOut(4), IndexOut(5),
                          IndexOut(6), IndexOut(7), IndexOut(8),
                          IndexOut(9)) +=
              (TIn(1.0) - coeff) * dataIn.Value(IndexIn(0), IndexIn(1), IndexIn(2),
                                                IndexIn(3), IndexIn(4), IndexIn(5),
                                                IndexIn(6), IndexIn(7), IndexIn(8),
                                                IndexIn(9));
            IndexIn(dim)++;
          }

        j = N - 1;
        while ((j >= 0) && ((IndexOut(j) == Length(j) - 1) || (j == dim)))
          {
            IndexOut(j) = 0;
            IndexIn(j) = 0;
            j--;
          }
        if (j != -1)
          {
            IndexOut(j)++;
            IndexIn(j)++;
          }

      }

  }

  // DIMENSION //
  ///////////////


  ///////////
  // POINT //

  //! Linear interpolation from data defined on regular grids to one point.
  /*!
    Linear interpolation from data defined on regular grids to
    data on one point.
    \param dataIn reference data.
    \param Coord coordinates of the interpolation point.
    \return data value on the interpolation point.
  */
  template < int N, class TIn, class TGIn,
             class TOut, class TGOut >
  void LinearInterpolationPoint(Data<TIn, N, TGIn>& dataIn,
                                Array<TGOut, 1>& Coord, TOut& dataOut)
  {

    if (Coord.extent(0) != N)
      throw WrongDim("LinearInterpolationPoint(Data&, Array&)",
                     "There are " + to_str(Coord.extent(0))
                     + " coordinates for output point "
                     + "and " + to_str(N) + " dimensions for input data. "
                     + "These numbers must be equal.");
    int i, k, l, m;
    Array<int, 1> IndexIn(10);
    Array<int, 1> LengthIn(10);
    Array<TIn, 1> Coeff(10);
    Array<bool, 1> Pos(10);
    TIn coeff;

    for (i = 0; i < 10; i++)
      {
        Pos(i) = 0;
        LengthIn(i) = dataIn.GetLength(i);
        IndexIn(i) = 0;
        if (i < N)
          {
            while ((IndexIn(i) < LengthIn(i))
                   && (dataIn[i](IndexIn(i)) < Coord(i)))
              IndexIn(i)++;
            if (IndexIn(i) == LengthIn(i))
              IndexIn(i) = LengthIn(i) - 1;
            else if (IndexIn(i) == 0)
              IndexIn(i) = 1;
            if (LengthIn(i) != 0)
              Coeff(i) = (Coord(i) - dataIn[i](IndexIn(i) - 1)) /
                (dataIn[i](IndexIn(i)) - dataIn[i](IndexIn(i) - 1));
            else
              Coeff(i) = TIn(0);
          }
      }

    dataOut = TOut(0);
    for (k = 0; k < int(pow(2.0, double(N))); k++)
      {
        l = k;
        coeff = TIn(1);
        for (m = 0; m < N; m++)
          {
            Pos(m) = l % 2;
            if (l % 2 == 1)
              coeff *= TIn(1) - Coeff(m);
            else
              coeff *= Coeff(m);
            l = l / 2;
          }

        dataOut += TOut(coeff *
                        dataIn.Value(IndexIn(0) - Pos(0),
                                     IndexIn(1) - Pos(1),
                                     IndexIn(2) - Pos(2),
                                     IndexIn(3) - Pos(3),
                                     IndexIn(4) - Pos(4),
                                     IndexIn(5) - Pos(5),
                                     IndexIn(6) - Pos(6),
                                     IndexIn(7) - Pos(7),
                                     IndexIn(8) - Pos(8),
                                     IndexIn(9) - Pos(9)));
      }
    // return dataOut;
  }

  // POINT //
  ///////////


}  // namespace SeldonData.

#define FILE_SELDONDATA_FUNCTIONS_CXX
#endif
