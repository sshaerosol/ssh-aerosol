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

#ifndef FILE_SELDONDATA_DATA_CXX

#include <cstdlib>
#include <limits>
#include "Data.hxx"

namespace SeldonData
{

  //! Default constructor.
  template<class T, int N, class TG>
  Data<T, N, TG>::Data()  throw():
    grids_(N)
  {
    for (int i = 0; i < N; i++)
      grids_(i) = NULL;
  }

  //! Constructor for one-dimensional data.
  /*!
    \param G0 grid for dimension #0.
  */
  template<class T, int N, class TG>
  Data<T, N, TG>::Data(Grid<TG>& G0)  throw():
    data_(G0.GetLength()), grids_(N)
  {

#ifdef SELDONDATA_DEBUG_CHECK_DIMENSIONS
    if (N != 1)
      {
        cout << "ERROR!" << endl;
        cout << "Wrong dimension in Data<T, N, TG>::Data(Grid<TG>, ...)."
             << endl << "   Required " + to_str(N)
          + " grids, but got 1 grid." << endl;
        abort();
      }
#endif

    grids_(0) = G0.Copy();
    SetVariables();
    InitData();

#ifdef SELDONDATA_DEBUG_CHECK_DIMENSIONS
    int i, j;
    for (i = 0; i < N; i++)
      for (j = 0; j < N; j++)
        if ((grids_(i)->GetLength(j) != 0)
            && (grids_(i)->GetLength(j) != grids_(j)->GetLength()))
          {
            cout << "ERROR!" << endl;
            cout << "Wrong dimension in Data<T, N, TG>::Data(Grid<TG>, ...)."
                 << endl << "   Length of grid #" << i << " along dimension #"
                 << j << " is " << grids_(i)->GetLength(j) << " but should be "
                 << grids_(j)->GetLength() << "." << endl;
            abort();
          }
#endif

  }

  //! Constructor for two-dimensional data.
  /*!
    \param G0 grid for dimension #0.
    \param G1 grid for dimension #1.
  */
  template<class T, int N, class TG>
  Data<T, N, TG>::Data(Grid<TG>& G0, Grid<TG>& G1)  throw():
    data_(G0.GetLength(), G1.GetLength()), grids_(N)
  {

#ifdef SELDONDATA_DEBUG_CHECK_DIMENSIONS
    if (N != 2)
      {
        cout << "ERROR!" << endl;
        cout << "Wrong dimension in Data<T, N, TG>::Data(Grid<TG>, ...)."
             << endl << "   Required " + to_str(N)
          + " grids, but got 2 grids." << endl;
        abort();
      }
#endif

    grids_(0) = G0.Copy();
    grids_(1) = G1.Copy();
    SetVariables();
    InitData();

#ifdef SELDONDATA_DEBUG_CHECK_DIMENSIONS
    int i, j;
    for (i = 0; i < N; i++)
      for (j = 0; j < N; j++)
        if ((grids_(i)->GetLength(j) != 0)
            && (grids_(i)->GetLength(j) != grids_(j)->GetLength()))
          {
            cout << "ERROR!" << endl;
            cout << "Wrong dimension in Data<T, N, TG>::Data(Grid<TG>, ...)."
                 << endl << "   Length of grid #" << i << " along dimension #"
                 << j << " is " << grids_(i)->GetLength(j) << " but should be "
                 << grids_(j)->GetLength() << "." << endl;
            abort();
          }
#endif

  }

  //! Constructor for three-dimensional data.
  /*!
    \param G0 grid for dimension #0.
    \param G1 grid for dimension #1.
    \param G2 grid for dimension #2.
  */
  template<class T, int N, class TG>
  Data<T, N, TG>::Data(Grid<TG>& G0, Grid<TG>& G1, Grid<TG>& G2)  throw():
    data_(G0.GetLength(), G1.GetLength(),
          G2.GetLength()), grids_(N)
  {

#ifdef SELDONDATA_DEBUG_CHECK_DIMENSIONS
    if (N != 3)
      {
        cout << "ERROR!" << endl;
        cout << "Wrong dimension in Data<T, N, TG>::Data(Grid<TG>, ...)."
             << endl << "   Required " + to_str(N)
          + " grids, but got 3 grids." << endl;
        abort();
      }
#endif

    grids_(0) = G0.Copy();
    grids_(1) = G1.Copy();
    grids_(2) = G2.Copy();
    SetVariables();
    InitData();


#ifdef SELDONDATA_DEBUG_CHECK_DIMENSIONS
    int i, j;
    for (i = 0; i < N; i++)
      for (j = 0; j < N; j++)
        if ((grids_(i)->GetLength(j) != 0)
            && (grids_(i)->GetLength(j) != grids_(j)->GetLength()))
          {
            cout << "ERROR!" << endl;
            cout << "Wrong dimension in Data<T, N, TG>::Data(Grid<TG>, ...)."
                 << endl << "   Length of grid #" << i << " along dimension #"
                 << j << " is " << grids_(i)->GetLength(j) << " but should be "
                 << grids_(j)->GetLength() << "." << endl;
            abort();
          }
#endif

  }

  //! Constructor for four-dimensional data.
  /*!
    \param G0 grid for dimension #0.
    \param G1 grid for dimension #1.
    \param G2 grid for dimension #2.
    \param G3 grid for dimension #3.
  */
  template<class T, int N, class TG>
  Data<T, N, TG>::Data(Grid<TG>& G0, Grid<TG>& G1, Grid<TG>& G2,
                       Grid<TG>& G3)  throw():
    data_(G0.GetLength(), G1.GetLength(),
          G2.GetLength(), G3.GetLength()), grids_(N)
  {

#ifdef SELDONDATA_DEBUG_CHECK_DIMENSIONS
    if (N != 4)
      {
        cout << "ERROR!" << endl;
        cout << "Wrong dimension in Data<T, N, TG>::Data(Grid<TG>, ...)."
             << endl << "   Required " + to_str(N)
          + " grids, but got 4 grids." << endl;
        abort();
      }
#endif

    grids_(0) = G0.Copy();
    grids_(1) = G1.Copy();
    grids_(2) = G2.Copy();
    grids_(3) = G3.Copy();
    SetVariables();
    InitData();


#ifdef SELDONDATA_DEBUG_CHECK_DIMENSIONS
    int i, j;
    for (i = 0; i < N; i++)
      for (j = 0; j < N; j++)
        if ((grids_(i)->GetLength(j) != 0)
            && (grids_(i)->GetLength(j) != grids_(j)->GetLength()))
          {
            cout << "ERROR!" << endl;
            cout << "Wrong dimension in Data<T, N, TG>::Data(Grid<TG>, ...)."
                 << endl << "   Length of grid #" << i << " along dimension #"
                 << j << " is " << grids_(i)->GetLength(j) << " but should be "
                 << grids_(j)->GetLength() << "." << endl;
            abort();
          }
#endif

  }

  //! Constructor for five-dimensional data.
  /*!
    \param G0 grid for dimension #0.
    \param G1 grid for dimension #1.
    \param G2 grid for dimension #2.
    \param G3 grid for dimension #3.
    \param G4 grid for dimension #4.
  */
  template<class T, int N, class TG>
  Data<T, N, TG>::Data(Grid<TG>& G0, Grid<TG>& G1, Grid<TG>& G2,
                       Grid<TG>& G3, Grid<TG>& G4)  throw():
    data_(G0.GetLength(), G1.GetLength(),
          G2.GetLength(), G3.GetLength(),
          G4.GetLength()), grids_(N)
  {

#ifdef SELDONDATA_DEBUG_CHECK_DIMENSIONS
    if (N != 5)
      {
        cout << "ERROR!" << endl;
        cout << "Wrong dimension in Data<T, N, TG>::Data(Grid<TG>, ...)."
             << endl << "   Required " + to_str(N)
          + " grids, but got 5 grids." << endl;
        abort();
      }
#endif

    grids_(0) = G0.Copy();
    grids_(1) = G1.Copy();
    grids_(2) = G2.Copy();
    grids_(3) = G3.Copy();
    grids_(4) = G4.Copy();
    SetVariables();
    InitData();


#ifdef SELDONDATA_DEBUG_CHECK_DIMENSIONS
    int i, j;
    for (i = 0; i < N; i++)
      for (j = 0; j < N; j++)
        if ((grids_(i)->GetLength(j) != 0)
            && (grids_(i)->GetLength(j) != grids_(j)->GetLength()))
          {
            cout << "ERROR!" << endl;
            cout << "Wrong dimension in Data<T, N, TG>::Data(Grid<TG>, ...)."
                 << endl << "   Length of grid #" << i << " along dimension #"
                 << j << " is " << grids_(i)->GetLength(j) << " but should be "
                 << grids_(j)->GetLength() << "." << endl;
            abort();
          }
#endif

  }

  //! Constructor for six-dimensional data.
  /*!
    \param G0 grid for dimension #0.
    \param G1 grid for dimension #1.
    \param G2 grid for dimension #2.
    \param G3 grid for dimension #3.
    \param G4 grid for dimension #4.
    \param G5 grid for dimension #5.
  */
  template<class T, int N, class TG>
  Data<T, N, TG>::Data(Grid<TG>& G0, Grid<TG>& G1, Grid<TG>& G2,
                       Grid<TG>& G3, Grid<TG>& G4, Grid<TG>& G5)  throw():
    data_(G0.GetLength(), G1.GetLength(),
          G2.GetLength(), G3.GetLength(), G4.GetLength(),
          G5.GetLength()), grids_(N)
  {

#ifdef SELDONDATA_DEBUG_CHECK_DIMENSIONS
    if (N != 6)
      {
        cout << "ERROR!" << endl;
        cout << "Wrong dimension in Data<T, N, TG>::Data(Grid<TG>, ...)."
             << endl << "   Required " + to_str(N)
          + " grids, but got 6 grids." << endl;
        abort();
      }
#endif

    grids_(0) = G0.Copy();
    grids_(1) = G1.Copy();
    grids_(2) = G2.Copy();
    grids_(3) = G3.Copy();
    grids_(4) = G4.Copy();
    grids_(5) = G5.Copy();
    SetVariables();
    InitData();


#ifdef SELDONDATA_DEBUG_CHECK_DIMENSIONS
    int i, j;
    for (i = 0; i < N; i++)
      for (j = 0; j < N; j++)
        if ((grids_(i)->GetLength(j) != 0)
            && (grids_(i)->GetLength(j) != grids_(j)->GetLength()))
          {
            cout << "ERROR!" << endl;
            cout << "Wrong dimension in Data<T, N, TG>::Data(Grid<TG>, ...)."
                 << endl << "   Length of grid #" << i << " along dimension #"
                 << j << " is " << grids_(i)->GetLength(j) << " but should be "
                 << grids_(j)->GetLength() << "." << endl;
            abort();
          }
#endif

  }

  //! Constructor for seven-dimensional data.
  /*!
    \param G0 grid for dimension #0.
    \param G1 grid for dimension #1.
    \param G2 grid for dimension #2.
    \param G3 grid for dimension #3.
    \param G4 grid for dimension #4.
    \param G5 grid for dimension #5.
    \param G6 grid for dimension #6.
  */
  template<class T, int N, class TG>
  Data<T, N, TG>::Data(Grid<TG>& G0, Grid<TG>& G1, Grid<TG>& G2,
                       Grid<TG>& G3, Grid<TG>& G4, Grid<TG>& G5,
                       Grid<TG>& G6)  throw():
    data_(G0.GetLength(), G1.GetLength(),
          G2.GetLength(), G3.GetLength(), G4.GetLength(),
          G5.GetLength(), G6.GetLength()), grids_(N)
  {

#ifdef SELDONDATA_DEBUG_CHECK_DIMENSIONS
    if (N != 7)
      {
        cout << "ERROR!" << endl;
        cout << "Wrong dimension in Data<T, N, TG>::Data(Grid<TG>, ...)."
             << endl << "   Required " + to_str(N)
          + " grids, but got 7 grids." << endl;
        abort();
      }
#endif

    grids_(0) = G0.Copy();
    grids_(1) = G1.Copy();
    grids_(2) = G2.Copy();
    grids_(3) = G3.Copy();
    grids_(4) = G4.Copy();
    grids_(5) = G5.Copy();
    grids_(6) = G6.Copy();
    SetVariables();
    InitData();


#ifdef SELDONDATA_DEBUG_CHECK_DIMENSIONS
    int i, j;
    for (i = 0; i < N; i++)
      for (j = 0; j < N; j++)
        if ((grids_(i)->GetLength(j) != 0)
            && (grids_(i)->GetLength(j) != grids_(j)->GetLength()))
          {
            cout << "ERROR!" << endl;
            cout << "Wrong dimension in Data<T, N, TG>::Data(Grid<TG>, ...)."
                 << endl << "   Length of grid #" << i << " along dimension #"
                 << j << " is " << grids_(i)->GetLength(j) << " but should be "
                 << grids_(j)->GetLength() << "." << endl;
            abort();
          }
#endif

  }

  //! Constructor for eight-dimensional data.
  /*!
    \param G0 grid for dimension #0.
    \param G1 grid for dimension #1.
    \param G2 grid for dimension #2.
    \param G3 grid for dimension #3.
    \param G4 grid for dimension #4.
    \param G5 grid for dimension #5.
    \param G6 grid for dimension #6.
    \param G7 grid for dimension #7.
  */
  template<class T, int N, class TG>
  Data<T, N, TG>::Data(Grid<TG>& G0, Grid<TG>& G1, Grid<TG>& G2,
                       Grid<TG>& G3, Grid<TG>& G4, Grid<TG>& G5,
                       Grid<TG>& G6, Grid<TG>& G7)  throw():
    data_(G0.GetLength(), G1.GetLength(),
          G2.GetLength(), G3.GetLength(), G4.GetLength(),
          G5.GetLength(), G6.GetLength(), G7.GetLength()),
    grids_(N)
  {

#ifdef SELDONDATA_DEBUG_CHECK_DIMENSIONS
    if (N != 8)
      {
        cout << "ERROR!" << endl;
        cout << "Wrong dimension in Data<T, N, TG>::Data(Grid<TG>, ...)."
             << endl << "   Required " + to_str(N)
          + " grids, but got 8 grids." << endl;
        abort();
      }
#endif

    grids_(0) = G0.Copy();
    grids_(1) = G1.Copy();
    grids_(2) = G2.Copy();
    grids_(3) = G3.Copy();
    grids_(4) = G4.Copy();
    grids_(5) = G5.Copy();
    grids_(6) = G6.Copy();
    grids_(7) = G7.Copy();
    SetVariables();
    InitData();


#ifdef SELDONDATA_DEBUG_CHECK_DIMENSIONS
    int i, j;
    for (i = 0; i < N; i++)
      for (j = 0; j < N; j++)
        if ((grids_(i)->GetLength(j) != 0)
            && (grids_(i)->GetLength(j) != grids_(j)->GetLength()))
          {
            cout << "ERROR!" << endl;
            cout << "Wrong dimension in Data<T, N, TG>::Data(Grid<TG>, ...)."
                 << endl << "   Length of grid #" << i << " along dimension #"
                 << j << " is " << grids_(i)->GetLength(j) << " but should be "
                 << grids_(j)->GetLength() << "." << endl;
            abort();
          }
#endif

  }

  //! Constructor for nine-dimensional data.
  /*!
    \param G0 grid for dimension #0.
    \param G1 grid for dimension #1.
    \param G2 grid for dimension #2.
    \param G3 grid for dimension #3.
    \param G4 grid for dimension #4.
    \param G5 grid for dimension #5.
    \param G6 grid for dimension #6.
    \param G7 grid for dimension #7.
    \param G8 grid for dimension #8.
  */
  template<class T, int N, class TG>
  Data<T, N, TG>::Data(Grid<TG>& G0, Grid<TG>& G1, Grid<TG>& G2,
                       Grid<TG>& G3, Grid<TG>& G4, Grid<TG>& G5,
                       Grid<TG>& G6, Grid<TG>& G7, Grid<TG>& G8)  throw():
    data_(G0.GetLength(), G1.GetLength(),
          G2.GetLength(), G3.GetLength(), G4.GetLength(),
          G5.GetLength(), G6.GetLength(), G7.GetLength(),
          G8.GetLength()), grids_(N)
  {

#ifdef SELDONDATA_DEBUG_CHECK_DIMENSIONS
    if (N != 9)
      {
        cout << "ERROR!" << endl;
        cout << "Wrong dimension in Data<T, N, TG>::Data(Grid<TG>, ...)."
             << endl << "   Required " + to_str(N)
          + " grids, but got 9 grids." << endl;
        abort();
      }
#endif

    grids_(0) = G0.Copy();
    grids_(1) = G1.Copy();
    grids_(2) = G2.Copy();
    grids_(3) = G3.Copy();
    grids_(4) = G4.Copy();
    grids_(5) = G5.Copy();
    grids_(6) = G6.Copy();
    grids_(7) = G7.Copy();
    grids_(8) = G8.Copy();
    SetVariables();
    InitData();


#ifdef SELDONDATA_DEBUG_CHECK_DIMENSIONS
    int i, j;
    for (i = 0; i < N; i++)
      for (j = 0; j < N; j++)
        if ((grids_(i)->GetLength(j) != 0)
            && (grids_(i)->GetLength(j) != grids_(j)->GetLength()))
          {
            cout << "ERROR!" << endl;
            cout << "Wrong dimension in Data<T, N, TG>::Data(Grid<TG>, ...)."
                 << endl << "   Length of grid #" << i << " along dimension #"
                 << j << " is " << grids_(i)->GetLength(j) << " but should be "
                 << grids_(j)->GetLength() << "." << endl;
            abort();
          }
#endif

  }

  //! Constructor for ten-dimensional data.
  /*!
    \param G0 grid for dimension #0.
    \param G1 grid for dimension #1.
    \param G2 grid for dimension #2.
    \param G3 grid for dimension #3.
    \param G4 grid for dimension #4.
    \param G5 grid for dimension #5.
    \param G6 grid for dimension #6.
    \param G7 grid for dimension #7.
    \param G8 grid for dimension #8.
    \param G9 grid for dimension #9.
  */
  template<class T, int N, class TG>
  Data<T, N, TG>::Data(Grid<TG>& G0, Grid<TG>& G1, Grid<TG>& G2,
                       Grid<TG>& G3, Grid<TG>& G4, Grid<TG>& G5,
                       Grid<TG>& G6, Grid<TG>& G7, Grid<TG>& G8,
                       Grid<TG>& G9)  throw():
    data_(G0.GetLength(), G1.GetLength(),
          G2.GetLength(), G3.GetLength(), G4.GetLength(),
          G5.GetLength(), G6.GetLength(), G7.GetLength(),
          G8.GetLength(), G9.GetLength()), grids_(N)
  {

#ifdef SELDONDATA_DEBUG_CHECK_DIMENSIONS
    if (N != 10)
      {
        cout << "ERROR!" << endl;
        cout << "Wrong dimension in Data<T, N, TG>::Data(Grid<TG>, ...)."
             << endl << "   Required " + to_str(N)
          + " grids, but got 10 grids." << endl;
        abort();
      }
#endif

    grids_(0) = G0.Copy();
    grids_(1) = G1.Copy();
    grids_(2) = G2.Copy();
    grids_(3) = G3.Copy();
    grids_(4) = G4.Copy();
    grids_(5) = G5.Copy();
    grids_(6) = G6.Copy();
    grids_(7) = G7.Copy();
    grids_(8) = G8.Copy();
    grids_(9) = G9.Copy();
    SetVariables();
    InitData();


#ifdef SELDONDATA_DEBUG_CHECK_DIMENSIONS
    int i, j;
    for (i = 0; i < N; i++)
      for (j = 0; j < N; j++)
        if ((grids_(i)->GetLength(j) != 0)
            && (grids_(i)->GetLength(j) != grids_(j)->GetLength()))
          {
            cout << "ERROR!" << endl;
            cout << "Wrong dimension in Data<T, N, TG>::Data(Grid<TG>, ...)."
                 << endl << "   Length of grid #" << i << " along dimension #"
                 << j << " is " << grids_(i)->GetLength(j) << " but should be "
                 << grids_(j)->GetLength() << "." << endl;
            abort();
          }
#endif

  }

  //! Constructor for one-dimensional data.
  /*!
    \param N0 grid length along dimension #0.
  */
  template<class T, int N, class TG>
  Data<T, N, TG>::Data(int N0)  throw():
    data_(N0), grids_(N)
  {

#ifdef SELDONDATA_DEBUG_CHECK_DIMENSIONS
    if (N != 1)
      {
        cout << "ERROR!" << endl;
        cout << "Wrong dimension in Data<T, N, TG>::Data(int, ...)."
             << endl << "   Required " + to_str(N)
          + " grids, but got 1 grid." << endl;
        abort();
      }
#endif

    RegularGrid<TG> G0(N0);

    grids_(0) = G0.Copy();
    SetVariables();
    InitData();

#ifdef SELDONDATA_DEBUG_CHECK_DIMENSIONS
    int i, j;
    for (i = 0; i < N; i++)
      for (j = 0; j < N; j++)
        if ((grids_(i)->GetLength(j) != 0)
            && (grids_(i)->GetLength(j) != grids_(j)->GetLength()))
          {
            cout << "ERROR!" << endl;
            cout << "Wrong dimension in Data<T, N, TG>::Data(int, ...)."
                 << endl << "   Length of grid #" << i << " along dimension #"
                 << j << " is " << grids_(i)->GetLength(j) << " but should be "
                 << grids_(j)->GetLength() << "." << endl;
            abort();
          }
#endif

  }

  //! Constructor for two-dimensional data.
  /*!
    \param N0 grid length along dimension #0.
    \param N1 grid length along dimension #1.
  */
  template<class T, int N, class TG>
  Data<T, N, TG>::Data(int N0, int N1)  throw():
    data_(N0, N1), grids_(N)
  {

#ifdef SELDONDATA_DEBUG_CHECK_DIMENSIONS
    if (N != 2)
      {
        cout << "ERROR!" << endl;
        cout << "Wrong dimension in Data<T, N, TG>::Data(int, ...)."
             << endl << "   Required " + to_str(N)
          + " grids, but got 2 grids." << endl;
        abort();
      }
#endif

    RegularGrid<TG> G0(N0);
    RegularGrid<TG> G1(N1);

    grids_(0) = G0.Copy();
    grids_(1) = G1.Copy();
    SetVariables();
    InitData();

#ifdef SELDONDATA_DEBUG_CHECK_DIMENSIONS
    int i, j;
    for (i = 0; i < N; i++)
      for (j = 0; j < N; j++)
        if ((grids_(i)->GetLength(j) != 0)
            && (grids_(i)->GetLength(j) != grids_(j)->GetLength()))
          {
            cout << "ERROR!" << endl;
            cout << "Wrong dimension in Data<T, N, TG>::Data(int, ...)."
                 << endl << "   Length of grid #" << i << " along dimension #"
                 << j << " is " << grids_(i)->GetLength(j) << " but should be "
                 << grids_(j)->GetLength() << "." << endl;
            abort();
          }
#endif

  }

  //! Constructor for three-dimensional data.
  /*!
    \param N0 grid length along dimension #0.
    \param N1 grid length along dimension #1.
    \param N2 grid length along dimension #2.
  */
  template<class T, int N, class TG>
  Data<T, N, TG>::Data(int N0, int N1, int N2)  throw():
    data_(N0, N1, N2), grids_(N)
  {

#ifdef SELDONDATA_DEBUG_CHECK_DIMENSIONS
    if (N != 3)
      {
        cout << "ERROR!" << endl;
        cout << "Wrong dimension in Data<T, N, TG>::Data(int, ...)."
             << endl << "   Required " + to_str(N)
          + " grids, but got 3 grids." << endl;
        abort();
      }
#endif

    RegularGrid<TG> G0(N0);
    RegularGrid<TG> G1(N1);
    RegularGrid<TG> G2(N2);

    grids_(0) = G0.Copy();
    grids_(1) = G1.Copy();
    grids_(2) = G2.Copy();
    SetVariables();
    InitData();


#ifdef SELDONDATA_DEBUG_CHECK_DIMENSIONS
    int i, j;
    for (i = 0; i < N; i++)
      for (j = 0; j < N; j++)
        if ((grids_(i)->GetLength(j) != 0)
            && (grids_(i)->GetLength(j) != grids_(j)->GetLength()))
          {
            cout << "ERROR!" << endl;
            cout << "Wrong dimension in Data<T, N, TG>::Data(int, ...)."
                 << endl << "   Length of grid #" << i << " along dimension #"
                 << j << " is " << grids_(i)->GetLength(j) << " but should be "
                 << grids_(j)->GetLength() << "." << endl;
            abort();
          }
#endif

  }

  //! Constructor for four-dimensional data.
  /*!
    \param N0 grid length along dimension #0.
    \param N1 grid length along dimension #1.
    \param N2 grid length along dimension #2.
    \param N3 grid length along dimension #3.
  */
  template<class T, int N, class TG>
  Data<T, N, TG>::Data(int N0, int N1, int N2,
                       int N3)  throw():
    data_(N0, N1, N2, N3), grids_(N)
  {

#ifdef SELDONDATA_DEBUG_CHECK_DIMENSIONS
    if (N != 4)
      {
        cout << "ERROR!" << endl;
        cout << "Wrong dimension in Data<T, N, TG>::Data(int, ...)."
             << endl << "   Required " + to_str(N)
          + " grids, but got 4 grids." << endl;
        abort();
      }
#endif

    RegularGrid<TG> G0(N0);
    RegularGrid<TG> G1(N1);
    RegularGrid<TG> G2(N2);
    RegularGrid<TG> G3(N3);

    grids_(0) = G0.Copy();
    grids_(1) = G1.Copy();
    grids_(2) = G2.Copy();
    grids_(3) = G3.Copy();
    SetVariables();
    InitData();


#ifdef SELDONDATA_DEBUG_CHECK_DIMENSIONS
    int i, j;
    for (i = 0; i < N; i++)
      for (j = 0; j < N; j++)
        if ((grids_(i)->GetLength(j) != 0)
            && (grids_(i)->GetLength(j) != grids_(j)->GetLength()))
          {
            cout << "ERROR!" << endl;
            cout << "Wrong dimension in Data<T, N, TG>::Data(int, ...)."
                 << endl << "   Length of grid #" << i << " along dimension #"
                 << j << " is " << grids_(i)->GetLength(j) << " but should be "
                 << grids_(j)->GetLength() << "." << endl;
            abort();
          }
#endif

  }

  //! Constructor for five-dimensional data.
  /*!
    \param N0 grid length along dimension #0.
    \param N1 grid length along dimension #1.
    \param N2 grid length along dimension #2.
    \param N3 grid length along dimension #3.
    \param N4 grid length along dimension #4.
  */
  template<class T, int N, class TG>
  Data<T, N, TG>::Data(int N0, int N1, int N2,
                       int N3, int N4)  throw():
    data_(N0, N1, N2, N3, N4), grids_(N)
  {

#ifdef SELDONDATA_DEBUG_CHECK_DIMENSIONS
    if (N != 5)
      {
        cout << "ERROR!" << endl;
        cout << "Wrong dimension in Data<T, N, TG>::Data(int, ...)."
             << endl << "   Required " + to_str(N)
          + " grids, but got 5 grids." << endl;
        abort();
      }
#endif

    RegularGrid<TG> G0(N0);
    RegularGrid<TG> G1(N1);
    RegularGrid<TG> G2(N2);
    RegularGrid<TG> G3(N3);
    RegularGrid<TG> G4(N4);

    grids_(0) = G0.Copy();
    grids_(1) = G1.Copy();
    grids_(2) = G2.Copy();
    grids_(3) = G3.Copy();
    grids_(4) = G4.Copy();
    SetVariables();
    InitData();


#ifdef SELDONDATA_DEBUG_CHECK_DIMENSIONS
    int i, j;
    for (i = 0; i < N; i++)
      for (j = 0; j < N; j++)
        if ((grids_(i)->GetLength(j) != 0)
            && (grids_(i)->GetLength(j) != grids_(j)->GetLength()))
          {
            cout << "ERROR!" << endl;
            cout << "Wrong dimension in Data<T, N, TG>::Data(int, ...)."
                 << endl << "   Length of grid #" << i << " along dimension #"
                 << j << " is " << grids_(i)->GetLength(j) << " but should be "
                 << grids_(j)->GetLength() << "." << endl;
            abort();
          }
#endif

  }

  //! Constructor for six-dimensional data.
  /*!
    \param N0 grid length along dimension #0.
    \param N1 grid length along dimension #1.
    \param N2 grid length along dimension #2.
    \param N3 grid length along dimension #3.
    \param N4 grid length along dimension #4.
    \param N5 grid length along dimension #5.
  */
  template<class T, int N, class TG>
  Data<T, N, TG>::Data(int N0, int N1, int N2,
                       int N3, int N4, int N5)  throw():
    data_(N0, N1, N2, N3, N4, N5), grids_(N)
  {

#ifdef SELDONDATA_DEBUG_CHECK_DIMENSIONS
    if (N != 6)
      {
        cout << "ERROR!" << endl;
        cout << "Wrong dimension in Data<T, N, TG>::Data(int, ...)."
             << endl << "   Required " + to_str(N)
          + " grids, but got 6 grids." << endl;
        abort();
      }
#endif

    RegularGrid<TG> G0(N0);
    RegularGrid<TG> G1(N1);
    RegularGrid<TG> G2(N2);
    RegularGrid<TG> G3(N3);
    RegularGrid<TG> G4(N4);
    RegularGrid<TG> G5(N5);

    grids_(0) = G0.Copy();
    grids_(1) = G1.Copy();
    grids_(2) = G2.Copy();
    grids_(3) = G3.Copy();
    grids_(4) = G4.Copy();
    grids_(5) = G5.Copy();
    SetVariables();
    InitData();


#ifdef SELDONDATA_DEBUG_CHECK_DIMENSIONS
    int i, j;
    for (i = 0; i < N; i++)
      for (j = 0; j < N; j++)
        if ((grids_(i)->GetLength(j) != 0)
            && (grids_(i)->GetLength(j) != grids_(j)->GetLength()))
          {
            cout << "ERROR!" << endl;
            cout << "Wrong dimension in Data<T, N, TG>::Data(int, ...)."
                 << endl << "   Length of grid #" << i << " along dimension #"
                 << j << " is " << grids_(i)->GetLength(j) << " but should be "
                 << grids_(j)->GetLength() << "." << endl;
            abort();
          }
#endif

  }

  //! Constructor for seven-dimensional data.
  /*!
    \param N0 grid length along dimension #0.
    \param N1 grid length along dimension #1.
    \param N2 grid length along dimension #2.
    \param N3 grid length along dimension #3.
    \param N4 grid length along dimension #4.
    \param N5 grid length along dimension #5.
    \param N6 grid length along dimension #6.
  */
  template<class T, int N, class TG>
  Data<T, N, TG>::Data(int N0, int N1, int N2,
                       int N3, int N4, int N5,
                       int N6)  throw():
    data_(N0, N1, N2, N3, N4, N5, N6), grids_(N)
  {

#ifdef SELDONDATA_DEBUG_CHECK_DIMENSIONS
    if (N != 7)
      {
        cout << "ERROR!" << endl;
        cout << "Wrong dimension in Data<T, N, TG>::Data(int, ...)."
             << endl << "   Required " + to_str(N)
          + " grids, but got 7 grids." << endl;
        abort();
      }
#endif

    RegularGrid<TG> G0(N0);
    RegularGrid<TG> G1(N1);
    RegularGrid<TG> G2(N2);
    RegularGrid<TG> G3(N3);
    RegularGrid<TG> G4(N4);
    RegularGrid<TG> G5(N5);
    RegularGrid<TG> G6(N6);

    grids_(0) = G0.Copy();
    grids_(1) = G1.Copy();
    grids_(2) = G2.Copy();
    grids_(3) = G3.Copy();
    grids_(4) = G4.Copy();
    grids_(5) = G5.Copy();
    grids_(6) = G6.Copy();
    SetVariables();
    InitData();


#ifdef SELDONDATA_DEBUG_CHECK_DIMENSIONS
    int i, j;
    for (i = 0; i < N; i++)
      for (j = 0; j < N; j++)
        if ((grids_(i)->GetLength(j) != 0)
            && (grids_(i)->GetLength(j) != grids_(j)->GetLength()))
          {
            cout << "ERROR!" << endl;
            cout << "Wrong dimension in Data<T, N, TG>::Data(int, ...)."
                 << endl << "   Length of grid #" << i << " along dimension #"
                 << j << " is " << grids_(i)->GetLength(j) << " but should be "
                 << grids_(j)->GetLength() << "." << endl;
            abort();
          }
#endif

  }

  //! Constructor for eight-dimensional data.
  /*!
    \param N0 grid length along dimension #0.
    \param N1 grid length along dimension #1.
    \param N2 grid length along dimension #2.
    \param N3 grid length along dimension #3.
    \param N4 grid length along dimension #4.
    \param N5 grid length along dimension #5.
    \param N6 grid length along dimension #6.
    \param N7 grid length along dimension #7.
  */
  template<class T, int N, class TG>
  Data<T, N, TG>::Data(int N0, int N1, int N2,
                       int N3, int N4, int N5,
                       int N6, int N7)  throw():
    data_(N0, N1, N2, N3, N4, N5, N6, N7), grids_(N)
  {

#ifdef SELDONDATA_DEBUG_CHECK_DIMENSIONS
    if (N != 8)
      {
        cout << "ERROR!" << endl;
        cout << "Wrong dimension in Data<T, N, TG>::Data(int, ...)."
             << endl << "   Required " + to_str(N)
          + " grids, but got 8 grids." << endl;
        abort();
      }
#endif

    RegularGrid<TG> G0(N0);
    RegularGrid<TG> G1(N1);
    RegularGrid<TG> G2(N2);
    RegularGrid<TG> G3(N3);
    RegularGrid<TG> G4(N4);
    RegularGrid<TG> G5(N5);
    RegularGrid<TG> G6(N6);
    RegularGrid<TG> G7(N7);

    grids_(0) = G0.Copy();
    grids_(1) = G1.Copy();
    grids_(2) = G2.Copy();
    grids_(3) = G3.Copy();
    grids_(4) = G4.Copy();
    grids_(5) = G5.Copy();
    grids_(6) = G6.Copy();
    grids_(7) = G7.Copy();
    SetVariables();
    InitData();


#ifdef SELDONDATA_DEBUG_CHECK_DIMENSIONS
    int i, j;
    for (i = 0; i < N; i++)
      for (j = 0; j < N; j++)
        if ((grids_(i)->GetLength(j) != 0)
            && (grids_(i)->GetLength(j) != grids_(j)->GetLength()))
          {
            cout << "ERROR!" << endl;
            cout << "Wrong dimension in Data<T, N, TG>::Data(int, ...)."
                 << endl << "   Length of grid #" << i << " along dimension #"
                 << j << " is " << grids_(i)->GetLength(j) << " but should be "
                 << grids_(j)->GetLength() << "." << endl;
            abort();
          }
#endif

  }

  //! Constructor for nine-dimensional data.
  /*!
    \param N0 grid length along dimension #0.
    \param N1 grid length along dimension #1.
    \param N2 grid length along dimension #2.
    \param N3 grid length along dimension #3.
    \param N4 grid length along dimension #4.
    \param N5 grid length along dimension #5.
    \param N6 grid length along dimension #6.
    \param N7 grid length along dimension #7.
    \param N8 grid length along dimension #8.
  */
  template<class T, int N, class TG>
  Data<T, N, TG>::Data(int N0, int N1, int N2,
                       int N3, int N4, int N5,
                       int N6, int N7, int N8)  throw():
    data_(N0, N1, N2, N3, N4, N5, N6, N7, N8), grids_(N)
  {

#ifdef SELDONDATA_DEBUG_CHECK_DIMENSIONS
    if (N != 9)
      {
        cout << "ERROR!" << endl;
        cout << "Wrong dimension in Data<T, N, TG>::Data(int, ...)."
             << endl << "   Required " + to_str(N)
          + " grids, but got 9 grids." << endl;
        abort();
      }
#endif

    RegularGrid<TG> G0(N0);
    RegularGrid<TG> G1(N1);
    RegularGrid<TG> G2(N2);
    RegularGrid<TG> G3(N3);
    RegularGrid<TG> G4(N4);
    RegularGrid<TG> G5(N5);
    RegularGrid<TG> G6(N6);
    RegularGrid<TG> G7(N7);
    RegularGrid<TG> G8(N8);

    grids_(0) = G0.Copy();
    grids_(1) = G1.Copy();
    grids_(2) = G2.Copy();
    grids_(3) = G3.Copy();
    grids_(4) = G4.Copy();
    grids_(5) = G5.Copy();
    grids_(6) = G6.Copy();
    grids_(7) = G7.Copy();
    grids_(8) = G8.Copy();
    SetVariables();
    InitData();


#ifdef SELDONDATA_DEBUG_CHECK_DIMENSIONS
    int i, j;
    for (i = 0; i < N; i++)
      for (j = 0; j < N; j++)
        if ((grids_(i)->GetLength(j) != 0)
            && (grids_(i)->GetLength(j) != grids_(j)->GetLength()))
          {
            cout << "ERROR!" << endl;
            cout << "Wrong dimension in Data<T, N, TG>::Data(int, ...)."
                 << endl << "   Length of grid #" << i << " along dimension #"
                 << j << " is " << grids_(i)->GetLength(j) << " but should be "
                 << grids_(j)->GetLength() << "." << endl;
            abort();
          }
#endif

  }

  //! Constructor for ten-dimensional data.
  /*!
    \param N0 grid length along dimension #0.
    \param N1 grid length along dimension #1.
    \param N2 grid length along dimension #2.
    \param N3 grid length along dimension #3.
    \param N4 grid length along dimension #4.
    \param N5 grid length along dimension #5.
    \param N6 grid length along dimension #6.
    \param N7 grid length along dimension #7.
    \param N8 grid length along dimension #8.
    \param N9 grid length along dimension #9.
  */
  template<class T, int N, class TG>
  Data<T, N, TG>::Data(int N0, int N1, int N2,
                       int N3, int N4, int N5,
                       int N6, int N7, int N8,
                       int N9)  throw():
    data_(N0, N1, N2, N3, N4, N5, N6, N7, N8, N9), grids_(N)
  {

#ifdef SELDONDATA_DEBUG_CHECK_DIMENSIONS
    if (N != 10)
      {
        cout << "ERROR!" << endl;
        cout << "Wrong dimension in Data<T, N, TG>::Data(int, ...)."
             << endl << "   Required " + to_str(N)
          + " grids, but got 10 grids." << endl;
        abort();
      }
#endif

    RegularGrid<TG> G0(N0);
    RegularGrid<TG> G1(N1);
    RegularGrid<TG> G2(N2);
    RegularGrid<TG> G3(N3);
    RegularGrid<TG> G4(N4);
    RegularGrid<TG> G5(N5);
    RegularGrid<TG> G6(N6);
    RegularGrid<TG> G7(N7);
    RegularGrid<TG> G8(N8);
    RegularGrid<TG> G9(N9);

    grids_(0) = G0.Copy();
    grids_(1) = G1.Copy();
    grids_(2) = G2.Copy();
    grids_(3) = G3.Copy();
    grids_(4) = G4.Copy();
    grids_(5) = G5.Copy();
    grids_(6) = G6.Copy();
    grids_(7) = G7.Copy();
    grids_(8) = G8.Copy();
    grids_(9) = G9.Copy();
    SetVariables();
    InitData();


#ifdef SELDONDATA_DEBUG_CHECK_DIMENSIONS
    int i, j;
    for (i = 0; i < N; i++)
      for (j = 0; j < N; j++)
        if ((grids_(i)->GetLength(j) != 0)
            && (grids_(i)->GetLength(j) != grids_(j)->GetLength()))
          {
            cout << "ERROR!" << endl;
            cout << "Wrong dimension in Data<T, N, TG>::Data(int, ...)."
                 << endl << "   Length of grid #" << i << " along dimension #"
                 << j << " is " << grids_(i)->GetLength(j) << " but should be "
                 << grids_(j)->GetLength() << "." << endl;
            abort();
          }
#endif

  }

  //! Constructor.
  template<class T, int N, class TG>
  Data<T, N, TG>::Data(const TinyVector<int, N>& shape)  throw():
    grids_(N)
  {
    for (int i = 0; i < N; i++)
      grids_(i) = NULL;

    if (N == 1)
      this->Resize(shape(0));
    else if (N == 2)
      this->Resize(shape(0), shape(1));
    else if (N == 3)
      this->Resize(shape(0), shape(1), shape(2));
    else if (N == 4)
      this->Resize(shape(0), shape(1), shape(2),
                   shape(3));
    else if (N == 5)
      this->Resize(shape(0), shape(1), shape(2),
                   shape(3), shape(4));
    else if (N == 6)
      this->Resize(shape(0), shape(1), shape(2),
                   shape(3), shape(4), shape(5));
    else if (N == 7)
      this->Resize(shape(0), shape(1), shape(2),
                   shape(3), shape(4), shape(5),
                   shape(6));
    else if (N == 8)
      this->Resize(shape(0), shape(1), shape(2),
                   shape(3), shape(4), shape(5),
                   shape(6), shape(7));
    else if (N == 9)
      this->Resize(shape(0), shape(1), shape(2),
                   shape(3), shape(4), shape(5),
                   shape(6), shape(7), shape(8));
    else if (N == 10)
      this->Resize(shape(0), shape(1), shape(2),
                   shape(3), shape(4), shape(5),
                   shape(6), shape(7), shape(8),
                   shape(9));
  }

  //! Constructor.
  /*!
    Builds a Data object from preexisting data.
  */
  template<class T, int N, class TG>
  Data<T, N, TG>::Data(T* data, const TinyVector<int, N>& shape,
                       preexistingMemoryPolicy policy)  throw():
    data_(data, shape, policy), grids_(N)
  {
    for (int i = 0; i < N; i++)
      {
        RegularGrid<TG> G(shape(i));
        grids_(i) = G.Copy();
      }
    SetVariables();
  }

  //! Copy constructor.
  template<class T, int N, class TG>
  template <class T0>
  Data<T, N, TG>::Data(Data<T0, N, TG>& data)  throw():
    grids_(N)
  {
    for (int i = 0; i < N; i++)
      if (data.GetGrid(i) != NULL)
        grids_(i) = data[i].Copy();
      else
        grids_(i) = NULL;

    data_.resize(data.GetArray().shape());
    data_ = data.GetArray();
  }

  //! Destructor.
  template<class T, int N, class TG>
  Data<T, N, TG>::~Data()  throw()
  {
    for (int i = 0; i < N; i++)
      if (grids_(i) != NULL)
        if (grids_(i)->GetPointers() == 1)
          delete grids_(i);
        else
          grids_(i)->SetPointers(grids_(i)->GetPointers() - 1);
  }

  //! Access operator for one-dimensional data.
  /*!
    \param i0 index for dimension #0.
    \return A reference to element #i0.
    \exception SeldonData::WrongDim data dimension is not 1.
    \exception SeldonData::WrongIndex index is out of range.
  */
  template<class T, int N, class TG>
  inline T& Data<T, N, TG>::operator()(int i0)
  {

#ifdef SELDONDATA_DEBUG_CHECK_DIMENSIONS
    if (N != 1)
      throw WrongDim("Data<T, " + to_str(N) + ">::operator()",
                     "operator() was called with 1 argument.");
#endif

#ifdef SELDONDATA_DEBUG_CHECK_INDICES
    if ((i0 < 0) || (i0 >= data_.extent(0)))
      throw WrongIndex("Data<T, " + to_str(N) + ">::operator()",
                       "Index along dimension #0 should be in [0, "
                       + to_str(data_.extent(0) - 1) + "], but is equal to "
                       + to_str(i0) + ".");
#endif

    return data_(i0);
  }

  //! Access operator for two-dimensional data.
  /*!
    \param i0 index for dimension #0.
    \param i1 index for dimension #1.
    \return A reference to element (i0, i1).
    \exception SeldonData::WrongDim data dimension is not 2.
    \exception SeldonData::WrongIndex an index is out of range.
  */
  template<class T, int N, class TG>
  inline T& Data<T, N, TG>::operator()(int i0, int i1)
  {

#ifdef SELDONDATA_DEBUG_CHECK_DIMENSIONS
    if (N != 2)
      throw WrongDim("Data<T, " + to_str(N) + ">::operator()",
                     "operator() was called with 2 arguments.");
#endif

#ifdef SELDONDATA_DEBUG_CHECK_INDICES
    if ((i0 < 0) || (i0 >= data_.extent(0)))
      throw WrongIndex("Data<T, " + to_str(N) + ">::operator()",
                       "Index along dimension #0 should be in [0, "
                       + to_str(data_.extent(0) - 1) + "], but is equal to "
                       + to_str(i0) + ".");
    if ((i1 < 0) || (i1 >= data_.extent(1)))
      throw WrongIndex("Data<T, " + to_str(N) + ">::operator()",
                       "Index along dimension #1 should be in [0, "
                       + to_str(data_.extent(1) - 1) + "], but is equal to "
                       + to_str(i1) + ".");
#endif

    return data_(i0, i1);
  }

  //! Access operator for three-dimensional data.
  /*!
    \param i0 index for dimension #0.
    \param i1 index for dimension #1.
    \param i2 index for dimension #2.
    \return A reference to element (i0, i1, i2).
    \exception SeldonData::WrongDim data dimension is not 3.
    \exception SeldonData::WrongIndex an index is out of range.
  */
  template<class T, int N, class TG>
  inline T& Data<T, N, TG>::operator()(int i0, int i1, int i2)
  {

#ifdef SELDONDATA_DEBUG_CHECK_DIMENSIONS
    if (N != 3)
      throw WrongDim("Data<T, " + to_str(N) + ">::operator()",
                     "operator() was called with 3 arguments.");
#endif

#ifdef SELDONDATA_DEBUG_CHECK_INDICES
    if ((i0 < 0) || (i0 >= data_.extent(0)))
      throw WrongIndex("Data<T, " + to_str(N) + ">::operator()",
                       "Index along dimension #0 should be in [0, "
                       + to_str(data_.extent(0) - 1) + "], but is equal to "
                       + to_str(i0) + ".");
    if ((i1 < 0) || (i1 >= data_.extent(1)))
      throw WrongIndex("Data<T, " + to_str(N) + ">::operator()",
                       "Index along dimension #1 should be in [0, "
                       + to_str(data_.extent(1) - 1) + "], but is equal to "
                       + to_str(i1) + ".");
    if ((i2 < 0) || (i2 >= data_.extent(2)))
      throw WrongIndex("Data<T, " + to_str(N) + ">::operator()",
                       "Index along dimension #2 should be in [0, "
                       + to_str(data_.extent(2) - 1) + "], but is equal to "
                       + to_str(i2) + ".");
#endif

    return data_(i0, i1, i2);
  }

  //! Access operator for four-dimensional data.
  /*!
    \param i0 index for dimension #0.
    \param i1 index for dimension #1.
    \param i2 index for dimension #2.
    \param i3 index for dimension #3.
    \return A reference to element (i0, i1, i2, i3).
    \exception SeldonData::WrongDim data dimension is not 4.
    \exception SeldonData::WrongIndex an index is out of range.
  */
  template<class T, int N, class TG>
  inline T& Data<T, N, TG>::operator()(int i0, int i1, int i2,
                                       int i3)

  {

#ifdef SELDONDATA_DEBUG_CHECK_DIMENSIONS
    if (N != 4)
      throw WrongDim("Data<T, " + to_str(N) + ">::operator()",
                     "operator() was called with 4 arguments.");
#endif

#ifdef SELDONDATA_DEBUG_CHECK_INDICES
    if ((i0 < 0) || (i0 >= data_.extent(0)))
      throw WrongIndex("Data<T, " + to_str(N) + ">::operator()",
                       "Index along dimension #0 should be in [0, "
                       + to_str(data_.extent(0) - 1) + "], but is equal to "
                       + to_str(i0) + ".");
    if ((i1 < 0) || (i1 >= data_.extent(1)))
      throw WrongIndex("Data<T, " + to_str(N) + ">::operator()",
                       "Index along dimension #1 should be in [0, "
                       + to_str(data_.extent(1) - 1) + "], but is equal to "
                       + to_str(i1) + ".");
    if ((i2 < 0) || (i2 >= data_.extent(2)))
      throw WrongIndex("Data<T, " + to_str(N) + ">::operator()",
                       "Index along dimension #2 should be in [0, "
                       + to_str(data_.extent(2) - 1) + "], but is equal to "
                       + to_str(i2) + ".");
    if ((i3 < 0) || (i3 >= data_.extent(3)))
      throw WrongIndex("Data<T, " + to_str(N) + ">::operator()",
                       "Index along dimension #3 should be in [0, "
                       + to_str(data_.extent(3) - 1) + "], but is equal to "
                       + to_str(i3) + ".");
#endif

    return data_(i0, i1, i2, i3);
  }

  //! Access operator for five-dimensional data.
  /*!
    \param i0 index for dimension #0.
    \param i1 index for dimension #1.
    \param i2 index for dimension #2.
    \param i3 index for dimension #3.
    \param i4 index for dimension #4.
    \return A reference to element (i0, i1, i2, i3, i4).
    \exception SeldonData::WrongDim data dimension is not 5.
    \exception SeldonData::WrongIndex an index is out of range.
  */
  template<class T, int N, class TG>
  inline T& Data<T, N, TG>::operator()(int i0, int i1, int i2,
                                       int i3, int i4)

  {

#ifdef SELDONDATA_DEBUG_CHECK_DIMENSIONS
    if (N != 5)
      throw WrongDim("Data<T, " + to_str(N) + ">::operator()",
                     "operator() was called with 5 arguments.");
#endif

#ifdef SELDONDATA_DEBUG_CHECK_INDICES
    if ((i0 < 0) || (i0 >= data_.extent(0)))
      throw WrongIndex("Data<T, " + to_str(N) + ">::operator()",
                       "Index along dimension #0 should be in [0, "
                       + to_str(data_.extent(0) - 1) + "], but is equal to "
                       + to_str(i0) + ".");
    if ((i1 < 0) || (i1 >= data_.extent(1)))
      throw WrongIndex("Data<T, " + to_str(N) + ">::operator()",
                       "Index along dimension #1 should be in [0, "
                       + to_str(data_.extent(1) - 1) + "], but is equal to "
                       + to_str(i1) + ".");
    if ((i2 < 0) || (i2 >= data_.extent(2)))
      throw WrongIndex("Data<T, " + to_str(N) + ">::operator()",
                       "Index along dimension #2 should be in [0, "
                       + to_str(data_.extent(2) - 1) + "], but is equal to "
                       + to_str(i2) + ".");
    if ((i3 < 0) || (i3 >= data_.extent(3)))
      throw WrongIndex("Data<T, " + to_str(N) + ">::operator()",
                       "Index along dimension #3 should be in [0, "
                       + to_str(data_.extent(3) - 1) + "], but is equal to "
                       + to_str(i3) + ".");
    if ((i4 < 0) || (i4 >= data_.extent(4)))
      throw WrongIndex("Data<T, " + to_str(N) + ">::operator()",
                       "Index along dimension #4 should be in [0, "
                       + to_str(data_.extent(4) - 1) + "], but is equal to "
                       + to_str(i4) + ".");
#endif

    return data_(i0, i1, i2, i3, i4);
  }

  //! Access operator for six-dimensional data.
  /*!
    \param i0 index for dimension #0.
    \param i1 index for dimension #1.
    \param i2 index for dimension #2.
    \param i3 index for dimension #3.
    \param i4 index for dimension #4.
    \param i5 index for dimension #5.
    \return A reference to element (i0, i1, i2, i3, i4, i5).
    \exception SeldonData::WrongDim data dimension is not 6.
    \exception SeldonData::WrongIndex an index is out of range.
  */
  template<class T, int N, class TG>
  inline T& Data<T, N, TG>::operator()(int i0, int i1, int i2,
                                       int i3, int i4, int i5)

  {

#ifdef SELDONDATA_DEBUG_CHECK_DIMENSIONS
    if (N != 6)
      throw WrongDim("Data<T, " + to_str(N) + ">::operator()",
                     "operator() was called with 6 arguments.");
#endif

#ifdef SELDONDATA_DEBUG_CHECK_INDICES
    if ((i0 < 0) || (i0 >= data_.extent(0)))
      throw WrongIndex("Data<T, " + to_str(N) + ">::operator()",
                       "Index along dimension #0 should be in [0, "
                       + to_str(data_.extent(0) - 1) + "], but is equal to "
                       + to_str(i0) + ".");
    if ((i1 < 0) || (i1 >= data_.extent(1)))
      throw WrongIndex("Data<T, " + to_str(N) + ">::operator()",
                       "Index along dimension #1 should be in [0, "
                       + to_str(data_.extent(1) - 1) + "], but is equal to "
                       + to_str(i1) + ".");
    if ((i2 < 0) || (i2 >= data_.extent(2)))
      throw WrongIndex("Data<T, " + to_str(N) + ">::operator()",
                       "Index along dimension #2 should be in [0, "
                       + to_str(data_.extent(2) - 1) + "], but is equal to "
                       + to_str(i2) + ".");
    if ((i3 < 0) || (i3 >= data_.extent(3)))
      throw WrongIndex("Data<T, " + to_str(N) + ">::operator()",
                       "Index along dimension #3 should be in [0, "
                       + to_str(data_.extent(3) - 1) + "], but is equal to "
                       + to_str(i3) + ".");
    if ((i4 < 0) || (i4 >= data_.extent(4)))
      throw WrongIndex("Data<T, " + to_str(N) + ">::operator()",
                       "Index along dimension #4 should be in [0, "
                       + to_str(data_.extent(4) - 1) + "], but is equal to "
                       + to_str(i4) + ".");
    if ((i5 < 0) || (i5 >= data_.extent(5)))
      throw WrongIndex("Data<T, " + to_str(N) + ">::operator()",
                       "Index along dimension #5 should be in [0, "
                       + to_str(data_.extent(5) - 1) + "], but is equal to "
                       + to_str(i5) + ".");
#endif

    return data_(i0, i1, i2, i3, i4, i5);
  }

  template<class T, int N, class TG>
  //! Access operator for seven-dimensional data.
  /*!
    \param i0 index for dimension #0.
    \param i1 index for dimension #1.
    \param i2 index for dimension #2.
    \param i3 index for dimension #3.
    \param i4 index for dimension #4.
    \param i5 index for dimension #5.
    \param i6 index for dimension #6.
    \return A reference to element (i0, i1, i2, i3,
    i4, i5, i6).
    \exception SeldonData::WrongDim data dimension is not 7.
    \exception SeldonData::WrongIndex an index is out of range.
  */
  inline T& Data<T, N, TG>::operator()(int i0, int i1, int i2,
                                       int i3, int i4, int i5,
                                       int i6)

  {

#ifdef SELDONDATA_DEBUG_CHECK_DIMENSIONS
    if (N != 7)
      throw WrongDim("Data<T, " + to_str(N) + ">::operator()",
                     "operator() was called with 7 arguments.");
#endif

#ifdef SELDONDATA_DEBUG_CHECK_INDICES
    if ((i0 < 0) || (i0 >= data_.extent(0)))
      throw WrongIndex("Data<T, " + to_str(N) + ">::operator()",
                       "Index along dimension #0 should be in [0, "
                       + to_str(data_.extent(0) - 1) + "], but is equal to "
                       + to_str(i0) + ".");
    if ((i1 < 0) || (i1 >= data_.extent(1)))
      throw WrongIndex("Data<T, " + to_str(N) + ">::operator()",
                       "Index along dimension #1 should be in [0, "
                       + to_str(data_.extent(1) - 1) + "], but is equal to "
                       + to_str(i1) + ".");
    if ((i2 < 0) || (i2 >= data_.extent(2)))
      throw WrongIndex("Data<T, " + to_str(N) + ">::operator()",
                       "Index along dimension #2 should be in [0, "
                       + to_str(data_.extent(2) - 1) + "], but is equal to "
                       + to_str(i2) + ".");
    if ((i3 < 0) || (i3 >= data_.extent(3)))
      throw WrongIndex("Data<T, " + to_str(N) + ">::operator()",
                       "Index along dimension #3 should be in [0, "
                       + to_str(data_.extent(3) - 1) + "], but is equal to "
                       + to_str(i3) + ".");
    if ((i4 < 0) || (i4 >= data_.extent(4)))
      throw WrongIndex("Data<T, " + to_str(N) + ">::operator()",
                       "Index along dimension #4 should be in [0, "
                       + to_str(data_.extent(4) - 1) + "], but is equal to "
                       + to_str(i4) + ".");
    if ((i5 < 0) || (i5 >= data_.extent(5)))
      throw WrongIndex("Data<T, " + to_str(N) + ">::operator()",
                       "Index along dimension #5 should be in [0, "
                       + to_str(data_.extent(5) - 1) + "], but is equal to "
                       + to_str(i5) + ".");
    if ((i6 < 0) || (i6 >= data_.extent(6)))
      throw WrongIndex("Data<T, " + to_str(N) + ">::operator()",
                       "Index along dimension #6 should be in [0, "
                       + to_str(data_.extent(6) - 1) + "], but is equal to "
                       + to_str(i6) + ".");
#endif

    return data_(i0, i1, i2, i3, i4, i5, i6);
  }

  template<class T, int N, class TG>
  //! Access operator for eight-dimensional data.
  /*!
    \param i0 index for dimension #0.
    \param i1 index for dimension #1.
    \param i2 index for dimension #2.
    \param i3 index for dimension #3.
    \param i4 index for dimension #4.
    \param i5 index for dimension #5.
    \param i6 index for dimension #6.
    \param i7 index for dimension #7.
    \return A reference to element (i0, i1, i2, i3, i4,
    i5, i6, i7).
    \exception SeldonData::WrongDim data dimension is not 8.
    \exception SeldonData::WrongIndex an index is out of range.
  */
  inline T& Data<T, N, TG>::operator()(int i0, int i1, int i2,
                                       int i3, int i4, int i5,
                                       int i6, int i7)

  {

#ifdef SELDONDATA_DEBUG_CHECK_DIMENSIONS
    if (N != 8)
      throw WrongDim("Data<T, " + to_str(N) + ">::operator()",
                     "operator() was called with 8 arguments.");
#endif

#ifdef SELDONDATA_DEBUG_CHECK_INDICES
    if ((i0 < 0) || (i0 >= data_.extent(0)))
      throw WrongIndex("Data<T, " + to_str(N) + ">::operator()",
                       "Index along dimension #0 should be in [0, "
                       + to_str(data_.extent(0) - 1) + "], but is equal to "
                       + to_str(i0) + ".");
    if ((i1 < 0) || (i1 >= data_.extent(1)))
      throw WrongIndex("Data<T, " + to_str(N) + ">::operator()",
                       "Index along dimension #1 should be in [0, "
                       + to_str(data_.extent(1) - 1) + "], but is equal to "
                       + to_str(i1) + ".");
    if ((i2 < 0) || (i2 >= data_.extent(2)))
      throw WrongIndex("Data<T, " + to_str(N) + ">::operator()",
                       "Index along dimension #2 should be in [0, "
                       + to_str(data_.extent(2) - 1) + "], but is equal to "
                       + to_str(i2) + ".");
    if ((i3 < 0) || (i3 >= data_.extent(3)))
      throw WrongIndex("Data<T, " + to_str(N) + ">::operator()",
                       "Index along dimension #3 should be in [0, "
                       + to_str(data_.extent(3) - 1) + "], but is equal to "
                       + to_str(i3) + ".");
    if ((i4 < 0) || (i4 >= data_.extent(4)))
      throw WrongIndex("Data<T, " + to_str(N) + ">::operator()",
                       "Index along dimension #4 should be in [0, "
                       + to_str(data_.extent(4) - 1) + "], but is equal to "
                       + to_str(i4) + ".");
    if ((i5 < 0) || (i5 >= data_.extent(5)))
      throw WrongIndex("Data<T, " + to_str(N) + ">::operator()",
                       "Index along dimension #5 should be in [0, "
                       + to_str(data_.extent(5) - 1) + "], but is equal to "
                       + to_str(i5) + ".");
    if ((i6 < 0) || (i6 >= data_.extent(6)))
      throw WrongIndex("Data<T, " + to_str(N) + ">::operator()",
                       "Index along dimension #6 should be in [0, "
                       + to_str(data_.extent(6) - 1) + "], but is equal to "
                       + to_str(i6) + ".");
    if ((i7 < 0) || (i7 >= data_.extent(7)))
      throw WrongIndex("Data<T, " + to_str(N) + ">::operator()",
                       "Index along dimension #7 should be in [0, "
                       + to_str(data_.extent(7) - 1) + "], but is equal to "
                       + to_str(i7) + ".");
#endif

    return data_(i0, i1, i2, i3, i4, i5, i6, i7);
  }

  template<class T, int N, class TG>
  //! Access operator for nine-dimensional data.
  /*!
    \param i0 index for dimension #0.
    \param i1 index for dimension #1.
    \param i2 index for dimension #2.
    \param i3 index for dimension #3.
    \param i4 index for dimension #4.
    \param i5 index for dimension #5.
    \param i6 index for dimension #6.
    \param i7 index for dimension #7.
    \param i8 index for dimension #8.
    \return A reference to element (i0, i1, i2, i3, i4,
    i5, i6, i7, i8).
    \exception SeldonData::WrongDim data dimension is not 9.
    \exception SeldonData::WrongIndex an index is out of range.
  */
  inline T& Data<T, N, TG>::operator()(int i0, int i1, int i2,
                                       int i3, int i4, int i5,
                                       int i6, int i7, int i8)

  {

#ifdef SELDONDATA_DEBUG_CHECK_DIMENSIONS
    if (N != 9)
      throw WrongDim("Data<T, " + to_str(N) + ">::operator()",
                     "operator() was called with 9 arguments.");
#endif

#ifdef SELDONDATA_DEBUG_CHECK_INDICES
    if ((i0 < 0) || (i0 >= data_.extent(0)))
      throw WrongIndex("Data<T, " + to_str(N) + ">::operator()",
                       "Index along dimension #0 should be in [0, "
                       + to_str(data_.extent(0) - 1) + "], but is equal to "
                       + to_str(i0) + ".");
    if ((i1 < 0) || (i1 >= data_.extent(1)))
      throw WrongIndex("Data<T, " + to_str(N) + ">::operator()",
                       "Index along dimension #1 should be in [0, "
                       + to_str(data_.extent(1) - 1) + "], but is equal to "
                       + to_str(i1) + ".");
    if ((i2 < 0) || (i2 >= data_.extent(2)))
      throw WrongIndex("Data<T, " + to_str(N) + ">::operator()",
                       "Index along dimension #2 should be in [0, "
                       + to_str(data_.extent(2) - 1) + "], but is equal to "
                       + to_str(i2) + ".");
    if ((i3 < 0) || (i3 >= data_.extent(3)))
      throw WrongIndex("Data<T, " + to_str(N) + ">::operator()",
                       "Index along dimension #3 should be in [0, "
                       + to_str(data_.extent(3) - 1) + "], but is equal to "
                       + to_str(i3) + ".");
    if ((i4 < 0) || (i4 >= data_.extent(4)))
      throw WrongIndex("Data<T, " + to_str(N) + ">::operator()",
                       "Index along dimension #4 should be in [0, "
                       + to_str(data_.extent(4) - 1) + "], but is equal to "
                       + to_str(i4) + ".");
    if ((i5 < 0) || (i5 >= data_.extent(5)))
      throw WrongIndex("Data<T, " + to_str(N) + ">::operator()",
                       "Index along dimension #5 should be in [0, "
                       + to_str(data_.extent(5) - 1) + "], but is equal to "
                       + to_str(i5) + ".");
    if ((i6 < 0) || (i6 >= data_.extent(6)))
      throw WrongIndex("Data<T, " + to_str(N) + ">::operator()",
                       "Index along dimension #6 should be in [0, "
                       + to_str(data_.extent(6) - 1) + "], but is equal to "
                       + to_str(i6) + ".");
    if ((i7 < 0) || (i7 >= data_.extent(7)))
      throw WrongIndex("Data<T, " + to_str(N) + ">::operator()",
                       "Index along dimension #7 should be in [0, "
                       + to_str(data_.extent(7) - 1) + "], but is equal to "
                       + to_str(i7) + ".");
    if ((i8 < 0) || (i8 >= data_.extent(8)))
      throw WrongIndex("Data<T, " + to_str(N) + ">::operator()",
                       "Index along dimension #8 should be in [0, "
                       + to_str(data_.extent(8) - 1) + "], but is equal to "
                       + to_str(i8) + ".");
#endif

    return data_(i0, i1, i2, i3, i4, i5, i6, i7, i8);
  }

  //! Access operator for ten-dimensional data.
  /*!
    \param i0 index for dimension #0.
    \param i1 index for dimension #1.
    \param i2 index for dimension #2.
    \param i3 index for dimension #3.
    \param i4 index for dimension #4.
    \param i5 index for dimension #5.
    \param i6 index for dimension #6.
    \param i7 index for dimension #7.
    \param i8 index for dimension #8.
    \param i9 index for dimension #9.
    \return A reference to element (i0, i1, i2, i3, i4,
    i5, i6, i7, i8, i9).
    \exception SeldonData::WrongDim data dimension is not 10.
    \exception SeldonData::WrongIndex an index is out of range.
  */
  template<class T, int N, class TG>
  inline T& Data<T, N, TG>::operator()(int i0, int i1, int i2,
                                       int i3, int i4, int i5,
                                       int i6, int i7, int i8,
                                       int i9)

  {

#ifdef SELDONDATA_DEBUG_CHECK_DIMENSIONS
    if (N != 10)
      throw WrongDim("Data<T, " + to_str(N) + ">::operator()",
                     "operator() was called with 10 arguments.");
#endif

#ifdef SELDONDATA_DEBUG_CHECK_INDICES
    if ((i0 < 0) || (i0 >= data_.extent(0)))
      throw WrongIndex("Data<T, " + to_str(N) + ">::operator()",
                       "Index along dimension #0 should be in [0, "
                       + to_str(data_.extent(0) - 1) + "], but is equal to "
                       + to_str(i0) + ".");
    if ((i1 < 0) || (i1 >= data_.extent(1)))
      throw WrongIndex("Data<T, " + to_str(N) + ">::operator()",
                       "Index along dimension #1 should be in [0, "
                       + to_str(data_.extent(1) - 1) + "], but is equal to "
                       + to_str(i1) + ".");
    if ((i2 < 0) || (i2 >= data_.extent(2)))
      throw WrongIndex("Data<T, " + to_str(N) + ">::operator()",
                       "Index along dimension #2 should be in [0, "
                       + to_str(data_.extent(2) - 1) + "], but is equal to "
                       + to_str(i2) + ".");
    if ((i3 < 0) || (i3 >= data_.extent(3)))
      throw WrongIndex("Data<T, " + to_str(N) + ">::operator()",
                       "Index along dimension #3 should be in [0, "
                       + to_str(data_.extent(3) - 1) + "], but is equal to "
                       + to_str(i3) + ".");
    if ((i4 < 0) || (i4 >= data_.extent(4)))
      throw WrongIndex("Data<T, " + to_str(N) + ">::operator()",
                       "Index along dimension #4 should be in [0, "
                       + to_str(data_.extent(4) - 1) + "], but is equal to "
                       + to_str(i4) + ".");
    if ((i5 < 0) || (i5 >= data_.extent(5)))
      throw WrongIndex("Data<T, " + to_str(N) + ">::operator()",
                       "Index along dimension #5 should be in [0, "
                       + to_str(data_.extent(5) - 1) + "], but is equal to "
                       + to_str(i5) + ".");
    if ((i6 < 0) || (i6 >= data_.extent(6)))
      throw WrongIndex("Data<T, " + to_str(N) + ">::operator()",
                       "Index along dimension #6 should be in [0, "
                       + to_str(data_.extent(6) - 1) + "], but is equal to "
                       + to_str(i6) + ".");
    if ((i7 < 0) || (i7 >= data_.extent(7)))
      throw WrongIndex("Data<T, " + to_str(N) + ">::operator()",
                       "Index along dimension #7 should be in [0, "
                       + to_str(data_.extent(7) - 1) + "], but is equal to "
                       + to_str(i7) + ".");
    if ((i8 < 0) || (i8 >= data_.extent(8)))
      throw WrongIndex("Data<T, " + to_str(N) + ">::operator()",
                       "Index along dimension #8 should be in [0, "
                       + to_str(data_.extent(8) - 1) + "], but is equal to "
                       + to_str(i8) + ".");
    if ((i9 < 0) || (i9 >= data_.extent(9)))
      throw WrongIndex("Data<T, " + to_str(N) + ">::operator()",
                       "Index along dimension #9 should be in [0, "
                       + to_str(data_.extent(9) - 1) + "], but is equal to "
                       + to_str(i9) + ".");
#endif

    return data_(i0, i1, i2, i3, i4, i5, i6, i7, i8, i9);
  }

  //! Access operator for one-dimensional data.
  /*!
    \param i0 index for dimension #0.
    \return A const reference to element #i0.
    \exception SeldonData::WrongDim data dimension is not 1.
    \exception SeldonData::WrongIndex index is out of range.
  */
  template<class T, int N, class TG>
  inline T Data<T, N, TG>::operator()(int i0)  const
  {

#ifdef SELDONDATA_DEBUG_CHECK_DIMENSIONS
    if (N != 1)
      throw WrongDim("Data<T, " + to_str(N) + ">::operator() const",
                     "operator() was called with 1 argument.");
#endif

#ifdef SELDONDATA_DEBUG_CHECK_INDICES
    if ((i0 < 0) || (i0 >= data_.extent(0)))
      throw WrongIndex("Data<T, " + to_str(N) + ">::operator() const",
                       "Index along dimension #0 should be in [0, "
                       + to_str(data_.extent(0) - 1) + "], but is equal to "
                       + to_str(i0) + ".");
#endif

    return data_(i0);
  }

  //! Access operator for two-dimensional data.
  /*!
    \param i0 index for dimension #0.
    \param i1 index for dimension #1.
    \return A const reference to element (i0, i1).
    \exception SeldonData::WrongDim data dimension is not 2.
    \exception SeldonData::WrongIndex an index is out of range.
  */
  template<class T, int N, class TG>
  inline T Data<T, N, TG>::operator()(int i0, int i1)  const
  {

#ifdef SELDONDATA_DEBUG_CHECK_DIMENSIONS
    if (N != 2)
      throw WrongDim("Data<T, " + to_str(N) + ">::operator() const",
                     "operator() was called with 2 arguments.");
#endif

#ifdef SELDONDATA_DEBUG_CHECK_INDICES
    if ((i0 < 0) || (i0 >= data_.extent(0)))
      throw WrongIndex("Data<T, " + to_str(N) + ">::operator() const",
                       "Index along dimension #0 should be in [0, "
                       + to_str(data_.extent(0) - 1) + "], but is equal to "
                       + to_str(i0) + ".");
    if ((i1 < 0) || (i1 >= data_.extent(1)))
      throw WrongIndex("Data<T, " + to_str(N) + ">::operator() const",
                       "Index along dimension #1 should be in [0, "
                       + to_str(data_.extent(1) - 1) + "], but is equal to "
                       + to_str(i1) + ".");
#endif

    return data_(i0, i1);
  }

  //! Access operator for three-dimensional data.
  /*!
    \param i0 index for dimension #0.
    \param i1 index for dimension #1.
    \param i2 index for dimension #2.
    \return A const reference to element (i0, i1, i2).
    \exception SeldonData::WrongDim data dimension is not 3.
    \exception SeldonData::WrongIndex an index is out of range.
  */
  template<class T, int N, class TG>
  inline T Data<T, N, TG>::operator()(int i0, int i1, int i2)  const
  {

#ifdef SELDONDATA_DEBUG_CHECK_DIMENSIONS
    if (N != 3)
      throw WrongDim("Data<T, " + to_str(N) + ">::operator() const",
                     "operator() was called with 3 arguments.");
#endif

#ifdef SELDONDATA_DEBUG_CHECK_INDICES
    if ((i0 < 0) || (i0 >= data_.extent(0)))
      throw WrongIndex("Data<T, " + to_str(N) + ">::operator() const",
                       "Index along dimension #0 should be in [0, "
                       + to_str(data_.extent(0) - 1) + "], but is equal to "
                       + to_str(i0) + ".");
    if ((i1 < 0) || (i1 >= data_.extent(1)))
      throw WrongIndex("Data<T, " + to_str(N) + ">::operator() const",
                       "Index along dimension #1 should be in [0, "
                       + to_str(data_.extent(1) - 1) + "], but is equal to "
                       + to_str(i1) + ".");
    if ((i2 < 0) || (i2 >= data_.extent(2)))
      throw WrongIndex("Data<T, " + to_str(N) + ">::operator() const",
                       "Index along dimension #2 should be in [0, "
                       + to_str(data_.extent(2) - 1) + "], but is equal to "
                       + to_str(i2) + ".");
#endif

    return data_(i0, i1, i2);
  }

  //! Access operator for four-dimensional data.
  /*!
    \param i0 index for dimension #0.
    \param i1 index for dimension #1.
    \param i2 index for dimension #2.
    \param i3 index for dimension #3.
    \return A const reference to element (i0, i1, i2, i3).
    \exception SeldonData::WrongDim data dimension is not 4.
    \exception SeldonData::WrongIndex an index is out of range.
  */
  template<class T, int N, class TG>
  inline T Data<T, N, TG>::operator()(int i0, int i1, int i2,
                                      int i3)  const

  {

#ifdef SELDONDATA_DEBUG_CHECK_DIMENSIONS
    if (N != 4)
      throw WrongDim("Data<T, " + to_str(N) + ">::operator() const",
                     "operator() was called with 4 arguments.");
#endif

#ifdef SELDONDATA_DEBUG_CHECK_INDICES
    if ((i0 < 0) || (i0 >= data_.extent(0)))
      throw WrongIndex("Data<T, " + to_str(N) + ">::operator() const",
                       "Index along dimension #0 should be in [0, "
                       + to_str(data_.extent(0) - 1) + "], but is equal to "
                       + to_str(i0) + ".");
    if ((i1 < 0) || (i1 >= data_.extent(1)))
      throw WrongIndex("Data<T, " + to_str(N) + ">::operator() const",
                       "Index along dimension #1 should be in [0, "
                       + to_str(data_.extent(1) - 1) + "], but is equal to "
                       + to_str(i1) + ".");
    if ((i2 < 0) || (i2 >= data_.extent(2)))
      throw WrongIndex("Data<T, " + to_str(N) + ">::operator() const",
                       "Index along dimension #2 should be in [0, "
                       + to_str(data_.extent(2) - 1) + "], but is equal to "
                       + to_str(i2) + ".");
    if ((i3 < 0) || (i3 >= data_.extent(3)))
      throw WrongIndex("Data<T, " + to_str(N) + ">::operator() const",
                       "Index along dimension #3 should be in [0, "
                       + to_str(data_.extent(3) - 1) + "], but is equal to "
                       + to_str(i3) + ".");
#endif

    return data_(i0, i1, i2, i3);
  }

  //! Access operator for five-dimensional data.
  /*!
    \param i0 index for dimension #0.
    \param i1 index for dimension #1.
    \param i2 index for dimension #2.
    \param i3 index for dimension #3.
    \param i4 index for dimension #4.
    \return A const reference to element (i0, i1, i2, i3, i4).
    \exception SeldonData::WrongDim data dimension is not 5.
    \exception SeldonData::WrongIndex an index is out of range.
  */
  template<class T, int N, class TG>
  inline T Data<T, N, TG>::operator()(int i0, int i1, int i2,
                                      int i3, int i4)  const

  {

#ifdef SELDONDATA_DEBUG_CHECK_DIMENSIONS
    if (N != 5)
      throw WrongDim("Data<T, " + to_str(N) + ">::operator() const",
                     "operator() was called with 5 arguments.");
#endif

#ifdef SELDONDATA_DEBUG_CHECK_INDICES
    if ((i0 < 0) || (i0 >= data_.extent(0)))
      throw WrongIndex("Data<T, " + to_str(N) + ">::operator() const",
                       "Index along dimension #0 should be in [0, "
                       + to_str(data_.extent(0) - 1) + "], but is equal to "
                       + to_str(i0) + ".");
    if ((i1 < 0) || (i1 >= data_.extent(1)))
      throw WrongIndex("Data<T, " + to_str(N) + ">::operator() const",
                       "Index along dimension #1 should be in [0, "
                       + to_str(data_.extent(1) - 1) + "], but is equal to "
                       + to_str(i1) + ".");
    if ((i2 < 0) || (i2 >= data_.extent(2)))
      throw WrongIndex("Data<T, " + to_str(N) + ">::operator() const",
                       "Index along dimension #2 should be in [0, "
                       + to_str(data_.extent(2) - 1) + "], but is equal to "
                       + to_str(i2) + ".");
    if ((i3 < 0) || (i3 >= data_.extent(3)))
      throw WrongIndex("Data<T, " + to_str(N) + ">::operator() const",
                       "Index along dimension #3 should be in [0, "
                       + to_str(data_.extent(3) - 1) + "], but is equal to "
                       + to_str(i3) + ".");
    if ((i4 < 0) || (i4 >= data_.extent(4)))
      throw WrongIndex("Data<T, " + to_str(N) + ">::operator() const",
                       "Index along dimension #4 should be in [0, "
                       + to_str(data_.extent(4) - 1) + "], but is equal to "
                       + to_str(i4) + ".");
#endif

    return data_(i0, i1, i2, i3, i4);
  }

  //! Access operator for six-dimensional data.
  /*!
    \param i0 index for dimension #0.
    \param i1 index for dimension #1.
    \param i2 index for dimension #2.
    \param i3 index for dimension #3.
    \param i4 index for dimension #4.
    \param i5 index for dimension #5.
    \return A const reference to element (i0, i1, i2, i3, i4, i5).
    \exception SeldonData::WrongDim data dimension is not 6.
    \exception SeldonData::WrongIndex an index is out of range.
  */
  template<class T, int N, class TG>
  inline T Data<T, N, TG>::operator()(int i0, int i1, int i2,
                                      int i3, int i4, int i5)  const

  {

#ifdef SELDONDATA_DEBUG_CHECK_DIMENSIONS
    if (N != 6)
      throw WrongDim("Data<T, " + to_str(N) + ">::operator() const",
                     "operator() was called with 6 arguments.");
#endif

#ifdef SELDONDATA_DEBUG_CHECK_INDICES
    if ((i0 < 0) || (i0 >= data_.extent(0)))
      throw WrongIndex("Data<T, " + to_str(N) + ">::operator() const",
                       "Index along dimension #0 should be in [0, "
                       + to_str(data_.extent(0) - 1) + "], but is equal to "
                       + to_str(i0) + ".");
    if ((i1 < 0) || (i1 >= data_.extent(1)))
      throw WrongIndex("Data<T, " + to_str(N) + ">::operator() const",
                       "Index along dimension #1 should be in [0, "
                       + to_str(data_.extent(1) - 1) + "], but is equal to "
                       + to_str(i1) + ".");
    if ((i2 < 0) || (i2 >= data_.extent(2)))
      throw WrongIndex("Data<T, " + to_str(N) + ">::operator() const",
                       "Index along dimension #2 should be in [0, "
                       + to_str(data_.extent(2) - 1) + "], but is equal to "
                       + to_str(i2) + ".");
    if ((i3 < 0) || (i3 >= data_.extent(3)))
      throw WrongIndex("Data<T, " + to_str(N) + ">::operator() const",
                       "Index along dimension #3 should be in [0, "
                       + to_str(data_.extent(3) - 1) + "], but is equal to "
                       + to_str(i3) + ".");
    if ((i4 < 0) || (i4 >= data_.extent(4)))
      throw WrongIndex("Data<T, " + to_str(N) + ">::operator() const",
                       "Index along dimension #4 should be in [0, "
                       + to_str(data_.extent(4) - 1) + "], but is equal to "
                       + to_str(i4) + ".");
    if ((i5 < 0) || (i5 >= data_.extent(5)))
      throw WrongIndex("Data<T, " + to_str(N) + ">::operator() const",
                       "Index along dimension #5 should be in [0, "
                       + to_str(data_.extent(5) - 1) + "], but is equal to "
                       + to_str(i5) + ".");
#endif

    return data_(i0, i1, i2, i3, i4, i5);
  }

  template<class T, int N, class TG>
  //! Access operator for seven-dimensional data.
  /*!
    \param i0 index for dimension #0.
    \param i1 index for dimension #1.
    \param i2 index for dimension #2.
    \param i3 index for dimension #3.
    \param i4 index for dimension #4.
    \param i5 index for dimension #5.
    \param i6 index for dimension #6.
    \return A const reference to element (i0, i1, i2, i3,
    i4, i5, i6).
    \exception SeldonData::WrongDim data dimension is not 7.
    \exception SeldonData::WrongIndex an index is out of range.
  */
  inline T Data<T, N, TG>::operator()(int i0, int i1, int i2,
                                      int i3, int i4, int i5,
                                      int i6)  const

  {

#ifdef SELDONDATA_DEBUG_CHECK_DIMENSIONS
    if (N != 7)
      throw WrongDim("Data<T, " + to_str(N) + ">::operator() const",
                     "operator() was called with 7 arguments.");
#endif

#ifdef SELDONDATA_DEBUG_CHECK_INDICES
    if ((i0 < 0) || (i0 >= data_.extent(0)))
      throw WrongIndex("Data<T, " + to_str(N) + ">::operator() const",
                       "Index along dimension #0 should be in [0, "
                       + to_str(data_.extent(0) - 1) + "], but is equal to "
                       + to_str(i0) + ".");
    if ((i1 < 0) || (i1 >= data_.extent(1)))
      throw WrongIndex("Data<T, " + to_str(N) + ">::operator() const",
                       "Index along dimension #1 should be in [0, "
                       + to_str(data_.extent(1) - 1) + "], but is equal to "
                       + to_str(i1) + ".");
    if ((i2 < 0) || (i2 >= data_.extent(2)))
      throw WrongIndex("Data<T, " + to_str(N) + ">::operator() const",
                       "Index along dimension #2 should be in [0, "
                       + to_str(data_.extent(2) - 1) + "], but is equal to "
                       + to_str(i2) + ".");
    if ((i3 < 0) || (i3 >= data_.extent(3)))
      throw WrongIndex("Data<T, " + to_str(N) + ">::operator() const",
                       "Index along dimension #3 should be in [0, "
                       + to_str(data_.extent(3) - 1) + "], but is equal to "
                       + to_str(i3) + ".");
    if ((i4 < 0) || (i4 >= data_.extent(4)))
      throw WrongIndex("Data<T, " + to_str(N) + ">::operator() const",
                       "Index along dimension #4 should be in [0, "
                       + to_str(data_.extent(4) - 1) + "], but is equal to "
                       + to_str(i4) + ".");
    if ((i5 < 0) || (i5 >= data_.extent(5)))
      throw WrongIndex("Data<T, " + to_str(N) + ">::operator() const",
                       "Index along dimension #5 should be in [0, "
                       + to_str(data_.extent(5) - 1) + "], but is equal to "
                       + to_str(i5) + ".");
    if ((i6 < 0) || (i6 >= data_.extent(6)))
      throw WrongIndex("Data<T, " + to_str(N) + ">::operator() const",
                       "Index along dimension #6 should be in [0, "
                       + to_str(data_.extent(6) - 1) + "], but is equal to "
                       + to_str(i6) + ".");
#endif

    return data_(i0, i1, i2, i3, i4, i5, i6);
  }

  template<class T, int N, class TG>
  //! Access operator for eight-dimensional data.
  /*!
    \param i0 index for dimension #0.
    \param i1 index for dimension #1.
    \param i2 index for dimension #2.
    \param i3 index for dimension #3.
    \param i4 index for dimension #4.
    \param i5 index for dimension #5.
    \param i6 index for dimension #6.
    \param i7 index for dimension #7.
    \return A const reference to element (i0, i1, i2, i3, i4,
    i5, i6, i7).
    \exception SeldonData::WrongDim data dimension is not 8.
    \exception SeldonData::WrongIndex an index is out of range.
  */
  inline T Data<T, N, TG>::operator()(int i0, int i1, int i2,
                                      int i3, int i4, int i5,
                                      int i6, int i7)  const

  {

#ifdef SELDONDATA_DEBUG_CHECK_DIMENSIONS
    if (N != 8)
      throw WrongDim("Data<T, " + to_str(N) + ">::operator() const",
                     "operator() was called with 8 arguments.");
#endif

#ifdef SELDONDATA_DEBUG_CHECK_INDICES
    if ((i0 < 0) || (i0 >= data_.extent(0)))
      throw WrongIndex("Data<T, " + to_str(N) + ">::operator() const",
                       "Index along dimension #0 should be in [0, "
                       + to_str(data_.extent(0) - 1) + "], but is equal to "
                       + to_str(i0) + ".");
    if ((i1 < 0) || (i1 >= data_.extent(1)))
      throw WrongIndex("Data<T, " + to_str(N) + ">::operator() const",
                       "Index along dimension #1 should be in [0, "
                       + to_str(data_.extent(1) - 1) + "], but is equal to "
                       + to_str(i1) + ".");
    if ((i2 < 0) || (i2 >= data_.extent(2)))
      throw WrongIndex("Data<T, " + to_str(N) + ">::operator() const",
                       "Index along dimension #2 should be in [0, "
                       + to_str(data_.extent(2) - 1) + "], but is equal to "
                       + to_str(i2) + ".");
    if ((i3 < 0) || (i3 >= data_.extent(3)))
      throw WrongIndex("Data<T, " + to_str(N) + ">::operator() const",
                       "Index along dimension #3 should be in [0, "
                       + to_str(data_.extent(3) - 1) + "], but is equal to "
                       + to_str(i3) + ".");
    if ((i4 < 0) || (i4 >= data_.extent(4)))
      throw WrongIndex("Data<T, " + to_str(N) + ">::operator() const",
                       "Index along dimension #4 should be in [0, "
                       + to_str(data_.extent(4) - 1) + "], but is equal to "
                       + to_str(i4) + ".");
    if ((i5 < 0) || (i5 >= data_.extent(5)))
      throw WrongIndex("Data<T, " + to_str(N) + ">::operator() const",
                       "Index along dimension #5 should be in [0, "
                       + to_str(data_.extent(5) - 1) + "], but is equal to "
                       + to_str(i5) + ".");
    if ((i6 < 0) || (i6 >= data_.extent(6)))
      throw WrongIndex("Data<T, " + to_str(N) + ">::operator() const",
                       "Index along dimension #6 should be in [0, "
                       + to_str(data_.extent(6) - 1) + "], but is equal to "
                       + to_str(i6) + ".");
    if ((i7 < 0) || (i7 >= data_.extent(7)))
      throw WrongIndex("Data<T, " + to_str(N) + ">::operator() const",
                       "Index along dimension #7 should be in [0, "
                       + to_str(data_.extent(7) - 1) + "], but is equal to "
                       + to_str(i7) + ".");
#endif

    return data_(i0, i1, i2, i3, i4, i5, i6, i7);
  }

  template<class T, int N, class TG>
  //! Access operator for nine-dimensional data.
  /*!
    \param i0 index for dimension #0.
    \param i1 index for dimension #1.
    \param i2 index for dimension #2.
    \param i3 index for dimension #3.
    \param i4 index for dimension #4.
    \param i5 index for dimension #5.
    \param i6 index for dimension #6.
    \param i7 index for dimension #7.
    \param i8 index for dimension #8.
    \return A const reference to element (i0, i1, i2, i3, i4,
    i5, i6, i7, i8).
    \exception SeldonData::WrongDim data dimension is not 9.
    \exception SeldonData::WrongIndex an index is out of range.
  */
  inline T Data<T, N, TG>::operator()(int i0, int i1, int i2,
                                      int i3, int i4, int i5,
                                      int i6, int i7, int i8)  const

  {

#ifdef SELDONDATA_DEBUG_CHECK_DIMENSIONS
    if (N != 9)
      throw WrongDim("Data<T, " + to_str(N) + ">::operator() const",
                     "operator() was called with 9 arguments.");
#endif

#ifdef SELDONDATA_DEBUG_CHECK_INDICES
    if ((i0 < 0) || (i0 >= data_.extent(0)))
      throw WrongIndex("Data<T, " + to_str(N) + ">::operator() const",
                       "Index along dimension #0 should be in [0, "
                       + to_str(data_.extent(0) - 1) + "], but is equal to "
                       + to_str(i0) + ".");
    if ((i1 < 0) || (i1 >= data_.extent(1)))
      throw WrongIndex("Data<T, " + to_str(N) + ">::operator() const",
                       "Index along dimension #1 should be in [0, "
                       + to_str(data_.extent(1) - 1) + "], but is equal to "
                       + to_str(i1) + ".");
    if ((i2 < 0) || (i2 >= data_.extent(2)))
      throw WrongIndex("Data<T, " + to_str(N) + ">::operator() const",
                       "Index along dimension #2 should be in [0, "
                       + to_str(data_.extent(2) - 1) + "], but is equal to "
                       + to_str(i2) + ".");
    if ((i3 < 0) || (i3 >= data_.extent(3)))
      throw WrongIndex("Data<T, " + to_str(N) + ">::operator() const",
                       "Index along dimension #3 should be in [0, "
                       + to_str(data_.extent(3) - 1) + "], but is equal to "
                       + to_str(i3) + ".");
    if ((i4 < 0) || (i4 >= data_.extent(4)))
      throw WrongIndex("Data<T, " + to_str(N) + ">::operator() const",
                       "Index along dimension #4 should be in [0, "
                       + to_str(data_.extent(4) - 1) + "], but is equal to "
                       + to_str(i4) + ".");
    if ((i5 < 0) || (i5 >= data_.extent(5)))
      throw WrongIndex("Data<T, " + to_str(N) + ">::operator() const",
                       "Index along dimension #5 should be in [0, "
                       + to_str(data_.extent(5) - 1) + "], but is equal to "
                       + to_str(i5) + ".");
    if ((i6 < 0) || (i6 >= data_.extent(6)))
      throw WrongIndex("Data<T, " + to_str(N) + ">::operator() const",
                       "Index along dimension #6 should be in [0, "
                       + to_str(data_.extent(6) - 1) + "], but is equal to "
                       + to_str(i6) + ".");
    if ((i7 < 0) || (i7 >= data_.extent(7)))
      throw WrongIndex("Data<T, " + to_str(N) + ">::operator() const",
                       "Index along dimension #7 should be in [0, "
                       + to_str(data_.extent(7) - 1) + "], but is equal to "
                       + to_str(i7) + ".");
    if ((i8 < 0) || (i8 >= data_.extent(8)))
      throw WrongIndex("Data<T, " + to_str(N) + ">::operator() const",
                       "Index along dimension #8 should be in [0, "
                       + to_str(data_.extent(8) - 1) + "], but is equal to "
                       + to_str(i8) + ".");
#endif

    return data_(i0, i1, i2, i3, i4, i5, i6, i7, i8);
  }

  //! Access operator for ten-dimensional data.
  /*!
    \param i0 index for dimension #0.
    \param i1 index for dimension #1.
    \param i2 index for dimension #2.
    \param i3 index for dimension #3.
    \param i4 index for dimension #4.
    \param i5 index for dimension #5.
    \param i6 index for dimension #6.
    \param i7 index for dimension #7.
    \param i8 index for dimension #8.
    \param i9 index for dimension #9.
    \return A const reference to element (i0, i1, i2, i3, i4,
    i5, i6, i7, i8, i9).
    \exception SeldonData::WrongDim data dimension is not 10.
    \exception SeldonData::WrongIndex an index is out of range.
  */
  template<class T, int N, class TG>
  inline T Data<T, N, TG>::operator()(int i0, int i1, int i2,
                                      int i3, int i4, int i5,
                                      int i6, int i7, int i8,
                                      int i9)  const

  {

#ifdef SELDONDATA_DEBUG_CHECK_DIMENSIONS
    if (N != 10)
      throw WrongDim("Data<T, " + to_str(N) + ">::operator() const",
                     "operator() was called with 10 arguments.");
#endif

#ifdef SELDONDATA_DEBUG_CHECK_INDICES
    if ((i0 < 0) || (i0 >= data_.extent(0)))
      throw WrongIndex("Data<T, " + to_str(N) + ">::operator() const",
                       "Index along dimension #0 should be in [0, "
                       + to_str(data_.extent(0) - 1) + "], but is equal to "
                       + to_str(i0) + ".");
    if ((i1 < 0) || (i1 >= data_.extent(1)))
      throw WrongIndex("Data<T, " + to_str(N) + ">::operator() const",
                       "Index along dimension #1 should be in [0, "
                       + to_str(data_.extent(1) - 1) + "], but is equal to "
                       + to_str(i1) + ".");
    if ((i2 < 0) || (i2 >= data_.extent(2)))
      throw WrongIndex("Data<T, " + to_str(N) + ">::operator() const",
                       "Index along dimension #2 should be in [0, "
                       + to_str(data_.extent(2) - 1) + "], but is equal to "
                       + to_str(i2) + ".");
    if ((i3 < 0) || (i3 >= data_.extent(3)))
      throw WrongIndex("Data<T, " + to_str(N) + ">::operator() const",
                       "Index along dimension #3 should be in [0, "
                       + to_str(data_.extent(3) - 1) + "], but is equal to "
                       + to_str(i3) + ".");
    if ((i4 < 0) || (i4 >= data_.extent(4)))
      throw WrongIndex("Data<T, " + to_str(N) + ">::operator() const",
                       "Index along dimension #4 should be in [0, "
                       + to_str(data_.extent(4) - 1) + "], but is equal to "
                       + to_str(i4) + ".");
    if ((i5 < 0) || (i5 >= data_.extent(5)))
      throw WrongIndex("Data<T, " + to_str(N) + ">::operator() const",
                       "Index along dimension #5 should be in [0, "
                       + to_str(data_.extent(5) - 1) + "], but is equal to "
                       + to_str(i5) + ".");
    if ((i6 < 0) || (i6 >= data_.extent(6)))
      throw WrongIndex("Data<T, " + to_str(N) + ">::operator() const",
                       "Index along dimension #6 should be in [0, "
                       + to_str(data_.extent(6) - 1) + "], but is equal to "
                       + to_str(i6) + ".");
    if ((i7 < 0) || (i7 >= data_.extent(7)))
      throw WrongIndex("Data<T, " + to_str(N) + ">::operator() const",
                       "Index along dimension #7 should be in [0, "
                       + to_str(data_.extent(7) - 1) + "], but is equal to "
                       + to_str(i7) + ".");
    if ((i8 < 0) || (i8 >= data_.extent(8)))
      throw WrongIndex("Data<T, " + to_str(N) + ">::operator() const",
                       "Index along dimension #8 should be in [0, "
                       + to_str(data_.extent(8) - 1) + "], but is equal to "
                       + to_str(i8) + ".");
    if ((i9 < 0) || (i9 >= data_.extent(9)))
      throw WrongIndex("Data<T, " + to_str(N) + ">::operator() const",
                       "Index along dimension #9 should be in [0, "
                       + to_str(data_.extent(9) - 1) + "], but is equal to "
                       + to_str(i9) + ".");
#endif

    return data_(i0, i1, i2, i3, i4, i5, i6, i7, i8, i9);
  }

  //! Access operator.
  /*!
    \param indices array of indices.
    \return The value of current data at (indices(0), indices(1), ...).
  */
  template<class T, int N, class TG>
  inline T& Data<T, N, TG>::operator()(const Array<int, 1>& indices)
  {

#ifdef SELDONDATA_DEBUG_CHECK_DIMENSIONS
    if (N != indices.numElements())
      throw WrongDim("Data<T, " + to_str(N) + ">::operator() (Array<int, 1>&)",
                     "The input array contains " + to_str(indices.numElements())
                     + " indice(s) instead of " + to_str(N) + ".");
#endif

    if (N == 1)
      return data_(indices(0));
    else if (N == 2)
      return data_(indices(0), indices(1));
    else if (N == 3)
      return data_(indices(0), indices(1), indices(2));
    else if (N == 4)
      return data_(indices(0), indices(1), indices(2),
                   indices(3));
    else if (N == 5)
      return data_(indices(0), indices(1), indices(2),
                   indices(3), indices(4));
    else if (N == 6)
      return data_(indices(0), indices(1), indices(2),
                   indices(3), indices(4), indices(5));
    else if (N == 7)
      return data_(indices(0), indices(1), indices(2),
                   indices(3), indices(4), indices(5),
                   indices(6));
    else if (N == 8)
      return data_(indices(0), indices(1), indices(2),
                   indices(3), indices(4), indices(5),
                   indices(6), indices(7));
    else if (N == 9)
      return data_(indices(0), indices(1), indices(2),
                   indices(3), indices(4), indices(5),
                   indices(6), indices(7), indices(8));
    else if (N == 10)
      return data_(indices(0), indices(1), indices(2),
                   indices(3), indices(4), indices(5),
                   indices(6), indices(7), indices(8),
                   indices(9));
  }

  //! General access function.
  /*!
    Only useful arguments are taken into account.
    \param i0 index for dimension #0.
    \param i1 index for dimension #1.
    \param i2 index for dimension #2.
    \param i3 index for dimension #3.
    \param i4 index for dimension #4.
    \param i5 index for dimension #5.
    \param i6 index for dimension #6.
    \param i7 index for dimension #7.
    \param i8 index for dimension #8.
    \param i9 index for dimension #9.
    \return A reference to the required element.
  */
  template<class T, int N, class TG>
  inline T& Data<T, N, TG>::Value(int i0, int i1, int i2,
                                  int i3, int i4, int i5,
                                  int i6, int i7, int i8,
                                  int i9)
  {
    if (N == 1)
      return (*this)(i0);
    if (N == 2)
      return (*this)(i0, i1);
    if (N == 3)
      return (*this)(i0, i1, i2);
    if (N == 4)
      return (*this)(i0, i1, i2, i3);
    if (N == 5)
      return (*this)(i0, i1, i2, i3, i4);
    if (N == 6)
      return (*this)(i0, i1, i2, i3, i4, i5);
    if (N == 7)
      return (*this)(i0, i1, i2, i3, i4, i5, i6);
    if (N == 8)
      return (*this)(i0, i1, i2, i3, i4, i5, i6, i7);
    if (N == 9)
      return (*this)(i0, i1, i2, i3, i4, i5, i6, i7, i8);
    if (N == 10)
      return (*this)(i0, i1, i2, i3, i4, i5, i6, i7, i8, i9);
  }

  //! General access function.
  /*!
    Only useful arguments are taken into account.
    \param i0 index for dimension #0.
    \param i1 index for dimension #1.
    \param i2 index for dimension #2.
    \param i3 index for dimension #3.
    \param i4 index for dimension #4.
    \param i5 index for dimension #5.
    \param i6 index for dimension #6.
    \param i7 index for dimension #7.
    \param i8 index for dimension #8.
    \param i9 index for dimension #9.
    \return Value of the required element.
  */
  template<class T, int N, class TG>
  inline T Data<T, N, TG>::Value(int i0, int i1, int i2,
                                 int i3, int i4, int i5,
                                 int i6, int i7, int i8,
                                 int i9) const
  {
    if (N == 1)
      return (*this)(i0);
    if (N == 2)
      return (*this)(i0, i1);
    if (N == 3)
      return (*this)(i0, i1, i2);
    if (N == 4)
      return (*this)(i0, i1, i2, i3);
    if (N == 5)
      return (*this)(i0, i1, i2, i3, i4);
    if (N == 6)
      return (*this)(i0, i1, i2, i3, i4, i5);
    if (N == 7)
      return (*this)(i0, i1, i2, i3, i4, i5, i6);
    if (N == 8)
      return (*this)(i0, i1, i2, i3, i4, i5, i6, i7);
    if (N == 9)
      return (*this)(i0, i1, i2, i3, i4, i5, i6, i7, i8);
    if (N == 10)
      return (*this)(i0, i1, i2, i3, i4, i5, i6, i7, i8, i9);
  }

  //! Performs a copy (duplication) of \a data.
  /*!
    \param data Data instance to be copied.
    \warning Grids may or may not be duplicated (in memory), depending on
    their own settings.
  */
  template<class T, int N, class TG>
  template <class T0>
  void Data<T, N, TG>::Copy(Data<T0, N, TG>& data)
  {
    for (int i = 0; i < N; i++)
      {
        if (grids_(i) != NULL)
          if (grids_(i)->GetPointers() == 1)
            delete grids_(i);
          else
            grids_(i)->SetPointers(grids_(i)->GetPointers() - 1);
        grids_(i) = data[i].Copy();
      }

    data_.resize(data.GetArray().shape());
    data_ = data.GetArray().copy();
  }

  //! Performs a reference copy of \a data.
  /*! On exit, current Data instance and \a data share the same data array.
    \param data Data instance to be copied.
    \warning Grids may or may not be duplicated (in memory), depending on
    their own settings.
  */
  template<class T, int N, class TG>
  template <class T0>
  void Data<T, N, TG>::ReferenceCopy(Data<T0, N, TG>& data)
  {
    for (int i = 0; i < N; i++)
      {
        if (grids_(i) != NULL)
          if (grids_(i)->GetPointers() == 1)
            delete grids_(i);
          else
            grids_(i)->SetPointers(grids_(i)->GetPointers() - 1);
        grids_(i) = data[i].Copy();
      }

    data_.resize(data.GetArray().shape());
    data_ = data.GetArray();
  }

  //! Returns the number of stored elements.
  /*!
    \return Number of elements.
  */
  template<class T, int N, class TG>
  inline int Data<T, N, TG>::GetNbElements()
  {
    return data_.numElements();
  }

  //! Returns the number of dimensions.
  /*!
    \return Number of dimensions.
  */
  template<class T, int N, class TG>
  inline int Data<T, N, TG>::GetNbDim()
  {
    return N;
  }

  //! Returns length along dimension #i.
  /*!
    \param i dimension number.
    \return Length along the i-th dimension.
  */
  template<class T, int N, class TG>
  inline int Data<T, N, TG>::GetLength(int i) const
  {
    if ((i >= 0) && (i < N))
      return data_.extent(i);
    else
      return 0;
  }

  //! Returns a pointer to the grid #i.
  /*!
    \param i grid to be returned.
    \return A pointer to the grid #i.
  */
  template<class T, int N, class TG>
  inline Grid<TG>* Data<T, N, TG>::GetGrid(int i)
  {
    return grids_(i);
  }

  //! Returns a reference to the grid #i.
  /*!
    \param i grid to be returned.
    \return A reference to the grid #i.
  */
  template<class T, int N, class TG>
  inline Grid<TG>& Data<T, N, TG>::operator[](int i)
  {
    return *grids_(i);
  }

  //! Returns a reference to the grid #i.
  /*!
    \param i grid to be returned.
    \return A const reference to the grid #i.
  */
  template<class T, int N, class TG>
  inline const Grid<TG>& Data<T, N, TG>::operator[](int i) const
  {
    return *grids_(i);
  }

  //! Returns the array of pointers to grids.
  /*!
    For every Data, an array gathers pointers to grids.
    The i-th element of the array is the grid related to
    dimension #i.
    \return A reference to the array of pointers to grids.
  */
  template<class T, int N, class TG>
  inline Array<Grid<TG>*, 1>& Data<T, N, TG>::GetGrids()
  {
    return grids_;
  }

  //! Returns the blitz::Array storing elements.
  /*!
    \return A reference to the blitz::Array storing data.
  */
  template<class T, int N, class TG>
  inline Array<T, N>& Data<T, N, TG>::GetArray()
  {
    return data_;
  }

  //! Returns the blitz::Array storing elements.
  /*!
    \return A reference to the blitz::Array storing data.
  */
  template<class T, int N, class TG>
  inline const Array<T, N>& Data<T, N, TG>::GetArray() const
  {
    return data_;
  }

  //! Returns the blitz::Array storing elements.
  /*!
    \return A reference to the blitz::Array storing data.
  */
  template<class T, int N, class TG>
  inline Array<T, N>& Data<T, N, TG>::operator()()
  {
    return data_;
  }

  //! Returns the array storing elements.
  /*!
    \return A pointer to the array storing data.
  */
  template<class T, int N, class TG>
  inline T* Data<T, N, TG>::GetData()
  {
    return data_.data();
  }

  //! Returns the array storing elements.
  /*!
    \return A constant pointer to the array storing data.
  */
  template<class T, int N, class TG>
  inline const T* Data<T, N, TG>::GetData() const
  {
    return data_.data();
  }

  //! Extracts data.
  /*!
    Current data (*this) is set to a sub-data of 'data'.
    \param data source data.
    \param r0 range for dimension #0 (source).
  */
  template<class T, int N, class TG>
  template<class DTG, class R0>
  void Data<T, N, TG>::SubData(Data<T, 1, DTG>& data, R0 r0)
  {
    Array<T, 1>& array = data.GetArray();
    this->data_ = array(r0);
  }

  //! Extracts data.
  /*!
    Current data (*this) is set to a sub-data of 'data'.
    \param data source data.
    \param r0 range for dimension #0 (source).
    \param r1 range for dimension #1 (source).
  */
  template<class T, int N, class TG>
  template<class DTG, class R0, class R1>
  void Data<T, N, TG>::SubData(Data<T, 2, DTG>& data, R0 r0, R1 r1)
  {
    Array<T, 2>& array = data.GetArray();
    this->data_ = array(r0, r1);
  }

  //! Extracts data.
  /*!
    Current data (*this) is set to a sub-data of 'data'.
    \param data source data.
    \param r0 range for dimension #0 (source).
    \param r1 range for dimension #1 (source).
    \param r2 range for dimension #2 (source).
  */
  template<class T, int N, class TG>
  template<class DTG, class R0, class R1, class R2>
  void Data<T, N, TG>::SubData(Data<T, 3, DTG>& data, R0 r0, R1 r1, R2 r2)
  {
    Array<T, 3>& array = data.GetArray();
    this->data_ = array(r0, r1, r2);
  }

  //! Extracts data.
  /*!
    Current data (*this) is set to a sub-data of 'data'.
    \param data source data.
    \param r0 range for dimension #0 (source).
    \param r1 range for dimension #1 (source).
    \param r2 range for dimension #2 (source).
    \param r3 range for dimension #3 (source).
  */
  template<class T, int N, class TG>
  template < class DTG, class R0, class R1, class R2,
             class R3 >
  void Data<T, N, TG>::SubData(Data<T, 4, DTG>& data, R0 r0, R1 r1, R2 r2,
                               R3 r3)
  {
    Array<T, 4>& array = data.GetArray();
    this->data_ = array(r0, r1, r2,
                        r3);
  }

  //! Extracts data.
  /*!
    Current data (*this) is set to a sub-data of 'data'.
    \param data source data.
    \param r0 range for dimension #0 (source).
    \param r1 range for dimension #1 (source).
    \param r2 range for dimension #2 (source).
    \param r3 range for dimension #3 (source).
    \param r4 range for dimension #4 (source).
  */
  template<class T, int N, class TG>
  template < class DTG, class R0, class R1, class R2,
             class R3, class R4 >
  void Data<T, N, TG>::SubData(Data<T, 5, DTG>& data, R0 r0, R1 r1, R2 r2,
                               R3 r3, R4 r4)
  {
    Array<T, 5>& array = data.GetArray();
    this->data_ = array(r0, r1, r2,
                        r3, r4);
  }

  //! Extracts data.
  /*!
    Current data (*this) is set to a sub-data of 'data'.
    \param data source data.
    \param r0 range for dimension #0 (source).
    \param r1 range for dimension #1 (source).
    \param r2 range for dimension #2 (source).
    \param r3 range for dimension #3 (source).
    \param r4 range for dimension #4 (source).
    \param r5 range for dimension #5 (source).
  */
  template<class T, int N, class TG>
  template < class DTG, class R0, class R1, class R2,
             class R3, class R4, class R5 >
  void Data<T, N, TG>::SubData(Data<T, 6, DTG>& data, R0 r0, R1 r1, R2 r2,
                               R3 r3, R4 r4, R5 r5)
  {
    Array<T, 6>& array = data.GetArray();
    this->data_ = array(r0, r1, r2,
                        r3, r4, r5);
  }

  //! Extracts data.
  /*!
    Current data (*this) is set to a sub-data of 'data'.
    \param data source data.
    \param r0 range for dimension #0 (source).
    \param r1 range for dimension #1 (source).
    \param r2 range for dimension #2 (source).
    \param r3 range for dimension #3 (source).
    \param r4 range for dimension #4 (source).
    \param r5 range for dimension #5 (source).
    \param r6 range for dimension #6 (source).
  */
  template<class T, int N, class TG>
  template < class DTG, class R0, class R1, class R2,
             class R3, class R4, class R5,
             class R6 >
  void Data<T, N, TG>::SubData(Data<T, 7, DTG>& data, R0 r0, R1 r1, R2 r2,
                               R3 r3, R4 r4, R5 r5, R6 r6)
  {
    Array<T, 7>& array = data.GetArray();
    this->data_ = array(r0, r1, r2,
                        r3, r4, r5, r6);
  }

  //! Extracts data.
  /*!
    Current data (*this) is set to a sub-data of 'data'.
    \param data source data.
    \param r0 range for dimension #0 (source).
    \param r1 range for dimension #1 (source).
    \param r2 range for dimension #2 (source).
    \param r3 range for dimension #3 (source).
    \param r4 range for dimension #4 (source).
    \param r5 range for dimension #5 (source).
    \param r6 range for dimension #6 (source).
    \param r7 range for dimension #7 (source).
  */
  template<class T, int N, class TG>
  template < class DTG, class R0, class R1, class R2,
             class R3, class R4, class R5,
             class R6, class R7 >
  void Data<T, N, TG>::SubData(Data<T, 8, DTG>& data, R0 r0, R1 r1, R2 r2,
                               R3 r3, R4 r4, R5 r5, R6 r6,
                               R7 r7)
  {
    Array<T, 8>& array = data.GetArray();
    this->data_ = array(r0, r1, r2,
                        r3, r4, r5, r6,
                        r7);
  }

  //! Extracts data.
  /*!
    Current data (*this) is set to a sub-data of 'data'.
    \param data source data.
    \param r0 range for dimension #0 (source).
    \param r1 range for dimension #1 (source).
    \param r2 range for dimension #2 (source).
    \param r3 range for dimension #3 (source).
    \param r4 range for dimension #4 (source).
    \param r5 range for dimension #5 (source).
    \param r6 range for dimension #6 (source).
    \param r7 range for dimension #7 (source).
    \param r8 range for dimension #8 (source).
  */
  template<class T, int N, class TG>
  template < class DTG, class R0, class R1, class R2,
             class R3, class R4, class R5,
             class R6, class R7, class R8 >
  void Data<T, N, TG>::SubData(Data<T, 9, DTG>& data, R0 r0, R1 r1, R2 r2,
                               R3 r3, R4 r4, R5 r5, R6 r6,
                               R7 r7, R8 r8)
  {
    Array<T, 9>& array = data.GetArray();
    this->data_ = array(r0, r1, r2,
                        r3, r4, r5, r6,
                        r7, r8);
  }

  //! Extracts data.
  /*!
    Current data (*this) is set to a sub-data of 'data'.
    \param data source data.
    \param r0 range for dimension #0 (source).
    \param r1 range for dimension #1 (source).
    \param r2 range for dimension #2 (source).
    \param r3 range for dimension #3 (source).
    \param r4 range for dimension #4 (source).
    \param r5 range for dimension #5 (source).
    \param r6 range for dimension #6 (source).
    \param r7 range for dimension #7 (source).
    \param r8 range for dimension #8 (source).
    \param r9 range for dimension #9 (source).
  */
  template<class T, int N, class TG>
  template < class DTG, class R0, class R1, class R2,
             class R3, class R4, class R5,
             class R6, class R7, class R8, class R9 >
  void Data<T, N, TG>::SubData(Data<T, 10, DTG>& data, R0 r0, R1 r1, R2 r2,
                               R3 r3, R4 r4, R5 r5, R6 r6,
                               R7 r7, R8 r8, R9 r9)
  {
    Array<T, 10>& array = data.GetArray();
    this->data_ = array(r0, r1, r2,
                        r3, r4, r5, r6,
                        r7, r8, r9);
  }

  //! Resize the grid.
  /*!
    \param G0 grid for dimension #0.
  */
  template<class T, int N, class TG>
  void Data<T, N, TG>::ResizeGrid(Grid<TG>& G0)
  {

#ifdef SELDONDATA_DEBUG_CHECK_DIMENSIONS
    if (N != 1)
      {
        cout << "ERROR!" << endl;
        cout << "Wrong dimension in Data<T, N, TG>::ResizeGrid(Grid<TG>, ...)."
             << endl << "   Required " + to_str(N)
          + " grids, but got 1 grid." << endl;
        abort();
      }
#endif

    this->ClearGrids();
    grids_(0) = G0.Copy();
    SetVariables();

#ifdef SELDONDATA_DEBUG_CHECK_DIMENSIONS
    int i, j;
    for (i = 0; i < N; i++)
      for (j = 0; j < N; j++)
        if ((grids_(i)->GetLength(j) != 0)
            && (grids_(i)->GetLength(j) != grids_(j)->GetLength()))
          {
            cout << "ERROR!" << endl;
            cout << "Wrong dimension in Data<T, N, TG>::ResizeGrid(Grid<TG>, ...)."
                 << endl << "   Length of grid #" << i << " along dimension #"
                 << j << " is " << grids_(i)->GetLength(j) << " but should be "
                 << grids_(j)->GetLength() << "." << endl;
            abort();
          }

    for (i = 0; i < N; i++)
      if (grids_(i)->GetLength(i) != data_.extent(i))
        throw WrongDim("Data<T, " + to_str(N) + ">::ResizeGrid(Grid<TG>, ...)",
                       "Length of grid #" + to_str(i) + " is "
                       + to_str(grids_(i)->GetLength(i)) + " but should be "
                       + to_str(data_.extent(i)) + ".");
#endif

  }

  //! Resizes grids.
  /*!
    \param G0 grid for dimension #0.
    \param G1 grid for dimension #1.
  */
  template<class T, int N, class TG>
  void Data<T, N, TG>::ResizeGrid(Grid<TG>& G0, Grid<TG>& G1)
  {

#ifdef SELDONDATA_DEBUG_CHECK_DIMENSIONS
    if (N != 2)
      {
        cout << "ERROR!" << endl;
        cout << "Wrong dimension in Data<T, N, TG>::ResizeGrid(Grid<TG>, ...)."
             << endl << "   Required " + to_str(N)
          + " grids, but got 2 grids." << endl;
        abort();
      }
#endif

    this->ClearGrids();
    grids_(0) = G0.Copy();
    grids_(1) = G1.Copy();
    SetVariables();

#ifdef SELDONDATA_DEBUG_CHECK_DIMENSIONS
    int i, j;
    for (i = 0; i < N; i++)
      for (j = 0; j < N; j++)
        if ((grids_(i)->GetLength(j) != 0)
            && (grids_(i)->GetLength(j) != grids_(j)->GetLength()))
          {
            cout << "ERROR!" << endl;
            cout << "Wrong dimension in Data<T, N, TG>::ResizeGrid(Grid<TG>, ...)."
                 << endl << "   Length of grid #" << i << " along dimension #"
                 << j << " is " << grids_(i)->GetLength(j) << " but should be "
                 << grids_(j)->GetLength() << "." << endl;
            abort();
          }

    for (i = 0; i < N; i++)
      if (grids_(i)->GetLength(i) != data_.extent(i))
        throw WrongDim("Data<T, " + to_str(N) + ">::ResizeGrid(Grid<TG>, ...)",
                       "Length of grid #" + to_str(i) + " is "
                       + to_str(grids_(i)->GetLength(i)) + " but should be "
                       + to_str(data_.extent(i)) + ".");
#endif

  }

  //! Resizes the grids.
  /*!
    \param G0 grid for dimension #0.
    \param G1 grid for dimension #1.
    \param G2 grid for dimension #2.
  */
  template<class T, int N, class TG>
  void Data<T, N, TG>::ResizeGrid(Grid<TG>& G0, Grid<TG>& G1, Grid<TG>& G2)
  {

#ifdef SELDONDATA_DEBUG_CHECK_DIMENSIONS
    if (N != 3)
      {
        cout << "ERROR!" << endl;
        cout << "Wrong dimension in Data<T, N, TG>::ResizeGrid(Grid<TG>, ...)."
             << endl << "   Required " + to_str(N)
          + " grids, but got 3 grids." << endl;
        abort();
      }
#endif

    this->ClearGrids();
    grids_(0) = G0.Copy();
    grids_(1) = G1.Copy();
    grids_(2) = G2.Copy();
    SetVariables();


#ifdef SELDONDATA_DEBUG_CHECK_DIMENSIONS
    int i, j;
    for (i = 0; i < N; i++)
      for (j = 0; j < N; j++)
        if ((grids_(i)->GetLength(j) != 0)
            && (grids_(i)->GetLength(j) != grids_(j)->GetLength()))
          {
            cout << "ERROR!" << endl;
            cout << "Wrong dimension in Data<T, N, TG>::ResizeGrid(Grid<TG>, ...)."
                 << endl << "   Length of grid #" << i << " along dimension #"
                 << j << " is " << grids_(i)->GetLength(j) << " but should be "
                 << grids_(j)->GetLength() << "." << endl;
            abort();
          }

    for (i = 0; i < N; i++)
      if (grids_(i)->GetLength(i) != data_.extent(i))
        throw WrongDim("Data<T, " + to_str(N) + ">::ResizeGrid(Grid<TG>, ...)",
                       "Length of grid #" + to_str(i) + " is "
                       + to_str(grids_(i)->GetLength(i)) + " but should be "
                       + to_str(data_.extent(i)) + ".");
#endif

  }

  //! Resizes the grids.
  /*!
    \param G0 grid for dimension #0.
    \param G1 grid for dimension #1.
    \param G2 grid for dimension #2.
    \param G3 grid for dimension #3.
  */
  template<class T, int N, class TG>
  void Data<T, N, TG>::ResizeGrid(Grid<TG>& G0, Grid<TG>& G1, Grid<TG>& G2,
                                  Grid<TG>& G3)
  {

#ifdef SELDONDATA_DEBUG_CHECK_DIMENSIONS
    if (N != 4)
      {
        cout << "ERROR!" << endl;
        cout << "Wrong dimension in Data<T, N, TG>::ResizeGrid(Grid<TG>, ...)."
             << endl << "   Required " + to_str(N)
          + " grids, but got 4 grids." << endl;
        abort();
      }
#endif

    this->ClearGrids();
    grids_(0) = G0.Copy();
    grids_(1) = G1.Copy();
    grids_(2) = G2.Copy();
    grids_(3) = G3.Copy();
    SetVariables();


#ifdef SELDONDATA_DEBUG_CHECK_DIMENSIONS
    int i, j;
    for (i = 0; i < N; i++)
      for (j = 0; j < N; j++)
        if ((grids_(i)->GetLength(j) != 0)
            && (grids_(i)->GetLength(j) != grids_(j)->GetLength()))
          {
            cout << "ERROR!" << endl;
            cout << "Wrong dimension in Data<T, N, TG>::ResizeGrid(Grid<TG>, ...)."
                 << endl << "   Length of grid #" << i << " along dimension #"
                 << j << " is " << grids_(i)->GetLength(j) << " but should be "
                 << grids_(j)->GetLength() << "." << endl;
            abort();
          }

    for (i = 0; i < N; i++)
      if (grids_(i)->GetLength(i) != data_.extent(i))
        throw WrongDim("Data<T, " + to_str(N) + ">::ResizeGrid(Grid<TG>, ...)",
                       "Length of grid #" + to_str(i) + " is "
                       + to_str(grids_(i)->GetLength(i)) + " but should be "
                       + to_str(data_.extent(i)) + ".");
#endif

  }

  //! Resizes the grids.
  /*!
    \param G0 grid for dimension #0.
    \param G1 grid for dimension #1.
    \param G2 grid for dimension #2.
    \param G3 grid for dimension #3.
    \param G4 grid for dimension #4.
  */
  template<class T, int N, class TG>
  void Data<T, N, TG>::ResizeGrid(Grid<TG>& G0, Grid<TG>& G1, Grid<TG>& G2,
                                  Grid<TG>& G3, Grid<TG>& G4)
  {

#ifdef SELDONDATA_DEBUG_CHECK_DIMENSIONS
    if (N != 5)
      {
        cout << "ERROR!" << endl;
        cout << "Wrong dimension in Data<T, N, TG>::ResizeGrid(Grid<TG>, ...)."
             << endl << "   Required " + to_str(N)
          + " grids, but got 5 grids." << endl;
        abort();
      }
#endif

    this->ClearGrids();
    grids_(0) = G0.Copy();
    grids_(1) = G1.Copy();
    grids_(2) = G2.Copy();
    grids_(3) = G3.Copy();
    grids_(4) = G4.Copy();
    SetVariables();


#ifdef SELDONDATA_DEBUG_CHECK_DIMENSIONS
    int i, j;
    for (i = 0; i < N; i++)
      for (j = 0; j < N; j++)
        if ((grids_(i)->GetLength(j) != 0)
            && (grids_(i)->GetLength(j) != grids_(j)->GetLength()))
          {
            cout << "ERROR!" << endl;
            cout << "Wrong dimension in Data<T, N, TG>::ResizeGrid(Grid<TG>, ...)."
                 << endl << "   Length of grid #" << i << " along dimension #"
                 << j << " is " << grids_(i)->GetLength(j) << " but should be "
                 << grids_(j)->GetLength() << "." << endl;
            abort();
          }

    for (i = 0; i < N; i++)
      if (grids_(i)->GetLength(i) != data_.extent(i))
        throw WrongDim("Data<T, " + to_str(N) + ">::ResizeGrid(Grid<TG>, ...)",
                       "Length of grid #" + to_str(i) + " is "
                       + to_str(grids_(i)->GetLength(i)) + " but should be "
                       + to_str(data_.extent(i)) + ".");
#endif

  }

  //! Resizes the grids.
  /*!
    \param G0 grid for dimension #0.
    \param G1 grid for dimension #1.
    \param G2 grid for dimension #2.
    \param G3 grid for dimension #3.
    \param G4 grid for dimension #4.
    \param G5 grid for dimension #5.
  */
  template<class T, int N, class TG>
  void Data<T, N, TG>::ResizeGrid(Grid<TG>& G0, Grid<TG>& G1, Grid<TG>& G2,
                                  Grid<TG>& G3, Grid<TG>& G4, Grid<TG>& G5)
  {

#ifdef SELDONDATA_DEBUG_CHECK_DIMENSIONS
    if (N != 6)
      {
        cout << "ERROR!" << endl;
        cout << "Wrong dimension in Data<T, N, TG>::ResizeGrid(Grid<TG>, ...)."
             << endl << "   Required " + to_str(N)
          + " grids, but got 6 grids." << endl;
        abort();
      }
#endif

    this->ClearGrids();
    grids_(0) = G0.Copy();
    grids_(1) = G1.Copy();
    grids_(2) = G2.Copy();
    grids_(3) = G3.Copy();
    grids_(4) = G4.Copy();
    grids_(5) = G5.Copy();
    SetVariables();


#ifdef SELDONDATA_DEBUG_CHECK_DIMENSIONS
    int i, j;
    for (i = 0; i < N; i++)
      for (j = 0; j < N; j++)
        if ((grids_(i)->GetLength(j) != 0)
            && (grids_(i)->GetLength(j) != grids_(j)->GetLength()))
          {
            cout << "ERROR!" << endl;
            cout << "Wrong dimension in Data<T, N, TG>::ResizeGrid(Grid<TG>, ...)."
                 << endl << "   Length of grid #" << i << " along dimension #"
                 << j << " is " << grids_(i)->GetLength(j) << " but should be "
                 << grids_(j)->GetLength() << "." << endl;
            abort();
          }

    for (i = 0; i < N; i++)
      if (grids_(i)->GetLength(i) != data_.extent(i))
        throw WrongDim("Data<T, " + to_str(N) + ">::ResizeGrid(Grid<TG>, ...)",
                       "Length of grid #" + to_str(i) + " is "
                       + to_str(grids_(i)->GetLength(i)) + " but should be "
                       + to_str(data_.extent(i)) + ".");
#endif

  }

  //! Resizes the grids.
  /*!
    \param G0 grid for dimension #0.
    \param G1 grid for dimension #1.
    \param G2 grid for dimension #2.
    \param G3 grid for dimension #3.
    \param G4 grid for dimension #4.
    \param G5 grid for dimension #5.
    \param G6 grid for dimension #6.
  */
  template<class T, int N, class TG>
  void Data<T, N, TG>::ResizeGrid(Grid<TG>& G0, Grid<TG>& G1, Grid<TG>& G2,
                                  Grid<TG>& G3, Grid<TG>& G4, Grid<TG>& G5,
                                  Grid<TG>& G6)
  {

#ifdef SELDONDATA_DEBUG_CHECK_DIMENSIONS
    if (N != 7)
      {
        cout << "ERROR!" << endl;
        cout << "Wrong dimension in Data<T, N, TG>::ResizeGrid(Grid<TG>, ...)."
             << endl << "   Required " + to_str(N)
          + " grids, but got 7 grids." << endl;
        abort();
      }
#endif

    this->ClearGrids();
    grids_(0) = G0.Copy();
    grids_(1) = G1.Copy();
    grids_(2) = G2.Copy();
    grids_(3) = G3.Copy();
    grids_(4) = G4.Copy();
    grids_(5) = G5.Copy();
    grids_(6) = G6.Copy();
    SetVariables();


#ifdef SELDONDATA_DEBUG_CHECK_DIMENSIONS
    int i, j;
    for (i = 0; i < N; i++)
      for (j = 0; j < N; j++)
        if ((grids_(i)->GetLength(j) != 0)
            && (grids_(i)->GetLength(j) != grids_(j)->GetLength()))
          {
            cout << "ERROR!" << endl;
            cout << "Wrong dimension in Data<T, N, TG>::ResizeGrid(Grid<TG>, ...)."
                 << endl << "   Length of grid #" << i << " along dimension #"
                 << j << " is " << grids_(i)->GetLength(j) << " but should be "
                 << grids_(j)->GetLength() << "." << endl;
            abort();
          }

    for (i = 0; i < N; i++)
      if (grids_(i)->GetLength(i) != data_.extent(i))
        throw WrongDim("Data<T, " + to_str(N) + ">::ResizeGrid(Grid<TG>, ...)",
                       "Length of grid #" + to_str(i) + " is "
                       + to_str(grids_(i)->GetLength(i)) + " but should be "
                       + to_str(data_.extent(i)) + ".");
#endif

  }

  //! Resizes the grids.
  /*!
    \param G0 grid for dimension #0.
    \param G1 grid for dimension #1.
    \param G2 grid for dimension #2.
    \param G3 grid for dimension #3.
    \param G4 grid for dimension #4.
    \param G5 grid for dimension #5.
    \param G6 grid for dimension #6.
    \param G7 grid for dimension #7.
  */
  template<class T, int N, class TG>
  void Data<T, N, TG>::ResizeGrid(Grid<TG>& G0, Grid<TG>& G1, Grid<TG>& G2,
                                  Grid<TG>& G3, Grid<TG>& G4, Grid<TG>& G5,
                                  Grid<TG>& G6, Grid<TG>& G7)
  {

#ifdef SELDONDATA_DEBUG_CHECK_DIMENSIONS
    if (N != 8)
      {
        cout << "ERROR!" << endl;
        cout << "Wrong dimension in Data<T, N, TG>::ResizeGrid(Grid<TG>, ...)."
             << endl << "   Required " + to_str(N)
          + " grids, but got 8 grids." << endl;
        abort();
      }
#endif

    this->ClearGrids();
    grids_(0) = G0.Copy();
    grids_(1) = G1.Copy();
    grids_(2) = G2.Copy();
    grids_(3) = G3.Copy();
    grids_(4) = G4.Copy();
    grids_(5) = G5.Copy();
    grids_(6) = G6.Copy();
    grids_(7) = G7.Copy();
    SetVariables();


#ifdef SELDONDATA_DEBUG_CHECK_DIMENSIONS
    int i, j;
    for (i = 0; i < N; i++)
      for (j = 0; j < N; j++)
        if ((grids_(i)->GetLength(j) != 0)
            && (grids_(i)->GetLength(j) != grids_(j)->GetLength()))
          {
            cout << "ERROR!" << endl;
            cout << "Wrong dimension in Data<T, N, TG>::ResizeGrid(Grid<TG>, ...)."
                 << endl << "   Length of grid #" << i << " along dimension #"
                 << j << " is " << grids_(i)->GetLength(j) << " but should be "
                 << grids_(j)->GetLength() << "." << endl;
            abort();
          }

    for (i = 0; i < N; i++)
      if (grids_(i)->GetLength(i) != data_.extent(i))
        throw WrongDim("Data<T, " + to_str(N) + ">::ResizeGrid(Grid<TG>, ...)",
                       "Length of grid #" + to_str(i) + " is "
                       + to_str(grids_(i)->GetLength(i)) + " but should be "
                       + to_str(data_.extent(i)) + ".");
#endif

  }

  //! Resizes the grids.
  /*!
    \param G0 grid for dimension #0.
    \param G1 grid for dimension #1.
    \param G2 grid for dimension #2.
    \param G3 grid for dimension #3.
    \param G4 grid for dimension #4.
    \param G5 grid for dimension #5.
    \param G6 grid for dimension #6.
    \param G7 grid for dimension #7.
    \param G8 grid for dimension #8.
  */
  template<class T, int N, class TG>
  void Data<T, N, TG>::ResizeGrid(Grid<TG>& G0, Grid<TG>& G1, Grid<TG>& G2,
                                  Grid<TG>& G3, Grid<TG>& G4, Grid<TG>& G5,
                                  Grid<TG>& G6, Grid<TG>& G7, Grid<TG>& G8)
  {

#ifdef SELDONDATA_DEBUG_CHECK_DIMENSIONS
    if (N != 9)
      {
        cout << "ERROR!" << endl;
        cout << "Wrong dimension in Data<T, N, TG>::ResizeGrid(Grid<TG>, ...)."
             << endl << "   Required " + to_str(N)
          + " grids, but got 9 grids." << endl;
        abort();
      }
#endif

    this->ClearGrids();
    grids_(0) = G0.Copy();
    grids_(1) = G1.Copy();
    grids_(2) = G2.Copy();
    grids_(3) = G3.Copy();
    grids_(4) = G4.Copy();
    grids_(5) = G5.Copy();
    grids_(6) = G6.Copy();
    grids_(7) = G7.Copy();
    grids_(8) = G8.Copy();
    SetVariables();


#ifdef SELDONDATA_DEBUG_CHECK_DIMENSIONS
    int i, j;
    for (i = 0; i < N; i++)
      for (j = 0; j < N; j++)
        if ((grids_(i)->GetLength(j) != 0)
            && (grids_(i)->GetLength(j) != grids_(j)->GetLength()))
          {
            cout << "ERROR!" << endl;
            cout << "Wrong dimension in Data<T, N, TG>::ResizeGrid(Grid<TG>, ...)."
                 << endl << "   Length of grid #" << i << " along dimension #"
                 << j << " is " << grids_(i)->GetLength(j) << " but should be "
                 << grids_(j)->GetLength() << "." << endl;
            abort();
          }

    for (i = 0; i < N; i++)
      if (grids_(i)->GetLength(i) != data_.extent(i))
        throw WrongDim("Data<T, " + to_str(N) + ">::ResizeGrid(Grid<TG>, ...)",
                       "Length of grid #" + to_str(i) + " is "
                       + to_str(grids_(i)->GetLength(i)) + " but should be "
                       + to_str(data_.extent(i)) + ".");
#endif

  }

  //! Resizes the grids.
  /*!
    \param G0 grid for dimension #0.
    \param G1 grid for dimension #1.
    \param G2 grid for dimension #2.
    \param G3 grid for dimension #3.
    \param G4 grid for dimension #4.
    \param G5 grid for dimension #5.
    \param G6 grid for dimension #6.
    \param G7 grid for dimension #7.
    \param G8 grid for dimension #8.
    \param G9 grid for dimension #9.
  */
  template<class T, int N, class TG>
  void Data<T, N, TG>::ResizeGrid(Grid<TG>& G0, Grid<TG>& G1, Grid<TG>& G2,
                                  Grid<TG>& G3, Grid<TG>& G4, Grid<TG>& G5,
                                  Grid<TG>& G6, Grid<TG>& G7, Grid<TG>& G8,
                                  Grid<TG>& G9)
  {

#ifdef SELDONDATA_DEBUG_CHECK_DIMENSIONS
    if (N != 10)
      {
        cout << "ERROR!" << endl;
        cout << "Wrong dimension in Data<T, N, TG>::ResizeGrid(Grid<TG>, ...)."
             << endl << "   Required " + to_str(N)
          + " grids, but got 10 grids." << endl;
        abort();
      }
#endif

    this->ClearGrids();
    grids_(0) = G0.Copy();
    grids_(1) = G1.Copy();
    grids_(2) = G2.Copy();
    grids_(3) = G3.Copy();
    grids_(4) = G4.Copy();
    grids_(5) = G5.Copy();
    grids_(6) = G6.Copy();
    grids_(7) = G7.Copy();
    grids_(8) = G8.Copy();
    grids_(9) = G9.Copy();
    SetVariables();


#ifdef SELDONDATA_DEBUG_CHECK_DIMENSIONS
    int i, j;
    for (i = 0; i < N; i++)
      for (j = 0; j < N; j++)
        if ((grids_(i)->GetLength(j) != 0)
            && (grids_(i)->GetLength(j) != grids_(j)->GetLength()))
          {
            cout << "ERROR!" << endl;
            cout << "Wrong dimension in Data<T, N, TG>::ResizeGrid(Grid<TG>, ...)."
                 << endl << "   Length of grid #" << i << " along dimension #"
                 << j << " is " << grids_(i)->GetLength(j) << " but should be "
                 << grids_(j)->GetLength() << "." << endl;
            abort();
          }

    for (i = 0; i < N; i++)
      if (grids_(i)->GetLength(i) != data_.extent(i))
        throw WrongDim("Data<T, " + to_str(N) + ">::ResizeGrid(Grid<TG>, ...)",
                       "Length of grid #" + to_str(i) + " is "
                       + to_str(grids_(i)->GetLength(i)) + " but should be "
                       + to_str(data_.extent(i)) + ".");
#endif

  }

  //! Resizes the grid.
  /*!
    \param N0 grid length along dimension #0.
  */
  template<class T, int N, class TG>
  void Data<T, N, TG>::ResizeGrid(int N0)
  {

    RegularGrid<TG> G0(N0);
    this->ResizeGrid(G0);

  }

  //! Resizes the grids.
  /*!
    \param N0 grid length along dimension #0.
    \param N1 grid length along dimension #1.
  */
  template<class T, int N, class TG>
  void Data<T, N, TG>::ResizeGrid(int N0, int N1)
  {

    RegularGrid<TG> G0(N0);
    RegularGrid<TG> G1(N1);
    this->ResizeGrid(G0, G1);

  }

  //! Resizes the grids.
  /*!
    \param N0 grid length along dimension #0.
    \param N1 grid length along dimension #1.
    \param N2 grid length along dimension #2.
  */
  template<class T, int N, class TG>
  void Data<T, N, TG>::ResizeGrid(int N0, int N1, int N2)
  {

    RegularGrid<TG> G0(N0);
    RegularGrid<TG> G1(N1);
    RegularGrid<TG> G2(N2);
    this->ResizeGrid(G0, G1, G2);

  }

  //! Resizes the grids.
  /*!
    \param N0 grid length along dimension #0.
    \param N1 grid length along dimension #1.
    \param N2 grid length along dimension #2.
    \param N3 grid length along dimension #3.
  */
  template<class T, int N, class TG>
  void Data<T, N, TG>::ResizeGrid(int N0, int N1, int N2,
                                  int N3)
  {

    RegularGrid<TG> G0(N0);
    RegularGrid<TG> G1(N1);
    RegularGrid<TG> G2(N2);
    RegularGrid<TG> G3(N3);
    this->ResizeGrid(G0, G1, G2, G3);

  }

  //! Resizes the grids.
  /*!
    \param N0 grid length along dimension #0.
    \param N1 grid length along dimension #1.
    \param N2 grid length along dimension #2.
    \param N3 grid length along dimension #3.
    \param N4 grid length along dimension #4.
  */
  template<class T, int N, class TG>
  void Data<T, N, TG>::ResizeGrid(int N0, int N1, int N2,
                                  int N3, int N4)
  {

    RegularGrid<TG> G0(N0);
    RegularGrid<TG> G1(N1);
    RegularGrid<TG> G2(N2);
    RegularGrid<TG> G3(N3);
    RegularGrid<TG> G4(N4);
    this->ResizeGrid(G0, G1, G2, G3, G4);

  }

  //! Resizes the grids.
  /*!
    \param N0 grid length along dimension #0.
    \param N1 grid length along dimension #1.
    \param N2 grid length along dimension #2.
    \param N3 grid length along dimension #3.
    \param N4 grid length along dimension #4.
    \param N5 grid length along dimension #5.
  */
  template<class T, int N, class TG>
  void Data<T, N, TG>::ResizeGrid(int N0, int N1, int N2,
                                  int N3, int N4, int N5)
  {

    RegularGrid<TG> G0(N0);
    RegularGrid<TG> G1(N1);
    RegularGrid<TG> G2(N2);
    RegularGrid<TG> G3(N3);
    RegularGrid<TG> G4(N4);
    RegularGrid<TG> G5(N5);
    this->ResizeGrid(G0, G1, G2, G3, G4, G5);

  }

  //! Resizes the grids.
  /*!
    \param N0 grid length along dimension #0.
    \param N1 grid length along dimension #1.
    \param N2 grid length along dimension #2.
    \param N3 grid length along dimension #3.
    \param N4 grid length along dimension #4.
    \param N5 grid length along dimension #5.
    \param N6 grid length along dimension #6.
  */
  template<class T, int N, class TG>
  void Data<T, N, TG>::ResizeGrid(int N0, int N1, int N2,
                                  int N3, int N4, int N5,
                                  int N6)
  {

    RegularGrid<TG> G0(N0);
    RegularGrid<TG> G1(N1);
    RegularGrid<TG> G2(N2);
    RegularGrid<TG> G3(N3);
    RegularGrid<TG> G4(N4);
    RegularGrid<TG> G5(N5);
    RegularGrid<TG> G6(N6);
    this->ResizeGrid(G0, G1, G2, G3, G4, G5, G6);

  }

  //! Resizes the grids.
  /*!
    \param N0 grid length along dimension #0.
    \param N1 grid length along dimension #1.
    \param N2 grid length along dimension #2.
    \param N3 grid length along dimension #3.
    \param N4 grid length along dimension #4.
    \param N5 grid length along dimension #5.
    \param N6 grid length along dimension #6.
    \param N7 grid length along dimension #7.
  */
  template<class T, int N, class TG>
  void Data<T, N, TG>::ResizeGrid(int N0, int N1, int N2,
                                  int N3, int N4, int N5,
                                  int N6, int N7)
  {

    RegularGrid<TG> G0(N0);
    RegularGrid<TG> G1(N1);
    RegularGrid<TG> G2(N2);
    RegularGrid<TG> G3(N3);
    RegularGrid<TG> G4(N4);
    RegularGrid<TG> G5(N5);
    RegularGrid<TG> G6(N6);
    RegularGrid<TG> G7(N7);
    this->ResizeGrid(G0, G1, G2, G3, G4, G5, G6, G7);

  }

  //! Resizes the grids.
  /*!
    \param N0 grid length along dimension #0.
    \param N1 grid length along dimension #1.
    \param N2 grid length along dimension #2.
    \param N3 grid length along dimension #3.
    \param N4 grid length along dimension #4.
    \param N5 grid length along dimension #5.
    \param N6 grid length along dimension #6.
    \param N7 grid length along dimension #7.
    \param N8 grid length along dimension #8.
  */
  template<class T, int N, class TG>
  void Data<T, N, TG>::ResizeGrid(int N0, int N1, int N2,
                                  int N3, int N4, int N5,
                                  int N6, int N7, int N8)
  {

    RegularGrid<TG> G0(N0);
    RegularGrid<TG> G1(N1);
    RegularGrid<TG> G2(N2);
    RegularGrid<TG> G3(N3);
    RegularGrid<TG> G4(N4);
    RegularGrid<TG> G5(N5);
    RegularGrid<TG> G6(N6);
    RegularGrid<TG> G7(N7);
    RegularGrid<TG> G8(N8);
    this->ResizeGrid(G0, G1, G2, G3, G4, G5, G6, G7, G8);

  }

  //! Resizes the grids.
  /*!
    \param N0 grid length along dimension #0.
    \param N1 grid length along dimension #1.
    \param N2 grid length along dimension #2.
    \param N3 grid length along dimension #3.
    \param N4 grid length along dimension #4.
    \param N5 grid length along dimension #5.
    \param N6 grid length along dimension #6.
    \param N7 grid length along dimension #7.
    \param N8 grid length along dimension #8.
    \param N9 grid length along dimension #9.
  */
  template<class T, int N, class TG>
  void Data<T, N, TG>::ResizeGrid(int N0, int N1, int N2,
                                  int N3, int N4, int N5,
                                  int N6, int N7, int N8,
                                  int N9)
  {

    RegularGrid<TG> G0(N0);
    RegularGrid<TG> G1(N1);
    RegularGrid<TG> G2(N2);
    RegularGrid<TG> G3(N3);
    RegularGrid<TG> G4(N4);
    RegularGrid<TG> G5(N5);
    RegularGrid<TG> G6(N6);
    RegularGrid<TG> G7(N7);
    RegularGrid<TG> G8(N8);
    RegularGrid<TG> G9(N9);
    this->ResizeGrid(G0, G1, G2, G3, G4, G5, G6, G7, G8, G9);

  }

  //! Clears the data.
  template<class T, int N, class TG>
  void Data<T, N, TG>::Resize()
  {
    this->ClearGrids();
    data_.free();
  }

  //! Resizes an one-dimensional data.
  /*!
    \param G0 grid for dimension #0.
  */
  template<class T, int N, class TG>
  void Data<T, N, TG>::Resize(Grid<TG>& G0)
  {

#ifdef SELDONDATA_DEBUG_CHECK_DIMENSIONS
    if (N != 1)
      {
        cout << "ERROR!" << endl;
        cout << "Wrong dimension in Data<T, N, TG>::Resize(Grid<TG>, ...)."
             << endl << "   Required " + to_str(N)
          + " grids, but got 1 grid." << endl;
        abort();
      }
#endif

    this->ClearGrids();
    grids_(0) = G0.Copy();
    SetVariables();

#ifdef SELDONDATA_DEBUG_CHECK_DIMENSIONS
    int i, j;
    for (i = 0; i < N; i++)
      for (j = 0; j < N; j++)
        if ((grids_(i)->GetLength(j) != 0)
            && (grids_(i)->GetLength(j) != grids_(j)->GetLength()))
          {
            cout << "ERROR!" << endl;
            cout << "Wrong dimension in Data<T, N, TG>::Resize(Grid<TG>, ...)."
                 << endl << "   Length of grid #" << i << " along dimension #"
                 << j << " is " << grids_(i)->GetLength(j) << " but should be "
                 << grids_(j)->GetLength() << "." << endl;
            abort();
          }
#endif

    data_.resize(G0.GetLength());
    InitData();
  }

  //! Resizes a two-dimensional data.
  /*!
    \param G0 grid for dimension #0.
    \param G1 grid for dimension #1.
  */
  template<class T, int N, class TG>
  void Data<T, N, TG>::Resize(Grid<TG>& G0, Grid<TG>& G1)
  {

#ifdef SELDONDATA_DEBUG_CHECK_DIMENSIONS
    if (N != 2)
      {
        cout << "ERROR!" << endl;
        cout << "Wrong dimension in Data<T, N, TG>::Resize(Grid<TG>, ...)."
             << endl << "   Required " + to_str(N)
          + " grids, but got 2 grids." << endl;
        abort();
      }
#endif

    this->ClearGrids();
    grids_(0) = G0.Copy();
    grids_(1) = G1.Copy();
    SetVariables();

#ifdef SELDONDATA_DEBUG_CHECK_DIMENSIONS
    int i, j;
    for (i = 0; i < N; i++)
      for (j = 0; j < N; j++)
        if ((grids_(i)->GetLength(j) != 0)
            && (grids_(i)->GetLength(j) != grids_(j)->GetLength()))
          {
            cout << "ERROR!" << endl;
            cout << "Wrong dimension in Data<T, N, TG>::Resize(Grid<TG>, ...)."
                 << endl << "   Length of grid #" << i << " along dimension #"
                 << j << " is " << grids_(i)->GetLength(j) << " but should be "
                 << grids_(j)->GetLength() << "." << endl;
            abort();
          }
#endif

    data_.resize(G0.GetLength(), G1.GetLength());
    InitData();

  }

  //! Resizes a three-dimensional data.
  /*!
    \param G0 grid for dimension #0.
    \param G1 grid for dimension #1.
    \param G2 grid for dimension #2.
  */
  template<class T, int N, class TG>
  void Data<T, N, TG>::Resize(Grid<TG>& G0, Grid<TG>& G1, Grid<TG>& G2)
  {

#ifdef SELDONDATA_DEBUG_CHECK_DIMENSIONS
    if (N != 3)
      {
        cout << "ERROR!" << endl;
        cout << "Wrong dimension in Data<T, N, TG>::Resize(Grid<TG>, ...)."
             << endl << "   Required " + to_str(N)
          + " grids, but got 3 grids." << endl;
        abort();
      }
#endif

    this->ClearGrids();
    grids_(0) = G0.Copy();
    grids_(1) = G1.Copy();
    grids_(2) = G2.Copy();
    SetVariables();


#ifdef SELDONDATA_DEBUG_CHECK_DIMENSIONS
    int i, j;
    for (i = 0; i < N; i++)
      for (j = 0; j < N; j++)
        if ((grids_(i)->GetLength(j) != 0)
            && (grids_(i)->GetLength(j) != grids_(j)->GetLength()))
          {
            cout << "ERROR!" << endl;
            cout << "Wrong dimension in Data<T, N, TG>::Resize(Grid<TG>, ...)."
                 << endl << "   Length of grid #" << i << " along dimension #"
                 << j << " is " << grids_(i)->GetLength(j) << " but should be "
                 << grids_(j)->GetLength() << "." << endl;
            abort();
          }
#endif

    data_.resize(G0.GetLength(), G1.GetLength(),
                 G2.GetLength());
    InitData();

  }

  //! Resizes a four-dimensional data.
  /*!
    \param G0 grid for dimension #0.
    \param G1 grid for dimension #1.
    \param G2 grid for dimension #2.
    \param G3 grid for dimension #3.
  */
  template<class T, int N, class TG>
  void Data<T, N, TG>::Resize(Grid<TG>& G0, Grid<TG>& G1, Grid<TG>& G2,
                              Grid<TG>& G3)
  {

#ifdef SELDONDATA_DEBUG_CHECK_DIMENSIONS
    if (N != 4)
      {
        cout << "ERROR!" << endl;
        cout << "Wrong dimension in Data<T, N, TG>::Resize(Grid<TG>, ...)."
             << endl << "   Required " + to_str(N)
          + " grids, but got 4 grids." << endl;
        abort();
      }
#endif

    this->ClearGrids();
    grids_(0) = G0.Copy();
    grids_(1) = G1.Copy();
    grids_(2) = G2.Copy();
    grids_(3) = G3.Copy();
    SetVariables();


#ifdef SELDONDATA_DEBUG_CHECK_DIMENSIONS
    int i, j;
    for (i = 0; i < N; i++)
      for (j = 0; j < N; j++)
        if ((grids_(i)->GetLength(j) != 0)
            && (grids_(i)->GetLength(j) != grids_(j)->GetLength()))
          {
            cout << "ERROR!" << endl;
            cout << "Wrong dimension in Data<T, N, TG>::Resize(Grid<TG>, ...)."
                 << endl << "   Length of grid #" << i << " along dimension #"
                 << j << " is " << grids_(i)->GetLength(j) << " but should be "
                 << grids_(j)->GetLength() << "." << endl;
            abort();
          }
#endif

    data_.resize(G0.GetLength(), G1.GetLength(),
                 G2.GetLength(), G3.GetLength());
    InitData();

  }

  //! Resizes a five-dimensional data.
  /*!
    \param G0 grid for dimension #0.
    \param G1 grid for dimension #1.
    \param G2 grid for dimension #2.
    \param G3 grid for dimension #3.
    \param G4 grid for dimension #4.
  */
  template<class T, int N, class TG>
  void Data<T, N, TG>::Resize(Grid<TG>& G0, Grid<TG>& G1, Grid<TG>& G2,
                              Grid<TG>& G3, Grid<TG>& G4)
  {

#ifdef SELDONDATA_DEBUG_CHECK_DIMENSIONS
    if (N != 5)
      {
        cout << "ERROR!" << endl;
        cout << "Wrong dimension in Data<T, N, TG>::Resize(Grid<TG>, ...)."
             << endl << "   Required " + to_str(N)
          + " grids, but got 5 grids." << endl;
        abort();
      }
#endif

    this->ClearGrids();
    grids_(0) = G0.Copy();
    grids_(1) = G1.Copy();
    grids_(2) = G2.Copy();
    grids_(3) = G3.Copy();
    grids_(4) = G4.Copy();
    SetVariables();


#ifdef SELDONDATA_DEBUG_CHECK_DIMENSIONS
    int i, j;
    for (i = 0; i < N; i++)
      for (j = 0; j < N; j++)
        if ((grids_(i)->GetLength(j) != 0)
            && (grids_(i)->GetLength(j) != grids_(j)->GetLength()))
          {
            cout << "ERROR!" << endl;
            cout << "Wrong dimension in Data<T, N, TG>::Resize(Grid<TG>, ...)."
                 << endl << "   Length of grid #" << i << " along dimension #"
                 << j << " is " << grids_(i)->GetLength(j) << " but should be "
                 << grids_(j)->GetLength() << "." << endl;
            abort();
          }
#endif

    data_.resize(G0.GetLength(), G1.GetLength(),
                 G2.GetLength(), G3.GetLength(),
                 G4.GetLength());
    InitData();

  }

  //! Resizes a six-dimensional data.
  /*!
    \param G0 grid for dimension #0.
    \param G1 grid for dimension #1.
    \param G2 grid for dimension #2.
    \param G3 grid for dimension #3.
    \param G4 grid for dimension #4.
    \param G5 grid for dimension #5.
  */
  template<class T, int N, class TG>
  void Data<T, N, TG>::Resize(Grid<TG>& G0, Grid<TG>& G1, Grid<TG>& G2,
                              Grid<TG>& G3, Grid<TG>& G4, Grid<TG>& G5)
  {

#ifdef SELDONDATA_DEBUG_CHECK_DIMENSIONS
    if (N != 6)
      {
        cout << "ERROR!" << endl;
        cout << "Wrong dimension in Data<T, N, TG>::Resize(Grid<TG>, ...)."
             << endl << "   Required " + to_str(N)
          + " grids, but got 6 grids." << endl;
        abort();
      }
#endif

    this->ClearGrids();
    grids_(0) = G0.Copy();
    grids_(1) = G1.Copy();
    grids_(2) = G2.Copy();
    grids_(3) = G3.Copy();
    grids_(4) = G4.Copy();
    grids_(5) = G5.Copy();
    SetVariables();


#ifdef SELDONDATA_DEBUG_CHECK_DIMENSIONS
    int i, j;
    for (i = 0; i < N; i++)
      for (j = 0; j < N; j++)
        if ((grids_(i)->GetLength(j) != 0)
            && (grids_(i)->GetLength(j) != grids_(j)->GetLength()))
          {
            cout << "ERROR!" << endl;
            cout << "Wrong dimension in Data<T, N, TG>::Resize(Grid<TG>, ...)."
                 << endl << "   Length of grid #" << i << " along dimension #"
                 << j << " is " << grids_(i)->GetLength(j) << " but should be "
                 << grids_(j)->GetLength() << "." << endl;
            abort();
          }
#endif

    data_.resize(G0.GetLength(), G1.GetLength(),
                 G2.GetLength(), G3.GetLength(), G4.GetLength(),
                 G5.GetLength());
    InitData();

  }

  //! Resizes a seven-dimensional data.
  /*!
    \param G0 grid for dimension #0.
    \param G1 grid for dimension #1.
    \param G2 grid for dimension #2.
    \param G3 grid for dimension #3.
    \param G4 grid for dimension #4.
    \param G5 grid for dimension #5.
    \param G6 grid for dimension #6.
  */
  template<class T, int N, class TG>
  void Data<T, N, TG>::Resize(Grid<TG>& G0, Grid<TG>& G1, Grid<TG>& G2,
                              Grid<TG>& G3, Grid<TG>& G4, Grid<TG>& G5,
                              Grid<TG>& G6)
  {

#ifdef SELDONDATA_DEBUG_CHECK_DIMENSIONS
    if (N != 7)
      {
        cout << "ERROR!" << endl;
        cout << "Wrong dimension in Data<T, N, TG>::Resize(Grid<TG>, ...)."
             << endl << "   Required " + to_str(N)
          + " grids, but got 7 grids." << endl;
        abort();
      }
#endif

    this->ClearGrids();
    grids_(0) = G0.Copy();
    grids_(1) = G1.Copy();
    grids_(2) = G2.Copy();
    grids_(3) = G3.Copy();
    grids_(4) = G4.Copy();
    grids_(5) = G5.Copy();
    grids_(6) = G6.Copy();
    SetVariables();


#ifdef SELDONDATA_DEBUG_CHECK_DIMENSIONS
    int i, j;
    for (i = 0; i < N; i++)
      for (j = 0; j < N; j++)
        if ((grids_(i)->GetLength(j) != 0)
            && (grids_(i)->GetLength(j) != grids_(j)->GetLength()))
          {
            cout << "ERROR!" << endl;
            cout << "Wrong dimension in Data<T, N, TG>::Resize(Grid<TG>, ...)."
                 << endl << "   Length of grid #" << i << " along dimension #"
                 << j << " is " << grids_(i)->GetLength(j) << " but should be "
                 << grids_(j)->GetLength() << "." << endl;
            abort();
          }
#endif

    data_.resize(G0.GetLength(), G1.GetLength(),
                 G2.GetLength(), G3.GetLength(), G4.GetLength(),
                 G5.GetLength(), G6.GetLength());
    InitData();

  }

  //! Resizes an eight-dimensional data.
  /*!
    \param G0 grid for dimension #0.
    \param G1 grid for dimension #1.
    \param G2 grid for dimension #2.
    \param G3 grid for dimension #3.
    \param G4 grid for dimension #4.
    \param G5 grid for dimension #5.
    \param G6 grid for dimension #6.
    \param G7 grid for dimension #7.
  */
  template<class T, int N, class TG>
  void Data<T, N, TG>::Resize(Grid<TG>& G0, Grid<TG>& G1, Grid<TG>& G2,
                              Grid<TG>& G3, Grid<TG>& G4, Grid<TG>& G5,
                              Grid<TG>& G6, Grid<TG>& G7)
  {

#ifdef SELDONDATA_DEBUG_CHECK_DIMENSIONS
    if (N != 8)
      {
        cout << "ERROR!" << endl;
        cout << "Wrong dimension in Data<T, N, TG>::Resize(Grid<TG>, ...)."
             << endl << "   Required " + to_str(N)
          + " grids, but got 8 grids." << endl;
        abort();
      }
#endif

    this->ClearGrids();
    grids_(0) = G0.Copy();
    grids_(1) = G1.Copy();
    grids_(2) = G2.Copy();
    grids_(3) = G3.Copy();
    grids_(4) = G4.Copy();
    grids_(5) = G5.Copy();
    grids_(6) = G6.Copy();
    grids_(7) = G7.Copy();
    SetVariables();


#ifdef SELDONDATA_DEBUG_CHECK_DIMENSIONS
    int i, j;
    for (i = 0; i < N; i++)
      for (j = 0; j < N; j++)
        if ((grids_(i)->GetLength(j) != 0)
            && (grids_(i)->GetLength(j) != grids_(j)->GetLength()))
          {
            cout << "ERROR!" << endl;
            cout << "Wrong dimension in Data<T, N, TG>::Resize(Grid<TG>, ...)."
                 << endl << "   Length of grid #" << i << " along dimension #"
                 << j << " is " << grids_(i)->GetLength(j) << " but should be "
                 << grids_(j)->GetLength() << "." << endl;
            abort();
          }
#endif

    data_.resize(G0.GetLength(), G1.GetLength(),
                 G2.GetLength(), G3.GetLength(), G4.GetLength(),
                 G5.GetLength(), G6.GetLength(), G7.GetLength());
    InitData();

  }

  //! Resizes a nine-dimensional data.
  /*!
    \param G0 grid for dimension #0.
    \param G1 grid for dimension #1.
    \param G2 grid for dimension #2.
    \param G3 grid for dimension #3.
    \param G4 grid for dimension #4.
    \param G5 grid for dimension #5.
    \param G6 grid for dimension #6.
    \param G7 grid for dimension #7.
    \param G8 grid for dimension #8.
  */
  template<class T, int N, class TG>
  void Data<T, N, TG>::Resize(Grid<TG>& G0, Grid<TG>& G1, Grid<TG>& G2,
                              Grid<TG>& G3, Grid<TG>& G4, Grid<TG>& G5,
                              Grid<TG>& G6, Grid<TG>& G7, Grid<TG>& G8)
  {

#ifdef SELDONDATA_DEBUG_CHECK_DIMENSIONS
    if (N != 9)
      {
        cout << "ERROR!" << endl;
        cout << "Wrong dimension in Data<T, N, TG>::Resize(Grid<TG>, ...)."
             << endl << "   Required " + to_str(N)
          + " grids, but got 9 grids." << endl;
        abort();
      }
#endif

    this->ClearGrids();
    grids_(0) = G0.Copy();
    grids_(1) = G1.Copy();
    grids_(2) = G2.Copy();
    grids_(3) = G3.Copy();
    grids_(4) = G4.Copy();
    grids_(5) = G5.Copy();
    grids_(6) = G6.Copy();
    grids_(7) = G7.Copy();
    grids_(8) = G8.Copy();
    SetVariables();


#ifdef SELDONDATA_DEBUG_CHECK_DIMENSIONS
    int i, j;
    for (i = 0; i < N; i++)
      for (j = 0; j < N; j++)
        if ((grids_(i)->GetLength(j) != 0)
            && (grids_(i)->GetLength(j) != grids_(j)->GetLength()))
          {
            cout << "ERROR!" << endl;
            cout << "Wrong dimension in Data<T, N, TG>::Resize(Grid<TG>, ...)."
                 << endl << "   Length of grid #" << i << " along dimension #"
                 << j << " is " << grids_(i)->GetLength(j) << " but should be "
                 << grids_(j)->GetLength() << "." << endl;
            abort();
          }
#endif

    data_.resize(G0.GetLength(), G1.GetLength(),
                 G2.GetLength(), G3.GetLength(), G4.GetLength(),
                 G5.GetLength(), G6.GetLength(), G7.GetLength(),
                 G8.GetLength());
    InitData();

  }

  //! Resizes a ten-dimensional data.
  /*!
    \param G0 grid for dimension #0.
    \param G1 grid for dimension #1.
    \param G2 grid for dimension #2.
    \param G3 grid for dimension #3.
    \param G4 grid for dimension #4.
    \param G5 grid for dimension #5.
    \param G6 grid for dimension #6.
    \param G7 grid for dimension #7.
    \param G8 grid for dimension #8.
    \param G9 grid for dimension #9.
  */
  template<class T, int N, class TG>
  void Data<T, N, TG>::Resize(Grid<TG>& G0, Grid<TG>& G1, Grid<TG>& G2,
                              Grid<TG>& G3, Grid<TG>& G4, Grid<TG>& G5,
                              Grid<TG>& G6, Grid<TG>& G7, Grid<TG>& G8,
                              Grid<TG>& G9)
  {

#ifdef SELDONDATA_DEBUG_CHECK_DIMENSIONS
    if (N != 10)
      {
        cout << "ERROR!" << endl;
        cout << "Wrong dimension in Data<T, N, TG>::Resize(Grid<TG>, ...)."
             << endl << "   Required " + to_str(N)
          + " grids, but got 10 grids." << endl;
        abort();
      }
#endif

    this->ClearGrids();
    grids_(0) = G0.Copy();
    grids_(1) = G1.Copy();
    grids_(2) = G2.Copy();
    grids_(3) = G3.Copy();
    grids_(4) = G4.Copy();
    grids_(5) = G5.Copy();
    grids_(6) = G6.Copy();
    grids_(7) = G7.Copy();
    grids_(8) = G8.Copy();
    grids_(9) = G9.Copy();
    SetVariables();


#ifdef SELDONDATA_DEBUG_CHECK_DIMENSIONS
    int i, j;
    for (i = 0; i < N; i++)
      for (j = 0; j < N; j++)
        if ((grids_(i)->GetLength(j) != 0)
            && (grids_(i)->GetLength(j) != grids_(j)->GetLength()))
          {
            cout << "ERROR!" << endl;
            cout << "Wrong dimension in Data<T, N, TG>::Resize(Grid<TG>, ...)."
                 << endl << "   Length of grid #" << i << " along dimension #"
                 << j << " is " << grids_(i)->GetLength(j) << " but should be "
                 << grids_(j)->GetLength() << "." << endl;
            abort();
          }
#endif

    data_.resize(G0.GetLength(), G1.GetLength(),
                 G2.GetLength(), G3.GetLength(), G4.GetLength(),
                 G5.GetLength(), G6.GetLength(), G7.GetLength(),
                 G8.GetLength(), G9.GetLength());
    InitData();

  }

  //! Resizes an one-dimensional data.
  /*!
    \param N0 grid length along dimension #0.
  */
  template<class T, int N, class TG>
  void Data<T, N, TG>::Resize(int N0)
  {

    RegularGrid<TG> G0(N0);
    this->Resize(G0);

  }

  //! Resizes a two-dimensional data.
  /*!
    \param N0 grid length along dimension #0.
    \param N1 grid length along dimension #1.
  */
  template<class T, int N, class TG>
  void Data<T, N, TG>::Resize(int N0, int N1)
  {

    RegularGrid<TG> G0(N0);
    RegularGrid<TG> G1(N1);
    this->Resize(G0, G1);

  }

  //! Resizes a three-dimensional data.
  /*!
    \param N0 grid length along dimension #0.
    \param N1 grid length along dimension #1.
    \param N2 grid length along dimension #2.
  */
  template<class T, int N, class TG>
  void Data<T, N, TG>::Resize(int N0, int N1, int N2)
  {

    RegularGrid<TG> G0(N0);
    RegularGrid<TG> G1(N1);
    RegularGrid<TG> G2(N2);
    this->Resize(G0, G1, G2);

  }

  //! Resizes a four-dimensional data.
  /*!
    \param N0 grid length along dimension #0.
    \param N1 grid length along dimension #1.
    \param N2 grid length along dimension #2.
    \param N3 grid length along dimension #3.
  */
  template<class T, int N, class TG>
  void Data<T, N, TG>::Resize(int N0, int N1, int N2,
                              int N3)
  {

    RegularGrid<TG> G0(N0);
    RegularGrid<TG> G1(N1);
    RegularGrid<TG> G2(N2);
    RegularGrid<TG> G3(N3);
    this->Resize(G0, G1, G2, G3);

  }

  //! Resizes a five-dimensional data.
  /*!
    \param N0 grid length along dimension #0.
    \param N1 grid length along dimension #1.
    \param N2 grid length along dimension #2.
    \param N3 grid length along dimension #3.
    \param N4 grid length along dimension #4.
  */
  template<class T, int N, class TG>
  void Data<T, N, TG>::Resize(int N0, int N1, int N2,
                              int N3, int N4)
  {

    RegularGrid<TG> G0(N0);
    RegularGrid<TG> G1(N1);
    RegularGrid<TG> G2(N2);
    RegularGrid<TG> G3(N3);
    RegularGrid<TG> G4(N4);
    this->Resize(G0, G1, G2, G3, G4);

  }

  //! Resizes a six-dimensional data.
  /*!
    \param N0 grid length along dimension #0.
    \param N1 grid length along dimension #1.
    \param N2 grid length along dimension #2.
    \param N3 grid length along dimension #3.
    \param N4 grid length along dimension #4.
    \param N5 grid length along dimension #5.
  */
  template<class T, int N, class TG>
  void Data<T, N, TG>::Resize(int N0, int N1, int N2,
                              int N3, int N4, int N5)
  {

    RegularGrid<TG> G0(N0);
    RegularGrid<TG> G1(N1);
    RegularGrid<TG> G2(N2);
    RegularGrid<TG> G3(N3);
    RegularGrid<TG> G4(N4);
    RegularGrid<TG> G5(N5);
    this->Resize(G0, G1, G2, G3, G4, G5);

  }

  //! Resizes a seven-dimensional data.
  /*!
    \param N0 grid length along dimension #0.
    \param N1 grid length along dimension #1.
    \param N2 grid length along dimension #2.
    \param N3 grid length along dimension #3.
    \param N4 grid length along dimension #4.
    \param N5 grid length along dimension #5.
    \param N6 grid length along dimension #6.
  */
  template<class T, int N, class TG>
  void Data<T, N, TG>::Resize(int N0, int N1, int N2,
                              int N3, int N4, int N5,
                              int N6)
  {

    RegularGrid<TG> G0(N0);
    RegularGrid<TG> G1(N1);
    RegularGrid<TG> G2(N2);
    RegularGrid<TG> G3(N3);
    RegularGrid<TG> G4(N4);
    RegularGrid<TG> G5(N5);
    RegularGrid<TG> G6(N6);
    this->Resize(G0, G1, G2, G3, G4, G5, G6);

  }

  //! Resizes an eight-dimensional data.
  /*!
    \param N0 grid length along dimension #0.
    \param N1 grid length along dimension #1.
    \param N2 grid length along dimension #2.
    \param N3 grid length along dimension #3.
    \param N4 grid length along dimension #4.
    \param N5 grid length along dimension #5.
    \param N6 grid length along dimension #6.
    \param N7 grid length along dimension #7.
  */
  template<class T, int N, class TG>
  void Data<T, N, TG>::Resize(int N0, int N1, int N2,
                              int N3, int N4, int N5,
                              int N6, int N7)
  {

    RegularGrid<TG> G0(N0);
    RegularGrid<TG> G1(N1);
    RegularGrid<TG> G2(N2);
    RegularGrid<TG> G3(N3);
    RegularGrid<TG> G4(N4);
    RegularGrid<TG> G5(N5);
    RegularGrid<TG> G6(N6);
    RegularGrid<TG> G7(N7);
    this->Resize(G0, G1, G2, G3, G4, G5, G6, G7);

  }

  //! Resizes a nine-dimensional data.
  /*!
    \param N0 grid length along dimension #0.
    \param N1 grid length along dimension #1.
    \param N2 grid length along dimension #2.
    \param N3 grid length along dimension #3.
    \param N4 grid length along dimension #4.
    \param N5 grid length along dimension #5.
    \param N6 grid length along dimension #6.
    \param N7 grid length along dimension #7.
    \param N8 grid length along dimension #8.
  */
  template<class T, int N, class TG>
  void Data<T, N, TG>::Resize(int N0, int N1, int N2,
                              int N3, int N4, int N5,
                              int N6, int N7, int N8)
  {

    RegularGrid<TG> G0(N0);
    RegularGrid<TG> G1(N1);
    RegularGrid<TG> G2(N2);
    RegularGrid<TG> G3(N3);
    RegularGrid<TG> G4(N4);
    RegularGrid<TG> G5(N5);
    RegularGrid<TG> G6(N6);
    RegularGrid<TG> G7(N7);
    RegularGrid<TG> G8(N8);
    this->Resize(G0, G1, G2, G3, G4, G5, G6, G7, G8);

  }

  //! Resizes a ten-dimensional data.
  /*!
    \param N0 grid length along dimension #0.
    \param N1 grid length along dimension #1.
    \param N2 grid length along dimension #2.
    \param N3 grid length along dimension #3.
    \param N4 grid length along dimension #4.
    \param N5 grid length along dimension #5.
    \param N6 grid length along dimension #6.
    \param N7 grid length along dimension #7.
    \param N8 grid length along dimension #8.
    \param N9 grid length along dimension #9.
  */
  template<class T, int N, class TG>
  void Data<T, N, TG>::Resize(int N0, int N1, int N2,
                              int N3, int N4, int N5,
                              int N6, int N7, int N8,
                              int N9)
  {

    RegularGrid<TG> G0(N0);
    RegularGrid<TG> G1(N1);
    RegularGrid<TG> G2(N2);
    RegularGrid<TG> G3(N3);
    RegularGrid<TG> G4(N4);
    RegularGrid<TG> G5(N5);
    RegularGrid<TG> G6(N6);
    RegularGrid<TG> G7(N7);
    RegularGrid<TG> G8(N8);
    RegularGrid<TG> G9(N9);
    this->Resize(G0, G1, G2, G3, G4, G5, G6, G7, G8, G9);

  }

  //! Resizes data.
  /*!
    \param shape the new shape.
    \note Data may be lost.
  */
  template<class T, int N, class TG>
  void Data<T, N, TG>::Resize(const TinyVector<int, N>& shape)
  {
    if (N == 1)
      this->Resize(shape(0));
    else if (N == 2)
      this->Resize(shape(0), shape(1));
    else if (N == 3)
      this->Resize(shape(0), shape(1), shape(2));
    else if (N == 4)
      this->Resize(shape(0), shape(1), shape(2),
                   shape(3));
    else if (N == 5)
      this->Resize(shape(0), shape(1), shape(2),
                   shape(3), shape(4));
    else if (N == 6)
      this->Resize(shape(0), shape(1), shape(2),
                   shape(3), shape(4), shape(5));
    else if (N == 7)
      this->Resize(shape(0), shape(1), shape(2),
                   shape(3), shape(4), shape(5),
                   shape(6));
    else if (N == 8)
      this->Resize(shape(0), shape(1), shape(2),
                   shape(3), shape(4), shape(5),
                   shape(6), shape(7));
    else if (N == 9)
      this->Resize(shape(0), shape(1), shape(2),
                   shape(3), shape(4), shape(5),
                   shape(6), shape(7), shape(8));
    else if (N == 10)
      this->Resize(shape(0), shape(1), shape(2),
                   shape(3), shape(4), shape(5),
                   shape(6), shape(7), shape(8),
                   shape(9));
  }

  //! Applies a given function on all elements.
  /*!
    \param function function to be applied.
  */
  template<class T, int N, class TG>
  void Data<T, N, TG>::Apply(void function(T&))
  {
    T* data = data_.data();
    int NbElements = data_.numElements();

    for (int i = 0; i < NbElements; i++)
      function(data[i]);
  }

  //! Multiplies all elements by a given value.
  /*!
    \param alpha number by which all elements are multiplied.
  */
  template<class T, int N, class TG>
  void Data<T, N, TG>::Mlt(T alpha)
  {
    T* data = data_.data();
    int NbElements = data_.numElements();

    for (int i = 0; i < NbElements; i++)
      data[i] *= alpha;
  }

  //! Adds a given value to all elements.
  /*!
    \param alpha number which is added to all elements.
  */
  template<class T, int N, class TG>
  void Data<T, N, TG>::Add(T alpha)
  {
    T* data = data_.data();
    int NbElements = data_.numElements();

    for (int i = 0; i < NbElements; i++)
      data[i] += alpha;
  }

  //! Applies a given function on all elements.
  /*!
    \param function function to be applied.
  */
  template<class T, int N, class TG>
  template<class F>
  void Data<T, N, TG>::Apply(F& function)
  {
    T* data = data_.data();
    int NbElements = data_.numElements();

    for (int i = 0; i < NbElements; i++)
      data[i] = function(data[i]);
  }

  //! Applies a given function on all elements.
  /*!
    \param function a pointer to a function to be applied.
  */
  template<class T, int N, class TG>
  void Data<T, N, TG>::Apply(T(function)(const T&))
  {
    T* data = data_.data();
    int NbElements = data_.numElements();

    for (int i = 0; i < NbElements; i++)
      data[i] = function(data[i]);
  }

  //! Applies a given function on all elements of a 'Data' instance
  //! and put the result in the current instance.
  /*!
    \param D data.
    \param function function to be applied on 'D'.
  */
  template<class T, int N, class TG>
  template<class T0, class TG0, class F>
  void Data<T, N, TG>::Apply(Data<T0, N, TG0>& D, F& function)
  {
    T* data = this->GetData();
    T0* data_in = D.GetData();

    int nb_elements = this->GetNbElements();

#ifdef SELDONDATA_DEBUG_CHECK_DIMENSIONS
    if (nb_elements != D.GetNbElements())
      throw WrongDim("Data<T, " + to_str(N) + ">::Apply(Data<T0, " + to_str(N)
                     + ", TG0>&, F& function)",
                     "Data sizes differ.");
#endif

    for (int i = 0; i < nb_elements; i++)
      data[i] = function(data_in[i]);
  }

  //! Returns the maximum.
  /*!
    \return The maximum.
  */
  template<class T, int N, class TG>
  T Data<T, N, TG>::GetMax() const
  {
    const T* data = data_.data();
    int NbElements = data_.numElements();

    T maximum = data[0];
    for (int i = 1; i < NbElements; i++)
      maximum = max(maximum, data[i]);

    return maximum;
  }

  //! Returns the maximum (absolute value).
  /*!
    \return The maximum (absolute value).
  */
  template<class T, int N, class TG>
  T Data<T, N, TG>::GetMaxAbs() const
  {
    const T* data = data_.data();
    int NbElements = data_.numElements();

    T maximum = abs(data[0]);
    for (int i = 1; i < NbElements; i++)
      maximum = max(maximum, abs(data[i]));

    return maximum;
  }

  //! Returns the (signed) value of the maximum in absolute value.
  /*!
    \return The (signed) value of the maximum in absolute value.
  */
  template<class T, int N, class TG>
  T Data<T, N, TG>::GetSignedMaxAbs() const
  {
    const T* data = data_.data();
    int NbElements = data_.numElements();

    T temp;
    T maximum = abs(data[0]);
    int max_index = 0;
    for (int i = 1; i < NbElements; i++)
      {
        temp = data[i];
        if (maximum < abs(temp))
          {
            maximum = abs(temp);
            max_index = i;
          }
      }

    return data[max_index];
  }

  //! Returns the minimum.
  /*!
    \return The minimum.
  */
  template<class T, int N, class TG>
  T Data<T, N, TG>::GetMin() const
  {
    const T* data = data_.data();
    int NbElements = data_.numElements();

    T minimum = data[0];
    for (int i = 1; i < NbElements; i++)
      minimum = min(minimum, data[i]);

    return minimum;
  }

  //! Returns indices of the maximum.
  /*!
    \return Array of indices of the maximum.
  */
  template<class T, int N, class TG>
  Array<int, 1> Data<T, N, TG>::GetMaxIndex() const
  {
    const T* data = data_.data();
    int NbElements = data_.numElements();

    T maximum = data[0];
    int i_max = 0;
    for (int i = 1; i < NbElements; i++)
      if (maximum < data[i])
        {
          maximum = data[i];
          i_max = i;
        }

    Array<int, 1> Index(N);
    for (int i = N - 1; i >= 0; i--)
      {
        Index(i) = i_max % this->GetLength(i);
        i_max = (i_max - Index(i)) / this->GetLength(i);
      }

    return Index;
  }

  //! Returns indices of the maximum in absolute value.
  /*!
    \return Array of indices of the maximum in absolute value.
  */
  template<class T, int N, class TG>
  Array<int, 1> Data<T, N, TG>::GetMaxAbsIndex() const
  {
    const T* data = data_.data();
    int NbElements = data_.numElements();

    T maximum = abs(data[0]);
    int i_max = 0;
    for (int i = 1; i < NbElements; i++)
      if (maximum < abs(data[i]))
        {
          maximum = abs(data[i]);
          i_max = i;
        }

    Array<int, 1> Index(N);
    for (int i = N - 1; i >= 0; i--)
      {
        Index(i) = i_max % this->GetLength(i);
        i_max = (i_max - Index(i)) / this->GetLength(i);
      }

    return Index;
  }

  //! Returns indices of the minimum.
  /*!
    \return Array of indices of the minimum.
  */
  template<class T, int N, class TG>
  Array<int, 1> Data<T, N, TG>::GetMinIndex() const
  {
    const T* data = data_.data();
    int NbElements = data_.numElements();

    T minimum = data[0];
    int i_min = 0;
    for (int i = 1; i < NbElements; i++)
      if (minimum > data[i])
        {
          minimum = data[i];
          i_min = i;
        }

    Array<int, 1> Index(N);
    for (int i = N - 1; i >= 0; i--)
      {
        Index(i) = i_min % this->GetLength(i);
        i_min = (i_min - Index(i)) / this->GetLength(i);
      }

    return Index;
  }

  //! Returns the sum of all values.
  /*!
    \return The sum of all values.
  */
  template<class T, int N, class TG>
  T Data<T, N, TG>::Sum() const
  {
    const T* data = data_.data();
    int NbElements = data_.numElements();

    T sum = T(0);
    for (int i = 0; i < NbElements; i++)
      sum += data[i];

    return sum;
  }

  //! Computes the sum of all values.
  /*!
    Computes the sum of all values with a type that can be different from the
    type of the values.
  */
  template<class T, int N, class TG>
  template <class Ts>
  void Data<T, N, TG>::Sum(Ts& sum) const
  {
    const T* data = data_.data();
    int NbElements = data_.numElements();

    sum = Ts(0);
    for (int i = 0; i < NbElements; i++)
      sum += data[i];
  }

  //! Returns the mean.
  /*!
    \return The mean.
  */
  template<class T, int N, class TG>
  T Data<T, N, TG>::Mean() const
  {
    return (this->Sum() / data_.numElements());
  }

  //! Computes the mean.
  /*!
    Computes the mean with a type that can be different from the type of
    the values.
  */
  template<class T, int N, class TG>
  template <class Ts>
  void Data<T, N, TG>::Mean(Ts& mean) const
  {
    this->Sum(mean);
    mean = mean / data_.numElements();
  }

  //! Returns the variance.
  /*!
    \return The variance.
  */
  template<class T, int N, class TG>
  T Data<T, N, TG>::Variance() const
  {
    const T* data = data_.data();
    int NbElements = data_.numElements();

    T mean = this->Mean();
    T var = T(0);
    for (int i = 0; i < NbElements; i++)
      var += (data[i] - mean) * (data[i] - mean);

    if (NbElements != 1)
      var = var / T(NbElements - 1);

    return var;
  }

  //! Computes the variance.
  /*!
    Computes the variance with a type that can be different from the type of
    the values.
  */
  template<class T, int N, class TG>
  template <class Ts>
  void Data<T, N, TG>::Variance(Ts& var) const
  {
    const T* data = data_.data();
    int NbElements = data_.numElements();

    Ts mean(0);
    this->Mean(mean);
    var = 0;
    for (int i = 0; i < NbElements; i++)
      var += (data[i] - mean) * (data[i] - mean);

    if (NbElements != 1)
      var = var / T(NbElements - 1);
  }

  //! Returns the standard deviation.
  /*!
    \return The standard deviation.
  */
  template<class T, int N, class TG>
  T Data<T, N, TG>::StandardDeviation() const
  {
    return sqrt(this->Variance());
  }

  //! Computes the standard deviation.
  /*!
    Computes the standard deviation with a type that can be different from
    the type of the values.
  */
  template<class T, int N, class TG>
  template <class Ts>
  void Data<T, N, TG>::StandardDeviation(Ts& std) const
  {
    this->Variance(std);
    std = sqrt(std);
  }

  //! Returns the norm 1.
  /*!
    \return The norm 1.
  */
  template<class T, int N, class TG>
  T Data<T, N, TG>::Norm1() const
  {
    const T* data = data_.data();
    int NbElements = data_.numElements();

    T norm = T(0);
    for (int i = 0; i < NbElements; i++)
      norm += abs(data[i]);

    return norm;
  }

  //! Returns the norm 2.
  /*!
    \return The norm 2.
  */
  template<class T, int N, class TG>
  T Data<T, N, TG>::Norm2() const
  {
    const T* data = data_.data();
    int NbElements = data_.numElements();

    T norm = T(0);
    for (int i = 0; i < NbElements; i++)
      norm += data[i] * data[i];

    return sqrt(norm);
  }

  //! Returns the norm p.
  /*!
    \param p norm parameter.
    \return The norm p.
  */
  template<class T, int N, class TG>
  T Data<T, N, TG>::Norm(T p) const
  {
    const T* data = data_.data();
    int NbElements = data_.numElements();

    T norm = T(0);
    for (int i = 0; i < NbElements; i++)
      norm += pow(data[i], p);

    return pow(norm, 1. / p);
  }

  //! Fills the data array.
  /*!
    Fills the array with integers 0, 1, 2, etc.
  */
  template<class T, int N, class TG>
  void Data<T, N, TG>::Fill()
  {
    T* data = data_.data();
    int NbElements = data_.numElements();

    for (int i = 0; i < NbElements; i++)
      data[i] = T(i);
  }

  //! Sets data to a constant value.
  /*!
    \param value value to which data is set.
  */
  template<class T, int N, class TG>
  void Data<T, N, TG>::Fill(T value)
  {
    T* data = data_.data();
    int NbElements = data_.numElements();

    for (int i = 0; i < NbElements; i++)
      data[i] = value;
  }

  //! Sets data to zero.
  template<class T, int N, class TG>
  void Data<T, N, TG>::SetZero()
  {
    T* data = data_.data();
    int NbElements = data_.numElements();

    for (int i = 0; i < NbElements; i++)
      data[i] = T(0);
  }

  //! Is data equal to zero?
  /*!
    \return true if data is zero, false otherwise.
  */
  template<class T, int N, class TG>
  bool Data<T, N, TG>::IsZero()
  {
    bool res = true;

    T* data = data_.data();
    int NbElements = data_.numElements();

    for (int i = 0; i < NbElements; i++)
      res = res && (data[i] == T(0));

    return res;
  }

  //! Sets data to a signaling NaN.
  template<class T, int N, class TG>
  void Data<T, N, TG>::SetNaN()
  {
    // The following line need to be revised (2017/10/12 YK).
    // data_ = std::numeric_limits<T>::signaling_NaN();
  }

  //! Thresholds data.
  /*! Any value less than 'threshold_min' is set to 'threshold_min'
    and any value more than 'threshold_max' is set to 'threshold_max'.
    \param threshold_min the lower threshold.
    \param threshold_max the upper threshold.
  */
  template<class T, int N, class TG>
  void Data<T, N, TG>::Threshold(T threshold_min, T threshold_max)
  {
    T* data = data_.data();
    int NbElements = data_.numElements();

    for (int i = 0; i < NbElements; i++)
      if (data[i] < threshold_min)
        data[i] = threshold_min;
      else if (data[i] > threshold_max)
        data[i] = threshold_max;
  }

  //! Threshold absolute value of data.
  /*! Any value less than 'threshold' (in absolute value)
    is set to 'threshold'.
    \param threshold the threshold.
  */
  template<class T, int N, class TG>
  void Data<T, N, TG>::ThresholdAbs(T threshold)
  {
    T* data = data_.data();
    int NbElements = data_.numElements();

    for (int i = 0; i < NbElements; i++)
      if (abs(data[i]) < threshold)
        data[i] = ((data[i] > 0) ? T(1.0) : T(-1.0)) * threshold;
  }

  //! Thresholds data.
  /*! Any value less than 'threshold' is set to 'threshold'.
    \param threshold the threshold.
  */
  template<class T, int N, class TG>
  void Data<T, N, TG>::ThresholdMin(T threshold)
  {
    T* data = data_.data();
    int NbElements = data_.numElements();

    for (int i = 0; i < NbElements; i++)
      if (data[i] < threshold)
        data[i] = threshold;
  }

  //! Thresholds data.
  /*! Any value more than 'threshold' is set to 'threshold'.
    \param threshold the threshold.
  */
  template<class T, int N, class TG>
  void Data<T, N, TG>::ThresholdMax(T threshold)
  {
    T* data = data_.data();
    int NbElements = data_.numElements();

    for (int i = 0; i < NbElements; i++)
      if (data[i] > threshold)
        data[i] = threshold;
  }

  //! Computes normalized gross error between two data sets.
  /*!
    Current data is the reference, and the input data is linearly
    interpolated.
    \param data data to be compared with current data.
    \param limit the NGE is computed only with values x such that
    |x| > limit.
    \return The normalized gross error.
  */
  template<class T, int N, class TG>
  template<class T0, class TG0>
  T Data<T, N, TG>::NGE_interpolation(Data<T0, N, TG0>& data, T limit)
  {
    Data<T, N, TG> data0(*this);
    LinearInterpolationGeneral(data, data0);

    return NGE(data0, limit);
  }

  //! Computes normalized gross error between two data sets.
  /*!
    Current data is the reference.
    \param data data to be compared with current data.
    \param limit the NGE is computed only with values x such that
    |x| > limit.
    \return The normalized gross error.
  */
  template<class T, int N, class TG>
  template<class T0, class TG0>
  T Data<T, N, TG>::NGE(Data<T0, N, TG0>& data, T limit)
  {
    T nge;

    T* data_arr = data_.data();
    T0* data0_arr = data.GetData();
    int NbElements = data_.numElements();

#ifdef SELDONDATA_DEBUG_CHECK_DIMENSIONS

    if (NbElements != data.GetArray().numElements())
      throw WrongDim("Data<T, " + to_str(N) + ">::NGE(Data<T, " + to_str(N) + ">&, T)",
                     "Data sizes differ.");

#endif

    int nb_elt = 0;
    nge = T(0);
    for (int i = 0; i < NbElements; i++)
      if (abs(data_arr[i]) > limit)
        {
          nb_elt++;
          nge += abs((data_arr[i] - data0_arr[i]) / data_arr[i]);
        }
    nge = nge / T(nb_elt);

    return nge;
  }

  //! Computes bias between two data sets.
  /*!
    Current data is the reference, and the input data is linearly
    interpolated.
    \param data data to be compared with current data.
    \return The bias.
  */
  template<class T, int N, class TG>
  template<class T0, class TG0>
  T Data<T, N, TG>::Bias_interpolation(Data<T0, N, TG0>& data)
  {
    Data<T, N, TG> data0(*this);
    LinearInterpolationGeneral(data, data0);

    return Bias(data0);
  }

  //! Computes bias between two data sets.
  /*!
    Current data is the reference.
    \param data data to be compared with current data.
    \return The bias.
  */
  template<class T, int N, class TG>
  template<class T0, class TG0>
  T Data<T, N, TG>::Bias(Data<T0, N, TG0>& data)
  {
    T bias = T(0);

    T* data_arr = data_.data();
    T0* data0_arr = data.GetData();
    int NbElements = data_.numElements();

#ifdef SELDONDATA_DEBUG_CHECK_DIMENSIONS

    if (NbElements != data.GetArray().numElements())
      throw WrongDim("Data<T, " + to_str(N) + ">::Bias(Data<T, " + to_str(N) + ">&, T)",
                     "Data sizes differ.");

#endif

    for (int i = 0; i < NbElements; i++)
      bias += data_arr[i] - data0_arr[i];
    bias = bias / T(NbElements);

    return bias;
  }

  //! Computes root mean square between two data sets.
  /*!
    Current data is the reference, and the input data is linearly
    interpolated.
    \param data data to be compared with current data.
    \return The root mean square.
  */
  template<class T, int N, class TG>
  template<class T0, class TG0>
  T Data<T, N, TG>::RMS_interpolation(Data<T0, N, TG0>& data)
  {
    Data<T, N, TG> data0(*this);
    LinearInterpolationGeneral(data, data0);

    return RMS(data0);
  }

  //! Computes root mean square between two data sets.
  /*!
    Current data is the reference.
    \param data data to be compared with current data.
    \return The root mean square.
  */
  template<class T, int N, class TG>
  template<class T0, class TG0>
  T Data<T, N, TG>::RMS(Data<T0, N, TG0>& data)
  {
    T rms;

    T* data_arr = data_.data();
    T0* data0_arr = data.GetData();
    int NbElements = data_.numElements();

#ifdef SELDONDATA_DEBUG_CHECK_DIMENSIONS

    if (NbElements != data.GetArray().numElements())
      throw WrongDim("Data<T, " + to_str(N) + ">::RMS(Data<T, " + to_str(N) + ">&)",
                     "Data sizes differ.");

#endif

    int nb_elt = 0;
    rms = T(0);
    for (int i = 0; i < NbElements; i++)
      rms += (data_arr[i] - data0_arr[i]) * (data_arr[i] - data0_arr[i]);

    rms = rms / T(NbElements);

    return sqrt(rms);
  }

  //! Computes the relative root mean square between two data sets.
  /*!
    Current data is the reference, and the input data is linearly
    interpolated.
    \param data data to be compared with current data.
    \return The relative root mean square.
  */
  template<class T, int N, class TG>
  template<class T0, class TG0>
  T Data<T, N, TG>::RelativeRMS_interpolation(Data<T0, N, TG0>& data)
  {
    Data<T, N, TG> data0(*this);
    LinearInterpolationGeneral(data, data0);

    return RelativeRMS(data0);
  }

  //! Computes the relative root mean square between two data sets.
  /*!
    Current data is the reference.
    \param data data to be compared with current data.
    \return The relative root mean square.
  */
  template<class T, int N, class TG>
  template<class T0, class TG0>
  T Data<T, N, TG>::RelativeRMS(Data<T0, N, TG0>& data)
  {
    T relative_rms;

    T* data_arr = data_.data();
    T0* data0_arr = data.GetData();
    int NbElements = data_.numElements();

#ifdef SELDONDATA_DEBUG_CHECK_DIMENSIONS

    if (NbElements != data.GetArray().numElements())
      throw WrongDim("Data<T, " + to_str(N) + ">::RelativeRMS(Data<T, " + to_str(N) + ">&)",
                     "Data sizes differ.");

#endif

    int nb_elt = 0;
    relative_rms = T(0);
    for (int i = 0; i < NbElements; i++)
      relative_rms += (data_arr[i] - data0_arr[i]) * (data_arr[i] - data0_arr[i])
        / (data_arr[i] * data_arr[i]);

    relative_rms = relative_rms / T(NbElements);

    return sqrt(relative_rms);
  }

  //! Computes the correlation between two data sets.
  /*!
    Current data is the reference, and the input data is linearly
    interpolated.
    \param data data to be compared with current data.
    \return The correlation.
  */
  template<class T, int N, class TG>
  template<class T0, class TG0>
  T Data<T, N, TG>::Corr_interpolation(Data<T0, N, TG0>& data)
  {
    Data<T, N, TG> data0(*this);
    LinearInterpolationGeneral(data, data0);

    return Corr(data0);
  }

  //! Computes the correlation between two data sets.
  /*!
    Current data is the reference.
    \param data data to be compared with current data.
    \return The correlation.
  */
  template<class T, int N, class TG>
  template<class T0, class TG0>
  T Data<T, N, TG>::Corr(Data<T0, N, TG0>& data)
  {
    T corr;

    T* data_arr = data_.data();
    T0* data0_arr = data.GetData();
    int NbElements = data_.numElements();

#ifdef SELDONDATA_DEBUG_CHECK_DIMENSIONS

    if (NbElements != data.GetArray().numElements())
      throw WrongDim("Data<T, " + to_str(N) + ">::Corr(Data<T, " + to_str(N) + ">&)",
                     "Data sizes differ.");

#endif

    int nb_elt = 0;
    corr = T(0);
    T mean = T(0);
    T0 mean0 = T0(0);
    T var = T(0);
    T0 var0 = T0(0);
    T covar = T(0);
    T temp;
    T0 temp0;

    // Means.
    for (int i = 0; i < NbElements; i++)
      {
        mean += data_arr[i];
        mean0 += data0_arr[i];
      }
    mean = mean / T(NbElements);
    mean0 = mean0 / T0(NbElements);

    // Co-variances.
    for (int i = 0; i < NbElements; i++)
      {
        temp = data_arr[i] - mean;
        temp0 = data0_arr[i] - mean0;
        covar += temp * temp0;
        var += temp * temp;
        var0 += temp0 * temp0;
      }

    corr = covar / sqrt(var * var0);

    return corr;
  }

  //! Returns the percentage of error less than 'threshold'.
  /*!
    Current data is the reference, and the input data is linearly
    interpolated.
    \param data data to be compared with current data.
    \return The percentage of error less than 'threshold'.
  */
  template<class T, int N, class TG>
  template<class T0, class TG0>
  T Data<T, N, TG>::ErrorLessThan_interpolation(Data<T0, N, TG0>& data, T threshold)
  {
    Data<T, N, TG> data0(*this);
    LinearInterpolationGeneral(data, data0);

    return ErrorLessThan(data0, threshold);
  }

  //! Returns the percentage of error less than 'threshold'.
  /*!
    \param data data to be compared with current data.
    \return The percentage of error less than 'threshold'.
  */
  template<class T, int N, class TG>
  template<class T0, class TG0>
  T Data<T, N, TG>::ErrorLessThan(Data<T0, N, TG0>& data, T threshold)
  {
    int nb_err = 0;

    T* data_arr = data_.data();
    T0* data0_arr = data.GetData();
    int NbElements = data_.numElements();

#ifdef SELDONDATA_DEBUG_CHECK_DIMENSIONS

    if (NbElements != data.GetArray().numElements())
      throw WrongDim("Data<T, " + to_str(N) + ">::ErrorLessThan(Data<T, " + to_str(N) + ">&, T)",
                     "Data sizes differ.");

#endif

    for (int i = 0; i < NbElements; i++)
      if (abs(data_arr[i] - data0_arr[i]) <= threshold)
        nb_err++;

    return (T(nb_err) / T(NbElements));
  }

  //! Switches dimensions of data.
  /*! Switches dimensions of data.
    \param NewDim vector of new dimensions indices.
    \param G0 new grid for dimension #0.
  */
  template<class T, int N, class TG>
  void Data<T, N, TG>::SwitchDimensions(TinyVector<int, N> NewDim, Grid<TG>& G0)
  {

#ifdef SELDONDATA_DEBUG_CHECK_DIMENSIONS
    if (N != 1)
      throw WrongDim(string("Data<T, ") + to_str(N) + ">::SwitchDimensions(TinyVector<int, "
                     + to_str(N) + "> NewDim, Grid<TG> G0)",
                     string("Required ") + to_str(N) + " grids, but got 1 grid.");

    // Checks if every dimension is given, and only once.
    TinyVector<bool, N> CheckNewDim;
    CheckNewDim = false;
    for (int i = 0; i < N; i++)
      {
        // Dimension index out of range.
        if (NewDim(i) > N || NewDim(i) < 0)
          throw WrongDim(string("Data<T, ") + to_str(N) + ">::SwitchDimensions(TinyVector<int, "
                         + to_str(N) + "> NewDim, Grid<TG> G0)",
                         string("NewDim(") + to_str(i) + ") is " + to_str(NewDim(i))
                         + " but should be in [0, " + to_str(N - 1) + "].");
        // One dimension is given twice.
        else if (CheckNewDim(NewDim(i)) == true)
          throw WrongDim(string("Data<T, ") + to_str(N) + ">::SwitchDimensions(TinyVector<int, "
                         + to_str(N) + "> NewDim, Grid<TG> G0)",
                         string("Dimension ") + to_str(NewDim(i))
                         + " is given twice or more in NewDim.");
        // First time this dimension is given.
        else
          CheckNewDim(NewDim(i)) = true;
      }

    // Checks if every dimensions is given.
    for (int i = 0; i < N; i++)
      if (CheckNewDim(i) == false)
        throw WrongDim(string("Data<T, ") + to_str(N) + ">::SwitchDimensions(TinyVector<int, "
                       + to_str(N) + "> NewDim, Grid<TG> G0)",
                       string("Dimension ") + to_str(i) + " lacks in NewDim.");
#endif

    int i, j;
    Array<int, 1> Index(10);

    // Gets a copy of the data array.
    Array<T, N> NewData(this->GetArray().data(), this->GetArray().shape(), duplicateData);

    // Re-shapes the data array.
    this->GetArray().resize(this->GetLength(NewDim(0)));
    // Resizes the grid of the data.
    this->ResizeGrid(G0);

    // Puts elements from old array to new array.
    Index = 0;
    for (i = 0; i < NewData.size(); i++)
      {
        this->Value(Index(NewDim(0))) = NewData(Index(0));
        j = N - 1;
        while ((j >= 0) && (Index(j) == NewData.extent(j) - 1))
          {
            Index(j) = 0;
            j--;
          }
        if (j != -1)
          Index(j)++;
      }
  }

  //! Switches dimensions of data.
  /*! Switches dimensions of data.
    \param NewDim vector of new dimensions indices.
    \param G0 new grid for dimension #0.
    \param G1 new grid for dimension #1.
  */
  template<class T, int N, class TG>
  void Data<T, N, TG>::SwitchDimensions(TinyVector<int, N> NewDim, Grid<TG>& G0,
                                        Grid<TG>& G1)
  {

#ifdef SELDONDATA_DEBUG_CHECK_DIMENSIONS
    if (N != 2)
      throw WrongDim(string("Data<T, ") + to_str(N) + ">::SwitchDimensions(TinyVector<int, "
                     + to_str(N) + "> NewDim, Grid<TG> G0, ...)",
                     string("Required ") + to_str(N) + " grids, but got 2 grids.");

    // Checks if every dimension is given, and only once.
    TinyVector<bool, N> CheckNewDim;
    CheckNewDim = false;
    for (int i = 0; i < N; i++)
      {
        // Dimension index out of range.
        if (NewDim(i) > N || NewDim(i) < 0)
          throw WrongDim(string("Data<T, ") + to_str(N) + ">::SwitchDimensions(TinyVector<int, "
                         + to_str(N) + "> NewDim, Grid<TG> G0, ...)",
                         string("NewDim(") + to_str(i) + ") is " + to_str(NewDim(i))
                         + " but should be in [0, " + to_str(N - 1) + "].");
        // One dimension is given twice.
        else if (CheckNewDim(NewDim(i)) == true)
          throw WrongDim(string("Data<T, ") + to_str(N) + ">::SwitchDimensions(TinyVector<int, "
                         + to_str(N) + "> NewDim, Grid<TG> G0, ...)",
                         string("Dimension ") + to_str(NewDim(i))
                         + " is given twice or more in NewDim.");
        // First time this dimension is given.
        else
          CheckNewDim(NewDim(i)) = true;
      }

    // Checks if every dimensions is given.
    for (int i = 0; i < N; i++)
      if (CheckNewDim(i) == false)
        throw WrongDim(string("Data<T, ") + to_str(N) + ">::SwitchDimensions(TinyVector<int, "
                       + to_str(N) + "> NewDim, Grid<TG> G0, ...)",
                       string("Dimension ") + to_str(i) + " lacks in NewDim.");
#endif

    int i, j;
    Array<int, 1> Index(10);

    // Gets a copy of the data array.
    Array<T, N> NewData(this->GetArray().data(), this->GetArray().shape(), duplicateData);

    // Re-shapes the data array.
    this->GetArray().resize(this->GetLength(NewDim(0)), this->GetLength(NewDim(1)));
    // Resizes the grid of the data.
    this->ResizeGrid(G0, G1);

    // Puts elements from old array to new array.
    Index = 0;
    for (i = 0; i < NewData.size(); i++)
      {
        this->Value(Index(NewDim(0)), Index(NewDim(1)))
          = NewData(Index(0), Index(1));
        j = N - 1;
        while ((j >= 0) && (Index(j) == NewData.extent(j) - 1))
          {
            Index(j) = 0;
            j--;
          }
        if (j != -1)
          Index(j)++;
      }
  }

  //! Switches dimensions of data.
  /*! Switches dimensions of data.
    \param NewDim vector of new dimensions indices.
    \param G0 new grid for dimension #0.
    \param G1 new grid for dimension #1.
    \param G2 new grid for dimension #2.
  */
  template<class T, int N, class TG>
  void Data<T, N, TG>::SwitchDimensions(TinyVector<int, N> NewDim, Grid<TG>& G0,
                                        Grid<TG>& G1, Grid<TG>& G2)
  {

#ifdef SELDONDATA_DEBUG_CHECK_DIMENSIONS
    if (N != 3)
      throw WrongDim(string("Data<T, ") + to_str(N) + ">::SwitchDimensions(TinyVector<int, "
                     + to_str(N) + "> NewDim, Grid<TG> G0, ...)",
                     string("Required ") + to_str(N) + " grids, but got 3 grids.");

    // Checks if every dimension is given, and only once.
    TinyVector<bool, N> CheckNewDim;
    CheckNewDim = false;
    for (int i = 0; i < N; i++)
      {
        // Dimension index out of range.
        if (NewDim(i) > N || NewDim(i) < 0)
          throw WrongDim(string("Data<T, ") + to_str(N) + ">::SwitchDimensions(TinyVector<int, "
                         + to_str(N) + "> NewDim, Grid<TG> G0, ...)",
                         string("NewDim(") + to_str(i) + ") is " + to_str(NewDim(i))
                         + " but should be in [0, " + to_str(N - 1) + "].");
        // One dimension is given twice.
        else if (CheckNewDim(NewDim(i)) == true)
          throw WrongDim(string("Data<T, ") + to_str(N) + ">::SwitchDimensions(TinyVector<int, "
                         + to_str(N) + "> NewDim, Grid<TG> G0, ...)",
                         string("Dimension ") + to_str(NewDim(i))
                         + " is given twice or more in NewDim.");
        // First time this dimension is given.
        else
          CheckNewDim(NewDim(i)) = true;
      }

    // Checks if every dimensions is given.
    for (int i = 0; i < N; i++)
      if (CheckNewDim(i) == false)
        throw WrongDim(string("Data<T, ") + to_str(N) + ">::SwitchDimensions(TinyVector<int, "
                       + to_str(N) + "> NewDim, Grid<TG> G0, ...)",
                       string("Dimension ") + to_str(i) + " lacks in NewDim.");
#endif

    int i, j;
    Array<int, 1> Index(10);

    // Gets a copy of the data array.
    Array<T, N> NewData(this->GetArray().data(), this->GetArray().shape(), duplicateData);

    // Re-shapes the data array.
    this->GetArray().resize(this->GetLength(NewDim(0)), this->GetLength(NewDim(1)),
                            this->GetLength(NewDim(2)));
    // Resizes the grid of the data.
    this->ResizeGrid(G0, G1, G2);

    // Puts elements from old array to new array.
    Index = 0;
    for (i = 0; i < (int) NewData.size(); i++)
      {
        this->Value(Index(NewDim(0)), Index(NewDim(1)), Index(NewDim(2)))
          = NewData(Index(0), Index(1), Index(2));
        j = N - 1;
        while ((j >= 0) && (Index(j) == NewData.extent(j) - 1))
          {
            Index(j) = 0;
            j--;
          }
        if (j != -1)
          Index(j)++;
      }
  }

  //! Switches dimensions of data.
  /*! Switches dimensions of data.
    \param NewDim vector of new dimensions indices.
    \param G0 new grid for dimension #0.
    \param G1 new grid for dimension #1.
    \param G2 new grid for dimension #2.
    \param G3 new grid for dimension #3.
  */
  template<class T, int N, class TG>
  void Data<T, N, TG>::SwitchDimensions(TinyVector<int, N> NewDim, Grid<TG>& G0,
                                        Grid<TG>& G1, Grid<TG>& G2, Grid<TG>& G3)
  {

#ifdef SELDONDATA_DEBUG_CHECK_DIMENSIONS
    if (N != 4)
      throw WrongDim(string("Data<T, ") + to_str(N) + ">::SwitchDimensions(TinyVector<int, "
                     + to_str(N) + "> NewDim, Grid<TG> G0, ...)",
                     string("Required ") + to_str(N) + " grids, but got 4 grids.");

    // Checks if every dimension is given, and only once.
    TinyVector<bool, N> CheckNewDim;
    CheckNewDim = false;
    for (int i = 0; i < N; i++)
      {
        // Dimension index out of range.
        if (NewDim(i) > N || NewDim(i) < 0)
          throw WrongDim(string("Data<T, ") + to_str(N) + ">::SwitchDimensions(TinyVector<int, "
                         + to_str(N) + "> NewDim, Grid<TG> G0, ...)",
                         string("NewDim(") + to_str(i) + ") is " + to_str(NewDim(i))
                         + " but should be in [0, " + to_str(N - 1) + "].");
        // One dimension is given twice.
        else if (CheckNewDim(NewDim(i)) == true)
          throw WrongDim(string("Data<T, ") + to_str(N) + ">::SwitchDimensions(TinyVector<int, "
                         + to_str(N) + "> NewDim, Grid<TG> G0, ...)",
                         string("Dimension ") + to_str(NewDim(i))
                         + " is given twice or more in NewDim.");
        // First time this dimension is given.
        else
          CheckNewDim(NewDim(i)) = true;
      }

    // Checks if every dimensions is given.
    for (int i = 0; i < N; i++)
      if (CheckNewDim(i) == false)
        throw WrongDim(string("Data<T, ") + to_str(N) + ">::SwitchDimensions(TinyVector<int, "
                       + to_str(N) + "> NewDim, Grid<TG> G0, ...)",
                       string("Dimension ") + to_str(i) + " lacks in NewDim.");
#endif

    int i, j;
    Array<int, 1> Index(10);

    // Gets a copy of the data array.
    Array<T, N> NewData(this->GetArray().data(), this->GetArray().shape(), duplicateData);

    // Re-shapes the data array.
    this->GetArray().resize(this->GetLength(NewDim(0)), this->GetLength(NewDim(1)),
                            this->GetLength(NewDim(2)), this->GetLength(NewDim(3)));
    // Resizes the grid of the data.
    this->ResizeGrid(G0, G1, G2, G3);

    // Puts elements from old array to new array.
    Index = 0;
    for (i = 0; i < (int) NewData.size(); i++)
      {
        this->Value(Index(NewDim(0)), Index(NewDim(1)),
                    Index(NewDim(2)), Index(NewDim(3)))
          = NewData(Index(0), Index(1), Index(2), Index(3));
        j = N - 1;
        while ((j >= 0) && (Index(j) == NewData.extent(j) - 1))
          {
            Index(j) = 0;
            j--;
          }
        if (j != -1)
          Index(j)++;
      }
  }

  //! Switches dimensions of data.
  /*! Switches dimensions of data.
    \param NewDim vector of new dimensions indices.
    \param G0 new grid for dimension #0.
    \param G1 new grid for dimension #1.
    \param G2 new grid for dimension #2.
    \param G3 new grid for dimension #3.
    \param G4 new grid for dimension #4.
  */
  template<class T, int N, class TG>
  void Data<T, N, TG>::SwitchDimensions(TinyVector<int, N> NewDim, Grid<TG>& G0,
                                        Grid<TG>& G1, Grid<TG>& G2, Grid<TG>& G3,
                                        Grid<TG>& G4)
  {

#ifdef SELDONDATA_DEBUG_CHECK_DIMENSIONS
    if (N != 5)
      throw WrongDim(string("Data<T, ") + to_str(N) + ">::SwitchDimensions(TinyVector<int, "
                     + to_str(N) + "> NewDim, Grid<TG> G0, ...)",
                     string("Required ") + to_str(N) + " grids, but got 5 grids.");

    // Checks if every dimension is given, and only once.
    TinyVector<bool, N> CheckNewDim;
    CheckNewDim = false;
    for (int i = 0; i < N; i++)
      {
        // Dimension index out of range.
        if (NewDim(i) > N || NewDim(i) < 0)
          throw WrongDim(string("Data<T, ") + to_str(N) + ">::SwitchDimensions(TinyVector<int, "
                         + to_str(N) + "> NewDim, Grid<TG> G0, ...)",
                         string("NewDim(") + to_str(i) + ") is " + to_str(NewDim(i))
                         + " but should be in [0, " + to_str(N - 1) + "].");
        // One dimension is given twice.
        else if (CheckNewDim(NewDim(i)) == true)
          throw WrongDim(string("Data<T, ") + to_str(N) + ">::SwitchDimensions(TinyVector<int, "
                         + to_str(N) + "> NewDim, Grid<TG> G0, ...)",
                         string("Dimension ") + to_str(NewDim(i))
                         + " is given twice or more in NewDim.");
        // First time this dimension is given.
        else
          CheckNewDim(NewDim(i)) = true;
      }

    // Checks if every dimensions is given.
    for (int i = 0; i < N; i++)
      if (CheckNewDim(i) == false)
        throw WrongDim(string("Data<T, ") + to_str(N) + ">::SwitchDimensions(TinyVector<int, "
                       + to_str(N) + "> NewDim, Grid<TG> G0, ...)",
                       string("Dimension ") + to_str(i) + " lacks in NewDim.");
#endif

    int i, j;
    Array<int, 1> Index(10);

    // Gets a copy of the data array.
    Array<T, N> NewData(this->GetArray().data(), this->GetArray().shape(), duplicateData);

    // Re-shapes the data array.
    this->GetArray().resize(this->GetLength(NewDim(0)), this->GetLength(NewDim(1)),
                            this->GetLength(NewDim(2)), this->GetLength(NewDim(3)),
                            this->GetLength(NewDim(4)));
    // Resizes the grid of the data.
    this->ResizeGrid(G0, G1, G2, G3, G4);

    // Puts elements from old array to new array.
    Index = 0;
    for (i = 0; i < NewData.size(); i++)
      {
        this->Value(Index(NewDim(0)), Index(NewDim(1)),
                    Index(NewDim(2)), Index(NewDim(3)),
                    Index(NewDim(4)))
          = NewData(Index(0), Index(1), Index(2), Index(3), Index(4));
        j = N - 1;
        while ((j >= 0) && (Index(j) == NewData.extent(j) - 1))
          {
            Index(j) = 0;
            j--;
          }
        if (j != -1)
          Index(j)++;
      }
  }

  //! Switches dimensions of data.
  /*! Switches dimensions of data.
    \param NewDim vector of new dimensions indices.
    \param G0 new grid for dimension #0.
    \param G1 new grid for dimension #1.
    \param G2 new grid for dimension #2.
    \param G3 new grid for dimension #3.
    \param G4 new grid for dimension #4.
    \param G5 new grid for dimension #5.
  */
  template<class T, int N, class TG>
  void Data<T, N, TG>::SwitchDimensions(TinyVector<int, N> NewDim, Grid<TG>& G0,
                                        Grid<TG>& G1, Grid<TG>& G2, Grid<TG>& G3,
                                        Grid<TG>& G4, Grid<TG>& G5)
  {

#ifdef SELDONDATA_DEBUG_CHECK_DIMENSIONS
    if (N != 6)
      throw WrongDim(string("Data<T, ") + to_str(N) + ">::SwitchDimensions(TinyVector<int, "
                     + to_str(N) + "> NewDim, Grid<TG> G0, ...)",
                     string("Required ") + to_str(N) + " grids, but got 6 grids.");

    // Checks if every dimension is given, and only once.
    TinyVector<bool, N> CheckNewDim;
    CheckNewDim = false;
    for (int i = 0; i < N; i++)
      {
        // Dimension index out of range.
        if (NewDim(i) > N || NewDim(i) < 0)
          throw WrongDim(string("Data<T, ") + to_str(N) + ">::SwitchDimensions(TinyVector<int, "
                         + to_str(N) + "> NewDim, Grid<TG> G0, ...)",
                         string("NewDim(") + to_str(i) + ") is " + to_str(NewDim(i))
                         + " but should be in [0, " + to_str(N - 1) + "].");
        // One dimension is given twice.
        else if (CheckNewDim(NewDim(i)) == true)
          throw WrongDim(string("Data<T, ") + to_str(N) + ">::SwitchDimensions(TinyVector<int, "
                         + to_str(N) + "> NewDim, Grid<TG> G0, ...)",
                         string("Dimension ") + to_str(NewDim(i))
                         + " is given twice or more in NewDim.");
        // First time this dimension is given.
        else
          CheckNewDim(NewDim(i)) = true;
      }

    // Checks if every dimensions is given.
    for (int i = 0; i < N; i++)
      if (CheckNewDim(i) == false)
        throw WrongDim(string("Data<T, ") + to_str(N) + ">::SwitchDimensions(TinyVector<int, "
                       + to_str(N) + "> NewDim, Grid<TG> G0, ...)",
                       string("Dimension ") + to_str(i) + " lacks in NewDim.");
#endif

    int i, j;
    Array<int, 1> Index(10);

    // Gets a copy of the data array.
    Array<T, N> NewData(this->GetArray().data(), this->GetArray().shape(), duplicateData);

    // Re-shapes the data array.
    this->GetArray().resize(this->GetLength(NewDim(0)), this->GetLength(NewDim(1)),
                            this->GetLength(NewDim(2)), this->GetLength(NewDim(3)),
                            this->GetLength(NewDim(4)), this->GetLength(NewDim(5)));
    // Resizes the grid of the data.
    this->ResizeGrid(G0, G1, G2, G3, G4, G5);

    // Puts elements from old array to new array.
    Index = 0;
    for (i = 0; i < NewData.size(); i++)
      {
        this->Value(Index(NewDim(0)), Index(NewDim(1)),
                    Index(NewDim(2)), Index(NewDim(3)),
                    Index(NewDim(4)), Index(NewDim(5)))
          = NewData(Index(0), Index(1), Index(2), Index(3), Index(4), Index(5));
        j = N - 1;
        while ((j >= 0) && (Index(j) == NewData.extent(j) - 1))
          {
            Index(j) = 0;
            j--;
          }
        if (j != -1)
          Index(j)++;
      }
  }

  //! Switches dimensions of data.
  /*! Switches dimensions of data.
    \param NewDim vector of new dimensions indices.
    \param G0 new grid for dimension #0.
    \param G1 new grid for dimension #1.
    \param G2 new grid for dimension #2.
    \param G3 new grid for dimension #3.
    \param G4 new grid for dimension #4.
    \param G5 new grid for dimension #5.
    \param G6 new grid for dimension #6.
  */
  template<class T, int N, class TG>
  void Data<T, N, TG>::SwitchDimensions(TinyVector<int, N> NewDim, Grid<TG>& G0,
                                        Grid<TG>& G1, Grid<TG>& G2, Grid<TG>& G3,
                                        Grid<TG>& G4, Grid<TG>& G5, Grid<TG>& G6)
  {

#ifdef SELDONDATA_DEBUG_CHECK_DIMENSIONS
    if (N != 7)
      throw WrongDim(string("Data<T, ") + to_str(N) + ">::SwitchDimensions(TinyVector<int, "
                     + to_str(N) + "> NewDim, Grid<TG> G0, ...)",
                     string("Required ") + to_str(N) + " grids, but got 7 grids.");

    // Checks if every dimension is given, and only once.
    TinyVector<bool, N> CheckNewDim;
    CheckNewDim = false;
    for (int i = 0; i < N; i++)
      {
        // Dimension index out of range.
        if (NewDim(i) > N || NewDim(i) < 0)
          throw WrongDim(string("Data<T, ") + to_str(N) + ">::SwitchDimensions(TinyVector<int, "
                         + to_str(N) + "> NewDim, Grid<TG> G0, ...)",
                         string("NewDim(") + to_str(i) + ") is " + to_str(NewDim(i))
                         + " but should be in [0, " + to_str(N - 1) + "].");
        // One dimension is given twice.
        else if (CheckNewDim(NewDim(i)) == true)
          throw WrongDim(string("Data<T, ") + to_str(N) + ">::SwitchDimensions(TinyVector<int, "
                         + to_str(N) + "> NewDim, Grid<TG> G0, ...)",
                         string("Dimension ") + to_str(NewDim(i))
                         + " is given twice or more in NewDim.");
        // First time this dimension is given.
        else
          CheckNewDim(NewDim(i)) = true;
      }

    // Checks if every dimensions is given.
    for (int i = 0; i < N; i++)
      if (CheckNewDim(i) == false)
        throw WrongDim(string("Data<T, ") + to_str(N) + ">::SwitchDimensions(TinyVector<int, "
                       + to_str(N) + "> NewDim, Grid<TG> G0, ...)",
                       string("Dimension ") + to_str(i) + " lacks in NewDim.");
#endif

    int i, j;
    Array<int, 1> Index(10);

    // Gets a copy of the data array.
    Array<T, N> NewData(this->GetArray().data(), this->GetArray().shape(), duplicateData);

    // Re-shapes the data array.
    this->GetArray().resize(this->GetLength(NewDim(0)), this->GetLength(NewDim(1)),
                            this->GetLength(NewDim(2)), this->GetLength(NewDim(3)),
                            this->GetLength(NewDim(4)), this->GetLength(NewDim(5)),
                            this->GetLength(NewDim(6)));
    // Resizes the grid of the data.
    this->ResizeGrid(G0, G1, G2, G3, G4, G5, G6);

    // Puts elements from old array to new array.
    Index = 0;
    for (i = 0; i < NewData.size(); i++)
      {
        this->Value(Index(NewDim(0)), Index(NewDim(1)),
                    Index(NewDim(2)), Index(NewDim(3)),
                    Index(NewDim(4)), Index(NewDim(5)),
                    Index(NewDim(6)))
          = NewData(Index(0), Index(1), Index(2), Index(3), Index(4), Index(5), Index(6));
        j = N - 1;
        while ((j >= 0) && (Index(j) == NewData.extent(j) - 1))
          {
            Index(j) = 0;
            j--;
          }
        if (j != -1)
          Index(j)++;
      }
  }

  //! Switches dimensions of data.
  /*! Switches dimensions of data.
    \param NewDim vector of new dimensions indices.
    \param G0 new grid for dimension #0.
    \param G1 new grid for dimension #1.
    \param G2 new grid for dimension #2.
    \param G3 new grid for dimension #3.
    \param G4 new grid for dimension #4.
    \param G5 new grid for dimension #5.
    \param G6 new grid for dimension #6.
    \param G7 new grid for dimension #7.
  */
  template<class T, int N, class TG>
  void Data<T, N, TG>::SwitchDimensions(TinyVector<int, N> NewDim, Grid<TG>& G0,
                                        Grid<TG>& G1, Grid<TG>& G2, Grid<TG>& G3,
                                        Grid<TG>& G4, Grid<TG>& G5, Grid<TG>& G6,
                                        Grid<TG>& G7)
  {

#ifdef SELDONDATA_DEBUG_CHECK_DIMENSIONS
    if (N != 8)
      throw WrongDim(string("Data<T, ") + to_str(N) + ">::SwitchDimensions(TinyVector<int, "
                     + to_str(N) + "> NewDim, Grid<TG> G0, ...)",
                     string("Required ") + to_str(N) + " grids, but got 8 grids.");

    // Checks if every dimension is given, and only once.
    TinyVector<bool, N> CheckNewDim;
    CheckNewDim = false;
    for (int i = 0; i < N; i++)
      {
        // Dimension index out of range.
        if (NewDim(i) > N || NewDim(i) < 0)
          throw WrongDim(string("Data<T, ") + to_str(N) + ">::SwitchDimensions(TinyVector<int, "
                         + to_str(N) + "> NewDim, Grid<TG> G0, ...)",
                         string("NewDim(") + to_str(i) + ") is " + to_str(NewDim(i))
                         + " but should be in [0, " + to_str(N - 1) + "].");
        // One dimension is given twice.
        else if (CheckNewDim(NewDim(i)) == true)
          throw WrongDim(string("Data<T, ") + to_str(N) + ">::SwitchDimensions(TinyVector<int, "
                         + to_str(N) + "> NewDim, Grid<TG> G0, ...)",
                         string("Dimension ") + to_str(NewDim(i))
                         + " is given twice or more in NewDim.");
        // First time this dimension is given.
        else
          CheckNewDim(NewDim(i)) = true;
      }

    // Checks if every dimensions is given.
    for (int i = 0; i < N; i++)
      if (CheckNewDim(i) == false)
        throw WrongDim(string("Data<T, ") + to_str(N) + ">::SwitchDimensions(TinyVector<int, "
                       + to_str(N) + "> NewDim, Grid<TG> G0, ...)",
                       string("Dimension ") + to_str(i) + " lacks in NewDim.");
#endif

    int i, j;
    Array<int, 1> Index(10);

    // Gets a copy of the data array.
    Array<T, N> NewData(this->GetArray().data(), this->GetArray().shape(), duplicateData);

    // Re-shapes the data array.
    this->GetArray().resize(this->GetLength(NewDim(0)), this->GetLength(NewDim(1)),
                            this->GetLength(NewDim(2)), this->GetLength(NewDim(3)),
                            this->GetLength(NewDim(4)), this->GetLength(NewDim(5)),
                            this->GetLength(NewDim(6)), this->GetLength(NewDim(7)));
    // Resizes the grid of the data.
    this->ResizeGrid(G0, G1, G2, G3, G4, G5, G6, G7);

    // Puts elements from old array to new array.
    Index = 0;
    for (i = 0; i < NewData.size(); i++)
      {
        this->Value(Index(NewDim(0)), Index(NewDim(1)),
                    Index(NewDim(2)), Index(NewDim(3)),
                    Index(NewDim(4)), Index(NewDim(5)),
                    Index(NewDim(6)), Index(NewDim(7)))
          = NewData(Index(0), Index(1), Index(2), Index(3),
                    Index(4), Index(5), Index(6), Index(7));
        j = N - 1;
        while ((j >= 0) && (Index(j) == NewData.extent(j) - 1))
          {
            Index(j) = 0;
            j--;
          }
        if (j != -1)
          Index(j)++;
      }
  }

  //! Switches dimensions of data.
  /*! Switches dimensions of data.
    \param NewDim vector of new dimensions indices.
    \param G0 new grid for dimension #0.
    \param G1 new grid for dimension #1.
    \param G2 new grid for dimension #2.
    \param G3 new grid for dimension #3.
    \param G4 new grid for dimension #4.
    \param G5 new grid for dimension #5.
    \param G6 new grid for dimension #6.
    \param G7 new grid for dimension #7.
    \param G8 new grid for dimension #8.
  */
  template<class T, int N, class TG>
  void Data<T, N, TG>::SwitchDimensions(TinyVector<int, N> NewDim, Grid<TG>& G0,
                                        Grid<TG>& G1, Grid<TG>& G2, Grid<TG>& G3,
                                        Grid<TG>& G4, Grid<TG>& G5, Grid<TG>& G6,
                                        Grid<TG>& G7, Grid<TG>& G8)
  {

#ifdef SELDONDATA_DEBUG_CHECK_DIMENSIONS
    if (N != 9)
      throw WrongDim(string("Data<T, ") + to_str(N) + ">::SwitchDimensions(TinyVector<int, "
                     + to_str(N) + "> NewDim, Grid<TG> G0, ...)",
                     string("Required ") + to_str(N) + " grids, but got 9 grids.");

    // Checks if every dimension is given, and only once.
    TinyVector<bool, N> CheckNewDim;
    CheckNewDim = false;
    for (int i = 0; i < N; i++)
      {
        // Dimension index out of range.
        if (NewDim(i) > N || NewDim(i) < 0)
          throw WrongDim(string("Data<T, ") + to_str(N) + ">::SwitchDimensions(TinyVector<int, "
                         + to_str(N) + "> NewDim, Grid<TG> G0, ...)",
                         string("NewDim(") + to_str(i) + ") is " + to_str(NewDim(i))
                         + " but should be in [0, " + to_str(N - 1) + "].");
        // One dimension is given twice.
        else if (CheckNewDim(NewDim(i)) == true)
          throw WrongDim(string("Data<T, ") + to_str(N) + ">::SwitchDimensions(TinyVector<int, "
                         + to_str(N) + "> NewDim, Grid<TG> G0, ...)",
                         string("Dimension ") + to_str(NewDim(i))
                         + " is given twice or more in NewDim.");
        // First time this dimension is given.
        else
          CheckNewDim(NewDim(i)) = true;
      }

    // Checks if every dimensions is given.
    for (int i = 0; i < N; i++)
      if (CheckNewDim(i) == false)
        throw WrongDim(string("Data<T, ") + to_str(N) + ">::SwitchDimensions(TinyVector<int, "
                       + to_str(N) + "> NewDim, Grid<TG> G0, ...)",
                       string("Dimension ") + to_str(i) + " lacks in NewDim.");
#endif

    int i, j;
    Array<int, 1> Index(10);

    // Gets a copy of the data array.
    Array<T, N> NewData(this->GetArray().data(), this->GetArray().shape(), duplicateData);

    // Re-shapes the data array.
    this->GetArray().resize(this->GetLength(NewDim(0)), this->GetLength(NewDim(1)),
                            this->GetLength(NewDim(2)), this->GetLength(NewDim(3)),
                            this->GetLength(NewDim(4)), this->GetLength(NewDim(5)),
                            this->GetLength(NewDim(6)), this->GetLength(NewDim(7)),
                            this->GetLength(NewDim(8)));
    // Resizes the grid of the data.
    this->ResizeGrid(G0, G1, G2, G3, G4, G5, G6, G7, G8);

    // Puts elements from old array to new array.
    Index = 0;
    for (i = 0; i < NewData.size(); i++)
      {
        this->Value(Index(NewDim(0)), Index(NewDim(1)),
                    Index(NewDim(2)), Index(NewDim(3)),
                    Index(NewDim(4)), Index(NewDim(5)),
                    Index(NewDim(6)), Index(NewDim(7)),
                    Index(NewDim(8)))
          = NewData(Index(0), Index(1), Index(2), Index(3),
                    Index(4), Index(5), Index(6), Index(7),
                    Index(8));
        j = N - 1;
        while ((j >= 0) && (Index(j) == NewData.extent(j) - 1))
          {
            Index(j) = 0;
            j--;
          }
        if (j != -1)
          Index(j)++;
      }
  }

  //! Switches dimensions of data.
  /*! Switches dimensions of data.
    \param NewDim vector of new dimensions indices.
    \param G0 new grid for dimension #0.
    \param G1 new grid for dimension #1.
    \param G2 new grid for dimension #2.
    \param G3 new grid for dimension #3.
    \param G4 new grid for dimension #4.
    \param G5 new grid for dimension #5.
    \param G6 new grid for dimension #6.
    \param G7 new grid for dimension #7.
    \param G8 new grid for dimension #8.
    \param G9 new grid for dimension #9.
  */
  template<class T, int N, class TG>
  void Data<T, N, TG>::SwitchDimensions(TinyVector<int, N> NewDim, Grid<TG>& G0,
                                        Grid<TG>& G1, Grid<TG>& G2, Grid<TG>& G3,
                                        Grid<TG>& G4, Grid<TG>& G5, Grid<TG>& G6,
                                        Grid<TG>& G7, Grid<TG>& G8, Grid<TG>& G9)
  {

#ifdef SELDONDATA_DEBUG_CHECK_DIMENSIONS
    if (N != 10)
      throw WrongDim(string("Data<T, ") + to_str(N) + ">::SwitchDimensions(TinyVector<int, "
                     + to_str(N) + "> NewDim, Grid<TG> G0, ...)",
                     string("Required ") + to_str(N) + " grids, but got 10 grids.");

    // Checks if every dimension is given, and only once.
    TinyVector<bool, N> CheckNewDim;
    CheckNewDim = false;
    for (int i = 0; i < N; i++)
      {
        // Dimension index out of range.
        if (NewDim(i) > N || NewDim(i) < 0)
          throw WrongDim(string("Data<T, ") + to_str(N) + ">::SwitchDimensions(TinyVector<int, "
                         + to_str(N) + "> NewDim, Grid<TG> G0, ...)",
                         string("NewDim(") + to_str(i) + ") is " + to_str(NewDim(i))
                         + " but should be in [0, " + to_str(N - 1) + "].");
        // One dimension is given twice.
        else if (CheckNewDim(NewDim(i)) == true)
          throw WrongDim(string("Data<T, ") + to_str(N) + ">::SwitchDimensions(TinyVector<int, "
                         + to_str(N) + "> NewDim, Grid<TG> G0, ...)",
                         string("Dimension ") + to_str(NewDim(i))
                         + " is given twice or more in NewDim.");
        // First time this dimension is given.
        else
          CheckNewDim(NewDim(i)) = true;
      }

    // Checks if every dimensions is given.
    for (int i = 0; i < N; i++)
      if (CheckNewDim(i) == false)
        throw WrongDim(string("Data<T, ") + to_str(N) + ">::SwitchDimensions(TinyVector<int, "
                       + to_str(N) + "> NewDim, Grid<TG> G0, ...)",
                       string("Dimension ") + to_str(i) + " lacks in NewDim.");
#endif

    int i, j;
    Array<int, 1> Index(10);

    // Gets a copy of the data array.
    Array<T, N> NewData(this->GetArray().data(), this->GetArray().shape(), duplicateData);

    // Re-shapes the data array.
    this->GetArray().resize(this->GetLength(NewDim(0)), this->GetLength(NewDim(1)),
                            this->GetLength(NewDim(2)), this->GetLength(NewDim(3)),
                            this->GetLength(NewDim(4)), this->GetLength(NewDim(5)),
                            this->GetLength(NewDim(6)), this->GetLength(NewDim(7)),
                            this->GetLength(NewDim(8)), this->GetLength(NewDim(9)));
    // Resizes the grid of the data.
    this->ResizeGrid(G0, G1, G2, G3, G4, G5, G6, G7, G8, G9);

    // Puts elements from old array to new array.
    Index = 0;
    for (i = 0; i < NewData.size(); i++)
      {
        this->Value(Index(NewDim(0)), Index(NewDim(1)),
                    Index(NewDim(2)), Index(NewDim(3)),
                    Index(NewDim(4)), Index(NewDim(5)),
                    Index(NewDim(6)), Index(NewDim(7)),
                    Index(NewDim(8)), Index(NewDim(9)))
          = NewData(Index(0), Index(1), Index(2), Index(3),
                    Index(4), Index(5), Index(6), Index(7),
                    Index(8), Index(9));
        j = N - 1;
        while ((j >= 0) && (Index(j) == NewData.extent(j) - 1))
          {
            Index(j) = 0;
            j--;
          }
        if (j != -1)
          Index(j)++;
      }
  }

  //! Reverses data along a given dimension.
  /*! Reverses data along dimension 'dim'. For example, if a three dimensional
    array A is reversed along dimension #1, on exit, A(i, j, k) = A(i, Ny-1-j, k),
    where Ny is the length along dimension #1.
    \param j dimension
  */
  template<class T, int N, class TG>
  void Data<T, N, TG>::ReverseData(int dim)
  {
    int i, j, k;

    int Ndim = this->GetLength(dim);

    Array<T, 1> tmp(Ndim);
    Array<int, 1> Index(10);

    for (int i = 0; i < 10; i++)
      Index(i) = 0;

    j = N - 1;
    for (i = 0; i < this->GetNbElements() / Ndim; i++)
      {

        for (k = 0; k < Ndim; k++)
          {
            Index(dim) = k;
            tmp(Ndim - 1 - k) = this->Value(Index(0), Index(1), Index(2),
                                            Index(3), Index(4), Index(5),
                                            Index(6), Index(7), Index(8),
                                            Index(9));
          }

        for (k = 0; k < Ndim; k++)
          {
            Index(dim) = k;
            this->Value(Index(0), Index(1), Index(2),
                        Index(3), Index(4), Index(5),
                        Index(6), Index(7), Index(8),
                        Index(9)) = tmp(k);
          }

        j = N - 1;
        while ((j == dim) || ((j >= 0) && (Index(j) == this->GetLength(j) - 1)))
          {
            Index(j) = 0;
            j--;
          }
        if (j != -1)
          Index(j)++;
      }
  }

  //! Change coordinates.
  /*!
    The coordinates transformation is provided by function
    f. f takes as input old coordinates and new coordinates (in this
    order). Old coordinates are assumed to be unchanged by f.
    \param f coordinates transformation. It must be an instance of
    FuncCoords_Base<TG> or of a class derived from FuncCoords_Base<TG>.
  */
  template<class T, int N, class TG>
  void Data<T, N, TG>::ChangeCoords(FuncCoords_Base<TG>& f)
  {

    int i;
    int i0, i1, i2, i3, i4, i5, i6, i7, i8, i9;

    Array<Grid<TG>*, 1> grids(N);
    for (i = 0; i < N; i++)
      grids(i) = grids_(i)->Duplicate();

    if (N == 1)
      for (i0 = 0; i0 < grids(0)->GetLength(); i0++)
        f(grids(0)->Value(i0),
          grids_(0)->Value(i0));
    else if (N == 2)
      for (i0 = 0; i0 < grids(0)->GetLength(); i0++)
        for (i1 = 0; i1 < grids(1)->GetLength(); i1++)
          f(grids(0)->Value(i0, i1),
            grids(1)->Value(i0, i1),
            grids_(0)->Value(i0, i1),
            grids_(1)->Value(i0, i1));
    else if (N == 3)
      for (i0 = 0; i0 < grids(0)->GetLength(); i0++)
        for (i1 = 0; i1 < grids(1)->GetLength(); i1++)
          for (i2 = 0; i2 < grids(2)->GetLength(); i2++)
            f(grids(0)->Value(i0, i1, i2),
              grids(1)->Value(i0, i1, i2),
              grids(2)->Value(i0, i1, i2),
              grids_(0)->Value(i0, i1, i2),
              grids_(1)->Value(i0, i1, i2),
              grids_(2)->Value(i0, i1, i2));
    else if (N == 4)
      for (i0 = 0; i0 < grids(0)->GetLength(); i0++)
        for (i1 = 0; i1 < grids(1)->GetLength(); i1++)
          for (i2 = 0; i2 < grids(2)->GetLength(); i2++)
            for (i3 = 0; i3 < grids(3)->GetLength(); i3++)
              f(grids(0)->Value(i0, i1, i2, i3),
                grids(1)->Value(i0, i1, i2, i3),
                grids(2)->Value(i0, i1, i2, i3),
                grids(3)->Value(i0, i1, i2, i3),
                grids_(0)->Value(i0, i1, i2, i3),
                grids_(1)->Value(i0, i1, i2, i3),
                grids_(2)->Value(i0, i1, i2, i3),
                grids_(3)->Value(i0, i1, i2, i3));
    else if (N == 5)
      for (i0 = 0; i0 < grids(0)->GetLength(); i0++)
        for (i1 = 0; i1 < grids(1)->GetLength(); i1++)
          for (i2 = 0; i2 < grids(2)->GetLength(); i2++)
            for (i3 = 0; i3 < grids(3)->GetLength(); i3++)
              for (i4 = 0; i4 < grids(4)->GetLength(); i4++)
                f(grids(0)->Value(i0, i1, i2, i3, i4),
                  grids(1)->Value(i0, i1, i2, i3, i4),
                  grids(2)->Value(i0, i1, i2, i3, i4),
                  grids(3)->Value(i0, i1, i2, i3, i4),
                  grids(4)->Value(i0, i1, i2, i3, i4),
                  grids_(0)->Value(i0, i1, i2, i3, i4),
                  grids_(1)->Value(i0, i1, i2, i3, i4),
                  grids_(2)->Value(i0, i1, i2, i3, i4),
                  grids_(3)->Value(i0, i1, i2, i3, i4),
                  grids_(4)->Value(i0, i1, i2, i3, i4));
    else if (N == 6)
      for (i0 = 0; i0 < grids(0)->GetLength(); i0++)
        for (i1 = 0; i1 < grids(1)->GetLength(); i1++)
          for (i2 = 0; i2 < grids(2)->GetLength(); i2++)
            for (i3 = 0; i3 < grids(3)->GetLength(); i3++)
              for (i4 = 0; i4 < grids(4)->GetLength(); i4++)
                for (i5 = 0; i5 < grids(5)->GetLength(); i5++)
                  f(grids(0)->Value(i0, i1, i2, i3,
                                    i4, i5),
                    grids(1)->Value(i0, i1, i2, i3,
                                    i4, i5),
                    grids(2)->Value(i0, i1, i2, i3,
                                    i4, i5),
                    grids(3)->Value(i0, i1, i2, i3,
                                    i4, i5),
                    grids(4)->Value(i0, i1, i2, i3,
                                    i4, i5),
                    grids(5)->Value(i0, i1, i2, i3,
                                    i4, i5),
                    grids_(0)->Value(i0, i1, i2, i3,
                                     i4, i5),
                    grids_(1)->Value(i0, i1, i2, i3,
                                     i4, i5),
                    grids_(2)->Value(i0, i1, i2, i3,
                                     i4, i5),
                    grids_(3)->Value(i0, i1, i2, i3,
                                     i4, i5),
                    grids_(4)->Value(i0, i1, i2, i3,
                                     i4, i5),
                    grids_(5)->Value(i0, i1, i2, i3,
                                     i4, i5));
    else if (N == 7)
      for (i0 = 0; i0 < grids(0)->GetLength(); i0++)
        for (i1 = 0; i1 < grids(1)->GetLength(); i1++)
          for (i2 = 0; i2 < grids(2)->GetLength(); i2++)
            for (i3 = 0; i3 < grids(3)->GetLength(); i3++)
              for (i4 = 0; i4 < grids(4)->GetLength(); i4++)
                for (i5 = 0; i5 < grids(5)->GetLength(); i5++)
                  for (i6 = 0; i6 < grids(6)->GetLength(); i6++)
                    f(grids(0)->Value(i0, i1, i2, i3,
                                      i4, i5, i6),
                      grids(1)->Value(i0, i1, i2, i3,
                                      i4, i5, i6),
                      grids(2)->Value(i0, i1, i2, i3,
                                      i4, i5, i6),
                      grids(3)->Value(i0, i1, i2, i3,
                                      i4, i5, i6),
                      grids(4)->Value(i0, i1, i2, i3,
                                      i4, i5, i6),
                      grids(5)->Value(i0, i1, i2, i3,
                                      i4, i5, i6),
                      grids(6)->Value(i0, i1, i2, i3,
                                      i4, i5, i6),
                      grids_(0)->Value(i0, i1, i2, i3,
                                       i4, i5, i6),
                      grids_(1)->Value(i0, i1, i2, i3,
                                       i4, i5, i6),
                      grids_(2)->Value(i0, i1, i2, i3,
                                       i4, i5, i6),
                      grids_(3)->Value(i0, i1, i2, i3,
                                       i4, i5, i6),
                      grids_(4)->Value(i0, i1, i2, i3,
                                       i4, i5, i6),
                      grids_(5)->Value(i0, i1, i2, i3,
                                       i4, i5, i6),
                      grids_(6)->Value(i0, i1, i2, i3,
                                       i4, i5, i6));
    else if (N == 8)
      for (i0 = 0; i0 < grids(0)->GetLength(); i0++)
        for (i1 = 0; i1 < grids(1)->GetLength(); i1++)
          for (i2 = 0; i2 < grids(2)->GetLength(); i2++)
            for (i3 = 0; i3 < grids(3)->GetLength(); i3++)
              for (i4 = 0; i4 < grids(4)->GetLength(); i4++)
                for (i5 = 0; i5 < grids(5)->GetLength(); i5++)
                  for (i6 = 0; i6 < grids(6)->GetLength(); i6++)
                    for (i7 = 0; i7 < grids(7)->GetLength(); i7++)
                      f(grids(0)->Value(i0, i1, i2, i3,
                                        i4, i5, i6, i7),
                        grids(1)->Value(i0, i1, i2, i3,
                                        i4, i5, i6, i7),
                        grids(2)->Value(i0, i1, i2, i3,
                                        i4, i5, i6, i7),
                        grids(3)->Value(i0, i1, i2, i3,
                                        i4, i5, i6, i7),
                        grids(4)->Value(i0, i1, i2, i3,
                                        i4, i5, i6, i7),
                        grids(5)->Value(i0, i1, i2, i3,
                                        i4, i5, i6, i7),
                        grids(6)->Value(i0, i1, i2, i3,
                                        i4, i5, i6, i7),
                        grids(7)->Value(i0, i1, i2, i3,
                                        i4, i5, i6, i7),
                        grids_(0)->Value(i0, i1, i2, i3,
                                         i4, i5, i6, i7),
                        grids_(1)->Value(i0, i1, i2, i3,
                                         i4, i5, i6, i7),
                        grids_(2)->Value(i0, i1, i2, i3,
                                         i4, i5, i6, i7),
                        grids_(3)->Value(i0, i1, i2, i3,
                                         i4, i5, i6, i7),
                        grids_(4)->Value(i0, i1, i2, i3,
                                         i4, i5, i6, i7),
                        grids_(5)->Value(i0, i1, i2, i3,
                                         i4, i5, i6, i7),
                        grids_(6)->Value(i0, i1, i2, i3,
                                         i4, i5, i6, i7),
                        grids_(7)->Value(i0, i1, i2, i3,
                                         i4, i5, i6, i7));
    else if (N == 9)
      for (i0 = 0; i0 < grids(0)->GetLength(); i0++)
        for (i1 = 0; i1 < grids(1)->GetLength(); i1++)
          for (i2 = 0; i2 < grids(2)->GetLength(); i2++)
            for (i3 = 0; i3 < grids(3)->GetLength(); i3++)
              for (i4 = 0; i4 < grids(4)->GetLength(); i4++)
                for (i5 = 0; i5 < grids(5)->GetLength(); i5++)
                  for (i6 = 0; i6 < grids(6)->GetLength(); i6++)
                    for (i7 = 0; i7 < grids(7)->GetLength(); i7++)
                      for (i8 = 0; i8 < grids(8)->GetLength(); i8++)
                        f(grids(0)->Value(i0, i1, i2, i3,
                                          i4, i5, i6, i7,
                                          i8),
                          grids(1)->Value(i0, i1, i2, i3,
                                          i4, i5, i6, i7,
                                          i8),
                          grids(2)->Value(i0, i1, i2, i3,
                                          i4, i5, i6, i7,
                                          i8),
                          grids(3)->Value(i0, i1, i2, i3,
                                          i4, i5, i6, i7,
                                          i8),
                          grids(4)->Value(i0, i1, i2, i3,
                                          i4, i5, i6, i7,
                                          i8),
                          grids(5)->Value(i0, i1, i2, i3,
                                          i4, i5, i6, i7,
                                          i8),
                          grids(6)->Value(i0, i1, i2, i3,
                                          i4, i5, i6, i7,
                                          i8),
                          grids(7)->Value(i0, i1, i2, i3,
                                          i4, i5, i6, i7,
                                          i8),
                          grids(8)->Value(i0, i1, i2, i3,
                                          i4, i5, i6, i7,
                                          i8),
                          grids_(0)->Value(i0, i1, i2, i3,
                                           i4, i5, i6, i7,
                                           i8),
                          grids_(1)->Value(i0, i1, i2, i3,
                                           i4, i5, i6, i7,
                                           i8),
                          grids_(2)->Value(i0, i1, i2, i3,
                                           i4, i5, i6, i7,
                                           i8),
                          grids_(3)->Value(i0, i1, i2, i3,
                                           i4, i5, i6, i7,
                                           i8),
                          grids_(4)->Value(i0, i1, i2, i3,
                                           i4, i5, i6, i7,
                                           i8),
                          grids_(5)->Value(i0, i1, i2, i3,
                                           i4, i5, i6, i7,
                                           i8),
                          grids_(6)->Value(i0, i1, i2, i3,
                                           i4, i5, i6, i7,
                                           i8),
                          grids_(7)->Value(i0, i1, i2, i3,
                                           i4, i5, i6, i7,
                                           i8),
                          grids_(8)->Value(i0, i1, i2, i3,
                                           i4, i5, i6, i7,
                                           i8));
    else if (N == 10)
      for (i0 = 0; i0 < grids(0)->GetLength(); i0++)
        for (i1 = 0; i1 < grids(1)->GetLength(); i1++)
          for (i2 = 0; i2 < grids(2)->GetLength(); i2++)
            for (i3 = 0; i3 < grids(3)->GetLength(); i3++)
              for (i4 = 0; i4 < grids(4)->GetLength(); i4++)
                for (i5 = 0; i5 < grids(5)->GetLength(); i5++)
                  for (i6 = 0; i6 < grids(6)->GetLength(); i6++)
                    for (i7 = 0; i7 < grids(7)->GetLength(); i7++)
                      for (i8 = 0; i8 < grids(8)->GetLength(); i8++)
                        for (i9 = 0; i9 < grids(9)->GetLength(); i9++)
                          f(grids(0)->Value(i0, i1, i2, i3,
                                            i4, i5, i6, i7,
                                            i8, i9),
                            grids(1)->Value(i0, i1, i2, i3,
                                            i4, i5, i6, i7,
                                            i8, i9),
                            grids(2)->Value(i0, i1, i2, i3,
                                            i4, i5, i6, i7,
                                            i8, i9),
                            grids(3)->Value(i0, i1, i2, i3,
                                            i4, i5, i6, i7,
                                            i8, i9),
                            grids(4)->Value(i0, i1, i2, i3,
                                            i4, i5, i6, i7,
                                            i8, i9),
                            grids(5)->Value(i0, i1, i2, i3,
                                            i4, i5, i6, i7,
                                            i8, i9),
                            grids(6)->Value(i0, i1, i2, i3,
                                            i4, i5, i6, i7,
                                            i8, i9),
                            grids(7)->Value(i0, i1, i2, i3,
                                            i4, i5, i6, i7,
                                            i8, i9),
                            grids(8)->Value(i0, i1, i2, i3,
                                            i4, i5, i6, i7,
                                            i8, i9),
                            grids(9)->Value(i0, i1, i2, i3,
                                            i4, i5, i6, i7,
                                            i8, i9),
                            grids_(0)->Value(i0, i1, i2, i3,
                                             i4, i5, i6, i7,
                                             i8, i9),
                            grids_(1)->Value(i0, i1, i2, i3,
                                             i4, i5, i6, i7,
                                             i8, i9),
                            grids_(2)->Value(i0, i1, i2, i3,
                                             i4, i5, i6, i7,
                                             i8, i9),
                            grids_(3)->Value(i0, i1, i2, i3,
                                             i4, i5, i6, i7,
                                             i8, i9),
                            grids_(4)->Value(i0, i1, i2, i3,
                                             i4, i5, i6, i7,
                                             i8, i9),
                            grids_(5)->Value(i0, i1, i2, i3,
                                             i4, i5, i6, i7,
                                             i8, i9),
                            grids_(6)->Value(i0, i1, i2, i3,
                                             i4, i5, i6, i7,
                                             i8, i9),
                            grids_(7)->Value(i0, i1, i2, i3,
                                             i4, i5, i6, i7,
                                             i8, i9),
                            grids_(8)->Value(i0, i1, i2, i3,
                                             i4, i5, i6, i7,
                                             i8, i9),
                            grids_(9)->Value(i0, i1, i2, i3,
                                             i4, i5, i6, i7,
                                             i8, i9));

    for (i = 0; i < N; i++)
      delete grids(i);

  }

  //! Coordinate transformation "in place".
  /*!
    Function f takes as inputs all coordinates and transforms those coordinates.
    This transformation is performed in place because function f works directly
    on its inputs.
    \param f coordinate transformation. It must be an instance of
    Function_Base<TG> or of a class derived from Function_Base<TG>.
    \warning One should use ChangeCoords instead of ChangeCoordsInPlace in order
    to save memory.
  */
  template<class T, int N, class TG>
  void Data<T, N, TG>::ChangeCoordsInPlace(Function_Base<TG>& f)
  {

    int i, j;
    int NbGrids = grids_.numElements();
    Array<Grid<TG>*, 1> grids_old(NbGrids),
      grids(NbGrids);

    Grid<TG>* old_grid;

    for (i = 0; i < NbGrids; i++)
      grids_old(i) = grids_(i)->Duplicate();

    for (i = 0; i < N; i++)
      {
        for (j = 0; j < NbGrids; j++)
          grids(j) = grids_old(j)->Duplicate();
        grids_(i)->ChangeCoordsInPlace(f, grids);
        for (j = 0; j < NbGrids; j++)
          delete grids(j);
      }

    for (int i = 0; i < N; i++)
      delete grids_old(i);

  }

  //! Displays data.
  template<class T, int N, class TG>
  void Data<T, N, TG>::Print()  const
  {
    cout << data_ << endl;
  }

  //! Displays minimum, mean, maximum and standard deviation of data.
  template<class T, int N, class TG>
  void Data<T, N, TG>::PrintInfo()  const
  {
    cout << "(min = " << this->GetMin()
         << ", max = " << this->GetMax()
         << ", mean = " << this->Mean()
         << ", std = " << this->StandardDeviation()
         << ")" << endl;
  }

  //! Returns a string with minimum, mean, maximum
  //! and standard deviation of data.
  template<class T, int N, class TG>
  string Data<T, N, TG>::InfoString()  const
  {
    return string("(min = ") + to_str(this->GetMin())
      + string(", max = ") + to_str(this->GetMax())
      + string(", mean = ") + to_str(this->Mean())
      + string(", std = ") + to_str(this->StandardDeviation())
      + string(")");
  }

  //! Clears grids.
  template<class T, int N, class TG>
  void Data<T, N, TG>::ClearGrids()
  {
    for (int i = 0; i < N; i++)
      if (grids_(i) != NULL)
        if (grids_(i)->GetPointers() == 1)
          delete grids_(i);
        else
          grids_(i)->SetPointers(grids_(i)->GetPointers() - 1);

    for (int i = 0; i < N; i++)
      grids_(i) = NULL;
  }

  //! For every grid, sets dimension to which the grid is related.
  template<class T, int N, class TG>
  void Data<T, N, TG>::SetVariables()
  {
    for (int i = 0; i < N; i++)
      grids_(i)->SetVariable(i);
  }

  //! For every grid, sets dimension to which the grid is related.
  template<class T, int N, class TG>
  void Data<T, N, TG>::InitData()
  {
#ifdef SELDONDATA_DEBUG_UNINITIALIZED_IS_NAN
    SetNaN();
#endif
  }

}  // namespace SeldonData.


#define FILE_SELDONDATA_DATA_CXX
#endif
