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

#ifndef FILE_SELDONDATA_GRID_CXX

#include "Grid.hxx"

namespace SeldonData
{

  //////////
  // GRID //
  //////////

  //! Default constructor.
  /*!
    All is set to zero (including grid length).
  */
  template<class T>
  Grid<T>::Grid()  throw()
  {
    length_ = 0;
    variable_ = 0;
    zero_ = value_type(0);

    duplicate_ = true;
    pointers_ = 1;
  }

  //! Constructor.
  /*!
    \param length grid length.
    \param variable dimension number related to the grid.
  */
  template<class T>
  Grid<T>::Grid(int length, int variable)  throw()
  {
    length_ = length;
    variable_ = variable;
    zero_ = value_type(0);

    duplicate_ = true;
    pointers_ = 1;
  }

  //! Copy constructor.
  /*!
    \param G grid to be copied.
  */
  template<class T>
  Grid<T>::Grid(const Grid<T>& G)  throw()
  {
    length_ = G.GetLength();
    variable_ = G.GetVariable();
    zero_ = value_type(0);

    duplicate_ = G.GetDuplicate();
    pointers_ = 1;
  }

  //! Destructor.
  template<class T>
  Grid<T>::~Grid()  throw()
  {
  }

  //! Affectation operator.
  /*!
    \param G grid to be copied.
    \return A reference to the current grid.
  */
  template<class T>
  Grid<T>& Grid<T>::operator= (const Grid<T>& G)
  {
    length_ = G.GetLength();
    variable_ = G.GetVariable();
    zero_ = value_type(0);

    duplicate_ = G.GetDuplicate();
    pointers_ = 1;

    return *this;
  }

  //! Returns grid length.
  /*!
    \return Length of the grid.
  */
  template<class T>
  inline int Grid<T>::GetLength() const
  {
    return length_;
  }

  //! Returns grid length along dimension #i.
  /*!
    \param i dimension number.
    \return Length of the grid along the i-th dimension.
  */
  template<class T>
  inline int Grid<T>::GetLength(int i) const
  {
    if (i == variable_)
      return length_;
    else
      return 0;
  }

  //! Returns dimension number related to the grid.
  /*!
    \return Dimension number related to the grid.
  */
  template<class T>
  inline int Grid<T>::GetVariable() const
  {
    return variable_;
  }

  //! Returns whether the grid depends on a given dimension.
  /*!
    \param i dimension number.
    \return true if the grid depends on dimension #i, false otherwise.
    \exception SeldonData::WrongDim Dimension number i is out of range.
  */
  template<class T>
  inline bool Grid<T>::IsDependent(int i) const
  {

#ifdef SELDONDATA_DEBUG_CHECK_DIMENSIONS
    if ((i < 0) || (i > 9))
      throw WrongDim("Grid<T>::IsDependent", "Dimension number is " + i);
#endif

    return (i == variable_);

  }

  //! Returns the number of elements in the grid.
  /*!
    \return Zero.
  */
  template<class T>
  inline int Grid<T>::GetNbElements() const
  {
    return 0;
  }

  //! Sets the dimension to which the grid is related.
  /*! A grid is used to store coordinates along a given dimension.
    The function sets this dimension.
    \param variable dimension to which the grid is related.
  */
  template<class T>
  inline void Grid<T>::SetVariable(int variable)
  {
#ifndef SELDONDATA_ALLOW_VARIABLE_CHANGE
    if (variable_ != variable && pointers_ > 1)
      throw Error("Grid::SetVariable(int)",
                  string("Unable to change grid main variable from ")
                  + to_str(variable_) + string(" to ") + to_str(variable)
                  + string(". Option disabled."));
#endif
    variable_ = variable;
  }

  //! Sets the number of pointers that point to the current grid.
  /*!
    \param pointers the number of pointers that point to the current grid.
  */
  template<class T>
  void Grid<T>::SetPointers(int pointers)
  {
    pointers_ = pointers;
  }

  //! Returns the number of pointers that point to the current grid.
  /*!
    \return The number of pointers that point to the current grid.
  */
  template<class T>
  int Grid<T>::GetPointers() const
  {
    return pointers_;
  }

  //! Sets whether the grid should be duplicated in certain cases.
  /*!
    \param duplicate true if the grid should
    be duplicated, false otherwise.
  */
  template<class T>
  void Grid<T>::SetDuplicate(bool duplicate)
  {
    duplicate_ = duplicate;
  }

  //! Should the grid be duplicated?
  /*! When an instance of 'Data' is created, one may want to
    duplicate the grid or not (if two instances of 'Data' should
    share a given grid).
    \return true if the grid should be duplicated, false otherwise.
  */
  template<class T>
  bool Grid<T>::GetDuplicate() const
  {
    return duplicate_;
  }

  //! Duplicates the grid and returns a pointer to the new copy.
  /*! After duplication, no memory is shared with the new grid.
    \return A pointer to a copy of the current grid.
    \exception SeldonData::NoMemory no more memory is available; duplication is impossible.
  */
  template<class T>
  Grid<T>* Grid<T>::Duplicate() const
  {

    Grid<T>* G = new Grid<T>(*this);

#ifdef SELDONDATA_DEBUG_CHECK_MEMORY
    if (G == NULL)
      throw NoMemory("Grid<T>::Duplicate");
#endif

    return G;

  }

  //! Returns a pointer to a copy of the grid or to the grid itself.
  /*! After copy, no memory is shared with the new grid if 'duplicate_'
    is set to true. Otherwise, the new grid is the same as the current
    grid, and the returned pointer is the 'this'.
    \return A pointer to a copy of the current grid, or to the current grid.
    \exception SeldonData::NoMemory no more memory is available; duplication is impossible.
  */
  template<class T>
  Grid<T>* Grid<T>::Copy()
  {

    Grid<T>* G;

    if (duplicate_)
      G = new Grid<T>(*this);
    else
      {
        G = this;
        pointers_++;
      }

#ifdef SELDONDATA_DEBUG_CHECK_MEMORY
    if (G == NULL)
      throw NoMemory("Grid<T>::Copy");
#endif

    return G;

  }

  //! Returns a reference to the i-th element of the grid.
  /*!
    \param i index of the element to be returned.
    \return Zero.
    \exception SeldonData::WrongIndex index is out of range (i.e. is not 0).
  */
  template<class T>
  typename Grid<T>::reference
  inline Grid<T>::operator()(int i)
  {

#ifdef SELDONDATA_DEBUG_CHECK_INDICES
    if (i != 0) throw WrongIndex("reference Grid<T>::operator() (int)",
                                 "Index is " + to_str(i) + " and should be 0.");
#endif

    return zero_;
  }

  //! Returns i-th element of the grid.
  /*!
    \param i index of the element to be returned.
    \return Zero.
    \exception SeldonData::WrongIndex index is out of range (i.e. is not 0).
  */
  template<class T>
  typename Grid<T>::value_type
  inline Grid<T>::operator()(int i) const
  {

#ifdef SELDONDATA_DEBUG_CHECK_INDICES
    if (i != 0) throw WrongIndex("value_type Grid<T>::operator() (int)",
                                 "Index is " + to_str(i) + " and should be 0.");
#endif

    return zero_;
  }

  //! Returns a reference to an element of the grid.
  /*!
    \param i0 index along dimension #0.
    \param i1 index along dimension #1.
    \param i2 index along dimension #2.
    \param i3 index along dimension #3.
    \param i4 index along dimension #4.
    \param i5 index along dimension #5.
    \param i6 index along dimension #6.
    \param i7 index along dimension #7.
    \param i8 index along dimension #8.
    \param i9 index along dimension #9.
    \return Zero.
    \exception SeldonData::WrongIndex index is out of range (i.e. is not 0).
  */
  template<class T>
  typename Grid<T>::reference
  inline Grid<T>::Value(int i0, int i1,
                        int i2, int i3,
                        int i4, int i5,
                        int i6, int i7,
                        int i8, int i9)
  {

#ifdef SELDONDATA_DEBUG_CHECK_INDICES
    bool out;
    out = (((variable_ == 0) && (i0 != 0)) ||
           ((variable_ == 1) && (i1 != 0)) ||
           ((variable_ == 2) && (i2 != 0)) ||
           ((variable_ == 3) && (i3 != 0)) ||
           ((variable_ == 4) && (i4 != 0)) ||
           ((variable_ == 5) && (i5 != 0)) ||
           ((variable_ == 6) && (i6 != 0)) ||
           ((variable_ == 7) && (i7 != 0)) ||
           ((variable_ == 8) && (i8 != 0)) ||
           ((variable_ == 9) && (i9 != 0)));
    if (out) throw WrongIndex("reference Grid<T>::Value",
                              "Index should be 0.");
#endif

    return zero_;
  }

  //! Returns an element of the grid.
  /*!
    \param i0 index along dimension #0.
    \param i1 index along dimension #1.
    \param i2 index along dimension #2.
    \param i3 index along dimension #3.
    \param i4 index along dimension #4.
    \param i5 index along dimension #5.
    \param i6 index along dimension #6.
    \param i7 index along dimension #7.
    \param i8 index along dimension #8.
    \param i9 index along dimension #9.
    \return Zero.
    \exception SeldonData::WrongIndex index is out of range (i.e. is not 0).
  */
  template<class T>
  typename Grid<T>::value_type
  inline Grid<T>::Value(int i0, int i1,
                        int i2, int i3,
                        int i4, int i5,
                        int i6, int i7,
                        int i8, int i9) const
  {

#ifdef SELDONDATA_DEBUG_CHECK_INDICES
    bool out;
    out = (((variable_ == 0) && (i0 != 0)) ||
           ((variable_ == 1) && (i1 != 0)) ||
           ((variable_ == 2) && (i2 != 0)) ||
           ((variable_ == 3) && (i3 != 0)) ||
           ((variable_ == 4) && (i4 != 0)) ||
           ((variable_ == 5) && (i5 != 0)) ||
           ((variable_ == 6) && (i6 != 0)) ||
           ((variable_ == 7) && (i7 != 0)) ||
           ((variable_ == 8) && (i8 != 0)) ||
           ((variable_ == 9) && (i9 != 0)));
    if (out) throw WrongIndex("value_type Grid<T>::Value",
                              "Index should be 0 (empty grid).");
#endif

    return zero_;
  }

  //! Coordinate transformation "in place".
  /*!
    Function f takes as inputs all coordinates and transforms those coordinates.
    This transformation is performed in place because function f works directly
    on its inputs.
    \param f coordinate transformation. It must be an instance of
    Function_Base<TG> or of a class derived from Function_Base<TG>.
    \param grids array of pointers to other grids (along other directions).
  */
  template<class T>
  void Grid<T>::ChangeCoordsInPlace(Function_Base<T>& f,
                                    Array<Grid<T>*, 1> grids)
  {

    int N = grids.numElements();
    int i0, i1, i2, i3, i4, i5, i6, i7, i8, i9;

    // Sets the variable_-th grid to current grid.
    if ((grids(variable_) != NULL) && (grids(variable_) != this))
      delete grids(variable_);
    grids(variable_) = this;

    // Calls f at all needed points.
    // Needed points are defined thanks to grid dependencies.
    // Therefore, loops are performed over dimensions upon which the grid depends.
    if (N == 1)
      for (i0 = 0; i0 < int(!IsDependent(0)) +
             int(IsDependent(0))*grids(0)->GetLength(); i0++)
        f(grids(0)->Value(i0));
    else if (N == 2)
      for (i0 = 0; i0 < int(!IsDependent(0)) +
             int(IsDependent(0))*grids(0)->GetLength(); i0++)
        for (i1 = 0; i1 < int(!IsDependent(1)) +
               int(IsDependent(1))*grids(1)->GetLength(); i1++)
          f(grids(0)->Value(i0, i1),
            grids(1)->Value(i0, i1));
    else if (N == 3)
      for (i0 = 0; i0 < int(!IsDependent(0)) +
             int(IsDependent(0))*grids(0)->GetLength(); i0++)
        for (i1 = 0; i1 < int(!IsDependent(1)) +
               int(IsDependent(1))*grids(1)->GetLength(); i1++)
          for (i2 = 0; i2 < int(!IsDependent(2)) +
                 int(IsDependent(2))*grids(2)->GetLength(); i2++)
            f(grids(0)->Value(i0, i1, i2),
              grids(1)->Value(i0, i1, i2),
              grids(2)->Value(i0, i1, i2));
    else if (N == 4)
      for (i0 = 0; i0 < int(!IsDependent(0)) +
             int(IsDependent(0))*grids(0)->GetLength(); i0++)
        for (i1 = 0; i1 < int(!IsDependent(1)) +
               int(IsDependent(1))*grids(1)->GetLength(); i1++)
          for (i2 = 0; i2 < int(!IsDependent(2)) +
                 int(IsDependent(2))*grids(2)->GetLength(); i2++)
            for (i3 = 0; i3 < int(!IsDependent(3)) +
                   int(IsDependent(3))*grids(3)->GetLength(); i3++)
              f(grids(0)->Value(i0, i1, i2, i3),
                grids(1)->Value(i0, i1, i2, i3),
                grids(2)->Value(i0, i1, i2, i3),
                grids(3)->Value(i0, i1, i2, i3));
    else if (N == 5)
      for (i0 = 0; i0 < int(!IsDependent(0)) +
             int(IsDependent(0))*grids(0)->GetLength(); i0++)
        for (i1 = 0; i1 < int(!IsDependent(1)) +
               int(IsDependent(1))*grids(1)->GetLength(); i1++)
          for (i2 = 0; i2 < int(!IsDependent(2)) +
                 int(IsDependent(2))*grids(2)->GetLength(); i2++)
            for (i3 = 0; i3 < int(!IsDependent(3)) +
                   int(IsDependent(3))*grids(3)->GetLength(); i3++)
              for (i4 = 0; i4 < int(!IsDependent(4)) +
                     int(IsDependent(4))*grids(4)->GetLength(); i4++)
                f(grids(0)->Value(i0, i1, i2, i3, i4),
                  grids(1)->Value(i0, i1, i2, i3, i4),
                  grids(2)->Value(i0, i1, i2, i3, i4),
                  grids(3)->Value(i0, i1, i2, i3, i4),
                  grids(4)->Value(i0, i1, i2, i3, i4));
    else if (N == 6)
      for (i0 = 0; i0 < int(!IsDependent(0)) +
             int(IsDependent(0))*grids(0)->GetLength(); i0++)
        for (i1 = 0; i1 < int(!IsDependent(1)) +
               int(IsDependent(1))*grids(1)->GetLength(); i1++)
          for (i2 = 0; i2 < int(!IsDependent(2)) +
                 int(IsDependent(2))*grids(2)->GetLength(); i2++)
            for (i3 = 0; i3 < int(!IsDependent(3)) +
                   int(IsDependent(3))*grids(3)->GetLength(); i3++)
              for (i4 = 0; i4 < int(!IsDependent(4)) +
                     int(IsDependent(4))*grids(4)->GetLength(); i4++)
                for (i5 = 0; i5 < int(!IsDependent(5)) +
                       int(IsDependent(5))*grids(5)->GetLength(); i5++)
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
                                    i4, i5));
    else if (N == 7)
      for (i0 = 0; i0 < int(!IsDependent(0)) +
             int(IsDependent(0))*grids(0)->GetLength(); i0++)
        for (i1 = 0; i1 < int(!IsDependent(1)) +
               int(IsDependent(1))*grids(1)->GetLength(); i1++)
          for (i2 = 0; i2 < int(!IsDependent(2)) +
                 int(IsDependent(2))*grids(2)->GetLength(); i2++)
            for (i3 = 0; i3 < int(!IsDependent(3)) +
                   int(IsDependent(3))*grids(3)->GetLength(); i3++)
              for (i4 = 0; i4 < int(!IsDependent(4)) +
                     int(IsDependent(4))*grids(4)->GetLength(); i4++)
                for (i5 = 0; i5 < int(!IsDependent(5)) +
                       int(IsDependent(5))*grids(5)->GetLength(); i5++)
                  for (i6 = 0; i6 < int(!IsDependent(6)) +
                         int(IsDependent(6))*grids(6)->GetLength(); i6++)
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
                                      i4, i5, i6));
    else if (N == 8)
      for (i0 = 0; i0 < int(!IsDependent(0)) +
             int(IsDependent(0))*grids(0)->GetLength(); i0++)
        for (i1 = 0; i1 < int(!IsDependent(1)) +
               int(IsDependent(1))*grids(1)->GetLength(); i1++)
          for (i2 = 0; i2 < int(!IsDependent(2)) +
                 int(IsDependent(2))*grids(2)->GetLength(); i2++)
            for (i3 = 0; i3 < int(!IsDependent(3)) +
                   int(IsDependent(3))*grids(3)->GetLength(); i3++)
              for (i4 = 0; i4 < int(!IsDependent(4)) +
                     int(IsDependent(4))*grids(4)->GetLength(); i4++)
                for (i5 = 0; i5 < int(!IsDependent(5)) +
                       int(IsDependent(5))*grids(5)->GetLength(); i5++)
                  for (i6 = 0; i6 < int(!IsDependent(6)) +
                         int(IsDependent(6))*grids(6)->GetLength(); i6++)
                    for (i7 = 0; i7 < int(!IsDependent(7)) +
                           int(IsDependent(7))*grids(7)->GetLength(); i7++)
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
                                        i4, i5, i6, i7));
    else if (N == 9)
      for (i0 = 0; i0 < int(!IsDependent(0)) +
             int(IsDependent(0))*grids(0)->GetLength(); i0++)
        for (i1 = 0; i1 < int(!IsDependent(1)) +
               int(IsDependent(1))*grids(1)->GetLength(); i1++)
          for (i2 = 0; i2 < int(!IsDependent(2)) +
                 int(IsDependent(2))*grids(2)->GetLength(); i2++)
            for (i3 = 0; i3 < int(!IsDependent(3)) +
                   int(IsDependent(3))*grids(3)->GetLength(); i3++)
              for (i4 = 0; i4 < int(!IsDependent(4)) +
                     int(IsDependent(4))*grids(4)->GetLength(); i4++)
                for (i5 = 0; i5 < int(!IsDependent(5)) +
                       int(IsDependent(5))*grids(5)->GetLength(); i5++)
                  for (i6 = 0; i6 < int(!IsDependent(6)) +
                         int(IsDependent(6))*grids(6)->GetLength(); i6++)
                    for (i7 = 0; i7 < int(!IsDependent(7)) +
                           int(IsDependent(7))*grids(7)->GetLength(); i7++)
                      for (i8 = 0; i8 < int(!IsDependent(8)) +
                             int(IsDependent(8))*grids(8)->GetLength(); i8++)
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
                                          i8));
    else if (N == 10)
      for (i0 = 0; i0 < int(!IsDependent(0)) +
             int(IsDependent(0))*grids(0)->GetLength(); i0++)
        for (i1 = 0; i1 < int(!IsDependent(1)) +
               int(IsDependent(1))*grids(1)->GetLength(); i1++)
          for (i2 = 0; i2 < int(!IsDependent(2)) +
                 int(IsDependent(2))*grids(2)->GetLength(); i2++)
            for (i3 = 0; i3 < int(!IsDependent(3)) +
                   int(IsDependent(3))*grids(3)->GetLength(); i3++)
              for (i4 = 0; i4 < int(!IsDependent(4)) +
                     int(IsDependent(4))*grids(4)->GetLength(); i4++)
                for (i5 = 0; i5 < int(!IsDependent(5)) +
                       int(IsDependent(5))*grids(5)->GetLength(); i5++)
                  for (i6 = 0; i6 < int(!IsDependent(6)) +
                         int(IsDependent(6))*grids(6)->GetLength(); i6++)
                    for (i7 = 0; i7 < int(!IsDependent(7)) +
                           int(IsDependent(7))*grids(7)->GetLength(); i7++)
                      for (i8 = 0; i8 < int(!IsDependent(8)) +
                             int(IsDependent(8))*grids(8)->GetLength(); i8++)
                        for (i9 = 0; i9 < int(!IsDependent(9)) +
                               int(IsDependent(9))*grids(9)->GetLength(); i9++)
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
                                            i8, i9));

  }

  //! Displays grid values.
  /*!
    Displays "Empty grid.".
  */
  template<class T>
  void Grid<T>::Print() const
  {
    cout << "Empty grid." << endl;
  }


  /////////////////
  // REGULARGRID //
  /////////////////

  //! Default constructor.
  /*!
    All is set to zero (including grid length).
  */
  template<class T>
  RegularGrid<T>::RegularGrid()  throw()
  {
  }

  //! Constructor.
  /*!
    Constructs a regular grid of length 'length', related to dimension
    'dimension', with zero as initial point, and with one as increment.
    \param length grid length.
    \param variable dimension number related to the grid.
  */
  template<class T>
  RegularGrid<T>::RegularGrid(int length, int variable)  throw():
    Grid<T>(length, variable), values_(this->length_)
  {
    for (int i = 0; i < this->length_; i++)
      values_(i) = value_type(i);
  }

  //! Constructor.
  /*!
    Constructs a regular grid of length 'length', related to dimension
    'variable', with 'start' as initial point, and 'inc' as increment.
    \param start initial-point coordinate.
    \param inc increment or step between two points.
    \param length grid length.
    \param variable dimension number related to the grid.
  */
  template<class T>
  RegularGrid<T>::RegularGrid(typename RegularGrid<T>::value_type start,
                              typename RegularGrid<T>::value_type inc,
                              int length, int variable)  throw():
    Grid<T>(length, variable), values_(this->length_)
  {
    for (int i = 0; i < this->length_; i++)
      values_(i) = start + value_type(i) * inc;
  }

  //! Constructor.
  /*!
    Constructs a regular grid related to dimension 'variable'.
    Points coordinates are provided through 'values'.
    \param values grid values (i.e. points coordinates).
    \param variable dimension number related to the grid.
  */
  template<class T>
  RegularGrid<T>::RegularGrid(const Array<typename RegularGrid<T>::value_type, 1>& values,
                              int variable)  throw():
    Grid<T>(values.numElements(), variable),
    // values_ is a copy of 'values'. Function copy()
    // is used to avoid memory overlap.
    values_(values.copy())
  {
  }

  //! Copy constructor.
  /*!
    \param G grid to be copied.
  */
  template<class T>
  RegularGrid<T>::RegularGrid(const Grid<T>& G)  throw():
    Grid<T>(G), values_(G.GetLength())
  {
    for (int i = 0; i < G.GetLength(); i++)
      values_(i) = G(i);
  }

  //! Copy constructor for regular grids.
  /*!
    \param G regular grid to be copied.
  */
  template<class T>
  RegularGrid<T>::RegularGrid(const RegularGrid<T>& G)  throw():
    Grid<T>(G),
    // values_ is a copy of G::values_. Function copy()
    // is used to avoid memory overlap.
    values_((G.GetArray()).copy())
  {
  }

  //! Destructor.
  template<class T>
  RegularGrid<T>::~RegularGrid()  throw()
  {
  }

  //! Affectation operator.
  /*!
    \param G grid to be copied.
    \return A reference to the current grid.
  */
  template<class T>
  RegularGrid<T>& RegularGrid<T>::operator= (const Grid<T>& G)
  {
    this->length_ = G.GetLength();
    this->variable_ = G.GetVariable();

    return (*this);
  }

  //! Affectation operator for regular grids.
  /*!
    \param G regular grid to be copied.
    \return A reference to the current grid.
  */
  template<class T>
  RegularGrid<T>&
  RegularGrid<T>::operator= (const RegularGrid<T>& G)
  {
#ifdef SELDONDATA_DEBUG_CHECK_DIMENSIONS
    if (this->pointers_ != 1 && this->length_ != G.GetLength())
      throw WrongDim("RegularGrid<T, n>::operator= (const RegularGrid<T, n>&)",
                     "Grids must have the same size.");
#endif

    this->length_ = G.GetLength();
    this->variable_ = G.GetVariable();

    values_.resize(G.GetArray().shape());
    values_ = G.GetArray().copy();

    return (*this);
  }

  //! Returns grid length.
  /*!
    \return Length of the grid.
  */
  template<class T>
  inline int RegularGrid<T>::GetLength() const
  {
    return values_.extent(0);
  }

  //! Returns grid length along dimension #i.
  /*!
    \param i dimension number.
    \return Length of the grid along the i-th dimension.
  */
  template<class T>
  inline int RegularGrid<T>::GetLength(int i) const
  {
    if (i == this->variable_)
      return values_.extent(0);
    else
      return 0;
  }

  //! Returns the number of elements in the grid.
  /*!
    \return The number of elements in the grid.
  */
  template<class T>
  inline int RegularGrid<T>::GetNbElements() const
  {
    return values_.extent(0);
  }

  //! Returns a reference to the array storing points coordinates.
  /*!
    \return A reference to the array storing points coordinates.
  */
  template<class T>
  inline Array<typename RegularGrid<T>::value_type, 1>&
  RegularGrid<T>::GetArray()
  {
    return values_;
  }

  //! Returns a reference to the array storing points coordinates.
  /*!
    \return A reference to the array storing points coordinates.
  */
  template<class T>
  inline const Array<typename RegularGrid<T>::value_type, 1>&
  RegularGrid<T>::GetArray() const
  {
    return values_;
  }

  //! Duplicates the grid and returns a pointer to the new copy.
  /*! After duplication, no memory is shared with the new grid.
    \return A pointer to a copy of the current grid.
    \exception SeldonData::NoMemory no more memory is available; duplication is impossible.
  */
  template<class T>
  Grid<T>* RegularGrid<T>::Duplicate() const
  {
    Grid<T>* G = new RegularGrid<T>(*this);

#ifdef SELDONDATA_DEBUG_CHECK_MEMORY
    if (G == NULL)
      throw NoMemory("RegularGrid<T>::Duplicate");
#endif

    return G;
  }


  //! Returns a pointer to a copy of the grid or to the grid itself.
  /*! After copy, no memory is shared with the new grid if 'duplicate_'
    is set to true. Otherwise, the new grid is the same as the current
    grid, and the returned pointer is the 'this'.
    \return A pointer to a copy of the current grid, or to the current grid.
    \exception SeldonData::NoMemory no more memory is available; duplication is impossible.
  */
  template<class T>
  Grid<T>* RegularGrid<T>::Copy()
  {

    Grid<T>* G;

    if (this->duplicate_)
      G = new RegularGrid<T>(*this);
    else
      {
        G = this;
        this->pointers_++;
      }

#ifdef SELDONDATA_DEBUG_CHECK_MEMORY
    if (G == NULL)
      throw NoMemory("RegularGrid<T>::Copy");
#endif

    return G;

  }

  //! Returns a reference to the i-th element of the grid.
  /*!
    \param i index of the element to be returned.
    \return A reference to the i-th element of the grid.
    \exception SeldonData::WrongIndex index is out of range.
  */
  template<class T>
  inline typename RegularGrid<T>::reference
  RegularGrid<T>::operator()(int i)
  {

#ifdef SELDONDATA_DEBUG_CHECK_INDICES
    if ((i < 0) || (i >= values_.extent(0)))
      throw WrongIndex("reference RegularGrid<T>::operator() (int)",
                       "Index along dimension #" + to_str(this->variable_)
                       + " should be in [0, " + to_str(values_.extent(0) - 1)
                       + "], but is equal to " + to_str(i) + ".");
#endif

    return values_(i);
  }

  //! Returns the i-th element of the grid.
  /*!
    \param i index of the element to be returned.
    \return The i-th element of the grid.
    \exception SeldonData::WrongIndex index is out of range.
  */
  template<class T>
  inline typename RegularGrid<T>::value_type
  RegularGrid<T>::operator()(int i) const
  {

#ifdef SELDONDATA_DEBUG_CHECK_INDICES
    if ((i < 0) || (i >= values_.extent(0)))
      throw WrongIndex("value_type RegularGrid<T>::operator() (int)",
                       "Index along dimension #" + to_str(this->variable_)
                       + " should be in [0, " + to_str(values_.extent(0) - 1)
                       + "], but is equal to " + to_str(i) + ".");
#endif

    return values_(i);
  }

  //! Returns a reference to an element of the grid.
  /*!
    \param i0 index along dimension #0.
    \param i1 index along dimension #1.
    \param i2 index along dimension #2.
    \param i3 index along dimension #3.
    \param i4 index along dimension #4.
    \param i5 index along dimension #5.
    \param i6 index along dimension #6.
    \param i7 index along dimension #7.
    \param i8 index along dimension #8.
    \param i9 index along dimension #9.
    \return A reference to the element to be returned.
    \exception SeldonData::WrongIndex index is out of range.
  */
  template<class T>
  inline typename RegularGrid<T>::reference
  RegularGrid<T>::Value(int i0, int i1,
                        int i2, int i3,
                        int i4, int i5,
                        int i6, int i7,
                        int i8, int i9)
  {

#ifdef SELDONDATA_DEBUG_CHECK_INDICES
    bool out;
    int j = -1;

    if (this->variable_ == 0)
      j = i0;
    else if (this->variable_ == 1)
      j = i1;
    else if (this->variable_ == 2)
      j = i2;
    else if (this->variable_ == 3)
      j = i3;
    else if (this->variable_ == 4)
      j = i4;
    else if (this->variable_ == 5)
      j = i5;
    else if (this->variable_ == 6)
      j = i6;
    else if (this->variable_ == 7)
      j = i7;
    else if (this->variable_ == 8)
      j = i8;
    else if (this->variable_ == 9)
      j = i9;

    out = ((j < 0) || (j >= values_.extent(0)));

    if (out) throw WrongIndex("reference RegularGrid<T>::Value",
                              "Index along dimension #" + to_str(this->variable_)
                              + " should be in [0, " + to_str(values_.extent(0) - 1)
                              + "], but is equal to " + to_str(j) + ".");

    return values_(j);
#endif

    if (this->variable_ == 0)
      return values_(i0);
    else if (this->variable_ == 1)
      return values_(i1);
    else if (this->variable_ == 2)
      return values_(i2);
    else if (this->variable_ == 3)
      return values_(i3);
    else if (this->variable_ == 4)
      return values_(i4);
    else if (this->variable_ == 5)
      return values_(i5);
    else if (this->variable_ == 6)
      return values_(i6);
    else if (this->variable_ == 7)
      return values_(i7);
    else if (this->variable_ == 8)
      return values_(i8);
    else if (this->variable_ == 9)
      return values_(i9);
    return values_(i0);
  }

  //! Returns an element of the grid.
  /*!
    \param i0 index along dimension #0.
    \param i1 index along dimension #1.
    \param i2 index along dimension #2.
    \param i3 index along dimension #3.
    \param i4 index along dimension #4.
    \param i5 index along dimension #5.
    \param i6 index along dimension #6.
    \param i7 index along dimension #7.
    \param i8 index along dimension #8.
    \param i9 index along dimension #9.
    \return The element to be returned.
    \exception SeldonData::WrongIndex index is out of range.
  */
  template<class T>
  inline typename RegularGrid<T>::value_type
  RegularGrid<T>::Value(int i0, int i1,
                        int i2, int i3,
                        int i4, int i5,
                        int i6, int i7,
                        int i8, int i9) const
  {

#ifdef SELDONDATA_DEBUG_CHECK_INDICES
    bool out;
    int j = -1;

    if (this->variable_ == 0)
      j = i0;
    else if (this->variable_ == 1)
      j = i1;
    else if (this->variable_ == 2)
      j = i2;
    else if (this->variable_ == 3)
      j = i3;
    else if (this->variable_ == 4)
      j = i4;
    else if (this->variable_ == 5)
      j = i5;
    else if (this->variable_ == 6)
      j = i6;
    else if (this->variable_ == 7)
      j = i7;
    else if (this->variable_ == 8)
      j = i8;
    else if (this->variable_ == 9)
      j = i9;

    out = ((j < 0) || (j >= values_.extent(0)));

    if (out) throw WrongIndex("value_type RegularGrid<T>::Value",
                              "Index along dimension #" + to_str(this->variable_)
                              + " should be in [0, " + to_str(values_.extent(0) - 1)
                              + "], but is equal to " + to_str(j) + ".");

    return values_(j);
#endif

    if (this->variable_ == 0)
      return values_(i0);
    else if (this->variable_ == 1)
      return values_(i1);
    else if (this->variable_ == 2)
      return values_(i2);
    else if (this->variable_ == 3)
      return values_(i3);
    else if (this->variable_ == 4)
      return values_(i4);
    else if (this->variable_ == 5)
      return values_(i5);
    else if (this->variable_ == 6)
      return values_(i6);
    else if (this->variable_ == 7)
      return values_(i7);
    else if (this->variable_ == 8)
      return values_(i8);
    else if (this->variable_ == 9)
      return values_(i9);
    return values_(i0);
  }

  //! Applies a given function on all elements.
  /*!
    \param function function to be applied.
  */
  template<class T>
  template<class F>
  void RegularGrid<T>::Apply(F& function)
  {
    T* data = values_.data();
    int nb_elements = this->GetLength();

    for (int i = 0; i < nb_elements; i++)
      data[i] = function(data[i]);
  }

  //! Applies a given function on all elements of a 'RegularGrid' instance
  //! and put the result in the current instance.
  /*!
    \param G grid.
    \param function function to be applied on 'G'.
  */
  template<class T>
  template<class T0, class F>
  void RegularGrid<T>::Apply(RegularGrid<T0>& G, F& function)
  {
    T* data = values_.data();
    T0* data_in = G.GetArray().data();

    int nb_elements = this->Length();

#ifdef SELDONDATA_DEBUG_CHECK_DIMENSIONS
    if (nb_elements != G.GetLength())
      throw WrongDim("RegularGrid<T>::Apply(RegularGrid<T0>&, F& function)",
                     "Grid lengths differ.");
#endif

    for (int i = 0; i < nb_elements; i++)
      data[i] = function(data_in[i]);
  }

  //! Displays grid values.
  /*!
    Displays "Regular grid:" followed by the dimension and
    the grid values.
  */
  template<class T>
  void RegularGrid<T>::Print() const
  {
    cout << "Regular grid:" << endl;
    cout << values_ << endl;
  }


  /////////////////
  // GENERALGRID //
  /////////////////

  //! Default constructor.
  /*!
    All is set to zero (including grid length).
  */
  template<class T, int n>
  GeneralGrid<T, n>::GeneralGrid()  throw()
  {
  }

  // variable_: variable stockee par Grid, soit variable de la grille.
  // main_variable_: dimension du tableau qui correspond a variable_.
  //! Constructor.
  /*!
    Constructs a general grid related to dimension 'variable'.
    Points coordinates are provided through 'values'. The number of
    dimensions of 'values' is the number 'n' of dimension upon which
    the grid depends. 'dependencies' gives upon which dimensions
    the grid depends. Notice that 'dependencies' must contain
    'variable'.
    \param values points coordinates.
    \param variable dimension number related to the grid.
    \param dependencies dimensions upon which the grid depends.
  */
  template<class T, int n>
  GeneralGrid<T, n>::GeneralGrid(Array<typename GeneralGrid<T, n>::value_type, n>& values,
                                 int variable,
                                 const TinyVector<int, n>& dependencies)  throw():
    // values_ is a copy of 'values'. Function copy()
    // is used to avoid memory overlap.
    values_(values.copy()),
    dependencies_(n)
  {
    int i;

    for (i = 0; i < n; i++)
      dependencies_(i) = dependencies(i);

    this->variable_ = variable;

    // Searching for the array dimension-number related to
    // the main dependency of the grid.
    i = 0;
    // To add: error if variable_ is not in dependencies_.
    while ((i < n) && (dependencies_(i) != this->variable_))
      i++;
    main_variable_ = i;

    this->length_ = values_.extent(main_variable_);
  }

  //! Main constructor.
  /*!
    Constructs a general grid related to dimension 'variable'.
    The number of coordinates to be stored is provided through 'values_shape'.
    The number of elements of 'values' is the number 'n' of dimension
    upon which the grid depends. 'dependencies' gives upon which dimensions
    the grid depends. Notice that 'dependencies' must contain
    'variable'.
    \param values_shape array of extents along dimensions.
    \param variable dimension number related to the grid.
    \param dependencies dimensions upon which the grid depends.
  */
  template<class T, int n>
  GeneralGrid<T, n>::GeneralGrid(const TinyVector<int, n>& values_shape,
                                 int variable,
                                 const TinyVector<int, n>& dependencies)  throw():
    // Allocates values_ according to its shape (i.e. its dimensions).
    values_(values_shape),
    dependencies_(n)
  {
    int i;

    for (i = 0; i < n; i++)
      dependencies_(i) = dependencies(i);

    this->variable_ = variable;

    // Searching for the array dimension-number related to
    // the main dependency of the grid.
    i = 0;
    // To add: error if variable_ is not in dependencies_.
    while ((i < n) && (dependencies_(i) != this->variable_))
      i++;
    main_variable_ = i;

    this->length_ = values_.extent(main_variable_);
  }

  //! Copy constructor.
  /*!
    \param G grid to be copied.
  */
  template<class T, int n>
  GeneralGrid<T, n>::GeneralGrid(const Grid<T>& G)  throw(): Grid<T>(G)
  {
  }

  //! Copy constructor for general grids.
  /*!
    \param G general grid to be copied.
  */
  template<class T, int n>
  GeneralGrid<T, n>::GeneralGrid(const GeneralGrid<T, n>& G)  throw():
    Grid<T>(G),
    // values_ and dependencies_ are copies of G::values_ and
    // G::dependencies respectively. Function copy()
    // is used to avoid memory overlap.
    values_((G.GetArray()).copy()),
    dependencies_((G.GetDependencies()).copy())
  {
    main_variable_ = G.GetMainVariable();
  }

  //! Destructor.
  template<class T, int n>
  GeneralGrid<T, n>::~GeneralGrid()  throw()
  {
  }

  //! Affectation operator.
  /*!
    \param G grid to be copied.
    \return A reference to the current grid.
  */
  template<class T, int n>
  GeneralGrid<T, n>& GeneralGrid<T, n>::operator= (const Grid<T>& G)
  {
    this->length_ = G.GetLength();
    this->variable_ = G.GetVariable();

    return (*this);
  }

  //! Affectation operator for general grids.
  /*!
    \param G general grid to be copied.
    \return A reference to the current grid.
  */
  template<class T, int n>
  GeneralGrid<T, n>&
  GeneralGrid<T, n>::operator= (const GeneralGrid<T, n>& G)
  {
#ifdef SELDONDATA_DEBUG_CHECK_DIMENSIONS
    bool dim = true;
    for (int i = 0; i < n; i++)
      dim = dim && (values_.extent(i) == (G.GetArray()).extent(i));
    if (this->pointers_ != 1 && !dim)
      throw WrongDim("GeneralGrid<T, n>::operator= (const GeneralGrid<T, n>&)",
                     "Grids must have the same size.");
#endif

    this->length_ = G.GetLength();
    this->variable_ = G.GetVariable();

    values_.resize(G.GetArray().shape());
    values_ = G.GetArray().copy();

    dependencies_.resize(G.GetDependencies().shape());
    dependencies_ = G.GetDependencies().copy();

    main_variable_ = G.GetMainVariable();

    return (*this);
  }

  //! Returns grid length.
  /*!
    \return Length of the grid.
  */
  template<class T, int n>
  inline int GeneralGrid<T, n>::GetLength() const
  {
    return values_.extent(main_variable_);
  }

  //! Returns grid length along dimension #i.
  /*!
    \param i dimension number.
    \return Length of the grid along the i-th dimension.
  */
  template<class T, int n>
  inline int GeneralGrid<T, n>::GetLength(int i) const
  {
    int dim = 0;

    while ((dim < n) && (dependencies_(dim) != i))
      dim++;

    if (dim == n)
      return 0;
    else
      return values_.extent(dim);
  }

  //! Returns the number of elements in the grid.
  /*!
    \return The number of elements in the grid.
  */
  template<class T, int n>
  inline int GeneralGrid<T, n>::GetNbElements() const
  {
    return values_.numElements();
  }

  //! Returns a reference to the array storing points coordinates.
  /*!
    \return A reference to the array storing points coordinates.
  */
  template<class T, int n>
  inline Array<typename GeneralGrid<T, n>::value_type, n>&
  GeneralGrid<T, n>::GetArray()
  {
    return values_;
  }

  //! Returns a reference to the array storing points coordinates.
  /*!
    \return A reference to the array storing points coordinates.
  */
  template<class T, int n>
  inline const Array<typename GeneralGrid<T, n>::value_type, n>&
  GeneralGrid<T, n>::GetArray() const
  {
    return values_;
  }

  //! Returns a reference to dependencies array.
  /*!
    \return A reference to dependencies array.
  */
  template<class T, int n>
  inline Array<int, 1>&
  GeneralGrid<T, n>::GetDependencies()
  {
    return dependencies_;
  }

  //! Returns a reference to dependencies array.
  /*!
    \return A reference to dependencies array.
  */
  template<class T, int n>
  inline const Array<int, 1>&
  GeneralGrid<T, n>::GetDependencies() const
  {
    return dependencies_;
  }

  //! Returns the dimension of the array (storing points coordinates)
  //! corresponding to the dimension related to the grid.
  /*!
    For instance, let the grid be related to direction z. If coordinate along z
    depends upon x and z, then the array storing coordinateds is a (Nx, Nz) matrix.
    If the "0-th dimension" of the matrix corresponds to x, then the "1-st
    dimension" is the dimension of the array corresponding to the dimension related
    to the grid (i.e. z). So, this function will return 1.
    \return Main dimension of the array.
  */
  template<class T, int n>
  inline int GeneralGrid<T, n>::GetMainVariable() const
  {
    return main_variable_;
  }

  //! Returns whether the grid depends on a given dimension.
  /*!
    \param i dimension number.
    \return true if the grid depends on dimension #i, false otherwise.
    \exception SeldonData::WrongDim Dimension number i is out of range.
  */
  template<class T, int n>
  inline bool GeneralGrid<T, n>::IsDependent(int i) const
  {

#ifdef SELDONDATA_DEBUG_CHECK_DIMENSIONS
    if ((i < 0) || (i > 9))
      throw WrongDim("GeneralGrid<T, n>::IsDependent", "Dimension number is " + i);
#endif

    bool res = false;
    for (int j = 0; j < (int) dependencies_.numElements(); j++)
      res = res || (dependencies_(j) == i);

    return res;
  }

  //! Duplicates the grid and returns a pointer to the new copy.
  /*! After duplication, no memory is shared with the new grid.
    \return A pointer to a copy of the current grid.
    \exception SeldonData::NoMemory no more memory is available; duplication is impossible.
  */
  template<class T, int n>
  Grid<T>* GeneralGrid<T, n>::Duplicate() const
  {
    Grid<T>* G = new GeneralGrid<T, n>(*this);

#ifdef SELDONDATA_DEBUG_CHECK_MEMORY
    if (G == NULL)
      throw NoMemory("GeneralGrid<T, n>::Duplicate");
#endif

    return G;
  }


  //! Returns a pointer to a copy of the grid or to the grid itself.
  /*! After copy, no memory is shared with the new grid if 'duplicate_'
    is set to true. Otherwise, the new grid is the same as the current
    grid, and the returned pointer is the 'this'.
    \return A pointer to a copy of the current grid, or to the current grid.
    \exception SeldonData::NoMemory no more memory is available; duplication is impossible.
  */
  template<class T, int n>
  Grid<T>* GeneralGrid<T, n>::Copy()
  {

    Grid<T>* G;

    if (this->duplicate_)
      G = new GeneralGrid<T, n>(*this);
    else
      {
        G = this;
        this->pointers_++;
      }

#ifdef SELDONDATA_DEBUG_CHECK_MEMORY
    if (G == NULL)
      throw NoMemory("GeneralGrid<T, n>::Copy");
#endif

    return G;

  }

  //! Not defined.
  template<class T, int n>
  inline typename GeneralGrid<T, n>::reference
  GeneralGrid<T, n>::operator()(int i)
  {
    throw Undefined("reference GeneralGrid<T, n>::operator() (int)",
                    "It has no meaning for general grids.");
  }

  //! Not defined.
  template<class T, int n>
  inline typename GeneralGrid<T, n>::value_type
  GeneralGrid<T, n>::operator()(int i) const
  {
    throw Undefined("value_type GeneralGrid<T, n>::operator() (int)",
                    "It has no meaning for general grids.");
  }

  //! Returns a reference to an element of the grid.
  /*!
    \param i0 index along dimension #0.
    \param i1 index along dimension #1.
    \param i2 index along dimension #2.
    \param i3 index along dimension #3.
    \param i4 index along dimension #4.
    \param i5 index along dimension #5.
    \param i6 index along dimension #6.
    \param i7 index along dimension #7.
    \param i8 index along dimension #8.
    \param i9 index along dimension #9.
    \return A reference to the element to be returned.
    \exception SeldonData::WrongIndex an index is out of range.
  */
  template<class T, int n>
  inline typename GeneralGrid<T, n>::reference
  GeneralGrid<T, n>::Value(int i0, int i1,
                           int i2, int i3,
                           int i4, int i5,
                           int i6, int i7,
                           int i8, int i9)
  {

    Array<int, 1> Indices(10);

    Indices(0) = i0;
    Indices(1) = i1;
    Indices(2) = i2;
    Indices(3) = i3;
    Indices(4) = i4;
    Indices(5) = i5;
    Indices(6) = i6;
    Indices(7) = i7;
    Indices(8) = i8;
    Indices(9) = i9;

#ifdef SELDONDATA_DEBUG_CHECK_INDICES
    bool out = false;
    int j(0);
    for (int i = 0; i < n; i++)
      if (!out && ((Indices(dependencies_(i)) < 0) ||
                   (Indices(dependencies_(i)) >= values_.extent(i))))
        {
          out = true;
          j = i;
        }
    if (out) throw WrongIndex("reference GeneralGrid<T, n>::Value",
                              "Index along dimension #" + to_str(dependencies_(j))
                              + " should be in [0, "
                              + to_str(values_.extent(j) - 1) + "], but is equal to "
                              + to_str(Indices(dependencies_(j))) + ".");
#endif

    if (n == 1)
      return values_(Indices(dependencies_(0)));
    if (n == 2)
      return values_(Indices(dependencies_(0)),
                     Indices(dependencies_(1)));
    if (n == 3)
      return values_(Indices(dependencies_(0)),
                     Indices(dependencies_(1)),
                     Indices(dependencies_(2)));
    if (n == 4)
      return values_(Indices(dependencies_(0)),
                     Indices(dependencies_(1)),
                     Indices(dependencies_(2)),
                     Indices(dependencies_(3)));
    if (n == 5)
      return values_(Indices(dependencies_(0)),
                     Indices(dependencies_(1)),
                     Indices(dependencies_(2)),
                     Indices(dependencies_(3)),
                     Indices(dependencies_(4)));
    if (n == 6)
      return values_(Indices(dependencies_(0)),
                     Indices(dependencies_(1)),
                     Indices(dependencies_(2)),
                     Indices(dependencies_(3)),
                     Indices(dependencies_(4)),
                     Indices(dependencies_(5)));
    if (n == 7)
      return values_(Indices(dependencies_(0)),
                     Indices(dependencies_(1)),
                     Indices(dependencies_(2)),
                     Indices(dependencies_(3)),
                     Indices(dependencies_(4)),
                     Indices(dependencies_(5)),
                     Indices(dependencies_(6)));
    if (n == 8)
      return values_(Indices(dependencies_(0)),
                     Indices(dependencies_(1)),
                     Indices(dependencies_(2)),
                     Indices(dependencies_(3)),
                     Indices(dependencies_(4)),
                     Indices(dependencies_(5)),
                     Indices(dependencies_(6)),
                     Indices(dependencies_(7)));
    if (n == 9)
      return values_(Indices(dependencies_(0)),
                     Indices(dependencies_(1)),
                     Indices(dependencies_(2)),
                     Indices(dependencies_(3)),
                     Indices(dependencies_(4)),
                     Indices(dependencies_(5)),
                     Indices(dependencies_(6)),
                     Indices(dependencies_(7)),
                     Indices(dependencies_(8)));
    if (n == 10)
      return values_(Indices(dependencies_(0)),
                     Indices(dependencies_(1)),
                     Indices(dependencies_(2)),
                     Indices(dependencies_(3)),
                     Indices(dependencies_(4)),
                     Indices(dependencies_(5)),
                     Indices(dependencies_(6)),
                     Indices(dependencies_(7)),
                     Indices(dependencies_(8)),
                     Indices(dependencies_(9)));

  }

  //! Returns an element of the grid.
  /*!
    \param i0 index along dimension #0.
    \param i1 index along dimension #1.
    \param i2 index along dimension #2.
    \param i3 index along dimension #3.
    \param i4 index along dimension #4.
    \param i5 index along dimension #5.
    \param i6 index along dimension #6.
    \param i7 index along dimension #7.
    \param i8 index along dimension #8.
    \param i9 index along dimension #9.
    \return The element to be returned.
    \exception SeldonData::WrongIndex an index is out of range.
  */
  template<class T, int n>
  inline typename GeneralGrid<T, n>::value_type
  GeneralGrid<T, n>::Value(int i0, int i1,
                           int i2, int i3,
                           int i4, int i5,
                           int i6, int i7,
                           int i8, int i9) const
  {
    Array<int, 1> Indices(10);

    Indices(0) = i0;
    Indices(1) = i1;
    Indices(2) = i2;
    Indices(3) = i3;
    Indices(4) = i4;
    Indices(5) = i5;
    Indices(6) = i6;
    Indices(7) = i7;
    Indices(8) = i8;
    Indices(9) = i9;

#ifdef SELDONDATA_DEBUG_CHECK_INDICES
    bool out = false;
    int j(0);
    for (int i = 0; i < n; i++)
      if (!out && ((Indices(dependencies_(i)) < 0) ||
                   (Indices(dependencies_(i)) >= values_.extent(i))))
        {
          out = true;
          j = i;
        }
    if (out) throw WrongIndex("value_type GeneralGrid<T, n>::Value",
                              "Index along dimension #" + to_str(dependencies_(j))
                              + " should be in [0, "
                              + to_str(values_.extent(j) - 1) + "], but is equal to "
                              + to_str(Indices(dependencies_(j))) + ".");
#endif

    if (n == 1)
      return values_(Indices(dependencies_(0)));
    if (n == 2)
      return values_(Indices(dependencies_(0)),
                     Indices(dependencies_(1)));
    if (n == 3)
      return values_(Indices(dependencies_(0)),
                     Indices(dependencies_(1)),
                     Indices(dependencies_(2)));
    if (n == 4)
      return values_(Indices(dependencies_(0)),
                     Indices(dependencies_(1)),
                     Indices(dependencies_(2)),
                     Indices(dependencies_(3)));
    if (n == 5)
      return values_(Indices(dependencies_(0)),
                     Indices(dependencies_(1)),
                     Indices(dependencies_(2)),
                     Indices(dependencies_(3)),
                     Indices(dependencies_(4)));
    if (n == 6)
      return values_(Indices(dependencies_(0)),
                     Indices(dependencies_(1)),
                     Indices(dependencies_(2)),
                     Indices(dependencies_(3)),
                     Indices(dependencies_(4)),
                     Indices(dependencies_(5)));
    if (n == 7)
      return values_(Indices(dependencies_(0)),
                     Indices(dependencies_(1)),
                     Indices(dependencies_(2)),
                     Indices(dependencies_(3)),
                     Indices(dependencies_(4)),
                     Indices(dependencies_(5)),
                     Indices(dependencies_(6)));
    if (n == 8)
      return values_(Indices(dependencies_(0)),
                     Indices(dependencies_(1)),
                     Indices(dependencies_(2)),
                     Indices(dependencies_(3)),
                     Indices(dependencies_(4)),
                     Indices(dependencies_(5)),
                     Indices(dependencies_(6)),
                     Indices(dependencies_(7)));
    if (n == 9)
      return values_(Indices(dependencies_(0)),
                     Indices(dependencies_(1)),
                     Indices(dependencies_(2)),
                     Indices(dependencies_(3)),
                     Indices(dependencies_(4)),
                     Indices(dependencies_(5)),
                     Indices(dependencies_(6)),
                     Indices(dependencies_(7)),
                     Indices(dependencies_(8)));
    if (n == 10)
      return values_(Indices(dependencies_(0)),
                     Indices(dependencies_(1)),
                     Indices(dependencies_(2)),
                     Indices(dependencies_(3)),
                     Indices(dependencies_(4)),
                     Indices(dependencies_(5)),
                     Indices(dependencies_(6)),
                     Indices(dependencies_(7)),
                     Indices(dependencies_(8)),
                     Indices(dependencies_(9)));

  }

  //! Not defined.
  template<class T, int n>
  void GeneralGrid<T, n>::ChangeCoordsInPlace(Function_Base<T>& f,
                                              Array<Grid<T>*, 1> grids_old)
  {
    throw Undefined("GeneralGrid<T, n>::ChangeCoordsInPlace",
                    "Not defined for general grids.");
  }

  //! Applies a given function on all elements.
  /*!
    \param function function to be applied.
  */
  template<class T, int n>
  template<class F>
  void GeneralGrid<T, n>::Apply(F& function)
  {
    T* data = values_.data();
    int nb_elements = values_.numElements();

    for (int i = 0; i < nb_elements; i++)
      data[i] = function(data[i]);
  }

  //! Applies a given function on all elements of a 'GeneralGrid' instance
  //! and put the result in the current instance.
  /*!
    \param G grid.
    \param function function to be applied on 'G'.
  */
  template<class T, int n>
  template<class T0, class F>
  void GeneralGrid<T, n>::Apply(GeneralGrid<T0, n>& G, F& function)
  {
    T* data = values_.data();
    T0* data_in = G.GetArray().data();

    int nb_elements = values_.numElements();

#ifdef SELDONDATA_DEBUG_CHECK_DIMENSIONS
    if (nb_elements != G.GetArray().numElements())
      throw WrongDim("GeneralGrid<T, " + to_str(n) + ">::Apply(GeneralGrid<T0, "
                     + to_str(n) + ">&, F& function)",
                     "Grid lengths differ.");
#endif

    for (int i = 0; i < nb_elements; i++)
      data[i] = function(data_in[i]);
  }

  //! Displays grid values.
  /*!
    Displays "General grid:" followed by dimensions and
    the grid values.
  */
  template<class T, int n>
  void GeneralGrid<T, n>::Print() const
  {
    cout << "General Grid:" << endl;
    cout << values_ << endl;
  }

}  // namespace SeldonData.


#define FILE_SELDONDATA_GRID_CXX
#endif
