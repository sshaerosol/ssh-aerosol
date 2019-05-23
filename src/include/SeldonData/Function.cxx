// Copyright (C) 2003-2011, INRIA
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

#ifndef FILE_SELDONDATA_FUNCTION_CXX

#include "Function.hxx"
#include <iostream>
using std::cout;
using std::endl;

namespace SeldonData
{


  ///////////////////
  // FUNCTION_BASE //
  ///////////////////


  //! Constructor.
  template<class T, class TOut>
  Function_Base<T, TOut>::Function_Base()
  {
  }

  //! Destructor.
  template<class T, class TOut>
  Function_Base<T, TOut>::~Function_Base()
  {
  }

  //! Undefined function.
  template<class T, class TOut>
  TOut Function_Base<T, TOut>::operator()(T& i0)
  {
    throw Undefined("Function_Base<T>::operator(1)");
  }

  //! Undefined function.
  template<class T, class TOut>
  TOut Function_Base<T, TOut>::operator()(T& i0, T& i1)
  {
    throw Undefined("Function_Base<T>::operator(2)");
  }

  //! Undefined function.
  template<class T, class TOut>
  TOut Function_Base<T, TOut>::operator()(T& i0, T& i1,
                                          T& i2)
  {
    throw Undefined("Function_Base<T>::operator(3)");
  }

  //! Undefined function.
  template<class T, class TOut>
  TOut Function_Base<T, TOut>::operator()(T& i0, T& i1,
                                          T& i2, T& i3)
  {
    throw Undefined("Function_Base<T>::operator(4)");
  }

  //! Undefined function.
  template<class T, class TOut>
  TOut Function_Base<T, TOut>::operator()(T& i0, T& i1,
                                          T& i2, T& i3,
                                          T& i4)
  {
    throw Undefined("Function_Base<T>::operator(5)");
  }

  //! Undefined function.
  template<class T, class TOut>
  TOut Function_Base<T, TOut>::operator()(T& i0, T& i1,
                                          T& i2, T& i3,
                                          T& i4, T& i5)
  {
    throw Undefined("Function_Base<T>::operator(6)");
  }

  //! Undefined function.
  template<class T, class TOut>
  TOut Function_Base<T, TOut>::operator()(T& i0, T& i1,
                                          T& i2, T& i3,
                                          T& i4, T& i5,
                                          T& i6)
  {
    throw Undefined("Function_Base<T>::operator(7)");
  }

  //! Undefined function.
  template<class T, class TOut>
  TOut Function_Base<T, TOut>::operator()(T& i0, T& i1,
                                          T& i2, T& i3,
                                          T& i4, T& i5,
                                          T& i6, T& i7)
  {
    throw Undefined("Function_Base<T>::operator(8)");
  }

  //! Undefined function.
  template<class T, class TOut>
  TOut Function_Base<T, TOut>::operator()(T& i0, T& i1,
                                          T& i2, T& i3,
                                          T& i4, T& i5,
                                          T& i6, T& i7,
                                          T& i8)
  {
    throw Undefined("Function_Base<T>::operator(9)");
  }

  //! Undefined function.
  template<class T, class TOut>
  TOut Function_Base<T, TOut>::operator()(T& i0, T& i1,
                                          T& i2, T& i3,
                                          T& i4, T& i5,
                                          T& i6, T& i7,
                                          T& i8, T& i9)
  {
    throw Undefined("Function_Base<T>::operator(10)");
  }

  //! Undefined function.
  template<class T, class TOut>
  TOut Function_Base<T, TOut>::operator()(T& i0, T& i1,
                                          T& i2, T& i3,
                                          T& i4, T& i5,
                                          T& i6, T& i7,
                                          T& i8, T& i9,
                                          T& i10)
  {
    throw Undefined("Function_Base<T>::operator(11)");
  }

  //! Undefined function.
  template<class T, class TOut>
  TOut Function_Base<T, TOut>::operator()(T& i0, T& i1,
                                          T& i2, T& i3,
                                          T& i4, T& i5,
                                          T& i6, T& i7,
                                          T& i8, T& i9,
                                          T& i10, T& i11)
  {
    throw Undefined("Function_Base<T>::operator(12)");
  }

  //! Undefined function.
  template<class T, class TOut>
  TOut Function_Base<T, TOut>::operator()(T& i0, T& i1,
                                          T& i2, T& i3,
                                          T& i4, T& i5,
                                          T& i6, T& i7,
                                          T& i8, T& i9,
                                          T& i10, T& i11,
                                          T& i12)
  {
    throw Undefined("Function_Base<T>::operator(13)");
  }

  //! Undefined function.
  template<class T, class TOut>
  TOut Function_Base<T, TOut>::operator()(T& i0, T& i1,
                                          T& i2, T& i3,
                                          T& i4, T& i5,
                                          T& i6, T& i7,
                                          T& i8, T& i9,
                                          T& i10, T& i11,
                                          T& i12, T& i13)
  {
    throw Undefined("Function_Base<T>::operator(14)");
  }

  //! Undefined function.
  template<class T, class TOut>
  TOut Function_Base<T, TOut>::operator()(T& i0, T& i1,
                                          T& i2, T& i3,
                                          T& i4, T& i5,
                                          T& i6, T& i7,
                                          T& i8, T& i9,
                                          T& i10, T& i11,
                                          T& i12, T& i13,
                                          T& i14)
  {
    throw Undefined("Function_Base<T>::operator(15)");
  }

  //! Undefined function.
  template<class T, class TOut>
  TOut Function_Base<T, TOut>::operator()(T& i0, T& i1,
                                          T& i2, T& i3,
                                          T& i4, T& i5,
                                          T& i6, T& i7,
                                          T& i8, T& i9,
                                          T& i10, T& i11,
                                          T& i12, T& i13,
                                          T& i14, T& i15)
  {
    throw Undefined("Function_Base<T>::operator(16)");
  }

  //! Undefined function.
  template<class T, class TOut>
  TOut Function_Base<T, TOut>::operator()(T& i0, T& i1,
                                          T& i2, T& i3,
                                          T& i4, T& i5,
                                          T& i6, T& i7,
                                          T& i8, T& i9,
                                          T& i10, T& i11,
                                          T& i12, T& i13,
                                          T& i14, T& i15,
                                          T& i16)
  {
    throw Undefined("Function_Base<T>::operator(17)");
  }

  //! Undefined function.
  template<class T, class TOut>
  TOut Function_Base<T, TOut>::operator()(T& i0, T& i1,
                                          T& i2, T& i3,
                                          T& i4, T& i5,
                                          T& i6, T& i7,
                                          T& i8, T& i9,
                                          T& i10, T& i11,
                                          T& i12, T& i13,
                                          T& i14, T& i15,
                                          T& i16, T& i17)
  {
    throw Undefined("Function_Base<T>::operator(18)");
  }

  //! Undefined function.
  template<class T, class TOut>
  TOut Function_Base<T, TOut>::operator()(T& i0, T& i1,
                                          T& i2, T& i3,
                                          T& i4, T& i5,
                                          T& i6, T& i7,
                                          T& i8, T& i9,
                                          T& i10, T& i11,
                                          T& i12, T& i13,
                                          T& i14, T& i15,
                                          T& i16, T& i17,
                                          T& i18)
  {
    throw Undefined("Function_Base<T>::operator(19)");
  }

  //! Undefined function.
  template<class T, class TOut>
  TOut Function_Base<T, TOut>::operator()(T& i0, T& i1,
                                          T& i2, T& i3,
                                          T& i4, T& i5,
                                          T& i6, T& i7,
                                          T& i8, T& i9,
                                          T& i10, T& i11,
                                          T& i12, T& i13,
                                          T& i14, T& i15,
                                          T& i16, T& i17,
                                          T& i18, T& i19)
  {
    throw Undefined("Function_Base<T>::operator(20)");
  }


  /////////////////////
  // FUNCCOORDS_BASE //
  /////////////////////


  //! Constructor.
  template<class T>
  FuncCoords_Base<T>::FuncCoords_Base()
  {
  }

  //! Destructor.
  template<class T>
  FuncCoords_Base<T>::~FuncCoords_Base()
  {
  }

  //! Undefined function.
  template<class T>
  void FuncCoords_Base<T>::operator()(T& i0)
  {
    throw Undefined("FuncCoords_Base<T>::operator(1)");
  }

  //! Undefined function.
  template<class T>
  void FuncCoords_Base<T>::operator()(const T& i0, T& i1)
  {
    throw Undefined("FuncCoords_Base<T>::operator(2)");
  }

  //! Undefined function.
  template<class T>
  void FuncCoords_Base<T>::operator()(const T& i0, const T& i1,
                                      T& i2, T& i3)
  {
    throw Undefined("FuncCoords_Base<T>::operator(4)");
  }

  //! Undefined function.
  template<class T>
  void FuncCoords_Base<T>::operator()(const T& i0, const T& i1,
                                      const T& i2, T& i3,
                                      T& i4, T& i5)
  {
    throw Undefined("FuncCoords_Base<T>::operator(6)");
  }

  //! Undefined function.
  template<class T>
  void FuncCoords_Base<T>::operator()(const T& i0, const T& i1,
                                      const T& i2, const T& i3,
                                      T& i4, T& i5,
                                      T& i6, T& i7)
  {
    throw Undefined("FuncCoords_Base<T>::operator(8)");
  }

  //! Undefined function.
  template<class T>
  void FuncCoords_Base<T>::operator()(const T& i0, const T& i1,
                                      const T& i2, const T& i3,
                                      const T& i4, const T& i5,
                                      T& i6, T& i7,
                                      T& i8, T& i9)
  {
    throw Undefined("FuncCoords_Base<T>::operator(10)");
  }

  //! Undefined function.
  template<class T>
  void FuncCoords_Base<T>::operator()(const T& i0, const T& i1,
                                      const T& i2, const T& i3,
                                      const T& i4, const T& i5,
                                      T& i6, T& i7,
                                      T& i8, T& i9,
                                      T& i10, T& i11)
  {
    throw Undefined("FuncCoords_Base<T>::operator(12)");
  }

  //! Undefined function.
  template<class T>
  void FuncCoords_Base<T>::operator()(const T& i0, const T& i1,
                                      const T& i2, const T& i3,
                                      const T& i4, const T& i5,
                                      const T& i6, T& i7,
                                      T& i8, T& i9,
                                      T& i10, T& i11,
                                      T& i12, T& i13)
  {
    throw Undefined("FuncCoords_Base<T>::operator(14)");
  }

  //! Undefined function.
  template<class T>
  void FuncCoords_Base<T>::operator()(const T& i0, const T& i1,
                                      const T& i2, const T& i3,
                                      const T& i4, const T& i5,
                                      const T& i6, const T& i7,
                                      T& i8, T& i9,
                                      T& i10, T& i11,
                                      T& i12, T& i13,
                                      T& i14, T& i15)
  {
    throw Undefined("FuncCoords_Base<T>::operator(16)");
  }

  //! Undefined function.
  template<class T>
  void FuncCoords_Base<T>::operator()(const T& i0, const T& i1,
                                      const T& i2, const T& i3,
                                      const T& i4, const T& i5,
                                      const T& i6, const T& i7,
                                      const T& i8, T& i9,
                                      T& i10, T& i11,
                                      T& i12, T& i13,
                                      T& i14, T& i15,
                                      T& i16, T& i17)
  {
    throw Undefined("FuncCoords_Base<T>::operator(18)");
  }

  //! Undefined function.
  template<class T>
  void FuncCoords_Base<T>::operator()(const T& i0, const T& i1,
                                      const T& i2, const T& i3,
                                      const T& i4, const T& i5,
                                      const T& i6, const T& i7,
                                      const T& i8, const T& i9,
                                      T& i10, T& i11,
                                      T& i12, T& i13,
                                      T& i14, T& i15,
                                      T& i16, T& i17,
                                      T& i18, T& i19)
  {
    throw Undefined("FuncCoords_Base<T>::operator(20)");
  }


}  // namespace SeldonData.

#define FILE_SELDONDATA_FUNCTION_CXX
#endif
