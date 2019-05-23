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

#ifndef FILE_SELDONDATA_FUNCTION_HXX

#include <iostream>
using std::cout;
using std::endl;

namespace SeldonData
{


  ///////////////////
  // FUNCTION_BASE //
  ///////////////////


  //! Based class for functions.
  /*!
    A function f should be a class derived from Function_Base.
    Then, operator() should be overloaded to define the function f.
  */
  template<class T, class TOut = void>
  class Function_Base
  {

  public:

    //! Constructor.
    Function_Base();
    //! Destructor.
    virtual ~Function_Base();

    virtual TOut operator()(T& i0);
    virtual TOut operator()(T& i0, T& i1);
    virtual TOut operator()(T& i0, T& i1,
                            T& i2);
    virtual TOut operator()(T& i0, T& i1,
                            T& i2, T& i3);
    virtual TOut operator()(T& i0, T& i1,
                            T& i2, T& i3,
                            T& i4);
    virtual TOut operator()(T& i0, T& i1,
                            T& i2, T& i3,
                            T& i4, T& i5);
    virtual TOut operator()(T& i0, T& i1,
                            T& i2, T& i3,
                            T& i4, T& i5,
                            T& i6);
    virtual TOut operator()(T& i0, T& i1,
                            T& i2, T& i3,
                            T& i4, T& i5,
                            T& i6, T& i7);
    virtual TOut operator()(T& i0, T& i1,
                            T& i2, T& i3,
                            T& i4, T& i5,
                            T& i6, T& i7,
                            T& i8);
    virtual TOut operator()(T& i0, T& i1,
                            T& i2, T& i3,
                            T& i4, T& i5,
                            T& i6, T& i7,
                            T& i8, T& i9);
    virtual TOut operator()(T& i0, T& i1,
                            T& i2, T& i3,
                            T& i4, T& i5,
                            T& i6, T& i7,
                            T& i8, T& i9,
                            T& i10);
    virtual TOut operator()(T& i0, T& i1,
                            T& i2, T& i3,
                            T& i4, T& i5,
                            T& i6, T& i7,
                            T& i8, T& i9,
                            T& i10, T& i11);
    virtual TOut operator()(T& i0, T& i1,
                            T& i2, T& i3,
                            T& i4, T& i5,
                            T& i6, T& i7,
                            T& i8, T& i9,
                            T& i10, T& i11,
                            T& i12);
    virtual TOut operator()(T& i0, T& i1,
                            T& i2, T& i3,
                            T& i4, T& i5,
                            T& i6, T& i7,
                            T& i8, T& i9,
                            T& i10, T& i11,
                            T& i12, T& i13);
    virtual TOut operator()(T& i0, T& i1,
                            T& i2, T& i3,
                            T& i4, T& i5,
                            T& i6, T& i7,
                            T& i8, T& i9,
                            T& i10, T& i11,
                            T& i12, T& i13,
                            T& i14);
    virtual TOut operator()(T& i0, T& i1,
                            T& i2, T& i3,
                            T& i4, T& i5,
                            T& i6, T& i7,
                            T& i8, T& i9,
                            T& i10, T& i11,
                            T& i12, T& i13,
                            T& i14, T& i15);
    virtual TOut operator()(T& i0, T& i1,
                            T& i2, T& i3,
                            T& i4, T& i5,
                            T& i6, T& i7,
                            T& i8, T& i9,
                            T& i10, T& i11,
                            T& i12, T& i13,
                            T& i14, T& i15,
                            T& i16);
    virtual TOut operator()(T& i0, T& i1,
                            T& i2, T& i3,
                            T& i4, T& i5,
                            T& i6, T& i7,
                            T& i8, T& i9,
                            T& i10, T& i11,
                            T& i12, T& i13,
                            T& i14, T& i15,
                            T& i16, T& i17);
    virtual TOut operator()(T& i0, T& i1,
                            T& i2, T& i3,
                            T& i4, T& i5,
                            T& i6, T& i7,
                            T& i8, T& i9,
                            T& i10, T& i11,
                            T& i12, T& i13,
                            T& i14, T& i15,
                            T& i16, T& i17,
                            T& i18);
    virtual TOut operator()(T& i0, T& i1,
                            T& i2, T& i3,
                            T& i4, T& i5,
                            T& i6, T& i7,
                            T& i8, T& i9,
                            T& i10, T& i11,
                            T& i12, T& i13,
                            T& i14, T& i15,
                            T& i16, T& i17,
                            T& i18, T& i19);
  };



  /////////////////////
  // FUNCCOORDS_BASE //
  /////////////////////


  //! Based class for coordinates transformations.
  /*!
    A function f for coordinates transformation should be a class derived
    from FuncCoords_Base. Then, operator() should be overloaded
    to define the function f.
  */
  template<class T>
  class FuncCoords_Base
  {

  public:

    //! Constructor.
    FuncCoords_Base();
    //! Destructor.
    ~FuncCoords_Base();

    virtual void operator()(T& i0);
    virtual void operator()(const T& i0, T& i1);
    virtual void operator()(const T& i0, const T& i1,
                            T& i2, T& i3);
    virtual void operator()(const T& i0, const T& i1,
                            const T& i2, T& i3,
                            T& i4, T& i5);
    virtual void operator()(const T& i0, const T& i1,
                            const T& i2, const T& i3,
                            T& i4, T& i5,
                            T& i6, T& i7);
    virtual void operator()(const T& i0, const T& i1,
                            const T& i2, const T& i3,
                            const T& i4, const T& i5,
                            T& i6, T& i7,
                            T& i8, T& i9);
    virtual void operator()(const T& i0, const T& i1,
                            const T& i2, const T& i3,
                            const T& i4, const T& i5,
                            T& i6, T& i7,
                            T& i8, T& i9,
                            T& i10, T& i11);
    virtual void operator()(const T& i0, const T& i1,
                            const T& i2, const T& i3,
                            const T& i4, const T& i5,
                            const T& i6, T& i7,
                            T& i8, T& i9,
                            T& i10, T& i11,
                            T& i12, T& i13);
    virtual void operator()(const T& i0, const T& i1,
                            const T& i2, const T& i3,
                            const T& i4, const T& i5,
                            const T& i6, const T& i7,
                            T& i8, T& i9,
                            T& i10, T& i11,
                            T& i12, T& i13,
                            T& i14, T& i15);
    virtual void operator()(const T& i0, const T& i1,
                            const T& i2, const T& i3,
                            const T& i4, const T& i5,
                            const T& i6, const T& i7,
                            const T& i8, T& i9,
                            T& i10, T& i11,
                            T& i12, T& i13,
                            T& i14, T& i15,
                            T& i16, T& i17);
    virtual void operator()(const T& i0, const T& i1,
                            const T& i2, const T& i3,
                            const T& i4, const T& i5,
                            const T& i6, const T& i7,
                            const T& i8, const T& i9,
                            T& i10, T& i11,
                            T& i12, T& i13,
                            T& i14, T& i15,
                            T& i16, T& i17,
                            T& i18, T& i19);

  };


}  // namespace SeldonData.

#define FILE_SELDONDATA_FUNCTION_HXX
#endif
