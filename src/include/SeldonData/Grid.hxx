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

#ifndef FILE_SELDONDATA_GRID_HXX

namespace SeldonData
{


  /********/
  /* Grid */
  /********/

  //! Base class for grids.
  template <class T>
  class Grid
  {

    // typedef declarations.
  public:
    typedef T value_type;
    typedef T* pointer;
    typedef const T* const_pointer;
    typedef T& reference;
    typedef const T& const_reference;

  protected:
    //! Length of the grid.
    int length_;
    //! Dimension related to the grid.
    int variable_;

    //! When copying the grid,
    //! should grid values be duplicated?
    bool duplicate_;
    //! How many pointers to the current grid
    // are there?
    int pointers_;

    //! Zero.
    value_type zero_;

  public:

    // Constructors.

    Grid()  throw();
    Grid(int length, int variable = 0)  throw();
    Grid(const Grid& G)  throw();

    // Destructor.

    virtual ~Grid()  throw();

    // Methods.

    virtual Grid<T>& operator= (const Grid<T>&);

    virtual int GetLength() const;
    virtual int GetLength(int i) const;
    int GetVariable() const;
    virtual bool IsDependent(int i) const;

    virtual int GetNbElements() const;

    void SetVariable(int variable);

    void SetPointers(int pointers);
    int GetPointers() const;

    void SetDuplicate(bool duplicate);
    bool GetDuplicate() const;
    virtual Grid<T>* Duplicate() const;
    virtual Grid<T>* Copy();

    virtual reference operator()(int i);
    virtual value_type operator()(int i) const;

    virtual reference Value(int i0, int i1 = -1,
                            int i2 = -1, int i3 = -1,
                            int i4 = -1, int i5 = -1,
                            int i6 = -1, int i7 = -1,
                            int i8 = -1, int i9 = -1);

    virtual value_type Value(int i0, int i1 = -1,
                             int i2 = -1, int i3 = -1,
                             int i4 = -1, int i5 = -1,
                             int i6 = -1, int i7 = -1,
                             int i8 = -1, int i9 = -1) const;

    virtual void ChangeCoordsInPlace(Function_Base<T>& f,
                                     Array<Grid<T>*, 1> grids);

    virtual void Print() const;

  };


  /***************/
  /* RegularGrid */
  /***************/

  //! Regular grids.
  template<class T>
  class RegularGrid: public Grid<T>
  {

    // typedef declarations.
  public:
    typedef T value_type;
    typedef T* pointer;
    typedef const T* const_pointer;
    typedef T& reference;
    typedef const T& const_reference;

  protected:
    //! Grid values.
    Array<value_type, 1> values_;

  public:

    // Constructors.

    RegularGrid()  throw();
    RegularGrid(int length, int variable = 0)  throw();
    RegularGrid(value_type start, value_type inc, int length, int variable = 0)  throw();
    RegularGrid(const Array<value_type, 1>& values, int variable = 0)  throw();
    RegularGrid(const Grid<T>& G)  throw();
    RegularGrid(const RegularGrid<T>& G)  throw();

    // Destructor.

    ~RegularGrid()  throw();

    // Methods.

    RegularGrid<T>& operator= (const Grid<T>&);
    RegularGrid<T>& operator= (const RegularGrid<T>&);

    int GetLength() const;
    int GetLength(int i) const;

    int GetNbElements() const;
    Array<value_type, 1>& GetArray();
    const Array<value_type, 1>& GetArray() const;

    Grid<T>* Duplicate() const;
    Grid<T>* Copy();

    reference operator()(int i);
    value_type operator()(int i) const;

    reference Value(int i0, int i1 = -1,
                    int i2 = -1, int i3 = -1,
                    int i4 = -1, int i5 = -1,
                    int i6 = -1, int i7 = -1,
                    int i8 = -1, int i9 = -1);

    value_type Value(int i0, int i1 = -1,
                     int i2 = -1, int i3 = -1,
                     int i4 = -1, int i5 = -1,
                     int i6 = -1, int i7 = -1,
                     int i8 = -1, int i9 = -1) const;

    template <class F>
    void Apply(F& function);
    template <class T0, class F>
    void Apply(RegularGrid<T0>&, F& function);

    void Print() const;

  };


  /***************/
  /* GeneralGrid */
  /***************/

  //! General grids.
  template<class T, int n>
  class GeneralGrid: public Grid<T>
  {

    // typedef declarations.
  public:
    typedef T value_type;
    typedef T* pointer;
    typedef const T* const_pointer;
    typedef T& reference;
    typedef const T& const_reference;

  protected:
    //! Grid values.
    Array<value_type, n> values_;
    //! Dimension upon which the grid depends.
    Array<int, 1> dependencies_;

  public:

    // Constructors.

    GeneralGrid()  throw();
    GeneralGrid(Array<value_type, n>& values,
                int variable,
                const TinyVector<int, n>& dependencies)  throw();
    GeneralGrid(const TinyVector<int, n>& values_shape,
                int variable,
                const TinyVector<int, n>& dependencies)  throw();
    GeneralGrid(const GeneralGrid<T, n>& G)  throw();
    GeneralGrid(const Grid<T>& G)  throw();

    // Destructors.

    ~GeneralGrid()  throw();

    // Methods.

    GeneralGrid<T, n>& operator= (const Grid<T>&);
    GeneralGrid<T, n>& operator= (const GeneralGrid<T, n>&);

    int GetLength() const;
    int GetLength(int i) const;

    int GetNbElements() const;
    Array<value_type, n>& GetArray();
    const Array<value_type, n>& GetArray() const;
    Array<int, 1>& GetDependencies();
    const Array<int, 1>& GetDependencies() const;
    int GetMainVariable() const;
    bool IsDependent(int i) const;

    Grid<T>* Duplicate() const;
    Grid<T>* Copy();

    reference operator()(int i);
    value_type operator()(int i) const;

    reference Value(int i0, int i1 = -1,
                    int i2 = -1, int i3 = -1,
                    int i4 = -1, int i5 = -1,
                    int i6 = -1, int i7 = -1,
                    int i8 = -1, int i9 = -1);

    value_type Value(int i0, int i1 = -1,
                     int i2 = -1, int i3 = -1,
                     int i4 = -1, int i5 = -1,
                     int i6 = -1, int i7 = -1,
                     int i8 = -1, int i9 = -1) const;

    void ChangeCoordsInPlace(Function_Base<T>& f, Array<Grid<T>*, 1> grids);

    template <class F>
    void Apply(F& function);
    template <class T0, class F>
    void Apply(GeneralGrid<T0, n>&, F& function);

    void Print() const;

  private:
    //! Dimension of 'values_' corresponding to the dimension
    //! to which the grid is related (i.e. 'variable').
    int main_variable_;

  };


}  // namespace Data.


#define FILE_SELDONDATA_GRID_HXX
#endif

