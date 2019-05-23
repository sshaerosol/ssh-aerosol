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

#ifndef FILE_SELDONDATA_DATA_HXX

#include <iostream>
using std::cout;
using std::endl;

namespace SeldonData
{

  //! Data class.
  template<class T, int N, class TG = T>
  class Data
  {

  protected:
    //! Array storing data.
    Array<T, N> data_;
    //! Array of pointers to grids.
    Array<Grid<TG>*, 1> grids_;

  public:

    // Constructors.

    Data()  throw();
    Data(Grid<TG>& G0)  throw();
    Data(Grid<TG>& G0, Grid<TG>& G1)  throw();
    Data(Grid<TG>& G0, Grid<TG>& G1,
         Grid<TG>& G2)  throw();
    Data(Grid<TG>& G0, Grid<TG>& G1,
         Grid<TG>& G2, Grid<TG>& G3)  throw();
    Data(Grid<TG>& G0, Grid<TG>& G1,
         Grid<TG>& G2, Grid<TG>& G3,
         Grid<TG>& G4)  throw();
    Data(Grid<TG>& G0, Grid<TG>& G1,
         Grid<TG>& G2, Grid<TG>& G3,
         Grid<TG>& G4, Grid<TG>& G5)  throw();
    Data(Grid<TG>& G0, Grid<TG>& G1,
         Grid<TG>& G2, Grid<TG>& G3,
         Grid<TG>& G4, Grid<TG>& G5,
         Grid<TG>& G6)  throw();
    Data(Grid<TG>& G0, Grid<TG>& G1,
         Grid<TG>& G2, Grid<TG>& G3,
         Grid<TG>& G4, Grid<TG>& G5,
         Grid<TG>& G6, Grid<TG>& G7)  throw();
    Data(Grid<TG>& G0, Grid<TG>& G1,
         Grid<TG>& G2, Grid<TG>& G3,
         Grid<TG>& G4, Grid<TG>& G5,
         Grid<TG>& G6, Grid<TG>& G7,
         Grid<TG>& G8)  throw();
    Data(Grid<TG>& G0, Grid<TG>& G1,
         Grid<TG>& G2, Grid<TG>& G3,
         Grid<TG>& G4, Grid<TG>& G5,
         Grid<TG>& G6, Grid<TG>& G7,
         Grid<TG>& G8, Grid<TG>& G9)  throw();
    Data(int N0)  throw();
    Data(int N0, int N1)  throw();
    Data(int N0, int N1,
         int N2)  throw();
    Data(int N0, int N1,
         int N2, int N3)  throw();
    Data(int N0, int N1,
         int N2, int N3,
         int N4)  throw();
    Data(int N0, int N1,
         int N2, int N3,
         int N4, int N5)  throw();
    Data(int N0, int N1,
         int N2, int N3,
         int N4, int N5,
         int N6)  throw();
    Data(int N0, int N1,
         int N2, int N3,
         int N4, int N5,
         int N6, int N7)  throw();
    Data(int N0, int N1,
         int N2, int N3,
         int N4, int N5,
         int N6, int N7,
         int N8)  throw();
    Data(int N0, int N1,
         int N2, int N3,
         int N4, int N5,
         int N6, int N7,
         int N8, int N9)  throw();
    Data(const TinyVector<int, N>& shape) throw();
    Data(T* data, const TinyVector<int, N>& shape,
         preexistingMemoryPolicy policy = neverDeleteData) throw();
    template <class T0>
    Data(Data<T0, N, TG>& data)  throw();

    // Destructor.

    ~Data()  throw();

    // Methods.

    T& operator()(int i0);
    T& operator()(int i0, int i1);
    T& operator()(int i0, int i1, int i2);
    T& operator()(int i0, int i1, int i2,
                  int i3);
    T& operator()(int i0, int i1, int i2,
                  int i3, int i4);
    T& operator()(int i0, int i1, int i2,
                  int i3, int i4, int i5);
    T& operator()(int i0, int i1, int i2,
                  int i3, int i4, int i5,
                  int i6);
    T& operator()(int i0, int i1, int i2,
                  int i3, int i4, int i5,
                  int i6, int i7);
    T& operator()(int i0, int i1, int i2,
                  int i3, int i4, int i5,
                  int i6, int i7, int i8);
    T& operator()(int i0, int i1, int i2,
                  int i3, int i4, int i5,
                  int i6, int i7, int i8,
                  int i9);

    T operator()(int i0) const;
    T operator()(int i0, int i1) const;
    T operator()(int i0, int i1, int i2) const;
    T operator()(int i0, int i1, int i2,
                 int i3) const;
    T operator()(int i0, int i1, int i2,
                 int i3, int i4) const;
    T operator()(int i0, int i1, int i2,
                 int i3, int i4, int i5) const;
    T operator()(int i0, int i1, int i2,
                 int i3, int i4, int i5,
                 int i6) const;
    T operator()(int i0, int i1, int i2,
                 int i3, int i4, int i5,
                 int i6, int i7) const;
    T operator()(int i0, int i1, int i2,
                 int i3, int i4, int i5,
                 int i6, int i7, int i8) const;
    T operator()(int i0, int i1, int i2,
                 int i3, int i4, int i5,
                 int i6, int i7, int i8,
                 int i9) const;

    T& operator()(const Array<int, 1>& indices);

    T& Value(int i0, int i1 = -1,
             int i2 = -1, int i3 = -1,
             int i4 = -1, int i5 = -1,
             int i6 = -1, int i7 = -1,
             int i8 = -1, int i9 = -1);

    T Value(int i0, int i1 = -1,
            int i2 = -1, int i3 = -1,
            int i4 = -1, int i5 = -1,
            int i6 = -1, int i7 = -1,
            int i8 = -1, int i9 = -1) const;

    template <class T0>
    void Copy(Data<T0, N, TG>& data);
    template <class T0>
    void ReferenceCopy(Data<T0, N, TG>& data);

    int GetNbElements();
    int GetNbDim();
    int GetLength(int dim) const;

    Grid<TG>* GetGrid(int i);
    Grid<TG>& operator [](int i);
    const Grid<TG>& operator [](int i) const;
    Array<Grid<TG>*, 1>& GetGrids();

    Array<T, N>& GetArray();
    const Array<T, N>& GetArray() const;
    Array<T, N>& operator()();
    T* GetData();
    const T* GetData() const;

    template<class DTG, class R0>
    void SubData(Data<T, 1, DTG>&, R0 r0);
    template<class DTG, class R0, class R1>
    void SubData(Data<T, 2, DTG>&, R0 r0, R1 r1);
    template<class DTG, class R0, class R1, class R2>
    void SubData(Data<T, 3, DTG>&, R0 r0, R1 r1, R2 r2);
    template < class DTG, class R0, class R1, class R2,
               class R3 >
    void SubData(Data<T, 4, DTG>&, R0 r0, R1 r1, R2 r2,
                 R3 r3);
    template < class DTG, class R0, class R1, class R2,
               class R3, class R4 >
    void SubData(Data<T, 5, DTG>&, R0 r0, R1 r1, R2 r2,
                 R3 r3, R4 r4);
    template < class DTG, class R0, class R1, class R2,
               class R3, class R4, class R5 >
    void SubData(Data<T, 6, DTG>&, R0 r0, R1 r1, R2 r2,
                 R3 r3, R4 r4, R5 r5);
    template < class DTG, class R0, class R1, class R2,
               class R3, class R4, class R5,
               class R6 >
    void SubData(Data<T, 7, DTG>&, R0 r0, R1 r1, R2 r2,
                 R3 r3, R4 r4, R5 r5, R6 r6);
    template < class DTG, class R0, class R1, class R2,
               class R3, class R4, class R5,
               class R6, class R7 >
    void SubData(Data<T, 8, DTG>&, R0 r0, R1 r1, R2 r2,
                 R3 r3, R4 r4, R5 r5, R6 r6,
                 R7 r7);
    template < class DTG, class R0, class R1, class R2,
               class R3, class R4, class R5,
               class R6, class R7, class R8 >
    void SubData(Data<T, 9, DTG>&, R0 r0, R1 r1, R2 r2,
                 R3 r3, R4 r4, R5 r5, R6 r6,
                 R7 r7, R8 r8);
    template < class DTG, class R0, class R1, class R2,
               class R3, class R4, class R5,
               class R6, class R7, class R8, class R9 >
    void SubData(Data<T, 10, DTG>&, R0 r0, R1 r1, R2 r2,
                 R3 r3, R4 r4, R5 r5, R6 r6,
                 R7 r7, R8 r8, R9 r9);

    void ResizeGrid();
    void ResizeGrid(Grid<TG>& G0);
    void ResizeGrid(Grid<TG>& G0, Grid<TG>& G1);
    void ResizeGrid(Grid<TG>& G0, Grid<TG>& G1,
                    Grid<TG>& G2);
    void ResizeGrid(Grid<TG>& G0, Grid<TG>& G1,
                    Grid<TG>& G2, Grid<TG>& G3);
    void ResizeGrid(Grid<TG>& G0, Grid<TG>& G1,
                    Grid<TG>& G2, Grid<TG>& G3,
                    Grid<TG>& G4);
    void ResizeGrid(Grid<TG>& G0, Grid<TG>& G1,
                    Grid<TG>& G2, Grid<TG>& G3,
                    Grid<TG>& G4, Grid<TG>& G5);
    void ResizeGrid(Grid<TG>& G0, Grid<TG>& G1,
                    Grid<TG>& G2, Grid<TG>& G3,
                    Grid<TG>& G4, Grid<TG>& G5,
                    Grid<TG>& G6);
    void ResizeGrid(Grid<TG>& G0, Grid<TG>& G1,
                    Grid<TG>& G2, Grid<TG>& G3,
                    Grid<TG>& G4, Grid<TG>& G5,
                    Grid<TG>& G6, Grid<TG>& G7);
    void ResizeGrid(Grid<TG>& G0, Grid<TG>& G1,
                    Grid<TG>& G2, Grid<TG>& G3,
                    Grid<TG>& G4, Grid<TG>& G5,
                    Grid<TG>& G6, Grid<TG>& G7,
                    Grid<TG>& G8);
    void ResizeGrid(Grid<TG>& G0, Grid<TG>& G1,
                    Grid<TG>& G2, Grid<TG>& G3,
                    Grid<TG>& G4, Grid<TG>& G5,
                    Grid<TG>& G6, Grid<TG>& G7,
                    Grid<TG>& G8, Grid<TG>& G9);

    void ResizeGrid(int N0);
    void ResizeGrid(int N0, int N1);
    void ResizeGrid(int N0, int N1,
                    int N2);
    void ResizeGrid(int N0, int N1,
                    int N2, int N3);
    void ResizeGrid(int N0, int N1,
                    int N2, int N3,
                    int N4);
    void ResizeGrid(int N0, int N1,
                    int N2, int N3,
                    int N4, int N5);
    void ResizeGrid(int N0, int N1,
                    int N2, int N3,
                    int N4, int N5,
                    int N6);
    void ResizeGrid(int N0, int N1,
                    int N2, int N3,
                    int N4, int N5,
                    int N6, int N7);
    void ResizeGrid(int N0, int N1,
                    int N2, int N3,
                    int N4, int N5,
                    int N6, int N7,
                    int N8);
    void ResizeGrid(int N0, int N1,
                    int N2, int N3,
                    int N4, int N5,
                    int N6, int N7,
                    int N8, int N9);

    void Resize();
    void Resize(Grid<TG>& G0);
    void Resize(Grid<TG>& G0, Grid<TG>& G1);
    void Resize(Grid<TG>& G0, Grid<TG>& G1,
                Grid<TG>& G2);
    void Resize(Grid<TG>& G0, Grid<TG>& G1,
                Grid<TG>& G2, Grid<TG>& G3);
    void Resize(Grid<TG>& G0, Grid<TG>& G1,
                Grid<TG>& G2, Grid<TG>& G3,
                Grid<TG>& G4);
    void Resize(Grid<TG>& G0, Grid<TG>& G1,
                Grid<TG>& G2, Grid<TG>& G3,
                Grid<TG>& G4, Grid<TG>& G5);
    void Resize(Grid<TG>& G0, Grid<TG>& G1,
                Grid<TG>& G2, Grid<TG>& G3,
                Grid<TG>& G4, Grid<TG>& G5,
                Grid<TG>& G6);
    void Resize(Grid<TG>& G0, Grid<TG>& G1,
                Grid<TG>& G2, Grid<TG>& G3,
                Grid<TG>& G4, Grid<TG>& G5,
                Grid<TG>& G6, Grid<TG>& G7);
    void Resize(Grid<TG>& G0, Grid<TG>& G1,
                Grid<TG>& G2, Grid<TG>& G3,
                Grid<TG>& G4, Grid<TG>& G5,
                Grid<TG>& G6, Grid<TG>& G7,
                Grid<TG>& G8);
    void Resize(Grid<TG>& G0, Grid<TG>& G1,
                Grid<TG>& G2, Grid<TG>& G3,
                Grid<TG>& G4, Grid<TG>& G5,
                Grid<TG>& G6, Grid<TG>& G7,
                Grid<TG>& G8, Grid<TG>& G9);

    void Resize(int N0);
    void Resize(int N0, int N1);
    void Resize(int N0, int N1,
                int N2);
    void Resize(int N0, int N1,
                int N2, int N3);
    void Resize(int N0, int N1,
                int N2, int N3,
                int N4);
    void Resize(int N0, int N1,
                int N2, int N3,
                int N4, int N5);
    void Resize(int N0, int N1,
                int N2, int N3,
                int N4, int N5,
                int N6);
    void Resize(int N0, int N1,
                int N2, int N3,
                int N4, int N5,
                int N6, int N7);
    void Resize(int N0, int N1,
                int N2, int N3,
                int N4, int N5,
                int N6, int N7,
                int N8);
    void Resize(int N0, int N1,
                int N2, int N3,
                int N4, int N5,
                int N6, int N7,
                int N8, int N9);

    void Resize(const TinyVector<int, N>& shape);

    // Calculus.
    void Mlt(T alpha);
    void Add(T alpha);
    void Apply(void function(T&));
    template <class F>
    void Apply(F& function);
    void Apply(T(function)(const T&));
    template <class T0, class TG0, class F>
    void Apply(Data<T0, N, TG0>&, F& function);

    T GetMax() const;
    T GetMaxAbs() const;
    T GetSignedMaxAbs() const;
    T GetMin() const;

    Array<int, 1> GetMaxIndex() const;
    Array<int, 1> GetMaxAbsIndex() const;
    Array<int, 1> GetMinIndex() const;

    T Sum() const;
    template <class Ts>
    void Sum(Ts& sum) const;
    T Mean() const;
    template <class Ts>
    void Mean(Ts& mean) const;
    T Variance() const;
    template <class Ts>
    void Variance(Ts& var) const;
    T StandardDeviation() const;
    template <class Ts>
    void StandardDeviation(Ts& std) const;

    T Norm1() const;
    T Norm2() const;
    T Norm(T p) const;

    void Fill();
    void Fill(T value);

    void SetZero();
    bool IsZero();
    void SetNaN();

    void Threshold(T threshold_min, T threshold_max);
    void ThresholdAbs(T threshold);
    void ThresholdMin(T threshold);
    void ThresholdMax(T threshold);

    template <class T0, class TG0>
    T NGE_interpolation(Data<T0, N, TG0>& data, T limit = T(0));
    template <class T0, class TG0>
    T NGE(Data<T0, N, TG0>& data, T limit = T(0));

    template <class T0, class TG0>
    T Bias_interpolation(Data<T0, N, TG0>& data);
    template <class T0, class TG0>
    T Bias(Data<T0, N, TG0>& data);

    template <class T0, class TG0>
    T RMS_interpolation(Data<T0, N, TG0>& data);
    template <class T0, class TG0>
    T RMS(Data<T0, N, TG0>& data);

    template <class T0, class TG0>
    T RelativeRMS_interpolation(Data<T0, N, TG0>& data);
    template <class T0, class TG0>
    T RelativeRMS(Data<T0, N, TG0>& data);

    template <class T0, class TG0>
    T Corr_interpolation(Data<T0, N, TG0>& data);
    template <class T0, class TG0>
    T Corr(Data<T0, N, TG0>& data);

    template <class T0, class TG0>
    T ErrorLessThan_interpolation(Data<T0, N, TG0>& data, T threshold);
    template <class T0, class TG0>
    T ErrorLessThan(Data<T0, N, TG0>& data, T threshold);

    void ReverseData(int dim = 0);
    void SwitchDimensions(TinyVector<int, N> NewDim, Grid<TG>& G0);
    void SwitchDimensions(TinyVector<int, N> NewDim, Grid<TG>& G0,
                          Grid<TG>& G1);
    void SwitchDimensions(TinyVector<int, N> NewDim, Grid<TG>& G0,
                          Grid<TG>& G1, Grid<TG>& G2);
    void SwitchDimensions(TinyVector<int, N> NewDim, Grid<TG>& G0,
                          Grid<TG>& G1, Grid<TG>& G2,
                          Grid<TG>& G3);
    void SwitchDimensions(TinyVector<int, N> NewDim, Grid<TG>& G0,
                          Grid<TG>& G1, Grid<TG>& G2,
                          Grid<TG>& G3, Grid<TG>& G4);
    void SwitchDimensions(TinyVector<int, N> NewDim, Grid<TG>& G0,
                          Grid<TG>& G1, Grid<TG>& G2,
                          Grid<TG>& G3, Grid<TG>& G4,
                          Grid<TG>& G5);
    void SwitchDimensions(TinyVector<int, N> NewDim, Grid<TG>& G0,
                          Grid<TG>& G1, Grid<TG>& G2,
                          Grid<TG>& G3, Grid<TG>& G4,
                          Grid<TG>& G5, Grid<TG>& G6);
    void SwitchDimensions(TinyVector<int, N> NewDim, Grid<TG>& G0,
                          Grid<TG>& G1, Grid<TG>& G2,
                          Grid<TG>& G3, Grid<TG>& G4,
                          Grid<TG>& G5, Grid<TG>& G6,
                          Grid<TG>& G7);
    void SwitchDimensions(TinyVector<int, N> NewDim, Grid<TG>& G0,
                          Grid<TG>& G1, Grid<TG>& G2,
                          Grid<TG>& G3, Grid<TG>& G4,
                          Grid<TG>& G5, Grid<TG>& G6,
                          Grid<TG>& G7, Grid<TG>& G8);
    void SwitchDimensions(TinyVector<int, N> NewDim, Grid<TG>& G0,
                          Grid<TG>& G1, Grid<TG>& G2,
                          Grid<TG>& G3, Grid<TG>& G4,
                          Grid<TG>& G5, Grid<TG>& G6,
                          Grid<TG>& G7, Grid<TG>& G8,
                          Grid<TG>& G9);

    void ChangeCoords(FuncCoords_Base<TG>& f);
    void ChangeCoordsInPlace(Function_Base<TG>& f);

    void Print() const;
    void PrintInfo() const;
    string InfoString() const;

  private:
    void ClearGrids();
    void SetVariables();
    void InitData();
  };


}  // namespace Data.


#define FILE_SELDONDATA_DATA_HXX
#endif
