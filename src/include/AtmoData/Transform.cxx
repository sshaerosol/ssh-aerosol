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


#ifndef ATMODATA_FILE_TRANSFORM_CXX

#include "Transform.hxx"

namespace AtmoData
{


  //! Decumulates data.
  /*!
    'data' stores values cumulated in time, where the time is the first
    dimension.
    \param data data to be decumulated.
    \param length number of time steps over which data is cumulated.
    \param start (optional) index of the first complete cycle. Default: 0.
  */
  template <class T, int N, class TG>
  void Decumulate(Data<T, N, TG>& data, int length,
                  int start)
  {
    unsigned int n = data.GetNbElements();
    unsigned int l = n / data.GetLength(0);

    for (unsigned int i = data.GetLength(0) - 1; i > 0; i--)
      if ((i - start) % length != 0)
        {
          T* data_arr = &data.GetData()[l * i];
          T* data_prev = &data.GetData()[l * (i - 1)];
          for (unsigned int j = 0; j < l; j++)
            data_arr[j] -= data_prev[j];
        }
  }


}  // namespace AtmoData.

#define ATMODATA_FILE_TRANSFORM_CXX
#endif
