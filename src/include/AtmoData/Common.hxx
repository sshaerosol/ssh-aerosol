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


#ifndef ATMODATA_FILE_COMMON_HXX

namespace AtmoData
{

#define SWAP_4(a) ((a >> 24) & 0xFF) | ((a >> 8) & 0xFF00) |	\
  ((a << 8) & 0x00FF0000) | ((a << 24) & 0xFF000000)

  inline float swap(float& x)
  {
    return (*(unsigned *)&x = SWAP_4(*(unsigned *)&x));
  }

  inline int swap(int& x)
  {
    return (*(unsigned *)&x = SWAP_4(*(unsigned *)&x));
  }

  inline unsigned int swap(unsigned int& x)
  {
    return (*(unsigned *)&x = SWAP_4(*(unsigned *)&x));
  }

  inline unsigned long swap(unsigned long& x)
  {
    return (*(unsigned *)&x = SWAP_4(*(unsigned *)&x));
  }

  template <class T, int N>
  inline void swap(Array<T, N>& A)
  {
    for (int i = 0; i < (int) A.size(); i++)
      swap(A.data()[i]);
  }

}  // namespace AtmoData.


#define ATMODATA_FILE_COMMON_HXX
#endif
