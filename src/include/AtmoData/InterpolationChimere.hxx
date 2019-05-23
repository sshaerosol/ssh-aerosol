// Copyright (C) 2006-2007, ENPC - INRIA - EDF R&D
// Author(s): Vivien Mallet
//
// This file is part of the air quality modeling system Polyphemus.
//
// Polyphemus is developed in the INRIA - ENPC joint project-team CLIME and in
// the ENPC - EDF R&D joint laboratory CEREA.
//
// Polyphemus is free software; you can redistribute it and/or modify it under
// the terms of the GNU General Public License as published by the Free
// Software Foundation; either version 2 of the License, or (at your option)
// any later version.
//
// Polyphemus is distributed in the hope that it will be useful, but WITHOUT
// ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
// FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for
// more details.
//
// For more information, visit the Polyphemus web site:
//      http://cerea.enpc.fr/polyphemus/

#ifndef ATMODATA_FILE_INTERPOLATIONCHIMERE_HXX

// This code is essentially based on the chemistry-transport model Chimere,
// distributed under GNU GPL -- copyright (C) 2005 Institut Pierre-Simon
// Laplace (CNRS), INERIS, LISA (CNRS).

namespace AtmoData
{

  template<class T>
    void
    InterpolationRegularInput (T x_min_in, T Delta_x_in,
                               const RegularGrid<T>& GridX_out,
                               Array<int, 1>& Index, Array<T, 1>& Weight);

  template<class T>
    void
    InterpolationRegularInput (const RegularGrid<T>& GridX_in,
                               const RegularGrid<T>& GridX_out,
                               Array<int, 1>& Index, Array<T, 1>& Weight);

  template<class T>
    void
    InterpolationRegularInput (T x_min_in, T Delta_x_in, T y_min_in,
                               T Delta_y_in, const Array<T, 2>& GridX_2D_out,
                               const Array<T, 2>& GridY_2D_out,
                               Array<int, 2>& Index_i_in,
                               Array<int, 2>& Index_j_in,
                               Array<T, 2>& Weight_i_in,
                               Array<T, 2>& Weight_j_in);

  template<class T>
    void
    InterpolationChimere (const Array<T, 2>& Longitude_in,
                          const Array<T, 2>& Latitude_in, T x_min_out,
                          T Delta_x_out, int Nx_out, T y_min_out, T Delta_y_out,
                          int Ny_out, Array<int, 2>& Index_i_in,
                          Array<int, 2>& Index_j_in, Array<T, 2>& Weight_i_in,
                          Array<T, 2>& Weight_j_in);

  template<class T>
    void
    HorizontalInterpolation (const Array<T, 2>& Field_in,
                             const Array<int, 1>& Index_x,
                             const Array<int, 1>& Index_y,
                             const Array<T, 1>& Weight_x,
                             const Array<T, 1>& Weight_y,
                             Array<T, 2>& Field_out);

  template<class T>
    void
    HorizontalInterpolation (const Array<T, 3>& Field_in,
                             const Array<int, 1>& Index_x,
                             const Array<int, 1>& Index_y,
                             const Array<T, 1>& Weight_x,
                             const Array<T, 1>& Weight_y,
                             Array<T, 3>& Field_out);

  template<class T>
    void
    HorizontalInterpolation (const Array<T, 3>& Field_in,
                             const Array<int, 2>& Index_x,
                             const Array<int, 2>& Index_y,
                             const Array<T, 2>& Weight_x,
                             const Array<T, 2>& Weight_y,
                             Array<T, 3>& Field_out);

  template<class T>
    void
    HorizontalInterpolation (Array<T, 4>& Field_in,
                             const Array<int, 2>& Index_x,
                             const Array<int, 2>& Index_y,
                             const Array<T, 2>& Weight_x,
                             const Array<T, 2>& Weight_y,
                             Array<T, 4>& Field_out);

  template<class T>
    void
    HorizontalInterpolation_ground (Array<T, 4>& Field_in,
                                    const Array<int, 2>& Index_x,
                                    const Array<int, 2>& Index_y,
                                    const Array<T, 2>& Weight_x,
                                    const Array<T, 2>& Weight_y,
                                    Array<T, 3>& GridZ_out,
                                    Array<T, 4>& Field_out, T ground_threshold =
                                        0.);

  template<class T>
    void
    VerticalAverage (const Array<T, 3>& GridZ_in, const Array<T, 4>& Field_in,
                     const Array<T, 4>& GridZ_out, Array<T, 4>& Field_out);

}  // namespace AtmoData.

#define ATMODATA_FILE_INTERPOLATIONCHIMERE_HXX
#endif
