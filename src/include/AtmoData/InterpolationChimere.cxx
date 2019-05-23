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

#ifndef ATMODATA_FILE_INTERPOLATIONCHIMERE_CXX

// This code is essentially based on the chemistry-transport model Chimere,
// distributed under GNU GPL -- copyright (C) 2005 Institut Pierre-Simon
// Laplace (CNRS), INERIS, LISA (CNRS).

namespace AtmoData
{

  template<class T>
    void
    InterpolationRegularInput (T x_min_in, T Delta_x_in,
                               const RegularGrid<T>& GridX_out,
                               Array<int, 1>& Index, Array<T, 1>& Weight)
    {
      int Nx_out = GridX_out.GetNbElements ();

      for (int i = 0; i < Nx_out; i++)
        {
          Index (i) = int ((GridX_out (i) - x_min_in) / Delta_x_in);
          Weight (i) = (GridX_out (i) - x_min_in - T (Index (i)) * Delta_x_in)
              / Delta_x_in;
        }
    }

  template<class T>
    void
    InterpolationRegularInput (const RegularGrid<T>& GridX_in,
                               const RegularGrid<T>& GridX_out,
                               Array<int, 1>& Index, Array<T, 1>& Weight)
    {
      T x_min_in = GridX_in.Value (0);
      T Delta_x_in = GridX_in.Value (1) - GridX_in.Value (0);
      InterpolationRegularInput (x_min_in, Delta_x_in, GridX_out, Index,
                                 Weight);
    }

  template<class T>
    void
    InterpolationRegularInput (T x_min_in, T Delta_x_in, T y_min_in,
                               T Delta_y_in, const Array<T, 2>& GridX_2D_out,
                               const Array<T, 2>& GridY_2D_out,
                               Array<int, 2>& Index_i_in,
                               Array<int, 2>& Index_j_in,
                               Array<T, 2>& Weight_i_in,
                               Array<T, 2>& Weight_j_in)
    {
      int j, i;
      int Ny_out = GridX_2D_out.extent (0);
      int Nx_out = GridX_2D_out.extent (1);

      for (j = 0; j < Ny_out; j++)
        for (i = 0; i < Nx_out; i++)
          {
            Index_i_in (j, i) = int (
                (GridX_2D_out (j, i) - x_min_in) / Delta_x_in);
            Index_j_in (j, i) = int (
                (GridY_2D_out (j, i) - y_min_in) / Delta_y_in);
            Weight_i_in (j, i) = (GridX_2D_out (j, i) - x_min_in
                - T (Index_i_in (j, i)) * Delta_x_in) / Delta_x_in;
            Weight_j_in (j, i) = (GridY_2D_out (j, i) - y_min_in
                - T (Index_j_in (j, i)) * Delta_y_in) / Delta_y_in;
          }
    }

  template<class T>
    void
    InterpolationChimere (const Array<T, 2>& Longitude_in,
                          const Array<T, 2>& Latitude_in, T x_min_out,
                          T Delta_x_out, int Nx_out, T y_min_out, T Delta_y_out,
                          int Ny_out, Array<int, 2>& Index_i_in,
                          Array<int, 2>& Index_j_in, Array<T, 2>& Weight_i_in,
                          Array<T, 2>& Weight_j_in)
    {
      int j, i;
      int j_in, i_in;
      int Nsec;
      T x, y, xx, yy;
      T d2, d2min, a (0.), a0, a1, b (0.), b0, b1;
      int ntry (100);

      int Nx_in = Longitude_in.extent (0);
      int Ny_in = Longitude_in.extent (1);

      for (j = 0; j < Ny_out; j++)
        for (i = 0; i < Nx_out; i++)
          {
            x = x_min_out + T (i) * Delta_x_out;
            y = y_min_out + T (j) * Delta_y_out;

            i_in = 0;
            j_in = 0;
            Nsec = 0;

            while (i_in != Nx_in - 1 && Nsec != 1)
              {
                Nsec = 0;

                if (Latitude_in (i_in + 1, j_in + 0)
                    != Latitude_in (i_in + 0, j_in + 0)
                    && y
                        <= max (Latitude_in (i_in + 1, j_in + 0),
                                Latitude_in (i_in + 0, j_in + 0))
                    && y
                        >= min (Latitude_in (i_in + 1, j_in + 0),
                                Latitude_in (i_in + 0, j_in + 0)))
                  {
                    a = (Latitude_in (i_in + 1, j_in + 0) - y)
                        / (Latitude_in (i_in + 1, j_in + 0)
                            - Latitude_in (i_in + 0, j_in + 0));
                    if (a >= 0 && a < 1.)
                      {
                        b = a * (Longitude_in (i_in + 0, j_in + 0) - x)
                            + (1. - a)
                                * (Longitude_in (i_in + 1, j_in + 0) - x);
                        if (b > 0.)
                          Nsec++;
                      }
                  }

                if (Latitude_in (i_in + 1, j_in + 1)
                    != Latitude_in (i_in + 1, j_in + 0)
                    && y
                        <= max (Latitude_in (i_in + 1, j_in + 1),
                                Latitude_in (i_in + 1, j_in + 0))
                    && y
                        >= min (Latitude_in (i_in + 1, j_in + 1),
                                Latitude_in (i_in + 1, j_in + 0)))
                  {
                    a = (Latitude_in (i_in + 1, j_in + 1) - y)
                        / (Latitude_in (i_in + 1, j_in + 1)
                            - Latitude_in (i_in + 1, j_in + 0));
                    if (a >= 0. && a < 1.)
                      {
                        b = a * (Longitude_in (i_in + 1, j_in + 0) - x)
                            + (1. - a)
                                * (Longitude_in (i_in + 1, j_in + 1) - x);
                        if (b > 0.)
                          Nsec++;
                      }
                  }

                if (Latitude_in (i_in + 0, j_in + 1)
                    != Latitude_in (i_in + 1, j_in + 1)
                    && y
                        <= max (Latitude_in (i_in + 0, j_in + 1),
                                Latitude_in (i_in + 1, j_in + 1))
                    && y
                        >= min (Latitude_in (i_in + 0, j_in + 1),
                                Latitude_in (i_in + 1, j_in + 1)))
                  {
                    a = (Latitude_in (i_in + 0, j_in + 1) - y)
                        / (Latitude_in (i_in + 0, j_in + 1)
                            - Latitude_in (i_in + 1, j_in + 1));
                    if (a >= 0. && a < 1.)
                      {
                        b = a * (Longitude_in (i_in + 1, j_in + 1) - x)
                            + (1. - a)
                                * (Longitude_in (i_in + 0, j_in + 1) - x);
                        if (b > 0.)
                          Nsec++;
                      }
                  }

                if (Latitude_in (i_in + 0, j_in + 0)
                    != Latitude_in (i_in + 0, j_in + 1)
                    && y
                        <= max (Latitude_in (i_in + 0, j_in + 0),
                                Latitude_in (i_in + 0, j_in + 1))
                    && y
                        >= min (Latitude_in (i_in + 0, j_in + 0),
                                Latitude_in (i_in + 0, j_in + 1)))
                  {
                    a = (Latitude_in (i_in + 0, j_in + 0) - y)
                        / (Latitude_in (i_in + 0, j_in + 0)
                            - Latitude_in (i_in + 0, j_in + 1));
                    if (a >= 0. && a < 1.)
                      {
                        b = a * (Longitude_in (i_in + 0, j_in + 1) - x)
                            + (1. - a)
                                * (Longitude_in (i_in + 0, j_in + 0) - x);
                        if (b > 0.)
                          Nsec++;
                      }
                  }

                if (x == Longitude_in (i_in + 0, j_in + 0)
                    && y == Latitude_in (i_in + 0, j_in + 0))
                  Nsec = 1;
                if (x == Longitude_in (i_in + 1, j_in + 0)
                    && y == Latitude_in (i_in + 1, j_in + 0))
                  Nsec = 1;
                if (x == Longitude_in (i_in + 1, j_in + 1)
                    && y == Latitude_in (i_in + 1, j_in + 1))
                  Nsec = 1;
                if (x == Longitude_in (i_in + 0, j_in + 1)
                    && y == Latitude_in (i_in + 0, j_in + 1))
                  Nsec = 1;

                if (Nsec == 1)
                  {
                    Index_i_in (j, i) = i_in;
                    Index_j_in (j, i) = j_in;
                  }

                j_in++;
                if (j_in == Ny_in - 2)
                  {
                    j_in = 0;
                    i_in++;
                  }
              }

            d2min = 1.e20;
            i_in = Index_i_in (j, i);
            j_in = Index_j_in (j, i);

            for (int ntx = 0; ntx < ntry; ntx++)
              {
                a0 = T (ntx) / T (ntry);
                a1 = 1. - a0;
                for (int nty = 0; nty < ntry; nty++)
                  {
                    b0 = T (nty) / T (ntry);
                    b1 = 1. - b0;
                    xx = a0 * b0 * Longitude_in (i_in + 1, j_in + 1)
                        + a1 * b0 * Longitude_in (i_in, j_in + 1)
                        + a0 * b1 * Longitude_in (i_in + 1, j_in)
                        + a1 * b1 * Longitude_in (i_in, j_in);
                    yy = a0 * b0 * Latitude_in (i_in + 1, j_in + 1)
                        + a1 * b0 * Latitude_in (i_in, j_in + 1)
                        + a0 * b1 * Latitude_in (i_in + 1, j_in)
                        + a1 * b1 * Latitude_in (i_in, j_in);
                    d2 = (x - xx) * (x - xx) + (y - yy) * (y - yy);
                    if (d2 < d2min)
                      {
                        d2min = d2;
                        a = a0;
                        b = b0;
                      }
                  }
              }
            Weight_i_in (j, i) = a;
            Weight_j_in (j, i) = b;
          }
    }

  template<class T>
    void
    HorizontalInterpolation (const Array<T, 2>& Field_in,
                             const Array<int, 1>& Index_x,
                             const Array<int, 1>& Index_y,
                             const Array<T, 1>& Weight_x,
                             const Array<T, 1>& Weight_y,
                             Array<T, 2>& Field_out)
    {
      int Ny = Field_out.extent (0);
      int Nx = Field_out.extent (1);

      int i, j, i_in, j_in;
      T weight_x, weight_y, one_weight_x, one_weight_y;

      for (j = 0; j < Ny; j++)
        for (i = 0; i < Nx; i++)
          {
            i_in = Index_x (i);
            j_in = Index_y (j);

            weight_x = Weight_x (i);
            weight_y = Weight_y (j);
            one_weight_x = 1. - weight_x;
            one_weight_y = 1. - weight_y;

            Field_out (j, i) = one_weight_x * one_weight_y
                * Field_in (j_in, i_in)
                + weight_x * one_weight_y * Field_in (j_in, i_in + 1)
                + one_weight_x * weight_y * Field_in (j_in + 1, i_in)
                + weight_x * weight_y * Field_in (j_in + 1, i_in + 1);
          }
    }

  template<class T>
    void
    HorizontalInterpolation (const Array<T, 3>& Field_in,
                             const Array<int, 1>& Index_x,
                             const Array<int, 1>& Index_y,
                             const Array<T, 1>& Weight_x,
                             const Array<T, 1>& Weight_y,
                             Array<T, 3>& Field_out)
    {
      int Nz = Field_out.extent (0);
      int Ny = Field_out.extent (1);
      int Nx = Field_out.extent (2);

      if (Nz != Field_in.extent (0))
        throw "Fields should be of the same size along z.";

      int i, j, k, i_in, j_in;
      T weight_x, weight_y, one_weight_x, one_weight_y;

      for (k = 0; k < Nz; k++)
        for (j = 0; j < Ny; j++)
          for (i = 0; i < Nx; i++)
            {
              i_in = Index_x (i);
              j_in = Index_y (j);

              weight_x = Weight_x (i);
              weight_y = Weight_y (j);
              one_weight_x = 1. - weight_x;
              one_weight_y = 1. - weight_y;

              Field_out (k, j, i) = one_weight_x * one_weight_y
                  * Field_in (k, j_in, i_in)
                  + weight_x * one_weight_y * Field_in (k, j_in, i_in + 1)
                  + one_weight_x * weight_y * Field_in (k, j_in + 1, i_in)
                  + weight_x * weight_y * Field_in (k, j_in + 1, i_in + 1);
            }
    }

  template<class T>
    void
    HorizontalInterpolation (const Array<T, 3>& Field_in,
                             const Array<int, 2>& Index_x,
                             const Array<int, 2>& Index_y,
                             const Array<T, 2>& Weight_x,
                             const Array<T, 2>& Weight_y,
                             Array<T, 3>& Field_out)
    {
      int Nz = Field_out.extent (0);
      int Ny = Field_out.extent (1);
      int Nx = Field_out.extent (2);

      if (Nz != Field_in.extent (0))
        throw "Fields should be of the same size along z.";

      int i, j, k, i_in, j_in;
      T weight_x, weight_y, one_weight_x, one_weight_y;

      for (k = 0; k < Nz; k++)
        for (j = 0; j < Ny; j++)
          for (i = 0; i < Nx; i++)
            {
              i_in = Index_x (j, i);
              j_in = Index_y (j, i);

              weight_x = Weight_x (j, i);
              weight_y = Weight_y (j, i);
              one_weight_x = 1. - weight_x;
              one_weight_y = 1. - weight_y;

              Field_out (k, j, i) = one_weight_x * one_weight_y
                  * Field_in (k, j_in, i_in)
                  + weight_x * one_weight_y * Field_in (k, j_in, i_in + 1)
                  + one_weight_x * weight_y * Field_in (k, j_in + 1, i_in)
                  + weight_x * weight_y * Field_in (k, j_in + 1, i_in + 1);
            }
    }

  template<class T>
    void
    HorizontalInterpolation (Array<T, 4>& Field_in,
                             const Array<int, 2>& Index_x,
                             const Array<int, 2>& Index_y,
                             const Array<T, 2>& Weight_x,
                             const Array<T, 2>& Weight_y,
                             Array<T, 4>& Field_out)
    {
      int Nt = Field_out.extent (0);

      if (Nt != Field_in.extent (0))
        throw "Fields should be of the same size along t.";

      TinyVector<int, 3> shape_in (Field_in.extent (1), Field_in.extent (2),
                                   Field_in.extent (3));
      TinyVector<int, 3> shape_out (Field_out.extent (1), Field_out.extent (2),
                                    Field_out.extent (3));

      for (int h = 0; h < Nt; h++)
        {
          Array<T, 3> Field_in_extract (&Field_in (h, 0, 0, 0), shape_in);
          Array<T, 3> Field_out_extract (&Field_out (h, 0, 0, 0), shape_out);
          HorizontalInterpolation (Field_in_extract, Index_x, Index_y, Weight_x,
                                   Weight_y, Field_out_extract);
        }
    }

  template<class T>
    void
    HorizontalInterpolation_ground (Array<T, 4>& Field_in,
                                    const Array<int, 2>& Index_x,
                                    const Array<int, 2>& Index_y,
                                    const Array<T, 2>& Weight_x,
                                    const Array<T, 2>& Weight_y,
                                    Array<T, 3>& GridZ_out,
                                    Array<T, 4>& Field_out, T ground_threshold)
    {
      int h, j, i;

      int Nt = Field_out.extent (0);
      int Nz = Field_out.extent (1);
      int Ny = Field_out.extent (2);
      int Nx = Field_out.extent (3);

      if (Nz != Field_in.extent (1) + 1)
        throw "Output field should have an extra level for ground.";

      // Horizontal interpolation at all levels except ground.
      // First excludes the ground.
      Array<T, 4> Field_out_extract (Field_out, Range::all (),
                                     Range (1, Nz - 1));
      HorizontalInterpolation (Field_in, Index_x, Index_y, Weight_x, Weight_y,
                               Field_out_extract);

      // Vertical interpolation at ground.
      for (h = 0; h < Nt; h++)
        for (j = 0; j < Ny; j++)
          for (i = 0; i < Nx; i++)
            {
              Field_out (h, 0, j, i) = Field_out (h, 1, j, i)
                  + (Field_out (h, 1, j, i) - Field_out (h, 2, j, i))
                      * GridZ_out (1, j, i)
                      / (GridZ_out (2, j, i) - GridZ_out (1, j, i));
              Field_out (h, 0, j, i) = max (ground_threshold,
                                            Field_out (h, 0, j, i));
            }
    }

  template<class T>
    void
    VerticalAverage (const Array<T, 3>& GridZ_in, const Array<T, 4>& Field_in,
                     const Array<T, 4>& GridZ_out, Array<T, 4>& Field_out)
    {
      int h, k, j, i, k_in;

      int Nt_out = Field_out.extent (0);
      int Nz_out = Field_out.extent (1);
      int Ny_out = Field_out.extent (2);
      int Nx_out = Field_out.extent (3);

      T cumulated, additional, distance, height;
      for (h = 0; h < Nt_out; h++)
        for (j = 0; j < Ny_out; j++)
          for (i = 0; i < Nx_out; i++)
            {
              cumulated = 0.;
              k_in = 1;
              for (k = 0; k < Nz_out; k++)
                {
                  while (GridZ_in (k_in, j, i) < GridZ_out (h, k + 1, j, i))
                    {
                      cumulated += .5
                          * (GridZ_in (k_in, j, i) - GridZ_in (k_in - 1, j, i))
                          * (Field_in (h, k_in, j, i)
                              + Field_in (h, k_in - 1, j, i));
                      k_in++;
                    }
                  distance = GridZ_out (h, k + 1, j, i)
                      - GridZ_in (k_in - 1, j, i);
                  height = GridZ_in (k_in, j, i) - GridZ_in (k_in - 1, j, i);
                  additional = Field_in (h, k_in - 1, j, i)
                      + .5 * distance / height
                          * (Field_in (h, k_in, j, i)
                              - Field_in (h, k_in - 1, j, i));
                  additional *= distance;
                  Field_out (h, k, j, i) = cumulated + additional;
                  Field_out (h, k, j, i) /= GridZ_out (h, k + 1, j, i)
                      - GridZ_out (h, k, j, i);
                  cumulated = -additional;
                }
            }
    }

}  // namespace AtmoData.

#define ATMODATA_FILE_INTERPOLATIONCHIMERE_CXX
#endif

