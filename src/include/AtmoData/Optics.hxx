// Copyright (C) 2003-2007, ENPC - INRIA - EDF R&D
// Author(s): Vivien Mallet, Elsa Real
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



#ifndef ATMODATA_FILE_OPTIC_HXX

namespace AtmoData
{

  template<class T>
  T sqr(T x);

  template<class T>
  T sqrt2(T a);

  template<class T>
  T arccos2(T a);

  template<class T>
  T sign2(T a);

  template<class T>
  void sqrt_for_complex(T a, T b, T& a_out, T& b_out);

  template<class T>
  void sqr_for_complex(T a, T b, T& a_out, T& b_out);

  template<class T>
  void quotient_for_complex(T divisor_real, T divisor_imaginary,
                            T dividend_real, T dividend_imaginary,
                            T& quotient_real, T& quotient_imaginary);

  template<class T>
  void compute_refractive_index
  (Data<T, 1> species_concentration,
   Data<T, 1> species_refractive_index_real,
   Data<T, 1> species_refractive_index_imaginary,
   T& global_refractive_index_real,
   T& global_refractive_index_imaginary);

  template<class T>
  void compute_refractive_index_Lorentz_Lorenz
  (Data<T, 1> species_concentration,
   Data<T, 1> species_refractive_index_real,
   Data<T, 1> species_refractive_index_imaginary,
   T& global_refractive_index_real,
   T& global_refractive_index_imaginary);

  template<class T>
  void compute_refractive_index_Maxwell_Garnet
  (T inclusion_concentration, T solution_concentration,
   T inclusion_refractive_index_real,
   T inclusion_refractive_index_imaginary,
   T solution_refractive_index_real,
   T solution_refractive_index_imaginary,
   T& global_refractive_index_real, T& global_refractive_index_imaginary);

  template<class T>
  void compute_Hanel_index(T dry_refractive_index_real,
                           T dry_refractive_index_imaginary,
                           T water_refractive_index_real,
                           T water_refractive_index_imaginary,
                           T dry_radius, T wet_radius,
                           T& wet_refractive_index_real,
                           T& wet_refractive_index_imaginary);



}


#define ATMODATA_FILE_OPTIC_HXX
#endif
