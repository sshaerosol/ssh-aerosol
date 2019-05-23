// Copyright (C) 2007, ENPC - INRIA - EDF R&D
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

#ifndef ATMODATA_FILE_OPTICS_CXX

#include "Optics.hxx"

#include <cmath>

namespace AtmoData
{
  //
  // Math fonctions.
  //

  template<class T>
  T sqr(T x)
  {
    return x * x;
  }

  template<class T>
  T sqrt2(T x)
  {
    return sqrt(double(x));
  }

  template<class T>
  T arccos2(T x)
  {
    return acos(double(x));
  }

  template<class T>
  T sign2(T x)
  {
    if (x != T(0.))
      return x / T(abs(double(x)));
    else
      return T(1.);
  }

  // Square root of a complex number.
  template<class T>
  void sqrt_for_complex(T re, T im, T& re_out, T& im_out)
  {
    T rho_in = sqrt2(sqr(re) + sqr(im));
    T phi_in = sign2(im) * arccos2(T(re / rho_in));
    T rho_out = sqrt2(rho_in);
    T phi_out = phi_in / 2.;
    re_out = rho_out * cos(phi_out);
    im_out = rho_out * sin(phi_out);
  }


  // Square of a complex number.
  template<class T>
  void sqr_for_complex(T re, T im, T& re_out, T& im_out)
  {
    re_out = sqr(re) - sqr(im);
    im_out = 2. * re * im;
  }


  // Quotient of two complex numbers.
  template<class T>
  void quotient_for_complex(T divisor_real, T divisor_imaginary,
                            T dividend_real, T dividend_imaginary,
                            T& quotient_real, T& quotient_imaginary)
  {
    quotient_real = (divisor_real * dividend_real
                     + divisor_imaginary * dividend_imaginary)
      / (sqr(dividend_real) + sqr(dividend_imaginary));

    quotient_imaginary = (divisor_imaginary * dividend_real
                          - divisor_real * dividend_imaginary)
      / (sqr(dividend_real) + sqr(dividend_imaginary));
  }

  // Computes the aerosol refractive index from the individual specie ones
  // in case of homogeneous mixing.
  /*!
    \param species_concentration concentrations of individual species.
    \param species_refractive_index_real real part of the refractive
    \        index of individual species.
    \param species_refractive_index_imaginary imaginary part of the
    \        refractive index of individual species.
    \param[out] global_refractive_index_real aerosol refractive index real part.
    \param[out] global_refractive_index_imaginary aerosol refractive index
    \             imaginary part.
  */
  template<class T>
  void compute_refractive_index
  (Data<T, 1> species_concentration,
   Data<T, 1> species_refractive_index_real,
   Data<T, 1> species_refractive_index_imaginary,
   T& global_refractive_index_real,
   T& global_refractive_index_imaginary)
  {
    int Nspecies = species_concentration[0].GetLength();

    T global_concentration(0);
    for (int i = 0; i < Nspecies; i++)
      global_concentration += species_concentration(i);

    global_refractive_index_real = 0;
    global_refractive_index_imaginary = 0;
    for (int i = 0; i < Nspecies; i++)
      {
        T frac = species_concentration(i) / global_concentration;
        global_refractive_index_real +=
          frac * species_refractive_index_real(i);
        global_refractive_index_imaginary +=
          frac * species_refractive_index_imaginary(i);
      }
  }

  // Computes the aerosol refractive index from the individual specie ones
  // (in case of homogeneous mixing) following the Lorentz_Lorenz
  // electromagnetism formula.
  /*!
    \param species_concentration concentrations of individual species.
    \param species_refractive_index_real real part of the refractive
    \        index of individual species.
    \param species_refractive_index_imaginary imaginary part of the
    \        refractive index of individual species.
    \param[out] global_refractive_index_real aerosol refractive index real part.
    \param[out] global_refractive_index_imaginary aerosol refractive index
    \             imaginary part.
  */
  template<class T>
  void compute_refractive_index_Lorentz_Lorenz
  (Data<T, 1> species_concentration,
   Data<T, 1> species_refractive_index_real,
   Data<T, 1> species_refractive_index_imaginary,
   T& global_refractive_index_real,
   T& global_refractive_index_imaginary)
  {
    int Nspecies = species_concentration[0].GetLength();
    T quotient_real, quotient_imaginary;
    global_refractive_index_real = 0.;
    global_refractive_index_imaginary = 0.;

    T global_concentration(0.);
    for (int i = 0; i < Nspecies; i++)
      global_concentration += species_concentration(i);

    T total_real(0.), total_imaginary(0.);
    for (int i = 0; i < Nspecies; i++)
      {
        T frac = species_concentration(i) / global_concentration;
        T square_real, square_imaginary;
        sqr_for_complex(species_refractive_index_real(i),
                        species_refractive_index_imaginary(i),
                        square_real, square_imaginary);
        quotient_for_complex(square_real - 1, square_imaginary,
                             square_real + 2, square_imaginary,
                             quotient_real, quotient_imaginary);
        total_real += quotient_real * frac;
        total_imaginary += quotient_imaginary * frac;
      }

    quotient_for_complex(2 * total_real + 1, 2 * total_imaginary,
                         1 - total_real, -total_imaginary,
                         quotient_real, quotient_imaginary);
    sqrt_for_complex(quotient_real, quotient_imaginary,
                     global_refractive_index_real,
                     global_refractive_index_imaginary);
  }

  //! Computes the aerosol refractive index from the individula ones
  //  in the case of a "core" mixing, using the Maxwell Garnet formula
  //  [Lesins et al., JGR 2002].
  /*!
    \param inclusion_concentration concentration of the core (black carbon).
    \param solution_concentration concentration of the solution
    \         (all species, including water, except black carbon).
    \param inclusion_refractive_index_real real part of the core
    \                                          refractive index.
    \param inclusion_refractive_index_imaginary imaginary part of the
    \                                         core refractive index.
    \param solution_refractive_index_real real part of the solution
    \                                              refractive index.
    \param solution_refractive_index_imaginary imaginary part of the
    \                                     solution refractive index.
    \param[out] global_refractive_index_real real part of the aerosol
    \                                           refractive index.
    \param[out] global_refractive_index_imaginary imaginary part of the
    \                                    aerosol refractive index.
  */
  template<class T>
  void compute_refractive_index_Maxwell_Garnet
  (T inclusion_concentration, T solution_concentration,
   T inclusion_refractive_index_real,
   T inclusion_refractive_index_imaginary,
   T solution_refractive_index_real,
   T solution_refractive_index_imaginary,
   T& global_refractive_index_real, T& global_refractive_index_imaginary)
  {
    T inclusion_ratio = inclusion_concentration /
      (inclusion_concentration + solution_concentration);

    // Dielectric constant is the square of the refractive index.
    T inclusion_dielectric_real = sqr(inclusion_refractive_index_real)
      - sqr(inclusion_refractive_index_imaginary);

    T inclusion_dielectric_imaginary = 2. * inclusion_refractive_index_real
      * inclusion_refractive_index_imaginary;

    T solution_dielectric_real = sqr(solution_refractive_index_real)
      - sqr(solution_refractive_index_imaginary);

    T solution_dielectric_imaginary = 2. * solution_refractive_index_real
      * solution_refractive_index_imaginary;

    // Computes dividend and divisor of the quotient.
    T dividend_real = inclusion_dielectric_real
      + 2. * solution_dielectric_real + 2. * inclusion_ratio
      * (inclusion_dielectric_real - solution_dielectric_real);

    T dividend_imaginary = inclusion_dielectric_imaginary
      + 2. * solution_dielectric_imaginary + 2. * inclusion_ratio
      * (inclusion_dielectric_imaginary - solution_dielectric_imaginary);

    T divisor_real = inclusion_dielectric_real
      + 2. * solution_dielectric_real
      - inclusion_ratio * (inclusion_dielectric_real
                           - solution_dielectric_real);

    T divisor_imaginary = inclusion_dielectric_imaginary
      + 2 * solution_dielectric_imaginary
      - inclusion_ratio
      * (inclusion_dielectric_imaginary - solution_dielectric_imaginary);

    T numerator_real = solution_dielectric_real * dividend_real
      - solution_dielectric_imaginary * dividend_imaginary;
    T numerator_imaginary = solution_dielectric_imaginary * dividend_real
      + solution_dielectric_real * dividend_imaginary;

    // Complex number with exponential.
    T rho_numerator = sqrt2(sqr(numerator_real) + sqr(numerator_imaginary));

    T theta_numerator = sign2(numerator_imaginary)
      * arccos2(numerator_real / rho_numerator);

    T rho_denominator = sqrt2(sqr(divisor_real) + sqr(divisor_imaginary));
    T theta_denominator = sign2(divisor_imaginary)
      * arccos2(divisor_real / rho_denominator);

    T rho = rho_numerator / rho_denominator;
    T theta = theta_numerator - theta_denominator;

    global_refractive_index_real = sqrt(rho) * cos(theta / 2.);
    global_refractive_index_imaginary = sqrt(rho) * sin(theta / 2.);
  }


  //! Computes the aerosol wet refractive index from the dry aerosol refractive
  // index and water refractive index, using the Hanel formula.
  /*!
    \param dry_refractive_index_real real part of the dry refractive index.
    \param dry_refractive_index_imaginary imaginary part of the dry
    \                                               refractive index.
    \param water_refractive_index_real real part of the water
    \                                        refractive index.
    \param water_refractive_index_imaginary imaginary part of
    \                             the water refractive index.
    \param dry_radius aerosol dry radius (in microg/m^3).
    \param wet_radius aerosol wet radius (in microg/m^3).
    \param[out] wet_refractive_index_real real part of the aerosol
    \                                    wet refractive index.
    \param[out] wet_refractive_index_imaginary imaginary part of the
    \                             aerosol wet refractive index.
  */
  template<class T>
  void compute_Hanel_index(T dry_refractive_index_real,
                           T dry_refractive_index_imaginary,
                           T water_refractive_index_real,
                           T water_refractive_index_imaginary,
                           T dry_radius, T wet_radius,
                           T& wet_refractive_index_real,
                           T& wet_refractive_index_imaginary)
  {
    wet_refractive_index_real = water_refractive_index_real
      + (dry_refractive_index_real - water_refractive_index_real)
      * wpower(dry_radius / wet_radius, T(3.));

    wet_refractive_index_imaginary = water_refractive_index_imaginary
      + (dry_refractive_index_imaginary - water_refractive_index_imaginary)
      * wpower(dry_radius / wet_radius, T(3.));
  }


}  // namespace AtmoData.

#define ATMODATA_FILE_DEPOSITION_CXX
#endif
