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


#ifndef ATMODATA_FILE_KZ_CXX

#include "Kz.hxx"

namespace AtmoData
{

  //! Computes vertical diffusion according to the Louis formula (1979).
  /*!
    \param U zonal wind.
    \param V meridional wind.
    \param Tp potential temperature.
    \param Kz (output) vertical diffusion coefficients at the interfaces.
    \param L0 scale parameter. Default: 100.
    \param B parameter. Default: 5.
    \param C parameter. Default: 5.
    \param D parameter. Default: 5.
    \param z0 scale parameter. Default: 1.
    \param a parameter. Default: 0.115.
    \param b parameter. Default: 0.175.
    \param delta_z0 parameter. Default: 0.01.
    \param Ka Von Karman constant. Default: 0.4.
    \note Kz is given at the interfaces. It is assumed to be 0 at the first
    interface and at the last one.
    \note The critical Richardson number is computed following Pielke (2002)
    and Nordeng (1986): Ric = a * (delta_z / delta_z0)^b.
  */
  template<class TU, class TV, class TTp, class T, class TG>
  void ComputeLouisKz(Data<TU, 4, TG>& U, Data<TV, 4, TG>& V,
                      Data<TTp, 4, TG>& Tp, Data<T, 4, TG>& Kz,
                      T L0, T B, T C, T D, T z0, T a, T b, T delta_z0, T Ka)
  {

    int h, i, j, k;

    int Nx = Kz.GetLength(3);
    int Ny = Kz.GetLength(2);
    int Nz = Kz.GetLength(1) - 1;
    int Nt = min(U.GetLength(0), V.GetLength(0));

    T l, R, F, L;
    T dWind_dz, derivative, dz;
    const T g(9.81);

    Grid<TG>& Levels = Kz[1];
    Grid<TG>& Nodes = U[1];

    Kz.SetZero();

    for (h = 0; h < Nt; h++)
      for (k = 1; k < Nz; k++)
        for (j = 0; j < Ny; j++)
          for (i = 0; i < Nx; i++)
            {

              /*********************/
              /* l = Ka * ---- ... */

              L = L0;

              l = Ka * (Levels(k) + z0) / (1.0 + Ka * (Levels(k) + z0) / L);

              /* l = Ka * ---- ... */
              /*********************/

              /**************/
              /* dWind / dz */

              dz = Nodes(k) - Nodes(k - 1);

              // dU/dz.

              derivative = (U(h, k, j, i + 1) + U(h, k, j, i)
                            - U(h, k - 1, j, i + 1) - U(h, k - 1, j, i))
                / dz * 0.5;

              derivative = derivative * derivative;
              dWind_dz = derivative;

              // dV/dz.

              derivative = (V(h, k, j + 1, i) + V(h, k, j, i)
                            - V(h, k - 1, j + 1, i) - V(h, k - 1, j, i))
                / dz * 0.5;

              derivative = derivative * derivative;
              dWind_dz += derivative;

              dWind_dz = sqrt(dWind_dz);

              /* dWind / dz */
              /**************/

              /***********/
              /* F(R, z) */

              derivative = (Tp(h, k, j, i) - Tp(h, k - 1, j, i))
                / (dz * Tp(h, k - 1, j, i));

              R = min(g * derivative / (dWind_dz * dWind_dz),
                      a * pow(dz / delta_z0, b));

              if (R >= 0)
                F = 1.0 / (1.0 + 3.0 * B * R * sqrt(1.0 + D * R));
              else
                {
                  F = 1.0 + 3.0 * B * C * sqrt(fabs(R) / 27.0) *
                    (l * l) / ((Levels(k) + z0) * (Levels(k) + z0));
                  F = 1.0 - 3.0 * B * R / F;
                }

              /* F(R, z) */
              /***********/

              Kz(h, k, j, i) = l * l * dWind_dz * F;

            }

  }


  //! Computes the Monin-Obukhov length.
  /*!
    \param FrictionModule friction module.
    \param SurfacePotentialTemperature surface potential temperature.
    \param PotentialTemperature potential temperature.
    \param SensibleHeat sensible heat.
    \param Evaporation evaporation.
    \param LMO (output) Monin-Obukhov length.
    \param Ka Von Karman constant. Default: 0.4.
  */
  template<class TU, class TTp, class T, class TG>
  void ComputeLMO(const Data<TU, 3, TG>& FrictionModule,
                  const Data<TTp, 3, TG>& SurfacePotentialTemperature,
                  const Data<TTp, 4, TG>& PotentialTemperature,
                  const Data<T, 3, TG>& SensibleHeat,
                  const Data<T, 3, TG>& Evaporation,
                  Data<T, 3, TG>& LMO, T Ka)
  {

    int h, i, j;

    int Nx = LMO.GetLength(2);
    int Ny = LMO.GetLength(1);
    int Nt = LMO.GetLength(0);

    const T g(9.81);

    for (h = 0; h < Nt; h++)
      for (j = 0; j < Ny; j++)
        for (i = 0; i < Nx; i++)
          {
            float u_star = FrictionModule(h, j, i);
            float theta_m = 0.5 * (SurfacePotentialTemperature(h, j, i)
                                   + PotentialTemperature(h, 0, j, i));
            LMO(h, j, i) = - u_star * u_star * u_star * theta_m
              / (Ka * g * (SensibleHeat(h, j, i)
                           + 0.608 * theta_m * Evaporation(h, j, i)));
          }
  }



  //! Computes the boundary layer height according to Troen & Mahrt.
  /*!
    \param SurfaceTemperature surface temperature.
    \param SurfacePotentialTemperature surface potential temperature.
    \param PotentialTemperature potential temperature.
    \param FrictionModule friction module.
    \param WindModule wind module.
    \param SensibleHeat sensible heat.
    \param LMO Monin-Obukhov length.
    \param GridZ_interf altitudes of interfaces (m).
    \param BoundaryHeight (output) boundary layer height.
    \param SBL Ratio between the SBL and the PBL. Default: 0.1.
    \param Ric critical Richardson number. Default: 0.21.
    \param C Troen & Mahrt coefficient. Default: 6.5.
    \param Ka Von Karman constant. Default: 0.4.
  */
  template<class TT, class TTp, class TU, class T, class TG>
  void ComputePBLH_TM(const Data<TT, 3, TG>& SurfaceTemperature,
                      const Data<TTp, 3, TG>& SurfacePotentialTemperature,
                      const Data<TTp, 4, TG>& PotentialTemperature,
                      const Data<TU, 3, TG>& FrictionModule,
                      const Data<TU, 4, TG>& WindModule,
                      const Data<T, 3, TG>& SensibleHeat,
                      const Data<T, 3, TG>& LMO,
                      const Grid<TG>& GridZ_interf,
                      Data<T, 3, TG>& BoundaryHeight,
                      T SBL, T Ric, T C, T Ka)
  {

    int h, i, j, k;

    int Nx = BoundaryHeight.GetLength(2);
    int Ny = BoundaryHeight.GetLength(1);
    int Nt = BoundaryHeight.GetLength(0);

    int Nz = WindModule.GetLength(1);
    const Grid<TG>& Levels = WindModule[1];

    T l_MO, u_star, w_star, ws, qv0;
    T buoyancy_temp,
      atm_pot_temp, prev_pot_temp,
      curr_pot_temp(0);

    const T g(9.81);

    for (h = 0; h < Nt; h++)
      for (j = 0; j < Ny; j++)
        for (i = 0; i < Nx; i++)
          {

            u_star = FrictionModule(h, j, i);
            l_MO = LMO(h, j, i);
            qv0 = SensibleHeat(h, j, i);

            if (l_MO < 0.)
              w_star = pow(Levels(1)
                           * g / SurfaceTemperature(h, j, i) * qv0,
                           T(1. / 3.));
            else
              w_star = 0.;

            ws = pow(u_star * u_star * u_star + 7. * SBL * Ka
                     * w_star * w_star * w_star, 1. / 3.);

            buoyancy_temp = (SurfacePotentialTemperature(h, j, i)
                             + PotentialTemperature(h, 0, j, i)) / 2.0;
            prev_pot_temp = atm_pot_temp
              = PotentialTemperature(h, 0, j, i) + C * qv0 / ws;

            // Boundary-layer height.
            k = 1;
            while ((k < Nz - 1)
                   && (PotentialTemperature(h, k, j, i)
                       < (curr_pot_temp = atm_pot_temp
                          + Ric * WindModule(h, k, j, i)
                          * WindModule(h, k, j, i)
                          / (Levels(k) * g) * buoyancy_temp)))
              {
                k++;
                prev_pot_temp = curr_pot_temp;
                if (l_MO < 0.)
                  w_star = pow(Levels(k)
                               * g / SurfaceTemperature(h, j, i) * qv0,
                               T(1. / 3.));
                else
                  w_star = 0.;
                ws = pow(u_star * u_star * u_star + 7. * SBL * Ka
                         * w_star * w_star * w_star, 1. / 3.);
                atm_pot_temp = PotentialTemperature(h, 0, j, i)
                  + C * qv0 / ws;
              }
            if (k == 1)
              BoundaryHeight(h, j, i) = GridZ_interf(1);
            else
              BoundaryHeight(h, j, i) = Levels(k - 1)
                + (PotentialTemperature(h, k - 1, j, i) - prev_pot_temp)
                / (curr_pot_temp - prev_pot_temp
                   + PotentialTemperature(h, k - 1, j, i)
                   - PotentialTemperature(h, k, j, i))
                * (Levels(k) - Levels(k - 1));
          }
  }


  //! Computes the boundary layer height based on Richardson number.
  /*!
    \param Richardson Richardson number.
    \param GridZ_interf altitudes of interfaces (m).
    \param BoundaryHeight (output) boundary layer height.
    \param Ric critical Richardson number. Default: 0.21.
  */
  template<class T, class TG>
  void ComputePBLH_Richardson(const Data<T, 4, TG>& Richardson,
                              const Grid<TG>& GridZ_interf,
                              Data<T, 3, TG>& BoundaryHeight,
                              T Ric)
  {

    int h, i, j, k;

    int Nx = BoundaryHeight.GetLength(2);
    int Ny = BoundaryHeight.GetLength(1);
    int Nt = BoundaryHeight.GetLength(0);

    int Nz = Richardson.GetLength(1);
    const Grid<TG>& Levels = Richardson[1];

    for (h = 0; h < Nt; h++)
      for (j = 0; j < Ny; j++)
        for (i = 0; i < Nx; i++)
          {
            k = 0;
            while ((k < Nz - 1) && (Richardson(h, k, j, i) < Ric))
              k++;
            if (k == 0)
              BoundaryHeight(h, j, i) = GridZ_interf(1);
            else if (k == Nz - 1)
              BoundaryHeight(h, j, i) = Levels(k);
            else
              BoundaryHeight(h, j, i) = Levels(k - 1)
                + (Ric - Richardson(h, k - 1, j, i))
                * (Levels(k) - Levels(k - 1))
                / (Richardson(h, k, j, i) - Richardson(h, k - 1, j, i));
          }
  }


  //! Computes vertical diffusion according to the Troen and Mahrt (1986).
  /*!
    \param SurfaceTemperature surface temperature.
    \param FrictionModule friction module.
    \param SensibleHeat sensible heat.
    \param LMO Monin-Obukhov length.
    \param BoundaryHeight boundary layer height.
    \param Kz (output) vertical diffusion coefficients at the interfaces.
    \param TM_stable if set to 'false', Troen and Mahrt parameterization only
    applied in unstable boundary layer. Default: true.
    \param SBL Ratio between the SBL and the PBL. Default: 0.1.
    \param p Troen & Mahrt coefficient. Default: 2.0.
    \param Ka Von Karman constant. Default: 0.4.

  */
  template<class TU, class TT, class T, class TG>
  void ComputeTM_Kz(const Data<TT, 3, TG>& SurfaceTemperature,
                    const Data<TU, 3, TG>& FrictionModule,
                    const Data<T, 3, TG>& SensibleHeat,
                    const Data<T, 3, TG>& LMO,
                    const Data<T, 3, TG>& BoundaryHeight,
                    Data<T, 4, TG>& Kz,
                    bool TM_stable, T SBL, T p, T Ka)
  {

    int h, i, j, k;

    int Nx = Kz.GetLength(3);
    int Ny = Kz.GetLength(2);
    int Nz = Kz.GetLength(1) - 1;
    int Nt = Kz.GetLength(0);

    Grid<TG>& GridZ_interf = Kz[1];

    const T g(9.81);

    T l_MO, u_star, w_star, ws, qv0, boundary_height, Phi_m;

    for (h = 0; h < Nt; h++)
      for (j = 0; j < Ny; j++)
        for (i = 0; i < Nx; i++)
          {

            u_star = FrictionModule(h, j, i);
            l_MO = LMO(h, j, i);
            qv0 = SensibleHeat(h, j, i);
            boundary_height = BoundaryHeight(h, j, i);

            // Computes ws according to the new boundary height.
            if (l_MO < 0.)
              w_star = pow(boundary_height * g / SurfaceTemperature(h, j, i) *
                           abs(qv0), T(1. / 3.));
            else
              w_star = 0.;
            ws = pow(u_star * u_star * u_star + 7. * SBL * Ka
                     * w_star * w_star * w_star, 1. / 3.);

            // Computes the nondimensional shear Phi_m and Kz.
            if (l_MO < 0)
              for (k = 1; k < Nz; k++)
                {
                  if (GridZ_interf(k) < SBL * boundary_height)
                    Phi_m = pow(1.0 - 7.0 * GridZ_interf(k) / l_MO, -1.0 / 3.0);
                  else
                    Phi_m = u_star / ws;
                  if (GridZ_interf(k) < boundary_height)
                    Kz(h, k, j, i) = u_star * Ka * GridZ_interf(k) / Phi_m
                      * pow(T(1.0) - GridZ_interf(k) / boundary_height, p);
                }
            else if (TM_stable)
              for (k = 1; k < Nz; k++)
                if (GridZ_interf(k) < boundary_height)
                  {
                    Phi_m = 1.0 + 4.7 * GridZ_interf(k) / l_MO;
                    Kz(h, k, j, i) = u_star * Ka * GridZ_interf(k) / Phi_m
                      * pow(T(1.0) - GridZ_interf(k) / boundary_height, p);
                  }

          }
  }

}  // namespace AtmoData.

#define ATMODATA_FILE_KZ_CXX
#endif
