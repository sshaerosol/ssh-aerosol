// Copyright (C) 2007, ENPC - INRIA - EDF R&D
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


#ifndef ATMODATA_FILE_PHOTOLYSIS_CXX

#include "Photolysis.hxx"

#include "Talos.hxx"
using namespace Talos;

namespace AtmoData
{


  // The following function is derived from the Fortran (77) subroutine
  // zenith found in TUV (Tropospheric Ultraviolet & Visible radiation model).
  // This subroutine was provided under the GNU General Public License
  // under the following copyright:
  // Copyright (C) 1994,95,96  University Corporation for Atmospheric
  // Research.
  /*! Calculates solar zenith angle for a given time and location.
    Calculation is based on equations given in:  Paltridge and Platt,
    Radiative Processes in Meteorology and Climatology, Elsevier,
    pp. 62,63, 1976.
    Fourier coefficients originally from:  Spencer, J.W., 1971, Fourier
    series representation of the position of the sun, Search, 2:172.
    Note:  This approximate program does not account for changes from year
    to year.
    \param lon longitude of location (degrees).
    \param lat latitude of location (degrees).
    \param date date in the form YYYYMMDD.
    \param ut local time in decimal UT (e.g., 16.25 means 15 minutes
    after 4 pm).
    \return solar zenith angle (degrees).
  */
  template<class T>
  T ZenithAngle(T lon, T lat, Date date, T ut)
  {

    T zen;
    // T azim, caz, raz;

    T lbut, lzut;
    T rlt;
    T rdecl, eqr, eqh, zpt;
    T csz, zr;

    const T pi = 3.1415926535898;
    const T dr = double(pi) / double(180.);

    // Converts latitude to radians.
    rlt = lat * dr;

    //Computes solar declination and equation of time (radians).
    ComputeDeclination(date, ut, rdecl, eqr);

    // Converts equation of time to hours.
    eqh = eqr * 24. / (2. * pi);

    // Calculates local hour angle (hours).
    lbut = 12. - eqh - lon * 24. / 360.;

    // Converts to angle from UT.
    lzut = 15. * (ut - lbut);
    zpt = lzut * dr;

    // Equation 2.4 for cosine of zenith angle.
    csz = sin(rlt) * sin(rdecl) + cos(rlt) * cos(rdecl) * cos(zpt);
    zr = acos(csz);
    zen = zr / dr;

    // Calculates local solar azimuth.
    // caz = (sin(rdecl) - sin(rlt) * cos(zr)) / (cos(rlt) * sin(zr));
    // raz = acos(caz);
    // azim = raz / dr;

    return zen;

  }


  //! Computes the cloud attenuation for photolysis rates.
  /*!
    \param Humidity relative humidity.
    \param CRH critical relative humidity.
    \param LiquidWaterContent liquid water content (kg/m^3).
    \param MediumCloudiness medium cloudiness (in [0, 1]).
    \param HighCloudiness high cloudiness (in [0, 1]).
    \param date_beg beginning date for the computation of attenuation.
    \param Delta_t time step for the computation of attenuation.
    \param Attenuation (output) cloud attenuation coefficient.
  */
  template <class TH, class TL, class TMC, class THC, class T, class TG>
  void ComputeAttenuation_LWC(Data<TH, 4, TG>& Humidity, Data<TH, 4, TG>& CRH,
                              Data<TL, 4, TG>& LiquidWaterContent,
                              Data<TMC, 3, TG>& MediumCloudiness,
                              Data<THC, 3, TG>& HighCloudiness,
                              Date date_beg, T Delta_t,
                              Data<T, 4, TG>& Attenuation)
  {
    int h, k, j, i;
    int Nt(Attenuation.GetLength(0));
    int Nz(Attenuation.GetLength(1));
    int Ny(Attenuation.GetLength(2));
    int Nx(Attenuation.GetLength(3));

    Date current_date = date_beg;

    // Index "0" and "1" refer to two contiguous levels.
    T rh0, rh1, rhc, dz, delta_z,
      lwc0, lwc1, lw, w, tau, tr;

    for (h = 0; h < Nt; h++)
      {
        for (j = 0; j < Ny; j++)
          for (i = 0; i < Nx; i++)
            {
              rh0 = Humidity(h, Nz - 1, j, i);
              // kg/m^3 to g/m^3.
              lwc0 = 1000. * LiquidWaterContent(h, Nz - 1, j, i);

              w = 0;

              for (k = Nz - 1; k >= 0; k--)
                {
                  if (k == Nz - 1)
                    dz = Attenuation[1].Value(h, Nz - 1, j, i)
                      - Attenuation[1].Value(h, Nz - 2, j, i);
                  else
                    dz = Attenuation[1].Value(h, k + 1, j, i)
                      - Attenuation[1].Value(h, k, j, i);

                  rh1 = Humidity(h, k, j, i);
                  // kg/m^3 to g/m^3.
                  lwc1 = 1000. * LiquidWaterContent(h, k, j, i);

                  // Critical relative humidity.
                  rhc = CRH(h, k, j, i);

                  if (rh0 > rhc && rh1 > rhc)  // In a cloud.
                    w += dz * (lwc0 + lwc1) / 2.0;
                  else if (rh0 > rhc && rh1 < rhc)  // Below a cloud.
                    {
                      delta_z = dz * (rh0 - rhc) / (rh0 - rh1);
                      w += lwc1 * delta_z + (lwc0 - lwc1) / dz * .5
                        * (2.0 * dz - delta_z) * delta_z;
                    }
                  else if (rh0 < rhc && rh1 > rhc)  // Above a cloud.
                    {
                      delta_z = dz * (rh1 - rhc) / (rh1 - rh0);
                      w += lwc1 * delta_z + (lwc0 - lwc1) / dz * .5
                        * delta_z * delta_z;
                    }
                  // For the next level.
                  rh0 = rh1;
                  lwc0 = lwc1;

                  // Computes liquid water path.
                  if (w > 0.)
                    lw = log10(w);
                  else
                    lw = 0.;

                  // Computes the cloud optical depth according
                  // to Stephens (1978).
                  if (lw <= 0.)
                    tau = 0.;
                  else
                    tau = pow(10., 0.2633 + 1.7095 * log(lw));

                  // Computes the cloud transmissivity.
                  if (tau < 5.)
                    tr = 1.;
                  else
                    tr = (5. - exp(-tau)) / (4. + 0.42 * tau);

                  // Zenith angle.
                  T cos_zenith_angle;
                  cos_zenith_angle =
                    cos(ZenithAngle(Attenuation[3].Value(h, k, j, i),
                                    Attenuation[2].Value(h, k, j, i),
                                    current_date,
                                    Attenuation[0].Value(h, k, j, i))
                        * 0.0174532925199433);
                  cos_zenith_angle = abs(cos_zenith_angle);

                  // Computes the attenuation coefficient.
                  if (tr == 1)
                    Attenuation(h, k, j, i) = 1.0
                      + (min(1.0, MediumCloudiness(h, j, i)
                             + HighCloudiness(h, j, i)))
                      * (1.6 * cos_zenith_angle - 1.0);
                  else
                    Attenuation(h, k, j, i) = 1.0
                      + (min(1.0, MediumCloudiness(h, j, i)
                             + HighCloudiness(h, j, i)))
                      * ((1. - tr) * cos_zenith_angle);
                }
            }
        current_date.AddHours(int(Delta_t));
      }
  }


  //! Computes the cloud attenuation for photolysis rates.
  /*!
    \param LiquidWaterContent liquid water content (kg/m^3).
    \param LowIndices vertical indices of base and top of low clouds.
    \param MediumIndices vertical indices of base and top of medium clouds.
    \param HighIndices vertical indices of base and top of high clouds.
    \param MediumCloudiness medium cloudiness (in [0, 1]).
    \param HighCloudiness high cloudiness (in [0, 1]).
    relative humidity as function of the altitude, the pressure and
    reference pressure.
    \param date_beg beginning date for the computation of attenuation.
    \param Delta_t time step for the computation of attenuation.
    \param Attenuation (output) cloud attenuation coefficient.
    \note Dimensions of LowIndices, MediumIndices and HighIndices are
    Nt x Ny x Nx x 2. Along the last dimension, those arrays store the index
    of the cloud base and the index of the cloud top (in this order). Those
    indices are indices of interfaces. E.g., if LowIndices(t, y, x, 0) equals
    2 and  LowIndices(t, y, x, 1) equals 4, then a cloud lies in layers 2
    and 3.
  */
  template <class TL, class TMC, class THC, class T, class TG>
  void ComputeAttenuation_LWC(Data<TL, 4, TG>& LiquidWaterContent,
                              Data<int, 4> LowIndices,
                              Data<int, 4> MediumIndices,
                              Data<int, 4> HighIndices,
                              Data<TMC, 3, TG>& MediumCloudiness,
                              Data<THC, 3, TG>& HighCloudiness,
                              Date date_beg, T Delta_t,
                              Data<T, 4, TG>& Attenuation)
  {
    int h, k, j, i;
    int Nt(Attenuation.GetLength(0));
    int Nz(Attenuation.GetLength(1));
    int Ny(Attenuation.GetLength(2));
    int Nx(Attenuation.GetLength(3));

    Date current_date = date_beg;

    Attenuation.Fill(1.);

    // Indices "0" and "1" refer to two contiguous levels.
    T dz, lw, w, tau, tr(0.), low_att, up_att, dist, href;
    // "l": low, "m": medium, "h": high.
    // "b": bottom, "t": top.
    int lb, lt, mb, mt, hb, ht, lower, upper;
    T cos_zenith_angle, alpha;

    for (h = 0; h < Nt; h++)
      {
        for (j = 0; j < Ny; j++)
          for (i = 0; i < Nx; i++)
            {

              // Low clouds are between lb and lt.
              lb = LowIndices(h, j, i, 0);
              lt = LowIndices(h, j, i, 1);
              // Medium clouds are between mb and mt.
              mb = MediumIndices(h, j, i, 0);
              mt = MediumIndices(h, j, i, 1);
              // High clouds are between hb and ht.
              hb = HighIndices(h, j, i, 0);
              ht = HighIndices(h, j, i, 1);

              // Searches for the cloud basis.
              lower = -1;
              if (lb == 0)
                if (mb == 0)
                  lower = hb;
                else
                  lower = mb;
              else
                lower = lb;

              // Searches for the cloud top.
              upper = max(lt, mt);
              upper = max(upper, ht);

              w = 0;

              // Starting from the top.
              for (k = Nz - 1; k >= 0; k--)
                {
                  if (k == Nz - 1)
                    dz = Attenuation[1].Value(h, Nz - 1, j, i)
                      - Attenuation[1].Value(h, Nz - 2, j, i);
                  else
                    dz = Attenuation[1].Value(h, k + 1, j, i)
                      - Attenuation[1].Value(h, k, j, i);

                  if ((k >= lb && k < lt) || (k >= mb && k < mt)
                      || (k >= hb && k < ht))  // In a cloud.
                    // kg/m^3 to g/m^3.
                    w += dz * 1000. * LiquidWaterContent(h, k, j, i);
                }

              // Computes liquid water path.
              if (w > 0.)
                lw = log10(w);
              else
                lw = 0.;

              // Computes the cloud optical depth according to Stephens
              // (1978).
              if (lw <= 0.)
                tau = 0.;
              else
                tau = pow(10., 0.2633 + 1.7095 * log(lw));

              // Computes the cloud transmissivity.
              if (tau > 5.)
                tr = (5. - exp(-tau)) / (4. + 0.42 * tau);

              /*** Below clouds ***/
              // If tau <= 5., nothing is done.
              for (k = 0; k < lower && tau > 5.; k++)
                {
                  // Zenith angle.
                  cos_zenith_angle =
                    cos(ZenithAngle(Attenuation[3].Value(h, k, j, i),
                                    Attenuation[2].Value(h, k, j, i),
                                    current_date,
                                    Attenuation[0].Value(h, k, j, i))
                        * 0.0174532925199433);
                  cos_zenith_angle = abs(cos_zenith_angle);

                  // Computes the attenuation coefficient.
                  Attenuation(h, k, j, i) = 1.0 + MediumCloudiness(h, j, i)
                    * (1.6 * tr * cos_zenith_angle - 1.0);
                }

              /*** Above clouds ***/
              // If tau <= 5., nothing is done.
              // Note: the loop is safe because if upper == -1, then tau == 0.
              for (k = upper; k < Nz && tau > 5.; k++)
                {
                  // Zenith angle.
                  cos_zenith_angle =
                    cos(ZenithAngle(Attenuation[3].Value(h, k, j, i),
                                    Attenuation[2].Value(h, k, j, i),
                                    current_date,
                                    Attenuation[0].Value(h, k, j, i))
                        * 0.0174532925199433);
                  cos_zenith_angle = abs(cos_zenith_angle);

                  // Computes the attenuation coefficient.
                  Attenuation(h, k, j, i) = 1.0
                    + MediumCloudiness(h, j, i) * (1. - tr)
                    * cos_zenith_angle;
                }

              /*** In cloud ***/
              // Computes the cloud thickness and the attenuation at both
              // ends.
              // Takes into account special cases in which the cloud
              // reaches the bottom or the top.
              if (upper == Nz && lower != 0)
                {
                  href = Attenuation[1].Value(h, lower - 1, j, i);
                  dist = 2. * Attenuation[1].Value(h, upper - 1, j, i)
                    - Attenuation[1].Value(h, upper - 2, j, i) - href;
                  up_att = 1.;
                  low_att = Attenuation(h, lower - 1, j, i);
                }
              else if (upper == Nz && lower == 0)
                {
                  href = Attenuation[1].Value(h, lower, j, i);
                  dist = 2. * Attenuation[1].Value(h, upper - 1, j, i)
                    - Attenuation[1].Value(h, upper - 2, j, i) - href;
                  up_att = 1.;
                  // Zenith angle.
                  cos_zenith_angle =
                    cos(ZenithAngle(Attenuation[3].Value(h, 0, j, i),
                                    Attenuation[2].Value(h, 0, j, i),
                                    current_date,
                                    Attenuation[0].Value(h, 0, j, i))
                        * 0.0174532925199433);
                  cos_zenith_angle = abs(cos_zenith_angle);
                  // Computes the attenuation coefficient.
                  low_att = 1.0 + MediumCloudiness(h, j, i)
                    * (1.6 * tr * cos_zenith_angle - 1.0);
                }
              else if (lower != 0)
                {
                  href = Attenuation[1].Value(h, lower - 1, j, i);
                  dist = Attenuation[1].Value(h, upper, j, i) - href;
                  up_att = Attenuation(h, upper, j, i);
                  low_att = Attenuation(h, lower - 1, j, i);
                }
              else
                {
                  href = Attenuation[1].Value(h, lower, j, i);
                  dist = Attenuation[1].Value(h, upper, j, i) - href;
                  up_att = Attenuation(h, upper, j, i);
                  // Zenith angle.
                  cos_zenith_angle =
                    cos(ZenithAngle(Attenuation[3].Value(h, 0, j, i),
                                    Attenuation[2].Value(h, 0, j, i),
                                    current_date,
                                    Attenuation[0].Value(h, 0, j, i))
                        * 0.0174532925199433);
                  cos_zenith_angle = abs(cos_zenith_angle);
                  // Computes the attenuation coefficient.
                  low_att = 1.0 + MediumCloudiness(h, j, i)
                    * (1.6 * tr * cos_zenith_angle - 1.0);
                }
              // If tau <= 5., nothing is done.
              for (k = lower; k < upper && tau > 5.; k++)
                {
                  // Computes the attenuation coefficient.
                  alpha = (Attenuation[1].Value(h, k, j, i) - href) / dist;
                  Attenuation(h, k, j, i) = alpha * up_att
                    + (1. - alpha) * low_att;
                }
            }
        current_date.AddHours(int(Delta_t));
      }
  }


  //! Computes the cloud attenuation for photolysis rates, following the
  // parameterization of the ESQUIF project (2001).
  /*!
    Formula : Attenuation = (1 - a * HighCloudiness) *
    (1 - b * MediumCloudiness) * exp(-c * B)
    \param MediumCloudiness medium cloudiness (in [0, 1]).
    \param HighCloudiness high cloudiness (in [0, 1]).
    \param RelativeHumidity relative humidity (kg/kg).
    \param Attenuation (output) cloud attenuation coefficient.
    \param a coefficient (see the formula). Default : 0.1.
    \param b coefficient (see the formula). Default : 0.3.
    \param c coefficient (see the formula). Default : 1.5.
  */
  template <class TMC, class THC, class TG, class TH, class T>
  void ComputeAttenuation_ESQUIF(Data<TMC, 3, TG>& MediumCloudiness,
                                 Data<THC, 3, TG>& HighCloudiness,
                                 Data<TH, 4, TG>& RelativeHumidity,
                                 Data<T, 4, TG>& Attenuation,
                                 T a, T b, T c)
  {
    int Nt(Attenuation.GetLength(0));
    int Nz(Attenuation.GetLength(1));
    int Ny(Attenuation.GetLength(2));
    int Nx(Attenuation.GetLength(3));
    int h, k, j, i;
    float dz;

    T B(0.), norm(0.);

    for (h = 0; h < Nt ; h++)
      for (j = 0; j < Ny ; j++)
        for (i = 0; i < Nx ; i++)
          {
            // Calculation of B.
            for (k = 0 ; k < Nz && RelativeHumidity[1].Value(h, k, j, i)
                   < 1500 ; k++)
              {
                dz = RelativeHumidity[1].Value(h, k + 1, j, i)
                  - RelativeHumidity[1].Value(h, k, j, i);

                if (RelativeHumidity(h, k, j, i) > 0.7)
                  B += (RelativeHumidity(h, k, j, i) - 0.7) * dz;

                norm += (1. - 0.7) * dz;
              }
            // Normalization.
            B /= norm;

            for (k = 0 ; k < Nz ; k++)
              Attenuation(h, k, j, i) = (1. - a * HighCloudiness(h, j, i))
                * (1. - b * MediumCloudiness(h, j, i))
                * exp(-c * B);
          }
  }

  // Computes the extinction for optical depth.
  // (OD is the integration of the extinction over the vertical).
  // The formula is based on  Pozzoli et al., 08 (JGR),
  // himself based on based on Rockel et al., 1991.
  /*!
    \param LiquidWaterContent liquid water content (kg/m^3).
    \param IceWaterContent Ice water content (kg/m^3).
    \param LiquidWaterExtinction (output) liquid cloud extinction .
    \param IceOpticalDepth (output) ice cloud extinction .
  */
  template <class TL, class T, class TG>
  void ComputeExtinction(const Data<TL, 4, TG>& LiquidWaterContent,
                         const Data<TL, 4, TG>& IceWaterContent,
                         Data<T, 4, TG>& LiquidWaterExtinction,
                         Data<T, 4, TG>& IceWaterExtinction)
  {
    double cf;

    int Nt(IceWaterContent.GetLength(0));
    int Nz(IceWaterContent.GetLength(1));
    int Ny(IceWaterContent.GetLength(2));
    int Nx(IceWaterContent.GetLength(3));

    LiquidWaterExtinction.Fill(0.);
    IceWaterExtinction.Fill(0.);

    T al = 1.488;
    T bl = -0.9374;
    T rleff = 12;
    T ai = 1.911;
    T bi = -1.0631;
    T rieff = 50;

    for (int h = 0; h < Nt; h++)
      for (int j = 0; j < Ny; j++)
        for (int i = 0; i < Nx; i++)
          // Starting from the top.
          for (int k = 0; k < Nz; k++)
            {
              // kg/m^3 to g/m^3.
              T tau = al * 1000. * LiquidWaterContent(h, k, j, i) * pow(rleff, bl);
              T icetau = ai * 1000. * IceWaterContent(h, k, j, i) * pow(rieff, bi);

              LiquidWaterExtinction(h, k, j, i) = tau;
              IceWaterExtinction(h, k, j, i) = icetau;
            }
  }

}  // namespace AtmoData.

#define ATMODATA_FILE_PHOTOLYSIS_CXX
#endif
