// Copyright (C) 2006-2007, ENPC - INRIA - EDF R&D
// Author(s): Irène Korsakissok
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


#ifndef ATMODATA_FILE_TIMEDIAGNOSIS_CXX


#include "TimeDiagnosis.hxx"


namespace AtmoData
{
  //! Computes the sun declination and the equation of time at a given date.
  /*!
    Calculates solar declination and equation of time for a given date.
    Calculation is based on equations given in:  Paltridge and Platt,
    Radiative Processes in Meteorology and Climatology, Elsevier,
    pp. 62,63, 1976.
    Fourier coefficients originally from:  Spencer, J.W., 1971, Fourier
    series representation of the position of the sun, Search, 2:172.
    Note:  This approximate program does not account for changes from year
    to year
    \param date date in the form YYYYMMDD.
    \param ut local time in decimal UT (e.g., 16.25 means 15 minutes
    after 4 pm).
    \param declination (output) sun declination (radians).
    \param time_equation (output) equation of time (radians).
  */
  template<class T>
  void ComputeDeclination(Date date, T ut, T& declination, T& time_equation)
  {
    T d, tz, sintz, costz, sin2tz, cos2tz, sin3tz, cos3tz;
    const T pi = 3.1415926535898;

    int j = date.GetNumberOfDays();

    // Calculates decimal Julian day from start of year.
    d = T(j) + ut / 24.;

    // Equation 3.8 for "day-angle" in radians.
    tz = 2.* pi * d / 365.;

    // Calculates sine and cosine from addition theoremes for
    // better performance;  the computation of sin2tz,
    // sin3tz, cos2tz and cos3tz is about 5-6 times faster
    // than the evaluation of the intrinsic functions.
    //
    // It is sin(x+y) = sin(x)*cos(y)+cos(x)*sin(y)
    // and   cos(x+y) = cos(x)*cos(y)-sin(x)*sin(y)
    //
    // sintz  = sin(tz)      costz  = cos(tz)
    // sin2tz = sin(2.*tz)   cos2tz = sin(2.*tz)
    // sin3tz = sin(3.*tz)   cos3tz = cos(3.*tz)

    sintz = sin(tz);
    costz = cos(tz);
    sin2tz = 2. * sintz * costz;
    cos2tz = costz * costz - sintz * sintz;
    sin3tz = sintz * cos2tz + costz * sin2tz;
    cos3tz = costz * cos2tz - sintz * sin2tz;

    // Equation 3.7 for declination in radians.
    declination = 0.006918 - 0.399912 * costz  + 0.070257 * sintz
      - 0.006758 * cos2tz + 0.000907 * sin2tz
      - 0.002697 * cos3tz + 0.001480 * sin3tz;

    // Equation 3.11 for Equation of time in radians.
    time_equation = 0.000075 + 0.001868 * costz  - 0.032077 * sintz
      - 0.014615 * cos2tz - 0.040849 * sin2tz;
  }


  //! Computes the sun declination and the equation of time at a given date.
  /*!
    Calculates solar declination and equation of time for a given date.
    Calculation is based on equations given in:  Paltridge and Platt,
    Radiative Processes in Meteorology and Climatology, Elsevier,
    pp. 62,63, 1976.
    Fourier coefficients originally from:  Spencer, J.W., 1971, Fourier
    series representation of the position of the sun, Search, 2:172.
    Note:  This approximate program does not account for changes from year
    to year
    \param idate date in the form YYYYMMDD.
    \param ut local time in decimal UT (e.g., 16.25 means 15 minutes
    after 4 pm).
    \param declination (output) sun declination (radians).
    \param time_equation (output) equation of time (radians).
  */
  template<class T>
  void ComputeDeclination(int idate, T ut, T& declination, T& time_equation)
  {
    ComputeDeclination(Date(idate), ut, declination, time_equation);
  }


  //! Computes the hours of sunrise and sunset at a given location and date.
  /*!
    \param lon longitude of location (degrees).
    \param lat latitude of location (degrees).
    \param idate date in the form YYYYMMDD.
    \param sunrise_hour (output) hour of sunrise in decimal UT
    (greenwich mean time).
    \param sunset_hour (output) hour of sunset in decimal UT
    (greenwich mean time).
  */
  template<class T>
  void ComputeSunHour(T lon, T lat, int idate,
                      T& sunrise_hour, T& sunset_hour)
  {
    T declination, time_equation, eqh;
    T csh,  h, lmtr, lmts;
    const T pi = 3.1415926535898;

    // Computes solar declination and equation of time (radians).
    ComputeDeclination(idate, T(0.), declination, time_equation);

    // Converts equation of time to hours.
    eqh = time_equation * 24. / (2. * pi);

    // Cosine of hour angle.
    csh = -tan(lat * pi / 180.) * tan(declination);
    if (csh > 1.)
      {
        sunrise_hour = 0.;
        sunset_hour = 0.;
      }
    else if (csh < -1.)
      {
        sunrise_hour = 1.;
        sunset_hour = 1.;
      }
    else
      {
        // Computes the hour angle in hours.
        h = acos(csh) * 24. / (2. * pi);

        // Computes the local mean time of sunrise and sunset in hours.
        lmtr = 12. - h - eqh;
        lmts = 12. + h - eqh;

        // Computes the universal time of sunrise and sunset in hours.
        sunrise_hour = lmtr - lon / 15.;
        sunset_hour = lmts - lon / 15.;
      }
  }


  //! Computes the hour of sunrise at a given location and date.
  /*!
    \param lon longitude of location (degrees).
    \param lat latitude of location (degrees).
    \param idate date in the form YYYYMMDD.
    \return the hour of sunrise in decimal UT.
  */
  template<class T>
  T ComputeSunriseHour(T lon, T lat, int idate)
  {
    T sunrise_hour, sunset_hour;
    ComputeSunHour(lon, lat, idate, sunrise_hour, sunset_hour);
    return sunrise_hour;
  }


  //! Computes the hour of sunset at a given location and date.
  /*!
    \param lon longitude of location (degrees).
    \param lat latitude of location (degrees).
    \param idate date in the form YYYYMMDD.
    \return the hour of sunset in decimal UT.
  */
  template<class T>
  T ComputeSunsetHour(T lon, T lat, int idate)
  {
    T sunrise_hour, sunset_hour;
    ComputeSunHour(lon, lat, idate, sunrise_hour, sunset_hour);
    return sunset_hour;
  }


  //! Checks whether it is nighttime or daytime at a given location and time.
  /*!
    \param lon longitude of location (degrees).
    \param lat latitude of location (degrees).
    \param idate date in the form YYYYMMDD.
    \param ut local time in decimal UT (e.g., 16.25 means 15 minutes
    after 4 pm).
    \return true if it is daytime, false otherwise.
  */
  template<class T>
  bool IsDay(T lon, T lat, int idate, T ut)
  {
    T sunrise_hour, sunset_hour;
    // Computes sunrise and sunset hours.
    ComputeSunHour(lon, lat, idate, sunrise_hour, sunset_hour);

    if (sunrise_hour != sunset_hour)
      return ut > sunrise_hour && ut < sunset_hour;
    else
      return sunrise_hour == 1.;
  }

  //! Checks whether it is nighttime or daytime at a given location and time.
  /*!
    \param lon longitude of location (degrees).
    \param lat latitude of location (degrees).
    \param date date and time (hour and minutes) in format Date.
    \return true if it is daytime, false otherwise.
  */
  template<class T>
  bool IsDay(T lon, T lat, Date date)
  {
    int idate, hr, mn;
    T sc, ut;

    // Gets current date in format YYYYMMDD.
    idate = date.GetDate();

    // Gets current hour in decimal UT.
    hr = date.GetHour();
    mn = date.GetMinutes();
    sc = date.GetSeconds();
    ut = T(hr) + T(mn) / 60. + T(sc) / 3600.;
    return IsDay(lon, lat, idate, ut);
  }

}  // namespace AtmoData.


#define ATMODATA_FILE_TIMEDIAGNOSIS_CXX
#endif
