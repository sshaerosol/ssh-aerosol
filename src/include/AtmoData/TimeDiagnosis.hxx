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


#ifndef ATMODATA_FILE_TIMEDIAGNOSIS_HXX


namespace AtmoData
{
  template<class T>
  void ComputeDeclination(Date date, T ut, T& declination, T& time_equation);

  template<class T>
  void ComputeDeclination(int idate, T ut, T& declination, T& time_equation);

  template<class T>
  void ComputeSunHour(T lon, T lat, int idate,
                      T& sunrise_hour, T& sunset_hour);

  template<class T>
  T ComputeSunriseHour(T lon, T lat, int idate);

  template<class T>
  T ComputeSunsetHour(T lon, T lat, int idate);

  template<class T>
  bool IsDay(T lon, T lat, int idate, T ut);

  template<class T>
  bool IsDay(T lon, T lat, Date date);


}  // namespace AtmoData.


#define ATMODATA_FILE_TIMEDIAGNOSIS_HXX
#endif
