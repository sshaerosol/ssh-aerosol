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


#ifndef ATMODATA_FILE_PHOTOLYSIS_HXX

namespace AtmoData
{

  template<class T>
  T ZenithAngle(T lon, T lat, Date date, T ut);

  template <class TH, class TL, class TMC, class THC, class T, class TG>
  void ComputeAttenuation_LWC(Data<TH, 4, TG>& Humidity, Data<TH, 4, TG>& CRH,
                              Data<TL, 4, TG>& LiquidWaterContent,
                              Data<TMC, 3, TG>& MediumCloudiness,
                              Data<THC, 3, TG>& HighCloudiness,
                              Date date_beg, T Delta_t,
                              Data<T, 4, TG>& Attenuation);

  template <class TL, class TMC, class THC, class T, class TG>
  void ComputeAttenuation_LWC(Data<TL, 4, TG>& LiquidWaterContent,
                              Data<int, 4> LowIndices,
                              Data<int, 4> MediumIndices,
                              Data<int, 4> HighIndices,
                              Data<TMC, 3, TG>& MediumCloudiness,
                              Data<THC, 3, TG>& HighCloudiness,
                              Date date_beg, T Delta_t,
                              Data<T, 4, TG>& Attenuation);

  template <class TMC, class THC, class TG, class TH, class T>
  void ComputeAttenuation_ESQUIF(Data<TMC, 3, TG>& MediumCloudiness,
                                 Data<THC, 3, TG>& HighCloudiness,
                                 Data<TH, 4, TG>& RelativeHumidity,
                                 Data<T, 4, TG>& Attenuation,
                                 T a = 0.1, T b = 0.3, T c = 1.5);

}  // namespace AtmoData.

#define ATMODATA_FILE_PHOTOLYSIS_HXX
#endif
