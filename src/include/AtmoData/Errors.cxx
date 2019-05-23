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


#ifndef ATMODATA_FILE_ERRORS_CXX

#include "Errors.hxx"

namespace AtmoData
{

  /*******
   * NGE *
   *******/

  //! Computes the normalized gross error between two data sets.
  /*!
    \param data_ref reference data.
    \param data_comp data to be compared with 'data_ref'.
    \param test boolean function with two parameters; the i-th
    elements of 'data_ref' and 'data_comp' are taken into account
    if 'test(data_ref(i), data_comp(i))' is true.
    \return The normalized gross error.
  */
  template < class T_ref, int N, class TG_ref,
             class T_comp, class TG_comp >
  T_ref NGE(Data<T_ref, N, TG_ref> data_ref,
            Data<T_comp, N, TG_comp>& data_comp,
            Function_Base<T_ref, bool>& test)
  {
    T_ref nge;

    T_ref* data_ref_arr = data_ref.GetData();
    T_comp* data_comp_arr = data_comp.GetData();
    int NbElements = data_ref.GetNbElements();

#ifdef SELDONDATA_DEBUG_CHECK_DIMENSIONS

    if (NbElements != data_comp.GetNbElements())
      throw WrongDim("AtmoData::NGE(Data<T_ref, " + to_str(N) +
                     ">&, Data<T_comp, " + to_str(N) +
                     ">&, Function_Base<T_ref, bool>&)",
                     "Data sizes differ.");

#endif

    int nb_elt = 0;
    nge = T_ref(0);
    for (int i = 0; i < NbElements; i++)
      if (test(data_ref_arr[i], data_comp_arr[i]))
        {
          nb_elt++;
          nge += abs((data_ref_arr[i] - data_comp_arr[i])
                     / data_ref_arr[i]);
        }
    nge = nge / T_ref(nb_elt);

    return nge;
  }


  /********
   * BIAS *
   ********/

  //! Computes the bias between two data sets.
  /*!
    \param data_ref reference data.
    \param data_comp data to be compared with 'data_ref'.
    \param test boolean function with two parameters; the i-th
    elements of 'data_ref' and 'data_comp' are taken into account
    if 'test(data_ref(i), data_comp(i))' is true.
    \return The bias.
  */
  template < class T_ref, int N, class TG_ref,
             class T_comp, class TG_comp >
  T_ref Bias(Data<T_ref, N, TG_ref> data_ref,
             Data<T_comp, N, TG_comp>& data_comp,
             Function_Base<T_ref, bool>& test)
  {
    T_ref bias;

    T_ref* data_ref_arr = data_ref.GetData();
    T_comp* data_comp_arr = data_comp.GetData();
    int NbElements = data_ref.GetNbElements();

#ifdef SELDONDATA_DEBUG_CHECK_DIMENSIONS

    if (NbElements != data_comp.GetNbElements())
      throw WrongDim("AtmoData::Bias(Data<T_ref, " + to_str(N) +
                     ">&, Data<T_comp, " + to_str(N) +
                     ">&, Function_Base<T_ref, bool>&)",
                     "Data sizes differ.");

#endif

    int nb_elt = 0;
    bias = T_ref(0);
    for (int i = 0; i < NbElements; i++)
      if (test(data_ref_arr[i], data_comp_arr[i]))
        {
          nb_elt++;
          bias += data_ref_arr[i] - data_comp_arr[i];
        }
    bias = bias / T_ref(nb_elt);

    return bias;
  }


  /*******
   * RMS *
   ********/

  //! Computes the root mean square between two data sets.
  /*!
    \param data_ref reference data.
    \param data_comp data to be compared with 'data_ref'.
    \param test boolean function with two parameters; the i-th
    elements of 'data_ref' and 'data_comp' are taken into account
    if 'test(data_ref(i), data_comp(i))' is true.
    \return The root mean square.
  */
  template < class T_ref, int N, class TG_ref,
             class T_comp, class TG_comp >
  T_ref RMS(Data<T_ref, N, TG_ref> data_ref,
            Data<T_comp, N, TG_comp>& data_comp,
            Function_Base<T_ref, bool>& test)
  {
    T_ref rms(0);

    T_ref* data_ref_arr = data_ref.GetData();
    T_comp* data_comp_arr = data_comp.GetData();
    int NbElements = data_ref.GetNbElements();

#ifdef SELDONDATA_DEBUG_CHECK_DIMENSIONS

    if (NbElements != data_comp.GetNbElements())
      throw WrongDim("AtmoData::RMS(Data<T_ref, " + to_str(N) +
                     ">&, Data<T_comp, " + to_str(N) +
                     ">&, Function_Base<T_ref, bool>&)",
                     "Data sizes differ.");

#endif

    int nb_elt = 0;
    for (int i = 0; i < NbElements; i++)
      if (test(data_ref_arr[i], data_comp_arr[i]))
        {
          nb_elt++;
          rms += (data_ref_arr[i] - data_comp_arr[i])
            * (data_ref_arr[i] - data_comp_arr[i]);
        }
    rms = sqrt(rms / T_ref(nb_elt));

    return rms;
  }


  /***************
   * RelativeRMS *
   ***************/

  //! Computes the relative root mean square between two data sets.
  /*!
    \param data_ref reference data.
    \param data_comp data to be compared with 'data_ref'.
    \param test boolean function with two parameters; the i-th
    elements of 'data_ref' and 'data_comp' are taken into account
    if 'test(data_ref(i), data_comp(i))' is true.
    \return The relative root mean square.
  */
  template < class T_ref, int N, class TG_ref,
             class T_comp, class TG_comp >
  T_ref RelativeRMS(Data<T_ref, N, TG_ref> data_ref,
                    Data<T_comp, N, TG_comp>& data_comp,
                    Function_Base<T_ref, bool>& test)
  {
    T_ref relative_rms(0);

    T_ref* data_ref_arr = data_ref.GetData();
    T_comp* data_comp_arr = data_comp.GetData();
    int NbElements = data_ref.GetNbElements();

#ifdef SELDONDATA_DEBUG_CHECK_DIMENSIONS

    if (NbElements != data_comp.GetNbElements())
      throw WrongDim("AtmoData::RelativeRMS(Data<T_ref, " + to_str(N) +
                     ">&, Data<T_comp, " + to_str(N) +
                     ">&, Function_Base<T_ref, bool>&)",
                     "Data sizes differ.");

#endif

    int nb_elt = 0;
    for (int i = 0; i < NbElements; i++)
      if (test(data_ref_arr[i], data_comp_arr[i]))
        {
          nb_elt++;
          relative_rms += (data_ref_arr[i] - data_comp_arr[i])
            * (data_ref_arr[i] - data_comp_arr[i])
            / (data_ref_arr[i] * data_ref_arr[i]);
        }
    relative_rms = sqrt(relative_rms / T_ref(nb_elt));

    return relative_rms;
  }


  /********
   * CORR *
   ********/

  //! Computes the correlation between two data sets.
  /*!
    \param data_ref reference data.
    \param data_comp data to be compared with 'data_ref'.
    \param test boolean function with two parameters; the i-th
    elements of 'data_ref' and 'data_comp' are taken into account
    if 'test(data_ref(i), data_comp(i))' is true.
    \return The correlation.
  */
  template < class T_ref, int N, class TG_ref,
             class T_comp, class TG_comp >
  T_ref Corr(Data<T_ref, N, TG_ref> data_ref,
             Data<T_comp, N, TG_comp>& data_comp,
             Function_Base<T_ref, bool>& test)
  {
    T_ref corr;

    T_ref* data_ref_arr = data_ref.GetData();
    T_comp* data_comp_arr = data_comp.GetData();
    int NbElements = data_ref.GetNbElements();

#ifdef SELDONDATA_DEBUG_CHECK_DIMENSIONS

    if (NbElements != data_comp.GetNbElements())
      throw WrongDim("AtmoData::Corr(Data<T_ref, " + to_str(N) +
                     ">&, Data<T_comp, " + to_str(N) +
                     ">&, Function_Base<T_ref, bool>&)",
                     "Data sizes differ.");

#endif

    int nb_elt = 0;
    corr = T_ref(0);
    T_ref mean_ref = T_ref(0);
    T_comp mean_comp = T_comp(0);
    T_ref var_ref = T_ref(0);
    T_comp var_comp = T_comp(0);
    T_ref covar = T_ref(0);
    T_ref temp_ref;
    T_comp temp_comp;

    // Means.
    for (int i = 0; i < NbElements; i++)
      if (test(data_ref_arr[i], data_comp_arr[i]))
        {
          nb_elt++;
          mean_ref += data_ref_arr[i];
          mean_comp += data_comp_arr[i];
        }
    mean_ref = mean_ref / T_ref(nb_elt);
    mean_comp = mean_comp / T_comp(nb_elt);

    // Co-variances.
    for (int i = 0; i < NbElements; i++)
      if (test(data_ref_arr[i], data_comp_arr[i]))
        {
          temp_ref = data_ref_arr[i] - mean_ref;
          temp_comp = data_comp_arr[i] - mean_comp;
          covar += temp_ref * temp_comp;
          var_ref += temp_ref * temp_ref;
          var_comp += temp_comp * temp_comp;
        }

    corr = covar / sqrt(var_ref * var_comp);

    return corr;
  }


}  // namespace AtmoData.


#define ATMODATA_FILE_ERRORS_CXX
#endif
