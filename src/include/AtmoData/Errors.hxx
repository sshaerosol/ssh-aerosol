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


#ifndef ATMODATA_FILE_ERRORS_HXX


namespace AtmoData
{

  // NGE.

  template < class T_ref, int N, class TG_ref,
             class T_comp, class TG_comp >
  T_ref NGE(Data<T_ref, N, TG_ref> data_ref,
            Data<T_comp, N, TG_comp>& data_comp,
            Function_Base<T_ref, bool>& test);

  // Bias.

  template < class T_ref, int N, class TG_ref,
             class T_comp, class TG_comp >
  T_ref Bias(Data<T_ref, N, TG_ref> data_ref,
             Data<T_comp, N, TG_comp>& data_comp,
             Function_Base<T_ref, bool>& test);

  // RMS.

  template < class T_ref, int N, class TG_ref,
             class T_comp, class TG_comp >
  T_ref RMS(Data<T_ref, N, TG_ref> data_ref,
            Data<T_comp, N, TG_comp>& data_comp,
            Function_Base<T_ref, bool>& test);

  // RelativeRMS.

  template < class T_ref, int N, class TG_ref,
             class T_comp, class TG_comp >
  T_ref RelativeRMS(Data<T_ref, N, TG_ref> data_ref,
                    Data<T_comp, N, TG_comp>& data_comp,
                    Function_Base<T_ref, bool>& test);

}  // namespace AtmoData.


#define ATMODATA_FILE_ERRORS_HXX
#endif
