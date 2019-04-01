// Copyright (C) 2007-2014, ENPC - INRIA - EDF R&D
// Author(s): Elsa Real
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


// C prototype for using FastJX from Polyphemus.
//
// For more information about FastJX, visit the Polyphemus web site:
//  http://www.ess.uci.edu/group/prather/scholar_software/fast-jx/7-1


#ifdef POLYPHEMUS_FASTJX

#ifdef POLYPHEMUS_SINGLE_UNDERSCORE
#undef POLYPHEMUS_DOUBLE_UNDERSCORE
#elif defined(__GNUG__) && __GNUG__ < 4 && !defined(__INTEL_COMPILER)
#define POLYPHEMUS_DOUBLE_UNDERSCORE
#endif

#ifdef POLYPHEMUS_DOUBLE_UNDERSCORE
#define _fastjx fastjx__
#else
#define _fastjx fastjx_
#endif

extern "C"
{
  void _fastjx(int* NxStep, int* NyStep, int* NzStep, float* x_min,
               float* delta_x, float* y_min, float* delta_y, float* GridZInterf,
               float* loudOpticalDepth, float* IceOpticalDepth, float* SurfacePressure,
               float* AerosolOpticalDepth, float* SingleScatteringAlbedo,
               float* ExtinctionEfficiencyFacor, float* PhaseFunction,
               float* PhotolysisRate, char* PhotolysisNameList,
               int*  nphot, int* year, int* iday, float* hour,
               const char* DirectoryParameter, int* DirectoryParameter_len);
}

#else

// A dummy implementation to generate an error when FastJX is required but not installed.
#define _fastjx(NxStep, NyStep, NzStep, x_min,				\
                delta_x, y_min, delta_y, GridZInterf,			\
                loudOpticalDepth, IceOpticalDepth, SurfacePressure,	\
                AerosolOpticalDepth, SingleScatteringAlbedo,		\
                ExtinctionEfficiencyFacor, PhaseFunction,		\
                PhotolysisRate, PhotolysisNameList,			\
                nphot, year, iday, hour,				\
                DirectoryParameter, DirectoryParameter_len)		\
  throw "The FastJX library is required but not installed or disabled.\n" \
  "\n"									\
  "The FastJX library have to be installed in the FastJX directory "	\
  "\"include/fastJX\". It can be installed by running \"scons\" in the " \
  "FastJX directory.\n"							\
  "\n"									\
  "Once the FastJX library is installed, just rebuild your application " \
  "with \"scons\".\n"							\
  "If \"scons\" is not used, please define manually the macro symbol "	\
  "\"POLYPHEMUS_FASTJX\" before rebuilding your application to get rid " \
  "of this message."
#endif
