// Copyright (C) 2003-2011, ENPC - INRIA - EDF R&D
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


#ifndef ATMODATA_FILE_ATMODATAHEADER_HXX

#include "SeldonDataHeader.hxx"
using namespace SeldonData;

//! AtmoData namespace.
namespace AtmoData
{
  using Talos::to_num;
  using Talos::to_str;
  using SeldonData::IOError;
  using SeldonData::WrongDim;
}

#include "Common.hxx"

#include "Meteorology.hxx"
#include "Kz.hxx"
#include "TimeDiagnosis.hxx"
#include "Photolysis.hxx"
#include "Emissions.hxx"
#include "Megan.hxx"
#include "Deposition.hxx"
#include "CoordTransform.hxx"
#include "Transform.hxx"
#include "Polair.hxx"
#include "Aerosol.hxx"
#include "Optics.hxx"
#include "InterpolationChimere.hxx"

#include "Format.hxx"

#include "Errors.hxx"

#define ATMODATA_FILE_ATMODATAHEADER_HXX
#endif
