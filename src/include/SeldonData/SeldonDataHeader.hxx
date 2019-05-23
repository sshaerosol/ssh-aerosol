// Copyright (C) 2003-2011, INRIA
// Author(s): Vivien Mallet
//
// This file is part of SeldonData library, used for data processing.
//
// SeldonData is free software; you can redistribute it and/or modify it under
// the terms of the GNU General Public License as published by the Free
// Software Foundation; either version 2 of the License, or (at your option)
// any later version.
//
// SeldonData is distributed in the hope that it will be useful, but WITHOUT
// ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
// FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for
// more details.
//
// For more information, visit the SeldonData home page:
//      http://vivienmallet.net/lib/seldondata/

#ifndef FILE_SELDONDATA_SELDONDATAHEADER_HXX

// Blitz++.
#include <blitz/array.h>
namespace blitz {}
using namespace blitz;

// Talos.
#include <TalosHeader.hxx>
using namespace Talos;

#include <string>
#include <sstream>
#include <exception>
#include <stdexcept>
#include <algorithm>


//! SeldonData namespace.
namespace SeldonData
{
  using Talos::to_num;
  using Talos::to_str;
}


//////////////////
// DEBUG LEVELS //
//////////////////

#ifdef SELDONDATA_DEBUG_LEVEL_4
#ifndef SELDONDATA_DEBUG_LEVEL_3
#define SELDONDATA_DEBUG_LEVEL_3
#endif
#endif

#ifdef SELDONDATA_DEBUG_LEVEL_3
#ifndef SELDONDATA_DEBUG_CHECK_INDICES
#define SELDONDATA_DEBUG_CHECK_INDICES
#endif
#ifndef SELDONDATA_DEBUG_LEVEL_2
#define SELDONDATA_DEBUG_LEVEL_2
#endif
#endif

#ifdef SELDONDATA_DEBUG_LEVEL_2
#ifndef SELDONDATA_DEBUG_CHECK_DIMENSIONS
#define SELDONDATA_DEBUG_CHECK_DIMENSIONS
#endif
#ifndef SELDONDATA_DEBUG_LEVEL_1
#define SELDONDATA_DEBUG_LEVEL_1
#endif
#endif

#ifdef SELDONDATA_DEBUG_LEVEL_1
#ifndef SELDONDATA_DEBUG_CHECK_MEMORY
#define SELDONDATA_DEBUG_CHECK_MEMORY
#endif
#ifndef SELDONDATA_DEBUG_CHECK_IO
#define SELDONDATA_DEBUG_CHECK_IO
#endif
#ifndef SELDONDATA_DEBUG_UNINITIALIZED_IS_NAN
#define SELDONDATA_DEBUG_UNINITIALIZED_IS_NAN
#endif
#ifndef SELDONDATA_DEBUG_LEVEL_0
#define SELDONDATA_DEBUG_LEVEL_0
#endif
#endif


#include "Errors.hxx"
#include "Function.hxx"
#include "Grid.hxx"
#include "Data.hxx"
#include "Format.hxx"

#define FILE_SELDONDATA_SELDONDATAHEADER_HXX
#endif
