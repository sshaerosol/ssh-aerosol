// Copyright (C) 2004-2008, INRIA
// Author(s): Vivien Mallet
//
// This file is part of Talos library, which provides miscellaneous tools to
// make up for C++ lacks and to ease C++ programming.
//
// Talos is free software; you can redistribute it and/or modify it under the
// terms of the GNU Lesser General Public License as published by the Free
// Software Foundation; either version 2.1 of the License, or (at your option)
// any later version.
//
// Talos is distributed in the hope that it will be useful, but WITHOUT ANY
// WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
// FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License for
// more details.
//
// You should have received a copy of the GNU Lesser General Public License
// along with Talos. If not, see http://www.gnu.org/licenses/.
//
// For more information, visit the Talos home page:
//     http://vivienmallet.net/lib/talos/


#ifndef TALOS_FILE_TALOSHEADER_HXX

namespace Talos
{
}

#include <iostream>

//! To display a message... call Hermes!
#ifndef ERR
#define ERR(x) std::cout << "Hermes - " #x << std::endl
#endif
//! To display a variable (with its name).
#ifndef DISP
#define DISP(x) std::cout << #x ": " << (x) << std::endl
#endif

#include "String.hxx"
#include "Date.hxx"
#include "Files.hxx"

#include "debug/debug.hxx"

#define TALOS_FILE_TALOSHEADER_HXX
#endif
