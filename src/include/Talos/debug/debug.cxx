// Copyright (C) 2015, ENPC
//    Author(s): Sylvain Doré
//
// This file is part of the air quality modeling system Polyphemus.
//
// Polyphemus is developed in the INRIA project-team CLIME and in
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
#ifndef TALOS_FILE_DEBUG_CXX

#include "exception.cxx"

#ifdef TALOS_DEBUG
#include "helpers.c"
#include "backtrace.c"
#include "signal.c"
#include "exception.cxx"
#include "terminate.cxx"

bool has_debug_signal = init_debug_signals();

// init_backtrace is only availalbe under gcc-4.8 or later
// And has_backtrace is not used (2017/09/21 Youngseob Kim)

//bool has_backtrace = (init_backtrace() == 0); 


namespace Talos
{
  terminate_handler std_terminate = std::set_terminate(debug_terminate);
}

#endif

#define TALOS_FILE_DEBUG_CXX
#endif
