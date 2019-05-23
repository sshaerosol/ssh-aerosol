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

#include <exception>
#include <string>
using namespace std;

#include "exception.hxx"

#include "helpers.h"
#include "signal.h"
#include "backtrace.h"

#include <cstdlib>
#include <iostream>
#include <typeinfo>

#ifndef WIN32
#include <pthread.h>
static pthread_mutex_t terminate_lock = PTHREAD_MUTEX_INITIALIZER;
#endif

namespace Talos
{
  void debug_terminate()
  {
#ifndef WIN32
    pthread_mutex_lock(&terminate_lock);
#endif
    // Other compilers that accept re-throw in the terminate handler
    // should be added at will.
#ifdef __GNUC__
    write_exception();
#else
    write_message(
                  "(No exception display because the compiler is not supported.)\n");
#endif

    write_backtrace(2);
    std::abort();
#ifndef WIN32
    pthread_mutex_unlock(&terminate_lock);
#endif
  }

}
