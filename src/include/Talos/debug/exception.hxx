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
#ifndef TALOS_FILE_EXCEPTION_HXX

namespace Talos
{

  //! Displays an error message corresponding to the caught exception.
  void write_exception();

}

#ifndef TALOS_DISABLE_TRY_END
  //! Convenient macros to catch exceptions.
  /*!
    Use TRY and END to catch exceptions:

    [...]
    TRY;

    Code that could throw an exception.

    END;

    [...]

    If an exception is caught, explanations are displayed to identify
    the problem (name of the involved function and comments).

    They are disabled with the 'TALOS_DEBUG' macros to allow for debugging
    exceptions.
  */
#if defined(TRY) or defined(END)
#  warning "'TRY'/'END' macros already defined, " \
  " the Talos library cannot define its own version."
#elif defined(TALOS_DEBUG)
  // A caught exception is considered as normal by the debugger. Only an uncaught
  // exception is considered as an error and make the debugger to stop on them.
  // This is why 'TRY'/'END' are disabled when 'TALOS_DEBUG' is defined:
#  define TRY {}
#  define END {}
#else
#  define TRY try {
#  define END }                 \
 catch (...)                    \
   {                            \
     Talos::write_exception();  \
     return 1;                  \
   }
#endif
#endif // ifndef TALOS_DISABLE_TRY_END

#define TALOS_FILE_EXCEPTION_HXX
#endif
