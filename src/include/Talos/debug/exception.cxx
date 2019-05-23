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
#ifndef TALOS_FILE_EXCEPTION_CXX

using namespace std;

#include "exception.hxx"

#include <exception>

namespace Talos
{

  //! Re-throws the current exception and re-catch it to display
  //! an appropriate error message.
  /** It is useful with compilers whose terminate handlers
   * accept re-throw of the uncaught exception.
   */
  void
  write_exception()
  {
    // TODO:
    // - improve message formating
    // - use typeid (ISO 2003) and cxxabi.h (gcc specific)
    //   to get the dynamic pretty typename.
    try
      {
        throw;
      }
    catch (std::exception& Err)
      {
        cerr << "[ERROR] 'std::exception' caught:" << endl;
        cerr << "    \"" << Err.what() << "\"" << endl;
      }
    catch (std::string& str)
      {
        cerr << "[ERROR] exception of type 'std::string' caught: ";
        cerr << "\"" << str << "\"" << endl;
      }
    catch (const char* str)
      {
        cerr << "[ERROR] exception of type 'char*' caught: ";
        cerr << "\"" << str << "\"" << endl;
      }
    catch (...)
      {
        cerr << "[ERROR] Unknown exception caught..." << endl;
      }
  }

}

#define TALOS_FILE_EXCEPTION_CXX
#endif
