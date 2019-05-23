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

// Displaying a back trace in a signal safe manner (which is a stronger
// assumption than multi-thread safe) is error-prone and can lead
// to dead-locks. The Gcc Backtrace library complies with this
// constraint provided that the BACKTRACE_USES_MALLOC macro is defined to 0
// (meaning mmap is used). Note also that multi-threading support is not
// necessary since a mutex is used to preserve a good ordering in the
// stack logging.
// See:
// https://github.com/gcc-mirror/gcc/blob/master/libbacktrace/
// Any attempt to port this functionality to other compiler/plate-form
// should be carried out with great care.

#include "backtrace.h"

#if __GNUC__ > 4 || ( __GNUC__ == 4 && __GNUC_MINOR__ >= 8 )
#  include <backtrace-supported.h>
#  if BACKTRACE_SUPPORTED && !BACKTRACE_USES_MALLOC
#    define HAVE_LIBBACKTRACE
#  endif
#endif


#ifdef HAVE_LIBBACKTRACE

#include <backtrace.h>
#include <stdlib.h>

#ifdef __cplusplus
#  include <cxxabi.h>
#endif

#ifndef WIN32
#include <pthread.h>
static pthread_mutex_t backtrace_lock = PTHREAD_MUTEX_INITIALIZER;
#endif

#include "helpers.h"

struct backtrace_state* state = 0;


static void error_callback(void *data, const char *msg, int errnum)
{
  static const char *last_msg = 0;
  if (msg == last_msg)
    return;
  last_msg = msg;

  char error_number[1024];
  kr_itoa(errnum, error_number, 10);
  write_message("[ERROR] Error while unwinding the call stack: ");
  write_message(msg);
  write_message(" (errnum=");
  write_message(error_number);
  write_message(")\n");
}


int init_backtrace()
{
  char executable_path[4096]; // Large enough for real life use.
  size_t executable_path_size = sizeof(executable_path) / sizeof(char);
  const char* filename = 0;

  if (get_executable_path(executable_path, executable_path_size) > 0)
    filename = executable_path;
  state = backtrace_create_state(filename,
                                 BACKTRACE_SUPPORTS_THREADS,
                                 error_callback, NULL);
  if (state == 0)
    return -1;
  return 0;
}


static int frame_callback(void *vdata, uintptr_t pc, const char *filename,
                          int lineno, const char *function)
{
  char line_number[64];
  char demangle_buffer[1024];

  if (filename == 0 && function == 0)
    return 0;

#ifdef __cplusplus
  if (function && function[0] == '_')
    {
      // abi::__cxa_demangle is requiring a malloc allocated buffer,
      // which is impossible since it has to be used in a signal handler.
      // However, malloc is required only for __cxa_demangle to call realloc
      // when needed, which won't happen with real life code and the following
      // buffer:
      int status;
      size_t demangle_buffer_size = sizeof(demangle_buffer) / sizeof(char);
      const char* demangled_function = abi::__cxa_demangle(function,
                                                           demangle_buffer,
                                                           &demangle_buffer_size,
                                                           &status);
      if (demangled_function)
        function = demangled_function;
    }
#endif

  if (function == 0)
    function = "<unknown function>";

  if (filename == 0)
    filename = "<unknown file>";

  kr_itoa(lineno, line_number, 10);

  write_message("  |  ");
  write_message(function);
  write_message("()\n  |      in ");
  write_message(filename);
  write_message(":");
  write_message(line_number);
  write_message("\n");

  return 0;
}


void write_backtrace(int skip)
{
#ifndef WIN32
  pthread_mutex_lock(&backtrace_lock);
#endif
  {
    if (state == 0)
      {
        write_message("\n[ERROR] Cannot display Backtrace!\n");
        return;
      }

    write_message("\nBacktrace:\n");
    backtrace_full(state, skip, frame_callback, error_callback, 0);
    write_message("\n");
  }
#ifndef WIN32
  pthread_mutex_unlock(&backtrace_lock);
#endif
}


#else


void write_backtrace(int)
{
  write_message(
                "(No backtrace display because the compiler is not supported.)\n");
}


#endif // HAVE_LIBBACKTRACE
