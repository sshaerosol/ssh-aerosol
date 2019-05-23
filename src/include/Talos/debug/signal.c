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

#include <fenv.h>
#include <signal.h>
#ifndef __cplusplus
#  include <stdbool.h>
#endif
#include <stdlib.h>
#include <sys/resource.h>
#include <unistd.h>

#include "helpers.h"

typedef void (*sighandler_t)(int);


///////////////////////////////////
// The purpose is:
// - enabling core dumps
// - enabling floating point exception
// - better error report in case of failure
///////////////////////////////////


// POSIX defines a maximum of 64 signal.
#ifdef _NSIG
#define MAX_SIGNAL_NUMBER _NSIG - 1
#else
#define MAX_SIGNAL_NUMBER 64
#endif


typedef void (signal_handler_t)(int, siginfo_t *, void *);

struct sigaction action_list[MAX_SIGNAL_NUMBER];

void register_signal_handler(int signal_id, signal_handler_t* handler)
{
  struct sigaction* action = &action_list[signal_id];

  sigfillset(&action->sa_mask);
  sigdelset(&action->sa_mask, signal_id);
  action->sa_flags = SA_SIGINFO;
  action->sa_flags |= SA_RESETHAND;

  if (handler)
    action->sa_sigaction = handler;
  else // Default handler is put back.
    action->sa_handler = SIG_DFL;

  sigaction(signal_id, action, 0);
}


static const char* fpe_description(siginfo_t* sinfo)
{
  const char* no_description = "Bad operation on float or integer.";
  if (sinfo == 0)
    return no_description;

  switch (sinfo->si_code)
    {
#ifdef FPE_NOOP   /* BUG: MacOSX uses this instead of INTDIV */
    case FPE_NOOP:
#endif
    case FPE_INTDIV:
      return "Integer divided by zero.";
    case FPE_INTOVF:
      return "Integer overflow.";
    case FPE_FLTDIV:
      return "Float divided by zero.";
    case FPE_FLTOVF:
      return "Float overflow.";
    case FPE_FLTUND:
      return "Float underflow.";
    case FPE_FLTRES:
      return "Float inexact result.";
    case FPE_FLTINV:
      return "Invalid operation on float.";
    case FPE_FLTSUB:
      return "Subscript out of range.";
    default:
      return no_description;
    }
}


static void write_signal_info(int sig, siginfo_t* sinfo)
{
  char signal_number[1024];
  write_message("Received signal ");
  kr_itoa(sig, signal_number, 10);
  write_message(signal_number);
  write_message(" ");

  const char* message = 0;
  switch (sig)
    {
    case SIGINT:
      write_message("(SIGINT) - Process interruption.\n"
                    "  That can occur when typing Ctrl-C (Unix/OSX) or "
                    "Ctrl-B (Windows) on a keyboard.");
      break;
    case SIGQUIT:
      write_message("(SIGQUIT) - Process quit.\n"
                    "  That can occur when typing Ctrl-\\ on a keyboard.");
      break;
    case SIGILL:
      write_message("(SIGILL) - Illegal instruction:\n"
                    "  A privileged or invalid instruction was encountered.");
      break;
    case SIGABRT:
      write_message("(SIGABRT) - The process execution is aborted.");
      break;
    case SIGBUS:
      write_message("(SIGBUS) - BUS error:\n"
                    "  the process requested an unreachable memory address\n"
                    "  that was was pointing out of memory or was unaligned.\n"
                    "  (can occur like segmentation fault on common memory"
                    " bugs...)");
      break;
    case SIGFPE:
      write_message("(SIGFPE) - Bad numerical operation:\n  ");
      write_message(fpe_description(sinfo));
      break;
    case SIGSEGV:
      write_message("(SIGSEGV) - Segmentation fault:\n"
                    "  a forbidden memory address was requested by the process\n"
                    "  (likely a memory bug).");
      break;
    default:
      // Nothing.
      break;
    }

  if (message != 0)
    {
      write_message(" - ");
      write_message(message);
    }
  write_message("\n");
}


void default_abort()
{
  // Puts back in place the default handler:
  register_signal_handler(SIGABRT, 0);
  // Triggers SIGABRT:
  ::abort();
}


static void error_signal_handler(int sig, siginfo_t* sinfo, void* ucontext)
{
  if (sig == SIGINT || sig == SIGQUIT)
    write_message("\n\n[WARNING] ");
  else
    write_message("\n\n[ERROR] ");
  write_signal_info(sig, sinfo);

#ifndef TALOS_NO_COREDUMP
  // Displaying the backtrace here does sometimes block the process...
  // So it is preferable to abort to get a core dump:
  write_message("\nHint: A core dump should be generated, you can open it with"
                " your debugger to get the backtrace and other informations.\n"
                "\n");
#endif

  default_abort();
  // Would have triggered the default signal hangler, which do not necessarily
  // result in a core dump:
  // raise(sig);
}


bool init_debug_signals()
{
  //=== Enables core dump.
#ifndef TALOS_NO_COREDUMP
  static struct rlimit core_limits;
  core_limits.rlim_cur = core_limits.rlim_max = RLIM_INFINITY;
  // Keep in mind that GCC address sanitizer is disabling
  // core dumps (because they would be enormous), and they can't be reenabled
  // even with this:
  setrlimit(RLIMIT_CORE, &core_limits);
#endif


  //=== Enabling floating point exception, aka FPE (or FE).

// With FPE signals, there is an abort since recovering from such
// a signal leads to an undefined behavior according to POSIX.
// Note that the GCC sanitizer is able to display an error without aborting
// the processus, so one might want to disable these FPE with the corresponding
// macros.
#ifndef TALOS_NO_FE_DIVBYZERO
  feenableexcept(FE_DIVBYZERO);
#endif

// FE_OVERFLOW is triggered when a float result was too large to be
// representable.
#ifndef TALOS_NO_FE_OVERFLOW
  feenableexcept(FE_OVERFLOW);
#endif

// FE_INVALID is triggered on float domain error.
// FE_INVALID does not work very well with the GCC sanitizer
// 'float-cast-overflow'. This signal is however useful, and can be preferred
// over 'float-cast-overflow' sanitizer for default debugging.
#ifndef TALOS_NO_FE_INVALID
  feenableexcept(FE_INVALID);
#endif

  // (FE_INEXACT and FE_UNDERFLOW should not worth enabling.)

  // TODO:
  // - put the comment on FPE limitations
  // - use the more portable fesetexcept

  // Installing debugging signal handlers.
  register_signal_handler(SIGINT,  error_signal_handler);
  register_signal_handler(SIGQUIT, error_signal_handler);
  register_signal_handler(SIGILL,  error_signal_handler);
  register_signal_handler(SIGBUS,  error_signal_handler);
  register_signal_handler(SIGFPE,  error_signal_handler);
  register_signal_handler(SIGSEGV, error_signal_handler);
  register_signal_handler(SIGABRT, error_signal_handler);

  return true;
}
