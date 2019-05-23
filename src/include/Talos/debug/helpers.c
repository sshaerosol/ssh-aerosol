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

#include "helpers.h"

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <sys/param.h>


#if defined(__linux__)
#define PROC_SELF "/proc/self/exe"
#elif defined(__sun)
#define PROC_SELF "/proc/self/path/a.out"
#elif defined(BSD)
#define PROC_SELF "/proc/curproc/file"
#endif

#ifdef PROC_SELF
ssize_t get_executable_path(char *path, size_t path_size)
{
  if (path_size <= 1)
    return 0;
  ssize_t written = readlink(PROC_SELF, path, path_size - 1);
  path[written + 1] = '\0';
  return written;
}
#else
ssize_t get_executable_path(char *path, size_t path_size)
{
  // Dummy implementation.
  return 0;
}
#endif

int write_message(const char *message)
{
  ssize_t to_write = strlen(message);
  ssize_t written = write(STDERR_FILENO, message, to_write);
  if (written != to_write)
    return -1;
  return 0;
}


/**
 * Ansi C "itoa" based on Kernighan & Ritchie's "Ansi C"
 * with slight modification to optimize for specific architecture:
 */

static void strreverse(char* begin, char* end)
{
  char aux;
  while (end > begin)
    {
      aux = *end;
      *end-- = *begin;
      *begin++ = aux;
    }
}


void kr_itoa(int value, char* str, int base)
{
  const char* num = "0123456789abcdefghijklmnopqrstuvwxyz";
  char* wstr = str;
  int sign;
  div_t res;

  // Validate base
  if (base < 2 || base > 35)
    {
      *wstr = '\0';
      return;
    }

  // Take care of sign
  if ((sign = value) < 0) value = -value;

  // Conversion. Number is reversed.
  do
    {
      res = div(value, base);
      *wstr++ = num[res.rem];
    }
  while ((value = res.quot));
  if (sign < 0) *wstr++ = '-';
  *wstr = '\0';

  // Reverse string
  strreverse(str, wstr - 1);
}
