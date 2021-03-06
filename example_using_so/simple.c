//
// How to compile and run this:
// ./compile && LD_LIBRARY_PATH=../src/:$LD_LIBRARY_PATH ./c_simple
//
/*

  !!-----------------------------------------------------------------------
  !!     Copyright (C) 2019 CEREA (ENPC) - INERIS
  !!     SSH-aerosol is distributed under the GNU General Public License v3
  !!-----------------------------------------------------------------------

  This file is also partly from Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2019 EDF S.A.

  This program is free software; you can redistribute it and/or modify it under
  the terms of the GNU General Public License as published by the Free Software
  Foundation; either version 2 of the License, or (at your option) any later
  version.

  This program is distributed in the hope that it will be useful, but WITHOUT
  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
  FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
  details.

  You should have received a copy of the GNU General Public License along with
  this program; if not, write to the Free Software Foundation, Inc., 51 Franklin
  Street, Fifth Floor, Boston, MA 02110-1301, USA.

*/

/*----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <stdbool.h>
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <math.h>
#include <dlfcn.h>

/*----------------------------------------------------------------------------
 * Get a shared library function pointer
 *
 * parameters:
 *   handle           <-- pointer to shared library (result of dlopen)
 *   name             <-- name of function symbol in library
 *   errors_are_fatal <-- abort if true, silently ignore if false
 *
 * returns:
 *   pointer to the given function in the shared library
 *----------------------------------------------------------------------------*/

static void *
_get_dl_function_pointer(void        *handle,
                         const char  *lib_path,
                         const char  *name,
                         bool         errors_are_fatal)
{

  void  *retval = NULL;
  char  *error = NULL;
  char  *name_ = NULL;

  dlerror();    /* Clear any existing error */

  retval = dlsym(handle, name);
  error = dlerror();

  if (error != NULL) { /* Try different symbol names */
    dlerror();    /* Clear any existing error */
    int _size_ = strlen(name) + strlen("_");
    name_ = (char*)malloc((_size_ + 1)*sizeof(char));
    strcpy(name_, name);
    strcat(name_, "_");
    retval = dlsym(handle, name_);
    error = dlerror();

    if (error != NULL) printf("Error on dlsym : %s\n",error);

    free(name_);
  }

  if (error != NULL && errors_are_fatal) {
    printf("Error, abort program\n");
    exit(1);
  }

  return retval;
}

void main(int argc, char** argv)
{

  // Number of gas species
  int ns;
  // Number of aerosol layers
  int nlayer;
  // Number of aerosols
  int nsize;
  // Time step
  double dtref;
  // Time
  double tt;
  // This is used to handle errors
  char *error = NULL;
  // This array will be used to store gas species concentrations
  double *gas = NULL;
  // This array will be used to store aerosol species concentrations
  double *aero = NULL;
  // This array will be used to store aerosol numbers
  double *aeronum = NULL;

  {

    printf("Start program.\n");

    void *handle;
    const char lib_path[] = "libssh-aerosol.so";

    dlerror();    /* Clear any existing error */

    // Try to load the shared library
    printf("dlopen\n");
    handle = dlopen(lib_path, RTLD_LAZY);
    error = dlerror();
    if (error != NULL) {
      printf("Error on dlopen : %s\n",error);
      return;
    }

    // Initialize SSH-aerosol with namelist.ssh
    {
      const char namelist_ssh[41] = "namelist.ssh";
      typedef void (*api_sshaerosol_initialize_t)(char*, int*, int*, int*);
      printf("get_function_pointer\n");
      api_sshaerosol_initialize_t api_sshaerosol_initialize =
      (api_sshaerosol_initialize_t) _get_dl_function_pointer(handle,
                                                            lib_path,
                                                           "api_sshaerosol_simple_initialize",
                                                            true);
      api_sshaerosol_initialize(&namelist_ssh, &ns, &nlayer, &nsize);
    }

    // Allocate memory
    gas = (double*) malloc(ns*sizeof(double));
    aero = (double*) malloc(nlayer*nsize*sizeof(double));
    aeronum = (double*) malloc(nsize*sizeof(double));
    if (gas == NULL || aero == NULL || aeronum == NULL) {
      printf("Memory allocation failed");
      exit(1);
    }

    dtref = 60;

    // Perform 60 time steps
    for (int idt = 0; idt < 60; idt++){

      // Call gaseous chemistry
      {
        typedef void (*api_sshaerosol_t)(double*);
        api_sshaerosol_t fct =
          (api_sshaerosol_t) _get_dl_function_pointer(handle,
                                                      lib_path,
                                                     "api_sshaerosol_simple_dt",
                                                      true);
        fct(&dtref);
      }

    }

    // Finalize SSH
    {
      typedef void (*api_sshaerosol_t)();
      api_sshaerosol_t fct =
        (api_sshaerosol_t) _get_dl_function_pointer(handle,
                                                    lib_path,
                                                   "api_sshaerosol_finalize",
                                                    true);
      fct();
    }

    // Close the shared library
    printf("dclose\n");
    dlclose(handle);
    error = dlerror();
    if (error != NULL) printf("Error on dlclose : %s\n",error);

  }

  // Free memory
  free(gas);
  free(aero);
  free(aeronum);

}

/*----------------------------------------------------------------------------*/
