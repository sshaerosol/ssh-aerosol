//
// How to compile and run this:
// ./compile && LD_LIBRARY_PATH=../src/:$LD_LIBRARY_PATH ./c_exemple
//
// Note:
//   - In this exemple, all the prototypes of the functions are
//     used and discarded on the fly. The prototypes should rather
//     be defined once with clearer names so they can be reused.
//   - The complete API can be found in ../src/include/Module/ModuleAPI.f90
//   - The API can be extended, feel free to contact us.
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
  // Number of reactions
  int nr;
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

    // Declare SSH-aerosol as not running standalone
    // This removes most of the output to stdout
    {
      typedef void (*api_set_sshaerosol_t)(bool*);
      api_set_sshaerosol_t fct =
        (api_set_sshaerosol_t) _get_dl_function_pointer(handle,
                                                       lib_path,
                                                      "api_sshaerosol_set_standalone",
                                                       true);
      bool flag = false;
      fct(&flag);
    }

    // Activate the logger
    // Only one rank should activate the logger in parallel
    {
      typedef void (*api_set_sshaerosol_t)(bool*);
      api_set_sshaerosol_t fct =
        (api_set_sshaerosol_t) _get_dl_function_pointer(handle,
                                                       lib_path,
                                                      "api_sshaerosol_set_logger",
                                                       true);
      bool flag = true;
      fct(&flag);
    }

    // Check that SSH-aerosol is not running standalone
    {
      typedef bool (*api_get_sshaerosol_t)();
      api_get_sshaerosol_t fct =
        (api_get_sshaerosol_t) _get_dl_function_pointer(handle,
                                                       lib_path,
                                                      "api_sshaerosol_get_standalone",
                                                       true);
      bool flag = fct();
      printf("Standalone is %s\n", flag ? "true" : "false");
    }

    // Check we have correctly activated the logger
    {
      typedef bool (*api_get_sshaerosol_t)();
      api_get_sshaerosol_t fct =
        (api_get_sshaerosol_t) _get_dl_function_pointer(handle,
                                                       lib_path,
                                                      "api_sshaerosol_get_logger",
                                                       true);
      bool flag = fct();
      printf("Logger is %s\n", flag ? "true" : "false");
    }

    // Initialize SSH-aerosol with namelist.ssh
    {
      const char namelist_ssh[41] = "namelist.ssh";
      typedef void (*api_sshaerosol_initialize_t)(char*);
      printf("get_function_pointer\n");
      api_sshaerosol_initialize_t api_sshaerosol_initialize =
      (api_sshaerosol_initialize_t) _get_dl_function_pointer(handle,
                                                            lib_path,
                                                           "api_sshaerosol_initialize",
                                                            true);
      api_sshaerosol_initialize(&namelist_ssh);
    }

    // Get the number of gas species
    {
      typedef int (*api_get_sshaerosol_t)();
      api_get_sshaerosol_t fct =
        (api_get_sshaerosol_t) _get_dl_function_pointer(handle,
                                                       lib_path,
                                                      "api_sshaerosol_get_ngas",
                                                       true);
      ns = fct();
      printf("N_gas : %d\n", ns);
    }

    // Get the number of aerosol layers
    {
      typedef int (*api_get_sshaerosol_t)();
      api_get_sshaerosol_t fct =
        (api_get_sshaerosol_t) _get_dl_function_pointer(handle,
                                                       lib_path,
                                                      "api_sshaerosol_get_nlayer",
                                                       true);
      nlayer = fct();
      printf("N_aerosol_layers : %d\n", nlayer);
    }

    // Get the number of aerosols
    {
      typedef int (*api_get_sshaerosol_t)();
      api_get_sshaerosol_t fct =
        (api_get_sshaerosol_t) _get_dl_function_pointer(handle,
                                                       lib_path,
                                                      "api_sshaerosol_get_nsize",
                                                       true);
      nsize = fct();
      printf("N_size : %d\n", nsize);
    }

    // Allocate memory
    gas = (double*) malloc(ns*sizeof(double));
    aero = (double*) malloc(nlayer*nsize*sizeof(double));
    aeronum = (double*) malloc(nsize*sizeof(double));
    if (gas == NULL || aero == NULL || aeronum == NULL) {
      printf("Memory allocation failed");
      exit(1);
    }

    // Get the name of a aerosol
    {
      typedef void (*api_get_sshaerosol_t)(int*, char*);
      api_get_sshaerosol_t fct =
        (api_get_sshaerosol_t) _get_dl_function_pointer(handle,
                                                        lib_path,
                                                       "api_sshaerosol_get_aero_name",
                                                        true);
      int iaero = 4;
      char name[81];
      fct(&iaero, &name);
      printf("iaero : %d, name : %s\n", iaero, name);
    }

    // Read the initial time
    {
      typedef double (*api_get_sshaerosol_dt_t)();
      api_get_sshaerosol_dt_t fct =
        (api_get_sshaerosol_dt_t) _get_dl_function_pointer(handle,
                                                          lib_path,
                                                         "api_sshaerosol_get_initial_t",
                                                          true);
      tt = fct();
      printf("Initial time : %f\n", tt);
    }

    // Read the time step used by SSH-aerosol
    {
      typedef double (*api_get_sshaerosol_dt_t)();
      api_get_sshaerosol_dt_t fct =
        (api_get_sshaerosol_dt_t) _get_dl_function_pointer(handle,
                                                          lib_path,
                                                         "api_sshaerosol_get_dt",
                                                          true);
      dtref = fct();
      printf("Delta t : %f\n", dtref);
    }

    // Set the time step of SSH-aerosol
    {
      typedef bool (*api_set_sshaerosol_dt_t)(double*);
      api_set_sshaerosol_dt_t fct =
        (api_set_sshaerosol_dt_t) _get_dl_function_pointer(handle,
                                                          lib_path,
                                                         "api_sshaerosol_set_dt",
                                                          true);
      double dt = 0.123456;
      bool check = fct(&dt);
      if (!check) printf("Warning: Time step was above DTAEROMIN");
    }

    // Read again the time step used by SSH-aerosol
    {
      typedef double (*api_get_sshaerosol_dt_t)();
      api_get_sshaerosol_dt_t fct =
        (api_get_sshaerosol_dt_t) _get_dl_function_pointer(handle,
                                                          lib_path,
                                                         "api_sshaerosol_get_dt",
                                                          true);
      double dt = fct();
      printf("Delta t : %f\n", dt);
    }

    // Restore the time step
    {
      typedef bool (*api_set_sshaerosol_dt_t)(double*);
      api_set_sshaerosol_dt_t fct =
        (api_set_sshaerosol_dt_t) _get_dl_function_pointer(handle,
                                                          lib_path,
                                                         "api_sshaerosol_set_dt",
                                                          true);
      bool check = fct(&dtref);
    }

    // This bypasses the API with direct access to the symbols of the shared library
    // Name of the symbols are compiler-dependent
    // In this exemple, the time step is read.
    // One can check the symbols in the library with the command nm (Linux):
    //   $ nm -D libssh-aerosol.so
    {
      void *nn = dlsym(handle,"__ainitialization_MOD_delta_t");
      error = dlerror();
      if (error != NULL) printf("Error on dlsym : %s\n",error);
      printf("Delta t : %f\n",*(double *)nn);
    }

    // Read the gas species concentrations in SSH-aerosol
    {
      typedef void* (*api_get_sshaerosol_t)(double*);
      api_get_sshaerosol_t fct =
        (api_get_sshaerosol_t) _get_dl_function_pointer(handle,
                                                       lib_path,
                                                      "api_sshaerosol_get_gas_",
                                                       true);
      fct(gas);
      printf("Get gaseous concentrations:");
      for (int i = 0; i < ns; i++) printf(" %f", gas[i]);
      printf("\n");
    }

    // Set the gas species concentrations in SSH-aerosol
    {
      typedef void (*api_set_sshaerosol_t)(double*);
      api_set_sshaerosol_t fct =
        (api_set_sshaerosol_t) _get_dl_function_pointer(handle,
                                                       lib_path,
                                                      "api_sshaerosol_set_gas_",
                                                       true);
      double data[ns];
      for (int i = 0; i < ns; i++) data[i] = i+1;
      fct(&data);
      printf("Set gaseous concentrations:");
      for (int i = 0; i < ns; i++) printf(" %f", data[i]);
      printf("\n");
    }

    // Read again the gas species concentrations
    {
      typedef void* (*api_get_sshaerosol_t)(double*);
      api_get_sshaerosol_t fct =
        (api_get_sshaerosol_t) _get_dl_function_pointer(handle,
                                                       lib_path,
                                                      "api_sshaerosol_get_gas_",
                                                       true);
      double data[ns];
      fct(&data);
      printf("Get gaseous concentrations:");
      for (int i = 0; i < ns; i++) printf(" %f", data[i]);
      printf("\n");
    }

    // Restore the original values
    {
      typedef void (*api_set_sshaerosol_t)(double*);
      api_set_sshaerosol_t fct =
        (api_set_sshaerosol_t) _get_dl_function_pointer(handle,
                                                       lib_path,
                                                      "api_sshaerosol_set_gas",
                                                       true);
      fct(gas);
    }

    // Read the aerosols concentrations
    {
      typedef void* (*api_get_sshaerosol_t)(double*);
      api_get_sshaerosol_t fct =
        (api_get_sshaerosol_t) _get_dl_function_pointer(handle,
                                                       lib_path,
                                                      "api_sshaerosol_get_aero_",
                                                       true);
      fct(aero);
      printf("Get aerosols concentrations:");
      for (int i = 0; i < nlayer*nsize; i++) printf(" %f", aero[i]);
      printf("\n");
    }

    // Set the aerosol concentrations
    {
      typedef void (*api_set_sshaerosol_t)(double*);
      api_set_sshaerosol_t fct =
        (api_set_sshaerosol_t) _get_dl_function_pointer(handle,
                                                       lib_path,
                                                      "api_sshaerosol_set_aero_",
                                                       true);
      double data[nlayer*nsize];
      for (int i = 0; i < nlayer*nsize; i++) data[i] = i+1;
      fct(&data);
      printf("Set aerosols concentrations:");
      for (int i = 0; i < nlayer*nsize; i++) printf(" %f", data[i]);
      printf("\n");
    }

    // Read again the aerosols concentration
    {
      typedef void* (*api_get_sshaerosol_t)(double*);
      api_get_sshaerosol_t fct =
        (api_get_sshaerosol_t) _get_dl_function_pointer(handle,
                                                       lib_path,
                                                      "api_sshaerosol_get_aero_",
                                                       true);
      double data[nlayer*nsize];
      fct(&data);
      printf("Get aerosols concentrations:");
      for (int i = 0; i < nlayer*nsize; i++) printf(" %f", data[i]);
      printf("\n");
    }

    // Restore the aerosols concentrations
    {
      typedef void (*api_set_sshaerosol_t)(double*);
      api_set_sshaerosol_t fct =
        (api_set_sshaerosol_t) _get_dl_function_pointer(handle,
                                                       lib_path,
                                                      "api_sshaerosol_set_aero_",
                                                       true);
      fct(aero);
    }

    // Read the aerosols numbers
    {
      typedef void* (*api_get_sshaerosol_t)(double*);
      api_get_sshaerosol_t fct =
        (api_get_sshaerosol_t) _get_dl_function_pointer(handle,
                                                       lib_path,
                                                      "api_sshaerosol_get_aero_num",
                                                       true);
      fct(aeronum);
      printf("Get aerosols numbers:");
      for (int i = 0; i < nsize; i++) printf(" %f", aeronum[i]);
      printf("\n");
    }

    // Set the aerosols numbers
    {
      typedef void (*api_set_sshaerosol_t)(double*);
      api_set_sshaerosol_t fct =
        (api_set_sshaerosol_t) _get_dl_function_pointer(handle,
                                                       lib_path,
                                                      "api_sshaerosol_set_aero_num",
                                                       true);
      double data[nsize];
      for (int i = 0; i < nsize; i++) data[i] = 0;
      fct(&data);
      printf("Set aerosols numbers:");
      for (int i = 0; i < nsize; i++) printf(" %f", data[i]);
      printf("\n");
    }

    // Read again the aerosols numbers
    {
      typedef void* (*api_get_sshaerosol_t)(double*);
      api_get_sshaerosol_t fct =
        (api_get_sshaerosol_t) _get_dl_function_pointer(handle,
                                                       lib_path,
                                                      "api_sshaerosol_get_aero_num",
                                                       true);
      double data[nsize];
      fct(&data);
      printf("Get aerosols numbers:");
      for (int i = 0; i < nsize; i++) printf(" %f", data[i]);
      printf("\n");
    }

    // Restore the aerosols numbers
    {
      typedef void (*api_set_sshaerosol_t)(double*);
      api_set_sshaerosol_t fct =
        (api_set_sshaerosol_t) _get_dl_function_pointer(handle,
                                                       lib_path,
                                                      "api_sshaerosol_set_aero_num",
                                                       true);
      fct(aeronum);
    }

    // Ask SSH-aerosol to perform some IO
    {
      typedef void (*api_sshaerosol_t)();
      api_sshaerosol_t fct =
        (api_sshaerosol_t) _get_dl_function_pointer(handle,
                                                    lib_path,
                                                   "api_sshaerosol_initoutput",
                                                    true);
      fct();
    }
    {
      typedef void (*api_sshaerosol_t)();
      api_sshaerosol_t fct =
        (api_sshaerosol_t) _get_dl_function_pointer(handle,
                                                    lib_path,
                                                   "api_sshaerosol_report",
                                                    true);
      fct();
    }
    {
      typedef void (*api_sshaerosol_t)();
      api_sshaerosol_t fct =
        (api_sshaerosol_t) _get_dl_function_pointer(handle,
                                                    lib_path,
                                                   "api_sshaerosol_output",
                                                    true);
      fct();
    }

    // Initialize photolysis
    {
      typedef void (*api_sshaerosol_t)();
      api_sshaerosol_t fct =
        (api_sshaerosol_t) _get_dl_function_pointer(handle,
                                                    lib_path,
                                                   "api_sshaerosol_initphoto",
                                                    true);
      fct();
    }

    // Perform a few time steps
    for (int idt = 0; idt < 60; idt++){

      //
      // A realistic external code would set the gas and aerosol concentrations here
      // It is also possible to change T, P, HR, pH, ...
      //

      // Update time in SSH-aerosol
      tt += dtref;
      {
        typedef void (*api_set_sshaerosol_dt_t)(double*);
        api_set_sshaerosol_dt_t fct =
          (api_set_sshaerosol_dt_t) _get_dl_function_pointer(handle,
                                                            lib_path,
                                                           "api_sshaerosol_set_current_t",
                                                            true);
        fct(&tt);
      }

      // Update photolysis
      {
        typedef void (*api_sshaerosol_t)();
        api_sshaerosol_t fct =
          (api_sshaerosol_t) _get_dl_function_pointer(handle,
                                                      lib_path,
                                                     "api_sshaerosol_updatephoto",
                                                      true);
        fct();
      }

      // Emissions
      {
        typedef void (*api_sshaerosol_t)();
        api_sshaerosol_t fct =
          (api_sshaerosol_t) _get_dl_function_pointer(handle,
                                                      lib_path,
                                                     "api_sshaerosol_emission",
                                                      true);
        fct();
      }

      // Call gaseous chemistry
      {
        typedef void (*api_sshaerosol_t)();
        api_sshaerosol_t fct =
          (api_sshaerosol_t) _get_dl_function_pointer(handle,
                                                      lib_path,
                                                     "api_sshaerosol_gaschemistry",
                                                      true);
        fct();
      }

      // Call aerosols dynamic
      {
        typedef void (*api_sshaerosol_t)();
        api_sshaerosol_t fct =
          (api_sshaerosol_t) _get_dl_function_pointer(handle,
                                                      lib_path,
                                                     "api_sshaerosol_aerodyn",
                                                      true);
        fct();
      }

      //
      // A realistic external code would read the gas and aerosol concentrations here
      //

      // Ask SSH-aerosol to perform some IO
      {
        typedef void (*api_sshaerosol_t)();
        api_sshaerosol_t fct =
          (api_sshaerosol_t) _get_dl_function_pointer(handle,
                                                      lib_path,
                                                     "api_sshaerosol_output",
                                                      true);
        fct();
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
