// Copyright (C) 2007, INRIA
// Author(s): Meryem Ahmed de Biasi
//
// This file is part of SeldonData library, used for data processing.
//
// SeldonData is free software; you can redistribute it and/or modify it under
// the terms of the GNU General Public License as published by the Free
// Software Foundation; either version 2 of the License, or (at your option)
// any later version.
//
// SeldonData is distributed in the hope that it will be useful, but WITHOUT
// ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
// FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for
// more details.
//
// For more information, visit the SeldonData home page:
// http://vivienmallet.net/lib/seldondata/

// Program to decode GRIB files.
// It is based on and makes calls to WGRIB, by Wesley Ebisuzaki.
// WGRIB home page : http://www.cpc.ncep.noaa.gov/products/wesley/wgrib.html


#include <float.h>
#include <cstdlib>
#include <stdio.h>
#include <string>

#include "wgrib_globals.h"
#include "pds4.h"
#include "gds.h"
#include "bms.h"
#include "bds.h"
#include "cnames.h"

#define MSEEK        1024
#define BUFF_ALLOC0  40000

extern "C"
{
  unsigned char *seek_grib(FILE *file, long *pos, long *len_grib,
                           unsigned char *buffer, unsigned int buf_len);
  double ibm2flt(unsigned char *ibm);
  int GDS_grid(unsigned char *gds, unsigned char *bds, int *nx, int *ny,
               long int *nxny);
  int read_grib(FILE *file, long pos, long len_grib, unsigned char *buffer);
  double ibm2flt(unsigned char *ibm);
  void BDS_unpack(float *flt, unsigned char *bds, unsigned char *bitmap,
                  int n_bits, int n, double ref, double scale);
  int BDS_NValues(unsigned char *bds);
  double int_power(double x, int y);
  int missing_points(unsigned char *bitmap, int n);
}

int decode_grib(FILE *grib_file, long *pos, int *param, float **array,
                int *nx, int *ny)
{
  unsigned char *buffer, *msg, *pds, *gds, *bms, *bds, *pointer;
  long len_grib, buffer_size, count = 1;
  int i, return_code = 0;
  long int nxny;
  double temp;

  if ((buffer = (unsigned char *) malloc(BUFF_ALLOC0)) == NULL)
    {
      fprintf(stderr, "not enough memory\n");
    }

  buffer_size = BUFF_ALLOC0;

  msg = seek_grib(grib_file, pos, &len_grib, buffer, MSEEK);
  if (msg == NULL)
    {
      fprintf(stderr, "missing GRIB record(s)\n");
      return 1;
    }

  /* read all whole grib record */
  if (len_grib + msg - buffer > buffer_size)
    {
      buffer_size = len_grib + msg - buffer + 1000;
      buffer = (unsigned char *) realloc((void *) buffer, buffer_size);
      if (buffer == NULL)
        {
          fprintf(stderr, "ran out of memory\n");
          return 1;
        }
    }

  if (read_grib(grib_file, *pos, len_grib, buffer) == 0)
    {
      fprintf(stderr, "error, could not read to end of record %ld\n", count);
      return 1;
    }

  msg = buffer;
  pds = (msg + 8);
  pointer = pds + PDS_LEN(pds);

  if (PDS_HAS_GDS(pds))
    {
      gds = pointer;
      pointer += GDS_LEN(gds);
    }
  else
    {
      gds = NULL;
    }

  if (PDS_HAS_BMS(pds))
    {
      bms = pointer;
      pointer += BMS_LEN(bms);
    }
  else
    {
      bms = NULL;
    }

  bds = pointer;
  pointer += BDS_LEN(bds);

  if (pointer[0] != 0x37 || pointer[1] != 0x37 ||
      pointer[2] != 0x37 || pointer[3] != 0x37)
    {
      fprintf(stderr, "\n\n    missing end section\n");
      fprintf(stderr, "%2x %2x %2x %2x\n", pointer[0], pointer[1],
              pointer[2], pointer[3]);
      return 1;
    }

  if (gds != NULL)
    GDS_grid(gds, bds, nx, ny, &nxny);
  else if (bms != NULL)
    {
      nxny = *nx = BMS_nxny(bms);
      *ny = 1;
    }
  else
    {
      if (BDS_NumBits(bds) == 0)
        {
          nxny = *nx = 1;
          fprintf(stderr, "Missing GDS, constant record .. cannot "
                  "determine number of data points\n");
        }
      else
        {
          nxny = *nx = BDS_NValues(bds);
        }
      *ny = 1;
    }

  if (gds && ! GDS_Harmonic(gds))
    {
      if (BDS_NumBits(bds) != 0)
        {
          i = BDS_NValues(bds);
          if (bms != NULL)
            {
              i += missing_points(BMS_bitmap(bms), nxny);
            }
          if (i != nxny)
            {
              fprintf(stderr, "grib header at record %ld: two values of nxny %ld %d\n",
                      count, nxny, i);
              fprintf(stderr, "   LEN %d DataStart %d UnusedBits %d #Bits %d nxny %ld\n",
                      BDS_LEN(bds), BDS_DataStart(bds), BDS_UnusedBits(bds),
                      BDS_NumBits(bds), nxny);
              return_code = 1;
              nxny = *nx = i;
              *ny = 1;
            }
        }
    }

  *param = PDS_PARAM(pds);


  if ((*array = (float *) malloc(sizeof(float) * nxny)) == NULL)
    {
      fprintf(stderr, "memory problems\n");
      exit(8);
    }
  temp = int_power(10.0, - PDS_DecimalScale(pds));
  BDS_unpack(*array, bds, BMS_bitmap(bms), BDS_NumBits(bds), nxny,
             temp * BDS_RefValue(bds),
             temp * int_power(2.0, BDS_BinScale(bds)));

  *pos += len_grib;
  count ++;

  free(buffer);

  return return_code;
}
