// Copyright (C) 2005-2007, ENPC - INRIA - EDF R&D
// Author(s): Edouard Debry
//
// This file is part of the air quality modeling system Polyphemus.
//
// Polyphemus is developed in the INRIA - ENPC joint project-team CLIME and in
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


// Generates boundary conditions for Polair3D based on Gocart concentrations
// for particulate matter.


#include "Talos.hxx"
using namespace Talos;

const int Nquad = 500;

// Class defining one log-normal density.
class LogNormal
{
private:
  float Ntot;
  float Dmean;
  float Sigma;

public:
  LogNormal();
  LogNormal(float, float, float);
  float NumLogNormal(float);
  float VolLogNormal(float);
};


// Constructor.
LogNormal::LogNormal()
{
  Ntot = 0.0;
  Dmean = 10.0;
  Sigma = 10.0;
}


// Constructor.
LogNormal::LogNormal(float Ntot0, float Dmean0, float Sigma0)
{
  if (Ntot0 * Dmean0 * Sigma0 <= 0.0)
    throw string("WARNING: strange parameters in LogNormal\n")
      + string("WARNING: Ntot Dmean Sigma: ") + to_str(Ntot0)
      + string("\t") + to_str(Dmean0) + string("\t") + to_str(Sigma0);

  Ntot = Ntot0;
  Dmean = Dmean0;
  Sigma = Sigma0;
}


// Computes the number concentration at a given diameter (µm) for a lognormal
// distribution.
float LogNormal::NumLogNormal(float diameter)
{
  // 1 / sqrt(2 * pi)
  const float sqrt_2pi = 0.3989422804020;

  float tmp1 = log(Sigma);
  float tmp2 = log(diameter / Dmean) / tmp1;

  return sqrt_2pi * Ntot / tmp1 * exp(-0.5 * tmp2 * tmp2);
}


// Computes the volume concentration at a given diameter (µm) for a lognormal
// distribution.
float LogNormal::VolLogNormal(float diameter)
{
  // pi / 6
  const float pi_6 = 0.523598775598;

  return pi_6 * diameter * diameter * diameter
    * this->NumLogNormal(diameter);
}



// Class for one modal aerosol distribution.
class ModalAerosol
{
private:
  float Dmin, Dmax;
  LogNormal Nuclei;
  LogNormal Accumulation;
  LogNormal Coarse;
  LogNormal NucleiSoluble;
  LogNormal AitkenSoluble;
  LogNormal AccumulationSoluble;
  LogNormal CoarseSoluble;
  LogNormal AitkenInsoluble;
  LogNormal AccumulationInsoluble;
  LogNormal CoarseInsoluble;

public:
  ModalAerosol();
  ModalAerosol(float, float, float*, float*, float*);
  ModalAerosol(float*, float*, float*);
  ModalAerosol(float, float, float, float, float);
  float NumDensity(float);
  float VolDensity(float);
  float NumQuantity(float, float);
  float VolQuantity(float, float);
  float NumDensity7(float);
  float VolDensity7(float);
  float NumQuantity7(float, float);
  float VolQuantity7(float, float);
  float VolDensity1(float);
  float VolQuantity1(float, float);
};


// Default constructor.
ModalAerosol::ModalAerosol():
  Nuclei(), Accumulation(), Coarse()
{
  Dmin = 0.001;
  Dmax = 10.0;
}


// Constructor.
ModalAerosol::ModalAerosol(float Dmin0, float Dmax0, float *Ntot0,
                           float *Dmean0, float *Sigma0):
  Nuclei(Ntot0[0], Dmean0[0], Sigma0[0]),
  Accumulation(Ntot0[1], Dmean0[1], Sigma0[1]),
  Coarse(Ntot0[2], Dmean0[2], Sigma0[2])
{
  if (Dmin0 * Dmax0 <= 0.000001 || Dmax0 <= Dmin0)
    throw string("WARNING: strange parameters in ModalAerosol\n")
      + string("WARNING: Dmin Dmax: ") + to_str(Dmin0)
      + string("\t") + to_str(Dmax0);

  Dmin = Dmin0;
  Dmax = Dmax0;
}

// Constructor.
ModalAerosol::ModalAerosol(float *Ntot0, float *Dmean0, float *Sigma0):
  NucleiSoluble(Ntot0[0], Dmean0[0], Sigma0[0]),
  AitkenSoluble(Ntot0[1], Dmean0[1], Sigma0[1]),
  AccumulationSoluble(Ntot0[2], Dmean0[2], Sigma0[2]),
  CoarseSoluble(Ntot0[3], Dmean0[3], Sigma0[3]),
  AitkenInsoluble(Ntot0[4], Dmean0[4], Sigma0[4]),
  AccumulationInsoluble(Ntot0[5], Dmean0[5], Sigma0[5]),
  CoarseInsoluble(Ntot0[6], Dmean0[6], Sigma0[6])
{
}

// Constructor.
ModalAerosol::ModalAerosol(float Dmin0, float Dmax0, float Ntot0, float Dmean0, float Sigma0):
  Nuclei(Ntot0, Dmean0, Sigma0)
{
  if (Dmin0 * Dmax0 <= 0.000001 || Dmax0 <= Dmin0)
    throw string("WARNING: strange parameters in ModalAerosol\n")
      + string("WARNING: Dmin Dmax: ") + to_str(Dmin0)
      + string("\t") + to_str(Dmax0);

  Dmin = Dmin0;
  Dmax = Dmax0;
}

// Computes the aerosol number concentration for a modal aerosol distribution
// at a given diameter (µm).
float ModalAerosol::NumDensity(float diameter)
{
  float num = 0.0;
  if (diameter >= Dmin && diameter <= Dmax)
    num = Nuclei.NumLogNormal(diameter) + Accumulation.NumLogNormal(diameter)
      + Coarse.NumLogNormal(diameter);
  return num;
}

// Computes the aerosol number concentration for a modal aerosol distribution
// at a given diameter (µm) for 7 modes.
float ModalAerosol::NumDensity7(float diameter)
{
  float num = 0.0;
  if (diameter >= Dmin && diameter <= Dmax)
    num = NucleiSoluble.NumLogNormal(diameter)
      + AccumulationSoluble.NumLogNormal(diameter)
      + CoarseSoluble.NumLogNormal(diameter)
      + AitkenInsoluble.NumLogNormal(diameter)
      + AccumulationInsoluble.NumLogNormal(diameter)
      + CoarseInsoluble.NumLogNormal(diameter);

  return num;
}

// Computes the aerosol volume concentration for a modal aerosol distribution
// at a given diameter (µm).
float ModalAerosol::VolDensity(float diameter)
{
  float vol = 0.0;
  if (diameter >= Dmin && diameter <= Dmax)
    vol = Nuclei.VolLogNormal(diameter) + Accumulation.VolLogNormal(diameter)
      + Coarse.VolLogNormal(diameter);
  return vol;
}

// Computes the aerosol volume concentration for a modal aerosol distribution
// at a given diameter (µm) for 1 mode.
float ModalAerosol::VolDensity1(float diameter)
{
  float vol = 0.0;
  if (diameter >= Dmin && diameter <= Dmax)
    vol = Nuclei.VolLogNormal(diameter);

  return vol;
}

// Computes the aerosol volume concentration for a modal aerosol distribution
// at a given diameter (µm) for 7 modes.
float ModalAerosol::VolDensity7(float diameter)
{
  float vol = 0.0;
  if (diameter >= Dmin && diameter <= Dmax)
    vol = NucleiSoluble.VolLogNormal(diameter)
      + AccumulationSoluble.VolLogNormal(diameter)
      + CoarseSoluble.VolLogNormal(diameter)
      + AitkenInsoluble.VolLogNormal(diameter)
      + AccumulationInsoluble.VolLogNormal(diameter)
      + CoarseInsoluble.VolLogNormal(diameter);
  return vol;
}


// Computes aerosol number quantity of modal
// distribution within an aerosol diameter range.
float ModalAerosol::NumQuantity(float diameter1, float diameter2)
{
  float log_diameter1 = log(diameter1);
  float log_diameter2 = log(diameter2);
  float log_difference = log_diameter2 - log_diameter1;
  float log_h = log_difference / (float) Nquad;

  float num = 0.0;
  if (log_difference > 0.0)
    {
      num = (this->NumDensity(diameter1) + this->NumDensity(diameter2)) * 0.5;
      for (int j = 1; j < Nquad; j++)
        {
          float log_diameter = log_diameter1 + log_h * (float) j;
          float diameter = exp(log_diameter);

          num += this->NumDensity(diameter);
        }
    }
  return num * log_h;
}

// Computes aerosol number quantity of modal
// distribution within an aerosol diameter range for 7 modes.
float ModalAerosol::NumQuantity7(float diameter1, float diameter2)
{
  float log_diameter1 = log(diameter1);
  float log_diameter2 = log(diameter2);
  float log_difference = log_diameter2 - log_diameter1;
  float log_h = log_difference / (float) Nquad;

  float num = 0.0;
  if (log_difference > 0.0)
    {
      num = (this->NumDensity7(diameter1) + this->NumDensity7(diameter2)) * 0.5;
      for (int j = 1; j < Nquad; j++)
        {
          float log_diameter = log_diameter1 + log_h * (float) j;
          float diameter = exp(log_diameter);

          num += this->NumDensity7(diameter);
        }
    }
  return num * log_h;
}

// Computes aerosol volume quantity of modal
// distribution within an aerosol diameter range
float ModalAerosol::VolQuantity(float diameter1, float diameter2)
{
  float log_diameter1 = log(diameter1);
  float log_diameter2 = log(diameter2);
  float log_difference = log_diameter2 - log_diameter1;
  float log_h = log_difference / (float) Nquad;

  float vol = 0.0;
  if (log_difference > 0.0)
    {
      vol = (this->VolDensity(diameter1) + this->VolDensity(diameter2)) * 0.5;
      for (int j = 1; j < Nquad; j++)
        {
          float log_diameter = log_diameter1 + log_h * (float) j;
          float diameter = exp(log_diameter);

          vol += this->VolDensity(diameter);
        }
    }
  return vol * log_h;
}

// Computes aerosol volume quantity of modal
// distribution within an aerosol diameter range
float ModalAerosol::VolQuantity1(float diameter1, float diameter2)
{
  float log_diameter1 = log(diameter1);
  float log_diameter2 = log(diameter2);
  float log_difference = log_diameter2 - log_diameter1;
  float log_h = log_difference / (float) Nquad;

  float vol = 0.0;
  if (log_difference > 0.0)
    {
      vol = (this->VolDensity1(diameter1) + this->VolDensity1(diameter2)) * 0.5;
      for (int j = 1; j < Nquad; j++)
        {
          float log_diameter = log_diameter1 + log_h * (float) j;
          float diameter = exp(log_diameter);

          vol += this->VolDensity1(diameter);
        }
    }
  return vol * log_h;
}

// Computes aerosol volume quantity of modal
// distribution within an aerosol diameter range
float ModalAerosol::VolQuantity7(float diameter1, float diameter2)
{
  float log_diameter1 = log(diameter1);
  float log_diameter2 = log(diameter2);
  float log_difference = log_diameter2 - log_diameter1;
  float log_h = log_difference / (float) Nquad;

  float vol = 0.0;
  if (log_difference > 0.0)
    {
      vol = (this->VolDensity7(diameter1) + this->VolDensity7(diameter2)) * 0.5;
      for (int j = 1; j < Nquad; j++)
        {
          float log_diameter = log_diameter1 + log_h * (float) j;
          float diameter = exp(log_diameter);

          vol += this->VolDensity7(diameter);
        }
    }
  return vol * log_h;
}

// Class of bin discretization.
class BinDist
{
private:
  int Nbin;
  float Dmin, Dmax;
  float* bound;

public:
  BinDist();
  BinDist(int, float, float);
  ~BinDist();
  float GetLowBound(int);
  float GetUpBound(int);
};


// Default constructor.
BinDist::BinDist()
{
  Nbin = 0;
  Dmin = 0.001;
  Dmax = 10.0;
  bound = NULL;
}


// Constructor.
BinDist::BinDist(int Nbin0, float Dmin0, float Dmax0)
{
  if (Dmin0 * Dmax0 <= 0.000001 || Dmax0 <= Dmin0 || Nbin0 <= 0)
    throw string("WARNING: strange parameters in BinDist\n")
      + string("WARNING: Dmin Dmax Nbin: ") + to_str(Dmin0) + string("\t")
      + to_str(Dmax0) + string("\t") + to_str(Nbin0);

  Dmin = Dmin0;
  Dmax = Dmax0;
  Nbin = Nbin0;

  bound = new float[Nbin + 1];

  float diameter_difference = log(Dmax / Dmin) / float(Nbin) ;

  for (int j = 0; j < Nbin + 1; j++)
    bound[j] = Dmin * exp(diameter_difference * float(j));
}


// Gets the lower bound of bin #i.
float BinDist::GetLowBound(int i)
{
  if (i >= 0 && i < Nbin + 1)
    return bound[i];
  else
    return Dmin;
}


// Gets the upper bound of bin #i.
float BinDist::GetUpBound(int i)
{
  if (i >= 0 && i < Nbin + 1)
    return bound[i + 1];
  else
    return Dmax;
}


// Destructor.
BinDist::~BinDist()
{
  if (bound != NULL)
    delete[] bound;
  Nbin = 0;
  Dmin = 0.0;
  Dmax = 0.0;
}
