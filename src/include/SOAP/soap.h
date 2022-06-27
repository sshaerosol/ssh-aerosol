//!!-----------------------------------------------------------------------
//!!     Copyright (C) 2019 CEREA (ENPC) - INERIS
//!!     SSH-aerosol is distributed under the GNU General Public License v3
//!!-----------------------------------------------------------------------

//////////////
// INCLUDES //

#include <vector>
#include "species.h"
#include "solving.cxx"
using namespace ssh_soap;

void soap_main_ssh(double LWC, double RH, double Temperature, 
		   double ionic, double chp, double& LWCorg,
		   double& deltat,
		   double DSD[], double csol[], double liquid[],
		   int ns_aer, int ns_aer_layers, int neq, double q[], double qaero[], double qaq[], 
		   double qgas[],
		   double lwc_Nsize[], double ionic_Nsize[], double chp_Nsize[],
		   double liquid_Nsize[], int nbin, int isoapdyn, 
                   int imethod, int soap_inorg,
		   char species_name[], int name_len, double molecular_weight_aer[],
		   double accomodation_coefficient[], int aerosol_type[], 
		   char partitioning[], char smiles[], double saturation_vapor_pressure[], double enthalpy_vaporization[],
                   double diffusion_coef[], int nlayer, int with_kelvin_effect,
		   double tequilibrium, double dtaeromin, double dorg,
		   int coupled_phases, int activity_model, double epser_soap, int i_hydrophilic,
		   int N_inert, int N_inorganic, int with_oligomerization);
