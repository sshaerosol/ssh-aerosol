//!!-----------------------------------------------------------------------
//!!     Copyright (C) 2019 CEREA (ENPC) - INERIS
//!!     SSH-aerosol is distributed under the GNU General Public License v3
//!!-----------------------------------------------------------------------

#ifndef SPECIES_CXX
#define SPECIES_CXX


#include <vector>
#include "species.h"
using namespace ssh_soap;
using namespace std;


void add_species_ssh( vector<species>& surrogate, species current_species, 
		      vector<string> species_list_aer,
		      double molecular_weight_aer[],
		      double accomodation_coefficient[],
		      double diffusion_coef[],
		      vector<string> species_part,
		      int nlayer, int i_hydrophilic,
	              int N_inert, int N_inorganic)
{

  int nsp = species_list_aer.size();
  int N_start = N_inert + N_inorganic;

  // Find the number in the aerosol species list
  current_species.soap_ind = -1;
  current_species.soap_ind_aero = -1;
  for (int i = 0; i < nsp; ++i)
    if (species_list_aer[i].substr(1,-1) == current_species.name)
      {
        current_species.soap_ind = i;
        current_species.soap_ind_aero = (i-N_start) * (nlayer-1+i_hydrophilic) + i;
        current_species.MM =  molecular_weight_aer[i] / 1.e6;
        current_species.accomodation_coefficient = accomodation_coefficient[i];
	current_species.is_generic=false;
	current_species.KDiffusion_air = 1.0e-5; //diffusion_coef[i];
	
	if (species_part[i] == "HPHO")
          {
            current_species.hydrophilic=false;
	    current_species.hydrophobic=true;
          }
        else if (species_part[i] == "HPHI")
          {
            current_species.hydrophilic=true;
	    current_species.hydrophobic=false;
          }  
        else if (species_part[i] == "BOTH")
          {
            current_species.hydrophilic=true;
            current_species.hydrophobic=true;
          }  
        else
          {
            cout << "WARNING: partitioning of species " << current_species.name << " not defined in species-list-aer file" << endl;
            cout << "Default partitioning processes for " << current_species.name << " is hydrophilic=" << current_species.hydrophilic  << " and hydrophobic="  << current_species.hydrophobic << endl;
          }
      }

  // if (current_species.soap_ind == -1)
  //   cout << "Warning: " << current_species.name << " is not found in " <<
  //     "the aerosol species list.\n";
  // else
  //   surrogate.push_back(current_species);

  if (current_species.soap_ind != -1)
    surrogate.push_back(current_species);
}

void add_generic_species_ssh(model_config &config, 
			     vector<species>& surrogate, 
			     vector<string> species_list_aer,
			     double molecular_weight_aer[],
			     double accomodation_coefficient[],
			     int aerosol_type[],
			     vector<string> species_smiles,
			     double saturation_vapor_pressure[],
			     double enthalpy_vaporization[],
			     vector<string> species_part,
			     int nlayer, int i_hydrophilic,
	                     int N_inert, int N_inorganic)
{

  int j,n,found;
  n=surrogate.size();
  int nsp = species_list_aer.size();
  int N_start = N_inert + N_inorganic;
  
  config.chemistry=false;
  // Find the number in the aerosol species list  
  for (int i = 0; i < nsp; ++i)
    if (aerosol_type[i]==4 and saturation_vapor_pressure[i]>0. and species_smiles[i]!="-")
      {
	found=0;
	for (j=0;j<n;j++)
	  if (species_list_aer[i].substr(1,-1) == surrogate[j].name)
	    {
	      found=1;
	    }      
      
	if (found==0)
	  {	  
	    species X;
	    X.name=species_list_aer[i].substr(1,-1);
	    X.smile=species_smiles[i];
	    X.is_inorganic_precursor=false;
	    X.Psat_ref=saturation_vapor_pressure[i]; // Saturation vapor pressure at Tref (torr)
	    X.kp_from_experiment=false;  // Use experimental partitioning constant at Tref?
	    X.Tref=298;         // Temperature of reference (K)
	    X.deltaH=enthalpy_vaporization[i];     // Enthalpy of vaporization (kJ/mol)
	    X.Henry=0.;     // Henry's law constant at Tref (M/atm)
	    X.aq_type="none"; // "none","diacid","monoacid" or "aldehyde"
	    
	    if (species_part[i] == "HPHO")
              {
                X.hydrophilic=false;
	        X.hydrophobic=true;
              }
            else if (species_part[i] == "HPHI")
              {
                X.hydrophilic=true;
	        X.hydrophobic=false;
              }  
            else if (species_part[i] == "BOTH")
              {
                X.hydrophilic=true;
	        X.hydrophobic=true;
              }  
            else
              {
                cout << "WARNING: partitioning of generic species " << species_list_aer[i].substr(1,-1) << " not defined." << endl;
                cout << species_list_aer[i].substr(1,-1) << " is considered both hydrophilic and hydrophobic" << endl;
                X.hydrophilic=true;
	        X.hydrophobic=true;
              }
              
	    X.nonvolatile=false; // Is the compound nonvolatile?
	    X.is_organic=true;  // Is the compound organic?
	    X.compute_gamma_org=true;  // Compute the activity coefficients of the organic phase for this compound?
	    X.compute_gamma_aq=true;  // Compute the activity coefficients of the aqueous phase for this compound?
	    X.Koligo_org=0.0;         //oligomeriation constant in the organic phase
	    X.rho=1300.0;
	    X.is_monomer=false;
	    X.rion=false;
	    X.KDiffusion_air=1.0e-5;
	    // BiA2D.accomodation_coefficient=alpha;
	    X.viscosity=1.68e12;
	    X.is_solid=false;	  
	    X.MM =  molecular_weight_aer[i] / 1.e6;
	    X.accomodation_coefficient = accomodation_coefficient[i];
	    X.soap_ind = i;
	    X.soap_ind_aero = (i-N_start) * (nlayer-1+i_hydrophilic) + i;
	    X.is_generic=true;
	    surrogate.push_back(X);	  
	  }
      }
    else if (species_list_aer[i].substr(1,5) == "Oligo")
      {
	int ij=-1;
	for (j=0;j<n;j++)
	  if (species_list_aer[i].substr(1,-1) == "Oligo"+surrogate[j].name)
	    {
	      ij=j;
	    }

	if (ij>0)
	  {
	    config.chemistry=true;
	    
	    surrogate[ij].is_monomer=true;
	    surrogate[ij].name_oligomer=species_list_aer[i].substr(1,-1);
	    surrogate[ij].moligo=config.moligo;
	    
	    species X;
	    X.name=species_list_aer[i].substr(1,-1);
	    X.is_inorganic_precursor=false;
	    X.smile="";
	    X.aq_type="none";
	    X.hydrophilic=surrogate[ij].hydrophilic;   // Does the species condense on the aqueous phase?
	    X.hydrophobic=surrogate[ij].hydrophobic;  // Does the species condense on the organic phase?
	    X.nonvolatile=true; // Is the compound nonvolatile?
	    X.kp_from_experiment=false;
	    X.is_organic=true;  // Is the compound organic?
	    X.compute_gamma_org=surrogate[ij].compute_gamma_org;  // Compute the activity coefficients of the organic phase for this compound?
	    X.compute_gamma_aq=surrogate[ij].compute_gamma_aq;  // Compute the activity coefficients of the aqueous phase for this compound?
	    X.Koligo_org=0.0;         //oligomeriation constant in the organic phase
	    X.rho=1300.0;
	    X.is_monomer=false;
	    X.rion=false;
	    X.KDiffusion_air=1.0e-5;
	    // BiA2D.accomodation_coefficient=alpha;
	    X.viscosity=1.68e12;
	    X.is_solid=false;	  
	    X.MM =  molecular_weight_aer[i] / 1.e6;
	    X.accomodation_coefficient = accomodation_coefficient[i];
	    X.soap_ind = i;
	    X.soap_ind_aero = (i-N_start) * (nlayer-1+i_hydrophilic) + i;
	    X.is_generic=false;
	    surrogate.push_back(X);	    

	  }
      }
}
  
void creation_species_ssh( model_config &config, vector<species>& surrogate, vector<string> species_list_aer,
			   double molecular_weight_aer[], double accomodation_coefficient[],
			   int aerosol_type[],
			   vector<string> species_smiles, double saturation_vapor_pressure[],
			   double enthalpy_vaporization[], 
			   double diffusion_coef[], vector<string> species_part,
			   int nlayer, int i_hydrophilic, bool compute_inorganic, int N_inert, int N_inorganic,
                           int with_oligomerization)
{
  int nsp = species_list_aer.size();
  // double alpha = 1.0; //0.01; // accommodation coefficient

  species BiA2D;
  BiA2D.name="BiA2D";
  BiA2D.is_inorganic_precursor=false;
  BiA2D.Psat_ref=1.43e-7; // Saturation vapor pressure at Tref (torr)
  BiA2D.kp_from_experiment=false;  // Use experimental partitioning constant at Tref?
  BiA2D.Tref=298;         // Temperature of reference (K)
  BiA2D.deltaH=109.0;     // Enthalpy of vaporization (kJ/mol)
  BiA2D.Henry=2.67e8;     // Henry's law constant at Tref (M/atm)
  BiA2D.aq_type="diacid"; // "none","diacid","monoacid" or "aldehyde"
  BiA2D.Kacidity1=3.95e-4;    // First acidity constant
  BiA2D.Kacidity2=7.70e-6;    // Second acidity constant
  BiA2D.hydrophilic=true;   // Does the species condense on the aqueous phase?
  BiA2D.hydrophobic=false;  // Does the species condense on the organic phase?
  BiA2D.nonvolatile=false; // Is the compound nonvolatile?
  BiA2D.is_organic=true;  // Is the compound organic?
  BiA2D.compute_gamma_org=true;  // Compute the activity coefficients of the organic phase for this compound?
  BiA2D.compute_gamma_aq=true;  // Compute the activity coefficients of the aqueous phase for this compound?
  BiA2D.Koligo_org=0.0;         //oligomeriation constant in the organic phase
  BiA2D.rho=1300.0;
  BiA2D.is_monomer=false;
  BiA2D.rion=false;
  //BiA2D.KDiffusion_air=1.0e-5;
  // BiA2D.accomodation_coefficient=alpha;
  BiA2D.viscosity=1.68e12;
  BiA2D.is_solid=false;
  BiA2D.is_generic=false;

  //Group: if no functionnal group in the species use the default species
  //for the computation of activity coefficients
  //
  double group_tmp_bia2d [] = {2.0,2.0,2.0,1.0, // group C
			       0.0,0.0,0.0,0.0, //group C[OH]
			       0.0,0.0,0.0,0.0, //group Calcohol
			       0.0,0.0,0.0,0.0, //group Calcohol-tail
			       0.0,0.0,0.0,0.0,0.0, //group C=C
			       0.0,0.0, //group aromatic carbon (AC)
			       0.0,0.0,0.0, // group //AC-C
			       0.0,  //group OH
			       0.0, //group H2O
			       0.0, //group ACOH
			       0.0,0.0, //group ketone
			       0.0,   //group aldehyde  
			       0.0,0.0, //group ester
			       0.0,0.0,0.0, //group ether 
			       2.0,  //group acid
			       0.0,   //group ACNO2
			       0.0,0.0,0.0, //group NO3
			       0.0,0.0,0.0, //group CO-OH
			       0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0, //group CO-OC
			       0.0,  //group PAN
			       0.0,  //group CO-OOH
			       0.0,  //group O=COC=O
			       0.0,0.0,0.0}; //group CHxNO2
  
  int size = sizeof(group_tmp_bia2d)/sizeof(double);
  assert(size == 60);
	
  for(int i = 0; i < size; ++i)
    BiA2D.groups[i] = group_tmp_bia2d[i];

  // Search the species name in the aerosol species list 
  // and add the species if its name matches with
  // the given list.
 
  add_species_ssh(surrogate, BiA2D, species_list_aer, molecular_weight_aer,
		  accomodation_coefficient,diffusion_coef,species_part,nlayer,i_hydrophilic,
		  N_inert, N_inorganic);

  /* ==== BiA1D ==== */ 

  species BiA1D;
  BiA1D.name="BiA1D";
  BiA1D.is_inorganic_precursor=false;
  BiA1D.Psat_ref=2.17e-7; // Saturation vapor pressure at Tref (torr)
  BiA1D.Tref=298;         // Temperature of reference (K)
  BiA1D.deltaH=50.0;     // Enthalpy of vaporization (kJ/mol)
  BiA1D.Henry=1.12e8;     // Henry's law constant at Tref (M/atm)
  BiA1D.aq_type="monoacid"; // "none","diacid","monoacid" or "aldehyde"
  BiA1D.Kacidity1=6.52e-4;    // First acidity constant
  BiA1D.hydrophilic=true;   // Does the species condense on the aqueous phase?
  BiA1D.hydrophobic=false;  // Does the species condense on the organic phase?
  BiA1D.nonvolatile=false; // Is the compound nonvolatile?
  BiA1D.kp_from_experiment=false;  // Use experimental partitioning constant at Tref?
  BiA1D.is_organic=true;  // Is the compound organic?
  BiA1D.compute_gamma_org=true;  // Compute the activity coefficients of the organic phase for this compound?
  BiA1D.compute_gamma_aq=true;  // Compute the activity coefficients of the aqueous phase for this compound?
  BiA1D.Koligo_org=0.0;      //oligomeriation constant in the organic phase
  BiA1D.rho=1300.0;  
  BiA1D.is_monomer=false;
  BiA1D.rion=false;
  //BiA1D.KDiffusion_air=1.0e-5;
  //  BiA1D.accomodation_coefficient=alpha;
  BiA1D.viscosity=1.68e12;
  BiA1D.is_solid=false;
  BiA1D.is_generic=false;

  //Group: if no functionnal group in the species use the default species
  //for the computation of activity coefficients
  //
  double group_tmp_bia1d [] = {2.0,1.0,2.0,1.0, // group C
			       0.0,0.0,0.0,0.0, //group C[OH]
			       0.0,0.0,0.0,0.0, //group Calcohol
			       0.0,0.0,0.0,0.0, //group Calcohol-tail
			       0.0,0.0,0.0,0.0,0.0, //group C=C
			       0.0,0.0, //group aromatic carbon (AC)
			       0.0,0.0,0.0, // group //AC-C
			       0.0,  //group OH
			       0.0, //group H2O
			       0.0, //group ACOH
			       1.0,0.0, //group ketone
			       0.0,   //group aldehyde  
			       0.0,0.0, //group ester
			       0.0,0.0,0.0, //group ether 
			       1.0,  //group acid
			       0.0,   //group ACNO2
			       0.0,0.0,0.0, //group NO3
			       0.0,0.0,0.0, //group CO-OH
  			       0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0, //group CO-OC
			       0.0,    //group PAN
			       0.0,  //group CO-OOH
			       0.0,  //group O=COC=O
			       0.0,0.0,0.0}; //group CHxNO2
  
  size = sizeof(group_tmp_bia1d)/sizeof(double);
  assert(size == 60);
	
  for(int i = 0; i < size; ++i)
    BiA1D.groups[i] = group_tmp_bia1d[i];

  // Search the species name in the aerosol species list 
  // and add the species if its name matches with
  // the given list.
  add_species_ssh(surrogate, BiA1D, species_list_aer, molecular_weight_aer,
              accomodation_coefficient,diffusion_coef,species_part,nlayer,i_hydrophilic,
		  N_inert, N_inorganic);

  species BiA0D;
  BiA0D.name="BiA0D";
  BiA0D.is_inorganic_precursor=false;
  BiA0D.Psat_ref=2.7e-4; // Saturation vapor pressure at Tref (torr)
  BiA0D.Tref=298;         // Temperature of reference (K)
  BiA0D.deltaH=50.0;     // Enthalpy of vaporization (kJ/mol)
  BiA0D.Henry=1.98e6;     // Henry's law constant at Tref (M/atm)
  BiA0D.Koligo_org=0.0;    //oligomeriation constant in the organic phase
  BiA0D.hydrophilic=true;   // Does the species condense on the aqueous phase?
  BiA0D.hydrophobic=false;  // Does the species condense on the organic phase?
  BiA0D.nonvolatile=false; // Is the compound nonvolatile?
  BiA0D.kp_from_experiment=false;  // Use experimental partitioning constant at Tref?
  BiA0D.is_organic=true;  // Is the compound organic?
  BiA0D.compute_gamma_org=true;  // Compute the activity coefficients of the organic phase for this compound?
  BiA0D.compute_gamma_aq=true;  // Compute the activity coefficients of the aqueous phase for this compound?
  //Parameters for the oligomerization of aldehyde in the aqueous phase:
  if(with_oligomerization == 1) {
         BiA0D.Koligo_aq=0.1;     
         BiA0D.aq_type="aldehyde";} // "none","diacid","monoacid" or "aldehyde"
  else {
         BiA0D.Koligo_aq=0.;     
         BiA0D.aq_type="none"; }// "none","diacid","monoacid" or "aldehyde"
  BiA0D.pHref=6.0;
  BiA0D.beta=1.91;
  BiA0D.rho=1300.0;
  BiA0D.is_monomer=false;
  BiA0D.rion=false;
  //BiA0D.KDiffusion_air=1.0e-5;
  //  BiA0D.accomodation_coefficient=alpha;
  BiA0D.viscosity=1.68e12;
  BiA0D.is_solid=false;
  BiA0D.is_generic=false;
  
  //Group: if no functionnal group in the species use the default species
  //for the computation of activity coefficients
  //
  double group_tmp_bia0d [] = {2.0,2.0,2.0,1.0, // group C
			       0.0,0.0,0.0,0.0, //group C[OH]
			       0.0,0.0,0.0,0.0, //group Calcohol
			       0.0,0.0,0.0,0.0, //group Calcohol-tail
			       0.0,0.0,0.0,0.0,0.0, //group C=C
			       0.0,0.0, //group aromatic carbon (AC)
			       0.0,0.0,0.0, // group //AC-C
			       0.0,  //group OH
			       0.0, //group H2O
			       0.0, //group ACOH
			       1.0,0.0, //group ketone
			       1.0,   //group aldehyde  
			       0.0,0.0, //group ester
			       0.0,0.0,0.0, //group ether 
			       0.0,  //group acid
			       0.0,   //group ACNO2
			       0.0,0.0,0.0, //group NO3
			       0.0,0.0,0.0, //group CO-OH
  			       0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0, //group CO-OC
			       0.0,  //group PAN
			       0.0,  //group CO-OOH
			       0.0,  //group O=COC=O
			       0.0,0.0,0.0}; //group CHxNO2

  size = sizeof(group_tmp_bia0d)/sizeof(double);
  assert(size == 60);
	
  for(int i = 0; i < size; ++i)
    BiA0D.groups[i] = group_tmp_bia0d[i];

  // Search the species name in the aerosol species list 
  // and add the species if its name matches with
  // the given list.
  add_species_ssh(surrogate, BiA0D, species_list_aer, molecular_weight_aer,
		  accomodation_coefficient,diffusion_coef,species_part,nlayer,i_hydrophilic,
		  N_inert, N_inorganic);
  
  species BiMT;
  BiMT.name="BiMT";
  BiMT.is_inorganic_precursor=false;
  BiMT.Psat_ref=1.45e-6; // Saturation vapor pressure at Tref (torr)
  BiMT.Tref=298;         // Temperature of reference (K)
  BiMT.deltaH=38.4;     // Enthalpy of vaporization (kJ/mol)
  BiMT.Henry=33.e9;     // Henry's law constant at Tref (M/atm)
  BiMT.aq_type="none"; // "none","diacid","monoacid" or "aldehyde"
  BiMT.hydrophilic=true;   // Does the species condense on the aqueous phase?
  BiMT.hydrophobic=false;  // Does the species condense on the organic phase?
  BiMT.nonvolatile=false; // Is the compound nonvolatile?
  BiMT.kp_from_experiment=false;  // Use experimental partitioning constant at Tref?
  BiMT.is_organic=true;  // Is the compound organic?
  BiMT.compute_gamma_org=true;  // Compute the activity coefficients of the organic phase for this compound?
  BiMT.compute_gamma_aq=true;  // Compute the activity coefficients of the aqueous phase for this compound?
  BiMT.Koligo_org=0.0;        //oligomeriation constant in the organic phase
  BiMT.rho=1300.0;
  BiMT.is_monomer=false;
  BiMT.rion=false;
  //BiMT.KDiffusion_air=1.0e-5;
  //  BiMT.accomodation_coefficient=alpha;
  BiMT.viscosity=1.68e12;
  BiMT.is_solid=false;
  BiMT.is_generic=false;
 
  //Group: if no functionnal group in the species use the default species
  //for the computation of activity coefficients
  //
  double group_tmp_bimt [] = {0.0,0.0,0.0,0.0, // group C
			      0.0,2.0,1.0,1.0, //group C[OH]
			      1.0,0.0,0.0,0.0, //group Calcohol
			      0.0,0.0,0.0,0.0, //group Calcohol-tail
			      0.0,0.0,0.0,0.0,0.0, //group C=C
			      0.0,0.0, //group aromatic carbon (AC)
			      0.0,0.0,0.0, // group //AC-C
			      4.0,  //group OH
			      0.0, //group H2O
			      0.0, //group ACOH
			      0.0,0.0, //group ketone
			      0.0,   //group aldehyde  
			      0.0,0.0, //group ester
			      0.0,0.0,0.0, //group ether 
			      0.0,  //group acid
			      0.0,   //group ACNO2
			      0.0,0.0,0.0, //group NO3
			      0.0,0.0,0.0, //group CO-OH
			      0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0, //group CO-OC
			      0.0,  //group PAN
			      0.0,  //group CO-OOH
			       0.0,  //group O=COC=O
			       0.0,0.0,0.0}; //group CHxNO2

  size = sizeof(group_tmp_bimt)/sizeof(double);
  assert(size == 60);
	
  for(int i = 0; i < size; ++i)
    BiMT.groups[i] = group_tmp_bimt[i];

  // Search the species name in the aerosol species list 
  // and add the species if its name matches with
  // the given list.  
  add_species_ssh(surrogate, BiMT, species_list_aer, molecular_weight_aer,
		  accomodation_coefficient,diffusion_coef,species_part,nlayer,i_hydrophilic,
		  N_inert, N_inorganic);

  species BiPER;
  BiPER.name="BiPER";
  BiPER.is_inorganic_precursor=false;
  BiPER.Psat_ref=2.61e-6; // Saturation vapor pressure at Tref (torr)
  BiPER.Tref=298;         // Temperature of reference (K)
  BiPER.deltaH=38.4;     // Enthalpy of vaporization (kJ/mol)
  BiPER.Henry=8.1e9;     // Henry's law constant at Tref (M/atm)
  BiPER.aq_type="none"; // "none","diacid","monoacid" or "aldehyde"
  BiPER.hydrophilic=true;   // Does the species condense on the aqueous phase?
  BiPER.hydrophobic=false;  // Does the species condense on the organic phase?
  BiPER.nonvolatile=false; // Is the compound nonvolatile?
  BiPER.kp_from_experiment=false;  // Use experimental partitioning constant at Tref?
  BiPER.is_organic=true;  // Is the compound organic?
  BiPER.compute_gamma_org=true;  // Compute the activity coefficients of the organic phase for this compound?
  BiPER.compute_gamma_aq=true;  // Compute the activity coefficients of the aqueous phase for this compound?
  BiPER.Koligo_org=0.0;         //oligomeriation constant in the organic phase
  BiPER.rho=1300.0;
  BiPER.is_monomer=false;
  BiPER.rion=false;
  //BiPER.KDiffusion_air=1.0e-5;
  //  BiPER.accomodation_coefficient=alpha;
  BiPER.viscosity=1.68e12;
  BiPER.is_solid=false;
  BiPER.is_generic=false;
  
  //Group: if no functionnal group in the species use the default species
  //for the computation of activity coefficients
  //
  double group_tmp_biper [] = {0.0,0.0,0.0,0.0, // group C
			       0.0,0.0,1.0,1.0, //group C[OH]
			       1.0,0.0,0.0,0.0, //group Calcohol
			       0.0,0.0,0.0,0.0, //group Calcohol-tail
			       0.0,0.0,0.0,0.0,0.0, //group C=C
			       0.0,0.0, //group aromatic carbon (AC)
			       0.0,0.0,0.0, // group //AC-C
			       2.0,  //group OH
			       0.0, //group H2O
			       0.0, //group ACOH
			       0.0,0.0, //group ketone
			       0.0,   //group aldehyde  
			       0.0,0.0, //group ester
			       0.0,0.0,0.0, //group ether 
			       0.0,  //group acid
			       0.0,   //group ACNO2
			       0.0,0.0,0.0, //group NO3
			       0.0,2.0,0.0, //group CO-OH
			       0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0, //group CO-OC
			       0.0,  //group PAN
			       0.0,  //group CO-OOH
			       0.0,  //group O=COC=O
			       0.0,0.0,0.0}; //group CHxNO2

  size = sizeof(group_tmp_biper)/sizeof(double);
  assert(size == 60);
	
  for(int i = 0; i < size; ++i)
    BiPER.groups[i] = group_tmp_biper[i];

  // Search the species name in the aerosol species list 
  // and add the species if its name matches with
  // the given list.
  add_species_ssh(surrogate, BiPER, species_list_aer, molecular_weight_aer,
		  accomodation_coefficient,diffusion_coef,species_part,nlayer,i_hydrophilic,
		  N_inert, N_inorganic);

  species BiDER;
  BiDER.name="BiDER";
  BiDER.is_inorganic_precursor=false;
  BiDER.Psat_ref=4.10e-7; // Saturation vapor pressure at Tref (torr)
  BiDER.Tref=298;         // Temperature of reference (K)
  BiDER.deltaH=38.4;     // Enthalpy of vaporization (kJ/mol)
  BiDER.Henry=89.1e9;     // Henry's law constant at Tref (M/atm)
  BiDER.aq_type="none"; // "none","diacid","monoacid" or "aldehyde"
  BiDER.hydrophilic=true;   // Does the species condense on the aqueous phase?
  BiDER.hydrophobic=false;  // Does the species condense on the organic phase?
  BiDER.nonvolatile=false; // Is the compound nonvolatile?
  BiDER.kp_from_experiment=false;  // Use experimental partitioning constant at Tref?
  BiDER.is_organic=true;  // Is the compound organic?
  BiDER.compute_gamma_org=true;  // Compute the activity coefficients of the organic phase for this compound?
  BiDER.compute_gamma_aq=true;  // Compute the activity coefficients of the aqueous phase for this compound?
  BiDER.Koligo_org=0.0;       //oligomeriation constant in the organic phase
  BiDER.rho=1300.0;
  BiDER.is_monomer=false;
  BiDER.rion=false;
  //BiDER.KDiffusion_air=1.0e-5;
  //  BiDER.accomodation_coefficient=alpha;
  BiDER.viscosity=1.68e12;
  BiDER.is_solid=false;
  BiDER.is_generic=false;
  
  //Group: if no functionnal group in the species use the default species
  //for the computation of activity coefficients
  //
  double group_tmp_bider [] = {0.0,0.0,0.0,0.0, // group C
			       0.0,2.0,1.0,1.0, //group C[OH]
			       1.0,0.0,0.0,0.0, //group Calcohol
			       0.0,0.0,0.0,0.0, //group Calcohol-tail
			       0.0,0.0,0.0,0.0,0.0, //group C=C
			       0.0,0.0, //group aromatic carbon (AC)
			       0.0,0.0,0.0, // group //AC-C
			       4.0,  //group OH
			       0.0, //group H2O
			       0.0, //group ACOH
			       0.0,0.0, //group ketone
			       0.0,   //group aldehyde  
			       0.0,0.0, //group ester
			       0.0,0.0,0.0, //group ether 
			       0.0,  //group acid
			       0.0,   //group ACNO2
			       0.0,0.0,0.0, //group NO3
			       0.0,0.0,0.0, //group CO-OH
  			       0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0, //group CO-OC
			       0.0,  //group PAN
			       0.0,  //group CO-OOH
			       0.0,  //group O=COC=O
			       0.0,0.0,0.0}; //group CHxNO2

  size = sizeof(group_tmp_bider)/sizeof(double);
  assert(size == 60);
	
  for(int i = 0; i < size; ++i)
    BiDER.groups[i] = group_tmp_bider[i];

  // Search the species name in the aerosol species list 
  // and add the species if its name matches with
  // the given list.
  add_species_ssh(surrogate, BiDER, species_list_aer, molecular_weight_aer,
              accomodation_coefficient,diffusion_coef,species_part,nlayer,i_hydrophilic,
	      N_inert, N_inorganic);

  species BiMGA;
  BiMGA.name="BiMGA";
  BiMGA.is_inorganic_precursor=false;
  BiMGA.Psat_ref=1.40e-5; // Saturation vapor pressure at Tref (torr)
  BiMGA.Tref=298;         // Temperature of reference (K)
  BiMGA.deltaH=43.2;     // Enthalpy of vaporization (kJ/mol)
  BiMGA.Henry=5.25e8;     // Henry's law constant at Tref (M/atm)
  BiMGA.aq_type="monoacid"; // "none","diacid","monoacid" or "aldehyde"
  BiMGA.Kacidity1=1.0e-4;    // First acidity constant
  BiMGA.Koligo_org=64.2;    //oligomeriation constant in the organic phase
  BiMGA.hydrophilic=true;   // Does the species condense on the aqueous phase?
  BiMGA.hydrophobic=false;  // Does the species condense on the organic phase?
  BiMGA.nonvolatile=false; // Is the compound nonvolatile?
  BiMGA.kp_from_experiment=false;  // Use experimental partitioning constant at Tref?
  BiMGA.is_organic=true;  // Is the compound organic?
  BiMGA.compute_gamma_org=true;  // Compute the activity coefficients of the organic phase for this compound?
  BiMGA.compute_gamma_aq=true;  // Compute the activity coefficients of the aqueous phase for this compound?
  BiMGA.rho=1300.0;
  BiMGA.is_monomer=false;
  BiMGA.rion=false;
  //BiMGA.KDiffusion_air=1.0e-5;
  //  BiMGA.accomodation_coefficient=alpha;
  BiMGA.viscosity=1.68e12;
  BiMGA.is_solid=false;
  BiMGA.is_generic=false;
  
  //Group: if no functionnal group in the species use the default species
  //for the computation of activity coefficients
  //
  double group_tmp_bimga [] = {0.0,0.0,0.0,0.0, // group C
			       0.0,0.0,1.0,1.0, //group C[OH]
			       1.0,0.0,0.0,0.0, //group Calcohol
			       0.0,0.0,0.0,0.0, //group Calcohol-tail
			       0.0,0.0,0.0,0.0,0.0, //group C=C
			       0.0,0.0, //group aromatic carbon (AC)
			       0.0,0.0,0.0, // group //AC-C
			       2.0,  //group OH
			       0.0, //group H2O
			       0.0, //group ACOH
			       0.0,0.0, //group ketone
			       0.0,   //group aldehyde  
			       0.0,0.0, //group ester
			       0.0,0.0,0.0, //group ether 
			       1.0,  //group acid
			       0.0,   //group ACNO2
			       0.0,0.0,0.0, //group NO3
			       0.0,0.0,0.0, //group CO-OH
			       0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0, //group CO-OC
			       0.0,  //group PAN
			       0.0,  //group CO-OOH
			       0.0,  //group O=COC=O
			       0.0,0.0,0.0}; //group CHxNO2

  size = sizeof(group_tmp_bimga)/sizeof(double);
  assert(size == 60);
	
  for(int i = 0; i < size; ++i)
    BiMGA.groups[i] = group_tmp_bimga[i];

  // Search the species name in the aerosol species list 
  // and add the species if its name matches with
  // the given list.
  add_species_ssh(surrogate, BiMGA, species_list_aer, molecular_weight_aer,
              accomodation_coefficient,diffusion_coef,species_part,nlayer,i_hydrophilic,
	      N_inert, N_inorganic);

  species AnBlP;
  AnBlP.name="AnBlP";
  AnBlP.is_inorganic_precursor=false;
  AnBlP.Psat_ref=6.8e-8; // Saturation vapor pressure at Tref (torr)
  AnBlP.Tref=298;         // Temperature of reference (K)
  AnBlP.deltaH=50.0;     // Enthalpy of vaporization (kJ/mol)
  AnBlP.aq_type="none";  // "none","diacid","monoacid" or "aldehyde"
  AnBlP.Henry=0.0;       //If the Henry's law constant is set to zero, the model compute the Henry's law constant from the saturation vapour pressure and the activity coefficients at infinite dilution given by unifac
  AnBlP.hydrophilic=false;  // Does the species condense on the aqueous phase?
  AnBlP.hydrophobic=true;  // Does the species condense on the organic phase?
  AnBlP.nonvolatile=false; // Is the compound nonvolatile?
  AnBlP.kp_from_experiment=false;  // Use experimental partitioning constant at Tref?
  AnBlP.is_organic=true;  // Is the compound organic?
  AnBlP.compute_gamma_org=true;  // Compute the activity coefficients of the organic phase for this compound?
  AnBlP.compute_gamma_aq=true;  // Compute the activity coefficients of the aqueous phase for this compound?
  AnBlP.Koligo_org=0.0;      //oligomeriation constant in the organic phase
  AnBlP.rho=1300.0;
  AnBlP.is_monomer=false;
  AnBlP.rion=false;
  //AnBlP.KDiffusion_air=1.0e-5;
  //  AnBlP.accomodation_coefficient=alpha;
  AnBlP.viscosity=1.68e12;
  AnBlP.is_solid=false;
  AnBlP.is_generic=false;
  
  //Group: if no functionnal group in the species use the default species
  //for the computation of activity coefficients
  //
  double group_tmp_anblp [] = {0.0,0.0,0.0,0.0, // group C
			       0.0,0.0,0.0,0.0, //group C[OH]
			       0.0,0.0,0.0,0.0, //group Calcohol
			       0.0,0.0,0.0,0.0, //group Calcohol-tail
			       0.0,0.0,0.0,0.0,0.0, //group C=C
			       3.0,1.0, //group aromatic carbon (AC)
			       1.0,0.0,0.0, // group //AC-C
			       0.0,  //group OH
			       0.0, //group H2O
			       0.0, //group ACOH
			       0.0,0.0, //group ketone
			       0.0,   //group aldehyde  
			       0.0,0.0, //group ester
			       0.0,0.0,0.0, //group ether 
			       1.0,  //group acid
			       1.0,   //group ACNO2
			       0.0,0.0,0.0, //group NO3
			       0.0,0.0,0.0, //group CO-OH
  			       0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0, //group CO-OC
			       0.0,  //group PAN
			       0.0,  //group CO-OOH
			       0.0,  //group O=COC=O
			       0.0,0.0,0.0}; //group CHxNO2

  size = sizeof(group_tmp_anblp)/sizeof(double);
  assert(size == 60);
	
  for(int i = 0; i < size; ++i)
    AnBlP.groups[i] = group_tmp_anblp[i];

  // Search the species name in the aerosol species list 
  // and add the species if its name matches with
  // the given list.
  add_species_ssh(surrogate, AnBlP, species_list_aer, molecular_weight_aer,
              accomodation_coefficient,diffusion_coef,species_part,nlayer,i_hydrophilic,
	      N_inert, N_inorganic);

  species AnBmP;
  AnBmP.name="AnBmP";
  AnBmP.is_inorganic_precursor=false;
  AnBmP.Psat_ref=8.4e-6; // Saturation vapor pressure at Tref (torr)
  AnBmP.Tref=298;         // Temperature of reference (K)
  AnBmP.deltaH=50.0;     // Enthalpy of vaporization (kJ/mol)
  AnBmP.aq_type="none";   // "none","diacid","monoacid" or "aldehyde"
  AnBmP.Henry=0.0;       //If the Henry's law constant is set to zero, the model compute the Henry's law constant from the saturation vapour pressure and the activity coefficients at infinite dilution given by unifac
  AnBmP.hydrophilic=false;   // Does the species condense on the aqueous phase?
  AnBmP.hydrophobic=true;  // Does the species condense on the organic phase?
  AnBmP.nonvolatile=false; // Is the compound nonvolatile?
  AnBmP.kp_from_experiment=false;  // Use experimental partitioning constant at Tref?
  AnBmP.is_organic=true;  // Is the compound organic?
  AnBmP.compute_gamma_org=true;  // Compute the activity coefficients of the organic phase for this compound?
  AnBmP.compute_gamma_aq=true;  // Compute the activity coefficients of the aqueous phase for this compound?
  AnBmP.Koligo_org=0.0;       //oligomeriation constant in the organic phase
  AnBmP.rho=1300.0;
  AnBmP.is_monomer=false;
  AnBmP.rion=false;
  //AnBmP.KDiffusion_air=1.0e-5;
  //  AnBmP.accomodation_coefficient=alpha;
  AnBmP.viscosity=1.68e12;
  AnBmP.is_solid=false;
  AnBmP.is_generic=false;
  
  //Group: if no functionnal group in the species use the default species
  //for the computation of activity coefficients
  //
  double group_tmp_anbmp [] = {0.0,0.0,0.0,0.0, // group C
			       0.0,0.0,0.0,0.0, //group C[OH]
			       0.0,0.0,0.0,0.0, //group Calcohol
			       0.0,0.0,0.0,0.0, //group Calcohol-tail
			       0.0,0.0,0.0,0.0,0.0, //group C=C
			       3.0,1.0, //group aromatic carbon (AC)
			       1.0,0.0,0.0, // group //AC-C
			       0.0,  //group OH
			       0.0, //group H2O
			       1.0, //group ACOH
			       0.0,0.0, //group ketone
			       0.0,   //group aldehyde  
			       0.0,0.0, //group ester
			       0.0,0.0,0.0, //group ether 
			       1.0,  //group acid
			       0.0,   //group ACNO2
			       0.0,0.0,0.0, //group NO3
			       0.0,0.0,0.0, //group CO-OH
			       0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0, //group CO-OC
			       0.0,  //group PAN
			       0.0,  //group CO-OOH
			       0.0,  //group O=COC=O
			       0.0,0.0,0.0}; //group CHxNO2

  size = sizeof(group_tmp_anbmp)/sizeof(double);
  assert(size == 60);
	
  for(int i = 0; i < size; ++i)
    AnBmP.groups[i] = group_tmp_anbmp[i];

  // Search the species name in the aerosol species list 
  // and add the species if its name matches with
  // the given list.
  add_species_ssh(surrogate, AnBmP, species_list_aer, molecular_weight_aer,
              accomodation_coefficient,diffusion_coef,species_part,nlayer,i_hydrophilic,
	      N_inert, N_inorganic);
  
  species BiBlP;
  BiBlP.name="BiBlP";
  BiBlP.is_inorganic_precursor=false;
  BiBlP.Psat_ref=0.60e-9; // Saturation vapor pressure at Tref (torr)
  BiBlP.Tref=298;         // Temperature of reference (K)
  BiBlP.deltaH=175.0;     // Enthalpy of vaporization (kJ/mol)
  BiBlP.aq_type="none";   // "none","diacid","monoacid" or "aldehyde"
  BiBlP.Henry=0.0;        //If the Henry's law constant is set to zero, the model compute the Henry's law constant from the saturation vapour pressure and the activity coefficients at infinite dilution given by unifac
  BiBlP.hydrophilic=false;   // Does the species condense on the aqueous phase?
  BiBlP.hydrophobic=true;  // Does the species condense on the organic phase?
  BiBlP.nonvolatile=false; // Is the compound nonvolatile?
  BiBlP.kp_from_experiment=false;  // Use experimental partitioning constant at Tref?
  BiBlP.is_organic=true;  // Is the compound organic?
  BiBlP.compute_gamma_org=true;  // Compute the activity coefficients of the organic phase for this compound?
  BiBlP.compute_gamma_aq=true;  // Compute the activity coefficients of the aqueous phase for this compound?
  BiBlP.Koligo_org=0.0;    //oligomeriation constant in the organic phase
  BiBlP.rho=1300.0;
  BiBlP.is_monomer=false;
  BiBlP.rion=false;
  //BiBlP.KDiffusion_air=1.0e-5;
  //  BiBlP.accomodation_coefficient=alpha;
  BiBlP.viscosity=1.68e12;  
  BiBlP.is_solid=false;
  BiBlP.is_generic=false;  
  
  //Group: if no functionnal group in the species use the default species
  //for the computation of activity coefficients
  //
  double group_tmp_biblp [] = {3.0,4.0,0.0,1.0, // group C
			       0.0,0.0,1.0,0.0, //group C[OH]
			       0.0,0.0,0.0,0.0, //group Calcohol
			       0.0,0.0,0.0,0.0, //group Calcohol-tail
			       0.0,1.0,0.0,1.0,0.0, //group C=C
			       0.0,0.0, //group aromatic carbon (AC)
			       0.0,0.0,0.0, // group //AC-C
			       1.0,  //group OH
			       0.0, //group H2O
			       0.0, //group ACOH
			       0.0,0.0, //group ketone
			       1.0,   //group aldehyde  
			       0.0,0.0, //group ester
			       0.0,0.0,0.0, //group ether 
			       0.0,  //group acid
			       0.0,   //group ACNO2
			       0.0,1.0,0.0, //group NO3
			       0.0,0.0,0.0, //group CO-OH
			       0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0, //group CO-OC
			       0.0,  //group PAN
			       0.0,  //group CO-OOH
			       0.0,  //group O=COC=O
			       0.0,0.0,0.0}; //group CHxNO2

  size = sizeof(group_tmp_biblp)/sizeof(double);
  assert(size == 60);
	
  for(int i = 0; i < size; ++i)
    BiBlP.groups[i] = group_tmp_biblp[i];

  // Search the species name in the aerosol species list 
  // and add the species if its name matches with
  // the given list.
  add_species_ssh(surrogate, BiBlP, species_list_aer, molecular_weight_aer,
              accomodation_coefficient,diffusion_coef,species_part,nlayer,i_hydrophilic,
	      N_inert, N_inorganic);
  
  species BiBmP;
  BiBmP.name="BiBmP";
  BiBmP.is_inorganic_precursor=false;
  BiBmP.Psat_ref=3.0e-7; // Saturation vapor pressure at Tref (torr)
  BiBmP.Tref=298;         // Temperature of reference (K)
  BiBmP.deltaH=175.0;     // Enthalpy of vaporization (kJ/mol)
  BiBmP.aq_type="none";   // "none","diacid","monoacid" or "aldehyde"
  BiBmP.Henry=0.0;      //If the Henry's law constant is set to zero, the model compute the Henry's law constant from the saturation vapour pressure and the activity coefficients at infinite dilution given by unifac
  BiBmP.hydrophilic=false;   // Does the species condense on the aqueous phase?
  BiBmP.hydrophobic=true;  // Does the species condense on the organic phase?
  BiBmP.nonvolatile=false; // Is the compound nonvolatile?
  BiBmP.kp_from_experiment=false;  // Use experimental partitioning constant at Tref?
  BiBmP.is_organic=true;  // Is the compound organic?
  BiBmP.compute_gamma_org=true;  // Compute the activity coefficients of the organic phase for this compound?
  BiBmP.compute_gamma_aq=true;  // Compute the activity coefficients of the aqueous phase for this compound?
  BiBmP.Koligo_org=0.0;     //oligomeriation constant in the organic phase
  BiBmP.rho=1300.0;
  BiBmP.is_monomer=false;
  BiBmP.rion=false;
  //BiBmP.KDiffusion_air=1.0e-5;
  //  BiBmP.accomodation_coefficient=alpha;
  BiBmP.viscosity=1.68e12;  
  BiBmP.is_solid=false;
  BiBmP.is_generic=false;
  
  //Group: if no functionnal group in the species use the default species
  //for the computation of activity coefficients
  //
  double group_tmp_bibmp [] = {3.0,4.0,0.0,1.0, // group C
			       0.0,0.0,0.0,0.0, //group C[OH]
			       0.0,0.0,0.0,0.0, //group Calcohol
			       0.0,0.0,0.0,0.0, //group Calcohol-tail
			       0.0,1.0,0.0,1.0,0.0, //group C=C
			       0.0,0.0, //group aromatic carbon (AC)
			       0.0,0.0,0.0, // group //AC-C
			       0.0,  //group OH
			       0.0, //group H2O
			       0.0, //group ACOH
			       0.0,1.0, //group ketone
			       1.0,   //group aldehyde  
			       0.0,0.0, //group ester
			       0.0,0.0,0.0, //group ether 
			       0.0,  //group acid
			       0.0,   //group ACNO2
			       0.0,0.0,0.0, //group NO3
			       0.0,0.0,0.0, //group CO-OH
			       0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0, //group CO-OC
			       0.0,  //group PAN
			       0.0,  //group CO-OOH
			       0.0,  //group O=COC=O
			       0.0,0.0,0.0}; //group CHxNO2

  size = sizeof(group_tmp_bibmp)/sizeof(double);
  assert(size == 60);
	
  for(int i = 0; i < size; ++i)
    BiBmP.groups[i] = group_tmp_bibmp[i];

  // Search the species name in the aerosol species list 
  // and add the species if its name matches with
  // the given list.
  add_species_ssh(surrogate, BiBmP, species_list_aer, molecular_weight_aer,
              accomodation_coefficient,diffusion_coef,species_part,nlayer,i_hydrophilic,
	      N_inert, N_inorganic);
  
  species AnClP;
  AnClP.name="AnClP";
  AnClP.is_inorganic_precursor=false;
  AnClP.nonvolatile=true;  // Is the compound nonvolatile?
  AnClP.hydrophilic=false; // Does the species condense on the aqueous phase?
  AnClP.hydrophobic=true;  // Does the species condense on the organic phase?
  AnClP.aq_type="none";  // "none","diacid","monoacid" or "aldehyde"
  //AnClP.Henry=0.0;
  AnClP.kp_from_experiment=false;  // Use experimental partitioning constant at Tref?
  AnClP.is_organic=true;  // Is the compound organic?
  AnClP.compute_gamma_org=true;  // Compute the activity coefficients of the organic phase for this compound?
  AnClP.compute_gamma_aq=true;  // Compute the activity coefficients of the aqueous phase for this compound?
  AnClP.rho=1300.0;
  AnClP.is_monomer=false;
  AnClP.rion=false;
  //AnClP.KDiffusion_air=1.0e-5;
  //  AnClP.accomodation_coefficient=alpha;
  AnClP.viscosity=1.68e12;  
  AnClP.is_solid=false;
  AnClP.is_generic=false;

  //Group: if no functionnal group in the species use the default species
  //for the computation of activity coefficients 
  double group_tmp_anclp [] = {0.0,0.0,0.0,0.0, // group C
			       0.0,0.0,0.0,0.0, //group C[OH]
			       0.0,0.0,0.0,0.0, //group Calcohol
			       0.0,0.0,0.0,0.0, //group Calcohol-tail
			       0.0,0.0,0.0,0.0,0.0, //group C=C
			       0.0,0.0, //group aromatic carbon (AC)
			       0.0,0.0,0.0, // group //AC-C
			       0.0,  //group OH
			       0.0, //group H2O
			       0.0, //group ACOH
			       0.0,0.0, //group ketone
			       0.0,   //group aldehyde  
			       0.0,0.0, //group ester
			       0.0,0.0,0.0, //group ether
			       0.0,  //group acid
			       0.0,   //group ACNO2
			       0.0,0.0,0.0, //group NO3
			       0.0,0.0,0.0, //group CO-OH
			       0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0, //group CO-OC
			       0.0,  //group PAN
			       0.0,  //group CO-OOH
			       0.0,  //group O=COC=O
			       0.0,0.0,0.0}; //group CHxNO2
  
  size = sizeof(group_tmp_anclp)/sizeof(double);
  assert(size == 60);
  
  for(int i = 0; i < size; ++i)
    AnClP.groups[i] = group_tmp_anclp[i];

  // Search the species name in the aerosol species list 
  // and add the species if its name matches with
  // the given list.
  add_species_ssh(surrogate, AnClP, species_list_aer, molecular_weight_aer,
              accomodation_coefficient,diffusion_coef,species_part,nlayer,i_hydrophilic,
	      N_inert, N_inorganic);

  species BiNGA;
  BiNGA.name="BiNGA";
  BiNGA.is_inorganic_precursor=false;
  BiNGA.Psat_ref=1.40e-5; // Saturation vapor pressure at Tref (torr)
  BiNGA.Tref=298;         // Temperature of reference (K)
  BiNGA.deltaH=43.2;     // Enthalpy of vaporization (kJ/mol)
  BiNGA.Henry=3.73e7;     // Henry's law constant at Tref (M/atm)
  BiNGA.aq_type="monoacid"; // "none","diacid","monoacid" or "aldehyde"
  BiNGA.Kacidity1=1.0e-4;    // First acidity constant
  BiNGA.Koligo_org=64.2;    //oligomeriation constant in the organic phase
  BiNGA.hydrophilic=false;   // Does the species condense on the aqueous phase?
  BiNGA.hydrophobic=true;  // Does the species condense on the organic phase?
  BiNGA.nonvolatile=false; // Is the compound nonvolatile?
  BiNGA.kp_from_experiment=false;  // Use experimental partitioning constant at Tref?
  BiNGA.is_organic=true;  // Is the compound organic?
  BiNGA.compute_gamma_org=true;  // Compute the activity coefficients of the organic phase for this compound?
  BiNGA.compute_gamma_aq=true;  // Compute the activity coefficients of the aqueous phase for this compound?
  BiNGA.rho=1300.0;
  BiNGA.is_monomer=false;
  BiNGA.rion=false;
  // BiNGA.KDiffusion_air=1.0e-5;
  //  BiNGA.accomodation_coefficient=alpha;
  BiNGA.viscosity=1.68e12;
  BiNGA.is_solid=false;
  BiNGA.is_generic=false;
  
  //Group: if no functionnal group in the species use the default species
  //for the computation of activity coefficients
  //
  double group_tmp_binga [] = {0.0,0.0,0.0,0.0, // group C
			       0.0,0.0,0.0,1.0, //group C[OH]
			       1.0,0.0,1.0,0.0, //group Calcohol
			       0.0,0.0,0.0,0.0, //group Calcohol-tail
			       0.0,0.0,0.0,0.0,0.0, //group C=C
			       0.0,0.0, //group aromatic carbon (AC)
			       0.0,0.0,0.0, // group //AC-C
			       1.0,  //group OH
			       0.0, //group H2O
			       0.0, //group ACOH
			       0.0,0.0, //group ketone
			       0.0,   //group aldehyde  
			       0.0,0.0, //group ester
			       0.0,0.0,0.0, //group ether 
			       1.0,  //group acid
			       0.0,   //group ACNO2
			       0.0,0.0,1.0, //group NO3
			       0.0,0.0,0.0, //group CO-OH
			       0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0, //group CO-OC
			       0.0,  //group PAN
			       0.0,  //group CO-OOH
			       0.0,  //group O=COC=O
			       0.0,0.0,0.0}; //group CHxNO2

  size = sizeof(group_tmp_binga)/sizeof(double);
  assert(size == 60);
	
  for(int i = 0; i < size; ++i)
    BiNGA.groups[i] = group_tmp_binga[i];

  // Search the species name in the aerosol species list 
  // and add the species if its name matches with
  // the given list.
  add_species_ssh(surrogate, BiNGA, species_list_aer, molecular_weight_aer,
              accomodation_coefficient,diffusion_coef,species_part,nlayer,i_hydrophilic,
	      N_inert, N_inorganic);

  species BiNIT3;
  BiNIT3.name="BiNIT3"; 
  BiNIT3.is_inorganic_precursor=false;  
  BiNIT3.Psat_ref=1.45e-6; // Saturation vapor pressure at Tref (torr)
  BiNIT3.Tref=298;         // Temperature of reference (K)
  BiNIT3.deltaH=38.4;     // Enthalpy of vaporization (kJ/mol)
  BiNIT3.aq_type="none";  // "none","diacid","monoacid" or "aldehyde"
  BiNIT3.Henry=0.0;       //If the Henry's law constant is set to zero, the model compute the Henry's law constant from the saturation vapour pressure and the activity coefficients at infinite dilution given by unifac
  BiNIT3.hydrophilic=false;   // Does the species condense on the aqueous phase?
  BiNIT3.hydrophobic=true;  // Does the species condense on the organic phase?
  BiNIT3.nonvolatile=false; // Is the compound nonvolatile?
  BiNIT3.kp_from_experiment=false;  // Use experimental partitioning constant at Tref?
  BiNIT3.is_organic=true;  // Is the compound organic?
  BiNIT3.compute_gamma_org=true;  // Compute the activity coefficients of the organic phase for this compound?
  BiNIT3.compute_gamma_aq=true;  // Compute the activity coefficients of the aqueous phase for this compound?
  BiNIT3.Koligo_org=0.0;      //oligomeriation constant in the organic phase
  BiNIT3.rho=1300.0;
  BiNIT3.is_monomer=false;
  BiNIT3.rion=false;
  // BiNIT3.KDiffusion_air=1.0e-5;
  //  BiNIT3.accomodation_coefficient=alpha;
  BiNIT3.viscosity=1.68e12;  
  BiNIT3.is_solid=false;
  BiNIT3.is_generic=false;
  
  //Group: if no functionnal group in the species use the default species
  //for the computation of activity coefficients
  //
  double group_tmp_binit3 [] = {0.0,0.0,0.0,0.0, // group C
				0.0,1.0,0.0,0.0, //group C[OH]
				1.0,0.0,0.0,0.0, //group Calcohol
				0.0,0.0,0.0,0.0, //group Calcohol-tail
				0.0,0.0,0.0,0.0,0.0, //group C=C
				0.0,0.0, //group aromatic carbon (AC)
				0.0,0.0,0.0, // group //AC-C
				1.0,  //group OH
				0.0, //group H2O
				0.0, //group ACOH
				0.0,0.0, //group ketone
				0.0,   //group aldehyde  
				0.0,0.0, //group ester
				0.0,0.0,0.0, //group ether 
				0.0,  //group acid
				0.0,   //group ACNO2
				1.0,2.0,0.0, //group NO3
				0.0,0.0,0.0, //group CO-OH
				0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0, //group CO-OC
				0.0,  //group PAN
				0.0,  //group CO-OOH
			       0.0,  //group O=COC=O
			       0.0,0.0,0.0}; //group CHxNO2

  size = sizeof(group_tmp_binit3)/sizeof(double);
  assert(size == 60);
	
  for(int i = 0; i < size; ++i)
    BiNIT3.groups[i] = group_tmp_binit3[i];

  // Search the species name in the aerosol species list 
  // and add the species if its name matches with
  // the given list.
  add_species_ssh(surrogate, BiNIT3, species_list_aer, molecular_weight_aer,
		  accomodation_coefficient,diffusion_coef,species_part,nlayer,i_hydrophilic,
		  N_inert, N_inorganic);
  
  species BiNIT;
  BiNIT.name="BiNIT";
  BiNIT.is_inorganic_precursor=false;
  BiNIT.Psat_ref=2.5e-6; // Saturation vapor pressure at Tref (torr)
  BiNIT.Tref=298;         // Temperature of reference (K)
  BiNIT.deltaH=50.0;     // Enthalpy of vaporization (kJ/mol)
  BiNIT.aq_type="none";  // "none","diacid","monoacid" or "aldehyde"
  BiNIT.Henry=0.0;       //If the Henry's law constant is set to zero, the model compute the Henry's law constant from the saturation vapour pressure and the activity coefficients at infinite dilution given by unifac
  BiNIT.hydrophilic=false;   // Does the species condense on the aqueous phase?
  BiNIT.hydrophobic=true;  // Does the species condense on the organic phase?
  BiNIT.nonvolatile=false; // Is the compound nonvolatile?
  BiNIT.kp_from_experiment=false;  // Use experimental partitioning constant at Tref?
  BiNIT.is_organic=true;  // Is the compound organic?
  BiNIT.compute_gamma_org=true;  // Compute the activity coefficients of the organic phase for this compound?
  BiNIT.compute_gamma_aq=true;  // Compute the activity coefficients of the aqueous phase for this compound?
  BiNIT.Koligo_org=0.0;
  BiNIT.rho=1300.0;
  BiNIT.is_monomer=false;
  BiNIT.rion=false;
  // BiNIT.KDiffusion_air=1.0e-5;
  //  BiNIT.accomodation_coefficient=alpha;
  BiNIT.viscosity=1.68e12;  
  BiNIT.is_solid=false;
  BiNIT.is_generic=false;
  
  //Group: if no functionnal group in the species use the default species
  //for the computation of activity coefficients
  //
  double group_tmp_binit [] = {3.0,1.0,2.0,2.0, // group C
			       0.0,0.0,1.0,0.0, //group C[OH]
			       0.0,0.0,0.0,0.0, //group Calcohol
			       0.0,0.0,0.0,0.0, //group Calcohol-tail   
			       0.0,0.0,0.0,0.0,0.0, //group C=C
			       0.0,0.0, //group aromatic carbon (AC)
			       0.0,0.0,0.0, // group //AC-C
			       1.0,  //group OH
			       0.0, //group H2O
			       0.0, //group ACOH
			       0.0,0.0, //group ketone
			       0.0,   //group aldehyde  
			       0.0,0.0, //group ester
			       0.0,0.0,0.0, //group ether 
			       0.0,  //group acid
			       0.0,   //group ACNO2
			       0.0,1.0,0.0, //group NO3
			       0.0,0.0,0.0, //group CO-OH
			       0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0, //group CO-OC
			       0.0,  //group PAN
			       0.0,  //group CO-OOH
			       0.0,  //group O=COC=O
			       0.0,0.0,0.0}; //group CHxNO2

  size = sizeof(group_tmp_binit)/sizeof(double);
  assert(size == 60);
	
  for(int i = 0; i < size; ++i)
    BiNIT.groups[i] = group_tmp_binit[i];

  // Search the species name in the aerosol species list 
  // and add the species if its name matches with
  // the given list.
  add_species_ssh(surrogate, BiNIT, species_list_aer, molecular_weight_aer,
		  accomodation_coefficient,diffusion_coef,species_part,nlayer,i_hydrophilic,
		  N_inert, N_inorganic);
  
  species POAlP;
  POAlP.name="POAlP";
  POAlP.is_inorganic_precursor=false;
  POAlP.nonvolatile=false;  // Is the compound nonvolatile?
  POAlP.hydrophilic=false; // Does the species condense on the aqueous phase?
  POAlP.hydrophobic=true;  // Does the species condense on the organic phase?
  POAlP.kp_from_experiment=true;  // Use experimental partitioning constant at Tref?
  //  POAlP.kp_experiment=0.0435;       
  //   POAlP.kp_experiment=0.0015;       
  POAlP.kp_experiment=1.1;       // Value of the experimental partitioning constant at Tref?
  POAlP.deltaH=106.0;     // Enthalpy of vaporization (kJ/mol)
  POAlP.Tref=298;
  POAlP.is_organic=true;  // Is the compound organic?
  POAlP.compute_gamma_org=true;  // Compute the activity coefficients of the organic phase for this compound?
  POAlP.compute_gamma_aq=true;  // Compute the activity coefficients of the aqueous phase for this compound
  POAlP.rho=1300.0;
  POAlP.is_monomer=false;
  POAlP.rion=false;
  // POAlP.KDiffusion_air=1.0e-5;
  //  POAlP.accomodation_coefficient=alpha;
  POAlP.viscosity=1.68e12;  
  POAlP.is_solid=false;
  POAlP.is_generic=false;

  //Group: if no functionnal group in the species use the default species
  //for the computation of activity coefficients 
  
  double group_tmp_poalp [] = {0.0,0.0,0.0,0.0, // group C
			       0.0,0.0,0.0,0.0, //group C[OH]
			       0.0,0.0,0.0,0.0, //group Calcohol
			       0.0,0.0,0.0,0.0, //group Calcohol-tail  
			       0.0,0.0,0.0,0.0,0.0, //group C=C
			       0.0,0.0, //group aromatic carbon (AC)
			       0.0,0.0,0.0, // group //AC-C
			       0.0,  //group OH
			       0.0, //group H2O
			       0.0, //group ACOH
			       0.0,0.0, //group ketone
			       0.0,   //group aldehyde  
			       0.0,0.0, //group ester
			       0.0,0.0,0.0, //group ether  
			       0.0,  //group acid
			       0.0,   //group ACNO2
			       0.0,0.0,0.0, //group NO3
			       0.0,0.0,0.0, //group CO-OH			       
			       0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0, //group CO-OC
			       0.0,  //group PAN
			       0.0,  //group CO-OOH
			       0.0,  //group O=COC=O
			       0.0,0.0,0.0}; //group CHxNO2
  
  size = sizeof(group_tmp_poalp)/sizeof(double);
  assert(size == 60);
  for(int i = 0; i < size; ++i)
    POAlP.groups[i] = group_tmp_poalp[i];

  // Search the species name in the aerosol species list 
  // and add the species if its name matches with
  // the given list.
  add_species_ssh(surrogate, POAlP, species_list_aer, molecular_weight_aer,
		  accomodation_coefficient,diffusion_coef,species_part,nlayer,i_hydrophilic,
		  N_inert, N_inorganic);
  
  species POAmP;
  POAmP.name="POAmP";
  POAmP.is_inorganic_precursor=false;
  POAmP.nonvolatile=false;  // Is the compound nonvolatile?
  POAmP.hydrophilic=false; // Does the species condense on the aqueous phase?
  POAmP.hydrophobic=true;  // Does the species condense on the organic phase?
  POAmP.kp_from_experiment=true;  // Use experimental partitioning constant at Tref?
  POAmP.kp_experiment=0.011;       // Value of the experimental partitioning constant at Tref?
  POAmP.deltaH=91.0;     // Enthalpy of vaporization (kJ/mol)
  POAmP.Tref=298;
  POAmP.is_organic=true;  // Is the compound organic?
  POAmP.compute_gamma_org=true;  // Compute the activity coefficients of the organic phase for this compound?
  POAmP.compute_gamma_aq=true;  // Compute the activity coefficients of the aqueous phase for this compound
  POAmP.rho=1300.0;
  POAmP.is_monomer=false;
  POAmP.rion=false;
  //POAmP.KDiffusion_air=1.0e-5;
  //  POAmP.accomodation_coefficient=alpha;
  POAmP.viscosity=1.68e12;  
  POAmP.is_solid=false;
  POAmP.is_generic=false;

  //Group: if no functionnal group in the species use the default species
  //for the computation of activity coefficients 
  
  double group_tmp_poamp [] = {0.0,0.0,0.0,0.0, // group C
			       0.0,0.0,0.0,0.0, //group C[OH]
			       0.0,0.0,0.0,0.0, //group Calcohol
			       0.0,0.0,0.0,0.0, //group Calcohol-tail  
			       0.0,0.0,0.0,0.0,0.0, //group C=C
			       0.0,0.0, //group aromatic carbon (AC)
			       0.0,0.0,0.0, // group //AC-C
			       0.0,  //group OH
			       0.0, //group H2O
			       0.0, //group ACOH
			       0.0,0.0, //group ketone
			       0.0,   //group aldehyde  
			       0.0,0.0, //group ester
			       0.0,0.0,0.0, //group ether  
			       0.0,  //group acid
			       0.0,   //group ACNO2
			       0.0,0.0,0.0, //group NO3
			       0.0,0.0,0.0, //group CO-OH
			       0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0, //group CO-OC
			       0.0,  //group PAN
			       0.0,  //group CO-OOH
			       0.0,  //group O=COC=O
			       0.0,0.0,0.0}; //group CHxNO2

  size = sizeof(group_tmp_poamp)/sizeof(double);
  assert(size == 60);
  for(int i = 0; i < size; ++i)
    POAmP.groups[i] = group_tmp_poamp[i];

  // Search the species name in the aerosol species list 
  // and add the species if its name matches with
  // the given list.
  add_species_ssh(surrogate, POAmP, species_list_aer, molecular_weight_aer,
              accomodation_coefficient,diffusion_coef,species_part,nlayer,i_hydrophilic,
	      N_inert, N_inorganic);
  
  species POAhP;
  POAhP.name="POAhP";
  POAhP.is_inorganic_precursor=false;
  POAhP.nonvolatile=false;  // Is the compound nonvolatile?
  POAhP.hydrophilic=false; // Does the species condense on the aqueous phase?
  POAhP.hydrophobic=true;  // Does the species condense on the organic phase?
  POAhP.kp_from_experiment=true;  // Use experimental partitioning constant at Tref?
  POAhP.kp_experiment=0.00031;       // Value of the experimental partitioning constant at Tref?
  POAhP.deltaH=79.0;     // Enthalpy of vaporization (kJ/mol)
  POAhP.Tref=298;
  POAhP.is_organic=true;  // Is the compound organic?
  POAhP.compute_gamma_org=true;  // Compute the activity coefficients of the organic phase for this compound?
  POAhP.compute_gamma_aq=true;  // Compute the activity coefficients of the aqueous phase for this compound
  POAhP.rho=1300.0;
  POAhP.is_monomer=false;
  POAhP.rion=false;
  //POAhP.KDiffusion_air=1.0e-5;
  //  POAhP.accomodation_coefficient=alpha;
  POAhP.viscosity=1.68e12;  
  POAhP.is_solid=false;
  POAhP.is_generic=false;

  //Group: if no functionnal group in the species use the default species
  //for the computation of activity coefficients 
  
  double group_tmp_poahp [] = {0.0,0.0,0.0,0.0, // group C
			       0.0,0.0,0.0,0.0, //group C[OH]
			       0.0,0.0,0.0,0.0, //group Calcohol
			       0.0,0.0,0.0,0.0, //group Calcohol-tail  
			       0.0,0.0,0.0,0.0,0.0, //group C=C
			       0.0,0.0, //group aromatic carbon (AC)
			       0.0,0.0,0.0, // group //AC-C
			       0.0,  //group OH
			       0.0, //group H2O
			       0.0, //group ACOH
			       0.0,0.0, //group ketone
			       0.0,   //group aldehyde  
			       0.0,0.0, //group ester
			       0.0,0.0,0.0, //group ether  
			       0.0,  //group acid
			       0.0,   //group ACNO2
			       0.0,0.0,0.0, //group NO3
			       0.0,0.0,0.0, //group CO-OH
			       0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0, //group CO-OC
			       0.0,  //group PAN
			       0.0,  //group CO-OOH
			       0.0,  //group O=COC=O
			       0.0,0.0,0.0}; //group CHxNO2

  size = sizeof(group_tmp_poahp)/sizeof(double);
  assert(size == 60);
  for(int i = 0; i < size; ++i)
    POAhP.groups[i] = group_tmp_poahp[i];

  // Search the species name in the aerosol species list 
  // and add the species if its name matches with
  // the given list.
  add_species_ssh(surrogate, POAhP, species_list_aer, molecular_weight_aer,
              accomodation_coefficient,diffusion_coef,species_part,nlayer,i_hydrophilic,
	      N_inert, N_inorganic);
  
  species SOAlP;
  SOAlP.name="SOAlP";
  SOAlP.is_inorganic_precursor=false;
  SOAlP.nonvolatile=false;  // Is the compound nonvolatile?
  SOAlP.hydrophilic=false; // Does the species condense on the aqueous phase?
  SOAlP.hydrophobic=true;  // Does the species condense on the organic phase?
  SOAlP.kp_from_experiment=true;  // Use experimental partitioning constant at Tref?
  SOAlP.kp_experiment=110.0;       // Value of the experimental partitioning constant at Tref?
  SOAlP.deltaH=106.0;     // Enthalpy of vaporization (kJ/mol)
  SOAlP.Tref=298;
  SOAlP.is_organic=true;  // Is the compound organic?
  SOAlP.compute_gamma_org=true;  // Compute the activity coefficients of the organic phase for this compound?
  SOAlP.compute_gamma_aq=true;  // Compute the activity coefficients of the aqueous phase for this compound
  SOAlP.rho=1300.0;
  SOAlP.is_monomer=false;
  SOAlP.rion=false;
  //SOAlP.KDiffusion_air=1.0e-5;
  //  SOAlP.accomodation_coefficient=alpha;
  SOAlP.viscosity=1.68e12;  
  SOAlP.is_solid=false;
  SOAlP.is_generic=false;
 
  //Group: if no functionnal group in the species use the default species
  //for the computation of activity coefficients 
  
  double group_tmp_soalp [] = {0.0,0.0,0.0,0.0, // group C
			       0.0,0.0,0.0,0.0, //group C[OH]
			       0.0,0.0,0.0,0.0, //group Calcohol
			       0.0,0.0,0.0,0.0, //group Calcohol-tail  
			       0.0,0.0,0.0,0.0,0.0, //group C=C
			       0.0,0.0, //group aromatic carbon (AC)
			       0.0,0.0,0.0, // group //AC-C
			       0.0,  //group OH
			       0.0, //group H2O
			       0.0, //group ACOH
			       0.0,0.0, //group ketone
			       0.0,   //group aldehyde  
			       0.0,0.0, //group ester
			       0.0,0.0,0.0, //group ether  
			       0.0,  //group acid
			       0.0,   //group ACNO2
			       0.0,0.0,0.0, //group NO3
			       0.0,0.0,0.0, //group CO-OH
			       0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0, //group CO-OC
			       0.0,  //group PAN
			       0.0,  //group CO-OOH
			       0.0,  //group O=COC=O
			       0.0,0.0,0.0}; //group CHxNO2

  size = sizeof(group_tmp_soalp)/sizeof(double);
  assert(size == 60);
  for(int i = 0; i < size; ++i)
    SOAlP.groups[i] = group_tmp_soalp[i];

  // Search the species name in the aerosol species list 
  // and add the species if its name matches with
  // the given list.
  add_species_ssh(surrogate, SOAlP, species_list_aer, molecular_weight_aer,
		  accomodation_coefficient,diffusion_coef,species_part,nlayer,i_hydrophilic,
		  N_inert, N_inorganic);
  
  species SOAmP;
  SOAmP.name="SOAmP";
  SOAmP.is_inorganic_precursor=false;
  SOAmP.nonvolatile=false;  // Is the compound nonvolatile?
  SOAmP.hydrophilic=false; // Does the species condense on the aqueous phase?
  SOAmP.hydrophobic=true;  // Does the species condense on the organic phase?
  SOAmP.kp_from_experiment=true;  // Use experimental partitioning constant at Tref?
  SOAmP.kp_experiment=1.10;       // Value of the experimental partitioning constant at Tref?
  SOAmP.deltaH=91.0;     // Enthalpy of vaporization (kJ/mol)
  SOAmP.Tref=298;
  SOAmP.is_organic=true;  // Is the compound organic?
  SOAmP.compute_gamma_org=true;  // Compute the activity coefficients of the organic phase for this compound?
  SOAmP.compute_gamma_aq=true;  // Compute the activity coefficients of the aqueous phase for this compound
  SOAmP.rho=1300.0;
  SOAmP.is_monomer=false;
  SOAmP.rion=false;
  //SOAmP.KDiffusion_air=1.0e-5;
  //  SOAmP.accomodation_coefficient=alpha;
  SOAmP.viscosity=1.68e12;
  SOAmP.is_solid=false;
  SOAmP.is_generic=false;

  //Group: if no functionnal group in the species use the default species
  //for the computation of activity coefficients 
  
  double group_tmp_soamp [] = {0.0,0.0,0.0,0.0, // group C
			       0.0,0.0,0.0,0.0, //group C[OH]
			       0.0,0.0,0.0,0.0, //group Calcohol
			       0.0,0.0,0.0,0.0, //group Calcohol-tail  
			       0.0,0.0,0.0,0.0,0.0, //group C=C
			       0.0,0.0, //group aromatic carbon (AC)
			       0.0,0.0,0.0, // group //AC-C
			       0.0,  //group OH
			       0.0, //group H2O
			       0.0, //group ACOH
			       0.0,0.0, //group ketone
			       0.0,   //group aldehyde  
			       0.0,0.0, //group ester
			       0.0,0.0,0.0, //group ether  
			       0.0,  //group acid
			       0.0,   //group ACNO2
			       0.0,0.0,0.0, //group NO3
			       0.0,0.0,0.0, //group CO-OH
			       0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0, //group CO-OC
			       0.0,  //group PAN
			       0.0,  //group CO-OOH
			       0.0,  //group O=COC=O
			       0.0,0.0,0.0}; //group CHxNO2

  size = sizeof(group_tmp_soamp)/sizeof(double);
  assert(size == 60);
  for(int i = 0; i < size; ++i)
    SOAmP.groups[i] = group_tmp_soamp[i];

  // Search the species name in the aerosol species list 
  // and add the species if its name matches with
  // the given list.
  add_species_ssh(surrogate, SOAmP, species_list_aer, molecular_weight_aer,
		  accomodation_coefficient,diffusion_coef,species_part,nlayer,i_hydrophilic,
		  N_inert, N_inorganic);
  
  species SOAhP;
  SOAhP.name="SOAhP";
  SOAhP.is_inorganic_precursor=false;
  SOAhP.nonvolatile=false;  // Is the compound nonvolatile?
  SOAhP.hydrophilic=false; // Does the species condense on the aqueous phase?
  SOAhP.hydrophobic=true;  // Does the species condense on the organic phase?
  SOAhP.kp_from_experiment=true;  // Use experimental partitioning constant at Tref?
  SOAhP.kp_experiment=0.031;       // Value of the experimental partitioning constant at Tref?
  SOAhP.deltaH=79.0;     // Enthalpy of vaporization (kJ/mol)
  SOAhP.Tref=298;
  SOAhP.is_organic=true;  // Is the compound organic?
  SOAhP.compute_gamma_org=true;  // Compute the activity coefficients of the organic phase for this compound?
  SOAhP.compute_gamma_aq=true;  // Compute the activity coefficients of the aqueous phase for this compound
  SOAhP.rho=1300.0;
  SOAhP.is_monomer=false;
  SOAhP.rion=false;
  //SOAhP.KDiffusion_air=1.0e-5;
  //  SOAhP.accomodation_coefficient=alpha;
  SOAhP.viscosity=1.68e12;  
  SOAhP.is_solid=false;
  SOAhP.is_generic=false;

  //Group: if no functionnal group in the species use the default species
  //for the computation of activity coefficients 
  
  double group_tmp_soahp [] = {0.0,0.0,0.0,0.0, // group C
			       0.0,0.0,0.0,0.0, //group C[OH]
			       0.0,0.0,0.0,0.0, //group Calcohol
			       0.0,0.0,0.0,0.0, //group Calcohol-tail  
			       0.0,0.0,0.0,0.0,0.0, //group C=C
			       0.0,0.0, //group aromatic carbon (AC)
			       0.0,0.0,0.0, // group //AC-C
			       0.0,  //group OH
			       0.0, //group H2O
			       0.0, //group ACOH
			       0.0,0.0, //group ketone
			       0.0,   //group aldehyde  
			       0.0,0.0, //group ester
			       0.0,0.0,0.0, //group ether  
			       0.0,  //group acid
			       0.0,   //group ACNO2
			       0.0,0.0,0.0, //group NO3
			       0.0,0.0,0.0, //group CO-OH
			       0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0, //group CO-OC
			       0.0,  //group PAN
			       0.0,  //group CO-OOH
			       0.0,  //group O=COC=O
			       0.0,0.0,0.0}; //group CHxNO2

  size = sizeof(group_tmp_soahp)/sizeof(double);
  assert(size == 60);
  for(int i = 0; i < size; ++i)
    SOAhP.groups[i] = group_tmp_soahp[i];

  // Search the species name in the aerosol species list 
  // and add the species if its name matches with
  // the given list.
  add_species_ssh(surrogate, SOAhP, species_list_aer, molecular_weight_aer,
		  accomodation_coefficient,diffusion_coef,species_part,nlayer,i_hydrophilic,
		  N_inert, N_inorganic);
  
  species Monomer;
  Monomer.name="Monomer";
  Monomer.is_inorganic_precursor=false;
  Monomer.Psat_ref=1.0e-14; // Saturation vapor pressure at Tref (torr)
  Monomer.Tref=298;         // Temperature of reference (K)
  Monomer.deltaH=50.0;     // Enthalpy of vaporization (kJ/mol)
  Monomer.aq_type="none";  // "none","diacid","monoacid" or "aldehyde"
  Monomer.Henry=0.0;       //If the Henry's law constant is set to zero, the model compute the Henry's law constant from the saturation vapour pressure and the activity coefficients at infinite dilution given by unifac
  Monomer.hydrophilic=false;  // Does the species condense on the aqueous phase?
  Monomer.hydrophobic=true;  // Does the species condense on the organic phase?
  Monomer.nonvolatile=true; // Is the compound nonvolatile?
  Monomer.kp_from_experiment=false;  // Use experimental partitioning constant at Tref?
  Monomer.is_organic=true;  // Is the compound organic?
  Monomer.compute_gamma_org=false;  // Compute the activity coefficients of the organic phase for this compound?
  Monomer.compute_gamma_aq=false;  // Compute the activity coefficients of the aqueous phase for this compound?
  Monomer.Koligo_org=0.0;      //oligomeriation constant in the organic phase
  Monomer.rho=1840.0;
  Monomer.is_monomer=false;
  Monomer.rion=false;
  //Monomer.KDiffusion_air=1.0e-5;
  //  Monomer.accomodation_coefficient=alpha;
  Monomer.viscosity=1.68e12;
  Monomer.is_solid=false;
  Monomer.is_generic=false;
  
  //Group: if no functionnal group in the species use the default species
  //for the computation of activity coefficients
  //
  double group_tmp_Monomer [] = {0.0,0.0,0.0,0.0, // group C
				 0.0,0.0,0.0,0.0, //group C[OH]
				 0.0,0.0,0.0,0.0, //group Calcohol
				 0.0,0.0,0.0,0.0, //group Calcohol-tail
				 0.0,0.0,0.0,0.0,0.0, //group C=C
				 3.0,1.0, //group aromatic carbon (AC)
				 1.0,0.0,0.0, // group //AC-C
				 0.0,  //group OH
				 0.0, //group H2O
				 0.0, //group ACOH
				 0.0,0.0, //group ketone
				 0.0,   //group aldehyde  
				 0.0,0.0, //group ester
				 0.0,0.0,0.0, //group ether 
				 1.0,  //group acid
				 1.0,   //group ACNO2
				 0.0,0.0,0.0, //group NO3
				 0.0,0.0,0.0, //group CO-OH
				 0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0, //group CO-OC
				 0.0,  //group PAN
				 0.0,  //group CO-OOH
			       0.0,  //group O=COC=O
			       0.0,0.0,0.0}; //group CHxNO2

  size = sizeof(group_tmp_Monomer)/sizeof(double);
  assert(size == 60);
	
  for(int i = 0; i < size; ++i)
    Monomer.groups[i] = group_tmp_Monomer[i];

  // Search the species name in the aerosol species list 
  // and add the species if its name matches with
  // the given list.
  add_species_ssh(surrogate, Monomer, species_list_aer, molecular_weight_aer,
		  accomodation_coefficient,diffusion_coef,species_part,nlayer,i_hydrophilic,
		  N_inert, N_inorganic);
  
  species Dimer;
  Dimer.name="Dimer";
  Dimer.is_inorganic_precursor=false;
  Dimer.Psat_ref=1.0e-14; // Saturation vapor pressure at Tref (torr)
  Dimer.Tref=298;         // Temperature of reference (K)
  Dimer.deltaH=50.0;     // Enthalpy of vaporization (kJ/mol)
  Dimer.aq_type="none";  // "none","diacid","monoacid" or "aldehyde"
  Dimer.Henry=0.0;       //If the Henry's law constant is set to zero, the model compute the Henry's law constant from the saturation vapour pressure and the activity coefficients at infinite dilution given by unifac
  Dimer.hydrophilic=false;  // Does the species condense on the aqueous phase?
  Dimer.hydrophobic=true;  // Does the species condense on the organic phase?
  Dimer.nonvolatile=false; // Is the compound nonvolatile?
  Dimer.kp_from_experiment=false;  // Use experimental partitioning constant at Tref?
  Dimer.is_organic=true;  // Is the compound organic?
  Dimer.compute_gamma_org=false;  // Compute the activity coefficients of the organic phase for this compound?
  Dimer.compute_gamma_aq=false;  // Compute the activity coefficients of the aqueous phase for this compound?
  Dimer.Koligo_org=0.0;      //oligomeriation constant in the organic phase
  Dimer.rho=1300.0;
  Dimer.is_monomer=false;
  Dimer.rion=false;
  // Dimer.KDiffusion_air=1.0e-5;
  //  Dimer.accomodation_coefficient=alpha;
  Dimer.viscosity=1.68e12;
  Dimer.is_solid=false;
  Dimer.is_generic=false;
  
  //Group: if no functionnal group in the species use the default species
  //for the computation of activity coefficients
  //
  double group_tmp_dimer [] = {0.0,0.0,0.0,0.0, // group C
			       0.0,0.0,0.0,0.0, //group C[OH]
			       0.0,0.0,0.0,0.0, //group Calcohol
			       0.0,0.0,0.0,0.0, //group Calcohol-tail
			       0.0,0.0,0.0,0.0,0.0, //group C=C
			       3.0,1.0, //group aromatic carbon (AC)
			       1.0,0.0,0.0, // group //AC-C
			       0.0,  //group OH
			       0.0, //group H2O
			       0.0, //group ACOH
			       0.0,0.0, //group ketone
			       0.0,   //group aldehyde  
			       0.0,0.0, //group ester
			       0.0,0.0,0.0, //group ether 
			       1.0,  //group acid
			       1.0,   //group ACNO2
			       0.0,0.0,0.0, //group NO3
			       0.0,0.0,0.0, //group CO-OH
			       0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0, //group CO-OC
			       0.0,  //group PAN
			       0.0,  //group CO-OOH
			       0.0,  //group O=COC=O
			       0.0,0.0,0.0}; //group CHxNO2

  size = sizeof(group_tmp_dimer)/sizeof(double);
  assert(size == 60);
	
  for(int i = 0; i < size; ++i)
    Dimer.groups[i] = group_tmp_dimer[i];

  // Search the species name in the aerosol species list 
  // and add the species if its name matches with
  // the given list.
  add_species_ssh(surrogate, Dimer, species_list_aer, molecular_weight_aer,
		  accomodation_coefficient,diffusion_coef,species_part,nlayer,i_hydrophilic,
		  N_inert, N_inorganic);

  species BiA3D;
  BiA3D.name="BiA3D";
  BiA3D.is_inorganic_precursor=false;
  BiA3D.Psat_ref=3.25e-7; // Saturation vapor pressure at Tref (torr)
  BiA3D.kp_from_experiment=false;  // Use experimental partitioning constant at Tref?
  BiA3D.Tref=298;         // Temperature of reference (K)
  BiA3D.deltaH=109.0;     // Enthalpy of vaporization (kJ/mol)
  BiA3D.Henry=2.67e8;     // Henry's law constant at Tref (M/atm)
  BiA3D.aq_type="diacid"; // "none","diacid","monoacid" or "aldehyde"
  BiA3D.Kacidity1=3.95e-4;    // First acidity constant
  BiA3D.Kacidity2=7.70e-6;    // Second acidity constant
  BiA3D.hydrophilic=true;   // Does the species condense on the aqueous phase?
  BiA3D.hydrophobic=false;  // Does the species condense on the organic phase?
  BiA3D.nonvolatile=true; // Is the compound nonvolatile?
  BiA3D.is_organic=true;  // Is the compound organic?
  BiA3D.compute_gamma_org=true;  // Compute the activity coefficients of the organic phase for this compound?
  BiA3D.compute_gamma_aq=true;  // Compute the activity coefficients of the aqueous phase for this compound?
  BiA3D.Koligo_org=0.0;         //oligomeriation constant in the organic phase
  BiA3D.rho=1300.0;
  BiA3D.is_monomer=false;
  BiA3D.rion=false;
  //  BiA3D.accomodation_coefficient=0.5;
  // BiA3D.KDiffusion_air=1.0e-5;
  BiA3D.viscosity=1.68e12;
  BiA3D.is_solid=false;
  BiA3D.is_generic=false;

  //Group: if no functionnal group in the species use the default species
  //for the computation of activity coefficients
  //
  double group_tmp_bia3d [] = {2.0,1.0,1.0,1.0, // group C
                               0.0,0.0,0.0,0.0, //group C[OH]
                               0.0,0.0,0.0,0.0, //group Calcohol
                               0.0,0.0,0.0,0.0, //group Calcohol-tail
                               0.0,0.0,0.0,0.0,0.0, //group C=C
                               0.0,0.0, //group aromatic carbon (AC)
                               0.0,0.0,0.0, // group //AC-C
                               0.0,  //group OH
                               0.0, //group H2O
                               0.0, //group ACOH
                               0.0,0.0, //group ketone
                               0.0,   //group aldehyde  
                               0.0,0.0, //group ester
                               0.0,0.0,0.0, //group ether 
                               3.0,  //group acid
                               0.0,   //group ACNO2
                               0.0,0.0,0.0, //group NO3
                               0.0,0.0,0.0, //group CO-OH
			       0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0, //group CO-OC
			       0.0,  //group PAN
			       0.0,  //group CO-OOH
			       0.0,  //group O=COC=O
			       0.0,0.0,0.0}; //group CHxNO2

  size = sizeof(group_tmp_bia3d)/sizeof(double);
  assert(size == 60);
	
  for(int i = 0; i < size; ++i)
    BiA3D.groups[i] = group_tmp_bia3d[i];

  // Search the species name in the aerosol species list 
  // and add the species if its name matches with
  // the given list.
  add_species_ssh(surrogate, BiA3D, species_list_aer, molecular_weight_aer,
		  accomodation_coefficient,diffusion_coef,species_part,nlayer,i_hydrophilic,
		  N_inert, N_inorganic);

  species ACIDMAL;
  ACIDMAL.name="ACIDMAL";
  ACIDMAL.is_inorganic_precursor=false;
  ACIDMAL.Psat_ref=4.59e-8; // Saturation vapor pressure at Tref (torr)
  ACIDMAL.Tref=298;
  ACIDMAL.deltaH= 50.0;     // Enthalpy of vaporization (kJ/mol)
  ACIDMAL.Henry=0.0;     // Henry's law constant at Tref (M/atm)
  ACIDMAL.aq_type="aldehyde"; // "none","diacid","monoacid" or "aldehyde"
  ACIDMAL.hydrophilic=true; // Does the species condense on the aqueous phase?
  ACIDMAL.hydrophobic=false;  // Does the species condense on the organic phase?
  ACIDMAL.nonvolatile=false;  // Is the compound nonvolatile?
  ACIDMAL.kp_from_experiment= false;  // Use experimental partitioning constant at Tref?
  ACIDMAL.kp_experiment= 2.56;       // Value of the experimental partitioning constant at Tref?
  ACIDMAL.is_organic=true;  // Is the compound organic?
  ACIDMAL.compute_gamma_org=true;  // Compute the activity coefficients of the organic phase for this compound?
  ACIDMAL.compute_gamma_aq=true;  // Compute the activity coefficients of the aqueous phase for this compound
  ACIDMAL.Koligo_org=0.0;
  //Parameters for the oligomerization of aldehyde in the aqueous phase as BiA0D:  
  if(with_oligomerization) ACIDMAL.Koligo_aq = 0.1;
  else ACIDMAL.Koligo_aq = 0.;
  ACIDMAL.pHref = 6.0;
  ACIDMAL.beta = 1.91;
  ACIDMAL.rho=1300.0;
  ACIDMAL.is_monomer=false;
  ACIDMAL.rion=false;
  //ACIDMAL.KDiffusion_air=1.0e-5;
  ACIDMAL.viscosity=1.68e12;
  ACIDMAL.is_solid = false;
  ACIDMAL.is_generic=false;

  //Group: if no functionnal group in the species use the default species
  //for the computation of activity coefficients 
    
  double group_tmp_acidmal [] = {0.0,0.0,0.0,0.0, // group C
                                 0.0,0.0,0.0,0.0, //group C[OH]
                                 0.0,0.0,0.0,0.0, //group Calcohol
                                 0.0,0.0,0.0,0.0, //group Calcohol-tail  
                                 0.0,1.0,0.0,0.0,1.0, //group C=C
                                 0.0,0.0, //group aromatic carbon (AC)
                                 0.0,0.0,0.0, // group //AC-C
                                 2.0,  //group OH
                                 0.0, //group H2O
                                 0.0, //group ACOH
                                 0.0,0.0, //group ketone
                                 1.0,   //group aldehyde  
                                 0.0,0.0, //group ester
                                 0.0,0.0,0.0, //group ether  
                                 1.0,  //group acid
                                 0.0,   //group ACNO2
                                 0.0,0.0,0.0, //group NO3
                                 0.0,0.0,0.0, //group CO-OH
				 0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0, //group CO-OC
				 0.0,  //group PAN
				 0.0,  //group CO-OOH
			       0.0,  //group O=COC=O
			       0.0,0.0,0.0}; //group CHxNO2

  size = sizeof(group_tmp_acidmal)/sizeof(double);
  assert(size == 60);
  for(int i = 0; i < size; ++i)
    ACIDMAL.groups[i] = group_tmp_acidmal[i];

  // Search the species name in the aerosol species list 
  // and add the species if its name matches with
  // the given list.
  add_species_ssh(surrogate, ACIDMAL, species_list_aer, molecular_weight_aer,
		  accomodation_coefficient,diffusion_coef,species_part,nlayer,i_hydrophilic,
		  N_inert, N_inorganic);

  species DHMB;
  DHMB.name="DHMB";
  DHMB.is_inorganic_precursor=false;
  DHMB.Psat_ref=3.52e-6; // Saturation vapor pressure at Tref (torr)
  DHMB.Tref=298;
  DHMB.MM=154;           // Molar mass (g/mol)
  DHMB.deltaH= 50.0;     // Enthalpy of vaporization (kJ/mol)
  DHMB.Henry=0.0;     // Henry's law constant at Tref (M/atm)
  DHMB.aq_type="none"; // "none","diacid","monoacid" or "aldehyde"
  DHMB.hydrophilic=true; // Does the species condense on the aqueous phase?
  DHMB.hydrophobic=false;  // Does the species condense on the organic phase?
  DHMB.nonvolatile=false;  // Is the compound nonvolatile?
  DHMB.kp_from_experiment=false;  // Use experimental partitioning constant at Tref?
  DHMB.kp_experiment= 0.034;      // Value of the experimental partitioning constant at Tref?
  DHMB.is_organic=true;  // Is the compound organic?
  DHMB.compute_gamma_org=true;  // Compute the activity coefficients of the organic phase for this compound?
  DHMB.compute_gamma_aq=true;  // Compute the activity coefficients of the aqueous phase for this compound
  DHMB.Koligo_org=0.0;
  DHMB.rho=1300.0;
  DHMB.is_monomer=false;
  DHMB.rion=false;
  //DHMB.KDiffusion_air=1.0e-5;
  DHMB.viscosity=1.68e12;  
  DHMB.is_solid=false;
  DHMB.is_generic=false;
  
  //Group: if no functionnal group in the species use the default species
  //for the computation of activity coefficients 
  
  double group_tmp_dhmb [] = {1.0,0.0,0.0,0.0, // group C
                              0.0,0.0,0.0,0.0, //group C[OH]
                              0.0,0.0,0.0,0.0, //group Calcohol
                              0.0,0.0,0.0,0.0, //group Calcohol-tail  
                              0.0,0.0,0.0,0.0,1.0, //group C=C
                              0.0,0.0, //group aromatic carbon (AC)
                              0.0,0.0,0.0, // group //AC-C
                              2.0,  //group OH
                              0.0, //group H2O
                              0.0, //group ACOH
                              0.0,2.0, //group ketone
                              0.0,   //group aldehyde  
                              0.0,0.0, //group ester
                              0.0,0.0,0.0, //group ether  
                              0.0,  //group acid
                              0.0,   //group ACNO2
                              0.0,0.0,0.0, //group NO3
                              0.0,0.0,0.0, //group CO-OH
			      0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0, //group CO-OC
			      0.0,  //group PAN
			      0.0,  //group CO-OOH
			       0.0,  //group O=COC=O
			       0.0,0.0,0.0}; //group CHxNO2

  size = sizeof(group_tmp_dhmb)/sizeof(double);
  assert(size == 60);
  for(int i = 0; i < size; ++i)
    DHMB.groups[i] = group_tmp_dhmb[i];

  // Search the species name in the aerosol species list 
  // and add the species if its name matches with
  // the given list.
  add_species_ssh(surrogate, DHMB, species_list_aer, molecular_weight_aer,
		  accomodation_coefficient,diffusion_coef,species_part,nlayer,i_hydrophilic,
		  N_inert, N_inorganic);

  species PAHlN;
  PAHlN.name="PAHlN";
  PAHlN.is_inorganic_precursor=false;
  PAHlN.Psat_ref=1e-12; // Saturation vapor pressure at Tref (torr)
  PAHlN.Tref=298;
  PAHlN.MM=198;           // Molar mass (g/mol)
  PAHlN.deltaH= 50.0;     // Enthalpy of vaporization (kJ/mol)
  PAHlN.Henry=0.0;     // Henry's law constant at Tref (M/atm)
  PAHlN.aq_type="diacid"; // "none","diacid","monoacid" or "aldehyde"
  PAHlN.Kacidity1=3.95e-4;    // First acidity constant as BiA3D
  PAHlN.Kacidity2=7.70e-6;    // Second acidity constant
  PAHlN.hydrophilic=true; // Does the species condense on the aqueous phase?
  PAHlN.hydrophobic=false;  // Does the species condense on the organic phase?
  PAHlN.nonvolatile=false;  // Is the compound nonvolatile?
  PAHlN.kp_from_experiment=false;  // Use experimental partitioning constant at Tref?
  PAHlN.kp_experiment=93817.62 ;       // Value of the experimental partitioning constant at Tref?
  PAHlN.is_organic=true;  // Is the compound organic?
  PAHlN.compute_gamma_org=true;  // Compute the activity coefficients of the organic phase for this compound?
  PAHlN.compute_gamma_aq=true;  // Compute the activity coefficients of the aqueous phase for this compound
  PAHlN.Koligo_org=0.0;
  PAHlN.rho=1300.0;
  PAHlN.is_monomer=false;
  PAHlN.rion=false;
  //PAHlN.KDiffusion_air=1.0e-5;
  PAHlN.viscosity=1.68e12;  
  PAHlN.is_solid=false;
  PAHlN.is_generic=false;
  
  //Group: if no functionnal group in the species use the default species
  //for the computation of activity coefficients 
  
  double group_tmp_pahln [] = {0.0,0.0,0.0,0.0, // group C
                               0.0,0.0,0.0,0.0, //group C[OH]
                               0.0,0.0,0.0,0.0, //group Calcohol
                               0.0,0.0,0.0,0.0, //group Calcohol-tail  
                               0.0,0.0,0.0,0.0,0.0, //group C=C
                               2.0,2.0, //group aromatic carbon (AC)
                               0.0,0.0,0.0, // group //AC-C
                               0.0,  //group OH
                               0.0, //group H2O
                               2.0, //group ACOH
                               0.0,0.0, //group ketone
                               0.0,   //group aldehyde  
                               0.0,0.0, //group ester
                               0.0,0.0,0.0, //group ether  
                               2.0,  //group acid
                               0.0,   //group ACNO2
                               0.0,0.0,0.0, //group NO3
                               0.0,0.0,0.0, //group CO-OH
			       0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0, //group CO-OC
			       0.0,  //group PAN
			       0.0,  //group CO-OOH
			       0.0,  //group O=COC=O
			       0.0,0.0,0.0}; //group CHxNO2

  size = sizeof(group_tmp_pahln)/sizeof(double);
  assert(size == 60);
  for(int i = 0; i < size; ++i)
    PAHlN.groups[i] = group_tmp_pahln[i];

  // Search the species name in the aerosol species list 
  // and add the species if its name matches with
  // the given list.
  add_species_ssh(surrogate, PAHlN, species_list_aer, molecular_weight_aer,
		  accomodation_coefficient,diffusion_coef,species_part,nlayer,i_hydrophilic,
		  N_inert, N_inorganic);
  
  species PAHhN;
  PAHhN.name="PAHhN";
  PAHhN.is_inorganic_precursor=false;
  PAHhN.Psat_ref=1e-6; // Saturation vapor pressure at Tref (torr)
  PAHhN.Tref=298;
  PAHhN.MM=166;           // Molar mass (g/mol)
  PAHhN.deltaH= 50.0;     // Enthalpy of vaporization (kJ/mol)
  PAHhN.Henry=0.0;     // Henry's law constant at Tref (M/atm)
  PAHhN.aq_type="diacid"; // "none","diacid","monoacid" or "aldehyde"
  PAHhN.Kacidity1=3.95e-4;    // First acidity constant as BiA3D
  PAHhN.Kacidity2=7.70e-6;    // Second acidity constant
  PAHhN.hydrophilic=true; // Does the species condense on the aqueous phase?
  PAHhN.hydrophobic=false;  // Does the species condense on the organic phase?
  PAHhN.nonvolatile=false;  // Is the compound nonvolatile?
  PAHhN.kp_from_experiment=false;  // Use experimental partitioning constant at Tref?
  PAHhN.kp_experiment=0.1119 ;       // Value of the experimental partitioning constant at Tref?
  PAHhN.is_organic=true;  // Is the compound organic?
  PAHhN.compute_gamma_org=true;  // Compute the activity coefficients of the organic phase for this compound?
  PAHhN.compute_gamma_aq=true;  // Compute the activity coefficients of the aqueous phase for this compound
  PAHhN.Koligo_org=0.0;
  PAHhN.rho=1300.0;
  PAHhN.is_monomer=false;
  PAHhN.rion=false;
  //PAHhN.KDiffusion_air=1.0e-5;
  PAHhN.viscosity=1.68e12;  
  PAHhN.is_solid=false;
  PAHhN.is_generic=false;
   
  //Group: if no functionnal group in the species use the default species
  //for the computation of activity coefficients 
  
  double group_tmp_pahhn [] = {0.0,0.0,0.0,0.0, // group C
                               0.0,0.0,0.0,0.0, //group C[OH]
                               0.0,0.0,0.0,0.0, //group Calcohol
                               0.0,0.0,0.0,0.0, //group Calcohol-tail  
                               0.0,0.0,0.0,0.0,0.0, //group C=C
                               4.0,2.0, //group aromatic carbon (AC)
                               0.0,0.0,0.0, // group //AC-C
                               0.0,  //group OH
                               0.0, //group H2O
                               0.0, //group ACOH
                               0.0,0.0, //group ketone
                               0.0,   //group aldehyde  
                               0.0,0.0, //group ester
                               0.0,0.0,0.0, //group ether  
                               2.0,  //group acid
                               0.0,   //group ACNO2
                               0.0,0.0,0.0, //group NO3
                               0.0,0.0,0.0, //group CO-OH
			       0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0, //group CO-OC
			       0.0,  //group PAN
			       0.0,  //group CO-OOH
			       0.0,  //group O=COC=O
			       0.0,0.0,0.0}; //group CHxNO2

  size = sizeof(group_tmp_pahhn)/sizeof(double);
  assert(size == 60);
  for(int i = 0; i < size; ++i)
    PAHhN.groups[i] = group_tmp_pahhn[i];

  // Search the species name in the aerosol species list 
  // and add the species if its name matches with
  // the given list.
  add_species_ssh(surrogate, PAHhN, species_list_aer, molecular_weight_aer,
		  accomodation_coefficient,diffusion_coef,species_part,nlayer,i_hydrophilic,
		  N_inert, N_inorganic);

  species PSYR;
  PSYR.name="PSYR";
  PSYR.is_inorganic_precursor=false;
  PSYR.Psat_ref=7.53e-6; // Saturation vapor pressure at Tref (torr)
  PSYR.Tref=298;
  PSYR.MM=186;           // Molar mass (g/mol)
  PSYR.deltaH= 50.0;     // Enthalpy of vaporization (kJ/mol)
  PSYR.Henry=0.0;     // Henry's law constant at Tref (M/atm)
  PSYR.aq_type="none"; // "none","diacid","monoacid" or "aldehyde"
  PSYR.hydrophilic=true; // Does the species condense on the aqueous phase?
  PSYR.hydrophobic=false;  // Does the species condense on the organic phase?
  PSYR.nonvolatile=false;  // Is the compound nonvolatile?
  PSYR.kp_from_experiment=false;  // Use experimental partitioning constant at Tref?
  PSYR.kp_experiment=1.294e-2 ;       // Value of the experimental partitioning constant at Tref?
  PSYR.is_organic=true;  // Is the compound organic?
  PSYR.compute_gamma_org=true;  // Compute the activity coefficients of the organic phase for this compound?
  PSYR.compute_gamma_aq=true;  // Compute the activity coefficients of the aqueous phase for this compound
  PSYR.Koligo_org=0.0;
  PSYR.rho=1300.0;
  PSYR.is_monomer=false;
  PSYR.rion=false;
  //PSYR.KDiffusion_air=1.0e-5;
  PSYR.viscosity=1.68e12; 
  PSYR.is_solid=false;
  PSYR.is_generic=false;

  //Group: if no functionnal group in the species use the default species
  //for the computation of activity coefficients 
  
  double group_tmp_psyr [] = {0.0,0.0,0.0,0.0, // group C
                              0.0,0.0,0.0,0.0, //group C[OH]
                              0.0,0.0,0.0,0.0, //group Calcohol
                              0.0,0.0,0.0,0.0, //group Calcohol-tail  
                              0.0,0.0,0.0,1.0,1.0, //group C=C
                              0.0,0.0, //group aromatic carbon (AC)
                              0.0,0.0,0.0, // group //AC-C
                              1.0,  //group OH
                              0.0, //group H2O
                              0.0, //group ACOH
                              0.0,0.0, //group ketone
                              2.0,   //group aldehyde  
                              0.0,0.0, //group ester
                              2.0,0.0,0.0, //group ether  
                              0.0,  //group acid
                              0.0,   //group ACNO2
                              0.0,0.0,0.0, //group NO3
                              0.0,0.0,0.0, //group CO-OH
			      0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0, //group CO-OC
			      0.0,  //group PAN
			      0.0,  //group CO-OOH
			       0.0,  //group O=COC=O
			       0.0,0.0,0.0}; //group CHxNO2

  size = sizeof(group_tmp_psyr)/sizeof(double);
  assert(size == 60);
  for(int i = 0; i < size; ++i)
    PSYR.groups[i] = group_tmp_psyr[i];

  // Search the species name in the aerosol species list 
  // and add the species if its name matches with
  // the given list.
  add_species_ssh(surrogate, PSYR, species_list_aer, molecular_weight_aer,
		  accomodation_coefficient,diffusion_coef,species_part,nlayer,i_hydrophilic,
		  N_inert, N_inorganic);

  species GHDPerox;
  GHDPerox.name="GHDPerox";
  GHDPerox.is_inorganic_precursor=false;
  GHDPerox.Psat_ref=5.41e-7; // Saturation vapor pressure at Tref (torr)
  GHDPerox.Tref=298;
  GHDPerox.MM=174;           // Molar mass (g/mol)
  GHDPerox.deltaH=50.0;     // Enthalpy of vaporization (kJ/mol)
  GHDPerox.Henry=0.0;     // Henry's law constant at Tref (M/atm)
  GHDPerox.aq_type="none"; // "none","diacid","monoacid" or "aldehyde"
  GHDPerox.hydrophilic=true; // Does the species condense on the aqueous phase?
  GHDPerox.hydrophobic=false;  // Does the species condense on the organic phase?
  GHDPerox.nonvolatile=false;  // Is the compound nonvolatile?
  GHDPerox.kp_from_experiment=false;  // Use experimental partitioning constant at Tref?
  GHDPerox.kp_experiment=0.1972 ;       // Value of the experimental partitioning constant at Tref?
  GHDPerox.is_organic=true;  // Is the compound organic?
  GHDPerox.compute_gamma_org=true;  // Compute the activity coefficients of the organic phase for this compound?
  GHDPerox.compute_gamma_aq=true;  // Compute the activity coefficients of the aqueous phase for this compound
  GHDPerox.Koligo_org=0.0;
  GHDPerox.rho=1300.0;
  GHDPerox.is_monomer=false;
  GHDPerox.rion=false;
  //GHDPerox.KDiffusion_air=1.0e-5;
  GHDPerox.viscosity=1.68e12; 
  GHDPerox.is_solid=false;
  GHDPerox.is_generic=false;
  
  //Group: if no functionnal group in the species use the default species
  //for the computation of activity coefficients 
  
  double group_tmp_ghdperox [] = {0.0,0.0,0.0,1.0, // group C
                                  0.0,0.0,1.0,0.0, //group C[OH]
                                  0.0,0.0,0.0,0.0, //group Calcohol
                                  0.0,0.0,0.0,0.0, //group Calcohol-tail  
                                  0.0,1.0,0.0,0.0,1.0, //group C=C
                                  0.0,0.0, //group aromatic carbon (AC)
                                  0.0,0.0,0.0, // group //AC-C
                                  2.0,  //group OH
                                  0.0, //group H2O
                                  0.0, //group ACOH
                                  0.0,0.0, //group ketone
                                  0.0,   //group aldehyde  
                                  0.0,0.0, //group ester
                                  1.0,0.0,0.0, //group ether  
                                  0.0,  //group acid
                                  0.0,   //group ACNO2
                                  0.0,0.0,0.0, //group NO3
                                  0.0,1.0,0.0, //group CO-OH
				  0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0, //group CO-OC
				  0.0,  //group PAN
				  0.0,  //group CO-OOH
			          0.0,  //group O=COC=O
			          0.0,0.0,0.0}; //group CHxNO2

  size = sizeof(group_tmp_ghdperox)/sizeof(double);
  assert(size == 60);
  for(int i = 0; i < size; ++i)
    GHDPerox.groups[i] = group_tmp_ghdperox[i];

  // Search the species name in the aerosol species list 
  // and add the species if its name matches with
  // the given list.
  add_species_ssh(surrogate, GHDPerox, species_list_aer, molecular_weight_aer,
		  accomodation_coefficient,diffusion_coef,species_part,nlayer,i_hydrophilic,
		  N_inert, N_inorganic);

  add_generic_species_ssh(config, surrogate, species_list_aer, molecular_weight_aer, accomodation_coefficient,
			  aerosol_type, species_smiles, saturation_vapor_pressure, enthalpy_vaporization, 
			  species_part,nlayer,i_hydrophilic,
		          N_inert, N_inorganic);

  
  species H2O;
  H2O.name="H2O";
  H2O.is_inorganic_precursor=false;
  H2O.is_organic=false;  // Is the compound organic?
  H2O.hydrophilic=true; // Does the species condense on the aqueous phase?
  H2O.hydrophobic=true;  // Does the species condense on the organic phase?
  H2O.nonvolatile=false;
  H2O.kp_from_experiment=false; 
  H2O.compute_gamma_org=true;  // Compute the activity coefficients of the organic phase for this compound?
  H2O.compute_gamma_aq=true;  // Compute the activity coefficients of the aqueous phase for this compound
  H2O.Koligo_org=0.0;
  H2O.rho=1000.0;
  //H2O.KDiffusion_air=1.0e-5;
  //  H2O.accomodation_coefficient=alpha;
  H2O.viscosity=1.0;  
  H2O.is_solid=false;
  H2O.is_monomer=false;
  H2O.rion=false;

  //Group: if no functionnal group in the species use the default species
  //for the computation of activity coefficients 
  double group_tmp_h2o [] = {0.0,0.0,0.0,0.0, // group C
			     0.0,0.0,0.0,0.0, //group C[OH]
			     0.0,0.0,0.0,0.0, //group Calcohol
			     0.0,0.0,0.0,0.0, //group Calcohol-tail  
			     0.0,0.0,0.0,0.0,0.0, //group C=C
			     0.0,0.0, //group aromatic carbon (AC)
			     0.0,0.0,0.0, // group //AC-C
			     0.0,  //group OH
			     1.0, //group H2O
			     0.0, //group ACOH
			     0.0,0.0, //group ketone
			     0.0,   //group aldehyde  
			     0.0,0.0, //group ester
			     0.0,0.0,0.0, //group ether 
			     0.0,  //group acid
			     0.0,   //group ACNO2
			     0.0,0.0,0.0, //group NO3
			     0.0,0.0,0.0, //group CO-OH
			     0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0, //group CO-OC
			     0.0,  //group PAN
			     0.0,  //group CO-OOH
			       0.0,  //group O=COC=O
			       0.0,0.0,0.0}; //group CHxNO2

  size = sizeof(group_tmp_h2o)/sizeof(double);
  assert(size == 60);
  
  for(int i = 0; i < size; ++i)
    H2O.groups[i] = group_tmp_h2o[i];

  // Search the species name in the aerosol species list 
  // and add the species if its name matches with
  // the given list.
  add_species_ssh(surrogate, H2O, species_list_aer, molecular_weight_aer,
		  accomodation_coefficient,diffusion_coef,species_part,nlayer,i_hydrophilic,
		  N_inert, N_inorganic);
  
  species SO4;
  SO4.name="SO4";
  SO4.is_inorganic_precursor=false;
  SO4.MM=96.0;           // Molar mass (g/mol)
  SO4.is_organic=false;  // Is the compound organic?
  SO4.hydrophilic=true; // Does the species condense on the aqueous phase?
  SO4.hydrophobic=false;  // Does the species condense on the organic phase?
  SO4.nonvolatile=false;
  SO4.kp_from_experiment=false;
  SO4.compute_gamma_org=false;  // Compute the activity coefficients of the organic phase for this compound?
  SO4.compute_gamma_aq=false;  // Compute the activity coefficients of the aqueous phase for this compound
  SO4.charge=-2;
  SO4.index_ion_aiomfac=11;
  SO4.rho=1300.0; //1400.0;    
  SO4.KDiffusion_air=1.0e-5;
  SO4.accomodation_coefficient=0.5;
  SO4.viscosity=1.0;
  SO4.soap_ind = -1;
  SO4.soap_ind_aero = -1;
  SO4.is_solid=false;
  SO4.is_monomer=false;
  SO4.rion=false;

  for (int i = 0; i < nsp; ++i)
    if (species_list_aer[i].substr(1,-1) == "SO4")
      {
        SO4.soap_ind = i;
        SO4.soap_ind_aero = i;
      }
  surrogate.push_back(SO4);


  species HSO4;
  HSO4.name="HSO4";
  HSO4.is_inorganic_precursor=false;
  HSO4.MM=97.0;           // Molar mass (g/mol)
  HSO4.is_organic=false;  // Is the compound organic?
  HSO4.hydrophilic=true; // Does the species condense on the aqueous phase?
  HSO4.hydrophobic=false;  // Does the species condense on the organic phase?
  HSO4.nonvolatile=false;
  HSO4.compute_gamma_org=false;  // Compute the activity coefficients of the organic phase for this compound?
  HSO4.compute_gamma_aq=false;  // Compute the activity coefficients of the aqueous phase for this compound
  HSO4.charge=-1;
  HSO4.index_ion_aiomfac=10;
  HSO4.rho=1300.0; //1400.0;
  HSO4.KDiffusion_air=1.0e-5;
  HSO4.accomodation_coefficient=0.5;
  HSO4.viscosity=1.0;
  HSO4.soap_ind = -1;
  HSO4.soap_ind_aero = -1;
  HSO4.is_solid=false;
  HSO4.is_monomer=false;
  HSO4.rion=false;

  for (int i = 0; i < nsp; ++i)
    if (species_list_aer[i].substr(1,-1) == "SO4")
      {
        HSO4.soap_ind = i;
        HSO4.soap_ind_aero = i;
      }
  surrogate.push_back(HSO4);

  species NO3;
  NO3.name="NO3";
  NO3.is_inorganic_precursor=false;
  NO3.MM=62.0;           // Molar mass (g/mol)
  NO3.is_organic=false;  // Is the compound organic?
  NO3.hydrophilic=true; // Does the species condense on the aqueous phase?
  NO3.hydrophobic=false;  // Does the species condense on the organic phase?
  NO3.nonvolatile=false;
  NO3.kp_from_experiment=false;
  NO3.compute_gamma_org=false;  // Compute the activity coefficients of the organic phase for this compound?
  NO3.compute_gamma_aq=false;  // Compute the activity coefficients of the aqueous phase for this compound
  NO3.charge=-1;
  NO3.index_ion_aiomfac=9;
  NO3.rho=1300.0; //1400.0;
  NO3.KDiffusion_air=1.0e-5;
  NO3.accomodation_coefficient=0.5;
  NO3.viscosity=1.0;
  NO3.soap_ind = -1;
  NO3.soap_ind_aero = -1;
  NO3.is_solid=false;
  NO3.is_monomer=false;
  NO3.rion=false;

  for (int i = 0; i < nsp; ++i)
    if (species_list_aer[i].substr(1,-1) == "NO3")
      {
        NO3.soap_ind = i;
        NO3.soap_ind_aero = i;
      }
  surrogate.push_back(NO3);

  species NH4;
  NH4.name="NH4";
  NH4.is_inorganic_precursor=false;
  NH4.MM=18.0;           // Molar mass (g/mol)
  NH4.is_organic=false;  // Is the compound organic?
  NH4.hydrophilic=true; // Does the species condense on the aqueous phase?
  NH4.hydrophobic=false;  // Does the species condense on the organic phase?
  NH4.nonvolatile=false;
  NH4.kp_from_experiment=false;
  NH4.compute_gamma_org=false;  // Compute the activity coefficients of the organic phase for this compound?
  NH4.compute_gamma_aq=false;  // Compute the activity coefficients of the aqueous phase for this compound
  NH4.charge=1.0;
  NH4.index_ion_aiomfac=4;
  NH4.rho=1300.0; //1400.0;
  NH4.KDiffusion_air=1.0e-5;
  NH4.accomodation_coefficient=0.5;
  NH4.viscosity=1.0;
  NH4.soap_ind = -1;
  NH4.soap_ind_aero = -1;
  NH4.is_solid=false;
  NH4.is_monomer=false;
  NH4.rion=false;

  for (int i = 0; i < nsp; ++i)
    if (species_list_aer[i].substr(1,-1) == "NH4")
      {
        NH4.soap_ind = i;
        NH4.soap_ind_aero = i;
      }
  surrogate.push_back(NH4);

  species H;
  H.name="H";
  H.is_inorganic_precursor=false;
  H.MM=1.0;           // Molar mass (g/mol)
  H.is_organic=false;  // Is the compound organic?
  H.hydrophilic=true; // Does the species condense on the aqueous phase?
  H.hydrophobic=false;  // Does the species condense on the organic phase?
  H.nonvolatile=false;
  H.kp_from_experiment=false;
  H.compute_gamma_org=false;  // Compute the activity coefficients of the organic phase for this compound?
  H.compute_gamma_aq=false;  // Compute the activity coefficients of the aqueous phase for this compound
  H.charge=1.0;
  H.rho=1300.0; //1400.0;
  H.index_ion_aiomfac=0;
  H.KDiffusion_air=1.0e-5;
  H.accomodation_coefficient=0.5;
  H.viscosity=1.0;
  H.soap_ind = -1;
  H.soap_ind_aero = -1;
  H.is_solid=false;
  H.is_monomer=false;
  H.rion=false;
  surrogate.push_back(H);

  species Na;
  Na.name="Na";
  Na.is_inorganic_precursor=false;
  Na.is_organic=false;  // Is the compound organic?
  Na.hydrophilic=true; // Does the species condense on the aqueous phase?
  Na.hydrophobic=false;  // Does the species condense on the organic phase?
  Na.nonvolatile=false;
  Na.kp_from_experiment=false;
  Na.compute_gamma_org=false;  // Compute the activity coefficients of the organic phase for this compound?
  Na.compute_gamma_aq=false;  // Compute the activity coefficients of the aqueous phase for this compound
  Na.charge=1.0;
  Na.index_ion_aiomfac=2;
  Na.rho=1300.0; //1400.0;
  Na.KDiffusion_air=1.0e-5;
  Na.viscosity=1.0;
  Na.is_solid=false;
  Na.is_monomer=false;
  Na.rion=false;

  // Find the number in the aerosol species list
  Na.soap_ind = -1;
  Na.soap_ind_aero = -1;
  for (int i = 0; i < nsp; ++i)
    if (species_list_aer[i].substr(1,-1) == "NA" or species_list_aer[i].substr(1,-1) == "Na")
      {
        Na.soap_ind = i;
        Na.soap_ind_aero = i;
        Na.MM =  molecular_weight_aer[i] / 1.e6;
        Na.accomodation_coefficient = accomodation_coefficient[i];
      }
  surrogate.push_back(Na);

  species Cl;
  Cl.name="Cl";
  Cl.is_inorganic_precursor=false;
  Cl.MM=35.5;           // Molar mass (g/mol)
  Cl.is_organic=false;  // Is the compound organic?
  Cl.hydrophilic=true; // Does the species condense on the aqueous phase?
  Cl.hydrophobic=false;  // Does the species condense on the organic phase?
  Cl.nonvolatile=false;
  Cl.kp_from_experiment=false;
  Cl.compute_gamma_org=false;  // Compute the activity coefficients of the organic phase for this compound?
  Cl.compute_gamma_aq=false;  // Compute the activity coefficients of the aqueous phase for this compound
  Cl.charge=-1.0;
  Cl.index_ion_aiomfac=7;
  Cl.rho=1300.0; // 1400.0;
  Cl.KDiffusion_air=1.0e-5;
  Cl.accomodation_coefficient=0.5;
  Cl.viscosity=1.0;
  Cl.soap_ind = -1;
  Cl.soap_ind_aero = -1;
  Cl.is_solid=false;
  Cl.is_monomer=false;
  Cl.rion=false;

  for (int i = 0; i < nsp; ++i)
    if (species_list_aer[i].substr(1,-1) == "HCL")
      {
        Cl.soap_ind = i;
        Cl.soap_ind_aero = i;
      }
  surrogate.push_back(Cl);

      species H2SO4;
      H2SO4.name="H2SO4";
      H2SO4.nonvolatile=true;
      H2SO4.is_organic=false;
      H2SO4.is_inorganic_precursor=true;
      H2SO4.Henry=1.0e13;     // Enthalpy of vaporization (kJ/mol)
      H2SO4.deltaH=0.0;     // Henry's law constant at Tref (M/atm)
      H2SO4.hydrophilic=true; // Does the species condense on the aqueous phase?
      H2SO4.hydrophobic=false;  // Does the species condense on the organic phase?
      H2SO4.compute_gamma_org=false;  // Compute the activity coefficients of the organic phase for this compound?
      H2SO4.compute_gamma_aq=false;  // Compute the activity coefficients of the aqueous phase for this compound
      //H2SO4.KDiffusion_air=1.07e-5;
      H2SO4.viscosity=1.0;
      H2SO4.is_solid=false;
  
      // Find the number in the aerosol species list
      H2SO4.soap_ind = -1;
      H2SO4.soap_ind_aero = -1;
      for (int i = 0; i < nsp; ++i)
	if (species_list_aer[i].substr(1,-1) == "SO4")
	  {
	    H2SO4.soap_ind = i;
	    H2SO4.soap_ind_aero = i;
	    H2SO4.MM =  molecular_weight_aer[i] / 1.e6;
	    H2SO4.accomodation_coefficient = accomodation_coefficient[i];
	    H2SO4.KDiffusion_air = diffusion_coef[i];
	  }
      surrogate.push_back(H2SO4);

      species NH3;
      NH3.name="NH3";
      //  NH3.MM=17.0;
      NH3.is_organic=false;
      NH3.nonvolatile=false;
      NH3.is_inorganic_precursor=true;
      NH3.Henry=57.64;     // Enthalpy of vaporization (kJ/mol)
      NH3.deltaH=0.0;     // Henry's law constant at Tref (M/atm)
      NH3.hydrophilic=true; // Does the species condense on the aqueous phase?
      NH3.hydrophobic=false;  // Does the species condense on the organic phase?
      NH3.compute_gamma_org=false;  // Compute the activity coefficients of the organic phase for this compound?
      NH3.compute_gamma_aq=false;  // Compute the activity coefficients of the aqueous phase for this compound
      //  NH3.accomodation_coefficient=0.5;
      NH3.KDiffusion_air=2.17e-5;
      NH3.viscosity=1.0;
      NH3.is_solid=false;
  
      // Find the number in the aerosol species list
      NH3.soap_ind = -1;
      NH3.soap_ind_aero = -1;
      for (int i = 0; i < nsp; ++i)
	if (species_list_aer[i].substr(1,-1) == "NH4")
	  {
	    NH3.soap_ind = i;
	    NH3.soap_ind_aero = i;
	    NH3.MM =  molecular_weight_aer[i] / 1.e6;
	    NH3.accomodation_coefficient = accomodation_coefficient[i];
	    NH3.KDiffusion_air = diffusion_coef[i];
	  }
      surrogate.push_back(NH3);

      species HNO3;
      HNO3.name="HNO3";
      //  HNO3.MM=63.0;
      HNO3.nonvolatile=false;
      HNO3.is_organic=false;
      HNO3.is_inorganic_precursor=true;
      HNO3.Henry=2.1e5;     // Enthalpy of vaporization (kJ/mol)
      HNO3.deltaH=0.0;     // Henry's law constant at Tref (M/atm)
      HNO3.hydrophilic=true; // Does the species condense on the aqueous phase?
      HNO3.hydrophobic=false;  // Does the species condense on the organic phase?
      HNO3.compute_gamma_org=false;  // Compute the activity coefficients of the organic phase for this compound?
      HNO3.compute_gamma_aq=false;  // Compute the activity coefficients of the aqueous phase for this compound
      //  HNO3.accomodation_coefficient=0.5;
      HNO3.KDiffusion_air=1.47e-5;
      HNO3.viscosity=1.0;
      HNO3.is_solid=false;
  
      // Find the number in the aerosol species list
      HNO3.soap_ind = -1;
      HNO3.soap_ind_aero = -1;
      for (int i = 0; i < nsp; ++i)
	if (species_list_aer[i].substr(1,-1) == "NO3")
	  {
	    HNO3.soap_ind = i;
	    HNO3.soap_ind_aero = i;
	    HNO3.MM =  molecular_weight_aer[i] / 1.e6;
	    HNO3.accomodation_coefficient = accomodation_coefficient[i];
	    HNO3.KDiffusion_air = diffusion_coef[i];
	  }
      surrogate.push_back(HNO3);

      species HCl;
      HCl.name="HCl";
      //  HCl.MM=36.5;
      HCl.nonvolatile=false;
      HCl.is_organic=false;
      HCl.is_inorganic_precursor=true;
      HCl.Henry=2.5e3;     // Enthalpy of vaporization (kJ/mol)
      HCl.deltaH=0.0;     // Henry's law constant at Tref (M/atm)
      HCl.hydrophilic=true; // Does the species condense on the aqueous phase?
      HCl.hydrophobic=false;  // Does the species condense on the organic phase?
      HCl.compute_gamma_org=false;  // Compute the activity coefficients of the organic phase for this compound?
      HCl.compute_gamma_aq=false;  // Compute the activity coefficients of the aqueous phase for this compound
      //  HCl.accomodation_coefficient=0.5;
      HCl.KDiffusion_air=1.72e-5;
      HCl.viscosity=1.0;
      HCl.is_solid=false;
  
      // Find the number in the aerosol species list
      HCl.soap_ind = -1;
      HCl.soap_ind_aero = -1;
      for (int i = 0; i < nsp; ++i)
	if (species_list_aer[i].substr(1,-1) == "HCL")
	  {
	    HCl.soap_ind = i;
	    HCl.soap_ind_aero = i;
	    HCl.MM =  molecular_weight_aer[i] / 1.e6;
	    HCl.accomodation_coefficient = accomodation_coefficient[i];
	    HCl.KDiffusion_air = diffusion_coef[i];
	  }
  
      surrogate.push_back(HCl);
    
}

#endif
