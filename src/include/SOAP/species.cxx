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
		      double saturation_vapor_pressure[],
		      double enthalpy_vaporization[],
		      double henry[],
		      double t_ref[], double mass_density[],
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
	if (current_species.is_organic or current_species.is_inorganic_precursor) 
	  current_species.KDiffusion_air = diffusion_coef[i];
	else
	  current_species.KDiffusion_air = 1.e-5;
	current_species.rho = mass_density[i]*1e9; 
	if (saturation_vapor_pressure[i]>0.) 
	  current_species.Psat_ref = saturation_vapor_pressure[i];

	if (henry[i]>0.) 
	  current_species.Henry = henry[i];

	if (t_ref[i]>0.) 
	  current_species.Tref = t_ref[i];
	else
	  current_species.Tref = 298.;

	if (enthalpy_vaporization[i]>0.)	 
	  current_species.deltaH = enthalpy_vaporization[i];	  
	
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
			     double diffusion_coef[],
			     int aerosol_type[],
			     vector<string> species_smiles,
			     double saturation_vapor_pressure[],
			     double enthalpy_vaporization[],
			     double henry[],
			     double t_ref[], double mass_density[],
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
    if (aerosol_type[i]==4)
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
	    if (species_smiles[i]=="-")
	      for(int igr = 0; igr < 60; ++igr)
		X.groups[igr] = 0.;		
	    
	    X.is_inorganic_precursor=false;
	    X.nonvolatile=false; // Is the compound nonvolatile?
	    if (saturation_vapor_pressure[i]>0.)
	      X.Psat_ref=saturation_vapor_pressure[i]; // Saturation vapor pressure at Tref (torr)
	    else
	      {
		X.Psat_ref=1e-15;
		X.nonvolatile=true;
	      }
	    X.kp_from_experiment=false;  // Use experimental partitioning constant at Tref?
	    if (t_ref[i]>0.)
	      X.Tref=t_ref[i];         // Temperature of reference (K)
	    else
	      X.Tref=298.;
	    
	    X.deltaH=enthalpy_vaporization[i];     // Enthalpy of vaporization (kJ/mol)
	    X.Henry=henry[i];     // Henry's law constant at Tref (M/atm)
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
              	    
	    X.is_organic=true;  // Is the compound organic?
	    X.compute_gamma_org=true;  // Compute the activity coefficients of the organic phase for this compound?
	    X.compute_gamma_aq=true;  // Compute the activity coefficients of the aqueous phase for this compound?
	    X.Koligo_org=0.0;         //oligomeriation constant in the organic phase
	    X.rho=mass_density[i]*1e9;
	    X.is_monomer=false;
	    X.rion=false;
	    X.KDiffusion_air=diffusion_coef[i];
	    // BiA2D.accomodation_coefficient=alpha;
	    X.viscosity=1.68e12;
	    X.is_solid=false;	  
	    X.MM =  molecular_weight_aer[i] / 1.e6;
	    X.accomodation_coefficient = accomodation_coefficient[i];
	    X.soap_ind = i;
	    X.soap_ind_aero = (i-N_start) * (nlayer-1+i_hydrophilic) + i;
	    X.is_generic=true;
	    X.is_ion=false;
	    X.is_solvent=true;
	    X.catalyzed_ph=false;	   
	    surrogate.push_back(X);
	  }
      }
  /*
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
	    //surrogate[ij].moligo=config.moligo;
	    
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
      }*/
}
  
void creation_species_ssh( model_config &config, vector<species>& surrogate, vector<string> species_list_aer,
			   double molecular_weight_aer[], double accomodation_coefficient[],
			   int aerosol_type[],
			   vector<string> species_smiles, double saturation_vapor_pressure[],
			   double enthalpy_vaporization[], 
			   double diffusion_coef[], double henry[], double t_ref[], double mass_density[],
			   vector<string> species_part,
			   int nlayer, int i_hydrophilic, bool compute_inorganic, int N_inert, int N_inorganic)
{
  int nsp = species_list_aer.size();
  // double alpha = 1.0; //0.01; // accommodation coefficient

  int with_ca=0;
  int with_co3=0;

  species BiA2D;
  BiA2D.name="BiA2D";
  BiA2D.is_inorganic_precursor=false;
  BiA2D.Psat_ref=4.43e-6; // Saturation vapor pressure at Tref (torr)
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
		  accomodation_coefficient,diffusion_coef,saturation_vapor_pressure,enthalpy_vaporization, henry, t_ref, mass_density, species_part,nlayer,i_hydrophilic,
		  N_inert, N_inorganic);

  // GENOA DELETE DETECTOR START //
  // NOTE: When running ./clean under genoa clean mode (-gcl=cleanall), the code block after "GENOA DELETE DETECTOR START" and before "GENOA DELETE DETECTOR END" will be removed.

  /* TOLexp SPECIES START */

  /* ==== A02000 ==== */ 

  species A02000;
  A02000.name="A02000";
  A02000.is_inorganic_precursor=false;
  A02000.Psat_ref=12.045; 		// Saturation vapor pressure at Tref (torr)
  A02000.Tref=298;         		// Temperature of reference (K)
  A02000.deltaH=45.00;     		// Enthalpy of vaporization (kJ/mol)
  A02000.Henry=0.0;     		// Henry's law constant at Tref (M/atm)
  A02000.aq_type="none"; 		// "none","diacid","monoacid" or "aldehyde"
  // A02000.Kacidity1=XXacidityXX;   		// First acidity constant
  A02000.hydrophilic=true;   	// Does the species condense on the aqueous phase?
  A02000.hydrophobic=true;  		// Does the species condense on the organic phase?
  A02000.nonvolatile=false; 		// Is the compound nonvolatile?
  A02000.kp_from_experiment=false; 	// Use experimental partitioning constant at Tref?
  A02000.is_organic=true;  		// Is the compound organic?
  A02000.compute_gamma_org=true;  	// Compute the activity coefficients of the organic phase for this compound?
  A02000.compute_gamma_aq=true;  	// Compute the activity coefficients of the aqueous phase for this compound?
  A02000.Koligo_org=0.0;      	//oligomeriation constant in the organic phase
  A02000.rho=1300.;  
  A02000.is_monomer=false;
  A02000.rion=false;
  A02000.KDiffusion_air=1.0e-5;
  //  A02000.accomodation_coefficient=alpha;
  A02000.viscosity=1.68e12;
  A02000.is_solid=false;
  A02000.is_generic=false;
  
  //Group: if no functionnal group in the species use the default species
  //for the computation of activity coefficients
  //
  double group_tmp_A02000 [] = {1.,0.,0.,0., // group C
				  0.,0.,0.,0., //group C[OH]
				  0.,0.,0.,0., //group Calcohol
				  0.,0.,0.,0., //group Calcohol-tail
				  0.,0.,0.,0.,0., //group C=C
				  0.,0., //group aromatic carbon (AC)
				  0.,0.,0., // group //AC-C
				  0.,  //group OH
				  0., //group H2O
				  0., //group ACOH
				  0.,0., //group ketone
				  0.,   //group aldehyde  
				  0.,0., //group ester
				  0.,0.,0., //group ether 
				  1.,  //group acid
				  0.,   //group ACNO2
				  0.,0.,0., //group NO3
				  0.,0.,0., //group CO-OH
  			          0.,0.,0.,0.,0.,0.,0.,0.,0., //group CO-OC
			          0.,  //group PAN
			          0.,  //group CO-OOH
			          0.,  //group O=COC=O
			          0.,0.,0.}; //group CHxNO2
  
  size = sizeof(group_tmp_A02000)/sizeof(double);
  assert(size = 60);
	
  for(int i = 0; i < size; ++i)
	A02000.groups[i] = group_tmp_A02000[i];

  // Search the species name in the aerosol species list 
  // and add the species if its name matches with
  // the given list.
  add_species_ssh(surrogate, A02000, species_list_aer, molecular_weight_aer,
              accomodation_coefficient,diffusion_coef,saturation_vapor_pressure,enthalpy_vaporization, henry, t_ref, mass_density, species_part,nlayer,i_hydrophilic,
		  N_inert, N_inorganic);



/* ==== AA2000 ==== */ 

  species AA2000;
  AA2000.name="AA2000";
  AA2000.is_inorganic_precursor=false;
  AA2000.Psat_ref=0.004795; 		// Saturation vapor pressure at Tref (torr)
  AA2000.Tref=298;         		// Temperature of reference (K)
  AA2000.deltaH=76.63;     		// Enthalpy of vaporization (kJ/mol)
  AA2000.Henry=0.0;     		// Henry's law constant at Tref (M/atm)
  AA2000.aq_type="none"; 		// "none","diacid","monoacid" or "aldehyde"
  // AA2000.Kacidity1=XXacidityXX;   		// First acidity constant
  AA2000.hydrophilic=true;   	// Does the species condense on the aqueous phase?
  AA2000.hydrophobic=true;  		// Does the species condense on the organic phase?
  AA2000.nonvolatile=false; 		// Is the compound nonvolatile?
  AA2000.kp_from_experiment=false; 	// Use experimental partitioning constant at Tref?
  AA2000.is_organic=true;  		// Is the compound organic?
  AA2000.compute_gamma_org=true;  	// Compute the activity coefficients of the organic phase for this compound?
  AA2000.compute_gamma_aq=true;  	// Compute the activity coefficients of the aqueous phase for this compound?
  AA2000.Koligo_org=0.0;      	//oligomeriation constant in the organic phase
  AA2000.rho=1300.;  
  AA2000.is_monomer=false;
  AA2000.rion=false;
  AA2000.KDiffusion_air=1.0e-5;
  //  AA2000.accomodation_coefficient=alpha;
  AA2000.viscosity=1.68e12;
  AA2000.is_solid=false;
  AA2000.is_generic=false;
  
  //Group: if no functionnal group in the species use the default species
  //for the computation of activity coefficients
  //
  double group_tmp_AA2000 [] = {0.,0.,0.,0., // group C
				  0.,0.,0.,0., //group C[OH]
				  0.,0.,0.,0., //group Calcohol
				  0.,0.,0.,0., //group Calcohol-tail
				  0.,0.,0.,0.,0., //group C=C
				  0.,0., //group aromatic carbon (AC)
				  0.,0.,0., // group //AC-C
				  0.,  //group OH
				  0., //group H2O
				  0., //group ACOH
				  0.,0., //group ketone
				  0.,   //group aldehyde  
				  0.,0., //group ester
				  0.,0.,0., //group ether 
				  2.,  //group acid
				  0.,   //group ACNO2
				  0.,0.,0., //group NO3
				  0.,0.,0., //group CO-OH
  			          0.,0.,0.,0.,0.,0.,0.,0.,0., //group CO-OC
			          0.,  //group PAN
			          0.,  //group CO-OOH
			          0.,  //group O=COC=O
			          0.,0.,0.}; //group CHxNO2
  
  size = sizeof(group_tmp_AA2000)/sizeof(double);
  assert(size = 60);
	
  for(int i = 0; i < size; ++i)
	AA2000.groups[i] = group_tmp_AA2000[i];

  // Search the species name in the aerosol species list 
  // and add the species if its name matches with
  // the given list.
  add_species_ssh(surrogate, AA2000, species_list_aer, molecular_weight_aer,
              accomodation_coefficient,diffusion_coef,saturation_vapor_pressure,enthalpy_vaporization, henry, t_ref, mass_density, species_part,nlayer,i_hydrophilic,
		  N_inert, N_inorganic);



/* ==== AA4000 ==== */ 

  species AA4000;
  AA4000.name="AA4000";
  AA4000.is_inorganic_precursor=false;
  AA4000.Psat_ref=3.025e-11; 		// Saturation vapor pressure at Tref (torr)
  AA4000.Tref=298;         		// Temperature of reference (K)
  AA4000.deltaH=165.93;     		// Enthalpy of vaporization (kJ/mol)
  AA4000.Henry=0.0;     		// Henry's law constant at Tref (M/atm)
  AA4000.aq_type="none"; 		// "none","diacid","monoacid" or "aldehyde"
  // AA4000.Kacidity1=XXacidityXX;   		// First acidity constant
  AA4000.hydrophilic=true;   	// Does the species condense on the aqueous phase?
  AA4000.hydrophobic=true;  		// Does the species condense on the organic phase?
  AA4000.nonvolatile=false; 		// Is the compound nonvolatile?
  AA4000.kp_from_experiment=false; 	// Use experimental partitioning constant at Tref?
  AA4000.is_organic=true;  		// Is the compound organic?
  AA4000.compute_gamma_org=true;  	// Compute the activity coefficients of the organic phase for this compound?
  AA4000.compute_gamma_aq=true;  	// Compute the activity coefficients of the aqueous phase for this compound?
  AA4000.Koligo_org=0.0;      	//oligomeriation constant in the organic phase
  AA4000.rho=1300.;  
  AA4000.is_monomer=false;
  AA4000.rion=false;
  AA4000.KDiffusion_air=1.0e-5;
  //  AA4000.accomodation_coefficient=alpha;
  AA4000.viscosity=1.68e12;
  AA4000.is_solid=false;
  AA4000.is_generic=false;
  
  //Group: if no functionnal group in the species use the default species
  //for the computation of activity coefficients
  //
  double group_tmp_AA4000 [] = {0.,0.,0.,0., // group C
				  0.,0.,0.,0., //group C[OH]
				  0.,0.,0.,0., //group Calcohol
				  0.,0.,0.,0., //group Calcohol-tail
				  0.,0.,0.,0.,1., //group C=C
				  0.,0., //group aromatic carbon (AC)
				  0.,0.,0., // group //AC-C
				  2.,  //group OH
				  0., //group H2O
				  0., //group ACOH
				  0.,0., //group ketone
				  0.,   //group aldehyde  
				  0.,0., //group ester
				  0.,0.,0., //group ether 
				  2.,  //group acid
				  0.,   //group ACNO2
				  0.,0.,0., //group NO3
				  0.,0.,0., //group CO-OH
  			          0.,0.,0.,0.,0.,0.,0.,0.,0., //group CO-OC
			          0.,  //group PAN
			          0.,  //group CO-OOH
			          0.,  //group O=COC=O
			          0.,0.,0.}; //group CHxNO2
  
  size = sizeof(group_tmp_AA4000)/sizeof(double);
  assert(size = 60);
	
  for(int i = 0; i < size; ++i)
	AA4000.groups[i] = group_tmp_AA4000[i];

  // Search the species name in the aerosol species list 
  // and add the species if its name matches with
  // the given list.
  add_species_ssh(surrogate, AA4000, species_list_aer, molecular_weight_aer,
              accomodation_coefficient,diffusion_coef,saturation_vapor_pressure,enthalpy_vaporization, henry, t_ref, mass_density, species_part,nlayer,i_hydrophilic,
		  N_inert, N_inorganic);



/* ==== AD2000 ==== */ 

  species AD2000;
  AD2000.name="AD2000";
  AD2000.is_inorganic_precursor=false;
  AD2000.Psat_ref=0.095678; 		// Saturation vapor pressure at Tref (torr)
  AD2000.Tref=298;         		// Temperature of reference (K)
  AD2000.deltaH=62.47;     		// Enthalpy of vaporization (kJ/mol)
  AD2000.Henry=0.0;     		// Henry's law constant at Tref (M/atm)
  AD2000.aq_type="none"; 		// "none","diacid","monoacid" or "aldehyde"
  // AD2000.Kacidity1=XXacidityXX;   		// First acidity constant
  AD2000.hydrophilic=true;   	// Does the species condense on the aqueous phase?
  AD2000.hydrophobic=true;  		// Does the species condense on the organic phase?
  AD2000.nonvolatile=false; 		// Is the compound nonvolatile?
  AD2000.kp_from_experiment=false; 	// Use experimental partitioning constant at Tref?
  AD2000.is_organic=true;  		// Is the compound organic?
  AD2000.compute_gamma_org=true;  	// Compute the activity coefficients of the organic phase for this compound?
  AD2000.compute_gamma_aq=true;  	// Compute the activity coefficients of the aqueous phase for this compound?
  AD2000.Koligo_org=0.0;      	//oligomeriation constant in the organic phase
  AD2000.rho=1300.;  
  AD2000.is_monomer=false;
  AD2000.rion=false;
  AD2000.KDiffusion_air=1.0e-5;
  //  AD2000.accomodation_coefficient=alpha;
  AD2000.viscosity=1.68e12;
  AD2000.is_solid=false;
  AD2000.is_generic=false;
  
  //Group: if no functionnal group in the species use the default species
  //for the computation of activity coefficients
  //
  double group_tmp_AD2000 [] = {0.,0.,0.,0., // group C
				  0.,0.,0.,0., //group C[OH]
				  0.,0.,0.,0., //group Calcohol
				  0.,0.,0.,0., //group Calcohol-tail
				  0.,0.,0.,0.,0., //group C=C
				  0.,0., //group aromatic carbon (AC)
				  0.,0.,0., // group //AC-C
				  0.,  //group OH
				  0., //group H2O
				  0., //group ACOH
				  0.,0., //group ketone
				  1.,   //group aldehyde  
				  0.,0., //group ester
				  0.,0.,0., //group ether 
				  1.,  //group acid
				  0.,   //group ACNO2
				  0.,0.,0., //group NO3
				  0.,0.,0., //group CO-OH
  			          0.,0.,0.,0.,0.,0.,0.,0.,0., //group CO-OC
			          0.,  //group PAN
			          0.,  //group CO-OOH
			          0.,  //group O=COC=O
			          0.,0.,0.}; //group CHxNO2
  
  size = sizeof(group_tmp_AD2000)/sizeof(double);
  assert(size = 60);
	
  for(int i = 0; i < size; ++i)
	AD2000.groups[i] = group_tmp_AD2000[i];

  // Search the species name in the aerosol species list 
  // and add the species if its name matches with
  // the given list.
  add_species_ssh(surrogate, AD2000, species_list_aer, molecular_weight_aer,
              accomodation_coefficient,diffusion_coef,saturation_vapor_pressure,enthalpy_vaporization, henry, t_ref, mass_density, species_part,nlayer,i_hydrophilic,
		  N_inert, N_inorganic);



/* ==== AD4000 ==== */ 

  species AD4000;
  AD4000.name="AD4000";
  AD4000.is_inorganic_precursor=false;
  AD4000.Psat_ref=7.6e-05; 		// Saturation vapor pressure at Tref (torr)
  AD4000.Tref=298;         		// Temperature of reference (K)
  AD4000.deltaH=96.23;     		// Enthalpy of vaporization (kJ/mol)
  AD4000.Henry=0.0;     		// Henry's law constant at Tref (M/atm)
  AD4000.aq_type="none"; 		// "none","diacid","monoacid" or "aldehyde"
  // AD4000.Kacidity1=XXacidityXX;   		// First acidity constant
  AD4000.hydrophilic=true;   	// Does the species condense on the aqueous phase?
  AD4000.hydrophobic=true;  		// Does the species condense on the organic phase?
  AD4000.nonvolatile=false; 		// Is the compound nonvolatile?
  AD4000.kp_from_experiment=false; 	// Use experimental partitioning constant at Tref?
  AD4000.is_organic=true;  		// Is the compound organic?
  AD4000.compute_gamma_org=true;  	// Compute the activity coefficients of the organic phase for this compound?
  AD4000.compute_gamma_aq=true;  	// Compute the activity coefficients of the aqueous phase for this compound?
  AD4000.Koligo_org=0.0;      	//oligomeriation constant in the organic phase
  AD4000.rho=1300.;  
  AD4000.is_monomer=false;
  AD4000.rion=false;
  AD4000.KDiffusion_air=1.0e-5;
  //  AD4000.accomodation_coefficient=alpha;
  AD4000.viscosity=1.68e12;
  AD4000.is_solid=false;
  AD4000.is_generic=false;
  
  //Group: if no functionnal group in the species use the default species
  //for the computation of activity coefficients
  //
  double group_tmp_AD4000 [] = {0.,0.,0.,0., // group C
				  0.,0.,0.,0., //group C[OH]
				  0.,0.,0.,0., //group Calcohol
				  0.,0.,0.,0., //group Calcohol-tail
				  0.,0.,0.,1.,0., //group C=C
				  0.,0., //group aromatic carbon (AC)
				  0.,0.,0., // group //AC-C
				  1.,  //group OH
				  0., //group H2O
				  0., //group ACOH
				  0.,0., //group ketone
				  1.,   //group aldehyde  
				  0.,0., //group ester
				  0.,0.,0., //group ether 
				  1.,  //group acid
				  0.,   //group ACNO2
				  0.,0.,0., //group NO3
				  0.,0.,0., //group CO-OH
  			          0.,0.,0.,0.,0.,0.,0.,0.,0., //group CO-OC
			          0.,  //group PAN
			          0.,  //group CO-OOH
			          0.,  //group O=COC=O
			          0.,0.,0.}; //group CHxNO2
  
  size = sizeof(group_tmp_AD4000)/sizeof(double);
  assert(size = 60);
	
  for(int i = 0; i < size; ++i)
	AD4000.groups[i] = group_tmp_AD4000[i];

  // Search the species name in the aerosol species list 
  // and add the species if its name matches with
  // the given list.
  add_species_ssh(surrogate, AD4000, species_list_aer, molecular_weight_aer,
              accomodation_coefficient,diffusion_coef,saturation_vapor_pressure,enthalpy_vaporization, henry, t_ref, mass_density, species_part,nlayer,i_hydrophilic,
		  N_inert, N_inorganic);



/* ==== AD4001 ==== */ 

  species AD4001;
  AD4001.name="AD4001";
  AD4001.is_inorganic_precursor=false;
  AD4001.Psat_ref=9.567e-8; 		// Saturation vapor pressure at Tref (torr)
  AD4001.Tref=298;         		// Temperature of reference (K)
  AD4001.deltaH=127.8;     		// Enthalpy of vaporization (kJ/mol)
  AD4001.Henry=0.0;     		// Henry's law constant at Tref (M/atm)
  AD4001.aq_type="none"; 		// "none","diacid","monoacid" or "aldehyde"
  // AD4001.Kacidity1=XXacidityXX;   		// First acidity constant
  AD4001.hydrophilic=true;   	// Does the species condense on the aqueous phase?
  AD4001.hydrophobic=true;  		// Does the species condense on the organic phase?
  AD4001.nonvolatile=false; 		// Is the compound nonvolatile?
  AD4001.kp_from_experiment=false; 	// Use experimental partitioning constant at Tref?
  AD4001.is_organic=true;  		// Is the compound organic?
  AD4001.compute_gamma_org=true;  	// Compute the activity coefficients of the organic phase for this compound?
  AD4001.compute_gamma_aq=true;  	// Compute the activity coefficients of the aqueous phase for this compound?
  AD4001.Koligo_org=0.0;      	//oligomeriation constant in the organic phase
  AD4001.rho=1300.;  
  AD4001.is_monomer=false;
  AD4001.rion=false;
  AD4001.KDiffusion_air=1.0e-5;
  //  AD4001.accomodation_coefficient=alpha;
  AD4001.viscosity=1.68e12;
  AD4001.is_solid=false;
  AD4001.is_generic=false;
  
  //Group: if no functionnal group in the species use the default species
  //for the computation of activity coefficients
  //
  double group_tmp_AD4001 [] = {0.,0.,0.,0., // group C
				  0.,0.,0.,0., //group C[OH]
				  0.,0.,0.,0., //group Calcohol
				  0.,0.,0.,0., //group Calcohol-tail
				  0.,0.,0.,0.,1., //group C=C
				  0.,0., //group aromatic carbon (AC)
				  0.,0.,0., // group //AC-C
				  2.,  //group OH
				  0., //group H2O
				  0., //group ACOH
				  0.,0., //group ketone
				  1.,   //group aldehyde  
				  0.,0., //group ester
				  0.,0.,0., //group ether 
				  1.,  //group acid
				  0.,   //group ACNO2
				  0.,0.,0., //group NO3
				  0.,0.,0., //group CO-OH
  			          0.,0.,0.,0.,0.,0.,0.,0.,0., //group CO-OC
			          0.,  //group PAN
			          0.,  //group CO-OOH
			          0.,  //group O=COC=O
			          0.,0.,0.}; //group CHxNO2
  
  size = sizeof(group_tmp_AD4001)/sizeof(double);
  assert(size = 60);
	
  for(int i = 0; i < size; ++i)
	AD4001.groups[i] = group_tmp_AD4001[i];

  // Search the species name in the aerosol species list 
  // and add the species if its name matches with
  // the given list.
  add_species_ssh(surrogate, AD4001, species_list_aer, molecular_weight_aer,
              accomodation_coefficient,diffusion_coef,saturation_vapor_pressure,enthalpy_vaporization, henry, t_ref, mass_density, species_part,nlayer,i_hydrophilic,
		  N_inert, N_inorganic);



/* ==== AD5000 ==== */ 

  species AD5000;
  AD5000.name="AD5000";
  AD5000.is_inorganic_precursor=false;
  AD5000.Psat_ref=0.000190; 		// Saturation vapor pressure at Tref (torr)
  AD5000.Tref=298;         		// Temperature of reference (K)
  AD5000.deltaH=91.87;     		// Enthalpy of vaporization (kJ/mol)
  AD5000.Henry=0.0;     		// Henry's law constant at Tref (M/atm)
  AD5000.aq_type="none"; 		// "none","diacid","monoacid" or "aldehyde"
  // AD5000.Kacidity1=XXacidityXX;   		// First acidity constant
  AD5000.hydrophilic=true;   	// Does the species condense on the aqueous phase?
  AD5000.hydrophobic=true;  		// Does the species condense on the organic phase?
  AD5000.nonvolatile=false; 		// Is the compound nonvolatile?
  AD5000.kp_from_experiment=false; 	// Use experimental partitioning constant at Tref?
  AD5000.is_organic=true;  		// Is the compound organic?
  AD5000.compute_gamma_org=true;  	// Compute the activity coefficients of the organic phase for this compound?
  AD5000.compute_gamma_aq=true;  	// Compute the activity coefficients of the aqueous phase for this compound?
  AD5000.Koligo_org=0.0;      	//oligomeriation constant in the organic phase
  AD5000.rho=1300.;  
  AD5000.is_monomer=false;
  AD5000.rion=false;
  AD5000.KDiffusion_air=1.0e-5;
  //  AD5000.accomodation_coefficient=alpha;
  AD5000.viscosity=1.68e12;
  AD5000.is_solid=false;
  AD5000.is_generic=false;
  
  //Group: if no functionnal group in the species use the default species
  //for the computation of activity coefficients
  //
  double group_tmp_AD5000 [] = {1.,0.,0.,0., // group C
				  0.,0.,0.,0., //group C[OH]
				  0.,0.,0.,0., //group Calcohol
				  0.,0.,0.,0., //group Calcohol-tail
				  0.,0.,0.,0.,1., //group C=C
				  0.,0., //group aromatic carbon (AC)
				  0.,0.,0., // group //AC-C
				  1.,  //group OH
				  0., //group H2O
				  0., //group ACOH
				  0.,0., //group ketone
				  1.,   //group aldehyde  
				  0.,0., //group ester
				  0.,0.,0., //group ether 
				  1.,  //group acid
				  0.,   //group ACNO2
				  0.,0.,0., //group NO3
				  0.,0.,0., //group CO-OH
  			          0.,0.,0.,0.,0.,0.,0.,0.,0., //group CO-OC
			          0.,  //group PAN
			          0.,  //group CO-OOH
			          0.,  //group O=COC=O
			          0.,0.,0.}; //group CHxNO2
  
  size = sizeof(group_tmp_AD5000)/sizeof(double);
  assert(size = 60);
	
  for(int i = 0; i < size; ++i)
	AD5000.groups[i] = group_tmp_AD5000[i];

  // Search the species name in the aerosol species list 
  // and add the species if its name matches with
  // the given list.
  add_species_ssh(surrogate, AD5000, species_list_aer, molecular_weight_aer,
              accomodation_coefficient,diffusion_coef,saturation_vapor_pressure,enthalpy_vaporization, henry, t_ref, mass_density, species_part,nlayer,i_hydrophilic,
		  N_inert, N_inorganic);



/* ==== AK3000 ==== */ 

  species AK3000;
  AK3000.name="AK3000";
  AK3000.is_inorganic_precursor=false;
  AK3000.Psat_ref=3.145e-1; 		// Saturation vapor pressure at Tref (torr)
  AK3000.Tref=298;         		// Temperature of reference (K)
  AK3000.deltaH=59.32;     		// Enthalpy of vaporization (kJ/mol)
  AK3000.Henry=0.0;     		// Henry's law constant at Tref (M/atm)
  AK3000.aq_type="none"; 		// "none","diacid","monoacid" or "aldehyde"
  // AK3000.Kacidity1=XXacidityXX;   		// First acidity constant
  AK3000.hydrophilic=true;   	// Does the species condense on the aqueous phase?
  AK3000.hydrophobic=true;  		// Does the species condense on the organic phase?
  AK3000.nonvolatile=false; 		// Is the compound nonvolatile?
  AK3000.kp_from_experiment=false; 	// Use experimental partitioning constant at Tref?
  AK3000.is_organic=true;  		// Is the compound organic?
  AK3000.compute_gamma_org=true;  	// Compute the activity coefficients of the organic phase for this compound?
  AK3000.compute_gamma_aq=true;  	// Compute the activity coefficients of the aqueous phase for this compound?
  AK3000.Koligo_org=0.0;      	//oligomeriation constant in the organic phase
  AK3000.rho=1300.;  
  AK3000.is_monomer=false;
  AK3000.rion=false;
  AK3000.KDiffusion_air=1.0e-5;
  //  AK3000.accomodation_coefficient=alpha;
  AK3000.viscosity=1.68e12;
  AK3000.is_solid=false;
  AK3000.is_generic=false;
  
  //Group: if no functionnal group in the species use the default species
  //for the computation of activity coefficients
  //
  double group_tmp_AK3000 [] = {0.,0.,0.,0., // group C
				  0.,0.,0.,0., //group C[OH]
				  0.,0.,0.,0., //group Calcohol
				  0.,0.,0.,0., //group Calcohol-tail
				  0.,0.,0.,0.,0., //group C=C
				  0.,0., //group aromatic carbon (AC)
				  0.,0.,0., // group //AC-C
				  0.,  //group OH
				  0., //group H2O
				  0., //group ACOH
				  1.,0., //group ketone
				  0.,   //group aldehyde  
				  0.,0., //group ester
				  0.,0.,0., //group ether 
				  1.,  //group acid
				  0.,   //group ACNO2
				  0.,0.,0., //group NO3
				  0.,0.,0., //group CO-OH
  			          0.,0.,0.,0.,0.,0.,0.,0.,0., //group CO-OC
			          0.,  //group PAN
			          0.,  //group CO-OOH
			          0.,  //group O=COC=O
			          0.,0.,0.}; //group CHxNO2
  
  size = sizeof(group_tmp_AK3000)/sizeof(double);
  assert(size = 60);
	
  for(int i = 0; i < size; ++i)
	AK3000.groups[i] = group_tmp_AK3000[i];

  // Search the species name in the aerosol species list 
  // and add the species if its name matches with
  // the given list.
  add_species_ssh(surrogate, AK3000, species_list_aer, molecular_weight_aer,
              accomodation_coefficient,diffusion_coef,saturation_vapor_pressure,enthalpy_vaporization, henry, t_ref, mass_density, species_part,nlayer,i_hydrophilic,
		  N_inert, N_inorganic);



/* ==== AK5000 ==== */ 

  species AK5000;
  AK5000.name="AK5000";
  AK5000.is_inorganic_precursor=false;
  AK5000.Psat_ref=6.036e-7; 		// Saturation vapor pressure at Tref (torr)
  AK5000.Tref=298;         		// Temperature of reference (K)
  AK5000.deltaH=119.1;     		// Enthalpy of vaporization (kJ/mol)
  AK5000.Henry=0.0;     		// Henry's law constant at Tref (M/atm)
  AK5000.aq_type="none"; 		// "none","diacid","monoacid" or "aldehyde"
  // AK5000.Kacidity1=XXacidityXX;   		// First acidity constant
  AK5000.hydrophilic=true;   	// Does the species condense on the aqueous phase?
  AK5000.hydrophobic=true;  		// Does the species condense on the organic phase?
  AK5000.nonvolatile=false; 		// Is the compound nonvolatile?
  AK5000.kp_from_experiment=false; 	// Use experimental partitioning constant at Tref?
  AK5000.is_organic=true;  		// Is the compound organic?
  AK5000.compute_gamma_org=true;  	// Compute the activity coefficients of the organic phase for this compound?
  AK5000.compute_gamma_aq=true;  	// Compute the activity coefficients of the aqueous phase for this compound?
  AK5000.Koligo_org=0.0;      	//oligomeriation constant in the organic phase
  AK5000.rho=1300.;  
  AK5000.is_monomer=false;
  AK5000.rion=false;
  AK5000.KDiffusion_air=1.0e-5;
  //  AK5000.accomodation_coefficient=alpha;
  AK5000.viscosity=1.68e12;
  AK5000.is_solid=false;
  AK5000.is_generic=false;
  
  //Group: if no functionnal group in the species use the default species
  //for the computation of activity coefficients
  //
  double group_tmp_AK5000 [] = {0.,0.,0.,0., // group C
				  0.,0.,0.,0., //group C[OH]
				  0.,0.,0.,0., //group Calcohol
				  0.,0.,0.,0., //group Calcohol-tail
				  0.,0.,0.,0.,1., //group C=C
				  0.,0., //group aromatic carbon (AC)
				  0.,0.,0., // group //AC-C
				  2.,  //group OH
				  0., //group H2O
				  0., //group ACOH
				  1.,0., //group ketone
				  0.,   //group aldehyde  
				  0.,0., //group ester
				  0.,0.,0., //group ether 
				  1.,  //group acid
				  0.,   //group ACNO2
				  0.,0.,0., //group NO3
				  0.,0.,0., //group CO-OH
  			          0.,0.,0.,0.,0.,0.,0.,0.,0., //group CO-OC
			          0.,  //group PAN
			          0.,  //group CO-OOH
			          0.,  //group O=COC=O
			          0.,0.,0.}; //group CHxNO2
  
  size = sizeof(group_tmp_AK5000)/sizeof(double);
  assert(size = 60);
	
  for(int i = 0; i < size; ++i)
	AK5000.groups[i] = group_tmp_AK5000[i];

  // Search the species name in the aerosol species list 
  // and add the species if its name matches with
  // the given list.
  add_species_ssh(surrogate, AK5000, species_list_aer, molecular_weight_aer,
              accomodation_coefficient,diffusion_coef,saturation_vapor_pressure,enthalpy_vaporization, henry, t_ref, mass_density, species_part,nlayer,i_hydrophilic,
		  N_inert, N_inorganic);



/* ==== AR0010 ==== */ 

  species AR0010;
  AR0010.name="AR0010";
  AR0010.is_inorganic_precursor=false;
  AR0010.Psat_ref=8.423; 		// Saturation vapor pressure at Tref (torr)
  AR0010.Tref=298;         		// Temperature of reference (K)
  AR0010.deltaH=45.88;     		// Enthalpy of vaporization (kJ/mol)
  AR0010.Henry=0.0;     		// Henry's law constant at Tref (M/atm)
  AR0010.aq_type="none"; 		// "none","diacid","monoacid" or "aldehyde"
  // AR0010.Kacidity1=XXacidityXX;   		// First acidity constant
  AR0010.hydrophilic=true;   	// Does the species condense on the aqueous phase?
  AR0010.hydrophobic=true;  		// Does the species condense on the organic phase?
  AR0010.nonvolatile=false; 		// Is the compound nonvolatile?
  AR0010.kp_from_experiment=false; 	// Use experimental partitioning constant at Tref?
  AR0010.is_organic=true;  		// Is the compound organic?
  AR0010.compute_gamma_org=true;  	// Compute the activity coefficients of the organic phase for this compound?
  AR0010.compute_gamma_aq=true;  	// Compute the activity coefficients of the aqueous phase for this compound?
  AR0010.Koligo_org=0.0;      	//oligomeriation constant in the organic phase
  AR0010.rho=1300.;  
  AR0010.is_monomer=false;
  AR0010.rion=false;
  AR0010.KDiffusion_air=1.0e-5;
  //  AR0010.accomodation_coefficient=alpha;
  AR0010.viscosity=1.68e12;
  AR0010.is_solid=false;
  AR0010.is_generic=false;
  
  //Group: if no functionnal group in the species use the default species
  //for the computation of activity coefficients
  //
  double group_tmp_AR0010 [] = {0.,0.,1.,0., // group C
				  0.,0.,0.,0., //group C[OH]
				  0.,0.,0.,0., //group Calcohol
				  0.,0.,0.,0., //group Calcohol-tail
				  0.,0.,0.,0.,0., //group C=C
				  0.,0., //group aromatic carbon (AC)
				  0.,0.,0., // group //AC-C
				  0.,  //group OH
				  0., //group H2O
				  0., //group ACOH
				  0.,0., //group ketone
				  2.,   //group aldehyde  
				  0.,0., //group ester
				  0.,0.,1., //group ether 
				  0.,  //group acid
				  0.,   //group ACNO2
				  0.,0.,0., //group NO3
				  0.,0.,0., //group CO-OH
  			          0.,0.,0.,0.,0.,0.,0.,0.,0., //group CO-OC
			          0.,  //group PAN
			          0.,  //group CO-OOH
			          0.,  //group O=COC=O
			          0.,0.,0.}; //group CHxNO2
  
  size = sizeof(group_tmp_AR0010)/sizeof(double);
  assert(size = 60);
	
  for(int i = 0; i < size; ++i)
	AR0010.groups[i] = group_tmp_AR0010[i];

  // Search the species name in the aerosol species list 
  // and add the species if its name matches with
  // the given list.
  add_species_ssh(surrogate, AR0010, species_list_aer, molecular_weight_aer,
              accomodation_coefficient,diffusion_coef,saturation_vapor_pressure,enthalpy_vaporization, henry, t_ref, mass_density, species_part,nlayer,i_hydrophilic,
		  N_inert, N_inorganic);



/* ==== AR0027 ==== */ 

  species AR0027;
  AR0027.name="AR0027";
  AR0027.is_inorganic_precursor=false;
  AR0027.Psat_ref=8.057e-3; 		// Saturation vapor pressure at Tref (torr)
  AR0027.Tref=298;         		// Temperature of reference (K)
  AR0027.deltaH=72.11;     		// Enthalpy of vaporization (kJ/mol)
  AR0027.Henry=0.0;     		// Henry's law constant at Tref (M/atm)
  AR0027.aq_type="none"; 		// "none","diacid","monoacid" or "aldehyde"
  // AR0027.Kacidity1=XXacidityXX;   		// First acidity constant
  AR0027.hydrophilic=true;   	// Does the species condense on the aqueous phase?
  AR0027.hydrophobic=true;  		// Does the species condense on the organic phase?
  AR0027.nonvolatile=false; 		// Is the compound nonvolatile?
  AR0027.kp_from_experiment=false; 	// Use experimental partitioning constant at Tref?
  AR0027.is_organic=true;  		// Is the compound organic?
  AR0027.compute_gamma_org=true;  	// Compute the activity coefficients of the organic phase for this compound?
  AR0027.compute_gamma_aq=true;  	// Compute the activity coefficients of the aqueous phase for this compound?
  AR0027.Koligo_org=0.0;      	//oligomeriation constant in the organic phase
  AR0027.rho=1300.;  
  AR0027.is_monomer=false;
  AR0027.rion=false;
  AR0027.KDiffusion_air=1.0e-5;
  //  AR0027.accomodation_coefficient=alpha;
  AR0027.viscosity=1.68e12;
  AR0027.is_solid=false;
  AR0027.is_generic=false;
  
  //Group: if no functionnal group in the species use the default species
  //for the computation of activity coefficients
  //
  double group_tmp_AR0027 [] = {0.,0.,0.,0., // group C
				  0.,0.,0.,0., //group C[OH]
				  0.,0.,0.,0., //group Calcohol
				  0.,0.,0.,0., //group Calcohol-tail
				  0.,0.,0.,0.,0., //group C=C
				  4.,0., //group aromatic carbon (AC)
				  0.,0.,0., // group //AC-C
				  0.,  //group OH
				  0., //group H2O
				  1., //group ACOH
				  0.,0., //group ketone
				  0.,   //group aldehyde  
				  0.,0., //group ester
				  0.,0.,0., //group ether 
				  0.,  //group acid
				  1.,   //group ACNO2
				  0.,0.,0., //group NO3
				  0.,0.,0., //group CO-OH
  			          0.,0.,0.,0.,0.,0.,0.,0.,0., //group CO-OC
			          0.,  //group PAN
			          0.,  //group CO-OOH
			          0.,  //group O=COC=O
			          0.,0.,0.}; //group CHxNO2
  
  size = sizeof(group_tmp_AR0027)/sizeof(double);
  assert(size = 60);
	
  for(int i = 0; i < size; ++i)
	AR0027.groups[i] = group_tmp_AR0027[i];

  // Search the species name in the aerosol species list 
  // and add the species if its name matches with
  // the given list.
  add_species_ssh(surrogate, AR0027, species_list_aer, molecular_weight_aer,
              accomodation_coefficient,diffusion_coef,saturation_vapor_pressure,enthalpy_vaporization, henry, t_ref, mass_density, species_part,nlayer,i_hydrophilic,
		  N_inert, N_inorganic);



/* ==== AR0043 ==== */ 

  species AR0043;
  AR0043.name="AR0043";
  AR0043.is_inorganic_precursor=false;
  AR0043.Psat_ref=5.299e-2; 		// Saturation vapor pressure at Tref (torr)
  AR0043.Tref=298;         		// Temperature of reference (K)
  AR0043.deltaH=66.03;     		// Enthalpy of vaporization (kJ/mol)
  AR0043.Henry=0.0;     		// Henry's law constant at Tref (M/atm)
  AR0043.aq_type="none"; 		// "none","diacid","monoacid" or "aldehyde"
  // AR0043.Kacidity1=XXacidityXX;   		// First acidity constant
  AR0043.hydrophilic=true;   	// Does the species condense on the aqueous phase?
  AR0043.hydrophobic=true;  		// Does the species condense on the organic phase?
  AR0043.nonvolatile=false; 		// Is the compound nonvolatile?
  AR0043.kp_from_experiment=false; 	// Use experimental partitioning constant at Tref?
  AR0043.is_organic=true;  		// Is the compound organic?
  AR0043.compute_gamma_org=true;  	// Compute the activity coefficients of the organic phase for this compound?
  AR0043.compute_gamma_aq=true;  	// Compute the activity coefficients of the aqueous phase for this compound?
  AR0043.Koligo_org=0.0;      	//oligomeriation constant in the organic phase
  AR0043.rho=1300.;  
  AR0043.is_monomer=false;
  AR0043.rion=false;
  AR0043.KDiffusion_air=1.0e-5;
  //  AR0043.accomodation_coefficient=alpha;
  AR0043.viscosity=1.68e12;
  AR0043.is_solid=false;
  AR0043.is_generic=false;
  
  //Group: if no functionnal group in the species use the default species
  //for the computation of activity coefficients
  //
  double group_tmp_AR0043 [] = {0.,0.,0.,0., // group C
				  0.,0.,0.,0., //group C[OH]
				  0.,0.,0.,0., //group Calcohol
				  0.,0.,0.,0., //group Calcohol-tail
				  0.,0.,0.,0.,0., //group C=C
				  5.,1., //group aromatic carbon (AC)
				  0.,0.,0., // group //AC-C
				  0.,  //group OH
				  0., //group H2O
				  0., //group ACOH
				  0.,0., //group ketone
				  0.,   //group aldehyde  
				  0.,0., //group ester
				  0.,0.,0., //group ether 
				  0.,  //group acid
				  0.,   //group ACNO2
				  0.,0.,0., //group NO3
				  1.,0.,0., //group CO-OH
  			          0.,0.,0.,0.,0.,0.,0.,0.,0., //group CO-OC
			          0.,  //group PAN
			          0.,  //group CO-OOH
			          0.,  //group O=COC=O
			          0.,0.,0.}; //group CHxNO2
  
  size = sizeof(group_tmp_AR0043)/sizeof(double);
  assert(size = 60);
	
  for(int i = 0; i < size; ++i)
	AR0043.groups[i] = group_tmp_AR0043[i];

  // Search the species name in the aerosol species list 
  // and add the species if its name matches with
  // the given list.
  add_species_ssh(surrogate, AR0043, species_list_aer, molecular_weight_aer,
              accomodation_coefficient,diffusion_coef,saturation_vapor_pressure,enthalpy_vaporization, henry, t_ref, mass_density, species_part,nlayer,i_hydrophilic,
		  N_inert, N_inorganic);



/* ==== AR0048 ==== */ 

  species AR0048;
  AR0048.name="AR0048";
  AR0048.is_inorganic_precursor=false;
  AR0048.Psat_ref=1.596e-5; 		// Saturation vapor pressure at Tref (torr)
  AR0048.Tref=298;         		// Temperature of reference (K)
  AR0048.deltaH=102.14;     		// Enthalpy of vaporization (kJ/mol)
  AR0048.Henry=0.0;     		// Henry's law constant at Tref (M/atm)
  AR0048.aq_type="none"; 		// "none","diacid","monoacid" or "aldehyde"
  // AR0048.Kacidity1=XXacidityXX;   		// First acidity constant
  AR0048.hydrophilic=true;   	// Does the species condense on the aqueous phase?
  AR0048.hydrophobic=true;  		// Does the species condense on the organic phase?
  AR0048.nonvolatile=false; 		// Is the compound nonvolatile?
  AR0048.kp_from_experiment=false; 	// Use experimental partitioning constant at Tref?
  AR0048.is_organic=true;  		// Is the compound organic?
  AR0048.compute_gamma_org=true;  	// Compute the activity coefficients of the organic phase for this compound?
  AR0048.compute_gamma_aq=true;  	// Compute the activity coefficients of the aqueous phase for this compound?
  AR0048.Koligo_org=0.0;      	//oligomeriation constant in the organic phase
  AR0048.rho=1300.;  
  AR0048.is_monomer=false;
  AR0048.rion=false;
  AR0048.KDiffusion_air=1.0e-5;
  //  AR0048.accomodation_coefficient=alpha;
  AR0048.viscosity=1.68e12;
  AR0048.is_solid=false;
  AR0048.is_generic=false;
  
  //Group: if no functionnal group in the species use the default species
  //for the computation of activity coefficients
  //
  double group_tmp_AR0048 [] = {0.,0.,0.,0., // group C
				  0.,0.,1.,0., //group C[OH]
				  0.,0.,0.,0., //group Calcohol
				  0.,0.,0.,0., //group Calcohol-tail
				  0.,0.,0.,0.,0., //group C=C
				  0.,0., //group aromatic carbon (AC)
				  0.,0.,0., // group //AC-C
				  1.,  //group OH
				  0., //group H2O
				  0., //group ACOH
				  0.,0., //group ketone
				  0.,   //group aldehyde  
				  0.,1., //group ester
				  0.,0.,0., //group ether 
				  0.,  //group acid
				  0.,   //group ACNO2
				  0.,0.,0., //group NO3
				  0.,1.,0., //group CO-OH
  			          0.,0.,0.,0.,0.,0.,0.,0.,0., //group CO-OC
			          0.,  //group PAN
			          0.,  //group CO-OOH
			          0.,  //group O=COC=O
			          0.,0.,0.}; //group CHxNO2
  
  size = sizeof(group_tmp_AR0048)/sizeof(double);
  assert(size = 60);
	
  for(int i = 0; i < size; ++i)
	AR0048.groups[i] = group_tmp_AR0048[i];

  // Search the species name in the aerosol species list 
  // and add the species if its name matches with
  // the given list.
  add_species_ssh(surrogate, AR0048, species_list_aer, molecular_weight_aer,
              accomodation_coefficient,diffusion_coef,saturation_vapor_pressure,enthalpy_vaporization, henry, t_ref, mass_density, species_part,nlayer,i_hydrophilic,
		  N_inert, N_inorganic);



/* ==== AR0087 ==== */ 

  species AR0087;
  AR0087.name="AR0087";
  AR0087.is_inorganic_precursor=false;
  AR0087.Psat_ref=3.221e-1; 		// Saturation vapor pressure at Tref (torr)
  AR0087.Tref=298;         		// Temperature of reference (K)
  AR0087.deltaH=58.10;     		// Enthalpy of vaporization (kJ/mol)
  AR0087.Henry=0.0;     		// Henry's law constant at Tref (M/atm)
  AR0087.aq_type="none"; 		// "none","diacid","monoacid" or "aldehyde"
  // AR0087.Kacidity1=XXacidityXX;   		// First acidity constant
  AR0087.hydrophilic=true;   	// Does the species condense on the aqueous phase?
  AR0087.hydrophobic=true;  		// Does the species condense on the organic phase?
  AR0087.nonvolatile=false; 		// Is the compound nonvolatile?
  AR0087.kp_from_experiment=false; 	// Use experimental partitioning constant at Tref?
  AR0087.is_organic=true;  		// Is the compound organic?
  AR0087.compute_gamma_org=true;  	// Compute the activity coefficients of the organic phase for this compound?
  AR0087.compute_gamma_aq=true;  	// Compute the activity coefficients of the aqueous phase for this compound?
  AR0087.Koligo_org=0.0;      	//oligomeriation constant in the organic phase
  AR0087.rho=1300.;  
  AR0087.is_monomer=false;
  AR0087.rion=false;
  AR0087.KDiffusion_air=1.0e-5;
  //  AR0087.accomodation_coefficient=alpha;
  AR0087.viscosity=1.68e12;
  AR0087.is_solid=false;
  AR0087.is_generic=false;
  
  //Group: if no functionnal group in the species use the default species
  //for the computation of activity coefficients
  //
  double group_tmp_AR0087 [] = {0.,0.,1.,0., // group C
				  0.,0.,0.,0., //group C[OH]
				  0.,0.,0.,0., //group Calcohol
				  0.,0.,0.,0., //group Calcohol-tail
				  0.,1.,0.,0.,0., //group C=C
				  0.,0., //group aromatic carbon (AC)
				  0.,0.,0., // group //AC-C
				  0.,  //group OH
				  0., //group H2O
				  0., //group ACOH
				  1.,0., //group ketone
				  1.,   //group aldehyde  
				  0.,0., //group ester
				  0.,0.,1., //group ether 
				  0.,  //group acid
				  0.,   //group ACNO2
				  0.,0.,0., //group NO3
				  0.,0.,0., //group CO-OH
  			          0.,0.,0.,0.,0.,0.,0.,0.,0., //group CO-OC
			          0.,  //group PAN
			          0.,  //group CO-OOH
			          0.,  //group O=COC=O
			          0.,0.,0.}; //group CHxNO2
  
  size = sizeof(group_tmp_AR0087)/sizeof(double);
  assert(size = 60);
	
  for(int i = 0; i < size; ++i)
	AR0087.groups[i] = group_tmp_AR0087[i];

  // Search the species name in the aerosol species list 
  // and add the species if its name matches with
  // the given list.
  add_species_ssh(surrogate, AR0087, species_list_aer, molecular_weight_aer,
              accomodation_coefficient,diffusion_coef,saturation_vapor_pressure,enthalpy_vaporization, henry, t_ref, mass_density, species_part,nlayer,i_hydrophilic,
		  N_inert, N_inorganic);



/* ==== AR0088 ==== */ 

  species AR0088;
  AR0088.name="AR0088";
  AR0088.is_inorganic_precursor=false;
  AR0088.Psat_ref=4.893e-1; 		// Saturation vapor pressure at Tref (torr)
  AR0088.Tref=298;         		// Temperature of reference (K)
  AR0088.deltaH=55.58;     		// Enthalpy of vaporization (kJ/mol)
  AR0088.Henry=0.0;     		// Henry's law constant at Tref (M/atm)
  AR0088.aq_type="none"; 		// "none","diacid","monoacid" or "aldehyde"
  // AR0088.Kacidity1=XXacidityXX;   		// First acidity constant
  AR0088.hydrophilic=true;   	// Does the species condense on the aqueous phase?
  AR0088.hydrophobic=true;  		// Does the species condense on the organic phase?
  AR0088.nonvolatile=false; 		// Is the compound nonvolatile?
  AR0088.kp_from_experiment=false; 	// Use experimental partitioning constant at Tref?
  AR0088.is_organic=true;  		// Is the compound organic?
  AR0088.compute_gamma_org=true;  	// Compute the activity coefficients of the organic phase for this compound?
  AR0088.compute_gamma_aq=true;  	// Compute the activity coefficients of the aqueous phase for this compound?
  AR0088.Koligo_org=0.0;      	//oligomeriation constant in the organic phase
  AR0088.rho=1300.;  
  AR0088.is_monomer=false;
  AR0088.rion=false;
  AR0088.KDiffusion_air=1.0e-5;
  //  AR0088.accomodation_coefficient=alpha;
  AR0088.viscosity=1.68e12;
  AR0088.is_solid=false;
  AR0088.is_generic=false;
  
  //Group: if no functionnal group in the species use the default species
  //for the computation of activity coefficients
  //
  double group_tmp_AR0088 [] = {0.,0.,0.,0., // group C
				  0.,0.,0.,0., //group C[OH]
				  0.,0.,0.,0., //group Calcohol
				  0.,0.,0.,0., //group Calcohol-tail
				  0.,0.,0.,0.,0., //group C=C
				  4.,0., //group aromatic carbon (AC)
				  1.,0.,0., // group //AC-C
				  0.,  //group OH
				  0., //group H2O
				  1., //group ACOH
				  0.,0., //group ketone
				  0.,   //group aldehyde  
				  0.,0., //group ester
				  0.,0.,0., //group ether 
				  0.,  //group acid
				  0.,   //group ACNO2
				  0.,0.,0., //group NO3
				  0.,0.,0., //group CO-OH
  			          0.,0.,0.,0.,0.,0.,0.,0.,0., //group CO-OC
			          0.,  //group PAN
			          0.,  //group CO-OOH
			          0.,  //group O=COC=O
			          0.,0.,0.}; //group CHxNO2
  
  size = sizeof(group_tmp_AR0088)/sizeof(double);
  assert(size = 60);
	
  for(int i = 0; i < size; ++i)
	AR0088.groups[i] = group_tmp_AR0088[i];

  // Search the species name in the aerosol species list 
  // and add the species if its name matches with
  // the given list.
  add_species_ssh(surrogate, AR0088, species_list_aer, molecular_weight_aer,
              accomodation_coefficient,diffusion_coef,saturation_vapor_pressure,enthalpy_vaporization, henry, t_ref, mass_density, species_part,nlayer,i_hydrophilic,
		  N_inert, N_inorganic);



/* ==== AR0090 ==== */ 

  species AR0090;
  AR0090.name="AR0090";
  AR0090.is_inorganic_precursor=false;
  AR0090.Psat_ref=1.173e-1; 		// Saturation vapor pressure at Tref (torr)
  AR0090.Tref=298;         		// Temperature of reference (K)
  AR0090.deltaH=62.08;     		// Enthalpy of vaporization (kJ/mol)
  AR0090.Henry=0.0;     		// Henry's law constant at Tref (M/atm)
  AR0090.aq_type="none"; 		// "none","diacid","monoacid" or "aldehyde"
  // AR0090.Kacidity1=XXacidityXX;   		// First acidity constant
  AR0090.hydrophilic=true;   	// Does the species condense on the aqueous phase?
  AR0090.hydrophobic=true;  		// Does the species condense on the organic phase?
  AR0090.nonvolatile=false; 		// Is the compound nonvolatile?
  AR0090.kp_from_experiment=false; 	// Use experimental partitioning constant at Tref?
  AR0090.is_organic=true;  		// Is the compound organic?
  AR0090.compute_gamma_org=true;  	// Compute the activity coefficients of the organic phase for this compound?
  AR0090.compute_gamma_aq=true;  	// Compute the activity coefficients of the aqueous phase for this compound?
  AR0090.Koligo_org=0.0;      	//oligomeriation constant in the organic phase
  AR0090.rho=1300.;  
  AR0090.is_monomer=false;
  AR0090.rion=false;
  AR0090.KDiffusion_air=1.0e-5;
  //  AR0090.accomodation_coefficient=alpha;
  AR0090.viscosity=1.68e12;
  AR0090.is_solid=false;
  AR0090.is_generic=false;
  
  //Group: if no functionnal group in the species use the default species
  //for the computation of activity coefficients
  //
  double group_tmp_AR0090 [] = {0.,0.,0.,0., // group C
				  0.,0.,0.,0., //group C[OH]
				  0.,0.,0.,0., //group Calcohol
				  0.,0.,0.,0., //group Calcohol-tail
				  0.,0.,0.,0.,0., //group C=C
				  5.,1., //group aromatic carbon (AC)
				  0.,0.,0., // group //AC-C
				  0.,  //group OH
				  0., //group H2O
				  0., //group ACOH
				  0.,0., //group ketone
				  0.,   //group aldehyde  
				  0.,0., //group ester
				  0.,0.,0., //group ether 
				  0.,  //group acid
				  0.,   //group ACNO2
				  1.,0.,0., //group NO3
				  0.,0.,0., //group CO-OH
  			          0.,0.,0.,0.,0.,0.,0.,0.,0., //group CO-OC
			          0.,  //group PAN
			          0.,  //group CO-OOH
			          0.,  //group O=COC=O
			          0.,0.,0.}; //group CHxNO2
  
  size = sizeof(group_tmp_AR0090)/sizeof(double);
  assert(size = 60);
	
  for(int i = 0; i < size; ++i)
	AR0090.groups[i] = group_tmp_AR0090[i];

  // Search the species name in the aerosol species list 
  // and add the species if its name matches with
  // the given list.
  add_species_ssh(surrogate, AR0090, species_list_aer, molecular_weight_aer,
              accomodation_coefficient,diffusion_coef,saturation_vapor_pressure,enthalpy_vaporization, henry, t_ref, mass_density, species_part,nlayer,i_hydrophilic,
		  N_inert, N_inorganic);



/* ==== AR0093 ==== */ 

  species AR0093;
  AR0093.name="AR0093";
  AR0093.is_inorganic_precursor=false;
  AR0093.Psat_ref=2.720e-3; 		// Saturation vapor pressure at Tref (torr)
  AR0093.Tref=298;         		// Temperature of reference (K)
  AR0093.deltaH=81.41;     		// Enthalpy of vaporization (kJ/mol)
  AR0093.Henry=0.0;     		// Henry's law constant at Tref (M/atm)
  AR0093.aq_type="none"; 		// "none","diacid","monoacid" or "aldehyde"
  // AR0093.Kacidity1=XXacidityXX;   		// First acidity constant
  AR0093.hydrophilic=true;   	// Does the species condense on the aqueous phase?
  AR0093.hydrophobic=true;  		// Does the species condense on the organic phase?
  AR0093.nonvolatile=false; 		// Is the compound nonvolatile?
  AR0093.kp_from_experiment=false; 	// Use experimental partitioning constant at Tref?
  AR0093.is_organic=true;  		// Is the compound organic?
  AR0093.compute_gamma_org=true;  	// Compute the activity coefficients of the organic phase for this compound?
  AR0093.compute_gamma_aq=true;  	// Compute the activity coefficients of the aqueous phase for this compound?
  AR0093.Koligo_org=0.0;      	//oligomeriation constant in the organic phase
  AR0093.rho=1300.;  
  AR0093.is_monomer=false;
  AR0093.rion=false;
  AR0093.KDiffusion_air=1.0e-5;
  //  AR0093.accomodation_coefficient=alpha;
  AR0093.viscosity=1.68e12;
  AR0093.is_solid=false;
  AR0093.is_generic=false;
  
  //Group: if no functionnal group in the species use the default species
  //for the computation of activity coefficients
  //
  double group_tmp_AR0093 [] = {1.,0.,0.,0., // group C
				  0.,0.,1.,0., //group C[OH]
				  0.,0.,0.,0., //group Calcohol
				  0.,0.,0.,0., //group Calcohol-tail
				  0.,1.,0.,0.,0., //group C=C
				  0.,0., //group aromatic carbon (AC)
				  0.,0.,0., // group //AC-C
				  1.,  //group OH
				  0., //group H2O
				  0., //group ACOH
				  0.,0., //group ketone
				  0.,   //group aldehyde  
				  0.,0., //group ester
				  0.,0.,0., //group ether 
				  0.,  //group acid
				  0.,   //group ACNO2
				  0.,1.,0., //group NO3
				  0.,0.,0., //group CO-OH
  			          0.,0.,0.,0.,0.,0.,0.,1.,0., //group CO-OC
			          0.,  //group PAN
			          0.,  //group CO-OOH
			          0.,  //group O=COC=O
			          0.,0.,0.}; //group CHxNO2
  
  size = sizeof(group_tmp_AR0093)/sizeof(double);
  assert(size = 60);
	
  for(int i = 0; i < size; ++i)
	AR0093.groups[i] = group_tmp_AR0093[i];

  // Search the species name in the aerosol species list 
  // and add the species if its name matches with
  // the given list.
  add_species_ssh(surrogate, AR0093, species_list_aer, molecular_weight_aer,
              accomodation_coefficient,diffusion_coef,saturation_vapor_pressure,enthalpy_vaporization, henry, t_ref, mass_density, species_part,nlayer,i_hydrophilic,
		  N_inert, N_inorganic);



/* ==== AR0094 ==== */ 

  species AR0094;
  AR0094.name="AR0094";
  AR0094.is_inorganic_precursor=false;
  AR0094.Psat_ref=4.529e-4; 		// Saturation vapor pressure at Tref (torr)
  AR0094.Tref=298;         		// Temperature of reference (K)
  AR0094.deltaH=89.49;     		// Enthalpy of vaporization (kJ/mol)
  AR0094.Henry=0.0;     		// Henry's law constant at Tref (M/atm)
  AR0094.aq_type="none"; 		// "none","diacid","monoacid" or "aldehyde"
  // AR0094.Kacidity1=XXacidityXX;   		// First acidity constant
  AR0094.hydrophilic=true;   	// Does the species condense on the aqueous phase?
  AR0094.hydrophobic=true;  		// Does the species condense on the organic phase?
  AR0094.nonvolatile=false; 		// Is the compound nonvolatile?
  AR0094.kp_from_experiment=false; 	// Use experimental partitioning constant at Tref?
  AR0094.is_organic=true;  		// Is the compound organic?
  AR0094.compute_gamma_org=true;  	// Compute the activity coefficients of the organic phase for this compound?
  AR0094.compute_gamma_aq=true;  	// Compute the activity coefficients of the aqueous phase for this compound?
  AR0094.Koligo_org=0.0;      	//oligomeriation constant in the organic phase
  AR0094.rho=1300.;  
  AR0094.is_monomer=false;
  AR0094.rion=false;
  AR0094.KDiffusion_air=1.0e-5;
  //  AR0094.accomodation_coefficient=alpha;
  AR0094.viscosity=1.68e12;
  AR0094.is_solid=false;
  AR0094.is_generic=false;
  
  //Group: if no functionnal group in the species use the default species
  //for the computation of activity coefficients
  //
  double group_tmp_AR0094 [] = {1.,0.,0.,0., // group C
				  0.,0.,1.,0., //group C[OH]
				  0.,0.,0.,0., //group Calcohol
				  0.,0.,0.,0., //group Calcohol-tail
				  0.,1.,0.,0.,0., //group C=C
				  0.,0., //group aromatic carbon (AC)
				  0.,0.,0., // group //AC-C
				  1.,  //group OH
				  0., //group H2O
				  0., //group ACOH
				  0.,0., //group ketone
				  0.,   //group aldehyde  
				  0.,0., //group ester
				  0.,0.,0., //group ether 
				  0.,  //group acid
				  0.,   //group ACNO2
				  0.,0.,0., //group NO3
				  0.,0.,1., //group CO-OH
  			          0.,0.,0.,0.,0.,0.,0.,1.,0., //group CO-OC
			          0.,  //group PAN
			          0.,  //group CO-OOH
			          0.,  //group O=COC=O
			          0.,0.,0.}; //group CHxNO2
  
  size = sizeof(group_tmp_AR0094)/sizeof(double);
  assert(size = 60);
	
  for(int i = 0; i < size; ++i)
	AR0094.groups[i] = group_tmp_AR0094[i];

  // Search the species name in the aerosol species list 
  // and add the species if its name matches with
  // the given list.
  add_species_ssh(surrogate, AR0094, species_list_aer, molecular_weight_aer,
              accomodation_coefficient,diffusion_coef,saturation_vapor_pressure,enthalpy_vaporization, henry, t_ref, mass_density, species_part,nlayer,i_hydrophilic,
		  N_inert, N_inorganic);



/* ==== AR0104 ==== */ 

  species AR0104;
  AR0104.name="AR0104";
  AR0104.is_inorganic_precursor=false;
  AR0104.Psat_ref=3.374e-2; 		// Saturation vapor pressure at Tref (torr)
  AR0104.Tref=298;         		// Temperature of reference (K)
  AR0104.deltaH=68.25;     		// Enthalpy of vaporization (kJ/mol)
  AR0104.Henry=0.0;     		// Henry's law constant at Tref (M/atm)
  AR0104.aq_type="none"; 		// "none","diacid","monoacid" or "aldehyde"
  // AR0104.Kacidity1=XXacidityXX;   		// First acidity constant
  AR0104.hydrophilic=true;   	// Does the species condense on the aqueous phase?
  AR0104.hydrophobic=true;  		// Does the species condense on the organic phase?
  AR0104.nonvolatile=false; 		// Is the compound nonvolatile?
  AR0104.kp_from_experiment=false; 	// Use experimental partitioning constant at Tref?
  AR0104.is_organic=true;  		// Is the compound organic?
  AR0104.compute_gamma_org=true;  	// Compute the activity coefficients of the organic phase for this compound?
  AR0104.compute_gamma_aq=true;  	// Compute the activity coefficients of the aqueous phase for this compound?
  AR0104.Koligo_org=0.0;      	//oligomeriation constant in the organic phase
  AR0104.rho=1300.;  
  AR0104.is_monomer=false;
  AR0104.rion=false;
  AR0104.KDiffusion_air=1.0e-5;
  //  AR0104.accomodation_coefficient=alpha;
  AR0104.viscosity=1.68e12;
  AR0104.is_solid=false;
  AR0104.is_generic=false;
  
  //Group: if no functionnal group in the species use the default species
  //for the computation of activity coefficients
  //
  double group_tmp_AR0104 [] = {0.,0.,0.,0., // group C
				  0.,0.,1.,0., //group C[OH]
				  0.,0.,0.,0., //group Calcohol
				  0.,0.,0.,0., //group Calcohol-tail
				  0.,1.,0.,0.,0., //group C=C
				  0.,0., //group aromatic carbon (AC)
				  0.,0.,0., // group //AC-C
				  1.,  //group OH
				  0., //group H2O
				  0., //group ACOH
				  1.,0., //group ketone
				  0.,   //group aldehyde  
				  0.,0., //group ester
				  0.,0.,0., //group ether 
				  0.,  //group acid
				  0.,   //group ACNO2
				  0.,0.,0., //group NO3
				  0.,0.,0., //group CO-OH
  			          0.,0.,0.,0.,0.,0.,0.,1.,0., //group CO-OC
			          0.,  //group PAN
			          0.,  //group CO-OOH
			          0.,  //group O=COC=O
			          0.,0.,0.}; //group CHxNO2
  
  size = sizeof(group_tmp_AR0104)/sizeof(double);
  assert(size = 60);
	
  for(int i = 0; i < size; ++i)
	AR0104.groups[i] = group_tmp_AR0104[i];

  // Search the species name in the aerosol species list 
  // and add the species if its name matches with
  // the given list.
  add_species_ssh(surrogate, AR0104, species_list_aer, molecular_weight_aer,
              accomodation_coefficient,diffusion_coef,saturation_vapor_pressure,enthalpy_vaporization, henry, t_ref, mass_density, species_part,nlayer,i_hydrophilic,
		  N_inert, N_inorganic);



/* ==== AR0105 ==== */ 

  species AR0105;
  AR0105.name="AR0105";
  AR0105.is_inorganic_precursor=false;
  AR0105.Psat_ref=8.381e-3; 		// Saturation vapor pressure at Tref (torr)
  AR0105.Tref=298;         		// Temperature of reference (K)
  AR0105.deltaH=74.46;     		// Enthalpy of vaporization (kJ/mol)
  AR0105.Henry=0.0;     		// Henry's law constant at Tref (M/atm)
  AR0105.aq_type="none"; 		// "none","diacid","monoacid" or "aldehyde"
  // AR0105.Kacidity1=XXacidityXX;   		// First acidity constant
  AR0105.hydrophilic=true;   	// Does the species condense on the aqueous phase?
  AR0105.hydrophobic=true;  		// Does the species condense on the organic phase?
  AR0105.nonvolatile=false; 		// Is the compound nonvolatile?
  AR0105.kp_from_experiment=false; 	// Use experimental partitioning constant at Tref?
  AR0105.is_organic=true;  		// Is the compound organic?
  AR0105.compute_gamma_org=true;  	// Compute the activity coefficients of the organic phase for this compound?
  AR0105.compute_gamma_aq=true;  	// Compute the activity coefficients of the aqueous phase for this compound?
  AR0105.Koligo_org=0.0;      	//oligomeriation constant in the organic phase
  AR0105.rho=1300.;  
  AR0105.is_monomer=false;
  AR0105.rion=false;
  AR0105.KDiffusion_air=1.0e-5;
  //  AR0105.accomodation_coefficient=alpha;
  AR0105.viscosity=1.68e12;
  AR0105.is_solid=false;
  AR0105.is_generic=false;
  
  //Group: if no functionnal group in the species use the default species
  //for the computation of activity coefficients
  //
  double group_tmp_AR0105 [] = {0.,0.,1.,0., // group C
				  0.,0.,0.,0., //group C[OH]
				  0.,0.,0.,0., //group Calcohol
				  0.,0.,0.,0., //group Calcohol-tail
				  0.,1.,0.,0.,0., //group C=C
				  0.,0., //group aromatic carbon (AC)
				  0.,0.,0., // group //AC-C
				  0.,  //group OH
				  0., //group H2O
				  0., //group ACOH
				  1.,0., //group ketone
				  0.,   //group aldehyde  
				  0.,0., //group ester
				  0.,0.,1., //group ether 
				  0.,  //group acid
				  0.,   //group ACNO2
				  0.,0.,0., //group NO3
				  0.,0.,0., //group CO-OH
  			          0.,0.,0.,0.,0.,0.,0.,0.,0., //group CO-OC
			          1.,  //group PAN
			          0.,  //group CO-OOH
			          0.,  //group O=COC=O
			          0.,0.,0.}; //group CHxNO2
  
  size = sizeof(group_tmp_AR0105)/sizeof(double);
  assert(size = 60);
	
  for(int i = 0; i < size; ++i)
	AR0105.groups[i] = group_tmp_AR0105[i];

  // Search the species name in the aerosol species list 
  // and add the species if its name matches with
  // the given list.
  add_species_ssh(surrogate, AR0105, species_list_aer, molecular_weight_aer,
              accomodation_coefficient,diffusion_coef,saturation_vapor_pressure,enthalpy_vaporization, henry, t_ref, mass_density, species_part,nlayer,i_hydrophilic,
		  N_inert, N_inorganic);



/* ==== AR0109 ==== */ 

  species AR0109;
  AR0109.name="AR0109";
  AR0109.is_inorganic_precursor=false;
  AR0109.Psat_ref=3.257e-5; 		// Saturation vapor pressure at Tref (torr)
  AR0109.Tref=298;         		// Temperature of reference (K)
  AR0109.deltaH=101.08;     		// Enthalpy of vaporization (kJ/mol)
  AR0109.Henry=0.0;     		// Henry's law constant at Tref (M/atm)
  AR0109.aq_type="none"; 		// "none","diacid","monoacid" or "aldehyde"
  // AR0109.Kacidity1=XXacidityXX;   		// First acidity constant
  AR0109.hydrophilic=true;   	// Does the species condense on the aqueous phase?
  AR0109.hydrophobic=true;  		// Does the species condense on the organic phase?
  AR0109.nonvolatile=false; 		// Is the compound nonvolatile?
  AR0109.kp_from_experiment=false; 	// Use experimental partitioning constant at Tref?
  AR0109.is_organic=true;  		// Is the compound organic?
  AR0109.compute_gamma_org=true;  	// Compute the activity coefficients of the organic phase for this compound?
  AR0109.compute_gamma_aq=true;  	// Compute the activity coefficients of the aqueous phase for this compound?
  AR0109.Koligo_org=0.0;      	//oligomeriation constant in the organic phase
  AR0109.rho=1300.;  
  AR0109.is_monomer=false;
  AR0109.rion=false;
  AR0109.KDiffusion_air=1.0e-5;
  //  AR0109.accomodation_coefficient=alpha;
  AR0109.viscosity=1.68e12;
  AR0109.is_solid=false;
  AR0109.is_generic=false;
  
  //Group: if no functionnal group in the species use the default species
  //for the computation of activity coefficients
  //
  double group_tmp_AR0109 [] = {0.,0.,1.,0., // group C
				  0.,0.,1.,0., //group C[OH]
				  0.,0.,0.,0., //group Calcohol
				  0.,0.,0.,0., //group Calcohol-tail
				  0.,0.,0.,0.,0., //group C=C
				  0.,0., //group aromatic carbon (AC)
				  0.,0.,0., // group //AC-C
				  1.,  //group OH
				  0., //group H2O
				  0., //group ACOH
				  1.,0., //group ketone
				  1.,   //group aldehyde  
				  0.,0., //group ester
				  0.,0.,1., //group ether 
				  0.,  //group acid
				  0.,   //group ACNO2
				  0.,1.,0., //group NO3
				  0.,0.,0., //group CO-OH
  			          0.,0.,0.,0.,0.,0.,0.,0.,0., //group CO-OC
			          0.,  //group PAN
			          0.,  //group CO-OOH
			          0.,  //group O=COC=O
			          0.,0.,0.}; //group CHxNO2
  
  size = sizeof(group_tmp_AR0109)/sizeof(double);
  assert(size = 60);
	
  for(int i = 0; i < size; ++i)
	AR0109.groups[i] = group_tmp_AR0109[i];

  // Search the species name in the aerosol species list 
  // and add the species if its name matches with
  // the given list.
  add_species_ssh(surrogate, AR0109, species_list_aer, molecular_weight_aer,
              accomodation_coefficient,diffusion_coef,saturation_vapor_pressure,enthalpy_vaporization, henry, t_ref, mass_density, species_part,nlayer,i_hydrophilic,
		  N_inert, N_inorganic);



/* ==== AR0110 ==== */ 

  species AR0110;
  AR0110.name="AR0110";
  AR0110.is_inorganic_precursor=false;
  AR0110.Psat_ref=4.446e-6; 		// Saturation vapor pressure at Tref (torr)
  AR0110.Tref=298;         		// Temperature of reference (K)
  AR0110.deltaH=110.53;     		// Enthalpy of vaporization (kJ/mol)
  AR0110.Henry=0.0;     		// Henry's law constant at Tref (M/atm)
  AR0110.aq_type="none"; 		// "none","diacid","monoacid" or "aldehyde"
  // AR0110.Kacidity1=XXacidityXX;   		// First acidity constant
  AR0110.hydrophilic=true;   	// Does the species condense on the aqueous phase?
  AR0110.hydrophobic=true;  		// Does the species condense on the organic phase?
  AR0110.nonvolatile=false; 		// Is the compound nonvolatile?
  AR0110.kp_from_experiment=false; 	// Use experimental partitioning constant at Tref?
  AR0110.is_organic=true;  		// Is the compound organic?
  AR0110.compute_gamma_org=true;  	// Compute the activity coefficients of the organic phase for this compound?
  AR0110.compute_gamma_aq=true;  	// Compute the activity coefficients of the aqueous phase for this compound?
  AR0110.Koligo_org=0.0;      	//oligomeriation constant in the organic phase
  AR0110.rho=1300.;  
  AR0110.is_monomer=false;
  AR0110.rion=false;
  AR0110.KDiffusion_air=1.0e-5;
  //  AR0110.accomodation_coefficient=alpha;
  AR0110.viscosity=1.68e12;
  AR0110.is_solid=false;
  AR0110.is_generic=false;
  
  //Group: if no functionnal group in the species use the default species
  //for the computation of activity coefficients
  //
  double group_tmp_AR0110 [] = {0.,0.,1.,0., // group C
				  0.,0.,1.,0., //group C[OH]
				  0.,0.,0.,0., //group Calcohol
				  0.,0.,0.,0., //group Calcohol-tail
				  0.,0.,0.,0.,0., //group C=C
				  0.,0., //group aromatic carbon (AC)
				  0.,0.,0., // group //AC-C
				  1.,  //group OH
				  0., //group H2O
				  0., //group ACOH
				  1.,0., //group ketone
				  1.,   //group aldehyde  
				  0.,0., //group ester
				  0.,0.,1., //group ether 
				  0.,  //group acid
				  0.,   //group ACNO2
				  0.,0.,0., //group NO3
				  0.,1.,0., //group CO-OH
  			          0.,0.,0.,0.,0.,0.,0.,0.,0., //group CO-OC
			          0.,  //group PAN
			          0.,  //group CO-OOH
			          0.,  //group O=COC=O
			          0.,0.,0.}; //group CHxNO2
  
  size = sizeof(group_tmp_AR0110)/sizeof(double);
  assert(size = 60);
	
  for(int i = 0; i < size; ++i)
	AR0110.groups[i] = group_tmp_AR0110[i];

  // Search the species name in the aerosol species list 
  // and add the species if its name matches with
  // the given list.
  add_species_ssh(surrogate, AR0110, species_list_aer, molecular_weight_aer,
              accomodation_coefficient,diffusion_coef,saturation_vapor_pressure,enthalpy_vaporization, henry, t_ref, mass_density, species_part,nlayer,i_hydrophilic,
		  N_inert, N_inorganic);



/* ==== AR0113 ==== */ 

  species AR0113;
  AR0113.name="AR0113";
  AR0113.is_inorganic_precursor=false;
  AR0113.Psat_ref=2.685e-3; 		// Saturation vapor pressure at Tref (torr)
  AR0113.Tref=298;         		// Temperature of reference (K)
  AR0113.deltaH=76.93;     		// Enthalpy of vaporization (kJ/mol)
  AR0113.Henry=0.0;     		// Henry's law constant at Tref (M/atm)
  AR0113.aq_type="none"; 		// "none","diacid","monoacid" or "aldehyde"
  // AR0113.Kacidity1=XXacidityXX;   		// First acidity constant
  AR0113.hydrophilic=true;   	// Does the species condense on the aqueous phase?
  AR0113.hydrophobic=true;  		// Does the species condense on the organic phase?
  AR0113.nonvolatile=false; 		// Is the compound nonvolatile?
  AR0113.kp_from_experiment=false; 	// Use experimental partitioning constant at Tref?
  AR0113.is_organic=true;  		// Is the compound organic?
  AR0113.compute_gamma_org=true;  	// Compute the activity coefficients of the organic phase for this compound?
  AR0113.compute_gamma_aq=true;  	// Compute the activity coefficients of the aqueous phase for this compound?
  AR0113.Koligo_org=0.0;      	//oligomeriation constant in the organic phase
  AR0113.rho=1300.;  
  AR0113.is_monomer=false;
  AR0113.rion=false;
  AR0113.KDiffusion_air=1.0e-5;
  //  AR0113.accomodation_coefficient=alpha;
  AR0113.viscosity=1.68e12;
  AR0113.is_solid=false;
  AR0113.is_generic=false;
  
  //Group: if no functionnal group in the species use the default species
  //for the computation of activity coefficients
  //
  double group_tmp_AR0113 [] = {0.,0.,0.,0., // group C
				  0.,0.,0.,0., //group C[OH]
				  0.,0.,0.,0., //group Calcohol
				  0.,0.,0.,0., //group Calcohol-tail
				  0.,0.,0.,0.,0., //group C=C
				  3.,0., //group aromatic carbon (AC)
				  1.,0.,0., // group //AC-C
				  0.,  //group OH
				  0., //group H2O
				  1., //group ACOH
				  0.,0., //group ketone
				  0.,   //group aldehyde  
				  0.,0., //group ester
				  0.,0.,0., //group ether 
				  0.,  //group acid
				  1.,   //group ACNO2
				  0.,0.,0., //group NO3
				  0.,0.,0., //group CO-OH
  			          0.,0.,0.,0.,0.,0.,0.,0.,0., //group CO-OC
			          0.,  //group PAN
			          0.,  //group CO-OOH
			          0.,  //group O=COC=O
			          0.,0.,0.}; //group CHxNO2
  
  size = sizeof(group_tmp_AR0113)/sizeof(double);
  assert(size = 60);
	
  for(int i = 0; i < size; ++i)
	AR0113.groups[i] = group_tmp_AR0113[i];

  // Search the species name in the aerosol species list 
  // and add the species if its name matches with
  // the given list.
  add_species_ssh(surrogate, AR0113, species_list_aer, molecular_weight_aer,
              accomodation_coefficient,diffusion_coef,saturation_vapor_pressure,enthalpy_vaporization, henry, t_ref, mass_density, species_part,nlayer,i_hydrophilic,
		  N_inert, N_inorganic);



/* ==== AR0115 ==== */ 

  species AR0115;
  AR0115.name="AR0115";
  AR0115.is_inorganic_precursor=false;
  AR0115.Psat_ref=8.105e-6; 		// Saturation vapor pressure at Tref (torr)
  AR0115.Tref=298;         		// Temperature of reference (K)
  AR0115.deltaH=109.69;     		// Enthalpy of vaporization (kJ/mol)
  AR0115.Henry=0.0;     		// Henry's law constant at Tref (M/atm)
  AR0115.aq_type="none"; 		// "none","diacid","monoacid" or "aldehyde"
  // AR0115.Kacidity1=XXacidityXX;   		// First acidity constant
  AR0115.hydrophilic=true;   	// Does the species condense on the aqueous phase?
  AR0115.hydrophobic=true;  		// Does the species condense on the organic phase?
  AR0115.nonvolatile=false; 		// Is the compound nonvolatile?
  AR0115.kp_from_experiment=false; 	// Use experimental partitioning constant at Tref?
  AR0115.is_organic=true;  		// Is the compound organic?
  AR0115.compute_gamma_org=true;  	// Compute the activity coefficients of the organic phase for this compound?
  AR0115.compute_gamma_aq=true;  	// Compute the activity coefficients of the aqueous phase for this compound?
  AR0115.Koligo_org=0.0;      	//oligomeriation constant in the organic phase
  AR0115.rho=1300.;  
  AR0115.is_monomer=false;
  AR0115.rion=false;
  AR0115.KDiffusion_air=1.0e-5;
  //  AR0115.accomodation_coefficient=alpha;
  AR0115.viscosity=1.68e12;
  AR0115.is_solid=false;
  AR0115.is_generic=false;
  
  //Group: if no functionnal group in the species use the default species
  //for the computation of activity coefficients
  //
  double group_tmp_AR0115 [] = {1.,0.,0.,0., // group C
				  0.,0.,1.,0., //group C[OH]
				  0.,0.,0.,0., //group Calcohol
				  0.,0.,0.,0., //group Calcohol-tail
				  0.,1.,0.,0.,0., //group C=C
				  0.,0., //group aromatic carbon (AC)
				  0.,0.,0., // group //AC-C
				  2.,  //group OH
				  0., //group H2O
				  0., //group ACOH
				  0.,0., //group ketone
				  0.,   //group aldehyde  
				  0.,0., //group ester
				  0.,0.,0., //group ether 
				  0.,  //group acid
				  0.,   //group ACNO2
				  0.,0.,0., //group NO3
				  0.,0.,1., //group CO-OH
  			          0.,0.,0.,0.,0.,0.,0.,1.,0., //group CO-OC
			          0.,  //group PAN
			          0.,  //group CO-OOH
			          0.,  //group O=COC=O
			          0.,0.,0.}; //group CHxNO2
  
  size = sizeof(group_tmp_AR0115)/sizeof(double);
  assert(size = 60);
	
  for(int i = 0; i < size; ++i)
	AR0115.groups[i] = group_tmp_AR0115[i];

  // Search the species name in the aerosol species list 
  // and add the species if its name matches with
  // the given list.
  add_species_ssh(surrogate, AR0115, species_list_aer, molecular_weight_aer,
              accomodation_coefficient,diffusion_coef,saturation_vapor_pressure,enthalpy_vaporization, henry, t_ref, mass_density, species_part,nlayer,i_hydrophilic,
		  N_inert, N_inorganic);



/* ==== AR0124 ==== */ 

  species AR0124;
  AR0124.name="AR0124";
  AR0124.is_inorganic_precursor=false;
  AR0124.Psat_ref=5.365e-4; 		// Saturation vapor pressure at Tref (torr)
  AR0124.Tref=298;         		// Temperature of reference (K)
  AR0124.deltaH=87.74;     		// Enthalpy of vaporization (kJ/mol)
  AR0124.Henry=0.0;     		// Henry's law constant at Tref (M/atm)
  AR0124.aq_type="none"; 		// "none","diacid","monoacid" or "aldehyde"
  // AR0124.Kacidity1=XXacidityXX;   		// First acidity constant
  AR0124.hydrophilic=true;   	// Does the species condense on the aqueous phase?
  AR0124.hydrophobic=true;  		// Does the species condense on the organic phase?
  AR0124.nonvolatile=false; 		// Is the compound nonvolatile?
  AR0124.kp_from_experiment=false; 	// Use experimental partitioning constant at Tref?
  AR0124.is_organic=true;  		// Is the compound organic?
  AR0124.compute_gamma_org=true;  	// Compute the activity coefficients of the organic phase for this compound?
  AR0124.compute_gamma_aq=true;  	// Compute the activity coefficients of the aqueous phase for this compound?
  AR0124.Koligo_org=0.0;      	//oligomeriation constant in the organic phase
  AR0124.rho=1300.;  
  AR0124.is_monomer=false;
  AR0124.rion=false;
  AR0124.KDiffusion_air=1.0e-5;
  //  AR0124.accomodation_coefficient=alpha;
  AR0124.viscosity=1.68e12;
  AR0124.is_solid=false;
  AR0124.is_generic=false;
  
  //Group: if no functionnal group in the species use the default species
  //for the computation of activity coefficients
  //
  double group_tmp_AR0124 [] = {0.,0.,0.,0., // group C
				  0.,0.,0.,1., //group C[OH]
				  0.,0.,0.,0., //group Calcohol
				  0.,0.,0.,0., //group Calcohol-tail
				  0.,0.,0.,0.,0., //group C=C
				  0.,0., //group aromatic carbon (AC)
				  0.,0.,0., // group //AC-C
				  1.,  //group OH
				  0., //group H2O
				  0., //group ACOH
				  1.,1., //group ketone
				  1.,   //group aldehyde  
				  0.,0., //group ester
				  0.,0.,1., //group ether 
				  0.,  //group acid
				  0.,   //group ACNO2
				  0.,0.,0., //group NO3
				  0.,0.,0., //group CO-OH
  			          0.,0.,0.,0.,0.,0.,0.,0.,0., //group CO-OC
			          0.,  //group PAN
			          0.,  //group CO-OOH
			          0.,  //group O=COC=O
			          0.,0.,0.}; //group CHxNO2
  
  size = sizeof(group_tmp_AR0124)/sizeof(double);
  assert(size = 60);
	
  for(int i = 0; i < size; ++i)
	AR0124.groups[i] = group_tmp_AR0124[i];

  // Search the species name in the aerosol species list 
  // and add the species if its name matches with
  // the given list.
  add_species_ssh(surrogate, AR0124, species_list_aer, molecular_weight_aer,
              accomodation_coefficient,diffusion_coef,saturation_vapor_pressure,enthalpy_vaporization, henry, t_ref, mass_density, species_part,nlayer,i_hydrophilic,
		  N_inert, N_inorganic);



/* ==== AR0127 ==== */ 

  species AR0127;
  AR0127.name="AR0127";
  AR0127.is_inorganic_precursor=false;
  AR0127.Psat_ref=5.111e-1; 		// Saturation vapor pressure at Tref (torr)
  AR0127.Tref=298;         		// Temperature of reference (K)
  AR0127.deltaH=53.60;     		// Enthalpy of vaporization (kJ/mol)
  AR0127.Henry=0.0;     		// Henry's law constant at Tref (M/atm)
  AR0127.aq_type="none"; 		// "none","diacid","monoacid" or "aldehyde"
  // AR0127.Kacidity1=XXacidityXX;   		// First acidity constant
  AR0127.hydrophilic=true;   	// Does the species condense on the aqueous phase?
  AR0127.hydrophobic=true;  		// Does the species condense on the organic phase?
  AR0127.nonvolatile=false; 		// Is the compound nonvolatile?
  AR0127.kp_from_experiment=false; 	// Use experimental partitioning constant at Tref?
  AR0127.is_organic=true;  		// Is the compound organic?
  AR0127.compute_gamma_org=true;  	// Compute the activity coefficients of the organic phase for this compound?
  AR0127.compute_gamma_aq=true;  	// Compute the activity coefficients of the aqueous phase for this compound?
  AR0127.Koligo_org=0.0;      	//oligomeriation constant in the organic phase
  AR0127.rho=1300.;  
  AR0127.is_monomer=false;
  AR0127.rion=false;
  AR0127.KDiffusion_air=1.0e-5;
  //  AR0127.accomodation_coefficient=alpha;
  AR0127.viscosity=1.68e12;
  AR0127.is_solid=false;
  AR0127.is_generic=false;
  
  //Group: if no functionnal group in the species use the default species
  //for the computation of activity coefficients
  //
  double group_tmp_AR0127 [] = {1.,0.,0.,0., // group C
				  0.,0.,0.,0., //group C[OH]
				  0.,0.,0.,0., //group Calcohol
				  0.,0.,0.,0., //group Calcohol-tail
				  0.,0.,0.,2.,0., //group C=C
				  0.,0., //group aromatic carbon (AC)
				  0.,0.,0., // group //AC-C
				  0.,  //group OH
				  0., //group H2O
				  0., //group ACOH
				  0.,2., //group ketone
				  0.,   //group aldehyde  
				  0.,0., //group ester
				  0.,0.,0., //group ether 
				  0.,  //group acid
				  0.,   //group ACNO2
				  0.,0.,0., //group NO3
				  0.,0.,0., //group CO-OH
  			          0.,0.,0.,0.,0.,0.,0.,0.,0., //group CO-OC
			          0.,  //group PAN
			          0.,  //group CO-OOH
			          0.,  //group O=COC=O
			          0.,0.,0.}; //group CHxNO2
  
  size = sizeof(group_tmp_AR0127)/sizeof(double);
  assert(size = 60);
	
  for(int i = 0; i < size; ++i)
	AR0127.groups[i] = group_tmp_AR0127[i];

  // Search the species name in the aerosol species list 
  // and add the species if its name matches with
  // the given list.
  add_species_ssh(surrogate, AR0127, species_list_aer, molecular_weight_aer,
              accomodation_coefficient,diffusion_coef,saturation_vapor_pressure,enthalpy_vaporization, henry, t_ref, mass_density, species_part,nlayer,i_hydrophilic,
		  N_inert, N_inorganic);



/* ==== AR0128 ==== */ 

  species AR0128;
  AR0128.name="AR0128";
  AR0128.is_inorganic_precursor=false;
  AR0128.Psat_ref=3.083e-5; 		// Saturation vapor pressure at Tref (torr)
  AR0128.Tref=298;         		// Temperature of reference (K)
  AR0128.deltaH=97.01;     		// Enthalpy of vaporization (kJ/mol)
  AR0128.Henry=0.0;     		// Henry's law constant at Tref (M/atm)
  AR0128.aq_type="none"; 		// "none","diacid","monoacid" or "aldehyde"
  // AR0128.Kacidity1=XXacidityXX;   		// First acidity constant
  AR0128.hydrophilic=true;   	// Does the species condense on the aqueous phase?
  AR0128.hydrophobic=true;  		// Does the species condense on the organic phase?
  AR0128.nonvolatile=false; 		// Is the compound nonvolatile?
  AR0128.kp_from_experiment=false; 	// Use experimental partitioning constant at Tref?
  AR0128.is_organic=true;  		// Is the compound organic?
  AR0128.compute_gamma_org=true;  	// Compute the activity coefficients of the organic phase for this compound?
  AR0128.compute_gamma_aq=true;  	// Compute the activity coefficients of the aqueous phase for this compound?
  AR0128.Koligo_org=0.0;      	//oligomeriation constant in the organic phase
  AR0128.rho=1300.;  
  AR0128.is_monomer=false;
  AR0128.rion=false;
  AR0128.KDiffusion_air=1.0e-5;
  //  AR0128.accomodation_coefficient=alpha;
  AR0128.viscosity=1.68e12;
  AR0128.is_solid=false;
  AR0128.is_generic=false;
  
  //Group: if no functionnal group in the species use the default species
  //for the computation of activity coefficients
  //
  double group_tmp_AR0128 [] = {0.,0.,0.,0., // group C
				  0.,0.,0.,0., //group C[OH]
				  0.,0.,0.,0., //group Calcohol
				  0.,0.,0.,0., //group Calcohol-tail
				  0.,0.,0.,0.,0., //group C=C
				  2.,0., //group aromatic carbon (AC)
				  1.,0.,0., // group //AC-C
				  0.,  //group OH
				  0., //group H2O
				  2., //group ACOH
				  0.,0., //group ketone
				  0.,   //group aldehyde  
				  0.,0., //group ester
				  0.,0.,0., //group ether 
				  0.,  //group acid
				  1.,   //group ACNO2
				  0.,0.,0., //group NO3
				  0.,0.,0., //group CO-OH
  			          0.,0.,0.,0.,0.,0.,0.,0.,0., //group CO-OC
			          0.,  //group PAN
			          0.,  //group CO-OOH
			          0.,  //group O=COC=O
			          0.,0.,0.}; //group CHxNO2
  
  size = sizeof(group_tmp_AR0128)/sizeof(double);
  assert(size = 60);
	
  for(int i = 0; i < size; ++i)
	AR0128.groups[i] = group_tmp_AR0128[i];

  // Search the species name in the aerosol species list 
  // and add the species if its name matches with
  // the given list.
  add_species_ssh(surrogate, AR0128, species_list_aer, molecular_weight_aer,
              accomodation_coefficient,diffusion_coef,saturation_vapor_pressure,enthalpy_vaporization, henry, t_ref, mass_density, species_part,nlayer,i_hydrophilic,
		  N_inert, N_inorganic);



/* ==== AR0130 ==== */ 

  species AR0130;
  AR0130.name="AR0130";
  AR0130.is_inorganic_precursor=false;
  AR0130.Psat_ref=6.917e-2; 		// Saturation vapor pressure at Tref (torr)
  AR0130.Tref=298;         		// Temperature of reference (K)
  AR0130.deltaH=62.99;     		// Enthalpy of vaporization (kJ/mol)
  AR0130.Henry=0.0;     		// Henry's law constant at Tref (M/atm)
  AR0130.aq_type="none"; 		// "none","diacid","monoacid" or "aldehyde"
  // AR0130.Kacidity1=XXacidityXX;   		// First acidity constant
  AR0130.hydrophilic=true;   	// Does the species condense on the aqueous phase?
  AR0130.hydrophobic=true;  		// Does the species condense on the organic phase?
  AR0130.nonvolatile=false; 		// Is the compound nonvolatile?
  AR0130.kp_from_experiment=false; 	// Use experimental partitioning constant at Tref?
  AR0130.is_organic=true;  		// Is the compound organic?
  AR0130.compute_gamma_org=true;  	// Compute the activity coefficients of the organic phase for this compound?
  AR0130.compute_gamma_aq=true;  	// Compute the activity coefficients of the aqueous phase for this compound?
  AR0130.Koligo_org=0.0;      	//oligomeriation constant in the organic phase
  AR0130.rho=1300.;  
  AR0130.is_monomer=false;
  AR0130.rion=false;
  AR0130.KDiffusion_air=1.0e-5;
  //  AR0130.accomodation_coefficient=alpha;
  AR0130.viscosity=1.68e12;
  AR0130.is_solid=false;
  AR0130.is_generic=false;
  
  //Group: if no functionnal group in the species use the default species
  //for the computation of activity coefficients
  //
  double group_tmp_AR0130 [] = {0.,0.,0.,0., // group C
				  0.,0.,0.,0., //group C[OH]
				  0.,0.,0.,0., //group Calcohol
				  0.,0.,0.,0., //group Calcohol-tail
				  0.,0.,0.,0.,0., //group C=C
				  5.,1., //group aromatic carbon (AC)
				  0.,0.,0., // group //AC-C
				  0.,  //group OH
				  0., //group H2O
				  0., //group ACOH
				  0.,0., //group ketone
				  0.,   //group aldehyde  
				  0.,0., //group ester
				  0.,0.,0., //group ether 
				  0.,  //group acid
				  0.,   //group ACNO2
				  0.,0.,0., //group NO3
				  0.,0.,0., //group CO-OH
  			          0.,0.,0.,0.,0.,0.,0.,0.,0., //group CO-OC
			          1.,  //group PAN
			          0.,  //group CO-OOH
			          0.,  //group O=COC=O
			          0.,0.,0.}; //group CHxNO2
  
  size = sizeof(group_tmp_AR0130)/sizeof(double);
  assert(size = 60);
	
  for(int i = 0; i < size; ++i)
	AR0130.groups[i] = group_tmp_AR0130[i];

  // Search the species name in the aerosol species list 
  // and add the species if its name matches with
  // the given list.
  add_species_ssh(surrogate, AR0130, species_list_aer, molecular_weight_aer,
              accomodation_coefficient,diffusion_coef,saturation_vapor_pressure,enthalpy_vaporization, henry, t_ref, mass_density, species_part,nlayer,i_hydrophilic,
		  N_inert, N_inorganic);



/* ==== AR0132 ==== */ 

  species AR0132;
  AR0132.name="AR0132";
  AR0132.is_inorganic_precursor=false;
  AR0132.Psat_ref=8.920e-3; 		// Saturation vapor pressure at Tref (torr)
  AR0132.Tref=298;         		// Temperature of reference (K)
  AR0132.deltaH=74.079;     		// Enthalpy of vaporization (kJ/mol)
  AR0132.Henry=0.0;     		// Henry's law constant at Tref (M/atm)
  AR0132.aq_type="none"; 		// "none","diacid","monoacid" or "aldehyde"
  // AR0132.Kacidity1=XXacidityXX;   		// First acidity constant
  AR0132.hydrophilic=true;   	// Does the species condense on the aqueous phase?
  AR0132.hydrophobic=true;  		// Does the species condense on the organic phase?
  AR0132.nonvolatile=false; 		// Is the compound nonvolatile?
  AR0132.kp_from_experiment=false; 	// Use experimental partitioning constant at Tref?
  AR0132.is_organic=true;  		// Is the compound organic?
  AR0132.compute_gamma_org=true;  	// Compute the activity coefficients of the organic phase for this compound?
  AR0132.compute_gamma_aq=true;  	// Compute the activity coefficients of the aqueous phase for this compound?
  AR0132.Koligo_org=0.0;      	//oligomeriation constant in the organic phase
  AR0132.rho=1300.;  
  AR0132.is_monomer=false;
  AR0132.rion=false;
  AR0132.KDiffusion_air=1.0e-5;
  //  AR0132.accomodation_coefficient=alpha;
  AR0132.viscosity=1.68e12;
  AR0132.is_solid=false;
  AR0132.is_generic=false;
  
  //Group: if no functionnal group in the species use the default species
  //for the computation of activity coefficients
  //
  double group_tmp_AR0132 [] = {0.,0.,0.,0., // group C
				  0.,0.,0.,0., //group C[OH]
				  0.,0.,0.,0., //group Calcohol
				  0.,0.,0.,0., //group Calcohol-tail
				  0.,0.,0.,0.,0., //group C=C
				  5.,1., //group aromatic carbon (AC)
				  0.,0.,0., // group //AC-C
				  0.,  //group OH
				  0., //group H2O
				  0., //group ACOH
				  0.,0., //group ketone
				  0.,   //group aldehyde  
				  0.,0., //group ester
				  0.,0.,0., //group ether 
				  1.,  //group acid
				  0.,   //group ACNO2
				  0.,0.,0., //group NO3
				  0.,0.,0., //group CO-OH
  			          0.,0.,0.,0.,0.,0.,0.,0.,0., //group CO-OC
			          0.,  //group PAN
			          0.,  //group CO-OOH
			          0.,  //group O=COC=O
			          0.,0.,0.}; //group CHxNO2
  
  size = sizeof(group_tmp_AR0132)/sizeof(double);
  assert(size = 60);
	
  for(int i = 0; i < size; ++i)
	AR0132.groups[i] = group_tmp_AR0132[i];

  // Search the species name in the aerosol species list 
  // and add the species if its name matches with
  // the given list.
  add_species_ssh(surrogate, AR0132, species_list_aer, molecular_weight_aer,
              accomodation_coefficient,diffusion_coef,saturation_vapor_pressure,enthalpy_vaporization, henry, t_ref, mass_density, species_part,nlayer,i_hydrophilic,
		  N_inert, N_inorganic);



/* ==== AR0134 ==== */ 

  species AR0134;
  AR0134.name="AR0134";
  AR0134.is_inorganic_precursor=false;
  AR0134.Psat_ref=1.392e-5; 		// Saturation vapor pressure at Tref (torr)
  AR0134.Tref=298;         		// Temperature of reference (K)
  AR0134.deltaH=102.92;     		// Enthalpy of vaporization (kJ/mol)
  AR0134.Henry=0.0;     		// Henry's law constant at Tref (M/atm)
  AR0134.aq_type="none"; 		// "none","diacid","monoacid" or "aldehyde"
  // AR0134.Kacidity1=XXacidityXX;   		// First acidity constant
  AR0134.hydrophilic=true;   	// Does the species condense on the aqueous phase?
  AR0134.hydrophobic=true;  		// Does the species condense on the organic phase?
  AR0134.nonvolatile=false; 		// Is the compound nonvolatile?
  AR0134.kp_from_experiment=false; 	// Use experimental partitioning constant at Tref?
  AR0134.is_organic=true;  		// Is the compound organic?
  AR0134.compute_gamma_org=true;  	// Compute the activity coefficients of the organic phase for this compound?
  AR0134.compute_gamma_aq=true;  	// Compute the activity coefficients of the aqueous phase for this compound?
  AR0134.Koligo_org=0.0;      	//oligomeriation constant in the organic phase
  AR0134.rho=1300.;  
  AR0134.is_monomer=false;
  AR0134.rion=false;
  AR0134.KDiffusion_air=1.0e-5;
  //  AR0134.accomodation_coefficient=alpha;
  AR0134.viscosity=1.68e12;
  AR0134.is_solid=false;
  AR0134.is_generic=false;
  
  //Group: if no functionnal group in the species use the default species
  //for the computation of activity coefficients
  //
  double group_tmp_AR0134 [] = {0.,0.,1.,0., // group C
				  0.,0.,1.,0., //group C[OH]
				  0.,0.,0.,0., //group Calcohol
				  0.,0.,0.,0., //group Calcohol-tail
				  0.,0.,0.,0.,0., //group C=C
				  0.,0., //group aromatic carbon (AC)
				  0.,0.,0., // group //AC-C
				  1.,  //group OH
				  0., //group H2O
				  0., //group ACOH
				  0.,0., //group ketone
				  0.,   //group aldehyde  
				  0.,1., //group ester
				  0.,0.,0., //group ether 
				  0.,  //group acid
				  0.,   //group ACNO2
				  0.,0.,0., //group NO3
				  1.,0.,0., //group CO-OH
  			          0.,0.,0.,0.,0.,0.,0.,0.,0., //group CO-OC
			          0.,  //group PAN
			          0.,  //group CO-OOH
			          0.,  //group O=COC=O
			          0.,0.,0.}; //group CHxNO2
  
  size = sizeof(group_tmp_AR0134)/sizeof(double);
  assert(size = 60);
	
  for(int i = 0; i < size; ++i)
	AR0134.groups[i] = group_tmp_AR0134[i];

  // Search the species name in the aerosol species list 
  // and add the species if its name matches with
  // the given list.
  add_species_ssh(surrogate, AR0134, species_list_aer, molecular_weight_aer,
              accomodation_coefficient,diffusion_coef,saturation_vapor_pressure,enthalpy_vaporization, henry, t_ref, mass_density, species_part,nlayer,i_hydrophilic,
		  N_inert, N_inorganic);



/* ==== AR0138 ==== */ 

  species AR0138;
  AR0138.name="AR0138";
  AR0138.is_inorganic_precursor=false;
  AR0138.Psat_ref=1.516e-5; 		// Saturation vapor pressure at Tref (torr)
  AR0138.Tref=298;         		// Temperature of reference (K)
  AR0138.deltaH=100.0;     		// Enthalpy of vaporization (kJ/mol)
  AR0138.Henry=0.0;     		// Henry's law constant at Tref (M/atm)
  AR0138.aq_type="none"; 		// "none","diacid","monoacid" or "aldehyde"
  // AR0138.Kacidity1=XXacidityXX;   		// First acidity constant
  AR0138.hydrophilic=true;   	// Does the species condense on the aqueous phase?
  AR0138.hydrophobic=true;  		// Does the species condense on the organic phase?
  AR0138.nonvolatile=false; 		// Is the compound nonvolatile?
  AR0138.kp_from_experiment=false; 	// Use experimental partitioning constant at Tref?
  AR0138.is_organic=true;  		// Is the compound organic?
  AR0138.compute_gamma_org=true;  	// Compute the activity coefficients of the organic phase for this compound?
  AR0138.compute_gamma_aq=true;  	// Compute the activity coefficients of the aqueous phase for this compound?
  AR0138.Koligo_org=0.0;      	//oligomeriation constant in the organic phase
  AR0138.rho=1300.;  
  AR0138.is_monomer=false;
  AR0138.rion=false;
  AR0138.KDiffusion_air=1.0e-5;
  //  AR0138.accomodation_coefficient=alpha;
  AR0138.viscosity=1.68e12;
  AR0138.is_solid=false;
  AR0138.is_generic=false;
  
  //Group: if no functionnal group in the species use the default species
  //for the computation of activity coefficients
  //
  double group_tmp_AR0138 [] = {0.,0.,0.,0., // group C
				  0.,0.,0.,0., //group C[OH]
				  0.,0.,0.,0., //group Calcohol
				  0.,0.,0.,0., //group Calcohol-tail
				  0.,0.,0.,0.,0., //group C=C
				  2.,0., //group aromatic carbon (AC)
				  1.,0.,0., // group //AC-C
				  0.,  //group OH
				  0., //group H2O
				  1., //group ACOH
				  0.,0., //group ketone
				  0.,   //group aldehyde  
				  0.,0., //group ester
				  0.,0.,0., //group ether 
				  0.,  //group acid
				  2.,   //group ACNO2
				  0.,0.,0., //group NO3
				  0.,0.,0., //group CO-OH
  			          0.,0.,0.,0.,0.,0.,0.,0.,0., //group CO-OC
			          0.,  //group PAN
			          0.,  //group CO-OOH
			          0.,  //group O=COC=O
			          0.,0.,0.}; //group CHxNO2
  
  size = sizeof(group_tmp_AR0138)/sizeof(double);
  assert(size = 60);
	
  for(int i = 0; i < size; ++i)
	AR0138.groups[i] = group_tmp_AR0138[i];

  // Search the species name in the aerosol species list 
  // and add the species if its name matches with
  // the given list.
  add_species_ssh(surrogate, AR0138, species_list_aer, molecular_weight_aer,
              accomodation_coefficient,diffusion_coef,saturation_vapor_pressure,enthalpy_vaporization, henry, t_ref, mass_density, species_part,nlayer,i_hydrophilic,
		  N_inert, N_inorganic);



/* ==== AR0140OOH ==== */ 

  species AR0140OOH;
  AR0140OOH.name="AR0140OOH";
  AR0140OOH.is_inorganic_precursor=false;
  AR0140OOH.Psat_ref=3.809e-6; 		// Saturation vapor pressure at Tref (torr)
  AR0140OOH.Tref=298;         		// Temperature of reference (K)
  AR0140OOH.deltaH=110.39;     		// Enthalpy of vaporization (kJ/mol)
  AR0140OOH.Henry=0.0;     		// Henry's law constant at Tref (M/atm)
  AR0140OOH.aq_type="none"; 		// "none","diacid","monoacid" or "aldehyde"
  // AR0140OOH.Kacidity1=XXacidityXX;   		// First acidity constant
  AR0140OOH.hydrophilic=true;   	// Does the species condense on the aqueous phase?
  AR0140OOH.hydrophobic=true;  		// Does the species condense on the organic phase?
  AR0140OOH.nonvolatile=false; 		// Is the compound nonvolatile?
  AR0140OOH.kp_from_experiment=false; 	// Use experimental partitioning constant at Tref?
  AR0140OOH.is_organic=true;  		// Is the compound organic?
  AR0140OOH.compute_gamma_org=true;  	// Compute the activity coefficients of the organic phase for this compound?
  AR0140OOH.compute_gamma_aq=true;  	// Compute the activity coefficients of the aqueous phase for this compound?
  AR0140OOH.Koligo_org=0.0;      	//oligomeriation constant in the organic phase
  AR0140OOH.rho=1300.;  
  AR0140OOH.is_monomer=false;
  AR0140OOH.rion=false;
  AR0140OOH.KDiffusion_air=1.0e-5;
  //  AR0140OOH.accomodation_coefficient=alpha;
  AR0140OOH.viscosity=1.68e12;
  AR0140OOH.is_solid=false;
  AR0140OOH.is_generic=false;
  
  //Group: if no functionnal group in the species use the default species
  //for the computation of activity coefficients
  //
  double group_tmp_AR0140OOH [] = {0.,0.,0.,0., // group C
				  0.,0.,0.,0., //group C[OH]
				  0.,0.,0.,0., //group Calcohol
				  0.,0.,0.,0., //group Calcohol-tail
				  0.,0.,0.,0.,1., //group C=C
				  0.,0., //group aromatic carbon (AC)
				  0.,0.,0., // group //AC-C
				  1.,  //group OH
				  0., //group H2O
				  0., //group ACOH
				  0.,2., //group ketone
				  0.,   //group aldehyde  
				  0.,0., //group ester
				  0.,0.,0., //group ether 
				  0.,  //group acid
				  0.,   //group ACNO2
				  0.,0.,0., //group NO3
				  1.,0.,0., //group CO-OH
  			          0.,0.,0.,0.,0.,0.,0.,0.,0., //group CO-OC
			          0.,  //group PAN
			          0.,  //group CO-OOH
			          0.,  //group O=COC=O
			          0.,0.,0.}; //group CHxNO2
  
  size = sizeof(group_tmp_AR0140OOH)/sizeof(double);
  assert(size = 60);
	
  for(int i = 0; i < size; ++i)
	AR0140OOH.groups[i] = group_tmp_AR0140OOH[i];

  // Search the species name in the aerosol species list 
  // and add the species if its name matches with
  // the given list.
  add_species_ssh(surrogate, AR0140OOH, species_list_aer, molecular_weight_aer,
              accomodation_coefficient,diffusion_coef,saturation_vapor_pressure,enthalpy_vaporization, henry, t_ref, mass_density, species_part,nlayer,i_hydrophilic,
		  N_inert, N_inorganic);



/* ==== AR0144 ==== */ 

  species AR0144;
  AR0144.name="AR0144";
  AR0144.is_inorganic_precursor=false;
  AR0144.Psat_ref=1.488e-4; 		// Saturation vapor pressure at Tref (torr)
  AR0144.Tref=298;         		// Temperature of reference (K)
  AR0144.deltaH=91.45;     		// Enthalpy of vaporization (kJ/mol)
  AR0144.Henry=0.0;     		// Henry's law constant at Tref (M/atm)
  AR0144.aq_type="none"; 		// "none","diacid","monoacid" or "aldehyde"
  // AR0144.Kacidity1=XXacidityXX;   		// First acidity constant
  AR0144.hydrophilic=true;   	// Does the species condense on the aqueous phase?
  AR0144.hydrophobic=true;  		// Does the species condense on the organic phase?
  AR0144.nonvolatile=false; 		// Is the compound nonvolatile?
  AR0144.kp_from_experiment=false; 	// Use experimental partitioning constant at Tref?
  AR0144.is_organic=true;  		// Is the compound organic?
  AR0144.compute_gamma_org=true;  	// Compute the activity coefficients of the organic phase for this compound?
  AR0144.compute_gamma_aq=true;  	// Compute the activity coefficients of the aqueous phase for this compound?
  AR0144.Koligo_org=0.0;      	//oligomeriation constant in the organic phase
  AR0144.rho=1300.;  
  AR0144.is_monomer=false;
  AR0144.rion=false;
  AR0144.KDiffusion_air=1.0e-5;
  //  AR0144.accomodation_coefficient=alpha;
  AR0144.viscosity=1.68e12;
  AR0144.is_solid=false;
  AR0144.is_generic=false;
  
  //Group: if no functionnal group in the species use the default species
  //for the computation of activity coefficients
  //
  double group_tmp_AR0144 [] = {0.,0.,0.,0., // group C
				  0.,0.,0.,0., //group C[OH]
				  0.,0.,0.,0., //group Calcohol
				  0.,0.,0.,0., //group Calcohol-tail
				  0.,0.,0.,0.,0., //group C=C
				  4.,1., //group aromatic carbon (AC)
				  0.,0.,0., // group //AC-C
				  0.,  //group OH
				  0., //group H2O
				  1., //group ACOH
				  0.,0., //group ketone
				  0.,   //group aldehyde  
				  0.,0., //group ester
				  0.,0.,0., //group ether 
				  0.,  //group acid
				  0.,   //group ACNO2
				  0.,0.,0., //group NO3
				  1.,0.,0., //group CO-OH
  			          0.,0.,0.,0.,0.,0.,0.,0.,0., //group CO-OC
			          0.,  //group PAN
			          0.,  //group CO-OOH
			          0.,  //group O=COC=O
			          0.,0.,0.}; //group CHxNO2
  
  size = sizeof(group_tmp_AR0144)/sizeof(double);
  assert(size = 60);
	
  for(int i = 0; i < size; ++i)
	AR0144.groups[i] = group_tmp_AR0144[i];

  // Search the species name in the aerosol species list 
  // and add the species if its name matches with
  // the given list.
  add_species_ssh(surrogate, AR0144, species_list_aer, molecular_weight_aer,
              accomodation_coefficient,diffusion_coef,saturation_vapor_pressure,enthalpy_vaporization, henry, t_ref, mass_density, species_part,nlayer,i_hydrophilic,
		  N_inert, N_inorganic);



/* ==== AR0153 ==== */ 

  species AR0153;
  AR0153.name="AR0153";
  AR0153.is_inorganic_precursor=false;
  AR0153.Psat_ref=2.00e-10; 		// Saturation vapor pressure at Tref (torr)
  AR0153.Tref=298;         		// Temperature of reference (K)
  AR0153.deltaH=159.39;     		// Enthalpy of vaporization (kJ/mol)
  AR0153.Henry=0.0;     		// Henry's law constant at Tref (M/atm)
  AR0153.aq_type="none"; 		// "none","diacid","monoacid" or "aldehyde"
  // AR0153.Kacidity1=XXacidityXX;   		// First acidity constant
  AR0153.hydrophilic=true;   	// Does the species condense on the aqueous phase?
  AR0153.hydrophobic=true;  		// Does the species condense on the organic phase?
  AR0153.nonvolatile=false; 		// Is the compound nonvolatile?
  AR0153.kp_from_experiment=false; 	// Use experimental partitioning constant at Tref?
  AR0153.is_organic=true;  		// Is the compound organic?
  AR0153.compute_gamma_org=true;  	// Compute the activity coefficients of the organic phase for this compound?
  AR0153.compute_gamma_aq=true;  	// Compute the activity coefficients of the aqueous phase for this compound?
  AR0153.Koligo_org=0.0;      	//oligomeriation constant in the organic phase
  AR0153.rho=1300.;  
  AR0153.is_monomer=false;
  AR0153.rion=false;
  AR0153.KDiffusion_air=1.0e-5;
  //  AR0153.accomodation_coefficient=alpha;
  AR0153.viscosity=1.68e12;
  AR0153.is_solid=false;
  AR0153.is_generic=false;
  
  //Group: if no functionnal group in the species use the default species
  //for the computation of activity coefficients
  //
  double group_tmp_AR0153 [] = {1.,0.,0.,0., // group C
				  0.,0.,0.,0., //group C[OH]
				  0.,0.,0.,0., //group Calcohol
				  0.,0.,0.,0., //group Calcohol-tail
				  0.,0.,0.,0.,1., //group C=C
				  0.,0., //group aromatic carbon (AC)
				  0.,0.,0., // group //AC-C
				  3.,  //group OH
				  0., //group H2O
				  0., //group ACOH
				  0.,0., //group ketone
				  0.,   //group aldehyde  
				  0.,0., //group ester
				  0.,0.,0., //group ether 
				  0.,  //group acid
				  0.,   //group ACNO2
				  0.,0.,0., //group NO3
				  0.,0.,1., //group CO-OH
  			          0.,0.,0.,0.,0.,0.,0.,1.,0., //group CO-OC
			          0.,  //group PAN
			          0.,  //group CO-OOH
			          0.,  //group O=COC=O
			          0.,0.,1.}; //group CHxNO2
  
  size = sizeof(group_tmp_AR0153)/sizeof(double);
  assert(size = 60);
	
  for(int i = 0; i < size; ++i)
	AR0153.groups[i] = group_tmp_AR0153[i];

  // Search the species name in the aerosol species list 
  // and add the species if its name matches with
  // the given list.
  add_species_ssh(surrogate, AR0153, species_list_aer, molecular_weight_aer,
              accomodation_coefficient,diffusion_coef,saturation_vapor_pressure,enthalpy_vaporization, henry, t_ref, mass_density, species_part,nlayer,i_hydrophilic,
		  N_inert, N_inorganic);



/* ==== AU4000 ==== */ 

  species AU4000;
  AU4000.name="AU4000";
  AU4000.is_inorganic_precursor=false;
  AU4000.Psat_ref=0.003809; 		// Saturation vapor pressure at Tref (torr)
  AU4000.Tref=298;         		// Temperature of reference (K)
  AU4000.deltaH=77.72;     		// Enthalpy of vaporization (kJ/mol)
  AU4000.Henry=0.0;     		// Henry's law constant at Tref (M/atm)
  AU4000.aq_type="none"; 		// "none","diacid","monoacid" or "aldehyde"
  // AU4000.Kacidity1=XXacidityXX;   		// First acidity constant
  AU4000.hydrophilic=true;   	// Does the species condense on the aqueous phase?
  AU4000.hydrophobic=true;  		// Does the species condense on the organic phase?
  AU4000.nonvolatile=false; 		// Is the compound nonvolatile?
  AU4000.kp_from_experiment=false; 	// Use experimental partitioning constant at Tref?
  AU4000.is_organic=true;  		// Is the compound organic?
  AU4000.compute_gamma_org=true;  	// Compute the activity coefficients of the organic phase for this compound?
  AU4000.compute_gamma_aq=true;  	// Compute the activity coefficients of the aqueous phase for this compound?
  AU4000.Koligo_org=0.0;      	//oligomeriation constant in the organic phase
  AU4000.rho=1300.;  
  AU4000.is_monomer=false;
  AU4000.rion=false;
  AU4000.KDiffusion_air=1.0e-5;
  //  AU4000.accomodation_coefficient=alpha;
  AU4000.viscosity=1.68e12;
  AU4000.is_solid=false;
  AU4000.is_generic=false;
  
  //Group: if no functionnal group in the species use the default species
  //for the computation of activity coefficients
  //
  double group_tmp_AU4000 [] = {0.,0.,0.,0., // group C
				  0.,0.,0.,0., //group C[OH]
				  0.,0.,0.,0., //group Calcohol
				  0.,0.,0.,0., //group Calcohol-tail
				  0.,1.,0.,0.,0., //group C=C
				  0.,0., //group aromatic carbon (AC)
				  0.,0.,0., // group //AC-C
				  0.,  //group OH
				  0., //group H2O
				  0., //group ACOH
				  0.,0., //group ketone
				  1.,   //group aldehyde  
				  0.,0., //group ester
				  0.,0.,0., //group ether 
				  1.,  //group acid
				  0.,   //group ACNO2
				  0.,0.,0., //group NO3
				  0.,0.,0., //group CO-OH
  			          0.,0.,0.,0.,0.,0.,0.,0.,0., //group CO-OC
			          0.,  //group PAN
			          0.,  //group CO-OOH
			          0.,  //group O=COC=O
			          0.,0.,0.}; //group CHxNO2
  
  size = sizeof(group_tmp_AU4000)/sizeof(double);
  assert(size = 60);
	
  for(int i = 0; i < size; ++i)
	AU4000.groups[i] = group_tmp_AU4000[i];

  // Search the species name in the aerosol species list 
  // and add the species if its name matches with
  // the given list.
  add_species_ssh(surrogate, AU4000, species_list_aer, molecular_weight_aer,
              accomodation_coefficient,diffusion_coef,saturation_vapor_pressure,enthalpy_vaporization, henry, t_ref, mass_density, species_part,nlayer,i_hydrophilic,
		  N_inert, N_inorganic);



/* ==== AU5000 ==== */ 

  species AU5000;
  AU5000.name="AU5000";
  AU5000.is_inorganic_precursor=false;
  AU5000.Psat_ref=0.001516; 		// Saturation vapor pressure at Tref (torr)
  AU5000.Tref=298;         		// Temperature of reference (K)
  AU5000.deltaH=82.07;     		// Enthalpy of vaporization (kJ/mol)
  AU5000.Henry=0.0;     		// Henry's law constant at Tref (M/atm)
  AU5000.aq_type="none"; 		// "none","diacid","monoacid" or "aldehyde"
  // AU5000.Kacidity1=XXacidityXX;   		// First acidity constant
  AU5000.hydrophilic=true;   	// Does the species condense on the aqueous phase?
  AU5000.hydrophobic=true;  		// Does the species condense on the organic phase?
  AU5000.nonvolatile=false; 		// Is the compound nonvolatile?
  AU5000.kp_from_experiment=false; 	// Use experimental partitioning constant at Tref?
  AU5000.is_organic=true;  		// Is the compound organic?
  AU5000.compute_gamma_org=true;  	// Compute the activity coefficients of the organic phase for this compound?
  AU5000.compute_gamma_aq=true;  	// Compute the activity coefficients of the aqueous phase for this compound?
  AU5000.Koligo_org=0.0;      	//oligomeriation constant in the organic phase
  AU5000.rho=1300.;  
  AU5000.is_monomer=false;
  AU5000.rion=false;
  AU5000.KDiffusion_air=1.0e-5;
  //  AU5000.accomodation_coefficient=alpha;
  AU5000.viscosity=1.68e12;
  AU5000.is_solid=false;
  AU5000.is_generic=false;
  
  //Group: if no functionnal group in the species use the default species
  //for the computation of activity coefficients
  //
  double group_tmp_AU5000 [] = {1.,0.,0.,0., // group C
				  0.,0.,0.,0., //group C[OH]
				  0.,0.,0.,0., //group Calcohol
				  0.,0.,0.,0., //group Calcohol-tail
				  0.,0.,0.,1.,0., //group C=C
				  0.,0., //group aromatic carbon (AC)
				  0.,0.,0., // group //AC-C
				  0.,  //group OH
				  0., //group H2O
				  0., //group ACOH
				  0.,0., //group ketone
				  1.,   //group aldehyde  
				  0.,0., //group ester
				  0.,0.,0., //group ether 
				  1.,  //group acid
				  0.,   //group ACNO2
				  0.,0.,0., //group NO3
				  0.,0.,0., //group CO-OH
  			          0.,0.,0.,0.,0.,0.,0.,0.,0., //group CO-OC
			          0.,  //group PAN
			          0.,  //group CO-OOH
			          0.,  //group O=COC=O
			          0.,0.,0.}; //group CHxNO2
  
  size = sizeof(group_tmp_AU5000)/sizeof(double);
  assert(size = 60);
	
  for(int i = 0; i < size; ++i)
	AU5000.groups[i] = group_tmp_AU5000[i];

  // Search the species name in the aerosol species list 
  // and add the species if its name matches with
  // the given list.
  add_species_ssh(surrogate, AU5000, species_list_aer, molecular_weight_aer,
              accomodation_coefficient,diffusion_coef,saturation_vapor_pressure,enthalpy_vaporization, henry, t_ref, mass_density, species_part,nlayer,i_hydrophilic,
		  N_inert, N_inorganic);



/* ==== AU5002 ==== */ 

  species AU5002;
  AU5002.name="AU5002";
  AU5002.is_inorganic_precursor=false;
  AU5002.Psat_ref=1.910e-2; 		// Saturation vapor pressure at Tref (torr)
  AU5002.Tref=298;         		// Temperature of reference (K)
  AU5002.deltaH=70.63;     		// Enthalpy of vaporization (kJ/mol)
  AU5002.Henry=0.0;     		// Henry's law constant at Tref (M/atm)
  AU5002.aq_type="none"; 		// "none","diacid","monoacid" or "aldehyde"
  // AU5002.Kacidity1=XXacidityXX;   		// First acidity constant
  AU5002.hydrophilic=true;   	// Does the species condense on the aqueous phase?
  AU5002.hydrophobic=true;  		// Does the species condense on the organic phase?
  AU5002.nonvolatile=false; 		// Is the compound nonvolatile?
  AU5002.kp_from_experiment=false; 	// Use experimental partitioning constant at Tref?
  AU5002.is_organic=true;  		// Is the compound organic?
  AU5002.compute_gamma_org=true;  	// Compute the activity coefficients of the organic phase for this compound?
  AU5002.compute_gamma_aq=true;  	// Compute the activity coefficients of the aqueous phase for this compound?
  AU5002.Koligo_org=0.0;      	//oligomeriation constant in the organic phase
  AU5002.rho=1300.;  
  AU5002.is_monomer=false;
  AU5002.rion=false;
  AU5002.KDiffusion_air=1.0e-5;
  //  AU5002.accomodation_coefficient=alpha;
  AU5002.viscosity=1.68e12;
  AU5002.is_solid=false;
  AU5002.is_generic=false;
  
  //Group: if no functionnal group in the species use the default species
  //for the computation of activity coefficients
  //
  double group_tmp_AU5002 [] = {0.,0.,0.,0., // group C
				  0.,0.,0.,0., //group C[OH]
				  0.,0.,0.,0., //group Calcohol
				  0.,0.,0.,0., //group Calcohol-tail
				  0.,1.,0.,0.,0., //group C=C
				  0.,0., //group aromatic carbon (AC)
				  0.,0.,0., // group //AC-C
				  0.,  //group OH
				  0., //group H2O
				  0., //group ACOH
				  1.,0., //group ketone
				  0.,   //group aldehyde  
				  0.,0., //group ester
				  0.,0.,0., //group ether 
				  1.,  //group acid
				  0.,   //group ACNO2
				  0.,0.,0., //group NO3
				  0.,0.,0., //group CO-OH
  			          0.,0.,0.,0.,0.,0.,0.,0.,0., //group CO-OC
			          0.,  //group PAN
			          0.,  //group CO-OOH
			          0.,  //group O=COC=O
			          0.,0.,0.}; //group CHxNO2
  
  size = sizeof(group_tmp_AU5002)/sizeof(double);
  assert(size = 60);
	
  for(int i = 0; i < size; ++i)
	AU5002.groups[i] = group_tmp_AU5002[i];

  // Search the species name in the aerosol species list 
  // and add the species if its name matches with
  // the given list.
  add_species_ssh(surrogate, AU5002, species_list_aer, molecular_weight_aer,
              accomodation_coefficient,diffusion_coef,saturation_vapor_pressure,enthalpy_vaporization, henry, t_ref, mass_density, species_part,nlayer,i_hydrophilic,
		  N_inert, N_inorganic);



/* ==== AU50DN ==== */ 

  species AU50DN;
  AU50DN.name="AU50DN";
  AU50DN.is_inorganic_precursor=false;
  AU50DN.Psat_ref=9.568e-5; 		// Saturation vapor pressure at Tref (torr)
  AU50DN.Tref=298;         		// Temperature of reference (K)
  AU50DN.deltaH=87.5;     		// Enthalpy of vaporization (kJ/mol)
  AU50DN.Henry=0.0;     		// Henry's law constant at Tref (M/atm)
  AU50DN.aq_type="none"; 		// "none","diacid","monoacid" or "aldehyde"
  // AU50DN.Kacidity1=XXacidityXX;   		// First acidity constant
  AU50DN.hydrophilic=true;   	// Does the species condense on the aqueous phase?
  AU50DN.hydrophobic=true;  		// Does the species condense on the organic phase?
  AU50DN.nonvolatile=false; 		// Is the compound nonvolatile?
  AU50DN.kp_from_experiment=false; 	// Use experimental partitioning constant at Tref?
  AU50DN.is_organic=true;  		// Is the compound organic?
  AU50DN.compute_gamma_org=true;  	// Compute the activity coefficients of the organic phase for this compound?
  AU50DN.compute_gamma_aq=true;  	// Compute the activity coefficients of the aqueous phase for this compound?
  AU50DN.Koligo_org=0.0;      	//oligomeriation constant in the organic phase
  AU50DN.rho=1300.;  
  AU50DN.is_monomer=false;
  AU50DN.rion=false;
  AU50DN.KDiffusion_air=1.0e-5;
  //  AU50DN.accomodation_coefficient=alpha;
  AU50DN.viscosity=1.68e12;
  AU50DN.is_solid=false;
  AU50DN.is_generic=false;
  
  //Group: if no functionnal group in the species use the default species
  //for the computation of activity coefficients
  //
  double group_tmp_AU50DN [] = {0.,0.,0.,0., // group C
				  0.,0.,0.,0., //group C[OH]
				  0.,0.,0.,0., //group Calcohol
				  0.,0.,0.,0., //group Calcohol-tail
				  0.,0.,0.,1.,0., //group C=C
				  0.,0., //group aromatic carbon (AC)
				  0.,0.,0., // group //AC-C
				  0.,  //group OH
				  0., //group H2O
				  0., //group ACOH
				  0.,0., //group ketone
				  1.,   //group aldehyde  
				  0.,0., //group ester
				  0.,0.,0., //group ether 
				  1.,  //group acid
				  0.,   //group ACNO2
				  0.,0.,0., //group NO3
				  0.,0.,0., //group CO-OH
  			          0.,0.,0.,0.,0.,0.,0.,0.,0., //group CO-OC
			          0.,  //group PAN
			          0.,  //group CO-OOH
			          0.,  //group O=COC=O
			          0.,1.,0.}; //group CHxNO2
  
  size = sizeof(group_tmp_AU50DN)/sizeof(double);
  assert(size = 60);
	
  for(int i = 0; i < size; ++i)
	AU50DN.groups[i] = group_tmp_AU50DN[i];

  // Search the species name in the aerosol species list 
  // and add the species if its name matches with
  // the given list.
  add_species_ssh(surrogate, AU50DN, species_list_aer, molecular_weight_aer,
              accomodation_coefficient,diffusion_coef,saturation_vapor_pressure,enthalpy_vaporization, henry, t_ref, mass_density, species_part,nlayer,i_hydrophilic,
		  N_inert, N_inorganic);



/* ==== AU6000 ==== */ 

  species AU6000;
  AU6000.name="AU6000";
  AU6000.is_inorganic_precursor=false;
  AU6000.Psat_ref=2.911e-4; 		// Saturation vapor pressure at Tref (torr)
  AU6000.Tref=298;         		// Temperature of reference (K)
  AU6000.deltaH=88.937;     		// Enthalpy of vaporization (kJ/mol)
  AU6000.Henry=0.0;     		// Henry's law constant at Tref (M/atm)
  AU6000.aq_type="none"; 		// "none","diacid","monoacid" or "aldehyde"
  // AU6000.Kacidity1=XXacidityXX;   		// First acidity constant
  AU6000.hydrophilic=true;   	// Does the species condense on the aqueous phase?
  AU6000.hydrophobic=true;  		// Does the species condense on the organic phase?
  AU6000.nonvolatile=false; 		// Is the compound nonvolatile?
  AU6000.kp_from_experiment=false; 	// Use experimental partitioning constant at Tref?
  AU6000.is_organic=true;  		// Is the compound organic?
  AU6000.compute_gamma_org=true;  	// Compute the activity coefficients of the organic phase for this compound?
  AU6000.compute_gamma_aq=true;  	// Compute the activity coefficients of the aqueous phase for this compound?
  AU6000.Koligo_org=0.0;      	//oligomeriation constant in the organic phase
  AU6000.rho=1300.;  
  AU6000.is_monomer=false;
  AU6000.rion=false;
  AU6000.KDiffusion_air=1.0e-5;
  //  AU6000.accomodation_coefficient=alpha;
  AU6000.viscosity=1.68e12;
  AU6000.is_solid=false;
  AU6000.is_generic=false;
  
  //Group: if no functionnal group in the species use the default species
  //for the computation of activity coefficients
  //
  double group_tmp_AU6000 [] = {1.,0.,0.,0., // group C
				  0.,0.,0.,0., //group C[OH]
				  0.,0.,0.,0., //group Calcohol
				  0.,0.,0.,0., //group Calcohol-tail
				  0.,2.,0.,0.,0., //group C=C
				  0.,0., //group aromatic carbon (AC)
				  0.,0.,0., // group //AC-C
				  0.,  //group OH
				  0., //group H2O
				  0., //group ACOH
				  0.,0., //group ketone
				  0.,   //group aldehyde  
				  0.,0., //group ester
				  0.,0.,0., //group ether 
				  1.,  //group acid
				  0.,   //group ACNO2
				  0.,0.,0., //group NO3
				  0.,0.,0., //group CO-OH
  			          0.,0.,0.,0.,0.,0.,0.,0.,0., //group CO-OC
			          0.,  //group PAN
			          0.,  //group CO-OOH
			          0.,  //group O=COC=O
			          0.,0.,0.}; //group CHxNO2
  
  size = sizeof(group_tmp_AU6000)/sizeof(double);
  assert(size = 60);
	
  for(int i = 0; i < size; ++i)
	AU6000.groups[i] = group_tmp_AU6000[i];

  // Search the species name in the aerosol species list 
  // and add the species if its name matches with
  // the given list.
  add_species_ssh(surrogate, AU6000, species_list_aer, molecular_weight_aer,
              accomodation_coefficient,diffusion_coef,saturation_vapor_pressure,enthalpy_vaporization, henry, t_ref, mass_density, species_part,nlayer,i_hydrophilic,
		  N_inert, N_inorganic);



/* ==== AU7000 ==== */ 

  species AU7000;
  AU7000.name="AU7000";
  AU7000.is_inorganic_precursor=false;
  AU7000.Psat_ref=3.969e-7; 		// Saturation vapor pressure at Tref (torr)
  AU7000.Tref=298;         		// Temperature of reference (K)
  AU7000.deltaH=120.38;     		// Enthalpy of vaporization (kJ/mol)
  AU7000.Henry=0.0;     		// Henry's law constant at Tref (M/atm)
  AU7000.aq_type="none"; 		// "none","diacid","monoacid" or "aldehyde"
  // AU7000.Kacidity1=XXacidityXX;   		// First acidity constant
  AU7000.hydrophilic=true;   	// Does the species condense on the aqueous phase?
  AU7000.hydrophobic=true;  		// Does the species condense on the organic phase?
  AU7000.nonvolatile=false; 		// Is the compound nonvolatile?
  AU7000.kp_from_experiment=false; 	// Use experimental partitioning constant at Tref?
  AU7000.is_organic=true;  		// Is the compound organic?
  AU7000.compute_gamma_org=true;  	// Compute the activity coefficients of the organic phase for this compound?
  AU7000.compute_gamma_aq=true;  	// Compute the activity coefficients of the aqueous phase for this compound?
  AU7000.Koligo_org=0.0;      	//oligomeriation constant in the organic phase
  AU7000.rho=1300.;  
  AU7000.is_monomer=false;
  AU7000.rion=false;
  AU7000.KDiffusion_air=1.0e-5;
  //  AU7000.accomodation_coefficient=alpha;
  AU7000.viscosity=1.68e12;
  AU7000.is_solid=false;
  AU7000.is_generic=false;
  
  //Group: if no functionnal group in the species use the default species
  //for the computation of activity coefficients
  //
  double group_tmp_AU7000 [] = {0.,0.,0.,0., // group C
				  0.,0.,0.,0., //group C[OH]
				  0.,0.,0.,0., //group Calcohol
				  0.,0.,0.,0., //group Calcohol-tail
				  0.,2.,0.,0.,0., //group C=C
				  0.,0., //group aromatic carbon (AC)
				  0.,0.,0., // group //AC-C
				  0.,  //group OH
				  0., //group H2O
				  0., //group ACOH
				  1.,0., //group ketone
				  0.,   //group aldehyde  
				  0.,0., //group ester
				  0.,0.,0., //group ether 
				  1.,  //group acid
				  0.,   //group ACNO2
				  0.,0.,0., //group NO3
				  0.,0.,0., //group CO-OH
  			          0.,0.,0.,0.,0.,0.,0.,0.,0., //group CO-OC
			          0.,  //group PAN
			          0.,  //group CO-OOH
			          0.,  //group O=COC=O
			          0.,0.,0.}; //group CHxNO2
  
  size = sizeof(group_tmp_AU7000)/sizeof(double);
  assert(size = 60);
	
  for(int i = 0; i < size; ++i)
	AU7000.groups[i] = group_tmp_AU7000[i];

  // Search the species name in the aerosol species list 
  // and add the species if its name matches with
  // the given list.
  add_species_ssh(surrogate, AU7000, species_list_aer, molecular_weight_aer,
              accomodation_coefficient,diffusion_coef,saturation_vapor_pressure,enthalpy_vaporization, henry, t_ref, mass_density, species_part,nlayer,i_hydrophilic,
		  N_inert, N_inorganic);



/* ==== BTOL3OHOOH ==== */ 

  species BTOL3OHOOH;
  BTOL3OHOOH.name="BTOL3OHOOH";
  BTOL3OHOOH.is_inorganic_precursor=false;
  BTOL3OHOOH.Psat_ref=3.809e-8; 		// Saturation vapor pressure at Tref (torr)
  BTOL3OHOOH.Tref=298;         		// Temperature of reference (K)
  BTOL3OHOOH.deltaH=132.17;     		// Enthalpy of vaporization (kJ/mol)
  BTOL3OHOOH.Henry=0.0;     		// Henry's law constant at Tref (M/atm)
  BTOL3OHOOH.aq_type="none"; 		// "none","diacid","monoacid" or "aldehyde"
  // BTOL3OHOOH.Kacidity1=XXacidityXX;   		// First acidity constant
  BTOL3OHOOH.hydrophilic=true;   	// Does the species condense on the aqueous phase?
  BTOL3OHOOH.hydrophobic=true;  		// Does the species condense on the organic phase?
  BTOL3OHOOH.nonvolatile=false; 		// Is the compound nonvolatile?
  BTOL3OHOOH.kp_from_experiment=false; 	// Use experimental partitioning constant at Tref?
  BTOL3OHOOH.is_organic=true;  		// Is the compound organic?
  BTOL3OHOOH.compute_gamma_org=true;  	// Compute the activity coefficients of the organic phase for this compound?
  BTOL3OHOOH.compute_gamma_aq=true;  	// Compute the activity coefficients of the aqueous phase for this compound?
  BTOL3OHOOH.Koligo_org=0.0;      	//oligomeriation constant in the organic phase
  BTOL3OHOOH.rho=1300.;  
  BTOL3OHOOH.is_monomer=false;
  BTOL3OHOOH.rion=false;
  BTOL3OHOOH.KDiffusion_air=1.0e-5;
  //  BTOL3OHOOH.accomodation_coefficient=alpha;
  BTOL3OHOOH.viscosity=1.68e12;
  BTOL3OHOOH.is_solid=false;
  BTOL3OHOOH.is_generic=false;
  
  //Group: if no functionnal group in the species use the default species
  //for the computation of activity coefficients
  //
  double group_tmp_BTOL3OHOOH [] = {1.,0.,0.,0., // group C
				  0.,0.,1.,0., //group C[OH]
				  0.,0.,0.,0., //group Calcohol
				  0.,0.,0.,0., //group Calcohol-tail
				  0.,0.,0.,0.,1., //group C=C
				  0.,0., //group aromatic carbon (AC)
				  0.,0.,0., // group //AC-C
				  3.,  //group OH
				  0., //group H2O
				  0., //group ACOH
				  0.,0., //group ketone
				  0.,   //group aldehyde  
				  0.,0., //group ester
				  0.,0.,0., //group ether 
				  0.,  //group acid
				  0.,   //group ACNO2
				  0.,0.,0., //group NO3
				  0.,1.,0., //group CO-OH
  			          0.,0.,0.,0.,0.,0.,0.,1.,0., //group CO-OC
			          0.,  //group PAN
			          0.,  //group CO-OOH
			          0.,  //group O=COC=O
			          0.,0.,0.}; //group CHxNO2
  
  size = sizeof(group_tmp_BTOL3OHOOH)/sizeof(double);
  assert(size = 60);
	
  for(int i = 0; i < size; ++i)
	BTOL3OHOOH.groups[i] = group_tmp_BTOL3OHOOH[i];

  // Search the species name in the aerosol species list 
  // and add the species if its name matches with
  // the given list.
  add_species_ssh(surrogate, BTOL3OHOOH, species_list_aer, molecular_weight_aer,
              accomodation_coefficient,diffusion_coef,saturation_vapor_pressure,enthalpy_vaporization, henry, t_ref, mass_density, species_part,nlayer,i_hydrophilic,
		  N_inert, N_inorganic);



/* ==== BTOL4OHOOH ==== */ 

  species BTOL4OHOOH;
  BTOL4OHOOH.name="BTOL4OHOOH";
  BTOL4OHOOH.is_inorganic_precursor=false;
  BTOL4OHOOH.Psat_ref=2.403e-10; 		// Saturation vapor pressure at Tref (torr)
  BTOL4OHOOH.Tref=298;         		// Temperature of reference (K)
  BTOL4OHOOH.deltaH=156.13;     		// Enthalpy of vaporization (kJ/mol)
  BTOL4OHOOH.Henry=0.0;     		// Henry's law constant at Tref (M/atm)
  BTOL4OHOOH.aq_type="none"; 		// "none","diacid","monoacid" or "aldehyde"
  // BTOL4OHOOH.Kacidity1=XXacidityXX;   		// First acidity constant
  BTOL4OHOOH.hydrophilic=true;   	// Does the species condense on the aqueous phase?
  BTOL4OHOOH.hydrophobic=true;  		// Does the species condense on the organic phase?
  BTOL4OHOOH.nonvolatile=false; 		// Is the compound nonvolatile?
  BTOL4OHOOH.kp_from_experiment=false; 	// Use experimental partitioning constant at Tref?
  BTOL4OHOOH.is_organic=true;  		// Is the compound organic?
  BTOL4OHOOH.compute_gamma_org=true;  	// Compute the activity coefficients of the organic phase for this compound?
  BTOL4OHOOH.compute_gamma_aq=true;  	// Compute the activity coefficients of the aqueous phase for this compound?
  BTOL4OHOOH.Koligo_org=0.0;      	//oligomeriation constant in the organic phase
  BTOL4OHOOH.rho=1300.;  
  BTOL4OHOOH.is_monomer=false;
  BTOL4OHOOH.rion=false;
  BTOL4OHOOH.KDiffusion_air=1.0e-5;
  //  BTOL4OHOOH.accomodation_coefficient=alpha;
  BTOL4OHOOH.viscosity=1.68e12;
  BTOL4OHOOH.is_solid=false;
  BTOL4OHOOH.is_generic=false;
  
  //Group: if no functionnal group in the species use the default species
  //for the computation of activity coefficients
  //
  double group_tmp_BTOL4OHOOH [] = {1.,0.,0.,0., // group C
				  0.,0.,1.,0., //group C[OH]
				  0.,0.,0.,0., //group Calcohol
				  0.,0.,0.,0., //group Calcohol-tail
				  0.,0.,0.,0.,1., //group C=C
				  0.,0., //group aromatic carbon (AC)
				  0.,0.,0., // group //AC-C
				  4.,  //group OH
				  0., //group H2O
				  0., //group ACOH
				  0.,0., //group ketone
				  0.,   //group aldehyde  
				  0.,0., //group ester
				  0.,0.,0., //group ether 
				  0.,  //group acid
				  0.,   //group ACNO2
				  0.,0.,0., //group NO3
				  0.,0.,1., //group CO-OH
  			          0.,0.,0.,0.,0.,0.,0.,1.,0., //group CO-OC
			          0.,  //group PAN
			          0.,  //group CO-OOH
			          0.,  //group O=COC=O
			          0.,0.,0.}; //group CHxNO2
  
  size = sizeof(group_tmp_BTOL4OHOOH)/sizeof(double);
  assert(size = 60);
	
  for(int i = 0; i < size; ++i)
	BTOL4OHOOH.groups[i] = group_tmp_BTOL4OHOOH[i];

  // Search the species name in the aerosol species list 
  // and add the species if its name matches with
  // the given list.
  add_species_ssh(surrogate, BTOL4OHOOH, species_list_aer, molecular_weight_aer,
              accomodation_coefficient,diffusion_coef,saturation_vapor_pressure,enthalpy_vaporization, henry, t_ref, mass_density, species_part,nlayer,i_hydrophilic,
		  N_inert, N_inorganic);



/* ==== BTOL5OHOOH ==== */ 

  species BTOL5OHOOH;
  BTOL5OHOOH.name="BTOL5OHOOH";
  BTOL5OHOOH.is_inorganic_precursor=false;
  BTOL5OHOOH.Psat_ref=6.036e-13; 		// Saturation vapor pressure at Tref (torr)
  BTOL5OHOOH.Tref=298;         		// Temperature of reference (K)
  BTOL5OHOOH.deltaH=184.44;     		// Enthalpy of vaporization (kJ/mol)
  BTOL5OHOOH.Henry=0.0;     		// Henry's law constant at Tref (M/atm)
  BTOL5OHOOH.aq_type="none"; 		// "none","diacid","monoacid" or "aldehyde"
  // BTOL5OHOOH.Kacidity1=XXacidityXX;   		// First acidity constant
  BTOL5OHOOH.hydrophilic=true;   	// Does the species condense on the aqueous phase?
  BTOL5OHOOH.hydrophobic=true;  		// Does the species condense on the organic phase?
  BTOL5OHOOH.nonvolatile=false; 		// Is the compound nonvolatile?
  BTOL5OHOOH.kp_from_experiment=false; 	// Use experimental partitioning constant at Tref?
  BTOL5OHOOH.is_organic=true;  		// Is the compound organic?
  BTOL5OHOOH.compute_gamma_org=true;  	// Compute the activity coefficients of the organic phase for this compound?
  BTOL5OHOOH.compute_gamma_aq=true;  	// Compute the activity coefficients of the aqueous phase for this compound?
  BTOL5OHOOH.Koligo_org=0.0;      	//oligomeriation constant in the organic phase
  BTOL5OHOOH.rho=1300.;  
  BTOL5OHOOH.is_monomer=false;
  BTOL5OHOOH.rion=false;
  BTOL5OHOOH.KDiffusion_air=1.0e-5;
  //  BTOL5OHOOH.accomodation_coefficient=alpha;
  BTOL5OHOOH.viscosity=1.68e12;
  BTOL5OHOOH.is_solid=false;
  BTOL5OHOOH.is_generic=false;
  
  //Group: if no functionnal group in the species use the default species
  //for the computation of activity coefficients
  //
  double group_tmp_BTOL5OHOOH [] = {1.,0.,0.,0., // group C
				  0.,0.,1.,0., //group C[OH]
				  0.,0.,0.,0., //group Calcohol
				  0.,0.,0.,0., //group Calcohol-tail
				  0.,0.,0.,0.,1., //group C=C
				  0.,0., //group aromatic carbon (AC)
				  0.,0.,0., // group //AC-C
				  5.,  //group OH
				  0., //group H2O
				  0., //group ACOH
				  0.,0., //group ketone
				  0.,   //group aldehyde  
				  0.,0., //group ester
				  0.,0.,0., //group ether 
				  0.,  //group acid
				  0.,   //group ACNO2
				  0.,0.,0., //group NO3
				  0.,0.,1., //group CO-OH
  			          0.,0.,0.,0.,0.,0.,0.,0.,1., //group CO-OC
			          0.,  //group PAN
			          0.,  //group CO-OOH
			          0.,  //group O=COC=O
			          0.,0.,0.}; //group CHxNO2
  
  size = sizeof(group_tmp_BTOL5OHOOH)/sizeof(double);
  assert(size = 60);
	
  for(int i = 0; i < size; ++i)
	BTOL5OHOOH.groups[i] = group_tmp_BTOL5OHOOH[i];

  // Search the species name in the aerosol species list 
  // and add the species if its name matches with
  // the given list.
  add_species_ssh(surrogate, BTOL5OHOOH, species_list_aer, molecular_weight_aer,
              accomodation_coefficient,diffusion_coef,saturation_vapor_pressure,enthalpy_vaporization, henry, t_ref, mass_density, species_part,nlayer,i_hydrophilic,
		  N_inert, N_inorganic);



/* ==== BZALDOH ==== */ 

  species BZALDOH;
  BZALDOH.name="BZALDOH";
  BZALDOH.is_inorganic_precursor=false;
  BZALDOH.Psat_ref=3.809e-2; 		// Saturation vapor pressure at Tref (torr)
  BZALDOH.Tref=298;         		// Temperature of reference (K)
  BZALDOH.deltaH=67.10;     		// Enthalpy of vaporization (kJ/mol)
  BZALDOH.Henry=0.0;     		// Henry's law constant at Tref (M/atm)
  BZALDOH.aq_type="none"; 		// "none","diacid","monoacid" or "aldehyde"
  // BZALDOH.Kacidity1=XXacidityXX;   		// First acidity constant
  BZALDOH.hydrophilic=true;   	// Does the species condense on the aqueous phase?
  BZALDOH.hydrophobic=true;  		// Does the species condense on the organic phase?
  BZALDOH.nonvolatile=false; 		// Is the compound nonvolatile?
  BZALDOH.kp_from_experiment=false; 	// Use experimental partitioning constant at Tref?
  BZALDOH.is_organic=true;  		// Is the compound organic?
  BZALDOH.compute_gamma_org=true;  	// Compute the activity coefficients of the organic phase for this compound?
  BZALDOH.compute_gamma_aq=true;  	// Compute the activity coefficients of the aqueous phase for this compound?
  BZALDOH.Koligo_org=0.0;      	//oligomeriation constant in the organic phase
  BZALDOH.rho=1300.;  
  BZALDOH.is_monomer=false;
  BZALDOH.rion=false;
  BZALDOH.KDiffusion_air=1.0e-5;
  //  BZALDOH.accomodation_coefficient=alpha;
  BZALDOH.viscosity=1.68e12;
  BZALDOH.is_solid=false;
  BZALDOH.is_generic=false;
  
  //Group: if no functionnal group in the species use the default species
  //for the computation of activity coefficients
  //
  double group_tmp_BZALDOH [] = {0.,0.,0.,0., // group C
				  0.,0.,0.,0., //group C[OH]
				  0.,0.,0.,0., //group Calcohol
				  0.,0.,0.,0., //group Calcohol-tail
				  0.,0.,0.,0.,0., //group C=C
				  4.,1., //group aromatic carbon (AC)
				  0.,0.,0., // group //AC-C
				  0.,  //group OH
				  0., //group H2O
				  1., //group ACOH
				  0.,0., //group ketone
				  1.,   //group aldehyde  
				  0.,0., //group ester
				  0.,0.,0., //group ether 
				  0.,  //group acid
				  0.,   //group ACNO2
				  0.,0.,0., //group NO3
				  0.,0.,0., //group CO-OH
  			          0.,0.,0.,0.,0.,0.,0.,0.,0., //group CO-OC
			          0.,  //group PAN
			          0.,  //group CO-OOH
			          0.,  //group O=COC=O
			          0.,0.,0.}; //group CHxNO2
  
  size = sizeof(group_tmp_BZALDOH)/sizeof(double);
  assert(size = 60);
	
  for(int i = 0; i < size; ++i)
	BZALDOH.groups[i] = group_tmp_BZALDOH[i];

  // Search the species name in the aerosol species list 
  // and add the species if its name matches with
  // the given list.
  add_species_ssh(surrogate, BZALDOH, species_list_aer, molecular_weight_aer,
              accomodation_coefficient,diffusion_coef,saturation_vapor_pressure,enthalpy_vaporization, henry, t_ref, mass_density, species_part,nlayer,i_hydrophilic,
		  N_inert, N_inorganic);



/* ==== C73K1OHOOH ==== */ 

  species C73K1OHOOH;
  C73K1OHOOH.name="C73K1OHOOH";
  C73K1OHOOH.is_inorganic_precursor=false;
  C73K1OHOOH.Psat_ref=4.795e-6; 		// Saturation vapor pressure at Tref (torr)
  C73K1OHOOH.Tref=298;         		// Temperature of reference (K)
  C73K1OHOOH.deltaH=109.30;     		// Enthalpy of vaporization (kJ/mol)
  C73K1OHOOH.Henry=0.0;     		// Henry's law constant at Tref (M/atm)
  C73K1OHOOH.aq_type="none"; 		// "none","diacid","monoacid" or "aldehyde"
  // C73K1OHOOH.Kacidity1=XXacidityXX;   		// First acidity constant
  C73K1OHOOH.hydrophilic=true;   	// Does the species condense on the aqueous phase?
  C73K1OHOOH.hydrophobic=true;  		// Does the species condense on the organic phase?
  C73K1OHOOH.nonvolatile=false; 		// Is the compound nonvolatile?
  C73K1OHOOH.kp_from_experiment=false; 	// Use experimental partitioning constant at Tref?
  C73K1OHOOH.is_organic=true;  		// Is the compound organic?
  C73K1OHOOH.compute_gamma_org=true;  	// Compute the activity coefficients of the organic phase for this compound?
  C73K1OHOOH.compute_gamma_aq=true;  	// Compute the activity coefficients of the aqueous phase for this compound?
  C73K1OHOOH.Koligo_org=0.0;      	//oligomeriation constant in the organic phase
  C73K1OHOOH.rho=1300.;  
  C73K1OHOOH.is_monomer=false;
  C73K1OHOOH.rion=false;
  C73K1OHOOH.KDiffusion_air=1.0e-5;
  //  C73K1OHOOH.accomodation_coefficient=alpha;
  C73K1OHOOH.viscosity=1.68e12;
  C73K1OHOOH.is_solid=false;
  C73K1OHOOH.is_generic=false;
  
  //Group: if no functionnal group in the species use the default species
  //for the computation of activity coefficients
  //
  double group_tmp_C73K1OHOOH [] = {0.,0.,0.,0., // group C
				  0.,0.,0.,0., //group C[OH]
				  0.,0.,0.,0., //group Calcohol
				  0.,0.,0.,0., //group Calcohol-tail
				  0.,0.,0.,0.,0., //group C=C
				  0.,0., //group aromatic carbon (AC)
				  0.,0.,0., // group //AC-C
				  1.,  //group OH
				  0., //group H2O
				  0., //group ACOH
				  0.,3., //group ketone
				  0.,   //group aldehyde  
				  0.,0., //group ester
				  0.,0.,0., //group ether 
				  0.,  //group acid
				  0.,   //group ACNO2
				  0.,0.,0., //group NO3
				  0.,0.,1., //group CO-OH
  			          0.,0.,0.,0.,0.,0.,0.,0.,0., //group CO-OC
			          0.,  //group PAN
			          0.,  //group CO-OOH
			          0.,  //group O=COC=O
			          0.,0.,0.}; //group CHxNO2
  
  size = sizeof(group_tmp_C73K1OHOOH)/sizeof(double);
  assert(size = 60);
	
  for(int i = 0; i < size; ++i)
	C73K1OHOOH.groups[i] = group_tmp_C73K1OHOOH[i];

  // Search the species name in the aerosol species list 
  // and add the species if its name matches with
  // the given list.
  add_species_ssh(surrogate, C73K1OHOOH, species_list_aer, molecular_weight_aer,
              accomodation_coefficient,diffusion_coef,saturation_vapor_pressure,enthalpy_vaporization, henry, t_ref, mass_density, species_part,nlayer,i_hydrophilic,
		  N_inert, N_inorganic);



/* ==== DD3000 ==== */ 

  species DD3000;
  DD3000.name="DD3000";
  DD3000.is_inorganic_precursor=false;
  DD3000.Psat_ref=1.288e-2; 		// Saturation vapor pressure at Tref (torr)
  DD3000.Tref=298;         		// Temperature of reference (K)
  DD3000.deltaH=74.475;     		// Enthalpy of vaporization (kJ/mol)
  DD3000.Henry=0.0;     		// Henry's law constant at Tref (M/atm)
  DD3000.aq_type="none"; 		// "none","diacid","monoacid" or "aldehyde"
  // DD3000.Kacidity1=XXacidityXX;   		// First acidity constant
  DD3000.hydrophilic=true;   	// Does the species condense on the aqueous phase?
  DD3000.hydrophobic=true;  		// Does the species condense on the organic phase?
  DD3000.nonvolatile=false; 		// Is the compound nonvolatile?
  DD3000.kp_from_experiment=false; 	// Use experimental partitioning constant at Tref?
  DD3000.is_organic=true;  		// Is the compound organic?
  DD3000.compute_gamma_org=true;  	// Compute the activity coefficients of the organic phase for this compound?
  DD3000.compute_gamma_aq=true;  	// Compute the activity coefficients of the aqueous phase for this compound?
  DD3000.Koligo_org=0.0;      	//oligomeriation constant in the organic phase
  DD3000.rho=1300.;  
  DD3000.is_monomer=false;
  DD3000.rion=false;
  DD3000.KDiffusion_air=1.0e-5;
  //  DD3000.accomodation_coefficient=alpha;
  DD3000.viscosity=1.68e12;
  DD3000.is_solid=false;
  DD3000.is_generic=false;
  
  //Group: if no functionnal group in the species use the default species
  //for the computation of activity coefficients
  //
  double group_tmp_DD3000 [] = {0.,0.,0.,0., // group C
				  0.,0.,1.,0., //group C[OH]
				  0.,0.,0.,0., //group Calcohol
				  0.,0.,0.,0., //group Calcohol-tail
				  0.,0.,0.,0.,0., //group C=C
				  0.,0., //group aromatic carbon (AC)
				  0.,0.,0., // group //AC-C
				  1.,  //group OH
				  0., //group H2O
				  0., //group ACOH
				  0.,0., //group ketone
				  2.,   //group aldehyde  
				  0.,0., //group ester
				  0.,0.,0., //group ether 
				  0.,  //group acid
				  0.,   //group ACNO2
				  0.,0.,0., //group NO3
				  0.,0.,0., //group CO-OH
  			          0.,0.,0.,0.,0.,0.,0.,0.,0., //group CO-OC
			          0.,  //group PAN
			          0.,  //group CO-OOH
			          0.,  //group O=COC=O
			          0.,0.,0.}; //group CHxNO2
  
  size = sizeof(group_tmp_DD3000)/sizeof(double);
  assert(size = 60);
	
  for(int i = 0; i < size; ++i)
	DD3000.groups[i] = group_tmp_DD3000[i];

  // Search the species name in the aerosol species list 
  // and add the species if its name matches with
  // the given list.
  add_species_ssh(surrogate, DD3000, species_list_aer, molecular_weight_aer,
              accomodation_coefficient,diffusion_coef,saturation_vapor_pressure,enthalpy_vaporization, henry, t_ref, mass_density, species_part,nlayer,i_hydrophilic,
		  N_inert, N_inorganic);



/* ==== DD3001 ==== */ 

  species DD3001;
  DD3001.name="DD3001";
  DD3001.is_inorganic_precursor=false;
  DD3001.Psat_ref=1.945; 		// Saturation vapor pressure at Tref (torr)
  DD3001.Tref=298;         		// Temperature of reference (K)
  DD3001.deltaH=51.14;     		// Enthalpy of vaporization (kJ/mol)
  DD3001.Henry=0.0;     		// Henry's law constant at Tref (M/atm)
  DD3001.aq_type="none"; 		// "none","diacid","monoacid" or "aldehyde"
  // DD3001.Kacidity1=XXacidityXX;   		// First acidity constant
  DD3001.hydrophilic=true;   	// Does the species condense on the aqueous phase?
  DD3001.hydrophobic=true;  		// Does the species condense on the organic phase?
  DD3001.nonvolatile=false; 		// Is the compound nonvolatile?
  DD3001.kp_from_experiment=false; 	// Use experimental partitioning constant at Tref?
  DD3001.is_organic=true;  		// Is the compound organic?
  DD3001.compute_gamma_org=true;  	// Compute the activity coefficients of the organic phase for this compound?
  DD3001.compute_gamma_aq=true;  	// Compute the activity coefficients of the aqueous phase for this compound?
  DD3001.Koligo_org=0.0;      	//oligomeriation constant in the organic phase
  DD3001.rho=1300.;  
  DD3001.is_monomer=false;
  DD3001.rion=false;
  DD3001.KDiffusion_air=1.0e-5;
  //  DD3001.accomodation_coefficient=alpha;
  DD3001.viscosity=1.68e12;
  DD3001.is_solid=false;
  DD3001.is_generic=false;
  
  //Group: if no functionnal group in the species use the default species
  //for the computation of activity coefficients
  //
  double group_tmp_DD3001 [] = {0.,0.,0.,0., // group C
				  0.,0.,0.,0., //group C[OH]
				  0.,0.,0.,0., //group Calcohol
				  0.,0.,0.,0., //group Calcohol-tail
				  0.,0.,0.,0.,0., //group C=C
				  0.,0., //group aromatic carbon (AC)
				  0.,0.,0., // group //AC-C
				  0.,  //group OH
				  0., //group H2O
				  0., //group ACOH
				  0.,1., //group ketone
				  2.,   //group aldehyde  
				  0.,0., //group ester
				  0.,0.,0., //group ether 
				  0.,  //group acid
				  0.,   //group ACNO2
				  0.,0.,0., //group NO3
				  0.,0.,0., //group CO-OH
  			          0.,0.,0.,0.,0.,0.,0.,0.,0., //group CO-OC
			          0.,  //group PAN
			          0.,  //group CO-OOH
			          0.,  //group O=COC=O
			          0.,0.,0.}; //group CHxNO2
  
  size = sizeof(group_tmp_DD3001)/sizeof(double);
  assert(size = 60);
	
  for(int i = 0; i < size; ++i)
	DD3001.groups[i] = group_tmp_DD3001[i];

  // Search the species name in the aerosol species list 
  // and add the species if its name matches with
  // the given list.
  add_species_ssh(surrogate, DD3001, species_list_aer, molecular_weight_aer,
              accomodation_coefficient,diffusion_coef,saturation_vapor_pressure,enthalpy_vaporization, henry, t_ref, mass_density, species_part,nlayer,i_hydrophilic,
		  N_inert, N_inorganic);



/* ==== DD5002 ==== */ 

  species DD5002;
  DD5002.name="DD5002";
  DD5002.is_inorganic_precursor=false;
  DD5002.Psat_ref=0.004795; 		// Saturation vapor pressure at Tref (torr)
  DD5002.Tref=298;         		// Temperature of reference (K)
  DD5002.deltaH=76.63;     		// Enthalpy of vaporization (kJ/mol)
  DD5002.Henry=0.0;     		// Henry's law constant at Tref (M/atm)
  DD5002.aq_type="none"; 		// "none","diacid","monoacid" or "aldehyde"
  // DD5002.Kacidity1=XXacidityXX;   		// First acidity constant
  DD5002.hydrophilic=true;   	// Does the species condense on the aqueous phase?
  DD5002.hydrophobic=true;  		// Does the species condense on the organic phase?
  DD5002.nonvolatile=false; 		// Is the compound nonvolatile?
  DD5002.kp_from_experiment=false; 	// Use experimental partitioning constant at Tref?
  DD5002.is_organic=true;  		// Is the compound organic?
  DD5002.compute_gamma_org=true;  	// Compute the activity coefficients of the organic phase for this compound?
  DD5002.compute_gamma_aq=true;  	// Compute the activity coefficients of the aqueous phase for this compound?
  DD5002.Koligo_org=0.0;      	//oligomeriation constant in the organic phase
  DD5002.rho=1300.;  
  DD5002.is_monomer=false;
  DD5002.rion=false;
  DD5002.KDiffusion_air=1.0e-5;
  //  DD5002.accomodation_coefficient=alpha;
  DD5002.viscosity=1.68e12;
  DD5002.is_solid=false;
  DD5002.is_generic=false;
  
  //Group: if no functionnal group in the species use the default species
  //for the computation of activity coefficients
  //
  double group_tmp_DD5002 [] = {0.,0.,0.,0., // group C
				  0.,0.,1.,0., //group C[OH]
				  0.,0.,0.,0., //group Calcohol
				  0.,0.,0.,0., //group Calcohol-tail
				  0.,0.,0.,0.,0., //group C=C
				  0.,0., //group aromatic carbon (AC)
				  0.,0.,0., // group //AC-C
				  1.,  //group OH
				  0., //group H2O
				  0., //group ACOH
				  1.,1., //group ketone
				  1.,   //group aldehyde  
				  0.,0., //group ester
				  0.,0.,0., //group ether 
				  0.,  //group acid
				  0.,   //group ACNO2
				  0.,0.,0., //group NO3
				  0.,0.,0., //group CO-OH
  			          0.,0.,0.,0.,0.,0.,0.,0.,0., //group CO-OC
			          0.,  //group PAN
			          0.,  //group CO-OOH
			          0.,  //group O=COC=O
			          0.,0.,0.}; //group CHxNO2
  
  size = sizeof(group_tmp_DD5002)/sizeof(double);
  assert(size = 60);
	
  for(int i = 0; i < size; ++i)
	DD5002.groups[i] = group_tmp_DD5002[i];

  // Search the species name in the aerosol species list 
  // and add the species if its name matches with
  // the given list.
  add_species_ssh(surrogate, DD5002, species_list_aer, molecular_weight_aer,
              accomodation_coefficient,diffusion_coef,saturation_vapor_pressure,enthalpy_vaporization, henry, t_ref, mass_density, species_part,nlayer,i_hydrophilic,
		  N_inert, N_inorganic);



/* ==== DK3000 ==== */ 

  species DK3000;
  DK3000.name="DK3000";
  DK3000.is_inorganic_precursor=false;
  DK3000.Psat_ref=1.131e+2; 		// Saturation vapor pressure at Tref (torr)
  DK3000.Tref=298;         		// Temperature of reference (K)
  DK3000.deltaH=35.19;     		// Enthalpy of vaporization (kJ/mol)
  DK3000.Henry=0.0;     		// Henry's law constant at Tref (M/atm)
  DK3000.aq_type="none"; 		// "none","diacid","monoacid" or "aldehyde"
  // DK3000.Kacidity1=XXacidityXX;   		// First acidity constant
  DK3000.hydrophilic=true;   	// Does the species condense on the aqueous phase?
  DK3000.hydrophobic=true;  		// Does the species condense on the organic phase?
  DK3000.nonvolatile=false; 		// Is the compound nonvolatile?
  DK3000.kp_from_experiment=false; 	// Use experimental partitioning constant at Tref?
  DK3000.is_organic=true;  		// Is the compound organic?
  DK3000.compute_gamma_org=true;  	// Compute the activity coefficients of the organic phase for this compound?
  DK3000.compute_gamma_aq=true;  	// Compute the activity coefficients of the aqueous phase for this compound?
  DK3000.Koligo_org=0.0;      	//oligomeriation constant in the organic phase
  DK3000.rho=1300.;  
  DK3000.is_monomer=false;
  DK3000.rion=false;
  DK3000.KDiffusion_air=1.0e-5;
  //  DK3000.accomodation_coefficient=alpha;
  DK3000.viscosity=1.68e12;
  DK3000.is_solid=false;
  DK3000.is_generic=false;
  
  //Group: if no functionnal group in the species use the default species
  //for the computation of activity coefficients
  //
  double group_tmp_DK3000 [] = {0.,0.,0.,0., // group C
				  0.,0.,0.,0., //group C[OH]
				  0.,0.,0.,0., //group Calcohol
				  0.,0.,0.,0., //group Calcohol-tail
				  0.,0.,0.,0.,0., //group C=C
				  0.,0., //group aromatic carbon (AC)
				  0.,0.,0., // group //AC-C
				  0.,  //group OH
				  0., //group H2O
				  0., //group ACOH
				  1.,0., //group ketone
				  1.,   //group aldehyde  
				  0.,0., //group ester
				  0.,0.,0., //group ether 
				  0.,  //group acid
				  0.,   //group ACNO2
				  0.,0.,0., //group NO3
				  0.,0.,0., //group CO-OH
  			          0.,0.,0.,0.,0.,0.,0.,0.,0., //group CO-OC
			          0.,  //group PAN
			          0.,  //group CO-OOH
			          0.,  //group O=COC=O
			          0.,0.,0.}; //group CHxNO2
  
  size = sizeof(group_tmp_DK3000)/sizeof(double);
  assert(size = 60);
	
  for(int i = 0; i < size; ++i)
	DK3000.groups[i] = group_tmp_DK3000[i];

  // Search the species name in the aerosol species list 
  // and add the species if its name matches with
  // the given list.
  add_species_ssh(surrogate, DK3000, species_list_aer, molecular_weight_aer,
              accomodation_coefficient,diffusion_coef,saturation_vapor_pressure,enthalpy_vaporization, henry, t_ref, mass_density, species_part,nlayer,i_hydrophilic,
		  N_inert, N_inorganic);

/* ==== IRDK3000 ==== */ 



  species IRDK3000;
  IRDK3000.name="IRDK3000";
  IRDK3000.is_inorganic_precursor=false;
  IRDK3000.Psat_ref=1.131e+2; 		// Saturation vapor pressure at Tref (torr)
  IRDK3000.Tref=298;         		// Temperature of reference (K)
  IRDK3000.deltaH=35.19;     		// Enthalpy of vaporization (kJ/mol)
  IRDK3000.Henry=0.0;     		// Henry's law constant at Tref (M/atm)
  IRDK3000.aq_type="none"; 		// "none","diacid","monoacid" or "aldehyde"
  // DK3000.Kacidity1=XXacidityXX;   		// First acidity constant
  IRDK3000.hydrophilic=true;   	// Does the species condense on the aqueous phase?
  IRDK3000.hydrophobic=false;  		// Does the species condense on the organic phase?
  IRDK3000.nonvolatile=true; 		// Is the compound nonvolatile?
  IRDK3000.kp_from_experiment=false; 	// Use experimental partitioning constant at Tref?
  IRDK3000.is_organic=true;  		// Is the compound organic?
  IRDK3000.compute_gamma_org=true;  	// Compute the activity coefficients of the organic phase for this compound?
  IRDK3000.compute_gamma_aq=true;  	// Compute the activity coefficients of the aqueous phase for this compound?
  IRDK3000.Koligo_org=0.0;      	//oligomeriation constant in the organic phase
  IRDK3000.rho=1300.;  
  IRDK3000.is_monomer=false;
  IRDK3000.rion=false;
  IRDK3000.KDiffusion_air=1.0e-5;
  //  DK3000.accomodation_coefficient=alpha;
  IRDK3000.viscosity=1.68e12;
  IRDK3000.is_solid=false;
  IRDK3000.is_generic=false;
  
  //Group: if no functionnal group in the species use the default species
  //for the computation of activity coefficients
  //
  double group_tmp_IRDK3000 [] = {0.,0.,0.,0., // group C
				  0.,0.,0.,0., //group C[OH]
				  0.,0.,0.,0., //group Calcohol
				  0.,0.,0.,0., //group Calcohol-tail
				  0.,0.,0.,0.,0., //group C=C
				  0.,0., //group aromatic carbon (AC)
				  0.,0.,0., // group //AC-C
				  0.,  //group OH
				  0., //group H2O
				  0., //group ACOH
				  1.,0., //group ketone
				  1.,   //group aldehyde  
				  0.,0., //group ester
				  0.,0.,0., //group ether 
				  0.,  //group acid
				  0.,   //group ACNO2
				  0.,0.,0., //group NO3
				  0.,0.,0., //group CO-OH
  			          0.,0.,0.,0.,0.,0.,0.,0.,0., //group CO-OC
			          0.,  //group PAN
			          0.,  //group CO-OOH
			          0.,  //group O=COC=O
			          0.,0.,0.}; //group CHxNO2
  
  size = sizeof(group_tmp_IRDK3000)/sizeof(double);
  assert(size = 60);
	
  for(int i = 0; i < size; ++i)
	IRDK3000.groups[i] = group_tmp_IRDK3000[i];

  // Search the species name in the aerosol species list 
  // and add the species if its name matches with
  // the given list.
  add_species_ssh(surrogate, IRDK3000, species_list_aer, molecular_weight_aer,
              accomodation_coefficient,diffusion_coef,saturation_vapor_pressure,enthalpy_vaporization, henry, t_ref, mass_density, species_part,nlayer,i_hydrophilic,
		  N_inert, N_inorganic);

/* ==== DK4000 ==== */ 

  species DK4000;
  DK4000.name="DK4000";
  DK4000.is_inorganic_precursor=false;
  DK4000.Psat_ref=1.112e-1; 		// Saturation vapor pressure at Tref (torr)
  DK4000.Tref=298;         		// Temperature of reference (K)
  DK4000.deltaH=63.91;     		// Enthalpy of vaporization (kJ/mol)
  DK4000.Henry=0.0;     		// Henry's law constant at Tref (M/atm)
  DK4000.aq_type="none"; 		// "none","diacid","monoacid" or "aldehyde"
  // DK4000.Kacidity1=XXacidityXX;   		// First acidity constant
  DK4000.hydrophilic=true;   	// Does the species condense on the aqueous phase?
  DK4000.hydrophobic=true;  		// Does the species condense on the organic phase?
  DK4000.nonvolatile=false; 		// Is the compound nonvolatile?
  DK4000.kp_from_experiment=false; 	// Use experimental partitioning constant at Tref?
  DK4000.is_organic=true;  		// Is the compound organic?
  DK4000.compute_gamma_org=true;  	// Compute the activity coefficients of the organic phase for this compound?
  DK4000.compute_gamma_aq=true;  	// Compute the activity coefficients of the aqueous phase for this compound?
  DK4000.Koligo_org=0.0;      	//oligomeriation constant in the organic phase
  DK4000.rho=1300.;  
  DK4000.is_monomer=false;
  DK4000.rion=false;
  DK4000.KDiffusion_air=1.0e-5;
  //  DK4000.accomodation_coefficient=alpha;
  DK4000.viscosity=1.68e12;
  DK4000.is_solid=false;
  DK4000.is_generic=false;
  
  //Group: if no functionnal group in the species use the default species
  //for the computation of activity coefficients
  //
  double group_tmp_DK4000 [] = {0.,0.,0.,0., // group C
				  0.,0.,1.,0., //group C[OH]
				  0.,0.,0.,0., //group Calcohol
				  0.,0.,0.,0., //group Calcohol-tail
				  0.,0.,0.,0.,0., //group C=C
				  0.,0., //group aromatic carbon (AC)
				  0.,0.,0., // group //AC-C
				  1.,  //group OH
				  0., //group H2O
				  0., //group ACOH
				  1.,0., //group ketone
				  1.,   //group aldehyde  
				  0.,0., //group ester
				  0.,0.,0., //group ether 
				  0.,  //group acid
				  0.,   //group ACNO2
				  0.,0.,0., //group NO3
				  0.,0.,0., //group CO-OH
  			          0.,0.,0.,0.,0.,0.,0.,0.,0., //group CO-OC
			          0.,  //group PAN
			          0.,  //group CO-OOH
			          0.,  //group O=COC=O
			          0.,0.,0.}; //group CHxNO2
  
  size = sizeof(group_tmp_DK4000)/sizeof(double);
  assert(size = 60);
	
  for(int i = 0; i < size; ++i)
	DK4000.groups[i] = group_tmp_DK4000[i];

  // Search the species name in the aerosol species list 
  // and add the species if its name matches with
  // the given list.
  add_species_ssh(surrogate, DK4000, species_list_aer, molecular_weight_aer,
              accomodation_coefficient,diffusion_coef,saturation_vapor_pressure,enthalpy_vaporization, henry, t_ref, mass_density, species_part,nlayer,i_hydrophilic,
		  N_inert, N_inorganic);



/* ==== DK4001 ==== */ 

  species DK4001;
  DK4001.name="DK4001";
  DK4001.is_inorganic_precursor=false;
  DK4001.Psat_ref=1.1814; 		// Saturation vapor pressure at Tref (torr)
  DK4001.Tref=298;         		// Temperature of reference (K)
  DK4001.deltaH=52.18;     		// Enthalpy of vaporization (kJ/mol)
  DK4001.Henry=0.0;     		// Henry's law constant at Tref (M/atm)
  DK4001.aq_type="none"; 		// "none","diacid","monoacid" or "aldehyde"
  // DK4001.Kacidity1=XXacidityXX;   		// First acidity constant
  DK4001.hydrophilic=true;   	// Does the species condense on the aqueous phase?
  DK4001.hydrophobic=true;  		// Does the species condense on the organic phase?
  DK4001.nonvolatile=false; 		// Is the compound nonvolatile?
  DK4001.kp_from_experiment=false; 	// Use experimental partitioning constant at Tref?
  DK4001.is_organic=true;  		// Is the compound organic?
  DK4001.compute_gamma_org=true;  	// Compute the activity coefficients of the organic phase for this compound?
  DK4001.compute_gamma_aq=true;  	// Compute the activity coefficients of the aqueous phase for this compound?
  DK4001.Koligo_org=0.0;      	//oligomeriation constant in the organic phase
  DK4001.rho=1300.;  
  DK4001.is_monomer=false;
  DK4001.rion=false;
  DK4001.KDiffusion_air=1.0e-5;
  //  DK4001.accomodation_coefficient=alpha;
  DK4001.viscosity=1.68e12;
  DK4001.is_solid=false;
  DK4001.is_generic=false;
  
  //Group: if no functionnal group in the species use the default species
  //for the computation of activity coefficients
  //
  double group_tmp_DK4001 [] = {0.,0.,0.,0., // group C
				  0.,0.,0.,0., //group C[OH]
				  0.,0.,0.,0., //group Calcohol
				  0.,0.,0.,0., //group Calcohol-tail
				  0.,0.,0.,0.,0., //group C=C
				  0.,0., //group aromatic carbon (AC)
				  0.,0.,0., // group //AC-C
				  0.,  //group OH
				  0., //group H2O
				  0., //group ACOH
				  1.,1., //group ketone
				  1.,   //group aldehyde  
				  0.,0., //group ester
				  0.,0.,0., //group ether 
				  0.,  //group acid
				  0.,   //group ACNO2
				  0.,0.,0., //group NO3
				  0.,0.,0., //group CO-OH
  			          0.,0.,0.,0.,0.,0.,0.,0.,0., //group CO-OC
			          0.,  //group PAN
			          0.,  //group CO-OOH
			          0.,  //group O=COC=O
			          0.,0.,0.}; //group CHxNO2
  
  size = sizeof(group_tmp_DK4001)/sizeof(double);
  assert(size = 60);
	
  for(int i = 0; i < size; ++i)
	DK4001.groups[i] = group_tmp_DK4001[i];

  // Search the species name in the aerosol species list 
  // and add the species if its name matches with
  // the given list.
  add_species_ssh(surrogate, DK4001, species_list_aer, molecular_weight_aer,
              accomodation_coefficient,diffusion_coef,saturation_vapor_pressure,enthalpy_vaporization, henry, t_ref, mass_density, species_part,nlayer,i_hydrophilic,
		  N_inert, N_inorganic);



/* ==== DK5000 ==== */ 

  species DK5000;
  DK5000.name="DK5000";
  DK5000.is_inorganic_precursor=false;
  DK5000.Psat_ref=0.00076; 		// Saturation vapor pressure at Tref (torr)
  DK5000.Tref=298;         		// Temperature of reference (K)
  DK5000.deltaH=85.34;     		// Enthalpy of vaporization (kJ/mol)
  DK5000.Henry=0.0;     		// Henry's law constant at Tref (M/atm)
  DK5000.aq_type="none"; 		// "none","diacid","monoacid" or "aldehyde"
  // DK5000.Kacidity1=XXacidityXX;   		// First acidity constant
  DK5000.hydrophilic=true;   	// Does the species condense on the aqueous phase?
  DK5000.hydrophobic=true;  		// Does the species condense on the organic phase?
  DK5000.nonvolatile=false; 		// Is the compound nonvolatile?
  DK5000.kp_from_experiment=false; 	// Use experimental partitioning constant at Tref?
  DK5000.is_organic=true;  		// Is the compound organic?
  DK5000.compute_gamma_org=true;  	// Compute the activity coefficients of the organic phase for this compound?
  DK5000.compute_gamma_aq=true;  	// Compute the activity coefficients of the aqueous phase for this compound?
  DK5000.Koligo_org=0.0;      	//oligomeriation constant in the organic phase
  DK5000.rho=1300.;  
  DK5000.is_monomer=false;
  DK5000.rion=false;
  DK5000.KDiffusion_air=1.0e-5;
  //  DK5000.accomodation_coefficient=alpha;
  DK5000.viscosity=1.68e12;
  DK5000.is_solid=false;
  DK5000.is_generic=false;
  
  //Group: if no functionnal group in the species use the default species
  //for the computation of activity coefficients
  //
  double group_tmp_DK5000 [] = {0.,0.,0.,0., // group C
				  0.,0.,0.,0., //group C[OH]
				  0.,0.,0.,0., //group Calcohol
				  0.,0.,0.,0., //group Calcohol-tail
				  0.,0.,0.,0.,1., //group C=C
				  0.,0., //group aromatic carbon (AC)
				  0.,0.,0., // group //AC-C
				  2.,  //group OH
				  0., //group H2O
				  0., //group ACOH
				  1.,0., //group ketone
				  1.,   //group aldehyde  
				  0.,0., //group ester
				  0.,0.,0., //group ether 
				  0.,  //group acid
				  0.,   //group ACNO2
				  0.,0.,0., //group NO3
				  0.,0.,0., //group CO-OH
  			          0.,0.,0.,0.,0.,0.,0.,0.,0., //group CO-OC
			          0.,  //group PAN
			          0.,  //group CO-OOH
			          0.,  //group O=COC=O
			          0.,0.,0.}; //group CHxNO2
  
  size = sizeof(group_tmp_DK5000)/sizeof(double);
  assert(size = 60);
	
  for(int i = 0; i < size; ++i)
	DK5000.groups[i] = group_tmp_DK5000[i];

  // Search the species name in the aerosol species list 
  // and add the species if its name matches with
  // the given list.
  add_species_ssh(surrogate, DK5000, species_list_aer, molecular_weight_aer,
              accomodation_coefficient,diffusion_coef,saturation_vapor_pressure,enthalpy_vaporization, henry, t_ref, mass_density, species_part,nlayer,i_hydrophilic,
		  N_inert, N_inorganic);



/* ==== ED4001 ==== */ 

  species ED4001;
  ED4001.name="ED4001";
  ED4001.is_inorganic_precursor=false;
  ED4001.Psat_ref=1.205; 		// Saturation vapor pressure at Tref (torr)
  ED4001.Tref=298;         		// Temperature of reference (K)
  ED4001.deltaH=50.0;     		// Enthalpy of vaporization (kJ/mol)
  ED4001.Henry=0.0;     		// Henry's law constant at Tref (M/atm)
  ED4001.aq_type="none"; 		// "none","diacid","monoacid" or "aldehyde"
  // ED4001.Kacidity1=XXacidityXX;   		// First acidity constant
  ED4001.hydrophilic=true;   	// Does the species condense on the aqueous phase?
  ED4001.hydrophobic=true;  		// Does the species condense on the organic phase?
  ED4001.nonvolatile=false; 		// Is the compound nonvolatile?
  ED4001.kp_from_experiment=false; 	// Use experimental partitioning constant at Tref?
  ED4001.is_organic=true;  		// Is the compound organic?
  ED4001.compute_gamma_org=true;  	// Compute the activity coefficients of the organic phase for this compound?
  ED4001.compute_gamma_aq=true;  	// Compute the activity coefficients of the aqueous phase for this compound?
  ED4001.Koligo_org=0.0;      	//oligomeriation constant in the organic phase
  ED4001.rho=1300.;  
  ED4001.is_monomer=false;
  ED4001.rion=false;
  ED4001.KDiffusion_air=1.0e-5;
  //  ED4001.accomodation_coefficient=alpha;
  ED4001.viscosity=1.68e12;
  ED4001.is_solid=false;
  ED4001.is_generic=false;
  
  //Group: if no functionnal group in the species use the default species
  //for the computation of activity coefficients
  //
  double group_tmp_ED4001 [] = {0.,0.,0.,0., // group C
				  0.,0.,0.,0., //group C[OH]
				  0.,0.,0.,0., //group Calcohol
				  0.,0.,0.,0., //group Calcohol-tail
				  0.,0.,0.,0.,0., //group C=C
				  0.,0., //group aromatic carbon (AC)
				  0.,0.,0., // group //AC-C
				  0.,  //group OH
				  0., //group H2O
				  0., //group ACOH
				  0.,0., //group ketone
				  2.,   //group aldehyde  
				  0.,1., //group ester
				  0.,0.,0., //group ether 
				  0.,  //group acid
				  0.,   //group ACNO2
				  0.,0.,0., //group NO3
				  0.,0.,0., //group CO-OH
  			          0.,0.,0.,0.,0.,0.,0.,0.,0., //group CO-OC
			          0.,  //group PAN
			          0.,  //group CO-OOH
			          0.,  //group O=COC=O
			          0.,0.,0.}; //group CHxNO2
  
  size = sizeof(group_tmp_ED4001)/sizeof(double);
  assert(size = 60);
	
  for(int i = 0; i < size; ++i)
	ED4001.groups[i] = group_tmp_ED4001[i];

  // Search the species name in the aerosol species list 
  // and add the species if its name matches with
  // the given list.
  add_species_ssh(surrogate, ED4001, species_list_aer, molecular_weight_aer,
              accomodation_coefficient,diffusion_coef,saturation_vapor_pressure,enthalpy_vaporization, henry, t_ref, mass_density, species_part,nlayer,i_hydrophilic,
		  N_inert, N_inorganic);



/* ==== ED5000 ==== */ 

  species ED5000;
  ED5000.name="ED5000";
  ED5000.is_inorganic_precursor=false;
  ED5000.Psat_ref=3.809e-1; 		// Saturation vapor pressure at Tref (torr)
  ED5000.Tref=298;         		// Temperature of reference (K)
  ED5000.deltaH=55.0;     		// Enthalpy of vaporization (kJ/mol)
  ED5000.Henry=0.0;     		// Henry's law constant at Tref (M/atm)
  ED5000.aq_type="none"; 		// "none","diacid","monoacid" or "aldehyde"
  // ED5000.Kacidity1=XXacidityXX;   		// First acidity constant
  ED5000.hydrophilic=true;   	// Does the species condense on the aqueous phase?
  ED5000.hydrophobic=true;  		// Does the species condense on the organic phase?
  ED5000.nonvolatile=false; 		// Is the compound nonvolatile?
  ED5000.kp_from_experiment=false; 	// Use experimental partitioning constant at Tref?
  ED5000.is_organic=true;  		// Is the compound organic?
  ED5000.compute_gamma_org=true;  	// Compute the activity coefficients of the organic phase for this compound?
  ED5000.compute_gamma_aq=true;  	// Compute the activity coefficients of the aqueous phase for this compound?
  ED5000.Koligo_org=0.0;      	//oligomeriation constant in the organic phase
  ED5000.rho=1300.;  
  ED5000.is_monomer=false;
  ED5000.rion=false;
  ED5000.KDiffusion_air=1.0e-5;
  //  ED5000.accomodation_coefficient=alpha;
  ED5000.viscosity=1.68e12;
  ED5000.is_solid=false;
  ED5000.is_generic=false;
  
  //Group: if no functionnal group in the species use the default species
  //for the computation of activity coefficients
  //
  double group_tmp_ED5000 [] = {0.,0.,0.,0., // group C
				  0.,0.,0.,0., //group C[OH]
				  0.,0.,0.,0., //group Calcohol
				  0.,0.,0.,0., //group Calcohol-tail
				  0.,0.,0.,0.,0., //group C=C
				  0.,0., //group aromatic carbon (AC)
				  0.,0.,0., // group //AC-C
				  0.,  //group OH
				  0., //group H2O
				  0., //group ACOH
				  0.,1., //group ketone
				  1.,   //group aldehyde  
				  1.,0., //group ester
				  0.,0.,0., //group ether 
				  0.,  //group acid
				  0.,   //group ACNO2
				  0.,0.,0., //group NO3
				  0.,0.,0., //group CO-OH
  			          0.,0.,0.,0.,0.,0.,0.,0.,0., //group CO-OC
			          0.,  //group PAN
			          0.,  //group CO-OOH
			          0.,  //group O=COC=O
			          0.,0.,0.}; //group CHxNO2
  
  size = sizeof(group_tmp_ED5000)/sizeof(double);
  assert(size = 60);
	
  for(int i = 0; i < size; ++i)
	ED5000.groups[i] = group_tmp_ED5000[i];

  // Search the species name in the aerosol species list 
  // and add the species if its name matches with
  // the given list.
  add_species_ssh(surrogate, ED5000, species_list_aer, molecular_weight_aer,
              accomodation_coefficient,diffusion_coef,saturation_vapor_pressure,enthalpy_vaporization, henry, t_ref, mass_density, species_part,nlayer,i_hydrophilic,
		  N_inert, N_inorganic);



/* ==== ED5000OOH ==== */ 

  species ED5000OOH;
  ED5000OOH.name="ED5000OOH";
  ED5000OOH.is_inorganic_precursor=false;
  ED5000OOH.Psat_ref=6.037e-3; 		// Saturation vapor pressure at Tref (torr)
  ED5000OOH.Tref=298;         		// Temperature of reference (K)
  ED5000OOH.deltaH=75.0;     		// Enthalpy of vaporization (kJ/mol)
  ED5000OOH.Henry=0.0;     		// Henry's law constant at Tref (M/atm)
  ED5000OOH.aq_type="none"; 		// "none","diacid","monoacid" or "aldehyde"
  // ED5000OOH.Kacidity1=XXacidityXX;   		// First acidity constant
  ED5000OOH.hydrophilic=true;   	// Does the species condense on the aqueous phase?
  ED5000OOH.hydrophobic=true;  		// Does the species condense on the organic phase?
  ED5000OOH.nonvolatile=false; 		// Is the compound nonvolatile?
  ED5000OOH.kp_from_experiment=false; 	// Use experimental partitioning constant at Tref?
  ED5000OOH.is_organic=true;  		// Is the compound organic?
  ED5000OOH.compute_gamma_org=true;  	// Compute the activity coefficients of the organic phase for this compound?
  ED5000OOH.compute_gamma_aq=true;  	// Compute the activity coefficients of the aqueous phase for this compound?
  ED5000OOH.Koligo_org=0.0;      	//oligomeriation constant in the organic phase
  ED5000OOH.rho=1300.;  
  ED5000OOH.is_monomer=false;
  ED5000OOH.rion=false;
  ED5000OOH.KDiffusion_air=1.0e-5;
  //  ED5000OOH.accomodation_coefficient=alpha;
  ED5000OOH.viscosity=1.68e12;
  ED5000OOH.is_solid=false;
  ED5000OOH.is_generic=false;
  
  //Group: if no functionnal group in the species use the default species
  //for the computation of activity coefficients
  //
  double group_tmp_ED5000OOH [] = {0.,0.,0.,0., // group C
				  0.,0.,0.,0., //group C[OH]
				  0.,0.,0.,0., //group Calcohol
				  0.,0.,0.,0., //group Calcohol-tail
				  0.,0.,0.,0.,0., //group C=C
				  0.,0., //group aromatic carbon (AC)
				  0.,0.,0., // group //AC-C
				  0.,  //group OH
				  0., //group H2O
				  0., //group ACOH
				  0.,1., //group ketone
				  0.,   //group aldehyde  
				  1.,0., //group ester
				  0.,0.,0., //group ether 
				  0.,  //group acid
				  0.,   //group ACNO2
				  0.,0.,0., //group NO3
				  0.,0.,0., //group CO-OH
  			          0.,0.,0.,0.,0.,0.,0.,0.,0., //group CO-OC
			          0.,  //group PAN
			          1.,  //group CO-OOH
			          0.,  //group O=COC=O
			          0.,0.,0.}; //group CHxNO2
  
  size = sizeof(group_tmp_ED5000OOH)/sizeof(double);
  assert(size = 60);
	
  for(int i = 0; i < size; ++i)
	ED5000OOH.groups[i] = group_tmp_ED5000OOH[i];

  // Search the species name in the aerosol species list 
  // and add the species if its name matches with
  // the given list.
  add_species_ssh(surrogate, ED5000OOH, species_list_aer, molecular_weight_aer,
              accomodation_coefficient,diffusion_coef,saturation_vapor_pressure,enthalpy_vaporization, henry, t_ref, mass_density, species_part,nlayer,i_hydrophilic,
		  N_inert, N_inorganic);



/* ==== ED5002OOH ==== */ 

  species ED5002OOH;
  ED5002OOH.name="ED5002OOH";
  ED5002OOH.is_inorganic_precursor=false;
  ED5002OOH.Psat_ref=1.205e-2; 		// Saturation vapor pressure at Tref (torr)
  ED5002OOH.Tref=298;         		// Temperature of reference (K)
  ED5002OOH.deltaH=70.0;     		// Enthalpy of vaporization (kJ/mol)
  ED5002OOH.Henry=0.0;     		// Henry's law constant at Tref (M/atm)
  ED5002OOH.aq_type="none"; 		// "none","diacid","monoacid" or "aldehyde"
  // ED5002OOH.Kacidity1=XXacidityXX;   		// First acidity constant
  ED5002OOH.hydrophilic=true;   	// Does the species condense on the aqueous phase?
  ED5002OOH.hydrophobic=true;  		// Does the species condense on the organic phase?
  ED5002OOH.nonvolatile=false; 		// Is the compound nonvolatile?
  ED5002OOH.kp_from_experiment=false; 	// Use experimental partitioning constant at Tref?
  ED5002OOH.is_organic=true;  		// Is the compound organic?
  ED5002OOH.compute_gamma_org=true;  	// Compute the activity coefficients of the organic phase for this compound?
  ED5002OOH.compute_gamma_aq=true;  	// Compute the activity coefficients of the aqueous phase for this compound?
  ED5002OOH.Koligo_org=0.0;      	//oligomeriation constant in the organic phase
  ED5002OOH.rho=1300.;  
  ED5002OOH.is_monomer=false;
  ED5002OOH.rion=false;
  ED5002OOH.KDiffusion_air=1.0e-5;
  //  ED5002OOH.accomodation_coefficient=alpha;
  ED5002OOH.viscosity=1.68e12;
  ED5002OOH.is_solid=false;
  ED5002OOH.is_generic=false;
  
  //Group: if no functionnal group in the species use the default species
  //for the computation of activity coefficients
  //
  double group_tmp_ED5002OOH [] = {0.,0.,0.,0., // group C
				  0.,0.,0.,0., //group C[OH]
				  0.,0.,0.,0., //group Calcohol
				  0.,0.,0.,0., //group Calcohol-tail
				  0.,0.,0.,0.,0., //group C=C
				  0.,0., //group aromatic carbon (AC)
				  0.,0.,0., // group //AC-C
				  0.,  //group OH
				  0., //group H2O
				  0., //group ACOH
				  0.,1., //group ketone
				  0.,   //group aldehyde  
				  1.,0., //group ester
				  0.,0.,0., //group ether 
				  0.,  //group acid
				  0.,   //group ACNO2
				  0.,0.,0., //group NO3
				  1.,0.,0., //group CO-OH
  			          0.,0.,0.,0.,0.,0.,0.,0.,0., //group CO-OC
			          0.,  //group PAN
			          0.,  //group CO-OOH
			          0.,  //group O=COC=O
			          0.,0.,0.}; //group CHxNO2
  
  size = sizeof(group_tmp_ED5002OOH)/sizeof(double);
  assert(size = 60);
	
  for(int i = 0; i < size; ++i)
	ED5002OOH.groups[i] = group_tmp_ED5002OOH[i];

  // Search the species name in the aerosol species list 
  // and add the species if its name matches with
  // the given list.
  add_species_ssh(surrogate, ED5002OOH, species_list_aer, molecular_weight_aer,
              accomodation_coefficient,diffusion_coef,saturation_vapor_pressure,enthalpy_vaporization, henry, t_ref, mass_density, species_part,nlayer,i_hydrophilic,
		  N_inert, N_inorganic);



/* ==== FUROH ==== */ 

  species FUROH;
  FUROH.name="FUROH";
  FUROH.is_inorganic_precursor=false;
  FUROH.Psat_ref=9.567e-2; 		// Saturation vapor pressure at Tref (torr)
  FUROH.Tref=298;         		// Temperature of reference (K)
  FUROH.deltaH=64.9;     		// Enthalpy of vaporization (kJ/mol)
  FUROH.Henry=0.0;     		// Henry's law constant at Tref (M/atm)
  FUROH.aq_type="none"; 		// "none","diacid","monoacid" or "aldehyde"
  // FUROH.Kacidity1=XXacidityXX;   		// First acidity constant
  FUROH.hydrophilic=true;   	// Does the species condense on the aqueous phase?
  FUROH.hydrophobic=true;  		// Does the species condense on the organic phase?
  FUROH.nonvolatile=false; 		// Is the compound nonvolatile?
  FUROH.kp_from_experiment=false; 	// Use experimental partitioning constant at Tref?
  FUROH.is_organic=true;  		// Is the compound organic?
  FUROH.compute_gamma_org=true;  	// Compute the activity coefficients of the organic phase for this compound?
  FUROH.compute_gamma_aq=true;  	// Compute the activity coefficients of the aqueous phase for this compound?
  FUROH.Koligo_org=0.0;      	//oligomeriation constant in the organic phase
  FUROH.rho=1300.;  
  FUROH.is_monomer=false;
  FUROH.rion=false;
  FUROH.KDiffusion_air=1.0e-5;
  //  FUROH.accomodation_coefficient=alpha;
  FUROH.viscosity=1.68e12;
  FUROH.is_solid=false;
  FUROH.is_generic=false;
  
  //Group: if no functionnal group in the species use the default species
  //for the computation of activity coefficients
  //
  double group_tmp_FUROH [] = {0.,0.,0.,0., // group C
				  0.,0.,0.,0., //group C[OH]
				  0.,0.,0.,0., //group Calcohol
				  0.,0.,0.,0., //group Calcohol-tail
				  0.,0.,0.,2.,0., //group C=C
				  0.,0., //group aromatic carbon (AC)
				  0.,0.,0., // group //AC-C
				  1.,  //group OH
				  0., //group H2O
				  0., //group ACOH
				  0.,0., //group ketone
				  1.,   //group aldehyde  
				  0.,0., //group ester
				  0.,0.,1., //group ether 
				  0.,  //group acid
				  0.,   //group ACNO2
				  0.,0.,0., //group NO3
				  0.,0.,0., //group CO-OH
  			          0.,0.,0.,0.,0.,0.,0.,0.,0., //group CO-OC
			          0.,  //group PAN
			          0.,  //group CO-OOH
			          0.,  //group O=COC=O
			          0.,0.,0.}; //group CHxNO2
  
  size = sizeof(group_tmp_FUROH)/sizeof(double);
  assert(size = 60);
	
  for(int i = 0; i < size; ++i)
	FUROH.groups[i] = group_tmp_FUROH[i];

  // Search the species name in the aerosol species list 
  // and add the species if its name matches with
  // the given list.
  add_species_ssh(surrogate, FUROH, species_list_aer, molecular_weight_aer,
              accomodation_coefficient,diffusion_coef,saturation_vapor_pressure,enthalpy_vaporization, henry, t_ref, mass_density, species_part,nlayer,i_hydrophilic,
		  N_inert, N_inorganic);



/* ==== FURON ==== */ 

  species FURON;
  FURON.name="FURON";
  FURON.is_inorganic_precursor=false;
  FURON.Psat_ref=9.180e-1; 		// Saturation vapor pressure at Tref (torr)
  FURON.Tref=298;         		// Temperature of reference (K)
  FURON.deltaH=51.15;     		// Enthalpy of vaporization (kJ/mol)
  FURON.Henry=0.0;     		// Henry's law constant at Tref (M/atm)
  FURON.aq_type="none"; 		// "none","diacid","monoacid" or "aldehyde"
  // FURON.Kacidity1=XXacidityXX;   		// First acidity constant
  FURON.hydrophilic=true;   	// Does the species condense on the aqueous phase?
  FURON.hydrophobic=true;  		// Does the species condense on the organic phase?
  FURON.nonvolatile=false; 		// Is the compound nonvolatile?
  FURON.kp_from_experiment=false; 	// Use experimental partitioning constant at Tref?
  FURON.is_organic=true;  		// Is the compound organic?
  FURON.compute_gamma_org=true;  	// Compute the activity coefficients of the organic phase for this compound?
  FURON.compute_gamma_aq=true;  	// Compute the activity coefficients of the aqueous phase for this compound?
  FURON.Koligo_org=0.0;      	//oligomeriation constant in the organic phase
  FURON.rho=1300.;  
  FURON.is_monomer=false;
  FURON.rion=false;
  FURON.KDiffusion_air=1.0e-5;
  //  FURON.accomodation_coefficient=alpha;
  FURON.viscosity=1.68e12;
  FURON.is_solid=false;
  FURON.is_generic=false;
  
  //Group: if no functionnal group in the species use the default species
  //for the computation of activity coefficients
  //
  double group_tmp_FURON [] = {0.,0.,0.,0., // group C
				  0.,0.,0.,0., //group C[OH]
				  0.,0.,0.,0., //group Calcohol
				  0.,0.,0.,0., //group Calcohol-tail
				  0.,0.,0.,1.,0., //group C=C
				  0.,0., //group aromatic carbon (AC)
				  0.,0.,0., // group //AC-C
				  0.,  //group OH
				  0., //group H2O
				  0., //group ACOH
				  0.,0., //group ketone
				  0.,   //group aldehyde  
				  0.,1., //group ester
				  0.,0.,0., //group ether 
				  0.,  //group acid
				  0.,   //group ACNO2
				  0.,0.,0., //group NO3
				  0.,0.,0., //group CO-OH
  			          0.,0.,0.,0.,0.,0.,0.,0.,0., //group CO-OC
			          0.,  //group PAN
			          0.,  //group CO-OOH
			          0.,  //group O=COC=O
			          0.,0.,0.}; //group CHxNO2
  
  size = sizeof(group_tmp_FURON)/sizeof(double);
  assert(size = 60);
	
  for(int i = 0; i < size; ++i)
	FURON.groups[i] = group_tmp_FURON[i];

  // Search the species name in the aerosol species list 
  // and add the species if its name matches with
  // the given list.
  add_species_ssh(surrogate, FURON, species_list_aer, molecular_weight_aer,
              accomodation_coefficient,diffusion_coef,saturation_vapor_pressure,enthalpy_vaporization, henry, t_ref, mass_density, species_part,nlayer,i_hydrophilic,
		  N_inert, N_inorganic);



/* ==== FURR3 ==== */ 

  species FURR3;
  FURR3.name="FURR3";
  FURR3.is_inorganic_precursor=false;
  FURR3.Psat_ref=1.516e1; 		// Saturation vapor pressure at Tref (torr)
  FURR3.Tref=298;         		// Temperature of reference (K)
  FURR3.deltaH=44.86;     		// Enthalpy of vaporization (kJ/mol)
  FURR3.Henry=0.0;     		// Henry's law constant at Tref (M/atm)
  FURR3.aq_type="none"; 		// "none","diacid","monoacid" or "aldehyde"
  // FURR3.Kacidity1=XXacidityXX;   		// First acidity constant
  FURR3.hydrophilic=true;   	// Does the species condense on the aqueous phase?
  FURR3.hydrophobic=true;  		// Does the species condense on the organic phase?
  FURR3.nonvolatile=false; 		// Is the compound nonvolatile?
  FURR3.kp_from_experiment=false; 	// Use experimental partitioning constant at Tref?
  FURR3.is_organic=true;  		// Is the compound organic?
  FURR3.compute_gamma_org=true;  	// Compute the activity coefficients of the organic phase for this compound?
  FURR3.compute_gamma_aq=true;  	// Compute the activity coefficients of the aqueous phase for this compound?
  FURR3.Koligo_org=0.0;      	//oligomeriation constant in the organic phase
  FURR3.rho=1300.;  
  FURR3.is_monomer=false;
  FURR3.rion=false;
  FURR3.KDiffusion_air=1.0e-5;
  //  FURR3.accomodation_coefficient=alpha;
  FURR3.viscosity=1.68e12;
  FURR3.is_solid=false;
  FURR3.is_generic=false;
  
  //Group: if no functionnal group in the species use the default species
  //for the computation of activity coefficients
  //
  double group_tmp_FURR3 [] = {0.,0.,0.,0., // group C
				  0.,0.,0.,0., //group C[OH]
				  0.,0.,0.,0., //group Calcohol
				  0.,0.,0.,0., //group Calcohol-tail
				  0.,1.,0.,1.,0., //group C=C
				  0.,0., //group aromatic carbon (AC)
				  0.,0.,0., // group //AC-C
				  0.,  //group OH
				  0., //group H2O
				  0., //group ACOH
				  0.,0., //group ketone
				  1.,   //group aldehyde  
				  0.,0., //group ester
				  0.,1.,0., //group ether 
				  0.,  //group acid
				  0.,   //group ACNO2
				  0.,0.,0., //group NO3
				  0.,0.,0., //group CO-OH
  			          0.,0.,0.,0.,0.,0.,0.,0.,0., //group CO-OC
			          0.,  //group PAN
			          0.,  //group CO-OOH
			          0.,  //group O=COC=O
			          0.,0.,0.}; //group CHxNO2
  
  size = sizeof(group_tmp_FURR3)/sizeof(double);
  assert(size = 60);
	
  for(int i = 0; i < size; ++i)
	FURR3.groups[i] = group_tmp_FURR3[i];

  // Search the species name in the aerosol species list 
  // and add the species if its name matches with
  // the given list.
  add_species_ssh(surrogate, FURR3, species_list_aer, molecular_weight_aer,
              accomodation_coefficient,diffusion_coef,saturation_vapor_pressure,enthalpy_vaporization, henry, t_ref, mass_density, species_part,nlayer,i_hydrophilic,
		  N_inert, N_inorganic);



/* ==== FURR5 ==== */ 

  species FURR5;
  FURR5.name="FURR5";
  FURR5.is_inorganic_precursor=false;
  FURR5.Psat_ref=8.468; 		// Saturation vapor pressure at Tref (torr)
  FURR5.Tref=298;         		// Temperature of reference (K)
  FURR5.deltaH=44.86;     		// Enthalpy of vaporization (kJ/mol)
  FURR5.Henry=0.0;     		// Henry's law constant at Tref (M/atm)
  FURR5.aq_type="none"; 		// "none","diacid","monoacid" or "aldehyde"
  // FURR5.Kacidity1=XXacidityXX;   		// First acidity constant
  FURR5.hydrophilic=true;   	// Does the species condense on the aqueous phase?
  FURR5.hydrophobic=true;  		// Does the species condense on the organic phase?
  FURR5.nonvolatile=false; 		// Is the compound nonvolatile?
  FURR5.kp_from_experiment=false; 	// Use experimental partitioning constant at Tref?
  FURR5.is_organic=true;  		// Is the compound organic?
  FURR5.compute_gamma_org=true;  	// Compute the activity coefficients of the organic phase for this compound?
  FURR5.compute_gamma_aq=true;  	// Compute the activity coefficients of the aqueous phase for this compound?
  FURR5.Koligo_org=0.0;      	//oligomeriation constant in the organic phase
  FURR5.rho=1300.;  
  FURR5.is_monomer=false;
  FURR5.rion=false;
  FURR5.KDiffusion_air=1.0e-5;
  //  FURR5.accomodation_coefficient=alpha;
  FURR5.viscosity=1.68e12;
  FURR5.is_solid=false;
  FURR5.is_generic=false;
  
  //Group: if no functionnal group in the species use the default species
  //for the computation of activity coefficients
  //
  double group_tmp_FURR5 [] = {0.,0.,0.,0., // group C
				  0.,0.,0.,0., //group C[OH]
				  0.,0.,0.,0., //group Calcohol
				  0.,0.,0.,0., //group Calcohol-tail
				  0.,1.,0.,1.,0., //group C=C
				  0.,0., //group aromatic carbon (AC)
				  0.,0.,0., // group //AC-C
				  0.,  //group OH
				  0., //group H2O
				  0., //group ACOH
				  0.,0., //group ketone
				  1.,   //group aldehyde  
				  0.,0., //group ester
				  0.,1.,0., //group ether 
				  0.,  //group acid
				  0.,   //group ACNO2
				  0.,0.,0., //group NO3
				  0.,0.,0., //group CO-OH
  			          0.,0.,0.,0.,0.,0.,0.,0.,0., //group CO-OC
			          0.,  //group PAN
			          0.,  //group CO-OOH
			          0.,  //group O=COC=O
			          0.,0.,0.}; //group CHxNO2
  
  size = sizeof(group_tmp_FURR5)/sizeof(double);
  assert(size = 60);
	
  for(int i = 0; i < size; ++i)
	FURR5.groups[i] = group_tmp_FURR5[i];

  // Search the species name in the aerosol species list 
  // and add the species if its name matches with
  // the given list.
  add_species_ssh(surrogate, FURR5, species_list_aer, molecular_weight_aer,
              accomodation_coefficient,diffusion_coef,saturation_vapor_pressure,enthalpy_vaporization, henry, t_ref, mass_density, species_part,nlayer,i_hydrophilic,
		  N_inert, N_inorganic);



/* ==== FURR6 ==== */ 

  species FURR6;
  FURR6.name="FURR6";
  FURR6.is_inorganic_precursor=false;
  FURR6.Psat_ref=1.909e-3; 		// Saturation vapor pressure at Tref (torr)
  FURR6.Tref=298;         		// Temperature of reference (K)
  FURR6.deltaH=80.6;     		// Enthalpy of vaporization (kJ/mol)
  FURR6.Henry=0.0;     		// Henry's law constant at Tref (M/atm)
  FURR6.aq_type="none"; 		// "none","diacid","monoacid" or "aldehyde"
  // FURR6.Kacidity1=XXacidityXX;   		// First acidity constant
  FURR6.hydrophilic=true;   	// Does the species condense on the aqueous phase?
  FURR6.hydrophobic=true;  		// Does the species condense on the organic phase?
  FURR6.nonvolatile=false; 		// Is the compound nonvolatile?
  FURR6.kp_from_experiment=false; 	// Use experimental partitioning constant at Tref?
  FURR6.is_organic=true;  		// Is the compound organic?
  FURR6.compute_gamma_org=true;  	// Compute the activity coefficients of the organic phase for this compound?
  FURR6.compute_gamma_aq=true;  	// Compute the activity coefficients of the aqueous phase for this compound?
  FURR6.Koligo_org=0.0;      	//oligomeriation constant in the organic phase
  FURR6.rho=1300.;  
  FURR6.is_monomer=false;
  FURR6.rion=false;
  FURR6.KDiffusion_air=1.0e-5;
  //  FURR6.accomodation_coefficient=alpha;
  FURR6.viscosity=1.68e12;
  FURR6.is_solid=false;
  FURR6.is_generic=false;
  
  //Group: if no functionnal group in the species use the default species
  //for the computation of activity coefficients
  //
  double group_tmp_FURR6 [] = {0.,0.,0.,0., // group C
				  0.,0.,0.,0., //group C[OH]
				  0.,0.,0.,0., //group Calcohol
				  0.,0.,0.,0., //group Calcohol-tail
				  0.,0.,0.,1.,0., //group C=C
				  0.,0., //group aromatic carbon (AC)
				  0.,0.,0., // group //AC-C
				  0.,  //group OH
				  0., //group H2O
				  0., //group ACOH
				  0.,1., //group ketone
				  0.,   //group aldehyde  
				  0.,1., //group ester
				  0.,0.,0., //group ether 
				  0.,  //group acid
				  0.,   //group ACNO2
				  0.,0.,0., //group NO3
				  0.,0.,0., //group CO-OH
  			          0.,0.,0.,0.,0.,0.,0.,0.,0., //group CO-OC
			          0.,  //group PAN
			          0.,  //group CO-OOH
			          0.,  //group O=COC=O
			          0.,0.,0.}; //group CHxNO2
  
  size = sizeof(group_tmp_FURR6)/sizeof(double);
  assert(size = 60);
	
  for(int i = 0; i < size; ++i)
	FURR6.groups[i] = group_tmp_FURR6[i];

  // Search the species name in the aerosol species list 
  // and add the species if its name matches with
  // the given list.
  add_species_ssh(surrogate, FURR6, species_list_aer, molecular_weight_aer,
              accomodation_coefficient,diffusion_coef,saturation_vapor_pressure,enthalpy_vaporization, henry, t_ref, mass_density, species_part,nlayer,i_hydrophilic,
		  N_inert, N_inorganic);



/* ==== FURR6OHOOH ==== */ 

  species FURR6OHOOH;
  FURR6OHOOH.name="FURR6OHOOH";
  FURR6OHOOH.is_inorganic_precursor=false;
  FURR6OHOOH.Psat_ref=4.795e-9; 		// Saturation vapor pressure at Tref (torr)
  FURR6OHOOH.Tref=298;         		// Temperature of reference (K)
  FURR6OHOOH.deltaH=125.0;     		// Enthalpy of vaporization (kJ/mol)
  FURR6OHOOH.Henry=0.0;     		// Henry's law constant at Tref (M/atm)
  FURR6OHOOH.aq_type="none"; 		// "none","diacid","monoacid" or "aldehyde"
  // FURR6OHOOH.Kacidity1=XXacidityXX;   		// First acidity constant
  FURR6OHOOH.hydrophilic=true;   	// Does the species condense on the aqueous phase?
  FURR6OHOOH.hydrophobic=true;  		// Does the species condense on the organic phase?
  FURR6OHOOH.nonvolatile=false; 		// Is the compound nonvolatile?
  FURR6OHOOH.kp_from_experiment=false; 	// Use experimental partitioning constant at Tref?
  FURR6OHOOH.is_organic=true;  		// Is the compound organic?
  FURR6OHOOH.compute_gamma_org=true;  	// Compute the activity coefficients of the organic phase for this compound?
  FURR6OHOOH.compute_gamma_aq=true;  	// Compute the activity coefficients of the aqueous phase for this compound?
  FURR6OHOOH.Koligo_org=0.0;      	//oligomeriation constant in the organic phase
  FURR6OHOOH.rho=1300.;  
  FURR6OHOOH.is_monomer=false;
  FURR6OHOOH.rion=false;
  FURR6OHOOH.KDiffusion_air=1.0e-5;
  //  FURR6OHOOH.accomodation_coefficient=alpha;
  FURR6OHOOH.viscosity=1.68e12;
  FURR6OHOOH.is_solid=false;
  FURR6OHOOH.is_generic=false;
  
  //Group: if no functionnal group in the species use the default species
  //for the computation of activity coefficients
  //
  double group_tmp_FURR6OHOOH [] = {0.,0.,0.,0., // group C
				  0.,0.,0.,0., //group C[OH]
				  0.,0.,0.,0., //group Calcohol
				  0.,0.,0.,0., //group Calcohol-tail
				  0.,0.,0.,0.,0., //group C=C
				  0.,0., //group aromatic carbon (AC)
				  0.,0.,0., // group //AC-C
				  1.,  //group OH
				  0., //group H2O
				  0., //group ACOH
				  0.,1., //group ketone
				  0.,   //group aldehyde  
				  0.,1., //group ester
				  0.,0.,0., //group ether 
				  0.,  //group acid
				  0.,   //group ACNO2
				  0.,0.,0., //group NO3
				  0.,0.,1., //group CO-OH
  			          0.,0.,0.,0.,0.,0.,0.,0.,0., //group CO-OC
			          0.,  //group PAN
			          0.,  //group CO-OOH
			          0.,  //group O=COC=O
			          0.,0.,0.}; //group CHxNO2
  
  size = sizeof(group_tmp_FURR6OHOOH)/sizeof(double);
  assert(size = 60);
	
  for(int i = 0; i < size; ++i)
	FURR6OHOOH.groups[i] = group_tmp_FURR6OHOOH[i];

  // Search the species name in the aerosol species list 
  // and add the species if its name matches with
  // the given list.
  add_species_ssh(surrogate, FURR6OHOOH, species_list_aer, molecular_weight_aer,
              accomodation_coefficient,diffusion_coef,saturation_vapor_pressure,enthalpy_vaporization, henry, t_ref, mass_density, species_part,nlayer,i_hydrophilic,
		  N_inert, N_inorganic);



/* ==== GH5002 ==== */ 

  species GH5002;
  GH5002.name="GH5002";
  GH5002.is_inorganic_precursor=false;
  GH5002.Psat_ref=1.491e-7; 		// Saturation vapor pressure at Tref (torr)
  GH5002.Tref=298;         		// Temperature of reference (K)
  GH5002.deltaH=126.79;     		// Enthalpy of vaporization (kJ/mol)
  GH5002.Henry=0.0;     		// Henry's law constant at Tref (M/atm)
  GH5002.aq_type="none"; 		// "none","diacid","monoacid" or "aldehyde"
  // GH5002.Kacidity1=XXacidityXX;   		// First acidity constant
  GH5002.hydrophilic=true;   	// Does the species condense on the aqueous phase?
  GH5002.hydrophobic=true;  		// Does the species condense on the organic phase?
  GH5002.nonvolatile=false; 		// Is the compound nonvolatile?
  GH5002.kp_from_experiment=false; 	// Use experimental partitioning constant at Tref?
  GH5002.is_organic=true;  		// Is the compound organic?
  GH5002.compute_gamma_org=true;  	// Compute the activity coefficients of the organic phase for this compound?
  GH5002.compute_gamma_aq=true;  	// Compute the activity coefficients of the aqueous phase for this compound?
  GH5002.Koligo_org=0.0;      	//oligomeriation constant in the organic phase
  GH5002.rho=1300.;  
  GH5002.is_monomer=false;
  GH5002.rion=false;
  GH5002.KDiffusion_air=1.0e-5;
  //  GH5002.accomodation_coefficient=alpha;
  GH5002.viscosity=1.68e12;
  GH5002.is_solid=false;
  GH5002.is_generic=false;
  
  //Group: if no functionnal group in the species use the default species
  //for the computation of activity coefficients
  //
  double group_tmp_GH5002 [] = {1.,0.,0.,0., // group C
				  0.,0.,0.,0., //group C[OH]
				  0.,0.,0.,0., //group Calcohol
				  0.,0.,0.,0., //group Calcohol-tail
				  0.,0.,0.,0.,0., //group C=C
				  0.,0., //group aromatic carbon (AC)
				  0.,0.,0., // group //AC-C
				  0.,  //group OH
				  0., //group H2O
				  0., //group ACOH
				  0.,0., //group ketone
				  1.,   //group aldehyde  
				  0.,0., //group ester
				  0.,0.,0., //group ether 
				  0.,  //group acid
				  0.,   //group ACNO2
				  0.,1.,0., //group NO3
				  0.,0.,1., //group CO-OH
  			          0.,0.,0.,0.,0.,0.,0.,0.,0., //group CO-OC
			          0.,  //group PAN
			          1.,  //group CO-OOH
			          0.,  //group O=COC=O
			          0.,0.,0.}; //group CHxNO2
  
  size = sizeof(group_tmp_GH5002)/sizeof(double);
  assert(size = 60);
	
  for(int i = 0; i < size; ++i)
	GH5002.groups[i] = group_tmp_GH5002[i];

  // Search the species name in the aerosol species list 
  // and add the species if its name matches with
  // the given list.
  add_species_ssh(surrogate, GH5002, species_list_aer, molecular_weight_aer,
              accomodation_coefficient,diffusion_coef,saturation_vapor_pressure,enthalpy_vaporization, henry, t_ref, mass_density, species_part,nlayer,i_hydrophilic,
		  N_inert, N_inorganic);



/* ==== HD2000 ==== */ 

  species HD2000;
  HD2000.name="HD2000";
  HD2000.is_inorganic_precursor=false;
  HD2000.Psat_ref=1.205; 		// Saturation vapor pressure at Tref (torr)
  HD2000.Tref=298;         		// Temperature of reference (K)
  HD2000.deltaH=50.0;     		// Enthalpy of vaporization (kJ/mol)
  HD2000.Henry=0.0;     		// Henry's law constant at Tref (M/atm)
  HD2000.aq_type="none"; 		// "none","diacid","monoacid" or "aldehyde"
  // HD2000.Kacidity1=XXacidityXX;   		// First acidity constant
  HD2000.hydrophilic=true;   	// Does the species condense on the aqueous phase?
  HD2000.hydrophobic=true;  		// Does the species condense on the organic phase?
  HD2000.nonvolatile=false; 		// Is the compound nonvolatile?
  HD2000.kp_from_experiment=false; 	// Use experimental partitioning constant at Tref?
  HD2000.is_organic=true;  		// Is the compound organic?
  HD2000.compute_gamma_org=true;  	// Compute the activity coefficients of the organic phase for this compound?
  HD2000.compute_gamma_aq=true;  	// Compute the activity coefficients of the aqueous phase for this compound?
  HD2000.Koligo_org=0.0;      	//oligomeriation constant in the organic phase
  HD2000.rho=1300.;  
  HD2000.is_monomer=false;
  HD2000.rion=false;
  HD2000.KDiffusion_air=1.0e-5;
  //  HD2000.accomodation_coefficient=alpha;
  HD2000.viscosity=1.68e12;
  HD2000.is_solid=false;
  HD2000.is_generic=false;
  
  //Group: if no functionnal group in the species use the default species
  //for the computation of activity coefficients
  //
  double group_tmp_HD2000 [] = {0.,0.,0.,0., // group C
				  0.,0.,0.,0., //group C[OH]
				  0.,0.,0.,0., //group Calcohol
				  0.,0.,0.,0., //group Calcohol-tail
				  0.,0.,0.,0.,0., //group C=C
				  0.,0., //group aromatic carbon (AC)
				  0.,0.,0., // group //AC-C
				  0.,  //group OH
				  0., //group H2O
				  0., //group ACOH
				  0.,0., //group ketone
				  1.,   //group aldehyde  
				  0.,0., //group ester
				  0.,0.,0., //group ether 
				  0.,  //group acid
				  0.,   //group ACNO2
				  0.,0.,0., //group NO3
				  1.,0.,0., //group CO-OH
  			          0.,0.,0.,0.,0.,0.,0.,0.,0., //group CO-OC
			          0.,  //group PAN
			          0.,  //group CO-OOH
			          0.,  //group O=COC=O
			          0.,0.,0.}; //group CHxNO2
  
  size = sizeof(group_tmp_HD2000)/sizeof(double);
  assert(size = 60);
	
  for(int i = 0; i < size; ++i)
	HD2000.groups[i] = group_tmp_HD2000[i];

  // Search the species name in the aerosol species list 
  // and add the species if its name matches with
  // the given list.
  add_species_ssh(surrogate, HD2000, species_list_aer, molecular_weight_aer,
              accomodation_coefficient,diffusion_coef,saturation_vapor_pressure,enthalpy_vaporization, henry, t_ref, mass_density, species_part,nlayer,i_hydrophilic,
		  N_inert, N_inorganic);



/* ==== HOM1O ==== */ 

  species HOM1O;
  HOM1O.name="HOM1O";
  HOM1O.is_inorganic_precursor=false;
  HOM1O.Psat_ref=3.256e-5; 		// Saturation vapor pressure at Tref (torr)
  HOM1O.Tref=298;         		// Temperature of reference (K)
  HOM1O.deltaH=100.0;     		// Enthalpy of vaporization (kJ/mol)
  HOM1O.Henry=0.0;     		// Henry's law constant at Tref (M/atm)
  HOM1O.aq_type="none"; 		// "none","diacid","monoacid" or "aldehyde"
  // HOM1O.Kacidity1=XXacidityXX;   		// First acidity constant
  HOM1O.hydrophilic=true;   	// Does the species condense on the aqueous phase?
  HOM1O.hydrophobic=true;  		// Does the species condense on the organic phase?
  HOM1O.nonvolatile=false; 		// Is the compound nonvolatile?
  HOM1O.kp_from_experiment=false; 	// Use experimental partitioning constant at Tref?
  HOM1O.is_organic=true;  		// Is the compound organic?
  HOM1O.compute_gamma_org=true;  	// Compute the activity coefficients of the organic phase for this compound?
  HOM1O.compute_gamma_aq=true;  	// Compute the activity coefficients of the aqueous phase for this compound?
  HOM1O.Koligo_org=0.0;      	//oligomeriation constant in the organic phase
  HOM1O.rho=1300.;  
  HOM1O.is_monomer=false;
  HOM1O.rion=false;
  HOM1O.KDiffusion_air=1.0e-5;
  //  HOM1O.accomodation_coefficient=alpha;
  HOM1O.viscosity=1.68e12;
  HOM1O.is_solid=false;
  HOM1O.is_generic=false;
  
  //Group: if no functionnal group in the species use the default species
  //for the computation of activity coefficients
  //
  double group_tmp_HOM1O [] = {0.,0.,0.,0., // group C
				  0.,0.,0.,0., //group C[OH]
				  0.,0.,0.,0., //group Calcohol
				  0.,0.,0.,0., //group Calcohol-tail
				  0.,0.,0.,1.,0., //group C=C
				  0.,0., //group aromatic carbon (AC)
				  0.,0.,0., // group //AC-C
				  1.,  //group OH
				  0., //group H2O
				  0., //group ACOH
				  0.,1., //group ketone
				  0.,   //group aldehyde  
				  0.,0., //group ester
				  0.,0.,0., //group ether 
				  0.,  //group acid
				  0.,   //group ACNO2
				  0.,0.,0., //group NO3
				  1.,0.,0., //group CO-OH
  			          0.,0.,0.,0.,1.,0.,0.,0.,0., //group CO-OC
			          0.,  //group PAN
			          0.,  //group CO-OOH
			          0.,  //group O=COC=O
			          0.,0.,0.}; //group CHxNO2
  
  size = sizeof(group_tmp_HOM1O)/sizeof(double);
  assert(size = 60);
	
  for(int i = 0; i < size; ++i)
	HOM1O.groups[i] = group_tmp_HOM1O[i];

  // Search the species name in the aerosol species list 
  // and add the species if its name matches with
  // the given list.
  add_species_ssh(surrogate, HOM1O, species_list_aer, molecular_weight_aer,
              accomodation_coefficient,diffusion_coef,saturation_vapor_pressure,enthalpy_vaporization, henry, t_ref, mass_density, species_part,nlayer,i_hydrophilic,
		  N_inert, N_inorganic);



/* ==== HOM1ONO2 ==== */ 

  species HOM1ONO2;
  HOM1ONO2.name="HOM1ONO2";
  HOM1ONO2.is_inorganic_precursor=false;
  HOM1ONO2.Psat_ref=1.516e-6; 		// Saturation vapor pressure at Tref (torr)
  HOM1ONO2.Tref=298;         		// Temperature of reference (K)
  HOM1ONO2.deltaH=105.0;     		// Enthalpy of vaporization (kJ/mol)
  HOM1ONO2.Henry=0.0;     		// Henry's law constant at Tref (M/atm)
  HOM1ONO2.aq_type="none"; 		// "none","diacid","monoacid" or "aldehyde"
  // HOM1ONO2.Kacidity1=XXacidityXX;   		// First acidity constant
  HOM1ONO2.hydrophilic=true;   	// Does the species condense on the aqueous phase?
  HOM1ONO2.hydrophobic=true;  		// Does the species condense on the organic phase?
  HOM1ONO2.nonvolatile=false; 		// Is the compound nonvolatile?
  HOM1ONO2.kp_from_experiment=false; 	// Use experimental partitioning constant at Tref?
  HOM1ONO2.is_organic=true;  		// Is the compound organic?
  HOM1ONO2.compute_gamma_org=true;  	// Compute the activity coefficients of the organic phase for this compound?
  HOM1ONO2.compute_gamma_aq=true;  	// Compute the activity coefficients of the aqueous phase for this compound?
  HOM1ONO2.Koligo_org=0.0;      	//oligomeriation constant in the organic phase
  HOM1ONO2.rho=1300.;  
  HOM1ONO2.is_monomer=false;
  HOM1ONO2.rion=false;
  HOM1ONO2.KDiffusion_air=1.0e-5;
  //  HOM1ONO2.accomodation_coefficient=alpha;
  HOM1ONO2.viscosity=1.68e12;
  HOM1ONO2.is_solid=false;
  HOM1ONO2.is_generic=false;
  
  //Group: if no functionnal group in the species use the default species
  //for the computation of activity coefficients
  //
  double group_tmp_HOM1ONO2 [] = {0.,0.,0.,0., // group C
				  0.,0.,1.,0., //group C[OH]
				  0.,0.,0.,0., //group Calcohol
				  0.,0.,0.,0., //group Calcohol-tail
				  0.,0.,0.,1.,0., //group C=C
				  0.,0., //group aromatic carbon (AC)
				  0.,0.,0., // group //AC-C
				  1.,  //group OH
				  0., //group H2O
				  0., //group ACOH
				  0.,0., //group ketone
				  0.,   //group aldehyde  
				  0.,0., //group ester
				  0.,0.,0., //group ether 
				  0.,  //group acid
				  0.,   //group ACNO2
				  0.,1.,0., //group NO3
				  1.,0.,0., //group CO-OH
  			          0.,0.,0.,0.,0.,0.,1.,0.,0., //group CO-OC
			          0.,  //group PAN
			          0.,  //group CO-OOH
			          0.,  //group O=COC=O
			          0.,0.,0.}; //group CHxNO2
  
  size = sizeof(group_tmp_HOM1ONO2)/sizeof(double);
  assert(size = 60);
	
  for(int i = 0; i < size; ++i)
	HOM1ONO2.groups[i] = group_tmp_HOM1ONO2[i];

  // Search the species name in the aerosol species list 
  // and add the species if its name matches with
  // the given list.
  add_species_ssh(surrogate, HOM1ONO2, species_list_aer, molecular_weight_aer,
              accomodation_coefficient,diffusion_coef,saturation_vapor_pressure,enthalpy_vaporization, henry, t_ref, mass_density, species_part,nlayer,i_hydrophilic,
		  N_inert, N_inorganic);



/* ==== HOM1OOH ==== */ 

  species HOM1OOH;
  HOM1OOH.name="HOM1OOH";
  HOM1OOH.is_inorganic_precursor=false;
  HOM1OOH.Psat_ref=1.516e-7; 		// Saturation vapor pressure at Tref (torr)
  HOM1OOH.Tref=298;         		// Temperature of reference (K)
  HOM1OOH.deltaH=122.0;     		// Enthalpy of vaporization (kJ/mol)
  HOM1OOH.Henry=0.0;     		// Henry's law constant at Tref (M/atm)
  HOM1OOH.aq_type="none"; 		// "none","diacid","monoacid" or "aldehyde"
  // HOM1OOH.Kacidity1=XXacidityXX;   		// First acidity constant
  HOM1OOH.hydrophilic=true;   	// Does the species condense on the aqueous phase?
  HOM1OOH.hydrophobic=true;  		// Does the species condense on the organic phase?
  HOM1OOH.nonvolatile=false; 		// Is the compound nonvolatile?
  HOM1OOH.kp_from_experiment=false; 	// Use experimental partitioning constant at Tref?
  HOM1OOH.is_organic=true;  		// Is the compound organic?
  HOM1OOH.compute_gamma_org=true;  	// Compute the activity coefficients of the organic phase for this compound?
  HOM1OOH.compute_gamma_aq=true;  	// Compute the activity coefficients of the aqueous phase for this compound?
  HOM1OOH.Koligo_org=0.0;      	//oligomeriation constant in the organic phase
  HOM1OOH.rho=1300.;  
  HOM1OOH.is_monomer=false;
  HOM1OOH.rion=false;
  HOM1OOH.KDiffusion_air=1.0e-5;
  //  HOM1OOH.accomodation_coefficient=alpha;
  HOM1OOH.viscosity=1.68e12;
  HOM1OOH.is_solid=false;
  HOM1OOH.is_generic=false;
  
  //Group: if no functionnal group in the species use the default species
  //for the computation of activity coefficients
  //
  double group_tmp_HOM1OOH [] = {0.,0.,0.,0., // group C
				  0.,0.,1.,0., //group C[OH]
				  0.,0.,0.,0., //group Calcohol
				  0.,0.,0.,0., //group Calcohol-tail
				  0.,0.,0.,1.,0., //group C=C
				  0.,0., //group aromatic carbon (AC)
				  0.,0.,0., // group //AC-C
				  1.,  //group OH
				  0., //group H2O
				  0., //group ACOH
				  0.,0., //group ketone
				  0.,   //group aldehyde  
				  0.,0., //group ester
				  0.,0.,0., //group ether 
				  0.,  //group acid
				  0.,   //group ACNO2
				  0.,0.,0., //group NO3
				  1.,1.,0., //group CO-OH
  			          0.,0.,0.,0.,0.,0.,1.,0.,0., //group CO-OC
			          0.,  //group PAN
			          0.,  //group CO-OOH
			          0.,  //group O=COC=O
			          0.,0.,0.}; //group CHxNO2
  
  size = sizeof(group_tmp_HOM1OOH)/sizeof(double);
  assert(size = 60);
	
  for(int i = 0; i < size; ++i)
	HOM1OOH.groups[i] = group_tmp_HOM1OOH[i];

  // Search the species name in the aerosol species list 
  // and add the species if its name matches with
  // the given list.
  add_species_ssh(surrogate, HOM1OOH, species_list_aer, molecular_weight_aer,
              accomodation_coefficient,diffusion_coef,saturation_vapor_pressure,enthalpy_vaporization, henry, t_ref, mass_density, species_part,nlayer,i_hydrophilic,
		  N_inert, N_inorganic);



/* ==== HOM2O ==== */ 

  species HOM2O;
  HOM2O.name="HOM2O";
  HOM2O.is_inorganic_precursor=false;
  HOM2O.Psat_ref=9.568e-9; 		// Saturation vapor pressure at Tref (torr)
  HOM2O.Tref=298;         		// Temperature of reference (K)
  HOM2O.deltaH=135.0;     		// Enthalpy of vaporization (kJ/mol)
  HOM2O.Henry=0.0;     		// Henry's law constant at Tref (M/atm)
  HOM2O.aq_type="none"; 		// "none","diacid","monoacid" or "aldehyde"
  // HOM2O.Kacidity1=XXacidityXX;   		// First acidity constant
  HOM2O.hydrophilic=true;   	// Does the species condense on the aqueous phase?
  HOM2O.hydrophobic=true;  		// Does the species condense on the organic phase?
  HOM2O.nonvolatile=false; 		// Is the compound nonvolatile?
  HOM2O.kp_from_experiment=false; 	// Use experimental partitioning constant at Tref?
  HOM2O.is_organic=true;  		// Is the compound organic?
  HOM2O.compute_gamma_org=true;  	// Compute the activity coefficients of the organic phase for this compound?
  HOM2O.compute_gamma_aq=true;  	// Compute the activity coefficients of the aqueous phase for this compound?
  HOM2O.Koligo_org=0.0;      	//oligomeriation constant in the organic phase
  HOM2O.rho=1300.;  
  HOM2O.is_monomer=false;
  HOM2O.rion=false;
  HOM2O.KDiffusion_air=1.0e-5;
  //  HOM2O.accomodation_coefficient=alpha;
  HOM2O.viscosity=1.68e12;
  HOM2O.is_solid=false;
  HOM2O.is_generic=false;
  
  //Group: if no functionnal group in the species use the default species
  //for the computation of activity coefficients
  //
  double group_tmp_HOM2O [] = {0.,0.,0.,0., // group C
				  0.,0.,0.,0., //group C[OH]
				  0.,0.,0.,0., //group Calcohol
				  0.,0.,0.,0., //group Calcohol-tail
				  0.,0.,0.,0.,1., //group C=C
				  0.,0., //group aromatic carbon (AC)
				  0.,0.,0., // group //AC-C
				  1.,  //group OH
				  0., //group H2O
				  0., //group ACOH
				  0.,1., //group ketone
				  0.,   //group aldehyde  
				  0.,0., //group ester
				  0.,0.,0., //group ether 
				  0.,  //group acid
				  0.,   //group ACNO2
				  0.,0.,0., //group NO3
				  1.,1.,0., //group CO-OH
  			          0.,0.,0.,0.,0.,0.,1.,0.,0., //group CO-OC
			          0.,  //group PAN
			          0.,  //group CO-OOH
			          0.,  //group O=COC=O
			          0.,0.,0.}; //group CHxNO2
  
  size = sizeof(group_tmp_HOM2O)/sizeof(double);
  assert(size = 60);
	
  for(int i = 0; i < size; ++i)
	HOM2O.groups[i] = group_tmp_HOM2O[i];

  // Search the species name in the aerosol species list 
  // and add the species if its name matches with
  // the given list.
  add_species_ssh(surrogate, HOM2O, species_list_aer, molecular_weight_aer,
              accomodation_coefficient,diffusion_coef,saturation_vapor_pressure,enthalpy_vaporization, henry, t_ref, mass_density, species_part,nlayer,i_hydrophilic,
		  N_inert, N_inorganic);



/* ==== HOM2ONO2 ==== */ 

  species HOM2ONO2;
  HOM2ONO2.name="HOM2ONO2";
  HOM2ONO2.is_inorganic_precursor=false;
  HOM2ONO2.Psat_ref=3.026e-10; 		// Saturation vapor pressure at Tref (torr)
  HOM2ONO2.Tref=298;         		// Temperature of reference (K)
  HOM2ONO2.deltaH=150.0;     		// Enthalpy of vaporization (kJ/mol)
  HOM2ONO2.Henry=0.0;     		// Henry's law constant at Tref (M/atm)
  HOM2ONO2.aq_type="none"; 		// "none","diacid","monoacid" or "aldehyde"
  // HOM2ONO2.Kacidity1=XXacidityXX;   		// First acidity constant
  HOM2ONO2.hydrophilic=true;   	// Does the species condense on the aqueous phase?
  HOM2ONO2.hydrophobic=true;  		// Does the species condense on the organic phase?
  HOM2ONO2.nonvolatile=false; 		// Is the compound nonvolatile?
  HOM2ONO2.kp_from_experiment=false; 	// Use experimental partitioning constant at Tref?
  HOM2ONO2.is_organic=true;  		// Is the compound organic?
  HOM2ONO2.compute_gamma_org=true;  	// Compute the activity coefficients of the organic phase for this compound?
  HOM2ONO2.compute_gamma_aq=true;  	// Compute the activity coefficients of the aqueous phase for this compound?
  HOM2ONO2.Koligo_org=0.0;      	//oligomeriation constant in the organic phase
  HOM2ONO2.rho=1300.;  
  HOM2ONO2.is_monomer=false;
  HOM2ONO2.rion=false;
  HOM2ONO2.KDiffusion_air=1.0e-5;
  //  HOM2ONO2.accomodation_coefficient=alpha;
  HOM2ONO2.viscosity=1.68e12;
  HOM2ONO2.is_solid=false;
  HOM2ONO2.is_generic=false;
  
  //Group: if no functionnal group in the species use the default species
  //for the computation of activity coefficients
  //
  double group_tmp_HOM2ONO2 [] = {0.,0.,0.,0., // group C
				  0.,0.,0.,0., //group C[OH]
				  0.,0.,0.,0., //group Calcohol
				  0.,0.,0.,0., //group Calcohol-tail
				  0.,0.,0.,0.,1., //group C=C
				  0.,0., //group aromatic carbon (AC)
				  0.,0.,0., // group //AC-C
				  1.,  //group OH
				  0., //group H2O
				  0., //group ACOH
				  0.,0., //group ketone
				  0.,   //group aldehyde  
				  0.,0., //group ester
				  0.,0.,0., //group ether 
				  0.,  //group acid
				  0.,   //group ACNO2
				  0.,1.,0., //group NO3
				  1.,1.,0., //group CO-OH
  			          0.,0.,0.,0.,0.,0.,1.,0.,0., //group CO-OC
			          0.,  //group PAN
			          0.,  //group CO-OOH
			          0.,  //group O=COC=O
			          0.,0.,0.}; //group CHxNO2
  
  size = sizeof(group_tmp_HOM2ONO2)/sizeof(double);
  assert(size = 60);
	
  for(int i = 0; i < size; ++i)
	HOM2ONO2.groups[i] = group_tmp_HOM2ONO2[i];

  // Search the species name in the aerosol species list 
  // and add the species if its name matches with
  // the given list.
  add_species_ssh(surrogate, HOM2ONO2, species_list_aer, molecular_weight_aer,
              accomodation_coefficient,diffusion_coef,saturation_vapor_pressure,enthalpy_vaporization, henry, t_ref, mass_density, species_part,nlayer,i_hydrophilic,
		  N_inert, N_inorganic);



/* ==== HOM2OOH ==== */ 

  species HOM2OOH;
  HOM2OOH.name="HOM2OOH";
  HOM2OOH.is_inorganic_precursor=false;
  HOM2OOH.Psat_ref=1.909e-11; 		// Saturation vapor pressure at Tref (torr)
  HOM2OOH.Tref=298;         		// Temperature of reference (K)
  HOM2OOH.deltaH=160.0;     		// Enthalpy of vaporization (kJ/mol)
  HOM2OOH.Henry=0.0;     		// Henry's law constant at Tref (M/atm)
  HOM2OOH.aq_type="none"; 		// "none","diacid","monoacid" or "aldehyde"
  // HOM2OOH.Kacidity1=XXacidityXX;   		// First acidity constant
  HOM2OOH.hydrophilic=true;   	// Does the species condense on the aqueous phase?
  HOM2OOH.hydrophobic=true;  		// Does the species condense on the organic phase?
  HOM2OOH.nonvolatile=false; 		// Is the compound nonvolatile?
  HOM2OOH.kp_from_experiment=false; 	// Use experimental partitioning constant at Tref?
  HOM2OOH.is_organic=true;  		// Is the compound organic?
  HOM2OOH.compute_gamma_org=true;  	// Compute the activity coefficients of the organic phase for this compound?
  HOM2OOH.compute_gamma_aq=true;  	// Compute the activity coefficients of the aqueous phase for this compound?
  HOM2OOH.Koligo_org=0.0;      	//oligomeriation constant in the organic phase
  HOM2OOH.rho=1300.;  
  HOM2OOH.is_monomer=false;
  HOM2OOH.rion=false;
  HOM2OOH.KDiffusion_air=1.0e-5;
  //  HOM2OOH.accomodation_coefficient=alpha;
  HOM2OOH.viscosity=1.68e12;
  HOM2OOH.is_solid=false;
  HOM2OOH.is_generic=false;
  
  //Group: if no functionnal group in the species use the default species
  //for the computation of activity coefficients
  //
  double group_tmp_HOM2OOH [] = {0.,0.,0.,0., // group C
				  0.,0.,0.,0., //group C[OH]
				  0.,0.,0.,0., //group Calcohol
				  0.,0.,0.,0., //group Calcohol-tail
				  0.,0.,0.,0.,1., //group C=C
				  0.,0., //group aromatic carbon (AC)
				  0.,0.,0., // group //AC-C
				  1.,  //group OH
				  0., //group H2O
				  0., //group ACOH
				  0.,0., //group ketone
				  0.,   //group aldehyde  
				  0.,0., //group ester
				  0.,0.,0., //group ether 
				  0.,  //group acid
				  0.,   //group ACNO2
				  0.,0.,0., //group NO3
				  1.,2.,0., //group CO-OH
  			          0.,0.,0.,0.,0.,0.,1.,0.,0., //group CO-OC
			          0.,  //group PAN
			          0.,  //group CO-OOH
			          0.,  //group O=COC=O
			          0.,0.,0.}; //group CHxNO2
  
  size = sizeof(group_tmp_HOM2OOH)/sizeof(double);
  assert(size = 60);
	
  for(int i = 0; i < size; ++i)
	HOM2OOH.groups[i] = group_tmp_HOM2OOH[i];

  // Search the species name in the aerosol species list 
  // and add the species if its name matches with
  // the given list.
  add_species_ssh(surrogate, HOM2OOH, species_list_aer, molecular_weight_aer,
              accomodation_coefficient,diffusion_coef,saturation_vapor_pressure,enthalpy_vaporization, henry, t_ref, mass_density, species_part,nlayer,i_hydrophilic,
		  N_inert, N_inorganic);



/* ==== MALAHY ==== */ 

  species MALAHY;
  MALAHY.name="MALAHY";
  MALAHY.is_inorganic_precursor=false;
  MALAHY.Psat_ref=2.563e-3; 		// Saturation vapor pressure at Tref (torr)
  MALAHY.Tref=298;         		// Temperature of reference (K)
  MALAHY.deltaH=75.42;     		// Enthalpy of vaporization (kJ/mol)
  MALAHY.Henry=0.0;     		// Henry's law constant at Tref (M/atm)
  MALAHY.aq_type="none"; 		// "none","diacid","monoacid" or "aldehyde"
  // MALAHY.Kacidity1=XXacidityXX;   		// First acidity constant
  MALAHY.hydrophilic=true;   	// Does the species condense on the aqueous phase?
  MALAHY.hydrophobic=true;  		// Does the species condense on the organic phase?
  MALAHY.nonvolatile=false; 		// Is the compound nonvolatile?
  MALAHY.kp_from_experiment=false; 	// Use experimental partitioning constant at Tref?
  MALAHY.is_organic=true;  		// Is the compound organic?
  MALAHY.compute_gamma_org=true;  	// Compute the activity coefficients of the organic phase for this compound?
  MALAHY.compute_gamma_aq=true;  	// Compute the activity coefficients of the aqueous phase for this compound?
  MALAHY.Koligo_org=0.0;      	//oligomeriation constant in the organic phase
  MALAHY.rho=1300.;  
  MALAHY.is_monomer=false;
  MALAHY.rion=false;
  MALAHY.KDiffusion_air=1.0e-5;
  //  MALAHY.accomodation_coefficient=alpha;
  MALAHY.viscosity=1.68e12;
  MALAHY.is_solid=false;
  MALAHY.is_generic=false;
  
  //Group: if no functionnal group in the species use the default species
  //for the computation of activity coefficients
  //
  double group_tmp_MALAHY [] = {0.,0.,0.,0., // group C
				  0.,0.,0.,0., //group C[OH]
				  0.,0.,0.,0., //group Calcohol
				  0.,0.,0.,0., //group Calcohol-tail
				  0.,1.,0.,0.,0., //group C=C
				  0.,0., //group aromatic carbon (AC)
				  0.,0.,0., // group //AC-C
				  0.,  //group OH
				  0., //group H2O
				  0., //group ACOH
				  0.,1., //group ketone
				  0.,   //group aldehyde  
				  0.,1., //group ester
				  0.,0.,0., //group ether 
				  0.,  //group acid
				  0.,   //group ACNO2
				  0.,0.,0., //group NO3
				  0.,0.,0., //group CO-OH
  			          0.,0.,0.,0.,0.,0.,0.,0.,0., //group CO-OC
			          0.,  //group PAN
			          0.,  //group CO-OOH
			          0.,  //group O=COC=O
			          0.,0.,0.}; //group CHxNO2
  
  size = sizeof(group_tmp_MALAHY)/sizeof(double);
  assert(size = 60);
	
  for(int i = 0; i < size; ++i)
	MALAHY.groups[i] = group_tmp_MALAHY[i];

  // Search the species name in the aerosol species list 
  // and add the species if its name matches with
  // the given list.
  add_species_ssh(surrogate, MALAHY, species_list_aer, molecular_weight_aer,
              accomodation_coefficient,diffusion_coef,saturation_vapor_pressure,enthalpy_vaporization, henry, t_ref, mass_density, species_part,nlayer,i_hydrophilic,
		  N_inert, N_inorganic);



/* ==== MBQN1K1OH ==== */ 

  species MBQN1K1OH;
  MBQN1K1OH.name="MBQN1K1OH";
  MBQN1K1OH.is_inorganic_precursor=false;
  MBQN1K1OH.Psat_ref=0.000191; 		// Saturation vapor pressure at Tref (torr)
  MBQN1K1OH.Tref=298;         		// Temperature of reference (K)
  MBQN1K1OH.deltaH=91.87;     		// Enthalpy of vaporization (kJ/mol)
  MBQN1K1OH.Henry=0.0;     		// Henry's law constant at Tref (M/atm)
  MBQN1K1OH.aq_type="none"; 		// "none","diacid","monoacid" or "aldehyde"
  // MBQN1K1OH.Kacidity1=XXacidityXX;   		// First acidity constant
  MBQN1K1OH.hydrophilic=true;   	// Does the species condense on the aqueous phase?
  MBQN1K1OH.hydrophobic=true;  		// Does the species condense on the organic phase?
  MBQN1K1OH.nonvolatile=false; 		// Is the compound nonvolatile?
  MBQN1K1OH.kp_from_experiment=false; 	// Use experimental partitioning constant at Tref?
  MBQN1K1OH.is_organic=true;  		// Is the compound organic?
  MBQN1K1OH.compute_gamma_org=true;  	// Compute the activity coefficients of the organic phase for this compound?
  MBQN1K1OH.compute_gamma_aq=true;  	// Compute the activity coefficients of the aqueous phase for this compound?
  MBQN1K1OH.Koligo_org=0.0;      	//oligomeriation constant in the organic phase
  MBQN1K1OH.rho=1300.;  
  MBQN1K1OH.is_monomer=false;
  MBQN1K1OH.rion=false;
  MBQN1K1OH.KDiffusion_air=1.0e-5;
  //  MBQN1K1OH.accomodation_coefficient=alpha;
  MBQN1K1OH.viscosity=1.68e12;
  MBQN1K1OH.is_solid=false;
  MBQN1K1OH.is_generic=false;
  
  //Group: if no functionnal group in the species use the default species
  //for the computation of activity coefficients
  //
  double group_tmp_MBQN1K1OH [] = {0.,0.,0.,0., // group C
				  0.,0.,0.,0., //group C[OH]
				  0.,0.,0.,0., //group Calcohol
				  0.,0.,0.,0., //group Calcohol-tail
				  0.,0.,0.,1.,0., //group C=C
				  0.,0., //group aromatic carbon (AC)
				  0.,0.,0., // group //AC-C
				  1.,  //group OH
				  0., //group H2O
				  0., //group ACOH
				  0.,3., //group ketone
				  0.,   //group aldehyde  
				  0.,0., //group ester
				  0.,0.,0., //group ether 
				  0.,  //group acid
				  0.,   //group ACNO2
				  0.,0.,0., //group NO3
				  0.,0.,0., //group CO-OH
  			          0.,0.,0.,0.,0.,0.,0.,0.,0., //group CO-OC
			          0.,  //group PAN
			          0.,  //group CO-OOH
			          0.,  //group O=COC=O
			          0.,0.,0.}; //group CHxNO2
  
  size = sizeof(group_tmp_MBQN1K1OH)/sizeof(double);
  assert(size = 60);
	
  for(int i = 0; i < size; ++i)
	MBQN1K1OH.groups[i] = group_tmp_MBQN1K1OH[i];

  // Search the species name in the aerosol species list 
  // and add the species if its name matches with
  // the given list.
  add_species_ssh(surrogate, MBQN1K1OH, species_list_aer, molecular_weight_aer,
              accomodation_coefficient,diffusion_coef,saturation_vapor_pressure,enthalpy_vaporization, henry, t_ref, mass_density, species_part,nlayer,i_hydrophilic,
		  N_inert, N_inorganic);



/* ==== MBQN1OH ==== */ 

  species MBQN1OH;
  MBQN1OH.name="MBQN1OH";
  MBQN1OH.is_inorganic_precursor=false;
  MBQN1OH.Psat_ref=0.003809; 		// Saturation vapor pressure at Tref (torr)
  MBQN1OH.Tref=298;         		// Temperature of reference (K)
  MBQN1OH.deltaH=77.72;     		// Enthalpy of vaporization (kJ/mol)
  MBQN1OH.Henry=0.0;     		// Henry's law constant at Tref (M/atm)
  MBQN1OH.aq_type="none"; 		// "none","diacid","monoacid" or "aldehyde"
  // MBQN1OH.Kacidity1=XXacidityXX;   		// First acidity constant
  MBQN1OH.hydrophilic=true;   	// Does the species condense on the aqueous phase?
  MBQN1OH.hydrophobic=true;  		// Does the species condense on the organic phase?
  MBQN1OH.nonvolatile=false; 		// Is the compound nonvolatile?
  MBQN1OH.kp_from_experiment=false; 	// Use experimental partitioning constant at Tref?
  MBQN1OH.is_organic=true;  		// Is the compound organic?
  MBQN1OH.compute_gamma_org=true;  	// Compute the activity coefficients of the organic phase for this compound?
  MBQN1OH.compute_gamma_aq=true;  	// Compute the activity coefficients of the aqueous phase for this compound?
  MBQN1OH.Koligo_org=0.0;      	//oligomeriation constant in the organic phase
  MBQN1OH.rho=1300.;  
  MBQN1OH.is_monomer=false;
  MBQN1OH.rion=false;
  MBQN1OH.KDiffusion_air=1.0e-5;
  //  MBQN1OH.accomodation_coefficient=alpha;
  MBQN1OH.viscosity=1.68e12;
  MBQN1OH.is_solid=false;
  MBQN1OH.is_generic=false;
  
  //Group: if no functionnal group in the species use the default species
  //for the computation of activity coefficients
  //
  double group_tmp_MBQN1OH [] = {0.,0.,0.,0., // group C
				  0.,0.,0.,0., //group C[OH]
				  0.,0.,0.,0., //group Calcohol
				  0.,0.,0.,0., //group Calcohol-tail
				  0.,1.,0.,1.,0., //group C=C
				  0.,0., //group aromatic carbon (AC)
				  0.,0.,0., // group //AC-C
				  1.,  //group OH
				  0., //group H2O
				  0., //group ACOH
				  0.,2., //group ketone
				  0.,   //group aldehyde  
				  0.,0., //group ester
				  0.,0.,0., //group ether 
				  0.,  //group acid
				  0.,   //group ACNO2
				  0.,0.,0., //group NO3
				  0.,0.,0., //group CO-OH
  			          0.,0.,0.,0.,0.,0.,0.,0.,0., //group CO-OC
			          0.,  //group PAN
			          0.,  //group CO-OOH
			          0.,  //group O=COC=O
			          0.,0.,0.}; //group CHxNO2
  
  size = sizeof(group_tmp_MBQN1OH)/sizeof(double);
  assert(size = 60);
	
  for(int i = 0; i < size; ++i)
	MBQN1OH.groups[i] = group_tmp_MBQN1OH[i];

  // Search the species name in the aerosol species list 
  // and add the species if its name matches with
  // the given list.
  add_species_ssh(surrogate, MBQN1OH, species_list_aer, molecular_weight_aer,
              accomodation_coefficient,diffusion_coef,saturation_vapor_pressure,enthalpy_vaporization, henry, t_ref, mass_density, species_part,nlayer,i_hydrophilic,
		  N_inert, N_inorganic);



/* ==== MBQN2OH ==== */ 

  species MBQN2OH;
  MBQN2OH.name="MBQN2OH";
  MBQN2OH.is_inorganic_precursor=false;
  MBQN2OH.Psat_ref=3.025e-5; 		// Saturation vapor pressure at Tref (torr)
  MBQN2OH.Tref=298;         		// Temperature of reference (K)
  MBQN2OH.deltaH=100.59;     		// Enthalpy of vaporization (kJ/mol)
  MBQN2OH.Henry=0.0;     		// Henry's law constant at Tref (M/atm)
  MBQN2OH.aq_type="none"; 		// "none","diacid","monoacid" or "aldehyde"
  // MBQN2OH.Kacidity1=XXacidityXX;   		// First acidity constant
  MBQN2OH.hydrophilic=true;   	// Does the species condense on the aqueous phase?
  MBQN2OH.hydrophobic=true;  		// Does the species condense on the organic phase?
  MBQN2OH.nonvolatile=false; 		// Is the compound nonvolatile?
  MBQN2OH.kp_from_experiment=false; 	// Use experimental partitioning constant at Tref?
  MBQN2OH.is_organic=true;  		// Is the compound organic?
  MBQN2OH.compute_gamma_org=true;  	// Compute the activity coefficients of the organic phase for this compound?
  MBQN2OH.compute_gamma_aq=true;  	// Compute the activity coefficients of the aqueous phase for this compound?
  MBQN2OH.Koligo_org=0.0;      	//oligomeriation constant in the organic phase
  MBQN2OH.rho=1300.;  
  MBQN2OH.is_monomer=false;
  MBQN2OH.rion=false;
  MBQN2OH.KDiffusion_air=1.0e-5;
  //  MBQN2OH.accomodation_coefficient=alpha;
  MBQN2OH.viscosity=1.68e12;
  MBQN2OH.is_solid=false;
  MBQN2OH.is_generic=false;
  
  //Group: if no functionnal group in the species use the default species
  //for the computation of activity coefficients
  //
  double group_tmp_MBQN2OH [] = {0.,0.,0.,0., // group C
				  0.,0.,0.,0., //group C[OH]
				  0.,0.,0.,0., //group Calcohol
				  0.,0.,0.,0., //group Calcohol-tail
				  0.,0.,0.,2.,0., //group C=C
				  0.,0., //group aromatic carbon (AC)
				  0.,0.,0., // group //AC-C
				  2.,  //group OH
				  0., //group H2O
				  0., //group ACOH
				  0.,2., //group ketone
				  0.,   //group aldehyde  
				  0.,0., //group ester
				  0.,0.,0., //group ether 
				  0.,  //group acid
				  0.,   //group ACNO2
				  0.,0.,0., //group NO3
				  0.,0.,0., //group CO-OH
  			          0.,0.,0.,0.,0.,0.,0.,0.,0., //group CO-OC
			          0.,  //group PAN
			          0.,  //group CO-OOH
			          0.,  //group O=COC=O
			          0.,0.,0.}; //group CHxNO2
  
  size = sizeof(group_tmp_MBQN2OH)/sizeof(double);
  assert(size = 60);
	
  for(int i = 0; i < size; ++i)
	MBQN2OH.groups[i] = group_tmp_MBQN2OH[i];

  // Search the species name in the aerosol species list 
  // and add the species if its name matches with
  // the given list.
  add_species_ssh(surrogate, MBQN2OH, species_list_aer, molecular_weight_aer,
              accomodation_coefficient,diffusion_coef,saturation_vapor_pressure,enthalpy_vaporization, henry, t_ref, mass_density, species_part,nlayer,i_hydrophilic,
		  N_inert, N_inorganic);



/* ==== MBQN3OH ==== */ 

  species MBQN3OH;
  MBQN3OH.name="MBQN3OH";
  MBQN3OH.is_inorganic_precursor=false;
  MBQN3OH.Psat_ref=3.025e-7; 		// Saturation vapor pressure at Tref (torr)
  MBQN3OH.Tref=298;         		// Temperature of reference (K)
  MBQN3OH.deltaH=122.37;     		// Enthalpy of vaporization (kJ/mol)
  MBQN3OH.Henry=0.0;     		// Henry's law constant at Tref (M/atm)
  MBQN3OH.aq_type="none"; 		// "none","diacid","monoacid" or "aldehyde"
  // MBQN3OH.Kacidity1=XXacidityXX;   		// First acidity constant
  MBQN3OH.hydrophilic=true;   	// Does the species condense on the aqueous phase?
  MBQN3OH.hydrophobic=true;  		// Does the species condense on the organic phase?
  MBQN3OH.nonvolatile=false; 		// Is the compound nonvolatile?
  MBQN3OH.kp_from_experiment=false; 	// Use experimental partitioning constant at Tref?
  MBQN3OH.is_organic=true;  		// Is the compound organic?
  MBQN3OH.compute_gamma_org=true;  	// Compute the activity coefficients of the organic phase for this compound?
  MBQN3OH.compute_gamma_aq=true;  	// Compute the activity coefficients of the aqueous phase for this compound?
  MBQN3OH.Koligo_org=0.0;      	//oligomeriation constant in the organic phase
  MBQN3OH.rho=1300.;  
  MBQN3OH.is_monomer=false;
  MBQN3OH.rion=false;
  MBQN3OH.KDiffusion_air=1.0e-5;
  //  MBQN3OH.accomodation_coefficient=alpha;
  MBQN3OH.viscosity=1.68e12;
  MBQN3OH.is_solid=false;
  MBQN3OH.is_generic=false;
  
  //Group: if no functionnal group in the species use the default species
  //for the computation of activity coefficients
  //
  double group_tmp_MBQN3OH [] = {0.,0.,0.,0., // group C
				  0.,0.,0.,0., //group C[OH]
				  0.,0.,0.,0., //group Calcohol
				  0.,0.,0.,0., //group Calcohol-tail
				  0.,0.,0.,1.,1., //group C=C
				  0.,0., //group aromatic carbon (AC)
				  0.,0.,0., // group //AC-C
				  3.,  //group OH
				  0., //group H2O
				  0., //group ACOH
				  0.,2., //group ketone
				  0.,   //group aldehyde  
				  0.,0., //group ester
				  0.,0.,0., //group ether 
				  0.,  //group acid
				  0.,   //group ACNO2
				  0.,0.,0., //group NO3
				  0.,0.,0., //group CO-OH
  			          0.,0.,0.,0.,0.,0.,0.,0.,0., //group CO-OC
			          0.,  //group PAN
			          0.,  //group CO-OOH
			          0.,  //group O=COC=O
			          0.,0.,0.}; //group CHxNO2
  
  size = sizeof(group_tmp_MBQN3OH)/sizeof(double);
  assert(size = 60);
	
  for(int i = 0; i < size; ++i)
	MBQN3OH.groups[i] = group_tmp_MBQN3OH[i];

  // Search the species name in the aerosol species list 
  // and add the species if its name matches with
  // the given list.
  add_species_ssh(surrogate, MBQN3OH, species_list_aer, molecular_weight_aer,
              accomodation_coefficient,diffusion_coef,saturation_vapor_pressure,enthalpy_vaporization, henry, t_ref, mass_density, species_part,nlayer,i_hydrophilic,
		  N_inert, N_inorganic);



/* ==== Me6Cy1U3K ==== */ 

  species Me6Cy1U3K;
  Me6Cy1U3K.name="Me6Cy1U3K";
  Me6Cy1U3K.is_inorganic_precursor=false;
  Me6Cy1U3K.Psat_ref=2.334e-2; 		// Saturation vapor pressure at Tref (torr)
  Me6Cy1U3K.Tref=298;         		// Temperature of reference (K)
  Me6Cy1U3K.deltaH=67.10;     		// Enthalpy of vaporization (kJ/mol)
  Me6Cy1U3K.Henry=0.0;     		// Henry's law constant at Tref (M/atm)
  Me6Cy1U3K.aq_type="none"; 		// "none","diacid","monoacid" or "aldehyde"
  // Me6Cy1U3K.Kacidity1=XXacidityXX;   		// First acidity constant
  Me6Cy1U3K.hydrophilic=true;   	// Does the species condense on the aqueous phase?
  Me6Cy1U3K.hydrophobic=true;  		// Does the species condense on the organic phase?
  Me6Cy1U3K.nonvolatile=false; 		// Is the compound nonvolatile?
  Me6Cy1U3K.kp_from_experiment=false; 	// Use experimental partitioning constant at Tref?
  Me6Cy1U3K.is_organic=true;  		// Is the compound organic?
  Me6Cy1U3K.compute_gamma_org=true;  	// Compute the activity coefficients of the organic phase for this compound?
  Me6Cy1U3K.compute_gamma_aq=true;  	// Compute the activity coefficients of the aqueous phase for this compound?
  Me6Cy1U3K.Koligo_org=0.0;      	//oligomeriation constant in the organic phase
  Me6Cy1U3K.rho=1300.;  
  Me6Cy1U3K.is_monomer=false;
  Me6Cy1U3K.rion=false;
  Me6Cy1U3K.KDiffusion_air=1.0e-5;
  //  Me6Cy1U3K.accomodation_coefficient=alpha;
  Me6Cy1U3K.viscosity=1.68e12;
  Me6Cy1U3K.is_solid=false;
  Me6Cy1U3K.is_generic=false;
  
  //Group: if no functionnal group in the species use the default species
  //for the computation of activity coefficients
  //
  double group_tmp_Me6Cy1U3K [] = {0.,0.,0.,0., // group C
				  0.,0.,0.,0., //group C[OH]
				  0.,0.,0.,0., //group Calcohol
				  0.,0.,0.,0., //group Calcohol-tail
				  0.,1.,0.,0.,0., //group C=C
				  0.,0., //group aromatic carbon (AC)
				  0.,0.,0., // group //AC-C
				  0.,  //group OH
				  0., //group H2O
				  0., //group ACOH
				  0.,3., //group ketone
				  0.,   //group aldehyde  
				  0.,0., //group ester
				  0.,0.,0., //group ether 
				  0.,  //group acid
				  0.,   //group ACNO2
				  0.,0.,0., //group NO3
				  0.,0.,0., //group CO-OH
  			          0.,0.,0.,0.,0.,0.,0.,0.,0., //group CO-OC
			          0.,  //group PAN
			          0.,  //group CO-OOH
			          0.,  //group O=COC=O
			          0.,0.,0.}; //group CHxNO2
  
  size = sizeof(group_tmp_Me6Cy1U3K)/sizeof(double);
  assert(size = 60);
	
  for(int i = 0; i < size; ++i)
	Me6Cy1U3K.groups[i] = group_tmp_Me6Cy1U3K[i];

  // Search the species name in the aerosol species list 
  // and add the species if its name matches with
  // the given list.
  add_species_ssh(surrogate, Me6Cy1U3K, species_list_aer, molecular_weight_aer,
              accomodation_coefficient,diffusion_coef,saturation_vapor_pressure,enthalpy_vaporization, henry, t_ref, mass_density, species_part,nlayer,i_hydrophilic,
		  N_inert, N_inorganic);



/* ==== MFUR ==== */ 

  species MFUR;
  MFUR.name="MFUR";
  MFUR.is_inorganic_precursor=false;
  MFUR.Psat_ref=6.043e-1; 		// Saturation vapor pressure at Tref (torr)
  MFUR.Tref=298;         		// Temperature of reference (K)
  MFUR.deltaH=52.77;     		// Enthalpy of vaporization (kJ/mol)
  MFUR.Henry=0.0;     		// Henry's law constant at Tref (M/atm)
  MFUR.aq_type="none"; 		// "none","diacid","monoacid" or "aldehyde"
  // MFUR.Kacidity1=XXacidityXX;   		// First acidity constant
  MFUR.hydrophilic=true;   	// Does the species condense on the aqueous phase?
  MFUR.hydrophobic=true;  		// Does the species condense on the organic phase?
  MFUR.nonvolatile=false; 		// Is the compound nonvolatile?
  MFUR.kp_from_experiment=false; 	// Use experimental partitioning constant at Tref?
  MFUR.is_organic=true;  		// Is the compound organic?
  MFUR.compute_gamma_org=true;  	// Compute the activity coefficients of the organic phase for this compound?
  MFUR.compute_gamma_aq=true;  	// Compute the activity coefficients of the aqueous phase for this compound?
  MFUR.Koligo_org=0.0;      	//oligomeriation constant in the organic phase
  MFUR.rho=1300.;  
  MFUR.is_monomer=false;
  MFUR.rion=false;
  MFUR.KDiffusion_air=1.0e-5;
  //  MFUR.accomodation_coefficient=alpha;
  MFUR.viscosity=1.68e12;
  MFUR.is_solid=false;
  MFUR.is_generic=false;
  
  //Group: if no functionnal group in the species use the default species
  //for the computation of activity coefficients
  //
  double group_tmp_MFUR [] = {1.,0.,0.,0., // group C
				  0.,0.,0.,0., //group C[OH]
				  0.,0.,0.,0., //group Calcohol
				  0.,0.,0.,0., //group Calcohol-tail
				  0.,0.,0.,1.,0., //group C=C
				  0.,0., //group aromatic carbon (AC)
				  0.,0.,0., // group //AC-C
				  0.,  //group OH
				  0., //group H2O
				  0., //group ACOH
				  0.,0., //group ketone
				  0.,   //group aldehyde  
				  0.,1., //group ester
				  0.,0.,0., //group ether 
				  0.,  //group acid
				  0.,   //group ACNO2
				  0.,0.,0., //group NO3
				  0.,0.,0., //group CO-OH
  			          0.,0.,0.,0.,0.,0.,0.,0.,0., //group CO-OC
			          0.,  //group PAN
			          0.,  //group CO-OOH
			          0.,  //group O=COC=O
			          0.,0.,0.}; //group CHxNO2
  
  size = sizeof(group_tmp_MFUR)/sizeof(double);
  assert(size = 60);
	
  for(int i = 0; i < size; ++i)
	MFUR.groups[i] = group_tmp_MFUR[i];

  // Search the species name in the aerosol species list 
  // and add the species if its name matches with
  // the given list.
  add_species_ssh(surrogate, MFUR, species_list_aer, molecular_weight_aer,
              accomodation_coefficient,diffusion_coef,saturation_vapor_pressure,enthalpy_vaporization, henry, t_ref, mass_density, species_part,nlayer,i_hydrophilic,
		  N_inert, N_inorganic);



/* ==== NTOL ==== */ 

  species NTOL;
  NTOL.name="NTOL";
  NTOL.is_inorganic_precursor=false;
  NTOL.Psat_ref=0.120452; 		// Saturation vapor pressure at Tref (torr)
  NTOL.Tref=298;         		// Temperature of reference (K)
  NTOL.deltaH=61.38;     		// Enthalpy of vaporization (kJ/mol)
  NTOL.Henry=0.0;     		// Henry's law constant at Tref (M/atm)
  NTOL.aq_type="none"; 		// "none","diacid","monoacid" or "aldehyde"
  // NTOL.Kacidity1=XXacidityXX;   		// First acidity constant
  NTOL.hydrophilic=true;   	// Does the species condense on the aqueous phase?
  NTOL.hydrophobic=true;  		// Does the species condense on the organic phase?
  NTOL.nonvolatile=false; 		// Is the compound nonvolatile?
  NTOL.kp_from_experiment=false; 	// Use experimental partitioning constant at Tref?
  NTOL.is_organic=true;  		// Is the compound organic?
  NTOL.compute_gamma_org=true;  	// Compute the activity coefficients of the organic phase for this compound?
  NTOL.compute_gamma_aq=true;  	// Compute the activity coefficients of the aqueous phase for this compound?
  NTOL.Koligo_org=0.0;      	//oligomeriation constant in the organic phase
  NTOL.rho=1300.;  
  NTOL.is_monomer=false;
  NTOL.rion=false;
  NTOL.KDiffusion_air=1.0e-5;
  //  NTOL.accomodation_coefficient=alpha;
  NTOL.viscosity=1.68e12;
  NTOL.is_solid=false;
  NTOL.is_generic=false;
  
  //Group: if no functionnal group in the species use the default species
  //for the computation of activity coefficients
  //
  double group_tmp_NTOL [] = {0.,0.,0.,0., // group C
				  0.,0.,0.,0., //group C[OH]
				  0.,0.,0.,0., //group Calcohol
				  0.,0.,0.,0., //group Calcohol-tail
				  0.,0.,0.,0.,0., //group C=C
				  4.,0., //group aromatic carbon (AC)
				  1.,0.,0., // group //AC-C
				  0.,  //group OH
				  0., //group H2O
				  0., //group ACOH
				  0.,0., //group ketone
				  0.,   //group aldehyde  
				  0.,0., //group ester
				  0.,0.,0., //group ether 
				  0.,  //group acid
				  1.,   //group ACNO2
				  0.,0.,0., //group NO3
				  0.,0.,0., //group CO-OH
  			          0.,0.,0.,0.,0.,0.,0.,0.,0., //group CO-OC
			          0.,  //group PAN
			          0.,  //group CO-OOH
			          0.,  //group O=COC=O
			          0.,0.,0.}; //group CHxNO2
  
  size = sizeof(group_tmp_NTOL)/sizeof(double);
  assert(size = 60);
	
  for(int i = 0; i < size; ++i)
	NTOL.groups[i] = group_tmp_NTOL[i];

  // Search the species name in the aerosol species list 
  // and add the species if its name matches with
  // the given list.
  add_species_ssh(surrogate, NTOL, species_list_aer, molecular_weight_aer,
              accomodation_coefficient,diffusion_coef,saturation_vapor_pressure,enthalpy_vaporization, henry, t_ref, mass_density, species_part,nlayer,i_hydrophilic,
		  N_inert, N_inorganic);



/* ==== P02000 ==== */ 

  species P02000;
  P02000.name="P02000";
  P02000.is_inorganic_precursor=false;
  P02000.Psat_ref=2.988e1; 		// Saturation vapor pressure at Tref (torr)
  P02000.Tref=298;         		// Temperature of reference (K)
  P02000.deltaH=39.74;     		// Enthalpy of vaporization (kJ/mol)
  P02000.Henry=0.0;     		// Henry's law constant at Tref (M/atm)
  P02000.aq_type="none"; 		// "none","diacid","monoacid" or "aldehyde"
  // P02000.Kacidity1=XXacidityXX;   		// First acidity constant
  P02000.hydrophilic=true;   	// Does the species condense on the aqueous phase?
  P02000.hydrophobic=true;  		// Does the species condense on the organic phase?
  P02000.nonvolatile=false; 		// Is the compound nonvolatile?
  P02000.kp_from_experiment=false; 	// Use experimental partitioning constant at Tref?
  P02000.is_organic=true;  		// Is the compound organic?
  P02000.compute_gamma_org=true;  	// Compute the activity coefficients of the organic phase for this compound?
  P02000.compute_gamma_aq=true;  	// Compute the activity coefficients of the aqueous phase for this compound?
  P02000.Koligo_org=0.0;      	//oligomeriation constant in the organic phase
  P02000.rho=1300.;  
  P02000.is_monomer=false;
  P02000.rion=false;
  P02000.KDiffusion_air=1.0e-5;
  //  P02000.accomodation_coefficient=alpha;
  P02000.viscosity=1.68e12;
  P02000.is_solid=false;
  P02000.is_generic=false;
  
  //Group: if no functionnal group in the species use the default species
  //for the computation of activity coefficients
  //
  double group_tmp_P02000 [] = {1.,0.,0.,0., // group C
				  0.,0.,0.,0., //group C[OH]
				  0.,0.,0.,0., //group Calcohol
				  0.,0.,0.,0., //group Calcohol-tail
				  0.,0.,0.,0.,0., //group C=C
				  0.,0., //group aromatic carbon (AC)
				  0.,0.,0., // group //AC-C
				  0.,  //group OH
				  0., //group H2O
				  0., //group ACOH
				  0.,0., //group ketone
				  0.,   //group aldehyde  
				  0.,0., //group ester
				  0.,0.,0., //group ether 
				  0.,  //group acid
				  0.,   //group ACNO2
				  0.,0.,0., //group NO3
				  0.,0.,0., //group CO-OH
  			          0.,0.,0.,0.,0.,0.,0.,0.,0., //group CO-OC
			          1.,  //group PAN
			          0.,  //group CO-OOH
			          0.,  //group O=COC=O
			          0.,0.,0.}; //group CHxNO2
  
  size = sizeof(group_tmp_P02000)/sizeof(double);
  assert(size = 60);
	
  for(int i = 0; i < size; ++i)
	P02000.groups[i] = group_tmp_P02000[i];

  // Search the species name in the aerosol species list 
  // and add the species if its name matches with
  // the given list.
  add_species_ssh(surrogate, P02000, species_list_aer, molecular_weight_aer,
              accomodation_coefficient,diffusion_coef,saturation_vapor_pressure,enthalpy_vaporization, henry, t_ref, mass_density, species_part,nlayer,i_hydrophilic,
		  N_inert, N_inorganic);



/* ==== PH5000 ==== */ 

  species PH5000;
  PH5000.name="PH5000";
  PH5000.is_inorganic_precursor=false;
  PH5000.Psat_ref=1.863e-6; 		// Saturation vapor pressure at Tref (torr)
  PH5000.Tref=298;         		// Temperature of reference (K)
  PH5000.deltaH=114.13;     		// Enthalpy of vaporization (kJ/mol)
  PH5000.Henry=0.0;     		// Henry's law constant at Tref (M/atm)
  PH5000.aq_type="none"; 		// "none","diacid","monoacid" or "aldehyde"
  // PH5000.Kacidity1=XXacidityXX;   		// First acidity constant
  PH5000.hydrophilic=true;   	// Does the species condense on the aqueous phase?
  PH5000.hydrophobic=true;  		// Does the species condense on the organic phase?
  PH5000.nonvolatile=false; 		// Is the compound nonvolatile?
  PH5000.kp_from_experiment=false; 	// Use experimental partitioning constant at Tref?
  PH5000.is_organic=true;  		// Is the compound organic?
  PH5000.compute_gamma_org=true;  	// Compute the activity coefficients of the organic phase for this compound?
  PH5000.compute_gamma_aq=true;  	// Compute the activity coefficients of the aqueous phase for this compound?
  PH5000.Koligo_org=0.0;      	//oligomeriation constant in the organic phase
  PH5000.rho=1300.;  
  PH5000.is_monomer=false;
  PH5000.rion=false;
  PH5000.KDiffusion_air=1.0e-5;
  //  PH5000.accomodation_coefficient=alpha;
  PH5000.viscosity=1.68e12;
  PH5000.is_solid=false;
  PH5000.is_generic=false;
  
  //Group: if no functionnal group in the species use the default species
  //for the computation of activity coefficients
  //
  double group_tmp_PH5000 [] = {1.,0.,0.,0., // group C
				  0.,0.,1.,0., //group C[OH]
				  0.,0.,0.,0., //group Calcohol
				  0.,0.,0.,0., //group Calcohol-tail
				  0.,0.,0.,0.,0., //group C=C
				  0.,0., //group aromatic carbon (AC)
				  0.,0.,0., // group //AC-C
				  1.,  //group OH
				  0., //group H2O
				  0., //group ACOH
				  0.,0., //group ketone
				  1.,   //group aldehyde  
				  0.,0., //group ester
				  0.,0.,0., //group ether 
				  0.,  //group acid
				  0.,   //group ACNO2
				  0.,0.,0., //group NO3
				  0.,0.,1., //group CO-OH
  			          0.,0.,0.,0.,0.,0.,0.,0.,0., //group CO-OC
			          1.,  //group PAN
			          0.,  //group CO-OOH
			          0.,  //group O=COC=O
			          0.,0.,0.}; //group CHxNO2
  
  size = sizeof(group_tmp_PH5000)/sizeof(double);
  assert(size = 60);
	
  for(int i = 0; i < size; ++i)
	PH5000.groups[i] = group_tmp_PH5000[i];

  // Search the species name in the aerosol species list 
  // and add the species if its name matches with
  // the given list.
  add_species_ssh(surrogate, PH5000, species_list_aer, molecular_weight_aer,
              accomodation_coefficient,diffusion_coef,saturation_vapor_pressure,enthalpy_vaporization, henry, t_ref, mass_density, species_part,nlayer,i_hydrophilic,
		  N_inert, N_inorganic);



/* ==== PH5002 ==== */ 

  species PH5002;
  PH5002.name="PH5002";
  PH5002.is_inorganic_precursor=false;
  PH5002.Psat_ref=1.289e-6; 		// Saturation vapor pressure at Tref (torr)
  PH5002.Tref=298;         		// Temperature of reference (K)
  PH5002.deltaH=115.66;     		// Enthalpy of vaporization (kJ/mol)
  PH5002.Henry=0.0;     		// Henry's law constant at Tref (M/atm)
  PH5002.aq_type="none"; 		// "none","diacid","monoacid" or "aldehyde"
  // PH5002.Kacidity1=XXacidityXX;   		// First acidity constant
  PH5002.hydrophilic=true;   	// Does the species condense on the aqueous phase?
  PH5002.hydrophobic=true;  		// Does the species condense on the organic phase?
  PH5002.nonvolatile=false; 		// Is the compound nonvolatile?
  PH5002.kp_from_experiment=false; 	// Use experimental partitioning constant at Tref?
  PH5002.is_organic=true;  		// Is the compound organic?
  PH5002.compute_gamma_org=true;  	// Compute the activity coefficients of the organic phase for this compound?
  PH5002.compute_gamma_aq=true;  	// Compute the activity coefficients of the aqueous phase for this compound?
  PH5002.Koligo_org=0.0;      	//oligomeriation constant in the organic phase
  PH5002.rho=1300.;  
  PH5002.is_monomer=false;
  PH5002.rion=false;
  PH5002.KDiffusion_air=1.0e-5;
  //  PH5002.accomodation_coefficient=alpha;
  PH5002.viscosity=1.68e12;
  PH5002.is_solid=false;
  PH5002.is_generic=false;
  
  //Group: if no functionnal group in the species use the default species
  //for the computation of activity coefficients
  //
  double group_tmp_PH5002 [] = {1.,0.,0.,0., // group C
				  0.,0.,0.,0., //group C[OH]
				  0.,0.,0.,0., //group Calcohol
				  0.,0.,0.,0., //group Calcohol-tail
				  0.,0.,0.,0.,0., //group C=C
				  0.,0., //group aromatic carbon (AC)
				  0.,0.,0., // group //AC-C
				  0.,  //group OH
				  0., //group H2O
				  0., //group ACOH
				  0.,0., //group ketone
				  1.,   //group aldehyde  
				  0.,0., //group ester
				  0.,0.,0., //group ether 
				  0.,  //group acid
				  0.,   //group ACNO2
				  0.,1.,0., //group NO3
				  0.,0.,1., //group CO-OH
  			          0.,0.,0.,0.,0.,0.,0.,0.,0., //group CO-OC
			          1.,  //group PAN
			          0.,  //group CO-OOH
			          0.,  //group O=COC=O
			          0.,0.,0.}; //group CHxNO2
  
  size = sizeof(group_tmp_PH5002)/sizeof(double);
  assert(size = 60);
	
  for(int i = 0; i < size; ++i)
	PH5002.groups[i] = group_tmp_PH5002[i];

  // Search the species name in the aerosol species list 
  // and add the species if its name matches with
  // the given list.
  add_species_ssh(surrogate, PH5002, species_list_aer, molecular_weight_aer,
              accomodation_coefficient,diffusion_coef,saturation_vapor_pressure,enthalpy_vaporization, henry, t_ref, mass_density, species_part,nlayer,i_hydrophilic,
		  N_inert, N_inorganic);



/* ==== PH5004 ==== */ 

  species PH5004;
  PH5004.name="PH5004";
  PH5004.is_inorganic_precursor=false;
  PH5004.Psat_ref=2.021e-9; 		// Saturation vapor pressure at Tref (torr)
  PH5004.Tref=298;         		// Temperature of reference (K)
  PH5004.deltaH=146.86;     		// Enthalpy of vaporization (kJ/mol)
  PH5004.Henry=0.0;     		// Henry's law constant at Tref (M/atm)
  PH5004.aq_type="none"; 		// "none","diacid","monoacid" or "aldehyde"
  // PH5004.Kacidity1=XXacidityXX;   		// First acidity constant
  PH5004.hydrophilic=true;   	// Does the species condense on the aqueous phase?
  PH5004.hydrophobic=true;  		// Does the species condense on the organic phase?
  PH5004.nonvolatile=false; 		// Is the compound nonvolatile?
  PH5004.kp_from_experiment=false; 	// Use experimental partitioning constant at Tref?
  PH5004.is_organic=true;  		// Is the compound organic?
  PH5004.compute_gamma_org=true;  	// Compute the activity coefficients of the organic phase for this compound?
  PH5004.compute_gamma_aq=true;  	// Compute the activity coefficients of the aqueous phase for this compound?
  PH5004.Koligo_org=0.0;      	//oligomeriation constant in the organic phase
  PH5004.rho=1300.;  
  PH5004.is_monomer=false;
  PH5004.rion=false;
  PH5004.KDiffusion_air=1.0e-5;
  //  PH5004.accomodation_coefficient=alpha;
  PH5004.viscosity=1.68e12;
  PH5004.is_solid=false;
  PH5004.is_generic=false;
  
  //Group: if no functionnal group in the species use the default species
  //for the computation of activity coefficients
  //
  double group_tmp_PH5004 [] = {1.,0.,0.,0., // group C
				  0.,0.,1.,0., //group C[OH]
				  0.,0.,0.,0., //group Calcohol
				  0.,0.,0.,0., //group Calcohol-tail
				  0.,0.,0.,0.,0., //group C=C
				  0.,0., //group aromatic carbon (AC)
				  0.,0.,0., // group //AC-C
				  1.,  //group OH
				  0., //group H2O
				  0., //group ACOH
				  0.,0., //group ketone
				  1.,   //group aldehyde  
				  0.,0., //group ester
				  0.,0.,0., //group ether 
				  0.,  //group acid
				  0.,   //group ACNO2
				  0.,0.,0., //group NO3
				  0.,1.,0., //group CO-OH
  			          0.,0.,0.,0.,0.,0.,0.,0.,0., //group CO-OC
			          1.,  //group PAN
			          0.,  //group CO-OOH
			          0.,  //group O=COC=O
			          0.,0.,0.}; //group CHxNO2
  
  size = sizeof(group_tmp_PH5004)/sizeof(double);
  assert(size = 60);
	
  for(int i = 0; i < size; ++i)
	PH5004.groups[i] = group_tmp_PH5004[i];

  // Search the species name in the aerosol species list 
  // and add the species if its name matches with
  // the given list.
  add_species_ssh(surrogate, PH5004, species_list_aer, molecular_weight_aer,
              accomodation_coefficient,diffusion_coef,saturation_vapor_pressure,enthalpy_vaporization, henry, t_ref, mass_density, species_part,nlayer,i_hydrophilic,
		  N_inert, N_inorganic);



/* ==== PK5001 ==== */ 

  species PK5001;
  PK5001.name="PK5001";
  PK5001.is_inorganic_precursor=false;
  PK5001.Psat_ref=1.836e-4; 		// Saturation vapor pressure at Tref (torr)
  PK5001.Tref=298;         		// Temperature of reference (K)
  PK5001.deltaH=91.76;     		// Enthalpy of vaporization (kJ/mol)
  PK5001.Henry=0.0;     		// Henry's law constant at Tref (M/atm)
  PK5001.aq_type="none"; 		// "none","diacid","monoacid" or "aldehyde"
  // PK5001.Kacidity1=XXacidityXX;   		// First acidity constant
  PK5001.hydrophilic=true;   	// Does the species condense on the aqueous phase?
  PK5001.hydrophobic=true;  		// Does the species condense on the organic phase?
  PK5001.nonvolatile=false; 		// Is the compound nonvolatile?
  PK5001.kp_from_experiment=false; 	// Use experimental partitioning constant at Tref?
  PK5001.is_organic=true;  		// Is the compound organic?
  PK5001.compute_gamma_org=true;  	// Compute the activity coefficients of the organic phase for this compound?
  PK5001.compute_gamma_aq=true;  	// Compute the activity coefficients of the aqueous phase for this compound?
  PK5001.Koligo_org=0.0;      	//oligomeriation constant in the organic phase
  PK5001.rho=1300.;  
  PK5001.is_monomer=false;
  PK5001.rion=false;
  PK5001.KDiffusion_air=1.0e-5;
  //  PK5001.accomodation_coefficient=alpha;
  PK5001.viscosity=1.68e12;
  PK5001.is_solid=false;
  PK5001.is_generic=false;
  
  //Group: if no functionnal group in the species use the default species
  //for the computation of activity coefficients
  //
  double group_tmp_PK5001 [] = {0.,0.,0.,0., // group C
				  0.,0.,0.,0., //group C[OH]
				  0.,0.,0.,0., //group Calcohol
				  0.,0.,0.,0., //group Calcohol-tail
				  0.,0.,0.,0.,0., //group C=C
				  0.,0., //group aromatic carbon (AC)
				  0.,0.,0., // group //AC-C
				  0.,  //group OH
				  1., //group H2O
				  0., //group ACOH
				  0.,2., //group ketone
				  0.,   //group aldehyde  
				  0.,0., //group ester
				  0.,0.,0., //group ether 
				  0.,  //group acid
				  0.,   //group ACNO2
				  0.,0.,0., //group NO3
				  0.,0.,0., //group CO-OH
  			          0.,0.,0.,0.,0.,0.,0.,0.,0., //group CO-OC
			          1.,  //group PAN
			          0.,  //group CO-OOH
			          0.,  //group O=COC=O
			          0.,0.,0.}; //group CHxNO2
  
  size = sizeof(group_tmp_PK5001)/sizeof(double);
  assert(size = 60);
	
  for(int i = 0; i < size; ++i)
	PK5001.groups[i] = group_tmp_PK5001[i];

  // Search the species name in the aerosol species list 
  // and add the species if its name matches with
  // the given list.
  add_species_ssh(surrogate, PK5001, species_list_aer, molecular_weight_aer,
              accomodation_coefficient,diffusion_coef,saturation_vapor_pressure,enthalpy_vaporization, henry, t_ref, mass_density, species_part,nlayer,i_hydrophilic,
		  N_inert, N_inorganic);



/* ==== PK5003 ==== */ 

  species PK5003;
  PK5003.name="PK5003";
  PK5003.is_inorganic_precursor=false;
  PK5003.Psat_ref=3.804e-7; 		// Saturation vapor pressure at Tref (torr)
  PK5003.Tref=298;         		// Temperature of reference (K)
  PK5003.deltaH=121.03;     		// Enthalpy of vaporization (kJ/mol)
  PK5003.Henry=0.0;     		// Henry's law constant at Tref (M/atm)
  PK5003.aq_type="none"; 		// "none","diacid","monoacid" or "aldehyde"
  // PK5003.Kacidity1=XXacidityXX;   		// First acidity constant
  PK5003.hydrophilic=true;   	// Does the species condense on the aqueous phase?
  PK5003.hydrophobic=true;  		// Does the species condense on the organic phase?
  PK5003.nonvolatile=false; 		// Is the compound nonvolatile?
  PK5003.kp_from_experiment=false; 	// Use experimental partitioning constant at Tref?
  PK5003.is_organic=true;  		// Is the compound organic?
  PK5003.compute_gamma_org=true;  	// Compute the activity coefficients of the organic phase for this compound?
  PK5003.compute_gamma_aq=true;  	// Compute the activity coefficients of the aqueous phase for this compound?
  PK5003.Koligo_org=0.0;      	//oligomeriation constant in the organic phase
  PK5003.rho=1300.;  
  PK5003.is_monomer=false;
  PK5003.rion=false;
  PK5003.KDiffusion_air=1.0e-5;
  //  PK5003.accomodation_coefficient=alpha;
  PK5003.viscosity=1.68e12;
  PK5003.is_solid=false;
  PK5003.is_generic=false;
  
  //Group: if no functionnal group in the species use the default species
  //for the computation of activity coefficients
  //
  double group_tmp_PK5003 [] = {0.,0.,0.,0., // group C
				  0.,0.,2.,0., //group C[OH]
				  0.,0.,0.,0., //group Calcohol
				  0.,0.,0.,0., //group Calcohol-tail
				  0.,0.,0.,0.,0., //group C=C
				  0.,0., //group aromatic carbon (AC)
				  0.,0.,0., // group //AC-C
				  3.,  //group OH
				  0., //group H2O
				  0., //group ACOH
				  0.,1., //group ketone
				  0.,   //group aldehyde  
				  0.,0., //group ester
				  0.,0.,0., //group ether 
				  0.,  //group acid
				  0.,   //group ACNO2
				  0.,0.,0., //group NO3
				  0.,0.,0., //group CO-OH
  			          0.,0.,0.,0.,0.,0.,0.,0.,0., //group CO-OC
			          1.,  //group PAN
			          0.,  //group CO-OOH
			          0.,  //group O=COC=O
			          0.,0.,0.}; //group CHxNO2
  
  size = sizeof(group_tmp_PK5003)/sizeof(double);
  assert(size = 60);
	
  for(int i = 0; i < size; ++i)
	PK5003.groups[i] = group_tmp_PK5003[i];

  // Search the species name in the aerosol species list 
  // and add the species if its name matches with
  // the given list.
  add_species_ssh(surrogate, PK5003, species_list_aer, molecular_weight_aer,
              accomodation_coefficient,diffusion_coef,saturation_vapor_pressure,enthalpy_vaporization, henry, t_ref, mass_density, species_part,nlayer,i_hydrophilic,
		  N_inert, N_inorganic);



/* ==== PP4000 ==== */ 

  species PP4000;
  PP4000.name="PP4000";
  PP4000.is_inorganic_precursor=false;
  PP4000.Psat_ref=3.008e-7; 		// Saturation vapor pressure at Tref (torr)
  PP4000.Tref=298;         		// Temperature of reference (K)
  PP4000.deltaH=122.40;     		// Enthalpy of vaporization (kJ/mol)
  PP4000.Henry=0.0;     		// Henry's law constant at Tref (M/atm)
  PP4000.aq_type="none"; 		// "none","diacid","monoacid" or "aldehyde"
  // PP4000.Kacidity1=XXacidityXX;   		// First acidity constant
  PP4000.hydrophilic=true;   	// Does the species condense on the aqueous phase?
  PP4000.hydrophobic=true;  		// Does the species condense on the organic phase?
  PP4000.nonvolatile=false; 		// Is the compound nonvolatile?
  PP4000.kp_from_experiment=false; 	// Use experimental partitioning constant at Tref?
  PP4000.is_organic=true;  		// Is the compound organic?
  PP4000.compute_gamma_org=true;  	// Compute the activity coefficients of the organic phase for this compound?
  PP4000.compute_gamma_aq=true;  	// Compute the activity coefficients of the aqueous phase for this compound?
  PP4000.Koligo_org=0.0;      	//oligomeriation constant in the organic phase
  PP4000.rho=1300.;  
  PP4000.is_monomer=false;
  PP4000.rion=false;
  PP4000.KDiffusion_air=1.0e-5;
  //  PP4000.accomodation_coefficient=alpha;
  PP4000.viscosity=1.68e12;
  PP4000.is_solid=false;
  PP4000.is_generic=false;
  
  //Group: if no functionnal group in the species use the default species
  //for the computation of activity coefficients
  //
  double group_tmp_PP4000 [] = {0.,0.,0.,0., // group C
				  0.,0.,2.,0., //group C[OH]
				  0.,0.,0.,0., //group Calcohol
				  0.,0.,0.,0., //group Calcohol-tail
				  0.,0.,0.,0.,0., //group C=C
				  0.,0., //group aromatic carbon (AC)
				  0.,0.,0., // group //AC-C
				  2.,  //group OH
				  0., //group H2O
				  0., //group ACOH
				  0.,0., //group ketone
				  0.,   //group aldehyde  
				  0.,0., //group ester
				  0.,0.,0., //group ether 
				  0.,  //group acid
				  0.,   //group ACNO2
				  0.,0.,0., //group NO3
				  0.,0.,0., //group CO-OH
  			          0.,0.,0.,0.,0.,0.,0.,0.,0., //group CO-OC
			          2.,  //group PAN
			          0.,  //group CO-OOH
			          0.,  //group O=COC=O
			          0.,0.,0.}; //group CHxNO2
  
  size = sizeof(group_tmp_PP4000)/sizeof(double);
  assert(size = 60);
	
  for(int i = 0; i < size; ++i)
	PP4000.groups[i] = group_tmp_PP4000[i];

  // Search the species name in the aerosol species list 
  // and add the species if its name matches with
  // the given list.
  add_species_ssh(surrogate, PP4000, species_list_aer, molecular_weight_aer,
              accomodation_coefficient,diffusion_coef,saturation_vapor_pressure,enthalpy_vaporization, henry, t_ref, mass_density, species_part,nlayer,i_hydrophilic,
		  N_inert, N_inorganic);



/* ==== PP4004 ==== */ 

  species PP4004;
  PP4004.name="PP4004";
  PP4004.is_inorganic_precursor=false;
  PP4004.Psat_ref=9.221e-6; 		// Saturation vapor pressure at Tref (torr)
  PP4004.Tref=298;         		// Temperature of reference (K)
  PP4004.deltaH=105.11;     		// Enthalpy of vaporization (kJ/mol)
  PP4004.Henry=0.0;     		// Henry's law constant at Tref (M/atm)
  PP4004.aq_type="none"; 		// "none","diacid","monoacid" or "aldehyde"
  // PP4004.Kacidity1=XXacidityXX;   		// First acidity constant
  PP4004.hydrophilic=true;   	// Does the species condense on the aqueous phase?
  PP4004.hydrophobic=true;  		// Does the species condense on the organic phase?
  PP4004.nonvolatile=false; 		// Is the compound nonvolatile?
  PP4004.kp_from_experiment=false; 	// Use experimental partitioning constant at Tref?
  PP4004.is_organic=true;  		// Is the compound organic?
  PP4004.compute_gamma_org=true;  	// Compute the activity coefficients of the organic phase for this compound?
  PP4004.compute_gamma_aq=true;  	// Compute the activity coefficients of the aqueous phase for this compound?
  PP4004.Koligo_org=0.0;      	//oligomeriation constant in the organic phase
  PP4004.rho=1300.;  
  PP4004.is_monomer=false;
  PP4004.rion=false;
  PP4004.KDiffusion_air=1.0e-5;
  //  PP4004.accomodation_coefficient=alpha;
  PP4004.viscosity=1.68e12;
  PP4004.is_solid=false;
  PP4004.is_generic=false;
  
  //Group: if no functionnal group in the species use the default species
  //for the computation of activity coefficients
  //
  double group_tmp_PP4004 [] = {0.,0.,0.,0., // group C
				  0.,0.,1.,0., //group C[OH]
				  0.,0.,0.,0., //group Calcohol
				  0.,0.,0.,0., //group Calcohol-tail
				  0.,0.,0.,0.,0., //group C=C
				  0.,0., //group aromatic carbon (AC)
				  0.,0.,0., // group //AC-C
				  1.,  //group OH
				  0., //group H2O
				  0., //group ACOH
				  0.,1., //group ketone
				  0.,   //group aldehyde  
				  0.,0., //group ester
				  0.,0.,0., //group ether 
				  0.,  //group acid
				  0.,   //group ACNO2
				  0.,0.,0., //group NO3
				  0.,0.,0., //group CO-OH
  			          0.,0.,0.,0.,0.,0.,0.,0.,0., //group CO-OC
			          2.,  //group PAN
			          0.,  //group CO-OOH
			          0.,  //group O=COC=O
			          0.,0.,0.}; //group CHxNO2
  
  size = sizeof(group_tmp_PP4004)/sizeof(double);
  assert(size = 60);
	
  for(int i = 0; i < size; ++i)
	PP4004.groups[i] = group_tmp_PP4004[i];

  // Search the species name in the aerosol species list 
  // and add the species if its name matches with
  // the given list.
  add_species_ssh(surrogate, PP4004, species_list_aer, molecular_weight_aer,
              accomodation_coefficient,diffusion_coef,saturation_vapor_pressure,enthalpy_vaporization, henry, t_ref, mass_density, species_part,nlayer,i_hydrophilic,
		  N_inert, N_inorganic);



/* ==== PU5000 ==== */ 

  species PU5000;
  PU5000.name="PU5000";
  PU5000.is_inorganic_precursor=false;
  PU5000.Psat_ref=1.493e-1; 		// Saturation vapor pressure at Tref (torr)
  PU5000.Tref=298;         		// Temperature of reference (K)
  PU5000.deltaH=60.71;     		// Enthalpy of vaporization (kJ/mol)
  PU5000.Henry=0.0;     		// Henry's law constant at Tref (M/atm)
  PU5000.aq_type="none"; 		// "none","diacid","monoacid" or "aldehyde"
  // PU5000.Kacidity1=XXacidityXX;   		// First acidity constant
  PU5000.hydrophilic=true;   	// Does the species condense on the aqueous phase?
  PU5000.hydrophobic=true;  		// Does the species condense on the organic phase?
  PU5000.nonvolatile=false; 		// Is the compound nonvolatile?
  PU5000.kp_from_experiment=false; 	// Use experimental partitioning constant at Tref?
  PU5000.is_organic=true;  		// Is the compound organic?
  PU5000.compute_gamma_org=true;  	// Compute the activity coefficients of the organic phase for this compound?
  PU5000.compute_gamma_aq=true;  	// Compute the activity coefficients of the aqueous phase for this compound?
  PU5000.Koligo_org=0.0;      	//oligomeriation constant in the organic phase
  PU5000.rho=1300.;  
  PU5000.is_monomer=false;
  PU5000.rion=false;
  PU5000.KDiffusion_air=1.0e-5;
  //  PU5000.accomodation_coefficient=alpha;
  PU5000.viscosity=1.68e12;
  PU5000.is_solid=false;
  PU5000.is_generic=false;
  
  //Group: if no functionnal group in the species use the default species
  //for the computation of activity coefficients
  //
  double group_tmp_PU5000 [] = {1.,0.,0.,0., // group C
				  0.,0.,0.,0., //group C[OH]
				  0.,0.,0.,0., //group Calcohol
				  0.,0.,0.,0., //group Calcohol-tail
				  0.,0.,0.,1.,0., //group C=C
				  0.,0., //group aromatic carbon (AC)
				  0.,0.,0., // group //AC-C
				  0.,  //group OH
				  0., //group H2O
				  0., //group ACOH
				  0.,0., //group ketone
				  1.,   //group aldehyde  
				  0.,0., //group ester
				  0.,0.,0., //group ether 
				  0.,  //group acid
				  0.,   //group ACNO2
				  0.,0.,0., //group NO3
				  0.,0.,0., //group CO-OH
  			          0.,0.,0.,0.,0.,0.,0.,0.,0., //group CO-OC
			          1.,  //group PAN
			          0.,  //group CO-OOH
			          0.,  //group O=COC=O
			          0.,0.,0.}; //group CHxNO2
  
  size = sizeof(group_tmp_PU5000)/sizeof(double);
  assert(size = 60);
	
  for(int i = 0; i < size; ++i)
	PU5000.groups[i] = group_tmp_PU5000[i];

  // Search the species name in the aerosol species list 
  // and add the species if its name matches with
  // the given list.
  add_species_ssh(surrogate, PU5000, species_list_aer, molecular_weight_aer,
              accomodation_coefficient,diffusion_coef,saturation_vapor_pressure,enthalpy_vaporization, henry, t_ref, mass_density, species_part,nlayer,i_hydrophilic,
		  N_inert, N_inorganic);



/* ==== PU5001 ==== */ 

  species PU5001;
  PU5001.name="PU5001";
  PU5001.is_inorganic_precursor=false;
  PU5001.Psat_ref=1.493e-1; 		// Saturation vapor pressure at Tref (torr)
  PU5001.Tref=298;         		// Temperature of reference (K)
  PU5001.deltaH=60.71;     		// Enthalpy of vaporization (kJ/mol)
  PU5001.Henry=0.0;     		// Henry's law constant at Tref (M/atm)
  PU5001.aq_type="none"; 		// "none","diacid","monoacid" or "aldehyde"
  // PU5001.Kacidity1=XXacidityXX;   		// First acidity constant
  PU5001.hydrophilic=true;   	// Does the species condense on the aqueous phase?
  PU5001.hydrophobic=true;  		// Does the species condense on the organic phase?
  PU5001.nonvolatile=false; 		// Is the compound nonvolatile?
  PU5001.kp_from_experiment=false; 	// Use experimental partitioning constant at Tref?
  PU5001.is_organic=true;  		// Is the compound organic?
  PU5001.compute_gamma_org=true;  	// Compute the activity coefficients of the organic phase for this compound?
  PU5001.compute_gamma_aq=true;  	// Compute the activity coefficients of the aqueous phase for this compound?
  PU5001.Koligo_org=0.0;      	//oligomeriation constant in the organic phase
  PU5001.rho=1300.;  
  PU5001.is_monomer=false;
  PU5001.rion=false;
  PU5001.KDiffusion_air=1.0e-5;
  //  PU5001.accomodation_coefficient=alpha;
  PU5001.viscosity=1.68e12;
  PU5001.is_solid=false;
  PU5001.is_generic=false;
  
  //Group: if no functionnal group in the species use the default species
  //for the computation of activity coefficients
  //
  double group_tmp_PU5001 [] = {1.,0.,0.,0., // group C
				  0.,0.,0.,0., //group C[OH]
				  0.,0.,0.,0., //group Calcohol
				  0.,0.,0.,0., //group Calcohol-tail
				  0.,0.,0.,1.,0., //group C=C
				  0.,0., //group aromatic carbon (AC)
				  0.,0.,0., // group //AC-C
				  0.,  //group OH
				  0., //group H2O
				  0., //group ACOH
				  0.,0., //group ketone
				  1.,   //group aldehyde  
				  0.,0., //group ester
				  0.,0.,0., //group ether 
				  0.,  //group acid
				  0.,   //group ACNO2
				  0.,0.,0., //group NO3
				  0.,0.,0., //group CO-OH
  			          0.,0.,0.,0.,0.,0.,0.,0.,0., //group CO-OC
			          1.,  //group PAN
			          0.,  //group CO-OOH
			          0.,  //group O=COC=O
			          0.,0.,0.}; //group CHxNO2
  
  size = sizeof(group_tmp_PU5001)/sizeof(double);
  assert(size = 60);
	
  for(int i = 0; i < size; ++i)
	PU5001.groups[i] = group_tmp_PU5001[i];

  // Search the species name in the aerosol species list 
  // and add the species if its name matches with
  // the given list.
  add_species_ssh(surrogate, PU5001, species_list_aer, molecular_weight_aer,
              accomodation_coefficient,diffusion_coef,saturation_vapor_pressure,enthalpy_vaporization, henry, t_ref, mass_density, species_part,nlayer,i_hydrophilic,
		  N_inert, N_inorganic);



/* ==== PU5002 ==== */ 

  species PU5002;
  PU5002.name="PU5002";
  PU5002.is_inorganic_precursor=false;
  PU5002.Psat_ref=9.097e-2; 		// Saturation vapor pressure at Tref (torr)
  PU5002.Tref=298;         		// Temperature of reference (K)
  PU5002.deltaH=62.59;     		// Enthalpy of vaporization (kJ/mol)
  PU5002.Henry=0.0;     		// Henry's law constant at Tref (M/atm)
  PU5002.aq_type="none"; 		// "none","diacid","monoacid" or "aldehyde"
  // PU5002.Kacidity1=XXacidityXX;   		// First acidity constant
  PU5002.hydrophilic=true;   	// Does the species condense on the aqueous phase?
  PU5002.hydrophobic=true;  		// Does the species condense on the organic phase?
  PU5002.nonvolatile=false; 		// Is the compound nonvolatile?
  PU5002.kp_from_experiment=false; 	// Use experimental partitioning constant at Tref?
  PU5002.is_organic=true;  		// Is the compound organic?
  PU5002.compute_gamma_org=true;  	// Compute the activity coefficients of the organic phase for this compound?
  PU5002.compute_gamma_aq=true;  	// Compute the activity coefficients of the aqueous phase for this compound?
  PU5002.Koligo_org=0.0;      	//oligomeriation constant in the organic phase
  PU5002.rho=1300.;  
  PU5002.is_monomer=false;
  PU5002.rion=false;
  PU5002.KDiffusion_air=1.0e-5;
  //  PU5002.accomodation_coefficient=alpha;
  PU5002.viscosity=1.68e12;
  PU5002.is_solid=false;
  PU5002.is_generic=false;
  
  //Group: if no functionnal group in the species use the default species
  //for the computation of activity coefficients
  //
  double group_tmp_PU5002 [] = {0.,0.,0.,0., // group C
				  0.,0.,0.,0., //group C[OH]
				  0.,0.,0.,0., //group Calcohol
				  0.,0.,0.,0., //group Calcohol-tail
				  0.,1.,0.,0.,0., //group C=C
				  0.,0., //group aromatic carbon (AC)
				  0.,0.,0., // group //AC-C
				  0.,  //group OH
				  0., //group H2O
				  0., //group ACOH
				  1.,0., //group ketone
				  0.,   //group aldehyde  
				  0.,0., //group ester
				  0.,0.,0., //group ether 
				  0.,  //group acid
				  0.,   //group ACNO2
				  0.,0.,0., //group NO3
				  0.,0.,0., //group CO-OH
  			          0.,0.,0.,0.,0.,0.,0.,0.,0., //group CO-OC
			          1.,  //group PAN
			          0.,  //group CO-OOH
			          0.,  //group O=COC=O
			          0.,0.,0.}; //group CHxNO2
  
  size = sizeof(group_tmp_PU5002)/sizeof(double);
  assert(size = 60);
	
  for(int i = 0; i < size; ++i)
	PU5002.groups[i] = group_tmp_PU5002[i];

  // Search the species name in the aerosol species list 
  // and add the species if its name matches with
  // the given list.
  add_species_ssh(surrogate, PU5002, species_list_aer, molecular_weight_aer,
              accomodation_coefficient,diffusion_coef,saturation_vapor_pressure,enthalpy_vaporization, henry, t_ref, mass_density, species_part,nlayer,i_hydrophilic,
		  N_inert, N_inorganic);



/* ==== TOL2OHOOH ==== */ 

  species TOL2OHOOH;
  TOL2OHOOH.name="TOL2OHOOH";
  TOL2OHOOH.is_inorganic_precursor=false;
  TOL2OHOOH.Psat_ref=1.204e-6; 		// Saturation vapor pressure at Tref (torr)
  TOL2OHOOH.Tref=298;         		// Temperature of reference (K)
  TOL2OHOOH.deltaH=115.83;     		// Enthalpy of vaporization (kJ/mol)
  TOL2OHOOH.Henry=0.0;     		// Henry's law constant at Tref (M/atm)
  TOL2OHOOH.aq_type="none"; 		// "none","diacid","monoacid" or "aldehyde"
  // TOL2OHOOH.Kacidity1=XXacidityXX;   		// First acidity constant
  TOL2OHOOH.hydrophilic=true;   	// Does the species condense on the aqueous phase?
  TOL2OHOOH.hydrophobic=true;  		// Does the species condense on the organic phase?
  TOL2OHOOH.nonvolatile=false; 		// Is the compound nonvolatile?
  TOL2OHOOH.kp_from_experiment=false; 	// Use experimental partitioning constant at Tref?
  TOL2OHOOH.is_organic=true;  		// Is the compound organic?
  TOL2OHOOH.compute_gamma_org=true;  	// Compute the activity coefficients of the organic phase for this compound?
  TOL2OHOOH.compute_gamma_aq=true;  	// Compute the activity coefficients of the aqueous phase for this compound?
  TOL2OHOOH.Koligo_org=0.0;      	//oligomeriation constant in the organic phase
  TOL2OHOOH.rho=1300.;  
  TOL2OHOOH.is_monomer=false;
  TOL2OHOOH.rion=false;
  TOL2OHOOH.KDiffusion_air=1.0e-5;
  //  TOL2OHOOH.accomodation_coefficient=alpha;
  TOL2OHOOH.viscosity=1.68e12;
  TOL2OHOOH.is_solid=false;
  TOL2OHOOH.is_generic=false;
  
  //Group: if no functionnal group in the species use the default species
  //for the computation of activity coefficients
  //
  double group_tmp_TOL2OHOOH [] = {0.,0.,0.,0., // group C
				  0.,0.,0.,0., //group C[OH]
				  0.,0.,0.,0., //group Calcohol
				  0.,0.,0.,0., //group Calcohol-tail
				  0.,0.,0.,0.,0., //group C=C
				  3.,1., //group aromatic carbon (AC)
				  0.,0.,0., // group //AC-C
				  0.,  //group OH
				  0., //group H2O
				  2., //group ACOH
				  0.,0., //group ketone
				  0.,   //group aldehyde  
				  0.,0., //group ester
				  0.,0.,0., //group ether 
				  0.,  //group acid
				  0.,   //group ACNO2
				  0.,0.,0., //group NO3
				  1.,0.,0., //group CO-OH
  			          0.,0.,0.,0.,0.,0.,0.,0.,0., //group CO-OC
			          0.,  //group PAN
			          0.,  //group CO-OOH
			          0.,  //group O=COC=O
			          0.,0.,0.}; //group CHxNO2
  
  size = sizeof(group_tmp_TOL2OHOOH)/sizeof(double);
  assert(size = 60);
	
  for(int i = 0; i < size; ++i)
	TOL2OHOOH.groups[i] = group_tmp_TOL2OHOOH[i];

  // Search the species name in the aerosol species list 
  // and add the species if its name matches with
  // the given list.
  add_species_ssh(surrogate, TOL2OHOOH, species_list_aer, molecular_weight_aer,
              accomodation_coefficient,diffusion_coef,saturation_vapor_pressure,enthalpy_vaporization, henry, t_ref, mass_density, species_part,nlayer,i_hydrophilic,
		  N_inert, N_inorganic);



/* ==== TOL3OH1NO2 ==== */ 

  species TOL3OH1NO2;
  TOL3OH1NO2.name="TOL3OH1NO2";
  TOL3OH1NO2.is_inorganic_precursor=false;
  TOL3OH1NO2.Psat_ref=3.0256e-7; 		// Saturation vapor pressure at Tref (torr)
  TOL3OH1NO2.Tref=298;         		// Temperature of reference (K)
  TOL3OH1NO2.deltaH=122.37;     		// Enthalpy of vaporization (kJ/mol)
  TOL3OH1NO2.Henry=0.0;     		// Henry's law constant at Tref (M/atm)
  TOL3OH1NO2.aq_type="none"; 		// "none","diacid","monoacid" or "aldehyde"
  // TOL3OH1NO2.Kacidity1=XXacidityXX;   		// First acidity constant
  TOL3OH1NO2.hydrophilic=true;   	// Does the species condense on the aqueous phase?
  TOL3OH1NO2.hydrophobic=true;  		// Does the species condense on the organic phase?
  TOL3OH1NO2.nonvolatile=false; 		// Is the compound nonvolatile?
  TOL3OH1NO2.kp_from_experiment=false; 	// Use experimental partitioning constant at Tref?
  TOL3OH1NO2.is_organic=true;  		// Is the compound organic?
  TOL3OH1NO2.compute_gamma_org=true;  	// Compute the activity coefficients of the organic phase for this compound?
  TOL3OH1NO2.compute_gamma_aq=true;  	// Compute the activity coefficients of the aqueous phase for this compound?
  TOL3OH1NO2.Koligo_org=0.0;      	//oligomeriation constant in the organic phase
  TOL3OH1NO2.rho=1300.;  
  TOL3OH1NO2.is_monomer=false;
  TOL3OH1NO2.rion=false;
  TOL3OH1NO2.KDiffusion_air=1.0e-5;
  //  TOL3OH1NO2.accomodation_coefficient=alpha;
  TOL3OH1NO2.viscosity=1.68e12;
  TOL3OH1NO2.is_solid=false;
  TOL3OH1NO2.is_generic=false;
  
  //Group: if no functionnal group in the species use the default species
  //for the computation of activity coefficients
  //
  double group_tmp_TOL3OH1NO2 [] = {0.,0.,0.,0., // group C
				  0.,0.,0.,0., //group C[OH]
				  0.,0.,0.,0., //group Calcohol
				  0.,0.,0.,0., //group Calcohol-tail
				  0.,0.,0.,0.,0., //group C=C
				  1.,0., //group aromatic carbon (AC)
				  1.,0.,0., // group //AC-C
				  0.,  //group OH
				  0., //group H2O
				  3., //group ACOH
				  0.,0., //group ketone
				  0.,   //group aldehyde  
				  0.,0., //group ester
				  0.,0.,0., //group ether 
				  0.,  //group acid
				  1.,   //group ACNO2
				  0.,0.,0., //group NO3
				  0.,0.,0., //group CO-OH
  			          0.,0.,0.,0.,0.,0.,0.,0.,0., //group CO-OC
			          0.,  //group PAN
			          0.,  //group CO-OOH
			          0.,  //group O=COC=O
			          0.,0.,0.}; //group CHxNO2
  
  size = sizeof(group_tmp_TOL3OH1NO2)/sizeof(double);
  assert(size = 60);
	
  for(int i = 0; i < size; ++i)
	TOL3OH1NO2.groups[i] = group_tmp_TOL3OH1NO2[i];

  // Search the species name in the aerosol species list 
  // and add the species if its name matches with
  // the given list.
  add_species_ssh(surrogate, TOL3OH1NO2, species_list_aer, molecular_weight_aer,
              accomodation_coefficient,diffusion_coef,saturation_vapor_pressure,enthalpy_vaporization, henry, t_ref, mass_density, species_part,nlayer,i_hydrophilic,
		  N_inert, N_inorganic);



/* ==== TOL3OH ==== */ 

  species TOL3OH;
  TOL3OH.name="TOL3OH";
  TOL3OH.is_inorganic_precursor=false;
  TOL3OH.Psat_ref=9.5678e-5; 		// Saturation vapor pressure at Tref (torr)
  TOL3OH.Tref=298;         		// Temperature of reference (K)
  TOL3OH.deltaH=95.14;     		// Enthalpy of vaporization (kJ/mol)
  TOL3OH.Henry=0.0;     		// Henry's law constant at Tref (M/atm)
  TOL3OH.aq_type="none"; 		// "none","diacid","monoacid" or "aldehyde"
  // TOL3OH.Kacidity1=XXacidityXX;   		// First acidity constant
  TOL3OH.hydrophilic=true;   	// Does the species condense on the aqueous phase?
  TOL3OH.hydrophobic=true;  		// Does the species condense on the organic phase?
  TOL3OH.nonvolatile=false; 		// Is the compound nonvolatile?
  TOL3OH.kp_from_experiment=false; 	// Use experimental partitioning constant at Tref?
  TOL3OH.is_organic=true;  		// Is the compound organic?
  TOL3OH.compute_gamma_org=true;  	// Compute the activity coefficients of the organic phase for this compound?
  TOL3OH.compute_gamma_aq=true;  	// Compute the activity coefficients of the aqueous phase for this compound?
  TOL3OH.Koligo_org=0.0;      	//oligomeriation constant in the organic phase
  TOL3OH.rho=1300.;  
  TOL3OH.is_monomer=false;
  TOL3OH.rion=false;
  TOL3OH.KDiffusion_air=1.0e-5;
  //  TOL3OH.accomodation_coefficient=alpha;
  TOL3OH.viscosity=1.68e12;
  TOL3OH.is_solid=false;
  TOL3OH.is_generic=false;
  
  //Group: if no functionnal group in the species use the default species
  //for the computation of activity coefficients
  //
  double group_tmp_TOL3OH [] = {0.,0.,0.,0., // group C
				  0.,0.,0.,0., //group C[OH]
				  0.,0.,0.,0., //group Calcohol
				  0.,0.,0.,0., //group Calcohol-tail
				  0.,0.,0.,0.,0., //group C=C
				  2.,0., //group aromatic carbon (AC)
				  1.,0.,0., // group //AC-C
				  0.,  //group OH
				  0., //group H2O
				  3., //group ACOH
				  0.,0., //group ketone
				  0.,   //group aldehyde  
				  0.,0., //group ester
				  0.,0.,0., //group ether 
				  0.,  //group acid
				  0.,   //group ACNO2
				  0.,0.,0., //group NO3
				  0.,0.,0., //group CO-OH
  			          0.,0.,0.,0.,0.,0.,0.,0.,0., //group CO-OC
			          0.,  //group PAN
			          0.,  //group CO-OOH
			          0.,  //group O=COC=O
			          0.,0.,0.}; //group CHxNO2
  
  size = sizeof(group_tmp_TOL3OH)/sizeof(double);
  assert(size = 60);
	
  for(int i = 0; i < size; ++i)
	TOL3OH.groups[i] = group_tmp_TOL3OH[i];

  // Search the species name in the aerosol species list 
  // and add the species if its name matches with
  // the given list.
  add_species_ssh(surrogate, TOL3OH, species_list_aer, molecular_weight_aer,
              accomodation_coefficient,diffusion_coef,saturation_vapor_pressure,enthalpy_vaporization, henry, t_ref, mass_density, species_part,nlayer,i_hydrophilic,
		  N_inert, N_inorganic);



/* ==== TOL3OHOOH ==== */ 

  species TOL3OHOOH;
  TOL3OHOOH.name="TOL3OHOOH";
  TOL3OHOOH.is_inorganic_precursor=false;
  TOL3OHOOH.Psat_ref=6.036e-9; 		// Saturation vapor pressure at Tref (torr)
  TOL3OHOOH.Tref=298;         		// Temperature of reference (K)
  TOL3OHOOH.deltaH=140.88;     		// Enthalpy of vaporization (kJ/mol)
  TOL3OHOOH.Henry=0.0;     		// Henry's law constant at Tref (M/atm)
  TOL3OHOOH.aq_type="none"; 		// "none","diacid","monoacid" or "aldehyde"
  // TOL3OHOOH.Kacidity1=XXacidityXX;   		// First acidity constant
  TOL3OHOOH.hydrophilic=true;   	// Does the species condense on the aqueous phase?
  TOL3OHOOH.hydrophobic=true;  		// Does the species condense on the organic phase?
  TOL3OHOOH.nonvolatile=false; 		// Is the compound nonvolatile?
  TOL3OHOOH.kp_from_experiment=false; 	// Use experimental partitioning constant at Tref?
  TOL3OHOOH.is_organic=true;  		// Is the compound organic?
  TOL3OHOOH.compute_gamma_org=true;  	// Compute the activity coefficients of the organic phase for this compound?
  TOL3OHOOH.compute_gamma_aq=true;  	// Compute the activity coefficients of the aqueous phase for this compound?
  TOL3OHOOH.Koligo_org=0.0;      	//oligomeriation constant in the organic phase
  TOL3OHOOH.rho=1300.;  
  TOL3OHOOH.is_monomer=false;
  TOL3OHOOH.rion=false;
  TOL3OHOOH.KDiffusion_air=1.0e-5;
  //  TOL3OHOOH.accomodation_coefficient=alpha;
  TOL3OHOOH.viscosity=1.68e12;
  TOL3OHOOH.is_solid=false;
  TOL3OHOOH.is_generic=false;
  
  //Group: if no functionnal group in the species use the default species
  //for the computation of activity coefficients
  //
  double group_tmp_TOL3OHOOH [] = {0.,0.,0.,0., // group C
				  0.,0.,0.,0., //group C[OH]
				  0.,0.,0.,0., //group Calcohol
				  0.,0.,0.,0., //group Calcohol-tail
				  0.,0.,0.,0.,0., //group C=C
				  2.,1., //group aromatic carbon (AC)
				  0.,0.,0., // group //AC-C
				  0.,  //group OH
				  0., //group H2O
				  3., //group ACOH
				  0.,0., //group ketone
				  0.,   //group aldehyde  
				  0.,0., //group ester
				  0.,0.,0., //group ether 
				  0.,  //group acid
				  0.,   //group ACNO2
				  0.,0.,0., //group NO3
				  1.,0.,0., //group CO-OH
  			          0.,0.,0.,0.,0.,0.,0.,0.,0., //group CO-OC
			          0.,  //group PAN
			          0.,  //group CO-OOH
			          0.,  //group O=COC=O
			          0.,0.,0.}; //group CHxNO2
  
  size = sizeof(group_tmp_TOL3OHOOH)/sizeof(double);
  assert(size = 60);
	
  for(int i = 0; i < size; ++i)
	TOL3OHOOH.groups[i] = group_tmp_TOL3OHOOH[i];

  // Search the species name in the aerosol species list 
  // and add the species if its name matches with
  // the given list.
  add_species_ssh(surrogate, TOL3OHOOH, species_list_aer, molecular_weight_aer,
              accomodation_coefficient,diffusion_coef,saturation_vapor_pressure,enthalpy_vaporization, henry, t_ref, mass_density, species_part,nlayer,i_hydrophilic,
		  N_inert, N_inorganic);



/* ==== TOL4OH1NO2 ==== */ 

  species TOL4OH1NO2;
  TOL4OH1NO2.name="TOL4OH1NO2";
  TOL4OH1NO2.is_inorganic_precursor=false;
  TOL4OH1NO2.Psat_ref=1.5164e-9; 		// Saturation vapor pressure at Tref (torr)
  TOL4OH1NO2.Tref=298;         		// Temperature of reference (K)
  TOL4OH1NO2.deltaH=147.41;     		// Enthalpy of vaporization (kJ/mol)
  TOL4OH1NO2.Henry=0.0;     		// Henry's law constant at Tref (M/atm)
  TOL4OH1NO2.aq_type="none"; 		// "none","diacid","monoacid" or "aldehyde"
  // TOL4OH1NO2.Kacidity1=XXacidityXX;   		// First acidity constant
  TOL4OH1NO2.hydrophilic=true;   	// Does the species condense on the aqueous phase?
  TOL4OH1NO2.hydrophobic=true;  		// Does the species condense on the organic phase?
  TOL4OH1NO2.nonvolatile=false; 		// Is the compound nonvolatile?
  TOL4OH1NO2.kp_from_experiment=false; 	// Use experimental partitioning constant at Tref?
  TOL4OH1NO2.is_organic=true;  		// Is the compound organic?
  TOL4OH1NO2.compute_gamma_org=true;  	// Compute the activity coefficients of the organic phase for this compound?
  TOL4OH1NO2.compute_gamma_aq=true;  	// Compute the activity coefficients of the aqueous phase for this compound?
  TOL4OH1NO2.Koligo_org=0.0;      	//oligomeriation constant in the organic phase
  TOL4OH1NO2.rho=1300.;  
  TOL4OH1NO2.is_monomer=false;
  TOL4OH1NO2.rion=false;
  TOL4OH1NO2.KDiffusion_air=1.0e-5;
  //  TOL4OH1NO2.accomodation_coefficient=alpha;
  TOL4OH1NO2.viscosity=1.68e12;
  TOL4OH1NO2.is_solid=false;
  TOL4OH1NO2.is_generic=false;
  
  //Group: if no functionnal group in the species use the default species
  //for the computation of activity coefficients
  //
  double group_tmp_TOL4OH1NO2 [] = {0.,0.,0.,0., // group C
				  0.,0.,0.,0., //group C[OH]
				  0.,0.,0.,0., //group Calcohol
				  0.,0.,0.,0., //group Calcohol-tail
				  0.,0.,0.,0.,0., //group C=C
				  0.,0., //group aromatic carbon (AC)
				  1.,0.,0., // group //AC-C
				  0.,  //group OH
				  0., //group H2O
				  4., //group ACOH
				  0.,0., //group ketone
				  0.,   //group aldehyde  
				  0.,0., //group ester
				  0.,0.,0., //group ether 
				  0.,  //group acid
				  1.,   //group ACNO2
				  0.,0.,0., //group NO3
				  0.,0.,0., //group CO-OH
  			          0.,0.,0.,0.,0.,0.,0.,0.,0., //group CO-OC
			          0.,  //group PAN
			          0.,  //group CO-OOH
			          0.,  //group O=COC=O
			          0.,0.,0.}; //group CHxNO2
  
  size = sizeof(group_tmp_TOL4OH1NO2)/sizeof(double);
  assert(size = 60);
	
  for(int i = 0; i < size; ++i)
	TOL4OH1NO2.groups[i] = group_tmp_TOL4OH1NO2[i];

  // Search the species name in the aerosol species list 
  // and add the species if its name matches with
  // the given list.
  add_species_ssh(surrogate, TOL4OH1NO2, species_list_aer, molecular_weight_aer,
              accomodation_coefficient,diffusion_coef,saturation_vapor_pressure,enthalpy_vaporization, henry, t_ref, mass_density, species_part,nlayer,i_hydrophilic,
		  N_inert, N_inorganic);



/* ==== TOL4OH ==== */ 

  species TOL4OH;
  TOL4OH.name="TOL4OH";
  TOL4OH.is_inorganic_precursor=false;
  TOL4OH.Psat_ref=6.0368e-7; 		// Saturation vapor pressure at Tref (torr)
  TOL4OH.Tref=298;         		// Temperature of reference (K)
  TOL4OH.deltaH=119.10;     		// Enthalpy of vaporization (kJ/mol)
  TOL4OH.Henry=0.0;     		// Henry's law constant at Tref (M/atm)
  TOL4OH.aq_type="none"; 		// "none","diacid","monoacid" or "aldehyde"
  // TOL4OH.Kacidity1=XXacidityXX;   		// First acidity constant
  TOL4OH.hydrophilic=true;   	// Does the species condense on the aqueous phase?
  TOL4OH.hydrophobic=true;  		// Does the species condense on the organic phase?
  TOL4OH.nonvolatile=false; 		// Is the compound nonvolatile?
  TOL4OH.kp_from_experiment=false; 	// Use experimental partitioning constant at Tref?
  TOL4OH.is_organic=true;  		// Is the compound organic?
  TOL4OH.compute_gamma_org=true;  	// Compute the activity coefficients of the organic phase for this compound?
  TOL4OH.compute_gamma_aq=true;  	// Compute the activity coefficients of the aqueous phase for this compound?
  TOL4OH.Koligo_org=0.0;      	//oligomeriation constant in the organic phase
  TOL4OH.rho=1300.;  
  TOL4OH.is_monomer=false;
  TOL4OH.rion=false;
  TOL4OH.KDiffusion_air=1.0e-5;
  //  TOL4OH.accomodation_coefficient=alpha;
  TOL4OH.viscosity=1.68e12;
  TOL4OH.is_solid=false;
  TOL4OH.is_generic=false;
  
  //Group: if no functionnal group in the species use the default species
  //for the computation of activity coefficients
  //
  double group_tmp_TOL4OH [] = {0.,0.,0.,0., // group C
				  0.,0.,0.,0., //group C[OH]
				  0.,0.,0.,0., //group Calcohol
				  0.,0.,0.,0., //group Calcohol-tail
				  0.,0.,0.,0.,0., //group C=C
				  1.,0., //group aromatic carbon (AC)
				  1.,0.,0., // group //AC-C
				  0.,  //group OH
				  0., //group H2O
				  4., //group ACOH
				  0.,0., //group ketone
				  0.,   //group aldehyde  
				  0.,0., //group ester
				  0.,0.,0., //group ether 
				  0.,  //group acid
				  0.,   //group ACNO2
				  0.,0.,0., //group NO3
				  0.,0.,0., //group CO-OH
  			          0.,0.,0.,0.,0.,0.,0.,0.,0., //group CO-OC
			          0.,  //group PAN
			          0.,  //group CO-OOH
			          0.,  //group O=COC=O
			          0.,0.,0.}; //group CHxNO2
  
  size = sizeof(group_tmp_TOL4OH)/sizeof(double);
  assert(size = 60);
	
  for(int i = 0; i < size; ++i)
	TOL4OH.groups[i] = group_tmp_TOL4OH[i];

  // Search the species name in the aerosol species list 
  // and add the species if its name matches with
  // the given list.
  add_species_ssh(surrogate, TOL4OH, species_list_aer, molecular_weight_aer,
              accomodation_coefficient,diffusion_coef,saturation_vapor_pressure,enthalpy_vaporization, henry, t_ref, mass_density, species_part,nlayer,i_hydrophilic,
		  N_inert, N_inorganic);



/* ==== TOL5OH ==== */ 

  species TOL5OH;
  TOL5OH.name="TOL5OH";
  TOL5OH.is_inorganic_precursor=false;
  TOL5OH.Psat_ref=3.0256e-9; 		// Saturation vapor pressure at Tref (torr)
  TOL5OH.Tref=298;         		// Temperature of reference (K)
  TOL5OH.deltaH=144.15;     		// Enthalpy of vaporization (kJ/mol)
  TOL5OH.Henry=0.0;     		// Henry's law constant at Tref (M/atm)
  TOL5OH.aq_type="none"; 		// "none","diacid","monoacid" or "aldehyde"
  // TOL5OH.Kacidity1=XXacidityXX;   		// First acidity constant
  TOL5OH.hydrophilic=true;   	// Does the species condense on the aqueous phase?
  TOL5OH.hydrophobic=true;  		// Does the species condense on the organic phase?
  TOL5OH.nonvolatile=false; 		// Is the compound nonvolatile?
  TOL5OH.kp_from_experiment=false; 	// Use experimental partitioning constant at Tref?
  TOL5OH.is_organic=true;  		// Is the compound organic?
  TOL5OH.compute_gamma_org=true;  	// Compute the activity coefficients of the organic phase for this compound?
  TOL5OH.compute_gamma_aq=true;  	// Compute the activity coefficients of the aqueous phase for this compound?
  TOL5OH.Koligo_org=0.0;      	//oligomeriation constant in the organic phase
  TOL5OH.rho=1300.;  
  TOL5OH.is_monomer=false;
  TOL5OH.rion=false;
  TOL5OH.KDiffusion_air=1.0e-5;
  //  TOL5OH.accomodation_coefficient=alpha;
  TOL5OH.viscosity=1.68e12;
  TOL5OH.is_solid=false;
  TOL5OH.is_generic=false;
  
  //Group: if no functionnal group in the species use the default species
  //for the computation of activity coefficients
  //
  double group_tmp_TOL5OH [] = {0.,0.,0.,0., // group C
				  0.,0.,0.,0., //group C[OH]
				  0.,0.,0.,0., //group Calcohol
				  0.,0.,0.,0., //group Calcohol-tail
				  0.,0.,0.,0.,0., //group C=C
				  0.,0., //group aromatic carbon (AC)
				  1.,0.,0., // group //AC-C
				  0.,  //group OH
				  0., //group H2O
				  5., //group ACOH
				  0.,0., //group ketone
				  0.,   //group aldehyde  
				  0.,0., //group ester
				  0.,0.,0., //group ether 
				  0.,  //group acid
				  0.,   //group ACNO2
				  0.,0.,0., //group NO3
				  0.,0.,0., //group CO-OH
  			          0.,0.,0.,0.,0.,0.,0.,0.,0., //group CO-OC
			          0.,  //group PAN
			          0.,  //group CO-OOH
			          0.,  //group O=COC=O
			          0.,0.,0.}; //group CHxNO2
  
  size = sizeof(group_tmp_TOL5OH)/sizeof(double);
  assert(size = 60);
	
  for(int i = 0; i < size; ++i)
	TOL5OH.groups[i] = group_tmp_TOL5OH[i];

  // Search the species name in the aerosol species list 
  // and add the species if its name matches with
  // the given list.
  add_species_ssh(surrogate, TOL5OH, species_list_aer, molecular_weight_aer,
              accomodation_coefficient,diffusion_coef,saturation_vapor_pressure,enthalpy_vaporization, henry, t_ref, mass_density, species_part,nlayer,i_hydrophilic,
		  N_inert, N_inorganic);



/* ==== UD4000 ==== */ 

  species UD4000;
  UD4000.name="UD4000";
  UD4000.is_inorganic_precursor=false;
  UD4000.Psat_ref=6.317e-1; 		// Saturation vapor pressure at Tref (torr)
  UD4000.Tref=298;         		// Temperature of reference (K)
  UD4000.deltaH=55.61;     		// Enthalpy of vaporization (kJ/mol)
  UD4000.Henry=0.0;     		// Henry's law constant at Tref (M/atm)
  UD4000.aq_type="none"; 		// "none","diacid","monoacid" or "aldehyde"
  // UD4000.Kacidity1=XXacidityXX;   		// First acidity constant
  UD4000.hydrophilic=true;   	// Does the species condense on the aqueous phase?
  UD4000.hydrophobic=true;  		// Does the species condense on the organic phase?
  UD4000.nonvolatile=false; 		// Is the compound nonvolatile?
  UD4000.kp_from_experiment=false; 	// Use experimental partitioning constant at Tref?
  UD4000.is_organic=true;  		// Is the compound organic?
  UD4000.compute_gamma_org=true;  	// Compute the activity coefficients of the organic phase for this compound?
  UD4000.compute_gamma_aq=true;  	// Compute the activity coefficients of the aqueous phase for this compound?
  UD4000.Koligo_org=0.0;      	//oligomeriation constant in the organic phase
  UD4000.rho=1300.;  
  UD4000.is_monomer=false;
  UD4000.rion=false;
  UD4000.KDiffusion_air=1.0e-5;
  //  UD4000.accomodation_coefficient=alpha;
  UD4000.viscosity=1.68e12;
  UD4000.is_solid=false;
  UD4000.is_generic=false;
  
  //Group: if no functionnal group in the species use the default species
  //for the computation of activity coefficients
  //
  double group_tmp_UD4000 [] = {0.,0.,0.,0., // group C
				  0.,0.,0.,0., //group C[OH]
				  0.,0.,0.,0., //group Calcohol
				  0.,0.,0.,0., //group Calcohol-tail
				  0.,1.,0.,0.,0., //group C=C
				  0.,0., //group aromatic carbon (AC)
				  0.,0.,0., // group //AC-C
				  0.,  //group OH
				  0., //group H2O
				  0., //group ACOH
				  0.,0., //group ketone
				  2.,   //group aldehyde  
				  0.,0., //group ester
				  0.,0.,0., //group ether 
				  0.,  //group acid
				  0.,   //group ACNO2
				  0.,0.,0., //group NO3
				  0.,0.,0., //group CO-OH
  			          0.,0.,0.,0.,0.,0.,0.,0.,0., //group CO-OC
			          0.,  //group PAN
			          0.,  //group CO-OOH
			          0.,  //group O=COC=O
			          0.,0.,0.}; //group CHxNO2
  
  size = sizeof(group_tmp_UD4000)/sizeof(double);
  assert(size = 60);
	
  for(int i = 0; i < size; ++i)
	UD4000.groups[i] = group_tmp_UD4000[i];

  // Search the species name in the aerosol species list 
  // and add the species if its name matches with
  // the given list.
  add_species_ssh(surrogate, UD4000, species_list_aer, molecular_weight_aer,
              accomodation_coefficient,diffusion_coef,saturation_vapor_pressure,enthalpy_vaporization, henry, t_ref, mass_density, species_part,nlayer,i_hydrophilic,
		  N_inert, N_inorganic);



/* ==== UD5000 ==== */ 

  species UD5000;
  UD5000.name="UD5000";
  UD5000.is_inorganic_precursor=false;
  UD5000.Psat_ref=5.127; 		// Saturation vapor pressure at Tref (torr)
  UD5000.Tref=298;         		// Temperature of reference (K)
  UD5000.deltaH=46.24;     		// Enthalpy of vaporization (kJ/mol)
  UD5000.Henry=0.0;     		// Henry's law constant at Tref (M/atm)
  UD5000.aq_type="none"; 		// "none","diacid","monoacid" or "aldehyde"
  // UD5000.Kacidity1=XXacidityXX;   		// First acidity constant
  UD5000.hydrophilic=true;   	// Does the species condense on the aqueous phase?
  UD5000.hydrophobic=true;  		// Does the species condense on the organic phase?
  UD5000.nonvolatile=false; 		// Is the compound nonvolatile?
  UD5000.kp_from_experiment=false; 	// Use experimental partitioning constant at Tref?
  UD5000.is_organic=true;  		// Is the compound organic?
  UD5000.compute_gamma_org=true;  	// Compute the activity coefficients of the organic phase for this compound?
  UD5000.compute_gamma_aq=true;  	// Compute the activity coefficients of the aqueous phase for this compound?
  UD5000.Koligo_org=0.0;      	//oligomeriation constant in the organic phase
  UD5000.rho=1300.;  
  UD5000.is_monomer=false;
  UD5000.rion=false;
  UD5000.KDiffusion_air=1.0e-5;
  //  UD5000.accomodation_coefficient=alpha;
  UD5000.viscosity=1.68e12;
  UD5000.is_solid=false;
  UD5000.is_generic=false;
  
  //Group: if no functionnal group in the species use the default species
  //for the computation of activity coefficients
  //
  double group_tmp_UD5000 [] = {0.,0.,0.,0., // group C
				  0.,0.,0.,0., //group C[OH]
				  0.,0.,0.,0., //group Calcohol
				  0.,0.,0.,0., //group Calcohol-tail
				  0.,1.,0.,0.,0., //group C=C
				  0.,0., //group aromatic carbon (AC)
				  0.,0.,0., // group //AC-C
				  0.,  //group OH
				  0., //group H2O
				  0., //group ACOH
				  1.,0., //group ketone
				  1.,   //group aldehyde  
				  0.,0., //group ester
				  0.,0.,0., //group ether 
				  0.,  //group acid
				  0.,   //group ACNO2
				  0.,0.,0., //group NO3
				  0.,0.,0., //group CO-OH
  			          0.,0.,0.,0.,0.,0.,0.,0.,0., //group CO-OC
			          0.,  //group PAN
			          0.,  //group CO-OOH
			          0.,  //group O=COC=O
			          0.,0.,0.}; //group CHxNO2
  
  size = sizeof(group_tmp_UD5000)/sizeof(double);
  assert(size = 60);
	
  for(int i = 0; i < size; ++i)
	UD5000.groups[i] = group_tmp_UD5000[i];

  // Search the species name in the aerosol species list 
  // and add the species if its name matches with
  // the given list.
  add_species_ssh(surrogate, UD5000, species_list_aer, molecular_weight_aer,
              accomodation_coefficient,diffusion_coef,saturation_vapor_pressure,enthalpy_vaporization, henry, t_ref, mass_density, species_part,nlayer,i_hydrophilic,
		  N_inert, N_inorganic);



/* ==== UD5001 ==== */ 

  species UD5001;
  UD5001.name="UD5001";
  UD5001.is_inorganic_precursor=false;
  UD5001.Psat_ref=5.405e-1; 		// Saturation vapor pressure at Tref (torr)
  UD5001.Tref=298;         		// Temperature of reference (K)
  UD5001.deltaH=55.99;     		// Enthalpy of vaporization (kJ/mol)
  UD5001.Henry=0.0;     		// Henry's law constant at Tref (M/atm)
  UD5001.aq_type="none"; 		// "none","diacid","monoacid" or "aldehyde"
  // UD5001.Kacidity1=XXacidityXX;   		// First acidity constant
  UD5001.hydrophilic=true;   	// Does the species condense on the aqueous phase?
  UD5001.hydrophobic=true;  		// Does the species condense on the organic phase?
  UD5001.nonvolatile=false; 		// Is the compound nonvolatile?
  UD5001.kp_from_experiment=false; 	// Use experimental partitioning constant at Tref?
  UD5001.is_organic=true;  		// Is the compound organic?
  UD5001.compute_gamma_org=true;  	// Compute the activity coefficients of the organic phase for this compound?
  UD5001.compute_gamma_aq=true;  	// Compute the activity coefficients of the aqueous phase for this compound?
  UD5001.Koligo_org=0.0;      	//oligomeriation constant in the organic phase
  UD5001.rho=1300.;  
  UD5001.is_monomer=false;
  UD5001.rion=false;
  UD5001.KDiffusion_air=1.0e-5;
  //  UD5001.accomodation_coefficient=alpha;
  UD5001.viscosity=1.68e12;
  UD5001.is_solid=false;
  UD5001.is_generic=false;
  
  //Group: if no functionnal group in the species use the default species
  //for the computation of activity coefficients
  //
  double group_tmp_UD5001 [] = {1.,0.,0.,0., // group C
				  0.,0.,0.,0., //group C[OH]
				  0.,0.,0.,0., //group Calcohol
				  0.,0.,0.,0., //group Calcohol-tail
				  0.,0.,0.,1.,0., //group C=C
				  0.,0., //group aromatic carbon (AC)
				  0.,0.,0., // group //AC-C
				  0.,  //group OH
				  0., //group H2O
				  0., //group ACOH
				  0.,0., //group ketone
				  2.,   //group aldehyde  
				  0.,0., //group ester
				  0.,0.,0., //group ether 
				  0.,  //group acid
				  0.,   //group ACNO2
				  0.,0.,0., //group NO3
				  0.,0.,0., //group CO-OH
  			          0.,0.,0.,0.,0.,0.,0.,0.,0., //group CO-OC
			          0.,  //group PAN
			          0.,  //group CO-OOH
			          0.,  //group O=COC=O
			          0.,0.,0.}; //group CHxNO2
  
  size = sizeof(group_tmp_UD5001)/sizeof(double);
  assert(size = 60);
	
  for(int i = 0; i < size; ++i)
	UD5001.groups[i] = group_tmp_UD5001[i];

  // Search the species name in the aerosol species list 
  // and add the species if its name matches with
  // the given list.
  add_species_ssh(surrogate, UD5001, species_list_aer, molecular_weight_aer,
              accomodation_coefficient,diffusion_coef,saturation_vapor_pressure,enthalpy_vaporization, henry, t_ref, mass_density, species_part,nlayer,i_hydrophilic,
		  N_inert, N_inorganic);



/* ==== UD5002 ==== */ 

  species UD5002;
  UD5002.name="UD5002";
  UD5002.is_inorganic_precursor=false;
  UD5002.Psat_ref=0.015164; 		// Saturation vapor pressure at Tref (torr)
  UD5002.Tref=298;         		// Temperature of reference (K)
  UD5002.deltaH=71.18;     		// Enthalpy of vaporization (kJ/mol)
  UD5002.Henry=0.0;     		// Henry's law constant at Tref (M/atm)
  UD5002.aq_type="none"; 		// "none","diacid","monoacid" or "aldehyde"
  // UD5002.Kacidity1=XXacidityXX;   		// First acidity constant
  UD5002.hydrophilic=true;   	// Does the species condense on the aqueous phase?
  UD5002.hydrophobic=true;  		// Does the species condense on the organic phase?
  UD5002.nonvolatile=false; 		// Is the compound nonvolatile?
  UD5002.kp_from_experiment=false; 	// Use experimental partitioning constant at Tref?
  UD5002.is_organic=true;  		// Is the compound organic?
  UD5002.compute_gamma_org=true;  	// Compute the activity coefficients of the organic phase for this compound?
  UD5002.compute_gamma_aq=true;  	// Compute the activity coefficients of the aqueous phase for this compound?
  UD5002.Koligo_org=0.0;      	//oligomeriation constant in the organic phase
  UD5002.rho=1300.;  
  UD5002.is_monomer=false;
  UD5002.rion=false;
  UD5002.KDiffusion_air=1.0e-5;
  //  UD5002.accomodation_coefficient=alpha;
  UD5002.viscosity=1.68e12;
  UD5002.is_solid=false;
  UD5002.is_generic=false;
  
  //Group: if no functionnal group in the species use the default species
  //for the computation of activity coefficients
  //
  double group_tmp_UD5002 [] = {1.,0.,0.,0., // group C
				  0.,0.,1.,0., //group C[OH]
				  0.,0.,0.,0., //group Calcohol
				  0.,0.,0.,0., //group Calcohol-tail
				  0.,1.,0.,0.,0., //group C=C
				  0.,0., //group aromatic carbon (AC)
				  0.,0.,0., // group //AC-C
				  1.,  //group OH
				  0., //group H2O
				  0., //group ACOH
				  0.,0., //group ketone
				  1.,   //group aldehyde  
				  0.,0., //group ester
				  0.,0.,0., //group ether 
				  0.,  //group acid
				  0.,   //group ACNO2
				  0.,0.,0., //group NO3
				  0.,0.,0., //group CO-OH
  			          0.,0.,0.,0.,0.,0.,0.,0.,0., //group CO-OC
			          0.,  //group PAN
			          0.,  //group CO-OOH
			          0.,  //group O=COC=O
			          0.,0.,0.}; //group CHxNO2
  
  size = sizeof(group_tmp_UD5002)/sizeof(double);
  assert(size = 60);
	
  for(int i = 0; i < size; ++i)
	UD5002.groups[i] = group_tmp_UD5002[i];

  // Search the species name in the aerosol species list 
  // and add the species if its name matches with
  // the given list.
  add_species_ssh(surrogate, UD5002, species_list_aer, molecular_weight_aer,
              accomodation_coefficient,diffusion_coef,saturation_vapor_pressure,enthalpy_vaporization, henry, t_ref, mass_density, species_part,nlayer,i_hydrophilic,
		  N_inert, N_inorganic);



/* ==== UD6000 ==== */ 

  species UD6000;
  UD6000.name="UD6000";
  UD6000.is_inorganic_precursor=false;
  UD6000.Psat_ref=0.0186; 		// Saturation vapor pressure at Tref (torr)
  UD6000.Tref=298;         		// Temperature of reference (K)
  UD6000.deltaH=71.224;     		// Enthalpy of vaporization (kJ/mol)
  UD6000.Henry=0.0;     		// Henry's law constant at Tref (M/atm)
  UD6000.aq_type="none"; 		// "none","diacid","monoacid" or "aldehyde"
  // UD6000.Kacidity1=XXacidityXX;   		// First acidity constant
  UD6000.hydrophilic=true;   	// Does the species condense on the aqueous phase?
  UD6000.hydrophobic=true;  		// Does the species condense on the organic phase?
  UD6000.nonvolatile=false; 		// Is the compound nonvolatile?
  UD6000.kp_from_experiment=false; 	// Use experimental partitioning constant at Tref?
  UD6000.is_organic=true;  		// Is the compound organic?
  UD6000.compute_gamma_org=true;  	// Compute the activity coefficients of the organic phase for this compound?
  UD6000.compute_gamma_aq=true;  	// Compute the activity coefficients of the aqueous phase for this compound?
  UD6000.Koligo_org=0.0;      	//oligomeriation constant in the organic phase
  UD6000.rho=1300.;  
  UD6000.is_monomer=false;
  UD6000.rion=false;
  UD6000.KDiffusion_air=1.0e-5;
  //  UD6000.accomodation_coefficient=alpha;
  UD6000.viscosity=1.68e12;
  UD6000.is_solid=false;
  UD6000.is_generic=false;
  
  //Group: if no functionnal group in the species use the default species
  //for the computation of activity coefficients
  //
  double group_tmp_UD6000 [] = {0.,0.,0.,0., // group C
				  0.,0.,1.,0., //group C[OH]
				  0.,0.,0.,0., //group Calcohol
				  0.,0.,0.,0., //group Calcohol-tail
				  0.,1.,0.,0.,0., //group C=C
				  0.,0., //group aromatic carbon (AC)
				  0.,0.,0., // group //AC-C
				  1.,  //group OH
				  0., //group H2O
				  0., //group ACOH
				  1.,0., //group ketone
				  1.,   //group aldehyde  
				  0.,0., //group ester
				  0.,0.,0., //group ether 
				  0.,  //group acid
				  0.,   //group ACNO2
				  0.,0.,0., //group NO3
				  0.,0.,0., //group CO-OH
  			          0.,0.,0.,0.,0.,0.,0.,0.,0., //group CO-OC
			          0.,  //group PAN
			          0.,  //group CO-OOH
			          0.,  //group O=COC=O
			          0.,0.,0.}; //group CHxNO2
  
  size = sizeof(group_tmp_UD6000)/sizeof(double);
  assert(size = 60);
	
  for(int i = 0; i < size; ++i)
	UD6000.groups[i] = group_tmp_UD6000[i];

  // Search the species name in the aerosol species list 
  // and add the species if its name matches with
  // the given list.
  add_species_ssh(surrogate, UD6000, species_list_aer, molecular_weight_aer,
              accomodation_coefficient,diffusion_coef,saturation_vapor_pressure,enthalpy_vaporization, henry, t_ref, mass_density, species_part,nlayer,i_hydrophilic,
		  N_inert, N_inorganic);



/* ==== UD7000 ==== */ 

  species UD7000;
  UD7000.name="UD7000";
  UD7000.is_inorganic_precursor=false;
  UD7000.Psat_ref=3.168e-4; 		// Saturation vapor pressure at Tref (torr)
  UD7000.Tref=298;         		// Temperature of reference (K)
  UD7000.deltaH=88.65;     		// Enthalpy of vaporization (kJ/mol)
  UD7000.Henry=0.0;     		// Henry's law constant at Tref (M/atm)
  UD7000.aq_type="none"; 		// "none","diacid","monoacid" or "aldehyde"
  // UD7000.Kacidity1=XXacidityXX;   		// First acidity constant
  UD7000.hydrophilic=true;   	// Does the species condense on the aqueous phase?
  UD7000.hydrophobic=true;  		// Does the species condense on the organic phase?
  UD7000.nonvolatile=false; 		// Is the compound nonvolatile?
  UD7000.kp_from_experiment=false; 	// Use experimental partitioning constant at Tref?
  UD7000.is_organic=true;  		// Is the compound organic?
  UD7000.compute_gamma_org=true;  	// Compute the activity coefficients of the organic phase for this compound?
  UD7000.compute_gamma_aq=true;  	// Compute the activity coefficients of the aqueous phase for this compound?
  UD7000.Koligo_org=0.0;      	//oligomeriation constant in the organic phase
  UD7000.rho=1300.;  
  UD7000.is_monomer=false;
  UD7000.rion=false;
  UD7000.KDiffusion_air=1.0e-5;
  //  UD7000.accomodation_coefficient=alpha;
  UD7000.viscosity=1.68e12;
  UD7000.is_solid=false;
  UD7000.is_generic=false;
  
  //Group: if no functionnal group in the species use the default species
  //for the computation of activity coefficients
  //
  double group_tmp_UD7000 [] = {0.,0.,0.,0., // group C
				  0.,0.,0.,0., //group C[OH]
				  0.,0.,0.,0., //group Calcohol
				  0.,0.,0.,0., //group Calcohol-tail
				  0.,0.,0.,1.,0., //group C=C
				  0.,0., //group aromatic carbon (AC)
				  0.,0.,0., // group //AC-C
				  1.,  //group OH
				  0., //group H2O
				  0., //group ACOH
				  1.,1., //group ketone
				  1.,   //group aldehyde  
				  0.,0., //group ester
				  0.,0.,0., //group ether 
				  0.,  //group acid
				  0.,   //group ACNO2
				  0.,0.,0., //group NO3
				  0.,0.,0., //group CO-OH
  			          0.,0.,0.,0.,0.,0.,0.,0.,0., //group CO-OC
			          0.,  //group PAN
			          0.,  //group CO-OOH
			          0.,  //group O=COC=O
			          0.,0.,0.}; //group CHxNO2
  
  size = sizeof(group_tmp_UD7000)/sizeof(double);
  assert(size = 60);
	
  for(int i = 0; i < size; ++i)
	UD7000.groups[i] = group_tmp_UD7000[i];

  add_species_ssh(surrogate, UD7000, species_list_aer, molecular_weight_aer,
                  accomodation_coefficient,diffusion_coef,saturation_vapor_pressure,enthalpy_vaporization, henry, t_ref, mass_density, species_part,nlayer,i_hydrophilic,
                  N_inert, N_inorganic);

/* ==== UU7000 ==== */ 

  species UU7000;
  UU7000.name="UU7000";
  UU7000.is_inorganic_precursor=false;
  UU7000.Psat_ref=6.618e-1; 		// Saturation vapor pressure at Tref (torr)
  UU7000.Tref=298;         		// Temperature of reference (K)
  UU7000.deltaH=54.228;     		// Enthalpy of vaporization (kJ/mol)
  UU7000.Henry=0.0;     		// Henry's law constant at Tref (M/atm)
  UU7000.aq_type="none"; 		// "none","diacid","monoacid" or "aldehyde"
  // UU7000.Kacidity1=XXacidityXX;   		// First acidity constant
  UU7000.hydrophilic=true;   	// Does the species condense on the aqueous phase?
  UU7000.hydrophobic=true;  		// Does the species condense on the organic phase?
  UU7000.nonvolatile=false; 		// Is the compound nonvolatile?
  UU7000.kp_from_experiment=false; 	// Use experimental partitioning constant at Tref?
  UU7000.is_organic=true;  		// Is the compound organic?
  UU7000.compute_gamma_org=true;  	// Compute the activity coefficients of the organic phase for this compound?
  UU7000.compute_gamma_aq=true;  	// Compute the activity coefficients of the aqueous phase for this compound?
  UU7000.Koligo_org=0.0;      	//oligomeriation constant in the organic phase
  UU7000.rho=1300.;  
  UU7000.is_monomer=false;
  UU7000.rion=false;
  UU7000.KDiffusion_air=1.0e-5;
  //  UU7000.accomodation_coefficient=alpha;
  UU7000.viscosity=1.68e12;
  UU7000.is_solid=false;
  UU7000.is_generic=false;
  
  //Group: if no functionnal group in the species use the default species
  //for the computation of activity coefficients
  //
  double group_tmp_UU7000 [] = {0.,0.,0.,0., // group C
				  0.,0.,0.,0., //group C[OH]
				  0.,0.,0.,0., //group Calcohol
				  0.,0.,0.,0., //group Calcohol-tail
				  0.,2.,0.,0.,0., //group C=C
				  0.,0., //group aromatic carbon (AC)
				  0.,0.,0., // group //AC-C
				  0.,  //group OH
				  0., //group H2O
				  0., //group ACOH
				  1.,0., //group ketone
				  1.,   //group aldehyde  
				  0.,0., //group ester
				  0.,0.,0., //group ether 
				  0.,  //group acid
				  0.,   //group ACNO2
				  0.,0.,0., //group NO3
				  0.,0.,0., //group CO-OH
  			          0.,0.,0.,0.,0.,0.,0.,0.,0., //group CO-OC
			          0.,  //group PAN
			          0.,  //group CO-OOH
			          0.,  //group O=COC=O
			          0.,0.,0.}; //group CHxNO2
  
  size = sizeof(group_tmp_UU7000)/sizeof(double);
  assert(size = 60);
	
  for(int i = 0; i < size; ++i)
	UU7000.groups[i] = group_tmp_UU7000[i];

  // Search the species name in the aerosol species list 
  // and add the species if its name matches with
  // the given list.
  add_species_ssh(surrogate, UU7000, species_list_aer, molecular_weight_aer,
              accomodation_coefficient,diffusion_coef,saturation_vapor_pressure,enthalpy_vaporization, henry, t_ref, mass_density, species_part,nlayer,i_hydrophilic,
		  N_inert, N_inorganic);
  
  /* TOLexp SPECIES END */

  // GENOA DELETE DETECTOR END //
  // NOTE: When running ./clean under genoa clean mode (-gcl=cleanall), the code block after "GENOA DELETE DETECTOR START" and before "GENOA DELETE DETECTOR END" will be removed.

  species GLY;
  GLY.name="GLY";
  GLY.is_inorganic_precursor=false;
  GLY.Psat_ref=17.2*100; // Saturation vapor pressure at Tref (torr)
  GLY.Tref=298;         // Temperature of reference (K)
  GLY.deltaH=50.0;     // Enthalpy of vaporization (kJ/mol)
  GLY.Henry=0.0;     // Henry's law constant at Tref (M/atm)
  //GLY.aq_type="hydrate"; // "none","diacid","monoacid" or "aldehyde"
  GLY.hydrophilic=true;   // Does the species condense on the aqueous phase?
  GLY.hydrophobic=false;  // Does the species condense on the organic phase?
  GLY.nonvolatile=false; // Is the compound nonvolatile?
  GLY.kp_from_experiment=false;  // Use experimental partitioning constant at Tref?
  GLY.is_organic=true;  // Is the compound organic?
  GLY.compute_gamma_org=true;  // Compute the activity coefficients of the organic phase for this compound?
  GLY.compute_gamma_aq=true;  // Compute the activity coefficients of the aqueous phase for this compound?
  GLY.Koligo_org=0.0;      //oligomeriation constant in the organic phase
  GLY.rho=1300.0;  
  GLY.is_monomer=false;
  GLY.rion=false;
  //GLY.KDiffusion_air=1.0e-5;
  //  GLY.accomodation_coefficient=alpha;
  GLY.viscosity=1.68e12;
  GLY.is_solid=false;
  GLY.is_generic=false;
  /*GLY.nion=2;
  GLY.rion_catalyzed=false;    
  GLY.ion.resize(GLY.nion);
  GLY.ion(0)="HSO4";
  GLY.ion(1)="NO3";
  GLY.kion.resize(GLY.nion);
  GLY.kion(0)=1.53e-7;
  GLY.kion(1)=1.53e-7;
  GLY.rion_product.resize(GLY.nion);
  GLY.rion_product(0)="GLY2";
  GLY.rion_product(1)="GLY2";
  GLY.rion_catalyzed(1)=false;
  GLY.rion_catalyzed(2)=false;*/
  
  //Group: if no functionnal group in the species use the default species
  //for the computation of activity coefficients
  //
  double group_tmp_gly [] = {0.0,0.0,0.0,0.0, // group C
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
			       2.0,   //group aldehyde  
			       0.0,0.0, //group ester
			       0.0,0.0,0.0, //group ether 
			       0.0,  //group acid
			       0.0,   //group ACNO2
			       0.0,0.0,0.0, //group NO3
			       0.0,0.0,0.0, //group CO-OH
  			       0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0, //group CO-OC
			       0.0,    //group PAN
			       0.0,  //group CO-OOH
			       0.0,  //group O=COC=O
			       0.0,0.0,0.0}; //group CHxNO2
  
  size = sizeof(group_tmp_gly)/sizeof(double);
  assert(size == 60);
	
  for(int i = 0; i < size; ++i)
    GLY.groups[i] = group_tmp_gly[i];

  // Search the species name in the aerosol species list 
  // and add the species if its name matches with
  // the given list.
  add_species_ssh(surrogate, GLY, species_list_aer, molecular_weight_aer,
              accomodation_coefficient,diffusion_coef,saturation_vapor_pressure,enthalpy_vaporization, henry, t_ref, mass_density, species_part,nlayer,i_hydrophilic,
		  N_inert, N_inorganic);


  
  species GLYOH;
  GLYOH.name="GLYOH";
  GLYOH.is_inorganic_precursor=false;
  GLYOH.Psat_ref=1.54e-2*100; // Saturation vapor pressure at Tref (torr)
  GLYOH.Tref=298;         // Temperature of reference (K)
  GLYOH.deltaH=50.0;     // Enthalpy of vaporization (kJ/mol)
  GLYOH.Henry=0.0;     // Henry's law constant at Tref (M/atm)
  //GLYOH.aq_type="hydrate"; // "none","diacid","monoacid" or "aldehyde"
  GLYOH.hydrophilic=true;   // Does the species condense on the aqueous phase?
  GLYOH.hydrophobic=false;  // Does the species condense on the organic phase?
  GLYOH.nonvolatile=false; // Is the compound nonvolatile?
  GLYOH.kp_from_experiment=false;  // Use experimental partitioning constant at Tref?
  GLYOH.is_organic=true;  // Is the compound organic?
  GLYOH.compute_gamma_org=true;  // Compute the activity coefficients of the organic phase for this compound?
  GLYOH.compute_gamma_aq=true;  // Compute the activity coefficients of the aqueous phase for this compound?
  GLYOH.Koligo_org=0.0;      //oligomeriation constant in the organic phase
  GLYOH.rho=1300.0;  
  GLYOH.is_monomer=false;
  GLYOH.rion=false;
  //GLYOH.KDiffusion_air=1.0e-5;
  //  GLYOH.accomodation_coefficient=alpha;
  GLYOH.viscosity=1.68e12;
  GLYOH.is_solid=false;
  GLYOH.is_generic=false;
  /*
  GLYOH.nion=2;
  GLYOH.rion_catalyzed=false;    
  GLYOH.ion.resize(GLYOH.nion);
  GLYOH.ion(0)="HSO4";
  GLYOH.ion(1)="NO3";
  GLYOH.kion.resize(GLYOH.nion);
  GLYOH.kion(0)=1.53e-7;
  GLYOH.kion(1)=1.53e-7;
  GLYOH.rion_product.resize(GLYOH.nion);
  GLYOH.rion_product(0)="GLYOH2";
  GLYOH.rion_product(1)="GLYOH2";
  GLYOH.rion_catalyzed(1)=false;
  GLYOH.rion_catalyzed(2)=false;*/
  
  //Group: if no functionnal group in the species use the default species
  //for the computation of activity coefficients
  //
  double group_tmp_glyoh [] = {0.0,0.0,0.0,0.0, // group C
			       0.0,0.0,1.0,0.0, //group C[OH]
			       0.0,0.0,0.0,0.0, //group Calcohol
			       0.0,0.0,0.0,0.0, //group Calcohol-tail
			       0.0,0.0,0.0,0.0,0.0, //group C=C
			       0.0,0.0, //group aromatic carbon (AC)
			       0.0,0.0,0.0, // group //AC-C
			       2.0,  //group OH
			       0.0, //group H2O
			       0.0, //group ACOH
			       0.0,0.0, //group ketone
			       1.0,   //group aldehyde  
			       0.0,0.0, //group ester
			       0.0,0.0,0.0, //group ether 
			       0.0,  //group acid
			       0.0,   //group ACNO2
			       0.0,0.0,0.0, //group NO3
			       0.0,0.0,0.0, //group CO-OH
  			       0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0, //group CO-OC
			       0.0,    //group PAN
			       0.0,  //group CO-OOH
			       0.0,  //group O=COC=O
			       0.0,0.0,0.0}; //group CHxNO2
  
  size = sizeof(group_tmp_glyoh)/sizeof(double);
  assert(size == 60);
	
  for(int i = 0; i < size; ++i)
    GLYOH.groups[i] = group_tmp_glyoh[i];

  // Search the species name in the aerosol species list 
  // and add the species if its name matches with
  // the given list.
  add_species_ssh(surrogate, GLYOH, species_list_aer, molecular_weight_aer,
              accomodation_coefficient,diffusion_coef,saturation_vapor_pressure,enthalpy_vaporization, henry, t_ref, mass_density, species_part,nlayer,i_hydrophilic,
		  N_inert, N_inorganic);


    
  species GLYOHOH;
  GLYOHOH.name="GLYOHOH";
  GLYOHOH.is_inorganic_precursor=false;
  GLYOHOH.Psat_ref=1.38e-5*100; // Saturation vapor pressure at Tref (torr)
  GLYOHOH.Tref=298;         // Temperature of reference (K)
  GLYOHOH.deltaH=50.0;     // Enthalpy of vaporization (kJ/mol)
  GLYOHOH.Henry=0.0;     // Henry's law constant at Tref (M/atm)
  GLYOHOH.aq_type="none"; // "none","diacid","monoacid" or "aldehyde"
  GLYOHOH.hydrophilic=true;   // Does the species condense on the aqueous phase?
  GLYOHOH.hydrophobic=false;  // Does the species condense on the organic phase?
  GLYOHOH.nonvolatile=false; // Is the compound nonvolatile?
  GLYOHOH.kp_from_experiment=false;  // Use experimental partitioning constant at Tref?
  GLYOHOH.is_organic=true;  // Is the compound organic?
  GLYOHOH.compute_gamma_org=true;  // Compute the activity coefficients of the organic phase for this compound?
  GLYOHOH.compute_gamma_aq=true;  // Compute the activity coefficients of the aqueous phase for this compound?
  GLYOHOH.Koligo_org=0.0;      //oligomeriation constant in the organic phase
  GLYOHOH.rho=1300.0;  
  GLYOHOH.is_monomer=false;
  GLYOHOH.rion=false;
  //GLYOHOH.KDiffusion_air=1.0e-5;
  //  GLYOHOH.accomodation_coefficient=alpha;
  GLYOHOH.viscosity=1.68e12;
  GLYOHOH.is_solid=false;
  GLYOHOH.is_generic=false;
  /*
  GLYOHOH.nion=2;
  GLYOHOH.rion_catalyzed=false;    
  GLYOHOH.ion.resize(GLYOHOH.nion);
  GLYOHOH.ion(0)="HSO4";
  GLYOHOH.ion(1)="NO3";
  GLYOHOH.kion.resize(GLYOHOH.nion);
  GLYOHOH.kion(0)=1.53e-7;
  GLYOHOH.kion(1)=1.53e-7;
  GLYOHOH.rion_product.resize(GLYOHOH.nion);
  GLYOHOH.rion_product(0)="GLYOHOH2";
  GLYOHOH.rion_product(1)="GLYOHOH2";
  GLYOHOH.rion_catalyzed(1)=false;
  GLYOHOH.rion_catalyzed(2)=false;*/
  
  //Group: if no functionnal group in the species use the default species
  //for the computation of activity coefficients
  //
  double group_tmp_glyohoh [] = {0.0,0.0,0.0,0.0, // group C
			       0.0,0.0,2.0,0.0, //group C[OH]
			       0.0,0.0,0.0,0.0, //group Calcohol
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
			       0.0,    //group PAN
			       0.0,  //group CO-OOH
			       0.0,  //group O=COC=O
			       0.0,0.0,0.0}; //group CHxNO2
  
  size = sizeof(group_tmp_glyohoh)/sizeof(double);
  assert(size == 60);
	
  for(int i = 0; i < size; ++i)
    GLYOHOH.groups[i] = group_tmp_glyohoh[i];


  // Search the species name in the aerosol species list 
  // and add the species if its name matches with
  // the given list.
  add_species_ssh(surrogate, GLYOHOH, species_list_aer, molecular_weight_aer,
                  accomodation_coefficient,diffusion_coef,saturation_vapor_pressure,enthalpy_vaporization, henry, t_ref, mass_density, species_part,nlayer,i_hydrophilic,
		  N_inert, N_inorganic);

  /* ==== BiA1D ==== */ 

  species BiA1D;
  BiA1D.name="BiA1D";
  BiA1D.is_inorganic_precursor=false;
  BiA1D.Psat_ref=2.98e-6; // Saturation vapor pressure at Tref (torr)
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
              accomodation_coefficient,diffusion_coef,saturation_vapor_pressure,enthalpy_vaporization, henry, t_ref, mass_density, species_part,nlayer,i_hydrophilic,
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
  BiA0D.aq_type="none"; // "none","diacid","monoacid" or "aldehyde"
  BiA0D.Koligo_org=0.0;
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
		  accomodation_coefficient,diffusion_coef,saturation_vapor_pressure,enthalpy_vaporization, henry, t_ref, mass_density, species_part,nlayer,i_hydrophilic,
		  N_inert, N_inorganic);
  
  species BiMT;
  BiMT.name="BiMT";
  BiMT.is_inorganic_precursor=false;
  BiMT.Psat_ref=1.45e-6; // Saturation vapor pressure at Tref (torr)
  BiMT.Tref=298;         // Temperature of reference (K)
  BiMT.deltaH=38.4;     // Enthalpy of vaporization (kJ/mol)
  BiMT.Henry=0.;     // Henry's law constant at Tref (M/atm)
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
		  accomodation_coefficient,diffusion_coef,saturation_vapor_pressure,enthalpy_vaporization, henry, t_ref, mass_density, species_part,nlayer,i_hydrophilic,
		  N_inert, N_inorganic);

  species BiPER;
  BiPER.name="BiPER";
  BiPER.is_inorganic_precursor=false;
  BiPER.Psat_ref=2.61e-6; // Saturation vapor pressure at Tref (torr)
  BiPER.Tref=298;         // Temperature of reference (K)
  BiPER.deltaH=38.4;     // Enthalpy of vaporization (kJ/mol)
  BiPER.Henry=0.;     // Henry's law constant at Tref (M/atm)
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
		  accomodation_coefficient,diffusion_coef,saturation_vapor_pressure,enthalpy_vaporization, henry, t_ref, mass_density, species_part,nlayer,i_hydrophilic,
		  N_inert, N_inorganic);

  species BiDER;
  BiDER.name="BiDER";
  BiDER.is_inorganic_precursor=false;
  BiDER.Psat_ref=4.10e-7; // Saturation vapor pressure at Tref (torr)
  BiDER.Tref=298;         // Temperature of reference (K)
  BiDER.deltaH=38.4;     // Enthalpy of vaporization (kJ/mol)
  BiDER.Henry=0.;     // Henry's law constant at Tref (M/atm)
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
              accomodation_coefficient,diffusion_coef,saturation_vapor_pressure,enthalpy_vaporization, henry, t_ref, mass_density, species_part,nlayer,i_hydrophilic,
	      N_inert, N_inorganic);

  species BiMGA;
  BiMGA.name="BiMGA";
  BiMGA.is_inorganic_precursor=false;
  BiMGA.Psat_ref=1.40e-5; // Saturation vapor pressure at Tref (torr)
  BiMGA.Tref=298;         // Temperature of reference (K)
  BiMGA.deltaH=43.2;     // Enthalpy of vaporization (kJ/mol)
  BiMGA.Henry=0.;     // Henry's law constant at Tref (M/atm)
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
  BiMGA.Koligo_org=0.0;
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
              accomodation_coefficient,diffusion_coef,saturation_vapor_pressure,enthalpy_vaporization, henry, t_ref, mass_density, species_part,nlayer,i_hydrophilic,
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
              accomodation_coefficient,diffusion_coef,saturation_vapor_pressure,enthalpy_vaporization, henry, t_ref, mass_density, species_part,nlayer,i_hydrophilic,
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
              accomodation_coefficient,diffusion_coef,saturation_vapor_pressure,enthalpy_vaporization, henry, t_ref, mass_density, species_part,nlayer,i_hydrophilic,
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
              accomodation_coefficient,diffusion_coef,saturation_vapor_pressure,enthalpy_vaporization, henry, t_ref, mass_density, species_part,nlayer,i_hydrophilic,
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
              accomodation_coefficient,diffusion_coef,saturation_vapor_pressure,enthalpy_vaporization, henry, t_ref, mass_density, species_part,nlayer,i_hydrophilic,
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
  AnClP.Koligo_org=0.0;  
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
              accomodation_coefficient,diffusion_coef,saturation_vapor_pressure,enthalpy_vaporization, henry, t_ref, mass_density, species_part,nlayer,i_hydrophilic,
	      N_inert, N_inorganic);

  species AnBOAhP;
  AnBOAhP.name="AnBOAhP";
  AnBOAhP.Psat_ref=8.63e-5;
  AnBOAhP.is_inorganic_precursor=false;
  AnBOAhP.MM=167.84;           // Molar mass (g/mol)
  AnBOAhP.nonvolatile=false;  // Is the compound nonvolatile?
  AnBOAhP.Henry=0.;
  AnBOAhP.hydrophilic=true; // Does the species condense on the aqueous phase?
  AnBOAhP.hydrophobic=true;  // Does the species condense on the organic phase?
  AnBOAhP.kp_from_experiment=false;  // Use experimental partitioning constant at Tref?
  AnBOAhP.kp_experiment=0.00023;       // Value of the experimental partitioning constant at Tref?
  AnBOAhP.deltaH=79.0;     // Enthalpy of vaporization (kJ/mol)
  AnBOAhP.Tref=298;
  AnBOAhP.is_organic=true;  // Is the compound organic?
  AnBOAhP.compute_gamma_org=true;  // Compute the activity coefficients of the organic phase for this compound?
  AnBOAhP.compute_gamma_aq=true;  // Compute the activity coefficients of the aqueous phase for this compound
  AnBOAhP.Koligo_org=0.0;         //oligomeriation constant in the organic phase
  AnBOAhP.rho=1300.0;
  AnBOAhP.is_monomer=false;
  AnBOAhP.rion=false;
  AnBOAhP.aq_type="none";
  AnBOAhP.viscosity=1.68e12;  
  AnBOAhP.is_solid=false;
  AnBOAhP.is_generic=false;

  //Group: if no functionnal group in the species use the default species
  //for the computation of activity coefficients 
  
  double group_tmp_anboahp [] = {0.33,3.3,0.0,0.0, // group C
                               0.0,0.0,0.0,0.0, //group C[OH]
                               0.0,0.0,0.0,0.0, //group Calcohol
                               0.0,0.0,0.0,0.0, //group Calcohol-tail
                               0.0,0.0,0.0,0.0,0.0, //group C=C
                               2.01,2.01, //group aromatic carbon (AC)
                               0.0,0.0,0.0, // group //AC-C
                               0.0,  //group OH
                               0.0, //group H2O
                               0.67, //group ACOH
                               0.0,0.0, //group ketone
                               0.67,   //group aldehyde  
                               0.0,0.0, //group ester
                               0.67,0.0,0.0, //group ether  
                               0.33,  //group acid
                               0.0,   //group ACNO2
                               0.0,0.0,0.0, //group NO3
                               0.0,0.0,0.0, //group CO-OH
			       0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0, //group CO-OC
			       0.0,  //group PAN
                               0.0,  //group CO-OOH
                               0.0,  //group O=COC=O
                               0.0,0.0,0.0}; //group CHxNO2


  size = sizeof(group_tmp_anboahp)/sizeof(double);
  assert(size == 60);
  for(int i = 0; i < size; ++i)
    AnBOAhP.groups[i] = group_tmp_anboahp[i];

  add_species_ssh(surrogate, AnBOAhP, species_list_aer, molecular_weight_aer,
                  accomodation_coefficient,diffusion_coef,saturation_vapor_pressure,enthalpy_vaporization, henry, t_ref, mass_density, species_part,nlayer,i_hydrophilic,
                  N_inert, N_inorganic);

  
  species AnBOAlP;
  AnBOAlP.name="AnBOAlP";
  AnBOAlP.Psat_ref=4.69e-9; 
  AnBOAlP.is_inorganic_precursor=false;
  AnBOAlP.MM=210.0;           // Molar mass (g/mol)
  AnBOAlP.nonvolatile=false;  // Is the compound nonvolatile?
  AnBOAlP.Henry=0.;
  AnBOAlP.hydrophilic=true; // Does the species condense on the aqueous phase?
  AnBOAlP.hydrophobic=true;  // Does the species condense on the organic phase?
  AnBOAlP.kp_from_experiment=false;  // Use experimental partitioning constant at Tref?
  AnBOAlP.kp_experiment=18.3;       // Value of the experimental partitioning constant at Tref?
  AnBOAlP.deltaH=106.0;     // Enthalpy of vaporization (kJ/mol)
  AnBOAlP.Tref=298;
  AnBOAlP.is_organic=true;  // Is the compound organic?
  AnBOAlP.compute_gamma_org=true;  // Compute the activity coefficients of the organic phase for this compound?
  AnBOAlP.compute_gamma_aq=true;  // Compute the activity coefficients of the aqueous phase for this compound
  AnBOAlP.Koligo_org=0.0;         //oligomeriation constant in the organic phase
  AnBOAlP.rho=1300.0;
  AnBOAlP.is_monomer=false;
  AnBOAlP.rion=false;
  AnBOAlP.aq_type="none";
  AnBOAlP.viscosity=1.68e12;  
  AnBOAlP.is_solid=false;
  AnBOAlP.is_generic=false;

  //Group: if no functionnal group in the species use the default species
  //for the computation of activity coefficients 

  double group_tmp_anboalp [] = {0.0,0.0,0.0,0.0, // group C
                               0.0,0.0,0.0,0.0, //group C[OH]
                               0.0,0.0,0.0,0.0, //group Calcohol
                               0.0,0.0,0.0,0.0, //group Calcohol-tail
                               0.0,0.0,0.0,0.0,0.0, //group C=C
                               3.0,3.0, //group aromatic carbon (AC)
                               0.0,0.0,0.0, // group //AC-C
                               0.0,  //group OH
                               0.0, //group H2O
                               0.0, //group ACOH
                               0.0,0.0, //group ketone
                               0.0,   //group aldehyde  
                               0.0,0.0, //group ester
                               0.0,0.0,0.0, //group ether  
                               3.,  //group acid
                               0.0,   //group ACNO2
                               0.0,0.0,0.0, //group NO3
                               0.0,0.0,0.0, //group CO-OH
			       0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0, //group CO-OC
			       0.0,  //group PAN
                               0.0,  //group CO-OOH
                               0.0,  //group O=COC=O
                               0.0,0.0,0.0}; //group CHxNO2


  size = sizeof(group_tmp_anboalp)/sizeof(double);
  assert(size == 60);
  for(int i = 0; i < size; ++i)
    AnBOAlP.groups[i] = group_tmp_anboalp[i];

  add_species_ssh(surrogate, AnBOAlP, species_list_aer, molecular_weight_aer,
                  accomodation_coefficient,diffusion_coef,saturation_vapor_pressure,enthalpy_vaporization, henry, t_ref, mass_density, species_part,nlayer,i_hydrophilic,
                  N_inert, N_inorganic);


  species AnBOAmP;
  AnBOAmP.name="AnBOAmP";
  AnBOAmP.Psat_ref=6.04e-7;
  AnBOAmP.is_inorganic_precursor=false;
  AnBOAmP.MM=180.98;           // Molar mass (g/mol)
  AnBOAmP.nonvolatile=false;  // Is the compound nonvolatile?
  AnBOAmP.Henry=0.;
  AnBOAmP.hydrophilic=true; // Does the species condense on the aqueous phase?
  AnBOAmP.hydrophobic=true;  // Does the species condense on the organic phase?
  AnBOAmP.kp_from_experiment=false;  // Use experimental partitioning constant at Tref?
  AnBOAmP.kp_experiment=0.040;       // Value of the experimental partitioning constant at Tref?
  AnBOAmP.deltaH=91.0;     // Enthalpy of vaporization (kJ/mol)
  AnBOAmP.Tref=298;
  AnBOAmP.is_organic=true;  // Is the compound organic?
  AnBOAmP.compute_gamma_org=true;  // Compute the activity coefficients of the organic phase for this compound?
  AnBOAmP.compute_gamma_aq=true;  // Compute the activity coefficients of the aqueous phase for this compound
  AnBOAmP.Koligo_org=0.0;         //oligomeriation constant in the organic phase
  AnBOAmP.rho=1300.0;
  AnBOAmP.is_monomer=false;
  AnBOAmP.rion=false;
  AnBOAmP.aq_type="none";
  AnBOAmP.viscosity=1.68e12;  
  AnBOAmP.is_solid=false;
  AnBOAmP.is_generic=false;

  //Group: if no functionnal group in the species use the default species
  //for the computation of activity coefficients 
  
  double group_tmp_anboamp [] = {0.17,2.38,0.0,0.0, // group C
                               0.0,0.0,0.99,0.0, //group C[OH]
                               0.0,0.0,0.33,0.0, //group Calcohol
                               0.0,0.0,0.0,0.0, //group Calcohol-tail
                               0.0,0.0,0.0,0.0,0.0, //group C=C
                               1.5,1.5, //group aromatic carbon (AC)
                               0.0,0.0,0.0, // group //AC-C
                               0.99,  //group OH
                               0.0, //group H2O
                               0.5, //group ACOH
                               0.0,0.0, //group ketone
                               0.0,   //group aldehyde  
                               0.0,0.0, //group ester
                               0.5,0.33,0.33, //group ether 
                               0.67,  //group acid
                               0.0,   //group ACNO2
                               0.0,0.0,0.0, //group NO3
                               0.0,0.0,0.0, //group CO-OH
			       0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0, //group CO-OC
			       0.0,  //group PAN
                               0.0,  //group CO-OOH
                               0.0,  //group O=COC=O
                               0.0,0.0,0.0}; //group CHxNO2


  size = sizeof(group_tmp_anboamp)/sizeof(double);
  assert(size == 60);
  for(int i = 0; i < size; ++i)
    AnBOAmP.groups[i] = group_tmp_anboamp[i];

  add_species_ssh(surrogate, AnBOAmP, species_list_aer, molecular_weight_aer,
                  accomodation_coefficient,diffusion_coef,saturation_vapor_pressure,enthalpy_vaporization, henry, t_ref, mass_density, species_part,nlayer,i_hydrophilic,
                  N_inert, N_inorganic);

  species AnBSOAhP;
  AnBSOAhP.name="AnBSOAhP";
  AnBSOAhP.Psat_ref=8.63e-7;
  AnBSOAhP.is_inorganic_precursor=false;
  AnBSOAhP.MM=234.976;           // Molar mass (g/mol)
  AnBSOAhP.nonvolatile=false;  // Is the compound nonvolatile?
  AnBSOAhP.Henry=0.;
  AnBSOAhP.hydrophilic=true; // Does the species condense on the aqueous phase?
  AnBSOAhP.hydrophobic=true;  // Does the species condense on the organic phase?
  AnBSOAhP.kp_from_experiment=false;  // Use experimental partitioning constant at Tref?
  AnBSOAhP.kp_experiment=0.00023;       // Value of the experimental partitioning constant at Tref?
  AnBSOAhP.deltaH=79.0;     // Enthalpy of vaporization (kJ/mol)
  AnBSOAhP.Tref=298;
  AnBSOAhP.is_organic=true;  // Is the compound organic?
  AnBSOAhP.compute_gamma_org=true;  // Compute the activity coefficients of the organic phase for this compound?
  AnBSOAhP.compute_gamma_aq=true;  // Compute the activity coefficients of the aqueous phase for this compound?
  AnBSOAhP.Koligo_org=0.0;         //oligomeriation constant in the organic phase
  AnBSOAhP.rho=1300.0;
  AnBSOAhP.is_monomer=false;
  AnBSOAhP.rion=false;
  AnBSOAhP.aq_type="none";
  AnBSOAhP.viscosity=1.68e12;  
  AnBSOAhP.is_solid=false;
  AnBSOAhP.is_generic=false;
  
  size = sizeof(group_tmp_anboahp)/sizeof(double);
  assert(size == 60);
  for(int i = 0; i < size; ++i)
    AnBSOAhP.groups[i] = group_tmp_anboahp[i];

  add_species_ssh(surrogate, AnBSOAhP, species_list_aer, molecular_weight_aer,
                  accomodation_coefficient,diffusion_coef,saturation_vapor_pressure,enthalpy_vaporization, henry, t_ref, mass_density, species_part,nlayer,i_hydrophilic,
                  N_inert, N_inorganic);

  species AnBSOAmP;
  AnBSOAmP.name="AnBSOAmP";
  AnBSOAmP.Psat_ref=6.04e-9;
  AnBSOAmP.is_inorganic_precursor=false;
  AnBSOAmP.MM=253.372;           // Molar mass (g/mol)
  AnBSOAmP.nonvolatile=false;  // Is the compound nonvolatile?
  AnBSOAmP.Henry=0.;
  AnBSOAmP.hydrophilic=true; // Does the species condense on the aqueous phase?
  AnBSOAmP.hydrophobic=true;  // Does the species condense on the organic phase?
  AnBSOAmP.kp_from_experiment=false;  // Use experimental partitioning constant at Tref?
  AnBSOAmP.kp_experiment=0.040;       // Value of the experimental partitioning constant at Tref?
  AnBSOAmP.deltaH=91.0;     // Enthalpy of vaporization (kJ/mol)
  AnBSOAmP.Tref=298;
  AnBSOAmP.is_organic=true;  // Is the compound organic?
  AnBSOAmP.compute_gamma_org=true;  // Compute the activity coefficients of the organic phase for this compound?
  AnBSOAmP.compute_gamma_aq=true;  // Compute the activity coefficients of the aqueous phase for this compound?
  AnBSOAmP.Koligo_org=0.0;         //oligomeriation constant in the organic phase
  AnBSOAmP.rho=1300.0;
  AnBSOAmP.is_monomer=false;
  AnBSOAmP.rion=false;
  AnBSOAmP.aq_type="none";
  AnBSOAmP.viscosity=1.68e12;  
  AnBSOAmP.is_solid=false;
  AnBSOAmP.is_generic=false;

  //Group: if no functionnal group in the species use the default species
  //for the computation of activity coefficients 

  size = sizeof(group_tmp_anboamp)/sizeof(double);
  assert(size == 60);
  for(int i = 0; i < size; ++i)
    AnBSOAmP.groups[i] = group_tmp_anboamp[i];

  add_species_ssh(surrogate, AnBSOAmP, species_list_aer, molecular_weight_aer,
                  accomodation_coefficient,diffusion_coef,saturation_vapor_pressure,enthalpy_vaporization, henry, t_ref, mass_density, species_part,nlayer,i_hydrophilic,
                  N_inert, N_inorganic);

  species BiNGA;
  BiNGA.name="BiNGA";
  BiNGA.is_inorganic_precursor=false;
  BiNGA.Psat_ref=1.40e-5; // Saturation vapor pressure at Tref (torr)
  BiNGA.Tref=298;         // Temperature of reference (K)
  BiNGA.deltaH=43.2;     // Enthalpy of vaporization (kJ/mol)
  BiNGA.Henry=0.;     // Henry's law constant at Tref (M/atm)
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
              accomodation_coefficient,diffusion_coef,saturation_vapor_pressure,enthalpy_vaporization, henry, t_ref, mass_density, species_part,nlayer,i_hydrophilic,
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
		  accomodation_coefficient,diffusion_coef,saturation_vapor_pressure,enthalpy_vaporization, henry, t_ref, mass_density, species_part,nlayer,i_hydrophilic,
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
		  accomodation_coefficient,diffusion_coef,saturation_vapor_pressure,enthalpy_vaporization, henry, t_ref, mass_density, species_part,nlayer,i_hydrophilic,
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
  POAlP.Koligo_org=0.0;  
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
		  accomodation_coefficient,diffusion_coef,saturation_vapor_pressure,enthalpy_vaporization, henry, t_ref, mass_density, species_part,nlayer,i_hydrophilic,
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
  POAmP.Koligo_org=0.0;  
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
              accomodation_coefficient,diffusion_coef,saturation_vapor_pressure,enthalpy_vaporization, henry, t_ref, mass_density, species_part,nlayer,i_hydrophilic,
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
  POAhP.Koligo_org=0.0;  
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
              accomodation_coefficient,diffusion_coef,saturation_vapor_pressure,enthalpy_vaporization, henry, t_ref, mass_density, species_part,nlayer,i_hydrophilic,
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
  SOAlP.Koligo_org=0.0;  
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
		  accomodation_coefficient,diffusion_coef,saturation_vapor_pressure,enthalpy_vaporization, henry, t_ref, mass_density, species_part,nlayer,i_hydrophilic,
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
  SOAmP.Koligo_org=0.0;  
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
		  accomodation_coefficient,diffusion_coef,saturation_vapor_pressure,enthalpy_vaporization, henry, t_ref, mass_density, species_part,nlayer,i_hydrophilic,
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
  SOAhP.Koligo_org=0.0;  
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
		  accomodation_coefficient,diffusion_coef,saturation_vapor_pressure,enthalpy_vaporization, henry, t_ref, mass_density, species_part,nlayer,i_hydrophilic,
	          N_inert, N_inorganic);

  //SOA species to account for non-reversible transformation into non-volatile species
  species irrSOA;
  irrSOA.name="irrSOA";
  irrSOA.is_inorganic_precursor=false;
  irrSOA.nonvolatile=true;  // Is the compound nonvolatile?
  irrSOA.hydrophilic=false; // Does the species condense on the aqueous phase?
  irrSOA.hydrophobic=true;  // Does the species condense on the organic phase?
  irrSOA.kp_from_experiment=false;  // Use experimental partitioning constant at Tref?
  irrSOA.kp_experiment=110.0;       // Value of the experimental partitioning constant at Tref?
  irrSOA.deltaH=106.0;     // Enthalpy of vaporization (kJ/mol)
  irrSOA.Tref=298;
  irrSOA.is_organic=true;  // Is the compound organic?
  irrSOA.compute_gamma_org=true;  // Compute the activity coefficients of the organic phase for this compound?
  irrSOA.compute_gamma_aq=true;  // Compute the activity coefficients of the aqueous phase for this compound
  irrSOA.rho=1300.0;
  irrSOA.Koligo_org=0.0;  
  irrSOA.is_monomer=false;
  irrSOA.rion=false;
  //irrSOA.KDiffusion_air=1.0e-5;
  //  irrSOA.accomodation_coefficient=alpha;
  irrSOA.viscosity=1.68e12;  
  irrSOA.is_solid=false;
  irrSOA.is_generic=false;
 
  //Group: if no functionnal group in the species use the default species
  //for the computation of activity coefficients 
 
  for(int i = 0; i < 60; ++i)
    irrSOA.groups[i] = 0.0;

  // Search the species name in the aerosol species list 
  // and add the species if its name matches with
  // the given list.
  add_species_ssh(surrogate, irrSOA, species_list_aer, molecular_weight_aer,
		  accomodation_coefficient,diffusion_coef,saturation_vapor_pressure,enthalpy_vaporization, henry, t_ref, mass_density, species_part,nlayer,i_hydrophilic,
		  N_inert, N_inorganic);

  
  species BOAlP;
  BOAlP.name="BOAlP";
  BOAlP.is_inorganic_precursor=false;
  BOAlP.nonvolatile=false;  // Is the compound nonvolatile?
  BOAlP.hydrophilic=false; // Does the species condense on the aqueous phase?
  BOAlP.hydrophobic=true;  // Does the species condense on the organic phase?
  BOAlP.kp_from_experiment=true;  // Use experimental partitioning constant at Tref?
  //  BOAlP.kp_experiment=0.0435;       
  //   BOAlP.kp_experiment=0.0015;       
  BOAlP.kp_experiment=18.3;       // Value of the experimental partitioning constant at Tref?
  BOAlP.deltaH=106.0;     // Enthalpy of vaporization (kJ/mol)
  BOAlP.Tref=298;
  BOAlP.is_organic=true;  // Is the compound organic?
  BOAlP.compute_gamma_org=true;  // Compute the activity coefficients of the organic phase for this compound?
  BOAlP.compute_gamma_aq=true;  // Compute the activity coefficients of the aqueous phase for this compound?
  BOAlP.Koligo_org=0.0;         //oligomeriation constant in the organic phase
  BOAlP.rho=1300.0;
  BOAlP.is_monomer=false;
  BOAlP.rion=false;
  // BOAlP.KDiffusion_air=1.0e-5;
  //  BOAlP.accomodation_coefficient=alpha;
  BOAlP.viscosity=1.68e12;  
  BOAlP.is_solid=false;
  BOAlP.is_generic=false;

  //Group: if no functionnal group in the species use the default species
  //for the computation of activity coefficients 
  
  double group_tmp_boalp [] = {0.0,0.0,0.0,0.0, // group C
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
  
  size = sizeof(group_tmp_boalp)/sizeof(double);
  assert(size == 60);
  for(int i = 0; i < size; ++i)
    BOAlP.groups[i] = group_tmp_boalp[i];

  // Search the species name in the aerosol species list 
  // and add the species if its name matches with
  // the given list.
  add_species_ssh(surrogate, BOAlP, species_list_aer, molecular_weight_aer,
		  accomodation_coefficient,diffusion_coef,saturation_vapor_pressure,enthalpy_vaporization, henry, t_ref, mass_density, species_part,nlayer,i_hydrophilic,
  	          N_inert, N_inorganic);
  
  species BOAmP;
  BOAmP.name="BOAmP";
  BOAmP.is_inorganic_precursor=false;
  BOAmP.nonvolatile=false;  // Is the compound nonvolatile?
  BOAmP.hydrophilic=false; // Does the species condense on the aqueous phase?
  BOAmP.hydrophobic=true;  // Does the species condense on the organic phase?
  BOAmP.kp_from_experiment=true;  // Use experimental partitioning constant at Tref?
  BOAmP.kp_experiment=0.040;       // Value of the experimental partitioning constant at Tref?
  BOAmP.deltaH=91.0;     // Enthalpy of vaporization (kJ/mol)
  BOAmP.Tref=298;
  BOAmP.is_organic=true;  // Is the compound organic?
  BOAmP.compute_gamma_org=true;  // Compute the activity coefficients of the organic phase for this compound?
  BOAmP.compute_gamma_aq=true;  // Compute the activity coefficients of the aqueous phase for this compound
  BOAmP.Koligo_org=0.0;         //oligomeriation constant in the organic phase
  BOAmP.rho=1300.0;
  BOAmP.is_monomer=false;
  BOAmP.rion=false;
  //BOAmP.KDiffusion_air=1.0e-5;
  //  BOAmP.accomodation_coefficient=alpha;
  BOAmP.viscosity=1.68e12;  
  BOAmP.is_solid=false;
  BOAmP.is_generic=false;

  //Group: if no functionnal group in the species use the default species
  //for the computation of activity coefficients 
  
  double group_tmp_boamp [] = {0.0,0.0,0.0,0.0, // group C
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


  size = sizeof(group_tmp_boamp)/sizeof(double);
  assert(size == 60);
  for(int i = 0; i < size; ++i)
    BOAmP.groups[i] = group_tmp_boamp[i];

  // Search the species name in the aerosol species list 
  // and add the species if its name matches with
  // the given list.
  add_species_ssh(surrogate, BOAmP, species_list_aer, molecular_weight_aer,
                  accomodation_coefficient,diffusion_coef,saturation_vapor_pressure,enthalpy_vaporization, henry, t_ref, mass_density, species_part,nlayer,i_hydrophilic,
                  N_inert, N_inorganic);
  
  species BOAhP;
  BOAhP.name="BOAhP";
  BOAhP.is_inorganic_precursor=false;
  BOAhP.nonvolatile=false;  // Is the compound nonvolatile?
  BOAhP.hydrophilic=false; // Does the species condense on the aqueous phase?
  BOAhP.hydrophobic=true;  // Does the species condense on the organic phase?
  BOAhP.kp_from_experiment=true;  // Use experimental partitioning constant at Tref?
  BOAhP.kp_experiment=0.00023;       // Value of the experimental partitioning constant at Tref?
  BOAhP.deltaH=79.0;     // Enthalpy of vaporization (kJ/mol)
  BOAhP.Tref=298;
  BOAhP.is_organic=true;  // Is the compound organic?
  BOAhP.compute_gamma_org=true;  // Compute the activity coefficients of the organic phase for this compound?
  BOAhP.compute_gamma_aq=true;  // Compute the activity coefficients of the aqueous phase for this compound?
  BOAhP.Koligo_org=0.0;         //oligomeriation constant in the organic phase
  BOAhP.rho=1300.0;
  BOAhP.is_monomer=false;
  BOAhP.rion=false;
  //BOAhP.KDiffusion_air=1.0e-5;
  //  BOAhP.accomodation_coefficient=alpha;
  BOAhP.viscosity=1.68e12;  
  BOAhP.is_solid=false;
  BOAhP.is_generic=false;

  //Group: if no functionnal group in the species use the default species
  //for the computation of activity coefficients 
  
  double group_tmp_boahp [] = {0.0,0.0,0.0,0.0, // group C
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

  size = sizeof(group_tmp_boahp)/sizeof(double);
  assert(size == 60);
  for(int i = 0; i < size; ++i)
    BOAhP.groups[i] = group_tmp_boahp[i];

  // Search the species name in the aerosol species list 
  // and add the species if its name matches with
  // the given list.
  add_species_ssh(surrogate, BOAhP, species_list_aer, molecular_weight_aer,
                  accomodation_coefficient,diffusion_coef,saturation_vapor_pressure,enthalpy_vaporization, henry, t_ref, mass_density, species_part,nlayer,i_hydrophilic,
                  N_inert, N_inorganic);
  
  species BSOAlP;
  BSOAlP.name="BSOAlP";
  BSOAlP.is_inorganic_precursor=false;
  BSOAlP.nonvolatile=false;  // Is the compound nonvolatile?
  BSOAlP.hydrophilic=false; // Does the species condense on the aqueous phase?
  BSOAlP.hydrophobic=true;  // Does the species condense on the organic phase?
  BSOAlP.kp_from_experiment=true;  // Use experimental partitioning constant at Tref?
  BSOAlP.kp_experiment=1830.0;       // Value of the experimental partitioning constant at Tref?
  BSOAlP.deltaH=106.0;     // Enthalpy of vaporization (kJ/mol)
  BSOAlP.Tref=298;
  BSOAlP.is_organic=true;  // Is the compound organic?
  BSOAlP.compute_gamma_org=true;  // Compute the activity coefficients of the organic phase for this compound?
  BSOAlP.compute_gamma_aq=true;  // Compute the activity coefficients of the aqueous phase for this compound?
  BSOAlP.Koligo_org=0.0;         //oligomeriation constant in the organic phase
  BSOAlP.rho=1300.0;
  BSOAlP.is_monomer=false;
  BSOAlP.rion=false;
  //BSOAlP.KDiffusion_air=1.0e-5;
  //  BSOAlP.accomodation_coefficient=alpha;
  BSOAlP.viscosity=1.68e12;  
  BSOAlP.is_solid=false;
  BSOAlP.is_generic=false;
 
  //Group: if no functionnal group in the species use the default species
  //for the computation of activity coefficients 
  
  double group_tmp_bsoalp [] = {0.0,0.0,0.0,0.0, // group C
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


  size = sizeof(group_tmp_bsoalp)/sizeof(double);
  assert(size == 60);
  for(int i = 0; i < size; ++i)
    BSOAlP.groups[i] = group_tmp_bsoalp[i];

  // Search the species name in the aerosol species list 
  // and add the species if its name matches with
  // the given list.
  add_species_ssh(surrogate, BSOAlP, species_list_aer, molecular_weight_aer,
		  accomodation_coefficient,diffusion_coef,saturation_vapor_pressure,enthalpy_vaporization, henry, t_ref, mass_density, species_part,nlayer,i_hydrophilic,
  	          N_inert, N_inorganic);
  
  species BSOAmP;
  BSOAmP.name="BSOAmP";
  BSOAmP.is_inorganic_precursor=false;
  BSOAmP.nonvolatile=false;  // Is the compound nonvolatile?
  BSOAmP.hydrophilic=false; // Does the species condense on the aqueous phase?
  BSOAmP.hydrophobic=true;  // Does the species condense on the organic phase?
  BSOAmP.kp_from_experiment=true;  // Use experimental partitioning constant at Tref?
  BSOAmP.kp_experiment=4.0;       // Value of the experimental partitioning constant at Tref?
  BSOAmP.deltaH=91.0;     // Enthalpy of vaporization (kJ/mol)
  BSOAmP.Tref=298;
  BSOAmP.is_organic=true;  // Is the compound organic?
  BSOAmP.compute_gamma_org=true;  // Compute the activity coefficients of the organic phase for this compound?
  BSOAmP.compute_gamma_aq=true;  // Compute the activity coefficients of the aqueous phase for this compound?
  BSOAmP.Koligo_org=0.0;         //oligomeriation constant in the organic phase
  BSOAmP.rho=1300.0;
  BSOAmP.is_monomer=false;
  BSOAmP.rion=false;
  //BSOAmP.KDiffusion_air=1.0e-5;
  //  BSOAmP.accomodation_coefficient=alpha;
  BSOAmP.viscosity=1.68e12;
  BSOAmP.is_solid=false;
  BSOAmP.is_generic=false;

  //Group: if no functionnal group in the species use the default species
  //for the computation of activity coefficients 
  
  double group_tmp_bsoamp [] = {0.0,0.0,0.0,0.0, // group C
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


  size = sizeof(group_tmp_bsoamp)/sizeof(double);
  assert(size == 60);
  for(int i = 0; i < size; ++i)
    BSOAmP.groups[i] = group_tmp_bsoamp[i];

  // Search the species name in the aerosol species list 
  // and add the species if its name matches with
  // the given list.
  add_species_ssh(surrogate, BSOAmP, species_list_aer, molecular_weight_aer,
		  accomodation_coefficient,diffusion_coef,saturation_vapor_pressure,enthalpy_vaporization, henry, t_ref, mass_density, species_part,nlayer,i_hydrophilic,
	          N_inert, N_inorganic);
  
  species BSOAhP;
  BSOAhP.name="BSOAhP";
  BSOAhP.is_inorganic_precursor=false;
  BSOAhP.nonvolatile=false;  // Is the compound nonvolatile?
  BSOAhP.hydrophilic=false; // Does the species condense on the aqueous phase?
  BSOAhP.hydrophobic=true;  // Does the species condense on the organic phase?
  BSOAhP.kp_from_experiment=true;  // Use experimental partitioning constant at Tref?
  BSOAhP.kp_experiment=0.023;       // Value of the experimental partitioning constant at Tref?
  BSOAhP.deltaH=79.0;     // Enthalpy of vaporization (kJ/mol)
  BSOAhP.Tref=298;
  BSOAhP.is_organic=true;  // Is the compound organic?
  BSOAhP.compute_gamma_org=true;  // Compute the activity coefficients of the organic phase for this compound?
  BSOAhP.compute_gamma_aq=true;  // Compute the activity coefficients of the aqueous phase for this compound?
  BSOAhP.Koligo_org=0.0;         //oligomeriation constant in the organic phase
  BSOAhP.rho=1300.0;
  BSOAhP.is_monomer=false;
  BSOAhP.rion=false;
  //BSOAhP.KDiffusion_air=1.0e-5;
  //  BSOAhP.accomodation_coefficient=alpha;
  BSOAhP.viscosity=1.68e12;  
  BSOAhP.is_solid=false;
  BSOAhP.is_generic=false;

  //Group: if no functionnal group in the species use the default species
  //for the computation of activity coefficients 
  
  double group_tmp_bsoahp [] = {0.0,0.0,0.0,0.0, // group C
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


  size = sizeof(group_tmp_bsoahp)/sizeof(double);
  assert(size == 60);
  for(int i = 0; i < size; ++i)
    BSOAhP.groups[i] = group_tmp_bsoahp[i];

  // Search the species name in the aerosol species list 
  // and add the species if its name matches with
  // the given list.
  add_species_ssh(surrogate, BSOAhP, species_list_aer, molecular_weight_aer,
		  accomodation_coefficient,diffusion_coef,saturation_vapor_pressure,enthalpy_vaporization, henry, t_ref, mass_density, species_part,nlayer,i_hydrophilic,
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
  Monomer.compute_gamma_org=true;  // Compute the activity coefficients of the organic phase for this compound?
  Monomer.compute_gamma_aq=true;  // Compute the activity coefficients of the aqueous phase for this compound?
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
  double group_tmp_Monomer [] = {0.0,1.0,2.0,0.0, // group C
				 0.0,0.0,0.0,0.0, //group C[OH]
				 0.0,0.0,0.0,0.0, //group Calcohol
				 0.0,0.0,0.0,0.0, //group Calcohol-tail
				 0.0,0.0,0.0,0.0,0.0, //group C=C
				 0.0,0.0, //group aromatic carbon (AC)
				 0.0,0.0,0.0, // group //AC-C
				 0.0,  //group OH
				 0.0, //group H2O
				 0.0, //group ACOH
				 0.0,2.0, //group ketone
				 0.0,   //group aldehyde  
				 0.0,0.0, //group ester
				 0.0,0.0,0.0, //group ether 
				 0.0,  //group acid
				 0.0,   //group ACNO2
				 0.0,0.0,0.0, //group NO3
				 0.0,2.0,0.0, //group CO-OH
				 0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0, //group CO-OC
				 0.0,  //group PAN
				 1.0,  //group CO-OOH
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
		  accomodation_coefficient,diffusion_coef,saturation_vapor_pressure,enthalpy_vaporization, henry, t_ref, mass_density, species_part,nlayer,i_hydrophilic,
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
  Dimer.compute_gamma_org=true;  // Compute the activity coefficients of the organic phase for this compound?
  Dimer.compute_gamma_aq=true;  // Compute the activity coefficients of the aqueous phase for this compound?
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
  double group_tmp_dimer [] = {1.0,4.0,4.0,0.0, // group C
			       0.0,0.0,0.0,0.0, //group C[OH]
			       0.0,0.0,0.0,0.0, //group Calcohol
			       0.0,0.0,0.0,0.0, //group Calcohol-tail
			       0.0,0.0,0.0,0.0,0.0, //group C=C
			       0.0,0.0, //group aromatic carbon (AC)
			       0.0,0.0,0.0, // group //AC-C
			       1.0,  //group OH
			       0.0, //group H2O
			       0.0, //group ACOH
			       0.0,2.0, //group ketone
			       1.0,   //group aldehyde  
			       0.0,0.0, //group ester
			       0.0,0.0,0.0, //group ether 
			       0.0,  //group acid
			       0.0,   //group ACNO2
			       0.0,0.0,0.0, //group NO3
			       0.0,1.0,0.0, //group CO-OH
			       0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0, //group CO-OC
			       0.0,  //group PAN
			       1.0,  //group CO-OOH
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
		  accomodation_coefficient,diffusion_coef,saturation_vapor_pressure,enthalpy_vaporization, henry, t_ref, mass_density, species_part,nlayer,i_hydrophilic,
		  N_inert, N_inorganic);

  species BiA3D;
  BiA3D.name="BiA3D";
  BiA3D.is_inorganic_precursor=false;
  BiA3D.Psat_ref=3.25e-7; // Saturation vapor pressure at Tref (torr)
  BiA3D.kp_from_experiment=false;  // Use experimental partitioning constant at Tref?
  BiA3D.Tref=298;         // Temperature of reference (K)
  BiA3D.deltaH=109.0;     // Enthalpy of vaporization (kJ/mol)
  BiA3D.Henry=0.;     // Henry's law constant at Tref (M/atm)
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
		  accomodation_coefficient,diffusion_coef,saturation_vapor_pressure,enthalpy_vaporization, henry, t_ref, mass_density, species_part,nlayer,i_hydrophilic,
		  N_inert, N_inorganic);

  species ACIDMAL;
  ACIDMAL.name="ACIDMAL";
  ACIDMAL.is_inorganic_precursor=false;
  ACIDMAL.Psat_ref=4.59e-8; // Saturation vapor pressure at Tref (torr)
  ACIDMAL.Tref=298;
  ACIDMAL.deltaH= 50.0;     // Enthalpy of vaporization (kJ/mol)
  ACIDMAL.Henry=0.0;     // Henry's law constant at Tref (M/atm)
  ACIDMAL.aq_type="none"; // "none","diacid","monoacid" or "aldehyde"
  ACIDMAL.hydrophilic=true; // Does the species condense on the aqueous phase?
  ACIDMAL.hydrophobic=false;  // Does the species condense on the organic phase?
  ACIDMAL.nonvolatile=false;  // Is the compound nonvolatile?
  ACIDMAL.kp_from_experiment= false;  // Use experimental partitioning constant at Tref?
  ACIDMAL.kp_experiment= 2.56;       // Value of the experimental partitioning constant at Tref?
  ACIDMAL.is_organic=true;  // Is the compound organic?
  ACIDMAL.compute_gamma_org=true;  // Compute the activity coefficients of the organic phase for this compound?
  ACIDMAL.compute_gamma_aq=true;  // Compute the activity coefficients of the aqueous phase for this compound
  ACIDMAL.aq_type="none"; // "none","diacid","monoacid" or "aldehyde"
  ACIDMAL.Koligo_org=0.0;
  //Parameters for the oligomerization of aldehyde in the aqueous phase as BiA0D:  
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
		  accomodation_coefficient,diffusion_coef,saturation_vapor_pressure,enthalpy_vaporization, henry, t_ref, mass_density, species_part,nlayer,i_hydrophilic,
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
		  accomodation_coefficient,diffusion_coef,saturation_vapor_pressure,enthalpy_vaporization, henry, t_ref, mass_density, species_part,nlayer,i_hydrophilic,
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
		  accomodation_coefficient,diffusion_coef,saturation_vapor_pressure,enthalpy_vaporization, henry, t_ref, mass_density, species_part,nlayer,i_hydrophilic,
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
		  accomodation_coefficient,diffusion_coef,saturation_vapor_pressure,enthalpy_vaporization, henry, t_ref, mass_density, species_part,nlayer,i_hydrophilic,
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
		  accomodation_coefficient,diffusion_coef,saturation_vapor_pressure,enthalpy_vaporization, henry, t_ref, mass_density, species_part,nlayer,i_hydrophilic,
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
		  accomodation_coefficient,diffusion_coef,saturation_vapor_pressure,enthalpy_vaporization, henry, t_ref, mass_density, species_part,nlayer,i_hydrophilic,
		  N_inert, N_inorganic);

  species C7H9O9;
  C7H9O9.name="C7H9O9";
  C7H9O9.is_inorganic_precursor=false;
  C7H9O9.Psat_ref=5.9e-13; // Saturation vapor pressure at Tref (torr)
  C7H9O9.Tref=298;         // Temperature of reference (K)
  C7H9O9.deltaH=100.0;     // Enthalpy of vaporization (kJ/mol)
//C7H9O9.aq_type="aldehyde";  // "none","diacid","monoacid" or "aldehyde"
  C7H9O9.aq_type="none";  // "none","diacid","monoacid" or "aldehyde"
  C7H9O9.Henry=0.0;       //If the Henry's law constant is set to zero, the model compute the Henry's law constant from the saturation vapour pressure and the activity coefficients at infinite dilution given by unifac
  C7H9O9.hydrophilic=true;  // Does the species condense on the aqueous phase?
  C7H9O9.hydrophobic=true;  // Does the species condense on the organic phase?
  C7H9O9.nonvolatile=false; // Is the compound nonvolatile?
  C7H9O9.kp_from_experiment=false;  // Use experimental partitioning constant at Tref?
  C7H9O9.is_organic=true;  // Is the compound organic?
  C7H9O9.compute_gamma_org=true;  // Compute the activity coefficients of the organic phase for this compound?
  C7H9O9.compute_gamma_aq=true;  // Compute the activity coefficients of the aqueous phase for this compound?
  C7H9O9.Koligo_org=0.0;      //oligomeriation constant in the organic phase
  C7H9O9.rho=1400.0;
  C7H9O9.is_monomer=false;
  C7H9O9.rion=false;
  //Monomer.KDiffusion_air=1.0e-5;
  //  Monomer.accomodation_coefficient=alpha;
  C7H9O9.viscosity=1.68e12;
  C7H9O9.is_solid=false;
  C7H9O9.is_generic=false;
  
  //Group: if no functionnal group in the species use the default species
  //for the computation of activity coefficients
  //
  double group_tmp_C7H9O9 [] = {0.0,0.0,0.0,0.0, // group C
				 0.0,0.0,0.0,0.0, //group C[OH]
				 0.0,0.0,0.0,0.0, //group Calcohol
				 0.0,0.0,0.0,0.0, //group Calcohol-tail
				 0.0,0.0,0.0,0.0,0.0, //group C=C
				 0.0,0.0, //group aromatic carbon (AC)
				 0.0,0.0,0.0, // group //AC-C
				 0.0,  //group OH
				 0.0, //group H2O
				 0.0, //group ACOH
				 1.0,1.0, //group ketone
				 1.0,   //group aldehyde  
				 0.0,0.0, //group ester
				 0.0,0.0,0.0, //group ether 
				 0.0,  //group acid
				 0.0,   //group ACNO2
				 0.0,0.0,0.0, //group NO3
				 0.0,3.0,0.0, //group CO-OH
				 0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0, //group CO-OC
				 0.0,  //group PAN
				 0.0,  //group CO-OOH
	  		         0.0,  //group O=COC=O
			         0.0,0.0,0.0}; //group CHxNO2

  size = sizeof(group_tmp_C7H9O9)/sizeof(double);
  assert(size == 60);
	
  for(int i = 0; i < size; ++i)
    C7H9O9.groups[i] = group_tmp_C7H9O9[i];

  // Search the species name in the aerosol species list 
  // and add the species if its name matches with
  // the given list.
  add_species_ssh(surrogate, C7H9O9, species_list_aer, molecular_weight_aer,
		  accomodation_coefficient,diffusion_coef,saturation_vapor_pressure,enthalpy_vaporization, henry, t_ref, mass_density, species_part,nlayer,i_hydrophilic,
		  N_inert, N_inorganic);
  

  add_generic_species_ssh(config, surrogate, species_list_aer, molecular_weight_aer, accomodation_coefficient, diffusion_coef,
			  aerosol_type, species_smiles, saturation_vapor_pressure, enthalpy_vaporization, henry, t_ref, mass_density,  
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
		  accomodation_coefficient,diffusion_coef,saturation_vapor_pressure,enthalpy_vaporization, henry, t_ref, mass_density, species_part,nlayer,i_hydrophilic,
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
  NO3.MM=63.0;           // Molar mass (g/mol)
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

  species OH;
  OH.name="OH";
  OH.is_inorganic_precursor=false;
  OH.MM=17.0;           // Molar mass (g/mol)
  OH.is_organic=false;  // Is the compound organic?
  OH.hydrophilic=true; // Does the species condense on the aqueous phase?
  OH.hydrophobic=false;  // Does the species condense on the organic phase?
  OH.nonvolatile=false;
  OH.kp_from_experiment=false;
  OH.compute_gamma_org=false;  // Compute the activity coefficients of the organic phase for this compound?
  OH.compute_gamma_aq=false;  // Compute the activity coefficients of the aqueous phase for this compound
  OH.charge=-1.0;
  OH.rho=1300.0; //1400.0;
  OH.index_ion_aiomfac=14;
  OH.KDiffusion_air=1.0e-5;
  OH.accomodation_coefficient=0.5;
  OH.viscosity=1.0;
  OH.soap_ind = -1;
  OH.soap_ind_aero = -1;
  OH.is_solid=false;
  OH.is_monomer=false;
  OH.rion=false;
  surrogate.push_back(OH);

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

  species HCO3;
  HCO3.name="HCO3";
  HCO3.is_inorganic_precursor=false;
  HCO3.MM=61.0;           // Molar mass (g/mol)
  HCO3.is_organic=false;  // Is the compound organic?
  HCO3.hydrophilic=true; // Does the species condense on the aqueous phase?
  HCO3.hydrophobic=false;  // Does the species condense on the organic phase?
  HCO3.nonvolatile=false;
  HCO3.kp_from_experiment=false;
  HCO3.compute_gamma_org=false;  // Compute the activity coefficients of the organic phase for this compound?
  HCO3.compute_gamma_aq=false;  // Compute the activity coefficients of the aqueous phase for this compound
  HCO3.charge=-1.0;
  HCO3.index_ion_aiomfac=15;
  HCO3.rho=1300.0; // 1400.0;
  HCO3.KDiffusion_air=1.0e-5;
  HCO3.accomodation_coefficient=0.5;
  HCO3.viscosity=1.0;
  HCO3.soap_ind = -1;
  HCO3.soap_ind_aero = -1;
  HCO3.is_solid=false;
  HCO3.is_monomer=false;
  HCO3.rion=false;
  for (int i = 0; i < nsp; ++i)
    if (species_list_aer[i].substr(1,-1) == "CO3")
      surrogate.push_back(HCO3);

  species CO3;
  CO3.name="CO3";
  CO3.is_inorganic_precursor=false;
  CO3.MM=60.0;           // Molar mass (g/mol)
  CO3.is_organic=false;  // Is the compound organic?
  CO3.hydrophilic=true; // Does the species condense on the aqueous phase?
  CO3.hydrophobic=false;  // Does the species condense on the organic phase?
  CO3.nonvolatile=false;
  CO3.kp_from_experiment=false;
  CO3.compute_gamma_org=false;  // Compute the activity coefficients of the organic phase for this compound?
  CO3.compute_gamma_aq=false;  // Compute the activity coefficients of the aqueous phase for this compound
  CO3.charge=-2.0;
  CO3.index_ion_aiomfac=16;
  CO3.rho=1300.0; // 1400.0;
  CO3.KDiffusion_air=1.0e-5;
  CO3.accomodation_coefficient=0.5;
  CO3.viscosity=1.0;
  CO3.soap_ind = -1;
  CO3.soap_ind_aero = -1;
  CO3.is_solid=false;
  CO3.is_monomer=false;
  CO3.rion=false;

  for (int i = 0; i < nsp; ++i)
    if (species_list_aer[i].substr(1,-1) == "CO3")
      {
        CO3.soap_ind = i;
        CO3.soap_ind_aero = i;
	surrogate.push_back(CO3);
	with_co3=1;
      }

  species Ca;
  Ca.name="Ca";
  Ca.is_inorganic_precursor=false;
  Ca.MM=40.078;           // Molar mass (g/mol)
  Ca.is_organic=false;  // Is the compound organic?
  Ca.hydrophilic=true; // Does the species condense on the aqueous phase?
  Ca.hydrophobic=false;  // Does the species condense on the organic phase?
  Ca.nonvolatile=false;
  Ca.kp_from_experiment=false;
  Ca.compute_gamma_org=false;  // Compute the activity coefficients of the organic phase for this compound?
  Ca.compute_gamma_aq=false;  // Compute the activity coefficients of the aqueous phase for this compound
  Ca.charge=2.0;
  Ca.index_ion_aiomfac=6;
  Ca.rho=1300.0; // 1400.0;
  Ca.KDiffusion_air=1.0e-5;
  Ca.accomodation_coefficient=0.5;
  Ca.viscosity=1.0;
  Ca.soap_ind = -1;
  Ca.soap_ind_aero = -1;
  Ca.is_solid=false;
  Ca.is_monomer=false;
  Ca.rion=false;

  for (int i = 0; i < nsp; ++i)
    if (species_list_aer[i].substr(1,-1) == "Ca")
      {
        Ca.soap_ind = i;
        Ca.soap_ind_aero = i;
	surrogate.push_back(Ca);
	with_ca=1;
      }

  species Mg;
  Mg.name="Mg";
  Mg.is_inorganic_precursor=false;
  Mg.MM=24.305;           // Molar mass (g/mol)
  Mg.is_organic=false;  // Is the compound organic?
  Mg.hydrophilic=true; // Does the species condense on the aqueous phase?
  Mg.hydrophobic=false;  // Does the species condense on the organic phase?
  Mg.nonvolatile=false;
  Mg.kp_from_experiment=false;
  Mg.compute_gamma_org=false;  // Compute the activity coefficients of the organic phase for this compound?
  Mg.compute_gamma_aq=false;  // Compute the activity coefficients of the aqueous phase for this compound
  Mg.charge=2.0;
  Mg.index_ion_aiomfac=5;
  Mg.rho=1300.0; // 1400.0;
  Mg.KDiffusion_air=1.0e-5;
  Mg.accomodation_coefficient=0.5;
  Mg.viscosity=1.0;
  Mg.soap_ind = -1;
  Mg.soap_ind_aero = -1;
  Mg.is_solid=false;
  Mg.is_monomer=false;
  Mg.rion=false;

  for (int i = 0; i < nsp; ++i)
    if (species_list_aer[i].substr(1,-1) == "Mg")
      {
        Mg.soap_ind = i;
        Mg.soap_ind_aero = i;
	surrogate.push_back(Mg);
      }

  species K;
  K.name="K";
  K.is_inorganic_precursor=false;
  K.MM=39.0983;           // Molar mass (g/mol)
  K.is_organic=false;  // Is the compound organic?
  K.hydrophilic=true; // Does the species condense on the aqueous phase?
  K.hydrophobic=false;  // Does the species condense on the organic phase?
  K.nonvolatile=false;
  K.kp_from_experiment=false;
  K.compute_gamma_org=false;  // Compute the activity coefficients of the organic phase for this compound?
  K.compute_gamma_aq=false;  // Compute the activity coefficients of the aqueous phase for this compound
  K.charge=1.0;
  K.index_ion_aiomfac=3;
  K.rho=1300.0; // 1400.0;
  K.KDiffusion_air=1.0e-5;
  K.accomodation_coefficient=0.5;
  K.viscosity=1.0;
  K.soap_ind = -1;
  K.soap_ind_aero = -1;
  K.is_solid=false;
  K.is_monomer=false;
  K.rion=false;

  for (int i = 0; i < nsp; ++i)
    if (species_list_aer[i].substr(1,-1) == "K")
      {
        K.soap_ind = i;
        K.soap_ind_aero = i;
	surrogate.push_back(K);
      }

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

  species CO2;
  CO2.name="CO2";
  CO2.MM=44.0;
  CO2.nonvolatile=false;
  CO2.is_organic=false;
  CO2.is_inorganic_precursor=true;
  CO2.Henry=2.5e3;     // Enthalpy of vaporization (kJ/mol)
  CO2.deltaH=0.0;     // Henry's law constant at Tref (M/atm)
  CO2.hydrophilic=true; // Does the species condense on the aqueous phase?
  CO2.hydrophobic=false;  // Does the species condense on the organic phase?
  CO2.compute_gamma_org=false;  // Compute the activity coefficients of the organic phase for this compound?
  CO2.compute_gamma_aq=false;  // Compute the activity coefficients of the aqueous phase for this compound
  //  CO2.accomodation_coefficient=0.5;
  CO2.KDiffusion_air=1.72e-5;
  CO2.viscosity=1.0;
  CO2.is_solid=false;
  
  // Find the number in the aerosol species list
  CO2.soap_ind = -1;
  CO2.soap_ind_aero = -1;       
  CO2.accomodation_coefficient = 0.1;
  /*
    for (int i = 0; i < nsp; ++i)
    if (species_list_aer[i].substr(1,-1) == "H2CO3")
    {	   
    CO2.accomodation_coefficient = accomodation_coefficient[i];
    CO2.KDiffusion_air = diffusion_coef[i];
    }*/
  for (int i = 0; i < nsp; ++i)
    if (species_list_aer[i].substr(1,-1) == "CO3")
      {
	CO2.soap_ind=i;
	surrogate.push_back(CO2);
      }

  if (with_ca==1 and with_co3==1)
    {
      if (config.compute_inorganic)
	config.solids=true;
      
      species CaCO3;
      CaCO3.name="CaCO3";
      CaCO3.MM=100.078;
      CaCO3.nonvolatile=false;
      CaCO3.is_organic=false;
      CaCO3.is_solid=true;
      CaCO3.is_inorganic_precursor=true;
      //  CaCO3.Henry=2.5e3;     // Enthalpy of vaporization (kJ/mol)
      //  CaCO3.deltaH=0.0;     // Henry's law constant at Tref (M/atm)
      CaCO3.hydrophilic=false; // Does the species condense on the aqueous phase?
      CaCO3.hydrophobic=false;  // Does the species condense on the organic phase?
      CaCO3.compute_gamma_org=false;  // Compute the activity coefficients of the organic phase for this compound?
      CaCO3.compute_gamma_aq=false;  // Compute the activity coefficients of the aqueous phase for this compound
      CaCO3.KDiffusion_air=1.72e-5;
      CaCO3.accomodation_coefficient=0.5;
      CaCO3.Ksol=2.25e-8; //11.73003492;   //4.199e−17*2.1e5*57.64*1.805e−5/1.01e-14;
      CaCO3.deltaH=0.; //8.16364555; //-74.735+29.17-1.5+13.79+22.52;
      CaCO3.dCp=0.; //-71.6059517; //6.025+16.83+26.92-5.39-26.92; //
      CaCO3.nion=2;
      CaCO3.ion1="Ca";
      CaCO3.pion1=1;
      CaCO3.ion2="CO3";
      CaCO3.pion2=1;
      CaCO3.viscosity=1.;
      surrogate.push_back(CaCO3);
    }

}

#endif
