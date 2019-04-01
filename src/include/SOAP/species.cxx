#ifndef SPECIES_CXX
#define SPECIES_CXX


#include <vector>
#include "species.h"
using namespace soap;
using namespace std;


void add_species( vector<species>& surrogate, species current_species, 
                  vector<string> species_list_aer)
{

  int nsp = species_list_aer.size();

  // Find the number in the aerosol species list
  current_species.soap_ind = -1;
  for (int i = 0; i < nsp; ++i)
    if (species_list_aer[i].substr(1,-1) == current_species.name)
      current_species.soap_ind = i;

  // if (current_species.soap_ind == -1)
  //   cout << "Warning: " << current_species.name << " is not found in " <<
  //     "the aerosol species list.\n";
  // else
  //   surrogate.push_back(current_species);

  if (current_species.soap_ind != -1)
    surrogate.push_back(current_species);

}
  
void creation_species( vector<species>& surrogate, vector<string> species_list_aer)
{
  int nsp = species_list_aer.size();
  double alpha = 0.5; //0.01; // accommodation coefficient

  species BiA2D;
  BiA2D.name="BiA2D";
  BiA2D.is_inorganic_precursor=false;
  BiA2D.Psat_ref=1.43e-7; // Saturation vapor pressure at Tref (torr)
  BiA2D.kp_from_experiment=false;  // Use experimental partitioning constant at Tref?
  BiA2D.Tref=298;         // Temperature of reference (K)
  BiA2D.MM=186;           // Molar mass (g/mol)
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
  BiA2D.KDiffusion_air=1.0e-5;
  BiA2D.accomodation_coefficient=alpha;
  BiA2D.viscosity=1.68e12;

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
							   0.0,0.0,0.0}; //group CO-OH
  
  int size = sizeof(group_tmp_bia2d)/sizeof(double);
  assert(size = 45);
	
  for(int i = 0; i < size; ++i)
	BiA2D.groups[i] = group_tmp_bia2d[i];

  // Search the species name in the aerosol species list 
  // and add the species if its name matches with
  // the given list.
  add_species(surrogate, BiA2D, species_list_aer);


  /* ==== BiA1D ==== */ 

  species BiA1D;
  BiA1D.name="BiA1D";
  BiA1D.is_inorganic_precursor=false;
  BiA1D.Psat_ref=2.17e-7; // Saturation vapor pressure at Tref (torr)
  BiA1D.Tref=298;         // Temperature of reference (K)
  BiA1D.MM=170.0;           // Molar mass (g/mol)
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
  BiA1D.KDiffusion_air=1.0e-5;
  BiA1D.accomodation_coefficient=alpha;
  BiA1D.viscosity=1.68e12;

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
							   0.0,0.0,0.0}; //group CO-OH
  
  size = sizeof(group_tmp_bia1d)/sizeof(double);
  assert(size = 45);
	
  for(int i = 0; i < size; ++i)
	BiA1D.groups[i] = group_tmp_bia1d[i];

  // Search the species name in the aerosol species list 
  // and add the species if its name matches with
  // the given list.
  add_species(surrogate, BiA1D, species_list_aer);

  species BiA0D;
  BiA0D.name="BiA0D";
  BiA0D.is_inorganic_precursor=false;
  BiA0D.Psat_ref=2.7e-4; // Saturation vapor pressure at Tref (torr)
  BiA0D.Tref=298;         // Temperature of reference (K)
  BiA0D.MM=168.0;           // Molar mass (g/mol)
  BiA0D.deltaH=50.0;     // Enthalpy of vaporization (kJ/mol)
  BiA0D.Henry=1.98e6;     // Henry's law constant at Tref (M/atm)
  BiA0D.aq_type="aldehyde"; // "none","diacid","monoacid" or "aldehyde"
  BiA0D.Koligo_org=0.0;    //oligomeriation constant in the organic phase
  BiA0D.hydrophilic=true;   // Does the species condense on the aqueous phase?
  BiA0D.hydrophobic=false;  // Does the species condense on the organic phase?
  BiA0D.nonvolatile=false; // Is the compound nonvolatile?
  BiA0D.kp_from_experiment=false;  // Use experimental partitioning constant at Tref?
  BiA0D.is_organic=true;  // Is the compound organic?
  BiA0D.compute_gamma_org=true;  // Compute the activity coefficients of the organic phase for this compound?
  BiA0D.compute_gamma_aq=true;  // Compute the activity coefficients of the aqueous phase for this compound?
  //Parameters for the oligomerization of aldehyde in the aqueous phase:
  BiA0D.Koligo_aq=0.1;     
  BiA0D.pHref=6.0;
  BiA0D.beta=1.91;
  BiA0D.rho=1300.0;
  BiA0D.KDiffusion_air=1.0e-5;
  BiA0D.accomodation_coefficient=alpha;
  BiA0D.viscosity=1.68e12;
  
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
							   0.0,0.0,0.0}; //group CO-OH

  size = sizeof(group_tmp_bia0d)/sizeof(double);
  assert(size = 45);
	
  for(int i = 0; i < size; ++i)
	BiA0D.groups[i] = group_tmp_bia0d[i];

  // Search the species name in the aerosol species list 
  // and add the species if its name matches with
  // the given list.
  add_species(surrogate, BiA0D, species_list_aer);
  
  
  species BiMT;
  BiMT.name="BiMT";
  BiMT.is_inorganic_precursor=false;
  BiMT.Psat_ref=1.45e-6; // Saturation vapor pressure at Tref (torr)
  BiMT.Tref=298;         // Temperature of reference (K)
  BiMT.MM=136;           // Molar mass (g/mol)
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
  BiMT.KDiffusion_air=1.0e-5;
  BiMT.accomodation_coefficient=alpha;
  BiMT.viscosity=1.68e12;
 
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
							  0.0,0.0,0.0}; //group CO-OH

  size = sizeof(group_tmp_bimt)/sizeof(double);
  assert(size == 45);
	
  for(int i = 0; i < size; ++i)
	BiMT.groups[i] = group_tmp_bimt[i];

  // Search the species name in the aerosol species list 
  // and add the species if its name matches with
  // the given list.
  add_species(surrogate, BiMT, species_list_aer);


  species BiPER;
  BiPER.name="BiPER";
  BiPER.is_inorganic_precursor=false;
  BiPER.Psat_ref=2.61e-6; // Saturation vapor pressure at Tref (torr)
  BiPER.Tref=298;         // Temperature of reference (K)
  BiPER.MM=168;           // Molar mass (g/mol)
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
  BiPER.KDiffusion_air=1.0e-5;
  BiPER.accomodation_coefficient=alpha;
  BiPER.viscosity=1.68e12;
  
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
							   0.0,2.0,0.0}; //group CO-OH

  size = sizeof(group_tmp_biper)/sizeof(double);
  assert(size == 45);
	
  for(int i = 0; i < size; ++i)
	BiPER.groups[i] = group_tmp_biper[i];

  // Search the species name in the aerosol species list 
  // and add the species if its name matches with
  // the given list.
  add_species(surrogate, BiPER, species_list_aer);


  species BiDER;
  BiDER.name="BiDER";
  BiDER.is_inorganic_precursor=false;
  BiDER.Psat_ref=4.10e-7; // Saturation vapor pressure at Tref (torr)
  BiDER.Tref=298;         // Temperature of reference (K)
  BiDER.MM=136;           // Molar mass (g/mol)
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
  BiDER.KDiffusion_air=1.0e-5;
  BiDER.accomodation_coefficient=alpha;
  BiDER.viscosity=1.68e12;
  
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
							   0.0,0.0,0.0}; //group CO-OH

  size = sizeof(group_tmp_bider)/sizeof(double);
  assert(size == 45);
	
  for(int i = 0; i < size; ++i)
	BiDER.groups[i] = group_tmp_bider[i];

  // Search the species name in the aerosol species list 
  // and add the species if its name matches with
  // the given list.
  add_species(surrogate, BiDER, species_list_aer);
  

  species BiMGA;
  BiMGA.name="BiMGA";
  BiMGA.is_inorganic_precursor=false;
  BiMGA.Psat_ref=1.40e-5; // Saturation vapor pressure at Tref (torr)
  BiMGA.Tref=298;         // Temperature of reference (K)
  BiMGA.MM=120;           // Molar mass (g/mol)
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
  BiMGA.KDiffusion_air=1.0e-5;
  BiMGA.accomodation_coefficient=alpha;
  BiMGA.viscosity=1.68e12;
  
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
							   0.0,0.0,0.0}; //group CO-OH

  size = sizeof(group_tmp_bimga)/sizeof(double);
  assert(size == 45);
	
  for(int i = 0; i < size; ++i)
	BiMGA.groups[i] = group_tmp_bimga[i];

  // Search the species name in the aerosol species list 
  // and add the species if its name matches with
  // the given list.
  add_species(surrogate, BiMGA, species_list_aer);


  species AnBlP;
  AnBlP.name="AnBlP";
  AnBlP.is_inorganic_precursor=false;
  AnBlP.Psat_ref=6.8e-8; // Saturation vapor pressure at Tref (torr)
  AnBlP.Tref=298;         // Temperature of reference (K)
  AnBlP.MM=167;           // Molar mass (g/mol)
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
  AnBlP.KDiffusion_air=1.0e-5;
  AnBlP.accomodation_coefficient=alpha;
  AnBlP.viscosity=1.68e12;
  
  
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
							   0.0,0.0,0.0}; //group CO-OH

  size = sizeof(group_tmp_anblp)/sizeof(double);
  assert(size == 45);
	
  for(int i = 0; i < size; ++i)
	AnBlP.groups[i] = group_tmp_anblp[i];

  // Search the species name in the aerosol species list 
  // and add the species if its name matches with
  // the given list.
  add_species(surrogate, AnBlP, species_list_aer);


  species AnBmP;
  AnBmP.name="AnBmP";
  AnBmP.is_inorganic_precursor=false;
  AnBmP.Psat_ref=8.4e-6; // Saturation vapor pressure at Tref (torr)
  AnBmP.Tref=298;         // Temperature of reference (K)
  AnBmP.MM=152;           // Molar mass (g/mol)
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
  AnBmP.KDiffusion_air=1.0e-5;
  AnBmP.accomodation_coefficient=alpha;
  AnBmP.viscosity=1.68e12;
  
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
							   0.0,0.0,0.0}; //group CO-OH

  size = sizeof(group_tmp_anbmp)/sizeof(double);
  assert(size == 45);
	
  for(int i = 0; i < size; ++i)
	AnBmP.groups[i] = group_tmp_anbmp[i];

  // Search the species name in the aerosol species list 
  // and add the species if its name matches with
  // the given list.
  add_species(surrogate, AnBmP, species_list_aer);


  species BiBlP;
  BiBlP.name="BiBlP";
  BiBlP.is_inorganic_precursor=false;
  BiBlP.Psat_ref=0.60e-9; // Saturation vapor pressure at Tref (torr)
  BiBlP.Tref=298;         // Temperature of reference (K)
  BiBlP.MM=298.0;           // Molar mass (g/mol)
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
  BiBlP.KDiffusion_air=1.0e-5;
  BiBlP.accomodation_coefficient=alpha;
  BiBlP.viscosity=1.68e12;  
  
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
							   0.0,0.0,0.0}; //group CO-OH

  size = sizeof(group_tmp_biblp)/sizeof(double);
  assert(size == 45);
	
  for(int i = 0; i < size; ++i)
	BiBlP.groups[i] = group_tmp_biblp[i];

  // Search the species name in the aerosol species list 
  // and add the species if its name matches with
  // the given list.
  add_species(surrogate, BiBlP, species_list_aer);

  species BiBmP;
  BiBmP.name="BiBmP";
  BiBmP.is_inorganic_precursor=false;
  BiBmP.Psat_ref=3.0e-7; // Saturation vapor pressure at Tref (torr)
  BiBmP.Tref=298;         // Temperature of reference (K)
  BiBmP.MM=236.0;           // Molar mass (g/mol)
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
  BiBmP.KDiffusion_air=1.0e-5;
  BiBmP.accomodation_coefficient=alpha;
  BiBmP.viscosity=1.68e12;  
  
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
							   0.0,0.0,0.0}; //group CO-OH

  size = sizeof(group_tmp_bibmp)/sizeof(double);
  assert(size == 45);
	
  for(int i = 0; i < size; ++i)
	BiBmP.groups[i] = group_tmp_bibmp[i];

  // Search the species name in the aerosol species list 
  // and add the species if its name matches with
  // the given list.
  add_species(surrogate, BiBmP, species_list_aer);
  
  species AnClP;
  AnClP.name="AnClP";
  AnClP.is_inorganic_precursor=false;
  AnClP.MM=167;           // Molar mass (g/mol)
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
  AnClP.KDiffusion_air=1.0e-5;
  AnClP.accomodation_coefficient=alpha;
  AnClP.viscosity=1.68e12;  

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
							   0.0,0.0,0.0}; //group CO-OH
  
  size = sizeof(group_tmp_anclp)/sizeof(double);
  assert(size == 45);
  
  for(int i = 0; i < size; ++i)
	AnClP.groups[i] = group_tmp_anclp[i];

  // Search the species name in the aerosol species list 
  // and add the species if its name matches with
  // the given list.
  add_species(surrogate, AnClP, species_list_aer);


  species BiNGA;
  BiNGA.name="BiNGA";
  BiNGA.is_inorganic_precursor=false;
  BiNGA.Psat_ref=1.40e-5; // Saturation vapor pressure at Tref (torr)
  BiNGA.Tref=298;         // Temperature of reference (K)
  BiNGA.MM=165;           // Molar mass (g/mol)
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
  BiNGA.KDiffusion_air=1.0e-5;
  BiNGA.accomodation_coefficient=alpha;
  BiNGA.viscosity=1.68e12;
  
  
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
							   0.0,0.0,0.0}; //group CO-OH

  size = sizeof(group_tmp_binga)/sizeof(double);
  assert(size == 45);
	
  for(int i = 0; i < size; ++i)
	BiNGA.groups[i] = group_tmp_binga[i];

  // Search the species name in the aerosol species list 
  // and add the species if its name matches with
  // the given list.
  add_species(surrogate, BiNGA, species_list_aer);


  species BiNIT3;
  BiNIT3.name="BiNIT3"; 
  BiNIT3.is_inorganic_precursor=false;  
  BiNIT3.Psat_ref=1.45e-6; // Saturation vapor pressure at Tref (torr)
  BiNIT3.Tref=298;         // Temperature of reference (K)
  BiNIT3.MM=272.0;           // Molar mass (g/mol)
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
  BiNIT3.KDiffusion_air=1.0e-5;
  BiNIT3.accomodation_coefficient=alpha;
  BiNIT3.viscosity=1.68e12;  
  
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
								0.0,0.0,0.0}; //group CO-OH

  size = sizeof(group_tmp_binit3)/sizeof(double);
  assert(size == 45);
	
  for(int i = 0; i < size; ++i)
	BiNIT3.groups[i] = group_tmp_binit3[i];

  // Search the species name in the aerosol species list 
  // and add the species if its name matches with
  // the given list.
  add_species(surrogate, BiNIT3, species_list_aer);

  species BiNIT;
  BiNIT.name="BiNIT";
  BiNIT.is_inorganic_precursor=false;
  BiNIT.Psat_ref=2.5e-6; // Saturation vapor pressure at Tref (torr)
  BiNIT.Tref=298;         // Temperature of reference (K)
  BiNIT.MM=215.0;           // Molar mass (g/mol)
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
  BiNIT.KDiffusion_air=1.0e-5;
  BiNIT.accomodation_coefficient=alpha;
  BiNIT.viscosity=1.68e12;  
  
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
							   0.0,0.0,0.0}; //group CO-OH

  size = sizeof(group_tmp_binit)/sizeof(double);
  assert(size == 45);
	
  for(int i = 0; i < size; ++i)
	BiNIT.groups[i] = group_tmp_binit[i];

  // Search the species name in the aerosol species list 
  // and add the species if its name matches with
  // the given list.
  add_species(surrogate, BiNIT, species_list_aer);
  
  species POAlP;
  POAlP.name="POAlP";
  POAlP.is_inorganic_precursor=false;
  POAlP.MM=280;           // Molar mass (g/mol)
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
  POAlP.KDiffusion_air=1.0e-5;
  POAlP.accomodation_coefficient=alpha;
  POAlP.viscosity=1.68e12;  

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
							   0.0,0.0,0.0}; //group CO-OH

  size = sizeof(group_tmp_poalp)/sizeof(double);
  assert(size == 45);
  for(int i = 0; i < size; ++i)
	POAlP.groups[i] = group_tmp_poalp[i];

  // Search the species name in the aerosol species list 
  // and add the species if its name matches with
  // the given list.
  add_species(surrogate, POAlP, species_list_aer);

  species POAmP;
  POAmP.name="POAmP";
  POAmP.is_inorganic_precursor=false;
  POAmP.MM=280;           // Molar mass (g/mol)
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
  POAmP.KDiffusion_air=1.0e-5;
  POAmP.accomodation_coefficient=alpha;
  POAmP.viscosity=1.68e12;  

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
							   0.0,0.0,0.0}; //group CO-OH

  size = sizeof(group_tmp_poamp)/sizeof(double);
  assert(size == 45);
  for(int i = 0; i < size; ++i)
	POAmP.groups[i] = group_tmp_poamp[i];

  // Search the species name in the aerosol species list 
  // and add the species if its name matches with
  // the given list.
  add_species(surrogate, POAmP, species_list_aer);

  species POAhP;
  POAhP.name="POAhP";
  POAhP.is_inorganic_precursor=false;
  POAhP.MM=280;           // Molar mass (g/mol)
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
  POAhP.KDiffusion_air=1.0e-5;
  POAhP.accomodation_coefficient=alpha;
  POAhP.viscosity=1.68e12;  

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
							   0.0,0.0,0.0}; //group CO-OH

  size = sizeof(group_tmp_poahp)/sizeof(double);
  assert(size == 45);
  for(int i = 0; i < size; ++i)
	POAhP.groups[i] = group_tmp_poahp[i];

  // Search the species name in the aerosol species list 
  // and add the species if its name matches with
  // the given list.
  add_species(surrogate, POAhP, species_list_aer);

  species SOAlP;
  SOAlP.name="SOAlP";
  SOAlP.is_inorganic_precursor=false;
  SOAlP.MM=392;           // Molar mass (g/mol)
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
  SOAlP.KDiffusion_air=1.0e-5;
  SOAlP.accomodation_coefficient=alpha;
  SOAlP.viscosity=1.68e12;  

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
							   0.0,0.0,0.0}; //group CO-OH

  size = sizeof(group_tmp_soalp)/sizeof(double);
  assert(size == 45);
  for(int i = 0; i < size; ++i)
	SOAlP.groups[i] = group_tmp_soalp[i];

  // Search the species name in the aerosol species list 
  // and add the species if its name matches with
  // the given list.
  add_species(surrogate, SOAlP, species_list_aer);

  species SOAmP;
  SOAmP.name="SOAmP";
  SOAmP.is_inorganic_precursor=false;
  SOAmP.MM=392;           // Molar mass (g/mol)
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
  SOAmP.KDiffusion_air=1.0e-5;
  SOAmP.accomodation_coefficient=alpha;
  SOAmP.viscosity=1.68e12;
  

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
							   0.0,0.0,0.0}; //group CO-OH

  size = sizeof(group_tmp_soamp)/sizeof(double);
  assert(size == 45);
  for(int i = 0; i < size; ++i)
	SOAmP.groups[i] = group_tmp_soamp[i];

  // Search the species name in the aerosol species list 
  // and add the species if its name matches with
  // the given list.
  add_species(surrogate, SOAmP, species_list_aer);

  species SOAhP;
  SOAhP.name="SOAhP";
  SOAhP.is_inorganic_precursor=false;
  SOAhP.MM=392;           // Molar mass (g/mol)
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
  SOAhP.KDiffusion_air=1.0e-5;
  SOAhP.accomodation_coefficient=alpha;
  SOAhP.viscosity=1.68e12;  

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
							   0.0,0.0,0.0}; //group CO-OH

  size = sizeof(group_tmp_soahp)/sizeof(double);
  assert(size == 45);
  for(int i = 0; i < size; ++i)
	SOAhP.groups[i] = group_tmp_soahp[i];

  // Search the species name in the aerosol species list 
  // and add the species if its name matches with
  // the given list.
  add_species(surrogate, SOAhP, species_list_aer);

  species Monomer;
  Monomer.name="Monomer";
  Monomer.is_inorganic_precursor=false;
  Monomer.Psat_ref=1.0e-14; // Saturation vapor pressure at Tref (torr)
  Monomer.Tref=298;         // Temperature of reference (K)
  Monomer.MM=278.0;           // Molar mass (g/mol)
  Monomer.deltaH=50.0;     // Enthalpy of vaporization (kJ/mol)
  Monomer.aq_type="none";  // "none","diacid","monoacid" or "aldehyde"
  Monomer.Henry=0.0;       //If the Henry's law constant is set to zero, the model compute the Henry's law constant from the saturation vapour pressure and the activity coefficients at infinite dilution given by unifac
  Monomer.hydrophilic=false;  // Does the species condense on the aqueous phase?
  Monomer.hydrophobic=true;  // Does the species condense on the organic phase?
  Monomer.nonvolatile=false; // Is the compound nonvolatile?
  Monomer.kp_from_experiment=false;  // Use experimental partitioning constant at Tref?
  Monomer.is_organic=true;  // Is the compound organic?
  Monomer.compute_gamma_org=false;  // Compute the activity coefficients of the organic phase for this compound?
  Monomer.compute_gamma_aq=false;  // Compute the activity coefficients of the aqueous phase for this compound?
  Monomer.Koligo_org=0.0;      //oligomeriation constant in the organic phase
  Monomer.rho=1300.0;
  Monomer.KDiffusion_air=1.0e-5;
  Monomer.accomodation_coefficient=alpha;
  Monomer.viscosity=1.68e12;
  
  
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
  							   0.0,0.0,0.0}; //group CO-OH

  size = sizeof(group_tmp_Monomer)/sizeof(double);
  assert(size == 45);
	
  for(int i = 0; i < size; ++i)
  	Monomer.groups[i] = group_tmp_Monomer[i];

  // Search the species name in the aerosol species list 
  // and add the species if its name matches with
  // the given list.
  add_species(surrogate, Monomer, species_list_aer);

  species Dimer;
  Dimer.name="Dimer";
  Dimer.is_inorganic_precursor=false;
  Dimer.Psat_ref=1.0e-14; // Saturation vapor pressure at Tref (torr)
  Dimer.Tref=298;         // Temperature of reference (K)
  Dimer.MM=432.0;           // Molar mass (g/mol)
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
  Dimer.KDiffusion_air=1.0e-5;
  Dimer.accomodation_coefficient=alpha;
  Dimer.viscosity=1.68e12;
  
  
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
  							   0.0,0.0,0.0}; //group CO-OH

  size = sizeof(group_tmp_dimer)/sizeof(double);
  assert(size == 45);
	
  for(int i = 0; i < size; ++i)
  	Dimer.groups[i] = group_tmp_dimer[i];

  // Search the species name in the aerosol species list 
  // and add the species if its name matches with
  // the given list.
  add_species(surrogate, Dimer, species_list_aer);


  species BiA3D;
  BiA3D.name="BiA3D";
  BiA3D.is_inorganic_precursor=false;
  BiA3D.Psat_ref=3.25e-7; // Saturation vapor pressure at Tref (torr)
  BiA3D.kp_from_experiment=false;  // Use experimental partitioning constant at Tref?
  BiA3D.Tref=298;         // Temperature of reference (K)
  BiA3D.MM=204;           // Molar mass (g/mol)
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
  BiA3D.accomodation_coefficient=0.5;
  BiA3D.KDiffusion_air=1.0e-5;
  BiA3D.viscosity=1.68e12;

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
                               0.0,0.0,0.0}; //group CO-OH
  


  size = sizeof(group_tmp_bia3d)/sizeof(double);
  assert(size = 45);
	
  for(int i = 0; i < size; ++i)
	BiA3D.groups[i] = group_tmp_bia3d[i];

  // Search the species name in the aerosol species list 
  // and add the species if its name matches with
  // the given list.
  add_species(surrogate, BiA3D, species_list_aer);

  species orgNIT ;
  orgNIT.name="orgNIT";
  orgNIT.is_inorganic_precursor=false;
  orgNIT.Psat_ref=5.0e-6; // Saturation vapor pressure at Tref (torr)
  orgNIT.Tref=298;         // Temperature of reference (K)
  orgNIT.MM=278.0;           // Molar mass (g/mol)
  orgNIT.deltaH=40.0;     // Enthalpy of vaporization (kJ/mol)
  orgNIT.aq_type="none";  // "none","diacid","monoacid" or "aldehyde"
  orgNIT.Henry=0.0;       //If the Henry's law constant is set to zero, the model compute the Henry's law constant from the saturation vapour pressure and the activity coefficients at infinite dilution given by unifac
  orgNIT.hydrophilic=false;  // Does the species condense on the aqueous phase?
  orgNIT.hydrophobic=true;  // Does the species condense on the organic phase?
  orgNIT.nonvolatile=false; // Is the compound nonvolatile?
  orgNIT.kp_from_experiment=false;  // Use experimental partitioning constant at Tref?
  orgNIT.is_organic=true;  // Is the compound organic?
  orgNIT.compute_gamma_org=false;  // Compute the activity coefficients of the organic phase for this compound?
  orgNIT.compute_gamma_aq=false;  // Compute the activity coefficients of the aqueous phase for this compound?
  orgNIT.Koligo_org=0.0;      //oligomeriation constant in the organic phase
  orgNIT.rho=1300.0;
  orgNIT.KDiffusion_air=1.0e-5;
  orgNIT.accomodation_coefficient=alpha;
  orgNIT.viscosity=1.68e12;
  
  
  //Group: if no functionnal group in the species use the default species
  //for the computation of activity coefficients
  //
  double group_tmp_orgNIT [] = {0.0,0.0,0.0,0.0, // group C
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
  							   0.0,0.0,0.0}; //group CO-OH

  size = sizeof(group_tmp_orgNIT)/sizeof(double);
  assert(size == 45);
	
  for(int i = 0; i < size; ++i)
  	orgNIT.groups[i] = group_tmp_orgNIT[i];

  // Search the species name in the aerosol species list 
  // and add the species if its name matches with
  // the given list.
  add_species(surrogate, orgNIT, species_list_aer);


  species H2O;
  H2O.name="H2O";
  H2O.is_inorganic_precursor=false;
  H2O.MM=18.0;           // Molar mass (g/mol)
  H2O.is_organic=false;  // Is the compound organic?
  H2O.hydrophilic=true; // Does the species condense on the aqueous phase?
  H2O.hydrophobic=true;  // Does the species condense on the organic phase?
  H2O.nonvolatile=false;
  H2O.kp_from_experiment=false; 
  H2O.compute_gamma_org=true;  // Compute the activity coefficients of the organic phase for this compound?
  H2O.compute_gamma_aq=true;  // Compute the activity coefficients of the aqueous phase for this compound
  H2O.Koligo_org=0.0;
  H2O.rho=1000.0;
  H2O.KDiffusion_air=1.0e-5;
  H2O.accomodation_coefficient=alpha;
  H2O.viscosity=1.0;  

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
							 0.0,0.0,0.0}; //group CO-OH

  size = sizeof(group_tmp_h2o)/sizeof(double);
  assert(size == 45);
  
  for(int i = 0; i < size; ++i)
    H2O.groups[i] = group_tmp_h2o[i];

  // Search the species name in the aerosol species list 
  // and add the species if its name matches with
  // the given list.
  add_species(surrogate, H2O, species_list_aer);

  species SO4;
  SO4.name="SO4";
  SO4.is_inorganic_precursor=false;
  SO4.MM=96.0;           // Molar mass (g/mol)
  SO4.is_organic=false;  // Is the compound organic?
  SO4.hydrophilic=true; // Does the species condense on the aqueous phase?
  SO4.hydrophobic=false;  // Does the species condense on the organic phase?
  SO4.compute_gamma_org=false;  // Compute the activity coefficients of the organic phase for this compound?
  SO4.compute_gamma_aq=false;  // Compute the activity coefficients of the aqueous phase for this compound
  SO4.charge=-2;
  SO4.index_ion_aiomfac=11;
  SO4.rho=1300.0; //1400.0;    
  SO4.KDiffusion_air=1.0e-5;
  SO4.accomodation_coefficient=0.5;
  SO4.viscosity=1.0;
  SO4.soap_ind = -1;
  surrogate.push_back(SO4);

  species HSO4;
  HSO4.name="HSO4";
  HSO4.is_inorganic_precursor=false;
  HSO4.MM=97.0;           // Molar mass (g/mol)
  HSO4.is_organic=false;  // Is the compound organic?
  HSO4.hydrophilic=true; // Does the species condense on the aqueous phase?
  HSO4.hydrophobic=false;  // Does the species condense on the organic phase?
  HSO4.compute_gamma_org=false;  // Compute the activity coefficients of the organic phase for this compound?
  HSO4.compute_gamma_aq=false;  // Compute the activity coefficients of the aqueous phase for this compound
  HSO4.charge=-1;
  HSO4.index_ion_aiomfac=10;
  HSO4.rho=1300.0; //1400.0;
  HSO4.KDiffusion_air=1.0e-5;
  HSO4.accomodation_coefficient=0.5;
  HSO4.viscosity=1.0;
  HSO4.soap_ind = -1;
  surrogate.push_back(HSO4);

  species NO3;
  NO3.name="NO3";
  NO3.is_inorganic_precursor=false;
  NO3.MM=62.0;           // Molar mass (g/mol)
  NO3.is_organic=false;  // Is the compound organic?
  NO3.hydrophilic=true; // Does the species condense on the aqueous phase?
  NO3.hydrophobic=false;  // Does the species condense on the organic phase?
  NO3.compute_gamma_org=false;  // Compute the activity coefficients of the organic phase for this compound?
  NO3.compute_gamma_aq=false;  // Compute the activity coefficients of the aqueous phase for this compound
  NO3.charge=-1;
  NO3.index_ion_aiomfac=9;
  NO3.rho=1300.0; //1400.0;
  NO3.KDiffusion_air=1.0e-5;
  NO3.accomodation_coefficient=0.5;
  NO3.viscosity=1.0;
  NO3.soap_ind = -1;
  surrogate.push_back(NO3);

  species NH4;
  NH4.name="NH4";
  NH4.is_inorganic_precursor=false;
  NH4.MM=18.0;           // Molar mass (g/mol)
  NH4.is_organic=false;  // Is the compound organic?
  NH4.hydrophilic=true; // Does the species condense on the aqueous phase?
  NH4.hydrophobic=false;  // Does the species condense on the organic phase?
  NH4.compute_gamma_org=false;  // Compute the activity coefficients of the organic phase for this compound?
  NH4.compute_gamma_aq=false;  // Compute the activity coefficients of the aqueous phase for this compound
  NH4.charge=1.0;
  NH4.index_ion_aiomfac=4;
  NH4.rho=1300.0; //1400.0;
  NH4.KDiffusion_air=1.0e-5;
  NH4.accomodation_coefficient=0.5;
  NH4.viscosity=1.0;
  NH4.soap_ind = -1;
  surrogate.push_back(NH4);

  species H;
  H.name="H";
  H.is_inorganic_precursor=false;
  H.MM=1.0;           // Molar mass (g/mol)
  H.is_organic=false;  // Is the compound organic?
  H.hydrophilic=true; // Does the species condense on the aqueous phase?
  H.hydrophobic=false;  // Does the species condense on the organic phase?
  H.compute_gamma_org=false;  // Compute the activity coefficients of the organic phase for this compound?
  H.compute_gamma_aq=false;  // Compute the activity coefficients of the aqueous phase for this compound
  H.charge=1.0;
  H.rho=1300.0; //1400.0;
  H.index_ion_aiomfac=0;
  H.KDiffusion_air=1.0e-5;
  H.accomodation_coefficient=0.5;
  H.viscosity=1.0;
  H.soap_ind = -1;
  surrogate.push_back(H);

  species Na;
  Na.name="Na";
  Na.is_inorganic_precursor=false;
  Na.MM=23.0;           // Molar mass (g/mol)
  Na.is_organic=false;  // Is the compound organic?
  Na.hydrophilic=true; // Does the species condense on the aqueous phase?
  Na.hydrophobic=false;  // Does the species condense on the organic phase?
  Na.compute_gamma_org=false;  // Compute the activity coefficients of the organic phase for this compound?
  Na.compute_gamma_aq=false;  // Compute the activity coefficients of the aqueous phase for this compound
  Na.charge=1.0;
  Na.index_ion_aiomfac=2;
  Na.rho=1300.0; //1400.0;
  Na.KDiffusion_air=1.0e-5;
  Na.accomodation_coefficient=0.5;
  Na.viscosity=1.0;

  // Find the number in the aerosol species list
  Na.soap_ind = -1;
  for (int i = 0; i < nsp; ++i)
    if (species_list_aer[i].substr(1,-1) == Na.name)
      Na.soap_ind = i;
  surrogate.push_back(Na);

  species Cl;
  Cl.name="Cl";
  Cl.is_inorganic_precursor=false;
  Cl.MM=35.5;           // Molar mass (g/mol)
  Cl.is_organic=false;  // Is the compound organic?
  Cl.hydrophilic=true; // Does the species condense on the aqueous phase?
  Cl.hydrophobic=false;  // Does the species condense on the organic phase?
  Cl.compute_gamma_org=false;  // Compute the activity coefficients of the organic phase for this compound?
  Cl.compute_gamma_aq=false;  // Compute the activity coefficients of the aqueous phase for this compound
  Cl.charge=-1.0;
  Cl.index_ion_aiomfac=7;
  Cl.rho=1300.0; // 1400.0;
  Cl.KDiffusion_air=1.0e-5;
  Cl.accomodation_coefficient=0.5;
  Cl.viscosity=1.0;
  Cl.soap_ind = -1;

  species H2SO4;
  H2SO4.name="H2SO4";
  H2SO4.MM=98.0;
  H2SO4.nonvolatile=true;
  H2SO4.is_organic=false;
  H2SO4.is_inorganic_precursor=true;
  H2SO4.Henry=1.0e13;     // Enthalpy of vaporization (kJ/mol)
  H2SO4.deltaH=0.0;     // Henry's law constant at Tref (M/atm)
  H2SO4.hydrophilic=true; // Does the species condense on the aqueous phase?
  H2SO4.hydrophobic=false;  // Does the species condense on the organic phase?
  H2SO4.compute_gamma_org=false;  // Compute the activity coefficients of the organic phase for this compound?
  H2SO4.compute_gamma_aq=false;  // Compute the activity coefficients of the aqueous phase for this compound
  H2SO4.accomodation_coefficient=0.5;  
  H2SO4.KDiffusion_air=1.07e-5;
  H2SO4.viscosity=1.0;

  // Find the number in the aerosol species list
  H2SO4.soap_ind = -1;
  for (int i = 0; i < nsp; ++i)
    if (species_list_aer[i].substr(1,-1) == "SO4")
      H2SO4.soap_ind = i;
  surrogate.push_back(H2SO4);

  species NH3;
  NH3.name="NH3";
  NH3.MM=17.0;
  NH3.is_organic=false;
  NH3.nonvolatile=false;
  NH3.is_inorganic_precursor=true;
  NH3.Henry=57.64;     // Enthalpy of vaporization (kJ/mol)
  NH3.deltaH=0.0;     // Henry's law constant at Tref (M/atm)
  NH3.hydrophilic=true; // Does the species condense on the aqueous phase?
  NH3.hydrophobic=false;  // Does the species condense on the organic phase?
  NH3.compute_gamma_org=false;  // Compute the activity coefficients of the organic phase for this compound?
  NH3.compute_gamma_aq=false;  // Compute the activity coefficients of the aqueous phase for this compound
  NH3.accomodation_coefficient=0.5;
  NH3.KDiffusion_air=2.17e-5;
  NH3.viscosity=1.0;

  // Find the number in the aerosol species list
  NH3.soap_ind = -1;
  for (int i = 0; i < nsp; ++i)
    if (species_list_aer[i].substr(1,-1) == "NH4")
      NH3.soap_ind = i;
  surrogate.push_back(NH3);

  species HNO3;
  HNO3.name="HNO3";
  HNO3.MM=63.0;
  HNO3.nonvolatile=false;
  HNO3.is_organic=false;
  HNO3.is_inorganic_precursor=true;
  HNO3.Henry=2.1e5;     // Enthalpy of vaporization (kJ/mol)
  HNO3.deltaH=0.0;     // Henry's law constant at Tref (M/atm)
  HNO3.hydrophilic=true; // Does the species condense on the aqueous phase?
  HNO3.hydrophobic=false;  // Does the species condense on the organic phase?
  HNO3.compute_gamma_org=false;  // Compute the activity coefficients of the organic phase for this compound?
  HNO3.compute_gamma_aq=false;  // Compute the activity coefficients of the aqueous phase for this compound
  HNO3.accomodation_coefficient=0.5;
  HNO3.KDiffusion_air=1.47e-5;
  HNO3.viscosity=1.0;

  // Find the number in the aerosol species list
  HNO3.soap_ind = -1;
  for (int i = 0; i < nsp; ++i)
    if (species_list_aer[i].substr(1,-1) == "NO3")
      HNO3.soap_ind = i;
  surrogate.push_back(HNO3);

  species HCl;
  HCl.name="HCl";
  HCl.MM=36.5;
  HCl.nonvolatile=false;
  HCl.is_organic=false;
  HCl.is_inorganic_precursor=true;
  HCl.Henry=2.5e3;     // Enthalpy of vaporization (kJ/mol)
  HCl.deltaH=0.0;     // Henry's law constant at Tref (M/atm)
  HCl.hydrophilic=true; // Does the species condense on the aqueous phase?
  HCl.hydrophobic=false;  // Does the species condense on the organic phase?
  HCl.compute_gamma_org=false;  // Compute the activity coefficients of the organic phase for this compound?
  HCl.compute_gamma_aq=false;  // Compute the activity coefficients of the aqueous phase for this compound
  HCl.accomodation_coefficient=0.5;
  HCl.KDiffusion_air=1.72e-5;
  HCl.viscosity=1.0;

  // Find the number in the aerosol species list
  HCl.soap_ind = -1;
  for (int i = 0; i < nsp; ++i)
    if (species_list_aer[i].substr(1,-1) == "HCL")
      HCl.soap_ind = i;

  surrogate.push_back(HCl);
}

#endif
