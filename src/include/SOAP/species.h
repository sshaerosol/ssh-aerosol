//!!-----------------------------------------------------------------------
//!!     Copyright (C) 2019 CEREA (ENPC) - INERIS
//!!     SSH-aerosol is distributed under the GNU General Public License v3
//!!-----------------------------------------------------------------------

#ifndef SPECIES_H
#define SPECIES_H

#include <blitz/array.h>
using namespace blitz;
using namespace std;

namespace ssh_soap
{  
  class model_config
  {
  public:
    bool equilibrium;        //True: system at equilibrium, False: Dynamic approach
	
	//common parameters:
	bool coupled_phases; //Is the system coupled (at least one compound condense on both phase)?
	bool hygroscopicity; //Does hygroscopicity has to be computed?
	bool compute_long_and_medium_range_interactions; // Take into account the influence of ions on
	                                                //    activity coefficients?
	bool compute_saturation; //Does phase separation has to be computed?
	bool initialized_saturation;
	bool compute_inorganic,compute_organic;
        bool coupling_organic_inorganic;
        bool compute_rho_aqueous;
        bool compute_aqueous_phase_properties;
	bool explicit_representation;	
        int number_of_org_inorg_cycles;
	int max_number_of_phases; //maximal number of organic phases
	string activity_model; // unifac or ideal
	double LWClimit; //LWC under which all the compounds there is no aqueous phase 
	double RHcoupling;
	int max_iter; //maximal number of iterations for the newton raphson method
	int iH2O,iHp; //indexes for water and H+
	int iHSO4m,iSO4mm;
	int iNH4p;
	int iNO3m;
	int iClm;        
        int iNa,iCa,iMg,iK;
        int iHCl,iHNO3,iNH3,iH2SO4;

	double rho_organic,rho_aqueous; //volumic masses of the organic phase and the aqueous phase
	double Ke, moligo;

	//equilibrium parameters:
	double precision; //absolute precision under which the system has been solved 
	
	//dynamic parameters:
	double relative_precision; //relative precision under which the system has been solved 
	bool use_global_dynamic_parameters;  //Assume the same composition over all bins and layers
	bool compute_kelvin_effect;      //Should the kelvin effect be computed?
	bool first_evaluation_of_saturation;  //Use initial concentrations to compute phase separation?
	bool first_evaluation_activity_coefficients; //Use initial concentrations to compute activity
                                                 // coefficients?
	int nlayer;  //number of layer into which divided the particle is divided   
	int nbins;   //number of bins for the diameter of particles
	Array<int, 2> nphase;  //number of phases in the particle
	double MOmin; //minimal concentrations of organics
	double tequilibrium;    //time under which equilibrium is assumed
	double EPSER;           //relative difference of ros2
	double deltatmin;       //minimal time step
	double surface_tension_org,surface_tension_aq; //surface tension of the organic
	                                               //and aqueous phases
	Array<double, 1> Vlayer; //Volume fraction of a layer (assumed equal to the mass fraction)
	Array<double, 1> alpha_layer; //ratio of the characteristic of the layer to the characteristic
	                              //time of diffusion
	Array<double, 2> Alayer;      //polynomial coefficients for the effect of the solid fraction fs
	                              //on the diffusion
	Array<double, 1> diameters;   //diameters of particles for each bin
	Array<double, 1> AQrho;       //density of the aqueous phase for each bin
        double kp_low_volatility;
	
	//Parameters of unifac and aiomfac:
	double Z;
	int nion_aiomfac, nion_unifac, ngroup_aiomfac, nmol_aiomfac;
	int nfunc_org, nfunc_aq, nmol_org, nmol_aq,nfunc_tot,nmol_tot;	
	Array<double, 2> Inter_org, Inter_aq, Inter2_aq, Inter2_org, Inter2_tot, Inter_tot, groups_org, groups_aq, groups_tot;

	Array<double, 2> InterB_org, InterB_aq, InterB_tot;
	Array<double, 2> InterC_org, InterC_aq, InterC_tot;
	Array<double, 1> RG_org, QG_org, RG_aq, QG_aq, RG_tot, QG_tot, RGions, QGions,Lions;
	Array<double, 2> groups_aiomfac,Mgroups_aiomfac,gMg_aiomfac;
	Array<double, 2> b1ki_aq,b2ki_aq,b1ca_aq,b2ca_aq,b3ca_aq,c1ca_aq,c2ca_aq,Rcc_aq;
	Array<double, 3> Qcca_aq;
	Array<double,2> dbound,Radius;
        Array<double,1> Rparam_org,Qparam_org,Lparam_org,Rparam_aq,Qparam_aq,Lparam_aq,Rparam_tot,Qparam_tot,Lparam_tot;
	Array<double,2> surface_fraction_moltot,surface_fraction_molorg,surface_fraction_molaq;
	Array<double,2> sum2mol_tot,sum2mol_org,sum2mol_aq;
	Array<double,2> group_activity_moltot,group_activity_molorg,group_activity_molaq;
	double dorg;
	int nh_max,nh_init,nh_aq_init,nh_org_init,nh_inorg_init;	
	int explicit_method;
	bool constant_dorg;	
        int ntemp; 
        bool SR_ions;
	bool temperature_dependancy;
        double koligo; //s-1
        double Keq_oligo;
        int nt;
        double dtchem_min;
	bool chemistry,solids;

        //For AIOMFAC:
        Array<double,1> molality,gamma_LR_ions,gamma_MR_ions, charges_ions,molar_mass_groups;
        Array<double,1> gamma_LR_solvents,gamma_MR_solvents,X_aiomfac,molar_mass_solvents;
	Array<double,1> gamma_ions_inf;

        int iiter;
	int imethod;

    double chpinit,ionicinit,initAQ;
    double molalmax;
    bool to_be_rejected;

  };

  class species
  {
  public:
	//Global parameters of species (see species.cxx)
	double Psat_ref,MM,deltaH,Tref,Koligo_org,Henry,GAMMAinf,kp_experiment;
	double Kacidity1,Kacidity2;
	double Koligo_aq,beta,pHref;
	double rho,viscosity;
	string aq_type,name;
	bool hydrophilic, hydrophobic,nonvolatile,kp_from_experiment,is_organic;
	bool compute_gamma_org,compute_gamma_aq,is_inorganic_precursor;
	int index_gamma_org, index_gamma_aq, index_gamma_tot,index_gamma_aiomfac;
	double groups[56],conc_mol;
	int index_ion_aiomfac,index_ion;
	int aqt;
	
	//local equilibrium values
	double Kaq_inorg;
	double Psat_ssh(double&);  
	double Kpart_org_ssh(double&, double&);
	double Kpart_aq_ssh(double&, double&);
	double Kp_eff_org_ssh(double&, double&);
	double Kp_exp_org_ssh(double&);
	double knudsen_function_ssh(double&, double, double);
	double Kp_eff_aq_ssh(model_config&, double&, double&, double&, double&, double&, double&, double&, double&);
        double Kp_eff_aqreal_ssh(model_config&, double&, double&, double&, double&, double&, double&, double&, double&);
        double Kp_eff_aqrealdyn_ssh(model_config&, double&, double&, double&, double&, double&, double&, double&, double&, int&);
	double Kequilibrium_ssh(double &);
	//Ap: concentrations in the organic phase
	//Ag: concentrations in the gas phase
	//Aaq: concentrations in the aqueous phase
	//Atot: total concentrations (Ap+Ag+Aaq)
	//Xaq: molar fraction in the aqueous phase
        //Waq: mass fraction in the aqueous phase
	//Xorg: molar fraction in the organic phase
	double kpi,keq,kaqi,fioni1,fioni2;
	double Ap,Ag,Aaq,Atot,Atot0,Atot1,Xaq,Xorg,Waq;
	double Aaq_save;
	double gamma_aq_old,Aaq_old,Xaq_old,Ag_old,Ap_old;
	double partial_pressure, partial_pressure_old;
	double gamma_org,gamma_org_old,gamma_aq,gamma_LR,gamma_SRMR,gamma_aq_tmp;
	double molality,charge;
	Array<double, 1> gamma_org_sat,Ap_sat,Ap_sat_save;  //for saturation 
	Array<double, 1> gamma_org_sat_old,Ap_sat_old,Xorg_sat,Xorg_sat_old; //for saturation
	int jmain_phase,jmain_phase_old; //for saturation

        //For solid species
        double Ksol,dCp;
        int iion1,iion2,iion3,pion1,pion2,pion3;
        string ion1,ion2,ion3;
        bool is_solid;
        double Ap2,Aaq2,Ag2,molality2;

	//local dynamic value:
	Array<double, 3> Ap_layer,Ap_layer_init,Ap_layer_init0,gamma_org_layer,gamma_org_layer0,Xinit;
	Array<double, 3> Kp;
	Array<double, 1> Aaq_bins,Aaq_bins_init,Aaq_bins_init0,time_aq,LR,SRMR,Kaq,dKaq;
	Array<double, 1> gamma_aq_bins,gamma_old;
	double KDiffusion_p,KDiffusion_air,accomodation_coefficient,Ag0,Aaq0;
	Array<double, 3> tau_diffusion,time;
	Array<double, 1> tau_air;
	Array <double, 4> k1,Jdn;
	Array<double, 2> k1_aq,Jdn_aq;
	Array<double, 2> dif_org;
        Array<double, 1> species_activity;
        Array<double, 1> Jdn_gas,flux_chem_tot;
        Array<double, 4> flux_chem;
	Array<double, 2> flux_chem_aq;
        Array<double, 1> flux_chem_gas;
	Array<double, 1> veckaqi,vecfioni1,vecfioni2;
	Array<double, 3> kprod,kloss,kloc;
	Array<double, 1> kprod_aq,kloss_aq,k1_gas;
	double kprod_gas,kloss_gas;
	double Ag1,Agt,fion1,fion2,ktot1,ktot2,Jdn_tot;	
	double moligo;
        bool is_monomer;
        bool is_ion,is_solvent;
        string name_oligomer;
        int ioligo;
        bool rion;           
        int nion;
        Array<string, 1> ion,rion_product;
        Array<bool, 1> rion_catalyzed;
        Array<double, 1> kion;
        Array<int, 1> iion,iproduct;
	double velocity,knui;

	double Aginit,Aaqinit,Apinit;
        double deltat_exp;
       
        int soap_ind; // Number in the aerosol species list
        int soap_ind_aero; // Number in the aerosol species list including layers
        string smile;
        bool is_generic;
	Array <double, 1> fac_corr_ph;
  };
  
}

#endif
