#ifndef SPECIES_H
#define SPECIES_H

#include <blitz/array.h>
using namespace blitz;
using namespace std;

namespace soap
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
        int iNa;
	double rho_organic,rho_aqueous; //volumic masses of the organic phase and the aqueous phase

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
	Array<double, 2> Inter_org, Inter_aq, Inter_tot, groups_org, groups_aq, groups_tot;
	Array<double, 2> InterB_org, InterB_aq, InterB_tot;
	Array<double, 2> InterC_org, InterC_aq, InterC_tot;
	Array<double, 1> RG_org, QG_org, RG_aq, QG_aq, RG_tot, QG_tot, RGions, QGions,Lions;
	Array<double, 2> groups_aiomfac,Mgroups_aiomfac;
	Array<double, 2> b1ki_aq,b2ki_aq,b1ca_aq,b2ca_aq,b3ca_aq,c1ca_aq,c2ca_aq,Rcc_aq;
	Array<double, 3> Qcca_aq;
	Array<double,2> dbound,Radius;
        Array<double,1> Rparam_org,Qparam_org,Lparam_org,Rparam_aq,Qparam_aq,Lparam_aq,Rparam_tot,Qparam_tot,Lparam_tot;
	double dorg;
	int nh_max,nh_init,nh_aq_init,nh_org_init,nh_inorg_init;	
	int explicit_method;
	bool constant_dorg;	
        int ntemp; 
        bool tabulation_unifac;
        bool SR_ions;
	bool temperature_dependancy;
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
	double groups[45];
	int index_ion_aiomfac,index_ion;
	
	//local equilibrium values
	double Kaq_inorg;
	double Psat(double&);  
	double Kpart_org(double&, double&);
	double Kpart_aq(double&, double&);
	double Kp_eff_org(double&, double&);
	double Kp_exp_org(double&);
	double knudsen_function(double&, double);
	double Kp_eff_aq(model_config&, double&, double&, double&, double&, double&, double&, double&, double&);
        double Kp_eff_aqreal(model_config&, double&, double&, double&, double&, double&, double&, double&, double&);
	double Kequilibrium(double &);
	//Ap: concentrations in the organic phase
	//Ag: concentrations in the gas phase
	//Aaq: concentrations in the aqueous phase
	//Atot: total concentrations (Ap+Ag+Aaq)
	//Xaq: molar fraction in the aqueous phase
        //Waq: mass fraction in the aqueous phase
	//Xorg: molar fraction in the organic phase
	double kpi,keq,kaqi,fioni1,fioni2;
	double Ap,Ag,Aaq,Atot,Xaq,Xorg,Waq;
	double Aaq_save;
	double gamma_aq_old,Aaq_old,Xaq_old,Ag_old,Ap_old;
	double partial_pressure, partial_pressure_old;
	double gamma_org,gamma_aq,gamma_LR,gamma_SRMR;
	double molality,charge;
	Array<double, 1> gamma_org_sat,Ap_sat,Ap_sat_save;  //for saturation 
	Array<double, 1> gamma_org_sat_old,Ap_sat_old,Xorg_sat,Xorg_sat_old; //for saturation
	int jmain_phase,jmain_phase_old; //for saturation

	//local dynamic value:
	Array<double, 3> Ap_layer,Ap_layer_init,Ap_layer_init0,gamma_org_layer,gamma_org_layer0,Xinit;
	Array<double, 3> Kp;
	Array<double, 1> Aaq_bins,Aaq_bins_init,Aaq_bins_init0,time_aq,LR,SRMR,Kaq;
	Array<double, 1> gamma_aq_bins,gamma_old;
	double KDiffusion_p,KDiffusion_air,accomodation_coefficient,Ag0;
	Array<double, 3> tau_diffusion,time;
	Array<double, 1> tau_air;
	Array <double, 4> k1,Jdn;
	Array<double, 2> k1_aq,Jdn_aq;
	Array<double, 2> dif_org;
        Array<double, 1> species_activity;
	double Ag1,Agt,fion1,fion2;

        int soap_ind; // Number in the aerosol species list
  };
  
}

#endif
