//////////////
// INCLUDES //

#include "soap.h"
#include "parameters.cxx"
#include <iostream>
#include <cmath>
#include <blitz/array.h>
#include "omp.h"

using namespace blitz;

// INCLUDES //
//////////////

extern "C" void soap_main_(double* LWC, double* RH, double* Temperature, 
                           double* ionic, double* chp, double& LWCorg,
                           double& deltat, double* DSD, double* csol, double* liquid,
                           int* ns_aer, int* neq, double* q, double* qaero, 
                           double* qgas, double* lwc_Nsize,
                           double* ionic_Nsize, double* chp_Nsize,
                           double* liquid_Nsize, int* nbin){

  return soap_main(*LWC, *RH, *Temperature, 
                   *ionic, *chp, LWCorg, 
                   deltat,DSD,csol, liquid,
                   *ns_aer, *neq, q, qaero, qgas,
                   lwc_Nsize,ionic_Nsize,chp_Nsize,
                   liquid_Nsize, *nbin);

}

/*! \brief Main function of SOAP
  \param LWC liquid water content absorbed by inorganic particles in the aqueous phase ([\mu g.m^-3]).
  \param LWCorg liquid water content absorbed by organic particles in the aqueous phase ([\mu g.m^-3]).
  \param psoap_config pointer of SOAP configuration 
  \param ionic ionic strength
  \param chp concentrations of H+ ions
  \param liquid liquid inorganic aerosols
*/
void soap_main(double LWC, double RH, double Temperature,
               double ionic, double chp, double& LWCorg,
               double& deltat, double DSD[], double csol[], double liquid[],
               int ns_aer, int neq, double q[], double qaero[], double qgas[],
               double lwc_Nsize[], double ionic_Nsize[], double chp_Nsize[],
               double liquid_Nsize[], int nbin)
{

  /*** General parameters and options ***/
  // before calling SOAP: if organics are dynamic, then the arrays
  // lwc_Nsize[], ionic_Nsize[], chp_Nsize[], liquid_Nsize need to be specified
  // if oganics are at equilibrium, then lwc, ionic, chp and liquid need to be specified.

  int i, b;

  // Make a reference with 'model_config' type to the pointer of SOAP configuration.
  model_config config;

  // Make a reference with 'surrogate' type to the pointer of SOAP species.
  vector<species> surrogate;


  // The order of species should not be changed (YK).
  vector<string> species_list_aer;
  string tmp[34] = {"PMD","PBC","PNA","PSO4",
		    "PNH4","PNO3","PHCL","PBiA2D","PBiA1D","PBiA0D","PAGLY","PAMGLY",
		    "PBiMT","PBiPER","PBiDER","PBiMGA","AnBlP","PAnBmP","PBiBlP","PBiBmP",
		    "PBiNGA","PBiNIT3","PBiNIT","PAnClP","PSOAlP","PSOAmP","PSOAhP",
                    "PPOAlP","PPOAmP","PPOAhP","PMonomer", "PDimer", "PBiA3D", "PH2O"}; 
  for (i = 0; i < 34; i++)
    species_list_aer.push_back(tmp[i]);  // YK
  cout << "Call soap:" << species_list_aer[0] << endl;
  config.nbins = nbin;
  config.equilibrium = true;
  config.nlayer=1; /// A MODIFIER

  parameters(config, surrogate, species_list_aer);

  // If Na and Cl are included.
  bool NaCl = false; 


  /*** Initialiation ***/

  // Initial mean molar mass in the organic phase.
  double MOW=200.0;

  double MOinit = 0.0;
  double AQinit = 0.0;
  int n = surrogate.size();
  
  // Compute the activity coefficients at infinite dilution 
  // and the Henry's law constant 
  compute_gamma_infini(config, surrogate);

  // Initialization for the equilibrium approach
  for (i = 0; i < n; ++i)
    {
      surrogate[i].Ap = 0.0;
      surrogate[i].Aaq = 0.0;
      surrogate[i].Ag = 0.0;
    }

  /*** Use the global equilibrium approach ***/
  
  if (config.equilibrium)
    {     
      for (i = 0; i < n; ++i)
	{
	  int iq = surrogate[i].soap_ind;
	  if (iq != -1)
	    {
	      if (surrogate[i].is_organic)
		if (surrogate[i].hydrophilic and LWC > config.LWClimit)
		  {
		    surrogate[i].Ag = qgas[iq];
		    surrogate[i].Ap = 0.0;
		    surrogate[i].Aaq = qaero[iq];
		  }
	      // Organic aersol is hydrophobic or LWC is equal or less than LWClimit
		else 
		  {
		    surrogate[i].Ag = qgas[iq];
		    surrogate[i].Ap = qaero[iq];
		    surrogate[i].Aaq = 0.0;
		  }

	      // Inorganic gas concentrations
	      if (surrogate[i].is_inorganic_precursor)
		{
		  surrogate[i].Ag = qgas[iq];
		  surrogate[i].Ap = 0.0;
		  surrogate[i].Aaq = 0.0;
		}
	    }
	  surrogate[i].Atot = surrogate[i].Ap + surrogate[i].Ag + surrogate[i].Aaq;
	}

      // Dissociated inorganic aerosols.
      // Assume that NH3 and HNO3 are totally dissociated.
      surrogate[config.iSO4mm].Aaq = liquid[4]*1.0e6*96.0; // SO4_2-
      surrogate[config.iHSO4m].Aaq = liquid[5]*1.0e6*97.0; // HSO4-
      surrogate[config.iNO3m].Aaq = liquid[6]*1.0e6*62.0; // NO3-
      surrogate[config.iNH4p].Aaq = liquid[2]*1.0e6*18.0; // NH4+
      surrogate[config.iHp].Aaq = liquid[0]*1.0e6; // H+

      if (NaCl)
	for (i = 0; i < n; ++i) // Na+ and Cl-
	  if ((surrogate[i].name == "Na") or (surrogate[i].name == "Cl"))
	    surrogate[i].Aaq = qaero[surrogate[i].soap_ind];


      // Set the total values in each phase. 
      if(LWC<config.LWClimit)
	{
	  LWC=0.0;
	  for(i=0;i<n;++i)
	    surrogate[i].Aaq=0.0;
	}
      else
	{
	  AQinit = LWC;
	  for (i = 0; i < n; ++i)      
	    {
	      AQinit += surrogate[i].Aaq;
	      MOinit += surrogate[i].Ap;
	    }
	}

      global_equilibrium(config,surrogate,
                         MOinit,MOW,
                         LWC, AQinit, ionic, chp,
                         Temperature, RH);

      // Give back the concentrations
      for (i = 0; i < n; ++i)
        {
          int iq = surrogate[i].soap_ind;
          if (iq != -1)
            {
              if (surrogate[i].is_organic)
                {
                  qgas[iq] = surrogate[i].Ag;
                  qaero[iq] = surrogate[i].Ap + surrogate[i].Aaq; 
                }
              if (i == config.iH2O)
                {
                  qaero[iq] = LWC + surrogate[i].Ap + surrogate[i].Aaq;
                  LWCorg = surrogate[i].Ap + surrogate[i].Aaq;
                }
            }
        }

    }

  /*** Use the dynamic approach ***/

  else
    {
      // Initialization for the dynamic approach	  
      int b,ilayer,iphase;
      Array<double,1> number, vsol;
      number.resize(config.nbins);
      vsol.resize(config.nbins);
      for (b = 0; b < config.nbins; ++b)
        {
	  config.diameters(b)=DSD[b];
	  number(b) = q[b]; // 1/m3
	  vsol(b) = csol[b] / 1400.0 * 1.0e-9; // in m3
	  // Compute fractions of sulfate ions.
	  surrogate[config.iSO4mm].Aaq_bins_init(b) = liquid_Nsize[4,b]*1.0e6*96.0; // SO4_2-
          surrogate[config.iHSO4m].Aaq_bins_init(b) = liquid_Nsize[5,b]*1.0e6*97.0; // HSO4-
          surrogate[config.iNO3m].Aaq_bins_init(b) = liquid_Nsize[6,b]*1.0e6*62.0; // NO3-
          surrogate[config.iNH4p].Aaq_bins_init(b) = liquid_Nsize[2,b]*1.0e6*18.0; // NH4+
          surrogate[config.iHp].Aaq_bins_init(b) = liquid_Nsize[0,b]*1.0e6; // H+
	}

      Array<double,1> total_so4_bins, frac_SO4mm_bins, frac_HSO4m_bins; 
      total_so4_bins.resize(config.nbins);
      frac_SO4mm_bins.resize(config.nbins);
      frac_HSO4m_bins.resize(config.nbins);
      for (b = 0; b < config.nbins; ++b)
        {
	  total_so4_bins(b) = surrogate[config.iSO4mm].Aaq_bins_init(b) + 
	    surrogate[config.iHSO4m].Aaq_bins_init(b);
	  if (total_so4_bins(b) == 0.0)
	    {
	      frac_SO4mm_bins(b) = 0.0;
	      frac_HSO4m_bins(b) = 0.0;
	    }  
	  else
	    {
	      frac_SO4mm_bins(b) = surrogate[config.iSO4mm].Aaq_bins_init(b) / total_so4_bins(b);
	      frac_HSO4m_bins(b) = 1.0 - frac_SO4mm_bins(b);
	    }
	  if (NaCl)
	    for (i = 0; i < n; ++i) // Na+ and Cl-
	      if ((surrogate[i].name == "Na") or (surrogate[i].name == "Cl"))
		{
		  int iq_aero = (surrogate[i].soap_ind + 1) * config.nbins; 
		  surrogate[i].Aaq_bins_init(b) = q[iq_aero + b];
		}
        }

      for (i = 0; i < n; ++i)
        {
          int iq = surrogate[i].soap_ind;
          if (iq != -1)
            {
              if (surrogate[i].is_organic)
                {             
		  surrogate[i].Ag = qgas[iq];
		  int iq_aero = (surrogate[i].soap_ind + 1) * config.nbins; 
                  for (b = 0; b < config.nbins; ++b)
                    {
                      if (surrogate[i].hydrophilic)
                        {
			  surrogate[i].Aaq_bins_init(b) = q[iq_aero + b];
			  surrogate[i].Aaq+=surrogate[i].Aaq_bins_init(b);
			}
		      else
                        surrogate[i].Aaq_bins_init(b) = 0.0;

                      for (ilayer=0;ilayer<config.nlayer;++ilayer)
                        for (iphase=0;iphase<config.max_number_of_phases;++iphase)
                          {
                            if (surrogate[i].hydrophobic and iphase==0)
                              {
				surrogate[i].Ap_layer_init(b,ilayer,iphase)=
				  config.Vlayer(ilayer) * q[iq_aero + b];
				surrogate[i].Ap+=surrogate[i].Ap_layer_init(b,ilayer,iphase);
                              }
                            else
                              surrogate[i].Ap_layer_init(b,ilayer,iphase)=0.0;
                          }
                    }
                }
              else if (NaCl and (surrogate[i].name == "Na" or surrogate[i].name == "Cl"))
                {
                  int iq_aero = (surrogate[i].soap_ind + 1) * config.nbins; 
                  surrogate[i].Aaq_bins_init(b) = q[iq_aero + b];

                }
              else
                {
                  surrogate[i].Aaq_bins_init = 0.0;
                  surrogate[i].Ap_layer_init = 0.0;
                }
             
	      // Inorganic gas concentrations
              if (surrogate[i].is_inorganic_precursor)
		{
		  surrogate[i].Ag = qgas[iq];
           
		}
            }
          else // liquid ions
            {
              for (b = 0; b < config.nbins; ++b)
                {
                  int ind;
                  // PMD    PBC     PNA    PSO4    PNH4    PNO3
                  // PHCL   PBiA2D  PBiA1D PBiA0D  PAGLY   PAMGLY ...
                  if (surrogate[i].name == "SO4")
                    {
                      ind = 3;
                      int iq_aero = (ind + 1) * config.nbins;                                   
                      surrogate[i].Aaq_bins_init(b) = frac_HSO4m_bins(b) * q[iq_aero + b];        
                    }
                  else if (surrogate[i].name == "HSO4")
                    {
                      ind = 3;
                      int iq_aero = (ind + 1) * config.nbins;                                   
                      surrogate[i].Aaq_bins_init(b) = frac_SO4mm_bins(b) * q[iq_aero + b];        
                    }
                  else if (surrogate[i].name == "NO3")
                    {
                      ind = 5;
                      int iq_aero = (ind + 1) * config.nbins; 
                      surrogate[i].Aaq_bins_init(b) = q[iq_aero + b];
                    }
                  else if (surrogate[i].name == "NH4")
                    {
                      ind = 4;
                      int iq_aero = (ind + 1) * config.nbins; 
                      surrogate[i].Aaq_bins_init(b) = q[iq_aero + b];
                    }
                  else
                    surrogate[i].Aaq_bins_init = 0.0;
                }
              surrogate[i].Ap_layer_init = 0.0;
            }
	  surrogate[i].Atot = surrogate[i].Ap + surrogate[i].Ag + surrogate[i].Aaq;

        }

      // Input data for the aqueous phase in each bin 
      Array<double,1> LWC_bins,AQinit_bins,chp_bins,ionic_bins;
      LWC_bins.resize(config.nbins);
      AQinit_bins.resize(config.nbins);
      chp_bins.resize(config.nbins);
      ionic_bins.resize(config.nbins);
      int iq_h2o = (surrogate[config.iH2O].soap_ind + 1) * config.nbins; 
      for (b=0;b<config.nbins;++b)
        {
          LWC_bins(b) = q[iq_h2o + b]; // or lwc_Nsize[b] - KS
          AQinit_bins(b) = LWC_bins(b);
          for (i=0;i<n;++i)
            AQinit_bins(b)+=surrogate[i].Aaq_bins_init(b);
          chp_bins(b)=chp_Nsize[b];
          ionic_bins(b)=ionic_Nsize[b];
        }
	  
      // Input data for the organic phase in each bin 
      Array<double,3> MOinit_layer,MOW_layer;
      MOinit_layer.resize(config.nbins,config.nlayer,config.max_number_of_phases);
      MOW_layer.resize(config.nbins,config.nlayer,config.max_number_of_phases);
      for (b=0;b<config.nbins;++b)
        for (ilayer=0;ilayer<config.nlayer;++ilayer)
          {
            for (iphase=0;iphase<config.max_number_of_phases;++iphase)
              MOinit_layer(b,ilayer,iphase)=0.0;
            for (i=0;i<n;++i)
              if(surrogate[i].hydrophobic)
                for (iphase=0;iphase<config.max_number_of_phases;++iphase)
		  MOinit_layer(b,ilayer,iphase)+=surrogate[i].Ap_layer_init(b,ilayer,iphase);
          }

      density_aqueous_phase(config, surrogate, LWC_bins, Temperature);

      dynamic_system(config,surrogate,
                     MOinit_layer, MOW_layer,number,vsol,
                     LWC_bins, AQinit_bins, ionic_bins, chp_bins,
                     Temperature, RH, deltat);

      // Give back the concentrations
      for (i = 0; i < n; ++i)
        {
          int iq = surrogate[i].soap_ind;
          if (iq != -1)
            {
              if (surrogate[i].is_organic)
                {              
                  int iq_gas = (ns_aer + 1) * config.nbins + surrogate[i].soap_ind;
                  int iq_aero = (surrogate[i].soap_ind + 1) * config.nbins; 
                  q[iq_gas] = surrogate[i].Ag;
                  for (b = 0; b < config.nbins; ++b)
                    {
                      q[iq_aero + b] = surrogate[i].Aaq_bins(b);
                      for (ilayer=0;ilayer<config.nlayer;++ilayer)
                        for (iphase=0;iphase<config.nphase(b,ilayer);++iphase)
                          q[iq_aero + b] += surrogate[i].Ap_layer(b,ilayer,iphase);
                    }
                }

              else if (i == config.iH2O)
		for (b = 0; b < config.nbins; ++b)
		  {
		    q[iq_h2o + b] = surrogate[i].Aaq_bins(b);
		    for (ilayer=0;ilayer<config.nlayer;++ilayer)
		      for (iphase=0;iphase<config.nphase(b,ilayer);++iphase)
			q[iq_h2o + b] += surrogate[i].Ap_layer(b,ilayer,iphase);
		  }
            }
        }

      for (b = 0; b < config.nbins; ++b)
        {
	  DSD[b] = config.diameters(b);
	  q[b] = number(b);
	  q[iq_h2o + b] += LWC_bins(b);           
        }
    }

  return;

}
