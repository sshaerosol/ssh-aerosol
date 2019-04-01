#include "properties.hxx"
#include "aiomfac.cxx"
#include "unifac.cxx"
using namespace soap;

void compute_gamma_infini(model_config &config, vector<species>& surrogate)
{
  // This routine compute the activity coefficients at infinite dilution GAMMAinf
  //     given by UNIFAC for each species and compute the Henry's law constant from
  //     the saturation vapour pressure if the Henry's law constant was not specified
  //
  //   H = 1000*760/(18.0*GAMMAinf*Psat)
  
  int n=surrogate.size();
  int i,j;
  Array<double, 1> X_unifac;
  Array<double, 1> gamma_unifac;
  double tiny=1.0e-100;
  X_unifac.resize(config.nmol_aq);
  gamma_unifac.resize(config.nmol_aq);
  for (i=0;i<n;++i)
    if(surrogate[i].hydrophilic and surrogate[i].compute_gamma_aq)
      if (surrogate[i].is_organic and surrogate[i].nonvolatile==false)
        {
          for (j=0;j<config.nmol_aq;++j)
            X_unifac(j)=0.0;
          X_unifac(surrogate[i].index_gamma_aq)=1.0e-10;
          X_unifac(surrogate[config.iH2O].index_gamma_aq)=
            1.0-X_unifac(surrogate[i].index_gamma_aq);
          unifac(config.nmol_aq,config.nfunc_aq,config.groups_aq,
                 X_unifac,config.Inter_aq,config.InterB_aq,config.InterC_aq,
		 config.RG_aq,
                 config.QG_aq,config.Rparam_aq,config.Qparam_aq,config.Lparam_aq,
                 config.Z,surrogate[i].Tref,gamma_unifac,config.temperature_dependancy);
          surrogate[i].GAMMAinf=gamma_unifac(surrogate[i].index_gamma_aq);
          
          if (surrogate[i].Henry <= tiny)
            surrogate[i].Henry=1000.0*760.0/(18.0*surrogate[i].GAMMAinf*surrogate[i].Psat_ref);
        }
      else
        {
          surrogate[i].GAMMAinf=1.0;
        }
    else
      surrogate[i].GAMMAinf=1.0;
}

void compute_density_aqueous_phase(model_config &config, vector<species>& surrogate, double &LWC, double &Temperature)
{
  //Compute density according to parameterisation in:
  //  Semmler, M., Luo, B. P. and Koop, T. Densities of liquid H+/NH4+/SO42-/NO3-/H2O solutions at tropospheric temperatures. Atmos. Environ., 40, 2006, 467-483.

  int n=surrogate.size();
  int i;
  double sumW;
  double temp,hcl;
  temp=0.0;
  sumW=0.0;
  
  for (i=0;i<n;++i)
    if (surrogate[i].hydrophilic and surrogate[i].is_organic)
      {
        surrogate[i].Waq=surrogate[i].Aaq;
        if (i==config.iH2O)
          surrogate[i].Waq+=LWC;
        sumW+=surrogate[i].Waq;
      }
    else if (surrogate[i].hydrophilic and i==config.iH2O)
      {
        surrogate[i].Waq=surrogate[i].Aaq+LWC;
        sumW+=surrogate[i].Waq;
      }
  
  double avail_nh4,avail_na;
  if (config.iNH4p>=0)
    {
      avail_nh4=surrogate[config.iNH4p].Aaq;
    }
  else 
    avail_nh4=0.0;

  if (config.iNa>=0)
    {
      avail_na=surrogate[config.iNa].Aaq;
    }
  else 
    avail_na=0.0;

  double waq_NaHSO4=0.0;
  double waq_NH4HSO4=0.0;
  double waq_H2SO4=0.0;
  double waq_Na2SO4=0.0;
  double waq_NH4NH4SO4=0.0;
  double waq_NaNO3=0.0;
  double waq_NH4NO3=0.0;
  double waq_HNO3=0.0;
  double waq_NaCl=0.0;
  double waq_NH4Cl=0.0;
  double waq_HCl=0.0;

  if (config.iHSO4m>=0)
    {
      double avail_hso4=surrogate[config.iHSO4m].Aaq;

      if (avail_hso4>0.0 and avail_na>0.0)
        {
          waq_NaHSO4=120.0*min(avail_hso4/surrogate[config.iHSO4m].MM,
                               avail_na/surrogate[config.iNa].MM);
          avail_hso4-=waq_NaHSO4/120.0*surrogate[config.iHSO4m].MM;
          avail_na-=waq_NaHSO4/120.0*surrogate[config.iNa].MM;
        }
      
      if (avail_hso4>0.0 and avail_nh4>0.0)
        {
          waq_NH4HSO4=115.0*min(avail_hso4/surrogate[config.iHSO4m].MM,
                                avail_nh4/surrogate[config.iNH4p].MM);
          avail_hso4-=waq_NH4HSO4/115.0*surrogate[config.iHSO4m].MM;
          avail_nh4-=waq_NH4HSO4/115.0*surrogate[config.iNH4p].MM;
        }

      if (avail_hso4>0.0)
        waq_H2SO4=98.0*avail_hso4/surrogate[config.iHSO4m].MM;  
    }

  if (config.iSO4mm>=0)
    {      
      double avail_so4=surrogate[config.iSO4mm].Aaq;

      if (avail_so4>0.0 and avail_na>0.0)
        {
          waq_Na2SO4=142.0*min(avail_so4/surrogate[config.iSO4mm].MM,
                               0.5*avail_na/surrogate[config.iNa].MM);
          avail_so4-=waq_Na2SO4/142.0*surrogate[config.iSO4mm].MM;
          avail_na-=2.0*waq_Na2SO4/142.0*surrogate[config.iNa].MM;
        }

      if (avail_so4>0.0 and avail_nh4>0.0)
        {
          waq_NH4NH4SO4=132.0*min(avail_so4/surrogate[config.iSO4mm].MM,
                                  0.5*avail_nh4/surrogate[config.iNH4p].MM);
          avail_so4-=waq_NH4NH4SO4/132.0*surrogate[config.iSO4mm].MM;
          avail_nh4-=2.0*waq_NH4NH4SO4/132.0*surrogate[config.iNH4p].MM;
        }

      if (avail_so4>0.0)
        waq_H2SO4+=98.0*avail_so4/surrogate[config.iSO4mm].MM;  

    }

  if (config.iNO3m>=0)
    {      
      double avail_no3=surrogate[config.iNO3m].Aaq;

      if (avail_no3>0.0 and avail_na>0.0)
        {
          waq_NaNO3=85.0*min(avail_no3/surrogate[config.iNO3m].MM,
                             avail_na/surrogate[config.iNa].MM);
          avail_no3-=waq_NaNO3/85.0*surrogate[config.iNO3m].MM;
          avail_na-=waq_NaNO3/85.0*surrogate[config.iNa].MM;
        }

      if (avail_no3>0.0 and avail_nh4)
        {
          waq_NH4NO3=80.0*min(avail_no3/surrogate[config.iNO3m].MM,
                              avail_nh4/surrogate[config.iNH4p].MM);
          avail_no3-=waq_NH4NO3/80.0*surrogate[config.iNO3m].MM;
          avail_nh4-=waq_NH4NO3/80.0*surrogate[config.iNH4p].MM;
        }

      if (avail_no3>0.0)
        waq_HNO3=63.0*avail_no3/surrogate[config.iNO3m].MM;  
    }

  if (config.iClm>=0)
    {      
      double avail_cl=surrogate[config.iClm].Aaq;

      if (avail_cl>0.0 and avail_na>0.0)
        {
          waq_NaCl=58.5*min(avail_cl/surrogate[config.iClm].MM,
                             avail_na/surrogate[config.iNa].MM);
          avail_cl-=waq_NaCl/58.5*surrogate[config.iClm].MM;
          avail_na-=waq_NaCl/58.5*surrogate[config.iNa].MM;
        }

      if (avail_cl>0.0 and avail_nh4)
        {
          waq_NH4Cl=53.5*min(avail_cl/surrogate[config.iClm].MM,
                              avail_nh4/surrogate[config.iNH4p].MM);
          avail_cl-=waq_NH4Cl/53.5*surrogate[config.iClm].MM;
          avail_nh4-=waq_NH4Cl/53.5*surrogate[config.iNH4p].MM;
        }

      if (avail_cl>0.0)
        waq_HCl=36.5*avail_cl/surrogate[config.iClm].MM;  
    }

  double waq_Na=avail_na;
  double waq_NH3=17.0*avail_nh4/surrogate[config.iNH4p].MM;

  sumW+=waq_NaHSO4+waq_NH4HSO4+waq_H2SO4+waq_Na2SO4+waq_NH4NH4SO4+waq_NaNO3
    +waq_NH4NO3+waq_HNO3+waq_NaCl+waq_NH4Cl+waq_HCl+waq_Na+waq_NH3+surrogate[config.iH2O].Aaq+LWC;

  double Tr=(Temperature-273.15)/273.15;
  double temp_rho;
  if (sumW>0.0)
    {
      for (i=0;i<n;++i)
        if (surrogate[i].hydrophilic and surrogate[i].is_organic)
          {
            surrogate[i].Waq/=sumW;
            temp+=surrogate[i].Waq/surrogate[i].rho;
          }
        else if (i==config.iH2O)
          {
            surrogate[i].Waq/=sumW;
            temp_rho=18.0/(17.89+1.815*Tr)*1000;
            temp+=surrogate[config.iH2O].Waq/temp_rho;
          }
      
      double w;
      //H2SO4
      w=waq_H2SO4/sumW;
      temp_rho=98.0/(36.2249+5.46891*w+8.73468*w*w+
                34.4556*Tr-65.0688*Tr*w+41.7174*Tr*w*w+
                13.6245*Tr*Tr-39.1920*Tr*Tr*w+25.2329*Tr*Tr*w*w)*1000;
      temp+=w/temp_rho;

      //(NH4)2SO4
      w=waq_NH4NH4SO4/sumW;
      temp_rho=132.0/(55.8085+40.6200*w-11.0797*w*w+
                      13.5837*Tr-15.8075*Tr*w)*1000;
      temp+=w/temp_rho;

      //HNO3
      w=waq_HNO3/sumW;
      temp_rho=63.0/(28.6889-1.35400*w+15.9166*w*w+
                     16.9281*Tr-1.36496*Tr*w)*1000;
      temp+=w/temp_rho;

      //NH4NO3
      w=waq_NH4NO3/sumW;
      temp_rho=80.0/(50.3779-6.91069*w+11.3024*w*w
                     -28.7432*Tr+155.419*Tr*w-157.049*Tr*w*w+
                     161.094*Tr*Tr-585.880*Tr*Tr*w+541.843*Tr*Tr*w*w)*1000;
      temp+=w/temp_rho;
                     
      //NH4HSO4
      w=waq_NH4HSO4/sumW;
      temp_rho=115.0/(47.6267+28.2931*w-8.90139*w*w
                      +32.3616*Tr-45.8649*Tr*w)*1000;
      temp+=w/temp_rho;
      
      //NaHSO4
      temp+=waq_NaHSO4/2740.0/sumW;
      
      //Na2SO4
      temp+=waq_Na2SO4/2700.0/sumW;
      
      //NaNO3
      temp+=waq_NaNO3/2260.0/sumW;

      //NaCl
      temp+=waq_NaCl/2165.0/sumW;

      //NH4Cl
      temp+=waq_NH4Cl/1530.0/sumW;

      //HCl
      temp+=waq_HCl/1150.0/sumW;

      //Na
      temp+=waq_Na/970.0/sumW;
      
      //NH3
      temp+=waq_NH3/910.0/sumW;

      //H2O
      temp+=(surrogate[i].Aaq+LWC)/1000.0/sumW;
    }
  

  if (temp>0.0)
    config.rho_aqueous=1.0/temp;
  else
    config.rho_aqueous=1000.0;
}

void density_aqueous_phase(model_config &config, vector<species>& surrogate,
                           Array<double, 1> &LWC, double Temperature)
{
  //compute the density of the aqueous phase based on the massic fraction of water
  //assumes that the volume of the aqueous phase is equal to the volume of water +
  //the volume of other compounds
  int n=surrogate.size();
  int b,i;
  
  if (config.compute_rho_aqueous)
    {
      for (b=0;b<config.nbins;++b)
        {
          for (i=0;i<n;++i)
            if (surrogate[i].hydrophilic)
              surrogate[i].Aaq=surrogate[i].Aaq_bins_init(b);
   
          compute_density_aqueous_phase(config, surrogate, LWC(b), Temperature);
          config.AQrho(b)=config.rho_aqueous;
        }
    }
  else
    for (b=0;b<config.nbins;++b)
      config.AQrho(b)=config.rho_aqueous;
}

void compute_viscosity(model_config &config, vector<species>& surrogate,
		       int &b, int &ilayer)
{
  int i,iphase;
  int n=surrogate.size();

  if (config.constant_dorg) //the organic phase coefficient diffusion is assumed constant
    {
      if (config.explicit_representation)
	{
	  for (i=0;i<n;++i)
	    if (surrogate[i].is_organic)
	      surrogate[i].dif_org(b,ilayer)=config.dorg;
	    else
	      surrogate[i].dif_org(b,ilayer)=1.0e-12;
	}
      else
	{
	  for (i=0;i<n;++i)
	    if (surrogate[i].is_organic)
	      surrogate[i].KDiffusion_p=config.dorg;
	    else
	      surrogate[i].KDiffusion_p=1.0e-12;
	}
    }
  else
    {
      //In this method the organic phase coefficient diffusion is computed from the viscosity of the mixture
      //The viscosity of the mixture is computed from the viscosity of each compound with a VBN method
      double molecular_radius=1.0e-10;
      double mo=0.0;
      for (i=0;i<n;++i)
	if (surrogate[i].hydrophobic)
	  for (iphase=0;iphase<config.nphase(b,ilayer);++iphase)
	    mo+=surrogate[i].Ap_layer_init(b,ilayer,iphase);	       
      
      double dorg;
      if (mo>config.MOmin*config.Vlayer(ilayer))
	{
	  double VBN=0.0;
	  double volume=0.0;
	  for (i=0;i<n;++i)
	    if (surrogate[i].hydrophobic)
	      for (iphase=0;iphase<config.nphase(b,ilayer);++iphase)
		{
		  VBN+=surrogate[i].Ap_layer_init(b,ilayer,iphase)/mo*(14.534*log(log(surrogate[i].viscosity+0.8))+10.975);
		  volume+=surrogate[i].Ap_layer_init(b,ilayer,iphase)/surrogate[i].rho;
		}

	  double rho=mo/volume;
	  double visco=(exp(exp((VBN-10.975)/14.534))-0.8)*1.0e-6*rho;
	  if (config.explicit_representation)
	    {
	      for (i=0;i<n;++i)
		if (surrogate[i].is_organic)
		  surrogate[i].dif_org(b,ilayer)=1.38e-23*298.0/(6.0*3.14159*molecular_radius*visco);
		else
		  surrogate[i].dif_org(b,ilayer)=1.0e-12;
	    }
	  else
	    {
	      for (i=0;i<n;++i)
		if (surrogate[i].is_organic)
		  surrogate[i].KDiffusion_p=1.38e-23*298.0/(6.0*3.14159*molecular_radius*visco);
		else
		  surrogate[i].KDiffusion_p=1.0e-12;
	    }
	}
      else
	if (config.explicit_representation)
	  {
	    for (i=0;i<n;++i)
	      surrogate[i].dif_org(b,ilayer)=1.0e-12;
	  }
	else
	  {
	    for (i=0;i<n;++i)
	      surrogate[i].KDiffusion_p=1.0e-12;
	  }
      
    }
}

void water_concentration(model_config &config, vector<species>& surrogate,
                         double &Temperature, double &RH)
{
  //compute the gas phase concentration of water
  double Pwater=surrogate[config.iH2O].Psat(Temperature)*RH;
  surrogate[config.iH2O].Atot=(Pwater/760.0*1.013e5)*surrogate[config.iH2O].MM*1.0e6/
    (8.314*Temperature);
  surrogate[config.iH2O].Ag=max(surrogate[config.iH2O].Atot
				-sum(surrogate[config.iH2O].Aaq_bins_init)
				-sum(surrogate[config.iH2O].Ap_layer_init)
				,0.0);
  
}

double species::knudsen_function(double &Temperature, double diameter)
{
  //compute the value of the knudsen function
  double velocity=pow(2.11714271498563e4*Temperature/MM,0.5);
  double knudsen=4.0*KDiffusion_air/(velocity*diameter);
  double value=0.75*accomodation_coefficient*(1.0+knudsen)/
    (pow(knudsen,2)+knudsen+0.283*knudsen*accomodation_coefficient+0.75*accomodation_coefficient);
  return value;
}


void tau_kmt(model_config &config,vector<species>& surrogate, double &temperature,
             Array<double, 1> &number)
{
  //compute the characteristic time for condensation
  int b,i;
  int n=surrogate.size();
  double pi=3.14159265358979323846;
  double R=8.314;

  for (b=0;b<config.nbins;++b)	  
    for (i=0;i<n;++i)
      if (number(b)>0.0)        
	surrogate[i].tau_air(b)=1.0/
	  (2.0*pi*surrogate[i].KDiffusion_air*
	   config.diameters(b)*1.0e-6*number(b)*
	   surrogate[i].knudsen_function(temperature, config.diameters(b)*1.0e-6));	         
      else
        surrogate[i].tau_air(b)=1.0e19; 
}

void tau_dif(model_config &config, vector<species>& surrogate,
             Array<double, 1> &number, Array<double, 1> &Vsol)
{
  //compute the characteristic time for diffusion in the organic phase
  int b,i,ilayer,iphase;
  int n=surrogate.size();
  double pi=3.14159265358979323846;
  double fs;
  double morphology_factor;
  //morphology_factor is a factor to take into account the morphology of the particle and that
  // particles are not entirely organic
  
  //here it is assumed that the center of the particle is solid and that the aqueous phase
  // does not impact the characteristic time for diffusion in the organic phase
  //the particle is still assumed spheric
  //fs: fraction of the volume which is solid
  
  for (b=0;b<config.nbins;++b)
    {
      if (number(b)>0.0)
        fs=(Vsol(b)/number(b))/(pi/6*pow(config.diameters(b)*1.0e-6,3));
      else
        fs = 0.0;
      for (ilayer=config.nlayer-1;ilayer>=0;--ilayer)
	for (iphase=0;iphase<config.nphase(b,ilayer);++iphase)
	  {
	    //determine the organic phase coefficient diffusion
            compute_viscosity(config, surrogate, b, ilayer);

	    morphology_factor=1.0+config.Alayer(ilayer,0)*pow(fs,4)+config.Alayer(ilayer,1)*pow(fs,3)
	      +config.Alayer(ilayer,2)*pow(fs,2)+config.Alayer(ilayer,3)*fs;
	    morphology_factor=max(morphology_factor,1.0e-3);
	    if (ilayer==config.nlayer-1)
	      for (i=0;i<n;++i)
		surrogate[i].tau_diffusion(b,ilayer,iphase)=pow(config.diameters(b)*1.0e-6,2)/
		  (4.0*pi*pi*surrogate[i].KDiffusion_p*config.alpha_layer(ilayer))
		  *morphology_factor;
	    else
	      {
		 double morphology_factor2=1.0+config.Alayer(ilayer+1,0)*pow(fs,4)+config.Alayer(ilayer+1,1)*pow(fs,3)
		   +config.Alayer(ilayer+1,2)*pow(fs,2)+config.Alayer(ilayer+1,3)*fs;
		 morphology_factor2=max(morphology_factor2,1.0e-3);

		for (i=0;i<n;++i)
		  surrogate[i].tau_diffusion(b,ilayer,iphase)=pow(config.diameters(b)*1.0e-6,2)/
		    (4.0*pi*pi*surrogate[i].KDiffusion_p*config.alpha_layer(ilayer))
		    *morphology_factor
		    -pow(config.diameters(b)*1.0e-6,2)/
		    (4.0*pi*pi*surrogate[i].KDiffusion_p*config.alpha_layer(ilayer+1))
		    *morphology_factor2
		    +surrogate[i].tau_diffusion(b,ilayer+1,iphase);
	      }	      
	  }
    }
}
  
double species::Psat(double &Temperature)
{
  //Compute the saturation vapor pressure at the specified temperature 
  double R=8.314; //ideal gas constant (J/K/mol)
  double value;
  if (name=="H2O")
    value=exp(13.7-5120.0/Temperature)*760.0;
  else
    value=Psat_ref*exp(-1000.0*deltaH/R*(1.0/Temperature-1.0/Tref));
  return value;
}

double species::Kpart_org(double &Temperature, double &MOW)
{
  //Compute the ideal partitioning constant between the gas phase and the organic phase
  //Ref: Odum et al. (Environ. Sci. Technol. 1996, vol=30)
  double R=8.206e-5; //ideal gas constant (m3.atm/mol/K)
  double value=760.0*R*Temperature/(MOW*1.0e6*Psat(Temperature));
  return value;
}

double species::Kp_eff_org(double &Temperature, double &MOW)
{
  //Compute the effective ideal partitioning constant between the gas phase
  //and the organic phase
  double value=Kpart_org(Temperature, MOW)*(1.0+Koligo_org);
  return value;
}


double species::Kp_exp_org(double &Temperature)
{
  //Compute the experimental partitioning constant between the gas phase and the organic phase
  double R=8.314; //ideal gas constant (J/K/mol)
  double value=kp_experiment*Temperature/Tref*exp(1000.0*deltaH/R*(1.0/Temperature-1.0/Tref));
  return value;
}

double species::Kpart_aq(double &Temperature, double &MMaq)
{
  //Compute the ideal partitioning constant between the gas phase and the aqueous phase
  double R=8.314; //ideal gas constant (J/K/mol)
  double rho_h2o=1000.0; //volumic mass of H2O
  double MH2O=18.0;
  double value=Henry*R*Temperature/(rho_h2o*1.0e6*1.013e5)*
    exp(1000.0*deltaH/R*(1.0/Temperature-1.0/Tref))*MH2O/MMaq;
  return value;
}

double species::Kp_eff_aq(model_config &config, double &Temperature, double &ionic, double &chp,
                          double &gammaH_LR, double &gammaH_SRMR,
                          double &MMaq, double &fion1, double &fion2)
{
  //Compute the effective ideal partitioning constant between the gas phase
  //and the organic phase
  // chp: concentrations of H+ ions in mol/L
  // ionic: ionic strenght
  // gammaH_LR: activity coefficients of H+ due to long range interactions
  // gammaH_SRMR: activity coefficients of H+ due to medium range and short range interactions
  // fion1: fraction of an acid H2A or HA which has been dissociated into HA- or A-
  // fion2: fraction of a diacid H2A which has been dissociated into A2-
  // MMaq: mean molar mass of the aqueous phase (g/mol)
  double value;
  double gamma;
  double XH2O_crit=0.2;
  fion1=0.0;
  fion2=0.0;
  if (aq_type=="diacid")
    {
      //Kp_effective=Kp_theoric*(1+HA-/H2A*(1+A2-/H2A))
      if (config.compute_aqueous_phase_properties) //config.compute_long_and_medium_range_interactions)
        {
          //For a species H2A:
          //gamma(H2A)=gamma_LR(H2A)*gamma_MR(H2A)*gamma_SR(H2A)
          //gamma(HA-)=gamma_LR(HA-)*gamma_MR(HA-)*gamma_SR(HA-)
          //gamma(A2-)=gamma_LR(A2-)*gamma_MR(A2-)*gamma_SR(A2-)
          //gamma(H+)=gamma_LR(H+)*gamma_MR(H+)*gamme_LR(H+)
          //Hypothesis:
          //gamma_SR(H2A)=gamma_SR(HA-)=gamma_SR(A2-)
          //gamma_MR(H2A)=gamma_MR(HA-)=gamma_MR(A2-)
          //gamma_LR(HA-)=gamma_LR(H+)
          //gamma_LR(A2-)=pow(gamma_LR(H+),2)
          double ratio_gamma1=pow(gammaH_LR,2.0)*gammaH_SRMR/gamma_LR;
          double ratio_gamma2=pow(gammaH_LR,2.0)*gammaH_SRMR;
          value=Kpart_aq(Temperature, MMaq)*
            (1.0+Kacidity1/(ratio_gamma1*chp)*(1.0+Kacidity2/(ratio_gamma2*chp)));
          fion1=(Kacidity1/(ratio_gamma1*chp))/
            (1.0+Kacidity1/(ratio_gamma1*chp)*(1.0+Kacidity2/(ratio_gamma2*chp)));
          fion2=(Kacidity1/(ratio_gamma1*chp))*(Kacidity2/(ratio_gamma2*chp))/
            (1.0+Kacidity1/(ratio_gamma1*chp)*(1.0+Kacidity2/(ratio_gamma2*chp)));
        }
      else
        {
          gamma=pow(10,-0.511*pow(298.0/Temperature,1.5)*pow(ionic,0.5)/(1.0+pow(ionic,0.5)));
          value=Kpart_aq(Temperature, MMaq)*
            (1.0+Kacidity1/(pow(gamma,2)*chp)*
             (1.0+Kacidity2/(pow(gamma,2)*chp)));
          fion1=(Kacidity1/(pow(gamma,2)*chp))/
            (1.0+Kacidity1/(pow(gamma,2)*chp)*(1.0+Kacidity2/(pow(gamma,2)*chp)));
          fion2=(Kacidity1/(pow(gamma,2)*chp))*(Kacidity2/(pow(gamma,2)*chp))/
            (1.0+Kacidity1/(pow(gamma,2)*chp)*(1.0+Kacidity2/(pow(gamma,2)*chp)));
        }
    }
  else if (aq_type=="monoacid")
    {
      //Kp_effective=Kp_theoric*(1+A-/HA)
      if (config.compute_aqueous_phase_properties) //config.compute_long_and_medium_range_interactions)
        {
          double ratio_gamma=pow(gammaH_LR,2.0)*gammaH_SRMR/gamma_LR;
          value=Kpart_aq(Temperature, MMaq)*(1.0+Kacidity1/(ratio_gamma*chp));
          fion1=(Kacidity1/(ratio_gamma*chp))
            /(1.0+Kacidity1/(ratio_gamma*chp));
        }
      else
        {
          gamma=pow(10,-0.511*pow(298.0/Temperature,1.5)*pow(ionic,0.5)/(1.0+pow(ionic,0.5)));
          value=Kpart_aq(Temperature, MMaq)*(1.0+Kacidity1/(pow(gamma,2)*chp));
          fion1=(Kacidity1/(pow(gamma,2)*chp))
            /(1.0+Kacidity1/(pow(gamma,2)*chp));
        }
    }
  else if (aq_type=="aldehyde")
    {
      //effective partitioning based on Pun and Seigneur (2007)
      if (config.compute_aqueous_phase_properties) //config.compute_long_and_medium_range_interactions)
        {
          value=Kpart_aq(Temperature, MMaq)*
            (1.0+Koligo_aq*pow(gammaH_LR*gammaH_SRMR*chp/pow(10,-pHref),beta));
        }
      else
        {
          gamma=pow(10,-0.511*pow(298.0/Temperature,1.5)*pow(ionic,0.5)/(1.0+pow(ionic,0.5)));
          value=Kpart_aq(Temperature, MMaq)*
            (1.0+Koligo_aq*pow(gamma*chp/pow(10,-pHref),beta));
        }
    }
  else if (aq_type=="none") //Kp_effective=Kp_theoric
    value=Kpart_aq(Temperature, MMaq);
  else
    {
      if (name!="H2O")
        cout << "WARNING: aq_type "+aq_type+" of species " +name+ " not defined." << endl;
      value=Kpart_aq(Temperature, MMaq);
    }
  
  return value;
}

double species::Kp_eff_aqreal(model_config &config,
                              double &Temperature, double &ionic, double &chp,
                              double &gammaH_LR, double &gammaH_SRMR,
                              double &MMaq, double &fion1, double &fion2)
{
  //Compute the effective ideal partitioning constant between the gas phase
  //and the organic phase
  // chp: concentrations of H+ ions in mol/L
  // ionic: ionic strenght
  // gammaH_LR: activity coefficients of H+ due to long range interactions
  // gammaH_SRMR: activity coefficients of H+ due to medium range and short interactions
  // fion1: fraction of an acid H2A or HA which has been dissociated into HA- or A-
  // fion2: fraction of a diacid H2A which has been dissociated into A2-
  // MMaq: mean molar mass of the aqueous phase (g/mol)
  double value;
  double gamma;
  double XH2O_crit=0.2;
  fion1=0.0;
  fion2=0.0;
  if (aq_type=="diacid")
    {
      //Kp_effective=Kp_theoric*(1+HA-/H2A*(1+A2-/H2A))
      if (config.compute_aqueous_phase_properties) //config.compute_long_and_medium_range_interactions)
        {
          //For a species H2A:
          //gamma(H2A)=gamma_LR(H2A)*gamma_MR(H2A)*gamma_SR(H2A)
          //gamma(HA-)=gamma_LR(HA-)*gamma_MR(HA-)*gamma_SR(HA-)
          //gamma(A2-)=gamma_LR(A2-)*gamma_MR(A2-)*gamma_SR(A2-)
          //gamma(H+)=gamma_LR(H+)*gamma_MR(H+)*gamma_SR(H+)
          //Hypothesis:
          //gamma_SR(H2A)=gamma_SR(HA-)=gamma_SR(A2-)
          //gamma_MR(H2A)=gamma_MR(HA-)=gamma_MR(A2-)
          //gamma_LR(HA-)=gamma_LR(H+)
          //gamma_LR(A2-)=pow(gamma_LR(H+),2)
          double ratio_gamma1=pow(gammaH_LR,2.0)*gammaH_SRMR/gamma_LR;
          double ratio_gamma2=pow(gammaH_LR,2.0)*gammaH_SRMR;
          value=kaqi/MMaq*
            (1.0+Kacidity1/(ratio_gamma1*chp)*(1.0+Kacidity2/(ratio_gamma2*chp)));
          fion1=(Kacidity1/(ratio_gamma1*chp))/
            (1.0+Kacidity1/(ratio_gamma1*chp)*(1.0+Kacidity2/(ratio_gamma2*chp)));
          fion2=(Kacidity1/(ratio_gamma1*chp))*(Kacidity2/(ratio_gamma2*chp))/
            (1.0+Kacidity1/(ratio_gamma1*chp)*(1.0+Kacidity2/(ratio_gamma2*chp)));
        }
      else
        {          
          value=kaqi/MMaq;
          fion1=fioni1;
          fion2=fioni2;
        }
    }
  else if (aq_type=="monoacid")
    {
      //Kp_effective=Kp_theoric*(1+A-/HA)
      if (config.compute_aqueous_phase_properties) //config.compute_long_and_medium_range_interactions)
        {
          double ratio_gamma=pow(gammaH_LR,2.0)*gammaH_SRMR/gamma_LR;
          value=kaqi/MMaq*(1.0+Kacidity1/(ratio_gamma*chp));
          fion1=(Kacidity1/(ratio_gamma*chp))
            /(1.0+Kacidity1/(ratio_gamma*chp));
        }
      else
        {          
          value=kaqi/MMaq;
          fion1=fioni1;
        }
    }
  else if (aq_type=="aldehyde")
    {
      //effective partitioning based on Pun and Seigneur (2007)
      if (config.compute_aqueous_phase_properties) //config.compute_long_and_medium_range_interactions)        
        value=kaqi/MMaq*(1.0+Koligo_aq*pow(gammaH_LR*gammaH_SRMR*chp/pow(10,-pHref),beta));        
      else        
        value=kaqi/MMaq;        
    }
  else if (aq_type=="none") //Kp_effective=Kp_theoric
    value=kaqi/MMaq;
  else
    {
      if (name!="H2O")
        cout << "WARNING: aq_type "+aq_type+" of species " +name+ " not defined." << endl;
      value=kaqi/MMaq;
    }
  
  return value;
}

void activity_coefficients_org(model_config &config, vector<species>& surrogate,
                               bool all_hydrophobic, double &Temperature,
                               double &MOW)
{
  //compute the activity coefficients with UNIFAC (short range interactions) for the organic phase
  //MOW: mean molar mass of the organic phase
  //all_hydrophobic: do all the compounds condense on the organic?(If low mass of particulate water)
  int n_unifac,i;
  int n=surrogate.size();
  double sum,sumX_unifac;
  Array<double, 1> X_unifac,gamma_unifac;

  //initialization of X_unifac
  
  if (all_hydrophobic==false)
    {
      X_unifac.resize(config.nmol_org);
      gamma_unifac.resize(config.nmol_org);
      n_unifac=config.nmol_org;
    }
  else
    {
      X_unifac.resize(config.nmol_tot);
      gamma_unifac.resize(config.nmol_tot);
      n_unifac=config.nmol_tot;
    }
  
  sum=0.0;
  sumX_unifac=0.0;
  
  for (i=0;i<n_unifac;++i)
    X_unifac(i)=0.0;

  //computation of X_unifac (molar fraction of compounds used in the computation of activity
  // coefficients)
  for (i=0;i<n;++i)
    if (surrogate[i].hydrophobic or
        (all_hydrophobic and (surrogate[i].is_organic or i==config.iH2O)))
      {
        surrogate[i].Xorg=surrogate[i].Ap/surrogate[i].MM;
        sum+=surrogate[i].Xorg;
        if (all_hydrophobic==false)
          {
            if (surrogate[i].index_gamma_org>=0)
              {
                sumX_unifac+=surrogate[i].Xorg;
                X_unifac(surrogate[i].index_gamma_org)+=surrogate[i].Xorg;
              }
          }
        else
          if (surrogate[i].index_gamma_tot>=0)
            {
              sumX_unifac+=surrogate[i].Xorg;
              X_unifac(surrogate[i].index_gamma_tot)+=surrogate[i].Xorg;
            }
      }
  
  if (sum>0.0)
    {
      MOW=0.0;
      for (i=0;i<n;++i)
        if (surrogate[i].hydrophobic or
            (all_hydrophobic and (surrogate[i].is_organic or i==config.iH2O)))
          {
            surrogate[i].Xorg=surrogate[i].Xorg/sum;
            MOW=MOW+surrogate[i].Xorg*surrogate[i].MM;
          }
    }
  else
    MOW=200.0;

  //Call of unifac 
  if (sumX_unifac>0.0 and config.activity_model=="unifac")
    {
      for (i=0;i<n_unifac;++i)
        {
          X_unifac(i)/=sumX_unifac;
          X_unifac(i)=max(1.0e-11, X_unifac(i));
        }

      if (config.tabulation_unifac)
        {
          if (all_hydrophobic==false)
            {
              Array <double, 1> spec_activity;
              spec_activity.resize(config.nmol_org);
              int itemp=int(Temperature-263.15);
              itemp=max(itemp,0);
              itemp=min(itemp,config.ntemp);
              double alpha;
              if (Temperature<=263.15) 
                alpha=0.0;            
              else if (Temperature>=263.15+1.0*config.ntemp)
                alpha=1.0;
              else
                alpha=Temperature-263.15-1.0*itemp;                           

              for (i=0;i<n;i++)
                if (surrogate[i].index_gamma_org>=0)                                     
                  spec_activity(surrogate[i].index_gamma_org)=
                    (1.0-alpha)*surrogate[i].species_activity(itemp)+alpha*surrogate[i].species_activity(itemp+1);

              tab_unifac(config.nmol_org,config.nfunc_org,config.groups_org,
                         X_unifac,config.Inter_org,config.InterB_org,config.InterC_org,config.RG_org,
                         config.QG_org,config.Rparam_org,config.Qparam_org,config.Lparam_org,
                         spec_activity,config.Z,Temperature,gamma_unifac,config.temperature_dependancy);

              for (i=0;i<n;++i)
                {
                  if (surrogate[i].index_gamma_org >=0)
                    surrogate[i].gamma_org=gamma_unifac(surrogate[i].index_gamma_org);
                  else
                    surrogate[i].gamma_org=1.0;
                }
            }
          else
            {
              Array <double, 1> spec_activity;
              spec_activity.resize(config.nmol_tot);
              int itemp=int(Temperature-263.15);
              itemp=max(itemp,0);
              itemp=min(itemp,config.ntemp);
              double alpha;
              if (Temperature<=263.15) 
                alpha=0.0;            
              else if (Temperature>=263.15+1.0*config.ntemp)
                alpha=1.0;
              else
                alpha=Temperature-263.15-itemp*1.0;              

              for (i=0;i<n;i++)
                if (surrogate[i].index_gamma_tot>=0)   
                  {                                  
                    spec_activity(surrogate[i].index_gamma_tot)=
                      (1.0-alpha)*surrogate[i].species_activity(itemp)+alpha*surrogate[i].species_activity(itemp+1);
                  }

              tab_unifac(config.nmol_tot,config.nfunc_tot,config.groups_tot,
                         X_unifac,config.Inter_tot,config.InterB_tot,config.InterC_tot,config.RG_tot,
                         config.QG_tot,config.Rparam_tot,config.Qparam_tot,config.Lparam_tot,         
                         spec_activity,config.Z,Temperature,gamma_unifac,
			 config.temperature_dependancy);

              for (i=0;i<n;++i)
                {
                  if (surrogate[i].index_gamma_tot >=0)
                    surrogate[i].gamma_org=gamma_unifac(surrogate[i].index_gamma_tot);
                  else
                    surrogate[i].gamma_org=1.0;
                }
            }
        }	  
      else
        {
          if (all_hydrophobic==false)
            {
              unifac(config.nmol_org,config.nfunc_org,config.groups_org,
                     X_unifac,config.Inter_org,config.InterB_org,config.InterC_org,config.RG_org,
                     config.QG_org,config.Rparam_org,config.Qparam_org,config.Lparam_org,                 
                     config.Z,Temperature,gamma_unifac,config.temperature_dependancy);

              for (i=0;i<n;++i)
                {
                  if (surrogate[i].index_gamma_org >=0)
                    surrogate[i].gamma_org=gamma_unifac(surrogate[i].index_gamma_org);
                  else
                    surrogate[i].gamma_org=1.0;
                }
            }
          else
            {
              unifac(config.nmol_tot,config.nfunc_tot,config.groups_tot,
                     X_unifac,config.Inter_tot,config.InterB_tot,config.InterC_tot,config.RG_tot,
                     config.QG_tot,config.Rparam_tot,config.Qparam_tot,config.Lparam_tot,         
                     config.Z,Temperature,gamma_unifac,config.temperature_dependancy);

              for (i=0;i<n;++i)
                {
                  if (surrogate[i].index_gamma_tot >=0)
                    surrogate[i].gamma_org=gamma_unifac(surrogate[i].index_gamma_tot);
                  else
                    surrogate[i].gamma_org=1.0;
                }
            }
        }
    }
  else
    for (i=0;i<n;++i)
      surrogate[i].gamma_org=1.0;
}


void activity_coefficients_dyn_org(model_config &config, vector<species>& surrogate,
                                   double &Temperature, Array<double, 3> &MOW)
{
  //compute the activity coefficients with UNIFAC (short range interactions) for the organic phase
  //MOW: mean molar mass of the organic phase
  //all_hydrophobic: do all the compounds condense on the organic?(If low mass of particulate water)
  int n_unifac,i;
  int n=surrogate.size();
  double sum,sumX_unifac;
  Array<double, 1> X_unifac,gamma_unifac;

  int b,ilayer,iphase;
  
  if (config.use_global_dynamic_parameters)
    {
      for (iphase=0;iphase<config.max_number_of_phases;++iphase)
        {
          double MOWtemp=0.0;
          for (i=0;i<n;++i)
            {
              surrogate[i].Ap=0.0;
              for (ilayer=0;ilayer<config.nlayer;++ilayer)
                for (b=0;b<config.nbins;++b)
                  surrogate[i].Ap+=surrogate[i].Ap_layer_init(b,ilayer,iphase);
            }
          activity_coefficients_org(config, surrogate, false, Temperature, MOWtemp);
		  
          for (i=0;i<n;++i)
            for (ilayer=0;ilayer<config.nlayer;++ilayer)
              for (b=0;b<config.nbins;++b)
                surrogate[i].gamma_org_layer(b,ilayer,iphase)=surrogate[i].gamma_org;
		  
          for (ilayer=0;ilayer<config.nlayer;++ilayer)
            for (b=0;b<config.nbins;++b)
              MOW(b,ilayer,iphase)=MOWtemp;
        }
    }
  else
    {
      for (ilayer=0;ilayer<config.nlayer;++ilayer)
        for (b=0;b<config.nbins;++b)
          for (iphase=0;iphase<config.nphase(b,ilayer);++iphase)
            {
              for (i=0;i<n;++i)
                surrogate[i].Ap=surrogate[i].Ap_layer_init(b,ilayer,iphase);
			  
              activity_coefficients_org(config, surrogate, false, Temperature,
                                        MOW(b,ilayer,iphase));
              for (i=0;i<n;++i)
                surrogate[i].gamma_org_layer(b,ilayer,iphase)=surrogate[i].gamma_org;
            }
    }
}


void activity_coefficients_dyn_sat(model_config &config, vector<species>& surrogate,
                                   double &Temperature, Array<double, 3> &MOW, int &b2, int &ilayer2)
{
  //compute the activity coefficients with UNIFAC (short range interactions) for the organic phase
  //MOW: mean molar mass of the organic phase
  //all_hydrophobic: do all the compounds condense on the organic?(If low mass of particulate water)
  int n_unifac,i;
  int n=surrogate.size();
  double sum,sumX_unifac;
  Array<double, 1> X_unifac,gamma_unifac;

  int b,ilayer,iphase;
  
  if (config.use_global_dynamic_parameters)
    {
      for (iphase=0;iphase<config.max_number_of_phases;++iphase)
        {
          double MOWtemp=0.0;
          for (i=0;i<n;++i)
            {
              surrogate[i].Ap=0.0;
              for (ilayer=0;ilayer<config.nlayer;++ilayer)
                for (b=0;b<config.nbins;++b)
                  surrogate[i].Ap+=surrogate[i].Ap_layer_init(b,ilayer,iphase);
            }
          activity_coefficients_org(config, surrogate, false, Temperature, MOWtemp);
		  
          for (i=0;i<n;++i)
            for (ilayer=0;ilayer<config.nlayer;++ilayer)
              for (b=0;b<config.nbins;++b)
                surrogate[i].gamma_org_layer(b,ilayer,iphase)=surrogate[i].gamma_org;
		  
          for (ilayer=0;ilayer<config.nlayer;++ilayer)
            for (b=0;b<config.nbins;++b)
              MOW(b,ilayer,iphase)=MOWtemp;
        }
    }
  else
    {
      for (iphase=0;iphase<config.nphase(b2,ilayer2);++iphase)
        {
          for (i=0;i<n;++i)
            surrogate[i].Ap=surrogate[i].Ap_layer_init(b2,ilayer2,iphase);
		  
          activity_coefficients_org(config, surrogate, false, Temperature,
                                    MOW(b2,ilayer2,iphase));
          for (i=0;i<n;++i)
            surrogate[i].gamma_org_layer(b2,ilayer2,iphase)=surrogate[i].gamma_org;
        }
    }
}

void compute_ionic_strenght(model_config &config, vector<species>& surrogate,
                            double &AQinit, double &conc_inorganic,
                            double &ionic, double &chp,
                            double &organion, double &ionic_organic, double factor)
{
  //computes the ionic strenth and recomputes pH
  int i;
  int n=surrogate.size();
  double conc_org=0.0;
  conc_org=max(AQinit-conc_inorganic,config.MOmin);

  if (config.compute_aqueous_phase_properties)
    {
      double Ke=1.0e-14; //Dissociation constant of water Ke=a(H+)*a(HO-)
      conc_inorganic=0.0;
      ionic=ionic_organic;      
      double chp_old=chp;
      
      //modify the pH of the aqueous phase to take into account organion
      double inorganion=0.0;
      double conc_org=0.0;
      double sumxt=0.0;
      for (i=0;i<n;++i)
        if (surrogate[i].is_organic==false and i!=config.iH2O and surrogate[i].is_inorganic_precursor==false)
          {
            surrogate[i].molality=surrogate[i].Aaq/surrogate[i].MM/conc_org*1000.0;           
            if (surrogate[i].name!="H")
              inorganion-=surrogate[i].molality*surrogate[i].charge;
          }
      
      chp=factor*0.5*(organion+inorganion+pow(pow(organion+inorganion,2)+4*Ke,0.5))+(1.0-factor)*chp;
      if (chp==0.0)
        chp=pow(10.0,-5.6);

      //compute the ionic strength
      for (i=0;i<n;++i)
        if (surrogate[i].is_organic==false and i!=config.iH2O and surrogate[i].is_inorganic_precursor==false)
          {
            ionic=ionic+0.5*surrogate[i].molality*pow(surrogate[i].charge,2);
            conc_inorganic+=surrogate[i].Aaq;
          }      
    }
  else
    {
      for (i=0;i<n;++i)
        if (surrogate[i].is_organic==false and i!=config.iH2O and surrogate[i].is_inorganic_precursor==false)
          {
            surrogate[i].molality=surrogate[i].Aaq/surrogate[i].MM/conc_org*1000.0;
            if (surrogate[i].name=="H")
              surrogate[i].molality=chp;
          }
    }
  ionic=min(ionic,1000.);
}

void compute_ionic_strenght2(model_config &config, vector<species>& surrogate,
                             double &AQinit, double &conc_inorganic,
                             double &ionic, double &chp,
                             double &organion, double &ionic_organic, double factor)
{
  //same as compute_ionic_strenth but does not recompute the pH
  int i;
  int n=surrogate.size();
  double conc_org=0.0;
  
  conc_org=max(AQinit-conc_inorganic,config.MOmin);
     
  if (config.compute_aqueous_phase_properties)
    {
      double Ke=1.0e-14; //Dissociation constant of water Ke=a(H+)*a(HO-)
      conc_inorganic=0.0;
      ionic=ionic_organic;      
      
      //modify the pH of the aqueous phase to take into account organion
      for (i=0;i<n;++i)
        if (surrogate[i].is_organic==false and i!=config.iH2O and surrogate[i].is_inorganic_precursor==false)          
	  surrogate[i].molality=surrogate[i].Aaq/surrogate[i].MM/conc_org*1000.0;                   

      //compute the ionic strength
      for (i=0;i<n;++i)
        if (surrogate[i].is_inorganic_precursor==false and surrogate[i].is_organic==false and i!=config.iH2O)
          {
            if (surrogate[i].name=="H")
              {                                
                surrogate[i].Aaq=chp/1000.*surrogate[i].MM*conc_org;
		surrogate[i].molality=chp; 
              }
            ionic=ionic+0.5*surrogate[i].molality*pow(surrogate[i].charge,2);
            conc_inorganic+=surrogate[i].Aaq;
	    
          }     	  
    }
  else
    {      
      for (i=0;i<n;++i)
        if (surrogate[i].is_organic==false and i!=config.iH2O and surrogate[i].is_inorganic_precursor==false)
          {
            surrogate[i].molality=surrogate[i].Aaq/surrogate[i].MM/conc_org*1000.0;
          }
    }
  ionic=min(ionic,1000.);
}


void activity_coefficients_saturation(model_config &config, vector<species>& surrogate,
                                      bool all_hydrophobic,double &Temperature,
                                      Array <double, 1> &MOW)
{
  //compute the activity coefficients with UNIFAC (short range interactions) for the organic phase
  //MOW: mean molar mass of the organic phase
  //all_hydrophobic: do all the compounds condense on the organic?(If low mass of particulate water)
  int n_unifac,i,j;
  int n=surrogate.size();
  int nphase=MOW.size();
  double sum,sumX_unifac;
  Array<double, 1> X_unifac,gamma_unifac;

  //initialization of X_unifac
  if (all_hydrophobic==false)
    {
      X_unifac.resize(config.nmol_org);
      gamma_unifac.resize(config.nmol_org);
      n_unifac=config.nmol_org;
    }
  else
    {
      X_unifac.resize(config.nmol_tot);
      gamma_unifac.resize(config.nmol_tot);
      n_unifac=config.nmol_tot;
    }

  for (j=0;j<nphase;++j)
    {
      sum=0.0;
      sumX_unifac=0.0;
      for (i=0;i<n_unifac;++i)
        X_unifac(i)=0.0;
	   
      //computation of X_unifac (molar fraction of compounds used in the computation of activity
      // coefficients)
      for (i=0;i<n;++i)
        if (surrogate[i].hydrophobic or
            (all_hydrophobic and (surrogate[i].is_organic or i==config.iH2O)))
          {
            surrogate[i].Xorg_sat(j)=surrogate[i].Ap_sat(j)/surrogate[i].MM;
            sum+=surrogate[i].Xorg_sat(j);
            if (all_hydrophobic==false)
              {
                if (surrogate[i].index_gamma_org>=0)
                  {
                    sumX_unifac+=surrogate[i].Xorg_sat(j);
                    X_unifac(surrogate[i].index_gamma_org)+=surrogate[i].Xorg_sat(j);
                  }
              }
            else
              if (surrogate[i].index_gamma_tot>=0)
                {
                  sumX_unifac+=surrogate[i].Xorg_sat(j);
                  X_unifac(surrogate[i].index_gamma_tot)+=surrogate[i].Xorg_sat(j);
                }
          }
  
      if (sum>0.0)
        {
          MOW(j)=0.0;
          for (i=0;i<n;++i)
            if (surrogate[i].hydrophobic or
                (all_hydrophobic and (surrogate[i].is_organic or i==config.iH2O)))
              {
                surrogate[i].Xorg_sat(j)=surrogate[i].Xorg_sat(j)/sum;
                MOW(j)+=surrogate[i].Xorg_sat(j)*surrogate[i].MM;
              }
        }
      else
        MOW(j)=200.0;

      //Call of unifac
      if (sumX_unifac>0.0 and config.activity_model=="unifac")
        {
          for (i=0;i<n_unifac;++i)
            {
              X_unifac(i)/=sumX_unifac;
              X_unifac(i)=max(1.0e-11, X_unifac(i));
            }
	  
          if (config.tabulation_unifac)
            {
              if (all_hydrophobic==false)
                {
                  Array <double, 1> spec_activity;
                  spec_activity.resize(config.nmol_org);
                  int itemp=int(Temperature-263.15);
                  itemp=max(itemp,0);
                  itemp=min(itemp,config.ntemp);
                  double alpha;
                  if (Temperature<=263.15) 
                    alpha=0.0;            
                  else if (Temperature>=263.15+1.0*config.ntemp)
                    alpha=1.0;
                  else
                    alpha=Temperature-263.15-1.0*itemp;                  

                  for (i=0;i<n;i++)
                    if (surrogate[i].index_gamma_org>=0)                                     
                      spec_activity(surrogate[i].index_gamma_org)=
                        (1.0-alpha)*surrogate[i].species_activity(itemp)+alpha*surrogate[i].species_activity(itemp+1);

                  tab_unifac(config.nmol_org,config.nfunc_org,config.groups_org,
                             X_unifac,config.Inter_org,config.InterB_org,config.InterC_org,config.RG_org,
                             config.QG_org,config.Rparam_org,config.Qparam_org,config.Lparam_org,
                             spec_activity,config.Z,Temperature,gamma_unifac,config.temperature_dependancy);

                  for (i=0;i<n;++i)
                    {
                      if (surrogate[i].index_gamma_org >=0)
                        surrogate[i].gamma_org_sat(j)=gamma_unifac(surrogate[i].index_gamma_org);
                      else
                        surrogate[i].gamma_org_sat(j)=1.0;
                    }
                }
              else
                {
                  Array <double, 1> spec_activity;
                  spec_activity.resize(config.nmol_tot);
                  int itemp=int(Temperature-263.15);
                  itemp=max(itemp,0);
                  itemp=min(itemp,config.ntemp);
                  double alpha;
                  if (Temperature<=263.15) 
                    alpha=0.0;            
                  else if (Temperature>=263.15+1.0*config.ntemp)
                    alpha=1.0;
                  else
                    alpha=Temperature-263.15-1.0*itemp;                  

                  for (i=0;i<n;i++)
                    if (surrogate[i].index_gamma_tot>=0)                                     
                      spec_activity(surrogate[i].index_gamma_tot)=
                        (1.0-alpha)*surrogate[i].species_activity(itemp)+alpha*surrogate[i].species_activity(itemp+1);

                  tab_unifac(config.nmol_tot,config.nfunc_tot,config.groups_tot,
                             X_unifac,config.Inter_tot,config.InterB_tot,config.InterC_tot,config.RG_tot,
                             config.QG_tot,config.Rparam_tot,config.Qparam_tot,config.Lparam_tot,         
                             spec_activity,config.Z,Temperature,gamma_unifac,config.temperature_dependancy);

                  for (i=0;i<n;++i)
                    {
                      if (surrogate[i].index_gamma_tot >=0)
                        surrogate[i].gamma_org_sat(j)=gamma_unifac(surrogate[i].index_gamma_tot);
                      else
                        surrogate[i].gamma_org_sat(j)=1.0;
                    }

                }
            }
          else
            {
              if (all_hydrophobic==false)
                {
                  unifac(config.nmol_org,config.nfunc_org,config.groups_org,
                         X_unifac,config.Inter_org,config.InterB_org,config.InterC_org,config.RG_org,
                         config.QG_org,config.Rparam_org,config.Qparam_org,config.Lparam_org,
                         config.Z,Temperature,gamma_unifac,config.temperature_dependancy);

                  for (i=0;i<n;++i)
                    {
                      if (surrogate[i].index_gamma_org >=0)
                        surrogate[i].gamma_org_sat(j)=gamma_unifac(surrogate[i].index_gamma_org);
                      else
                        surrogate[i].gamma_org_sat(j)=1.0;
                    }
                }
              else
                {
                  unifac(config.nmol_tot,config.nfunc_tot,config.groups_tot,
                         X_unifac,config.Inter_tot,config.InterB_tot,config.InterC_tot,config.RG_tot,
                         config.QG_tot,config.Rparam_tot,config.Qparam_tot,config.Lparam_tot,
                         config.Z,Temperature,gamma_unifac,config.temperature_dependancy);
              
                  for (i=0;i<n;++i)
                    {
                      if (surrogate[i].index_gamma_tot >=0)
                        surrogate[i].gamma_org_sat(j)=gamma_unifac(surrogate[i].index_gamma_tot);
                      else
                        surrogate[i].gamma_org_sat(j)=1.0;
                    }
                }
            }
        }
      else
        for (i=0;i<n;++i)
          surrogate[i].gamma_org_sat(j)=1.0;
    }
}

void hygroscopicity_coupled(model_config &config, vector<species>& surrogate, 
			    double &Temperature, double &MOW, double &MMaq,
			    double &RH, double &LWC, double &conc_inorganic,
			    double &MOinit, double &MO, 
			    double &AQinit, double &AQ,
			    double &deriv_error1_MO, 
			    double &deriv_error1_AQ,
			    double &deriv_error2_MO,
			    double &deriv_error2_AQ,
			    double factor)
{
  //Compute the absorption of water by on organic phase with a concentration of MOinit
  //MOW: mean molar mass
  //RH: relative humidity
  //surrogate[config.iH2O].Ap : mass of water in the organic phase
  //at equilibrium according to Raoult's law:
  // Xwater = RH/gamma_water
  // with Xwater = moles water/(sum of moles) = Ap/MMwater * MOW/MOinit  

  if (config.activity_model=="unifac" or config.compute_inorganic)
    {     
      double Pwater=exp(13.7-5120.0/Temperature)*1.013e5;
      double alpha_org;
      if (surrogate[config.iH2O].hydrophobic)
	alpha_org=surrogate[config.iH2O].MM/MOW/surrogate[config.iH2O].gamma_org;
      else 
	alpha_org=0.0;
      double beta_org=alpha_org*8.314*Temperature/(surrogate[config.iH2O].MM*Pwater)*1.0e-6;
      double alpha_aq,beta_aq;
      if (surrogate[config.iH2O].hydrophilic)
	{
	  alpha_aq=surrogate[config.iH2O].MM/MMaq/surrogate[config.iH2O].gamma_aq;
	  beta_aq=alpha_aq*8.314*Temperature/(surrogate[config.iH2O].MM*Pwater)*1.0e-6; 
	}
      else
	{
	  alpha_aq=0.0;
	  beta_aq=0.0;
	}
      
      surrogate[config.iH2O].Ap=factor*alpha_org*RH*MOinit/(1.0+beta_org*MOinit+beta_aq*AQinit)+(1.0-factor)*surrogate[config.iH2O].Ap;
      if (config.compute_inorganic)
	surrogate[config.iH2O].Aaq=factor*alpha_aq*RH*AQinit/(1.0+beta_org*MOinit+beta_aq*AQinit)+(1.0-factor)*surrogate[config.iH2O].Aaq;
      else
	surrogate[config.iH2O].Aaq=factor*max(alpha_aq*RH*AQinit/(1.0+beta_org*MOinit+beta_aq*AQinit)-LWC,0.0)+(1.0-factor)*surrogate[config.iH2O].Aaq;
      MO+=surrogate[config.iH2O].Ap;
      AQ+=surrogate[config.iH2O].Aaq;
      deriv_error1_MO-=alpha_org*RH/(1.0+beta_org*MOinit+beta_aq*AQinit)-alpha_org*RH*beta_org*MOinit/pow(1.0+beta_org*MOinit+beta_aq*AQinit,2.0);
      deriv_error1_AQ+=alpha_org*RH*beta_aq*MOinit/pow(1.0+beta_org*MOinit+beta_aq*AQinit,2.0);

      if (surrogate[config.iH2O].Aaq>0.0)
	{
	  deriv_error2_AQ-=alpha_aq*RH/(1.0+beta_aq*AQinit+beta_org*MOinit)-alpha_aq*RH*beta_aq*AQinit/pow(1.0+beta_aq*AQinit+beta_org*MOinit,2.0);
	  deriv_error2_MO+=alpha_aq*RH*beta_org*AQinit/pow(1.0+beta_aq*AQinit+beta_org*MOinit,2.0);
	}

    }
  else
    {
      double Mads=AQinit-LWC-conc_inorganic;
      double Pwater=exp(13.7-5120.0/Temperature)*1.013e5;
      double alpha_org=0.0;
      if (surrogate[config.iH2O].hydrophobic)
	alpha_org=surrogate[config.iH2O].MM/MOW;      
      else
	alpha_org=0.0;
      double beta_org=alpha_org*8.314*Temperature/(surrogate[config.iH2O].MM*Pwater)*1.0e-6;
      double alpha_aq,beta_aq;
      if (surrogate[config.iH2O].hydrophilic)
	{
	  alpha_aq=surrogate[config.iH2O].MM/MMaq;
	  beta_aq=alpha_aq*8.314*Temperature/(surrogate[config.iH2O].MM*Pwater)*1.0e-6; 
	}      
      else
	{
	  alpha_aq=0.0;
	  beta_aq=0.0;
	}          
      double fdispo=1.0-LWC*beta_aq/(alpha_aq*RH);
 
      surrogate[config.iH2O].Ap=factor*alpha_org*RH*MOinit/(1.0+beta_org*MOinit+beta_aq*Mads)*fdispo+(1.0-factor)*surrogate[config.iH2O].Ap;
      surrogate[config.iH2O].Aaq=factor*alpha_aq*RH*Mads/(1.0+beta_org*MOinit+beta_aq*Mads)*fdispo+(1.0-factor)*surrogate[config.iH2O].Aaq;
      MO+=surrogate[config.iH2O].Ap;
      AQ+=surrogate[config.iH2O].Aaq;

      deriv_error1_MO-=(alpha_org*RH/(1.0+beta_org*MOinit+beta_aq*Mads)-alpha_org*RH*beta_org*MOinit/pow(1.0+beta_org*MOinit+beta_aq*Mads,2.0))*fdispo;
      deriv_error1_AQ+=(alpha_org*RH*beta_aq*MOinit/pow(1.0+beta_org*MOinit+beta_aq*Mads,2.0))*fdispo;

      deriv_error2_AQ-=(alpha_aq*RH/(1.0+beta_aq*Mads+beta_org*MOinit)-alpha_aq*RH*beta_aq*Mads/pow(1.0+beta_aq*Mads+beta_org*MOinit,2.0))*fdispo;
      deriv_error2_MO+=(alpha_aq*RH*beta_org*Mads/pow(1.0+beta_aq*Mads+beta_org*MOinit,2.0))*fdispo;
      
    }

}

void hygroscopicity_org(model_config &config, vector<species>& surrogate, double &Temperature, double &MOW,
                        double &RH, double &MOinit, double &MO, double &derivative, double factor)
{
  //Compute the absorption of water by on organic phase with a concentration of MOinit
  //MOW: mean molar mass
  //RH: relative humidity
  //surrogate[config.iH2O].Ap : mass of water in the organic phase
  //at equilibrium according to Raoult's law:
  // Xwater = RH/gamma_water
  // with Xwater = moles water/(sum of moles) = Ap/MMwater * MOW/MOinit  

  double alpha=surrogate[config.iH2O].MM/MOW/surrogate[config.iH2O].gamma_org;
  double Pwater=exp(13.7-5120.0/Temperature)*1.013e5;
  double beta=alpha*8.314*Temperature/(surrogate[config.iH2O].MM*Pwater)*1.0e-6;

  surrogate[config.iH2O].Ap=factor*alpha*RH*MOinit/(1.0+beta*MOinit)+(1.0-factor)*surrogate[config.iH2O].Ap;
  MO+=surrogate[config.iH2O].Ap;
  derivative-=alpha*RH/(1.0+beta*MOinit)-alpha*RH*beta*MOinit/pow(1.0+beta*MOinit,2.0);
}

void hygroscopicity_org_sat(model_config &config, vector<species> &surrogate, double &Temperature,
                            Array <double, 1> &MOW, double &MMaq,
                            double &RH, double LWC, double conc_inorganic,
			    Array <double, 1> &MOinit, Array <double,1> &MO, 
			    double &AQinit, double &AQ, 
			    Array <double, 2> &Jacobian,			    
			    bool compute_aq, double factor)
{
  //Compute the absorption of water by on organic phase with a concentration of MOinit
  //MOW: mean molar mass
  //RH: relative humidity
  //surrogate[config.iH2O].Ap : mass of water in the organic phase
  //at equilibrium according to Raoult's law:
  // Xwater = RH/gamma_water
  // with Xwater = moles water/(sum of moles) = Ap/MMwater * MOW/MOinit

  int nphase=MOinit.size();
  int iphase,jphase;

  double Pwater=exp(13.7-5120.0/Temperature)*1.013e5;
  Array <double,1> alpha_org,beta_org;
  alpha_org.resize(nphase);
  beta_org.resize(nphase);
  double alpha_aq=0.0;
  double beta_aq=0.0;
  double Mads=AQinit-LWC-conc_inorganic;
  if (compute_aq)
    if (config.activity_model=="unifac")
      {
	alpha_aq=surrogate[config.iH2O].MM/MMaq/surrogate[config.iH2O].gamma_aq;
	beta_aq=alpha_aq*8.314*Temperature/(surrogate[config.iH2O].MM*Pwater)*1.0e-6;
      }
    else
      {	
	alpha_aq=surrogate[config.iH2O].MM/MMaq;
	beta_aq=alpha_aq*8.314*Temperature/(surrogate[config.iH2O].MM*Pwater)*1.0e-6; 
      }

  double fdispo=1.0;
  if (config.activity_model!="unifac" and compute_aq)
    fdispo=(1.0-LWC*beta_aq/(alpha_aq*RH));
  double sum=1.0;
  if (config.activity_model=="unifac")
    sum+=beta_aq*AQinit;
  else
    sum+=beta_aq*Mads;
  for (iphase=0;iphase<nphase;iphase++)
    {
      if (surrogate[config.iH2O].hydrophobic)
	alpha_org(iphase)=surrogate[config.iH2O].MM/MOW(iphase)/surrogate[config.iH2O].gamma_org_sat(iphase);
      else 
	alpha_org(iphase)=0.0;
      beta_org(iphase)=alpha_org(iphase)*8.314*Temperature/(surrogate[config.iH2O].MM*Pwater)*1.0e-6;
      sum+=beta_org(iphase)*MOinit(iphase);
    }
     
  if (compute_aq)
    {
      if (config.activity_model=="unifac")
	{
	  if (config.compute_inorganic)
	    surrogate[config.iH2O].Aaq=factor*max(alpha_aq*RH*AQinit/sum-LWC,0.0)+(1.0-factor)*surrogate[config.iH2O].Aaq;      
	  else
	    surrogate[config.iH2O].Aaq=factor*alpha_aq*RH*AQinit/sum+(1.0-factor)*surrogate[config.iH2O].Aaq;  

	  AQ+=surrogate[config.iH2O].Aaq;    
	  if (surrogate[config.iH2O].Aaq>0.0)
	    {
	      Jacobian(nphase,nphase)-=alpha_aq*RH/sum-alpha_aq*RH*beta_aq*AQinit/pow(sum,2.0);
	      for (iphase=0;iphase<nphase;iphase++)
		Jacobian(nphase,iphase)+=alpha_aq*RH*beta_org(iphase)*AQinit/pow(sum,2.0);
	    }
	}
      else
	{
	  surrogate[config.iH2O].Aaq=factor*alpha_aq*RH*Mads/sum*fdispo+(1.0-factor)*surrogate[config.iH2O].Aaq;      
	  AQ+=surrogate[config.iH2O].Aaq;    	  
	  Jacobian(nphase,nphase)-=(alpha_aq*RH/sum-alpha_aq*RH*beta_aq*Mads/pow(sum,2.0))*fdispo;
	  for (iphase=0;iphase<nphase;iphase++)
	    Jacobian(nphase,iphase)+=alpha_aq*RH*beta_org(iphase)*Mads/pow(sum,2.0)*fdispo;	    
	}
    }

  for (iphase=0;iphase<nphase;iphase++)
    {      
      surrogate[config.iH2O].Ap_sat(iphase)=factor*alpha_org(iphase)*RH*MOinit(iphase)/sum*fdispo
	+(1.0-factor)*surrogate[config.iH2O].Ap_sat(iphase);
      MO(iphase)+=surrogate[config.iH2O].Ap_sat(iphase);
      for (jphase=0;jphase<nphase;jphase++)
	{
	  if (iphase!=jphase)
	    Jacobian(iphase,jphase)+=alpha_org(iphase)*RH*beta_aq*MOinit(iphase)/pow(sum,2.0)*fdispo;
	  else
	    Jacobian(iphase,jphase)-=(alpha_org(iphase)*RH/sum-alpha_org(iphase)*RH*beta_org(iphase)*MOinit(iphase)/pow(sum,2.0))*fdispo;
	  
	  if (compute_aq)
	    if (config.activity_model=="unifac")
	      {
		if (surrogate[config.iH2O].Aaq>0.0)
		  Jacobian(iphase,nphase)+=alpha_org(iphase)*RH*beta_aq*AQinit/pow(sum,2.0)*fdispo;
	      }
	    else
	      Jacobian(iphase,nphase)+=alpha_org(iphase)*RH*beta_aq*Mads/pow(sum,2.0)*fdispo;
	}
    }




}

void activity_coefficients_aq(model_config config, vector<species>& surrogate,
                              double Temperature, double LWC,
                              double &MMaq, double &XH2O)
{
  //compute the activity coefficients of the aqueous phase
  //LWC: liquid water content due only to inorganic compounds
  //MMaq: mean molar mass of the aqueous phase
  Array<double, 1> X_unifac, gamma_unifac,Xions,gamma_ions;
  int n=surrogate.size();
  int i,iion;
  X_unifac.resize(config.nmol_aq);
  Xions.resize(config.nion_unifac);
  gamma_unifac.resize(config.nmol_aq);
  gamma_ions.resize(config.nion_unifac);
  double sumX_unifac=0.0;
  double sum=0.0;

  //initialization 
  X_unifac=0.0;
  Xions=0.0;

  MMaq=0.0;
  XH2O=0.0;
  
  //computation of X_unifac (molar fraction of compounds used in the computation of activity
  // coefficients) and MMaq 
  double conc_org=LWC;
  double sum_org=0.0;
  for (i=0;i<n;++i)
    if (surrogate[i].hydrophilic)
      {
        surrogate[i].Xaq=surrogate[i].Aaq/surrogate[i].MM;
        if (i==config.iH2O)
          surrogate[i].Xaq+=LWC/surrogate[i].MM;
        sum+=surrogate[i].Xaq;
        if (surrogate[i].index_gamma_aq>=0)
          {	   
            sumX_unifac+=surrogate[i].Xaq;
            X_unifac(surrogate[i].index_gamma_aq)+=surrogate[i].Xaq;
	    conc_org+=surrogate[i].Aaq;
	    sum_org+=surrogate[i].Xaq;
          }
      }

  conc_org=max(conc_org,config.MOmin);

  double xsms=0.0;
  double summolal=0.0;        
  if (sum_org>0.0)
    for (i=0;i<n;++i)     
      {
        if (surrogate[i].is_organic or i==config.iH2O)	          
          {
            if (surrogate[i].hydrophilic)
              xsms+=surrogate[i].Xaq/sum_org*surrogate[i].MM*1.0e-3;	
          }
        else if (surrogate[i].is_inorganic_precursor==false)
          summolal+=surrogate[i].Aaq/surrogate[i].MM/conc_org*1000.0;             
      } 

  
  double gamma_molal=1.0;
  if (xsms>0.0)
    gamma_molal=1.0/(0.018/xsms+0.018*summolal);   

  for (i=0;i<n;++i)
    if (surrogate[i].hydrophilic)
      {
        surrogate[i].Xaq/=sum;
        MMaq+=surrogate[i].Xaq*surrogate[i].MM;
      }

  if (config.iH2O>=0)
    XH2O=surrogate[config.iH2O].Xaq;

  //Call of unifac
  if (config.compute_organic and sumX_unifac>0.0 
      and config.activity_model=="unifac")
    {
      if (config.SR_ions)
        {
          iion=0;
          for (i=0;i<n;++i)
            if (surrogate[i].is_organic==false and i!=config.iH2O and surrogate[i].is_inorganic_precursor==false)
              {
                sumX_unifac+=surrogate[i].Xaq;
                Xions(iion)=surrogate[i].Xaq;
                iion++;
              }

          for (i=0;i<config.nion_unifac;i++)
            Xions(i)/=sumX_unifac;
        }

      for (i=0;i<config.nmol_aq;++i)
        {
          X_unifac(i)/=sumX_unifac;
          X_unifac(i)=max(1.0e-11, X_unifac(i));
        }

      if (config.SR_ions)
        {
          if (config.tabulation_unifac)
            {
              Array <double, 1> spec_activity;
              spec_activity.resize(config.nmol_aq);
              int itemp=int(Temperature-263.15);
              itemp=max(itemp,0);
              itemp=min(itemp,config.ntemp);
              double alpha;
              if (Temperature<=263.15) 
                alpha=0.0;            
              else if (Temperature>=263.15+1.0*config.ntemp)
                alpha=1.0;
              else
                alpha=Temperature-263.15-1.0*itemp;          

              for (i=0;i<n;i++)
                if (surrogate[i].index_gamma_aq>=0)                                     
                  spec_activity(surrogate[i].index_gamma_aq)=
                    (1.0-alpha)*surrogate[i].species_activity(itemp)+alpha*surrogate[i].species_activity(itemp+1);

              tab_unifac_aq(config.nmol_aq,config.nion_unifac,config.nfunc_aq,config.groups_aq,
                            X_unifac,Xions,config.Inter_aq,config.InterB_aq,config.InterC_aq,
			    config.RG_aq,
                            config.QG_aq,config.Rparam_aq,config.Qparam_aq,config.Lparam_aq,
                            config.RGions,config.QGions, config.Lions,
                            spec_activity,config.Z,Temperature,gamma_unifac,gamma_ions,
			    config.temperature_dependancy);
            }
          else
            unifac_aq(config.nmol_aq,config.nion_unifac,config.nfunc_aq,config.groups_aq,
                      X_unifac,Xions,config.Inter_aq,config.InterB_aq,config.InterC_aq,config.RG_aq,
                      config.QG_aq,config.Rparam_aq,config.Qparam_aq,config.Lparam_aq,
                      config.RGions,config.QGions, config.Lions,
                      config.Z,Temperature,gamma_unifac,gamma_ions,config.temperature_dependancy);
        
	  iion=0;
          for (i=0;i<n;++i)
            {
              if (surrogate[i].is_organic==false and i!=config.iH2O and surrogate[i].is_inorganic_precursor==false)
                {
                  surrogate[i].gamma_aq=gamma_ions(iion)*gamma_molal;
		    
		  if (i==config.iHp)
		    surrogate[i].gamma_aq=gamma_ions(iion)*gamma_molal;
                  iion++;
                }
              else
                {
                  if (surrogate[i].index_gamma_aq >=0)
                    surrogate[i].gamma_aq=gamma_unifac(surrogate[i].index_gamma_aq)
                      /surrogate[i].GAMMAinf;
                  else
                    surrogate[i].gamma_aq=1.0;
                }
            }
        }
      else
        {
          if (config.tabulation_unifac)
            {
              Array <double, 1> spec_activity;
              spec_activity.resize(config.nmol_aq);
              int itemp=int(Temperature-263.15);
              itemp=max(itemp,0);
              itemp=min(itemp,config.ntemp);
              double alpha;
              if (Temperature<=263.15) 
                alpha=0.0;            
              else if (Temperature>=263.15+1.0*config.ntemp)
                alpha=1.0;
              else
                alpha=Temperature-263.15-1.0*itemp;          

              for (i=0;i<n;i++)
                if (surrogate[i].index_gamma_aq>=0)                                     
                  spec_activity(surrogate[i].index_gamma_aq)=
                    (1.0-alpha)*surrogate[i].species_activity(itemp)+alpha*surrogate[i].species_activity(itemp+1);

              tab_unifac(config.nmol_aq,config.nfunc_aq,config.groups_aq,
                         X_unifac,config.Inter_aq,config.InterB_aq,config.InterC_aq,config.RG_aq,
                         config.QG_aq,config.Rparam_aq,config.Qparam_aq,config.Lparam_aq,                        
                         spec_activity,config.Z,Temperature,gamma_unifac,config.temperature_dependancy);
            }
          else
            unifac(config.nmol_aq,config.nfunc_aq,config.groups_aq,
                   X_unifac,config.Inter_aq,config.InterB_aq,config.InterC_aq,
		   config.RG_aq,
                   config.QG_aq,config.Rparam_aq,config.Qparam_aq,config.Lparam_aq,                  
                   config.Z,Temperature,gamma_unifac,config.temperature_dependancy);
        
          for (i=0;i<n;++i)
            {
              if (surrogate[i].index_gamma_aq >=0)
                surrogate[i].gamma_aq=gamma_unifac(surrogate[i].index_gamma_aq)
                  /surrogate[i].GAMMAinf;
              else
                surrogate[i].gamma_aq=1.0;
            }
        }

    }
  else
    for (i=0;i<n;++i)
      surrogate[i].gamma_aq=1.0;
}

void hygroscopicity_aq(model_config &config, vector<species>& surrogate,
                       double &Temperature, double &RH,
                       double &AQinit, double &LWC,
                       double &conc_inorganic, double &MMaq,
                       double &AQ, double &derivative, double factor)
{
  //LWC liquid water content due to inorganic coumpounds (µg/m3)
  //Compute the absorption of water by the organic compounds in the aqueous phase

  if (config.activity_model=="unifac")
    {
      double alpha=surrogate[config.iH2O].MM/MMaq/surrogate[config.iH2O].gamma_aq;
      double Pwater=exp(13.7-5120.0/Temperature)*1.013e5;
      double beta=alpha*8.314*Temperature/(surrogate[config.iH2O].MM*Pwater)*1.0e-6;     

      surrogate[config.iH2O].Aaq=factor*
        max(alpha*RH*AQinit/(1.0+beta*AQinit)-LWC,0.0)
        +(1.0-factor)*surrogate[config.iH2O].Aaq;

      AQ+=surrogate[config.iH2O].Aaq;

      if (surrogate[config.iH2O].Aaq>0.0)
	derivative-=alpha*RH/(1.0+beta*AQinit)-alpha*RH*beta*AQinit/pow(1.0+beta*AQinit,2.0);
    }
  else
    {
      double Mads=AQinit-LWC-conc_inorganic; //mass of the organic compounds in the aqueous phase
      double alpha=surrogate[config.iH2O].MM/MMaq/surrogate[config.iH2O].gamma_aq;
      double Pwater=exp(13.7-5120.0/Temperature)*1.013e5;
      double beta=alpha*8.314*Temperature/(surrogate[config.iH2O].MM*Pwater)*1.0e-6;

      surrogate[config.iH2O].Aaq=factor*alpha*RH*Mads/(1.0+beta*AQinit)
        +(1.0-factor)*surrogate[config.iH2O].Aaq;

      AQ+=surrogate[config.iH2O].Aaq;

      if (surrogate[config.iH2O].Aaq>0.0)
	derivative-=alpha*RH/(1.0+beta*AQinit)-alpha*RH*beta*Mads/pow(1.0+beta*AQinit,2.0);

    }
}

void hygroscopicity_tot(model_config &config, vector<species>& surrogate,
			double &Temperature, double &RH,
			double &AQinit, 
			double &conc_inorganic, double &MMaq,
			double &AQ, double &derivative, double factor)
{
  //LWC liquid water content due to inorganic coumpounds (µg/m3)
  //Compute the absorption of water by the inorganic and organic compounds in the aqueous phase
  /*surrogate[config.iH2O].Aaq=factor*surrogate[config.iH2O].MM/MMaq*AQinit*RH/surrogate[config.iH2O].gamma_aq
    +(1.0-factor)*surrogate[config.iH2O].Aaq; //mass of water absorbed by the organic and inorganic compounds
  //in the aqueous phase
  AQ+=surrogate[config.iH2O].Aaq;
  derivative-=surrogate[config.iH2O].MM/MMaq*RH/surrogate[config.iH2O].gamma_aq;*/

  double alpha=surrogate[config.iH2O].MM/MMaq/surrogate[config.iH2O].gamma_aq;
  double Pwater=exp(13.7-5120.0/Temperature)*1.013e5;
  double beta=alpha*8.314*Temperature/(surrogate[config.iH2O].MM*Pwater)*1.0e-6;

  surrogate[config.iH2O].Aaq=factor*alpha*RH*AQinit/(1.0+beta*AQinit)
    +(1.0-factor)*surrogate[config.iH2O].Aaq;

  AQ+=surrogate[config.iH2O].Aaq;
  derivative-=alpha*RH/(1.0+beta*AQinit)-alpha*RH*beta*AQinit/pow(1.0+beta*AQinit,2.0);

}

void activity_coefficients_LR_MR(model_config &config, vector<species>& surrogate,
                                 double &Temperature, double LWC, double &ionic)
{
  //routine to call aiomfac to compute activity coefficients due to long range and medium range
  // interactions
  // gamma = gamma_SR (short range computed previously with UNIFAC)
  //         * gamma_MR (medium range)
  //         * gamma_LR (long range)
  // gamma_MR = gamma_LR = 1.0 without ions (used therefore only for the aqueous phase)
  //ionic = ionic strength (mol/kg)
  
  double Xtot=0.0;
  int n=surrogate.size();
  int i,j;
  Array<double,1> molality,gamma_LR_ions,gamma_MR_ions, charges_ions,molar_mass_groups;
  Array<double,1> gamma_LR_solvents,gamma_MR_solvents,X_aiomfac,molar_mass_solvents;
  molality.resize(config.nion_aiomfac);
  gamma_LR_ions.resize(config.nion_aiomfac);
  gamma_MR_ions.resize(config.nion_aiomfac);
  gamma_LR_solvents.resize(config.nmol_aiomfac);
  gamma_MR_solvents.resize(config.nmol_aiomfac);
  X_aiomfac.resize(config.nmol_aiomfac);
  charges_ions.resize(config.nion_aiomfac);
  molar_mass_solvents.resize(config.nmol_aiomfac);
  molar_mass_groups.resize(config.ngroup_aiomfac);
 
  for (i=0;i<n;++i) //molality (mol/kg) and charge of inorganic ions
    if (surrogate[i].is_organic==false and i!=config.iH2O and surrogate[i].is_inorganic_precursor==false)
      {
        molality(surrogate[i].index_ion)=surrogate[i].molality;	
        charges_ions(surrogate[i].index_ion)=surrogate[i].charge;      
      }
    
  //compute the molar fraction of solvent species used in aiomfac
  // and their molar masses (in kg/mol)
  for (i=0;i<n;++i)
    if ((surrogate[i].is_organic or i==config.iH2O) and surrogate[i].hydrophilic
        and surrogate[i].index_gamma_aiomfac>=0)
      {
        surrogate[i].Xaq=surrogate[i].Aaq/surrogate[i].MM;
        if (i==config.iH2O)
          surrogate[i].Xaq+=LWC/surrogate[i].MM;
        Xtot+=surrogate[i].Xaq;
      }

  if (Xtot==0.0) 
    {
      surrogate[config.iH2O].Xaq=1.0;
      Xtot=1.0;
    }
 
  for (i=0;i<n;++i)
    if ((surrogate[i].is_organic or i==config.iH2O) and surrogate[i].hydrophilic
        and surrogate[i].index_gamma_aiomfac>=0)
      {
        surrogate[i].Xaq/=Xtot;
        X_aiomfac(surrogate[i].index_gamma_aiomfac)=surrogate[i].Xaq;
        molar_mass_solvents(surrogate[i].index_gamma_aiomfac)=surrogate[i].MM*0.001;
      }
  
  //Compute the mean molar masses of aiomfac groups in the aqueous phase according to:
  //   the number of aiomfac groups in the molecule (config.groups_aiomfac)
  //   and the molar mass (in kg/mol) of the group in the molecule (config.Mgroups_aiomfac)
  
  for (j=0;j<config.ngroup_aiomfac;++j)
    {
      Xtot=0.0;
      molar_mass_groups(j)=0.0;
      for (i=0;i<config.nmol_aiomfac;++i)
        Xtot+=config.groups_aiomfac(i,j)*X_aiomfac(i);

      if (Xtot>0.0)
        for (i=0;i<config.nmol_aiomfac;++i)
          molar_mass_groups(j)+=config.Mgroups_aiomfac(i,j)*
            config.groups_aiomfac(i,j)*X_aiomfac(i)/Xtot*0.001; 
      else
	{
	  for (i=0;i<config.nmol_aiomfac;++i)
	    if (config.groups_aiomfac(i,j)>0.0)
	      {
		molar_mass_groups(j)+=config.groups_aiomfac(i,j)*config.Mgroups_aiomfac(i,j)*0.001;
		Xtot+=config.groups_aiomfac(i,j);
	      }	  	  
	  molar_mass_groups(j)/=Xtot;
	}
    }

  //Call of aiomfac
  aiomfac(X_aiomfac, gamma_LR_solvents,gamma_MR_solvents,molar_mass_solvents,molar_mass_groups,
          molality,gamma_LR_ions,gamma_MR_ions,charges_ions,
          Temperature,ionic,config.ngroup_aiomfac,config.groups_aiomfac,
          config.b1ki_aq,config.b2ki_aq,config.b1ca_aq,
          config.b2ca_aq, config.b3ca_aq, config.c1ca_aq,
          config.c2ca_aq, config.Rcc_aq, config.Qcca_aq);
  
  for (i=0;i<n;++i)
    if (surrogate[i].index_gamma_aiomfac>=0)
      {
        surrogate[i].gamma_LR=gamma_LR_solvents(surrogate[i].index_gamma_aiomfac);
        surrogate[i].gamma_SRMR=surrogate[i].gamma_aq*gamma_MR_solvents(surrogate[i].index_gamma_aiomfac);       
        surrogate[i].gamma_aq*=gamma_LR_solvents(surrogate[i].index_gamma_aiomfac)*
          gamma_MR_solvents(surrogate[i].index_gamma_aiomfac);
      }
  
  for (i=0;i<n;++i)
    if (surrogate[i].is_organic==false and i!=config.iH2O and surrogate[i].is_inorganic_precursor==false)
      {
        surrogate[i].gamma_LR=gamma_LR_ions(surrogate[i].index_ion);
        surrogate[i].gamma_SRMR=surrogate[i].gamma_aq*gamma_MR_ions(surrogate[i].index_ion);        
        surrogate[i].gamma_aq=surrogate[i].gamma_LR*surrogate[i].gamma_SRMR;
      }
}

void compute_organion(model_config &config, vector<species>& surrogate,
                      double &Temperature,  Array <double, 1> &MMaq,
                      Array <double, 1> &AQinit, Array<double, 1> &ionic,
                      Array<double, 1> &chp, Array <double, 1> &ionic_organic,
                      Array <double, 1> &organion, Array <double, 1> &LWC)
{
  double fion1,fion2;
  double molality1,molality2;
  int i,b;
  int n=surrogate.size();
  double Kp;
  
  for (b=0;b<config.nbins;++b)
    {
      double conc_org=LWC(b);
      for (i=0;i<n;i++)
	if (surrogate[i].hydrophilic)
	  conc_org+=surrogate[i].Aaq_bins_init(b);
      conc_org=max(conc_org,config.MOmin);

      ionic_organic(b)=0.0;
      organion(b)=0.0;      
      for (i=0;i<n;++i)
	if (surrogate[i].hydrophilic and surrogate[i].is_organic)
	  if(surrogate[i].nonvolatile==false)
	    {
	      Kp=surrogate[i].Kp_eff_aq(config, Temperature, ionic(b), chp(b),surrogate[config.iHp].LR(b),
					surrogate[config.iHp].SRMR(b),MMaq(b), fion1, fion2);
	      //molality1: molality of ions HA- or A-  
	      molality1=surrogate[i].Aaq_bins_init(b)*fion1/surrogate[i].MM/conc_org*1000.0;
	      //molality2: molality of ions A2-
	      molality2=surrogate[i].Aaq_bins_init(b)*fion2/surrogate[i].MM/conc_org*1000.0;
	      //compute ionic_organic and organion
	      ionic_organic(b)+=0.5*molality1+0.5*molality2*4;
	      organion(b)+=molality1+2*molality2;
	    }
	  
    }
  
}

void compute_organion2(model_config &config, vector<species>& surrogate,
		       double &Temperature,  Array <double, 1> &MMaq,
		       Array <double, 1> &AQ, Array<double, 1> &ionic,
		       Array<double, 1> &chp, Array <double, 1> &ionic_organic,
		       Array <double, 1> &organion, Array <double, 1> &LWC)
{
  double fion1,fion2;
  double molality1,molality2;
  int i,b;
  int n=surrogate.size();
  double Kp;
  
  for (b=0;b<config.nbins;++b)
    {
      ionic_organic(b)=0.0;
      organion(b)=0.0;

      double conc_org=LWC(b);
      for (i=0;i<n;i++)
	if (surrogate[i].hydrophilic)
	  conc_org+=surrogate[i].Aaq_bins_init(b);
      conc_org=max(conc_org,config.MOmin);
      
      for (i=0;i<n;++i)
	if (surrogate[i].hydrophilic and surrogate[i].is_organic)
	  if(surrogate[i].nonvolatile==false)
	    {
	      Kp=surrogate[i].Kp_eff_aq(config, Temperature, ionic(b), chp(b),surrogate[config.iHp].LR(b),
					surrogate[config.iHp].SRMR(b),MMaq(b), fion1, fion2);
	      //molality1: molality of ions HA- or A-  
	      molality1=surrogate[i].Aaq_bins(b)*fion1/surrogate[i].MM/conc_org*1000.0;
	      //molality2: molality of ions A2-
	      molality2=surrogate[i].Aaq_bins(b)*fion2/surrogate[i].MM/conc_org*1000.0;
	      //compute ionic_organic and organion
	      ionic_organic(b)+=0.5*molality1+0.5*molality2*4;
	      organion(b)+=molality1+2*molality2;
	    }	 
    }
  
}

double species::Kequilibrium(double &Temperature)
{
  double K0,deltaH_over_RT0,deltaCp0_over_R,value;
  double T0=298.15;
  K0=0.0;
  deltaH_over_RT0=0.0;
  deltaCp0_over_R=0.0;
  if (name=="H2SO4")
    {
      K0=1.015e-2;  //mol/kg
      deltaH_over_RT0=8.85;
      deltaCp0_over_R=25.14;
    }
  else if (name=="NH3")
    {
      double Keau=1.010e-14*exp(22.52*(T0/Temperature-1.0)-26.92*(1+log(T0/Temperature)-T0/Temperature));
      K0=1.805e-5/Keau; //mol/kg
      deltaH_over_RT0=-1.5;
      deltaCp0_over_R=26.92;
    }
  else if (name=="HNO3")
    {
      K0=12.0; //mol/kg
      deltaH_over_RT0=0.0; //29.17;
      deltaCp0_over_R=0.0; //16.83;
    }
  else if (name=="HCl")
    {
      K0=788; //mol/kg
      deltaH_over_RT0=0.0;
      deltaCp0_over_R=0.0;
    }
  else 
    cout << "error: inorganic precursor " << name << " not defined. " << endl;

  value=K0*exp(-deltaH_over_RT0*(T0/Temperature-1.0)-deltaCp0_over_R*(1+log(T0/Temperature)-T0/Temperature));
  return value;
}

void Kp_inorganic(model_config config, vector<species>& surrogate, double& Temperature, double& chp, double& MMaq)
{
  int n=surrogate.size();
  int i;
  double R=8.314; //ideal gas constant (J/K/mol)
  double rho_h2o=1000.0; //volumic mass of H2O
  double MH2O=18.0;
  double deltaH_over_RT0,deltaCp0_over_R;
  double T0=298.15;

  for (i=0;i<n;i++)
    if (surrogate[i].is_inorganic_precursor)
      {
        if (surrogate[i].name=="H2SO4")
          surrogate[i].Kaq_inorg=1.0e10;		
        else if (surrogate[i].name=="NH3")
          {
            deltaH_over_RT0=13.79;
            deltaCp0_over_R=-5.39;
            surrogate[i].Kaq_inorg=
              surrogate[i].Henry*exp(-deltaH_over_RT0*(T0/Temperature-1.0)-deltaCp0_over_R*(1+log(T0/Temperature)-T0/Temperature))
              *MH2O/MMaq*R*Temperature/(rho_h2o*1.0e6*1.013e5)
              *(1.0+surrogate[i].Kequilibrium(Temperature)*chp*surrogate[config.iHp].gamma_aq/surrogate[config.iNH4p].gamma_aq);
          }
        else if (surrogate[i].name=="HNO3")
          {
            deltaH_over_RT0=29.17;
            deltaCp0_over_R=16.83;
            surrogate[i].Kaq_inorg=
              surrogate[i].Henry*exp(-deltaH_over_RT0*(T0/Temperature-1.0)-deltaCp0_over_R*(1+log(T0/Temperature)-T0/Temperature))
              *R*Temperature/(rho_h2o*1.0e6*1.013e5)*MH2O/MMaq
              *(1.0+surrogate[i].Kequilibrium(Temperature)/(chp*surrogate[config.iHp].gamma_aq*surrogate[config.iNO3m].gamma_aq));
          }
        else if (surrogate[i].name=="HCl")
          {
            deltaH_over_RT0=30.20;
            deltaCp0_over_R=19.91;
            surrogate[i].Kaq_inorg=surrogate[i].Henry*exp(-deltaH_over_RT0*(T0/Temperature-1.0)-deltaCp0_over_R*(1+log(T0/Temperature)-T0/Temperature))
              *R*Temperature/(rho_h2o*1.0e6*1.013e5)*MH2O/MMaq
              *(1.0+surrogate[i].Kequilibrium(Temperature)/(chp*surrogate[config.iHp].gamma_aq*surrogate[config.iClm].gamma_aq));
          }
	    
      }
}

void Kpideal_inorganic(model_config config, vector<species>& surrogate, double& Temperature)
{
  int n=surrogate.size();
  int i;
  double R=8.314; //ideal gas constant (J/K/mol)
  double rho_h2o=1000.0; //volumic mass of H2O
  double MH2O=18.0;
  double deltaH_over_RT0,deltaCp0_over_R;
  double T0=298.15;

  for (i=0;i<n;i++)
    if (surrogate[i].is_inorganic_precursor)
      {
        if (surrogate[i].name=="H2SO4")
          surrogate[i].kpi=1.0e10;		
        else if (surrogate[i].name=="NH3")
          {
            deltaH_over_RT0=13.79;
            deltaCp0_over_R=-5.39;
            surrogate[i].kpi=
              surrogate[i].Henry*exp(-deltaH_over_RT0*(T0/Temperature-1.0)-deltaCp0_over_R*(1+log(T0/Temperature)-T0/Temperature))
              *MH2O*R*Temperature/(rho_h2o*1.0e6*1.013e5);
          }
        else if (surrogate[i].name=="HNO3")
          {
            deltaH_over_RT0=29.17;
            deltaCp0_over_R=16.83;
            surrogate[i].kpi=
              surrogate[i].Henry*exp(-deltaH_over_RT0*(T0/Temperature-1.0)-deltaCp0_over_R*(1+log(T0/Temperature)-T0/Temperature))
              *R*Temperature/(rho_h2o*1.0e6*1.013e5)*MH2O;
          }
        else if (surrogate[i].name=="HCl")
          {
            deltaH_over_RT0=30.20;
            deltaCp0_over_R=19.91;
            surrogate[i].kpi=surrogate[i].Henry*exp(-deltaH_over_RT0*(T0/Temperature-1.0)-deltaCp0_over_R*(1+log(T0/Temperature)-T0/Temperature))
              *R*Temperature/(rho_h2o*1.0e6*1.013e5)*MH2O;
          }
	    
      }
}
 
void Kpreal_inorganic(model_config config, vector<species>& surrogate, double& chp, double& MMaq)
{
  int n=surrogate.size();
  int i;
  double R=8.314; //ideal gas constant (J/K/mol)
  double rho_h2o=1000.0; //volumic mass of H2O
  double MH2O=18.0;
  double deltaH_over_RT0,deltaCp0_over_R;
  double T0=298.15;

  for (i=0;i<n;i++)
    if (surrogate[i].is_inorganic_precursor)
      {
        if (surrogate[i].name=="H2SO4")
          surrogate[i].Kaq_inorg=1.0e10;		
        else if (surrogate[i].name=="NH3")
          {
            deltaH_over_RT0=13.79;
            deltaCp0_over_R=-5.39;
            surrogate[i].Kaq_inorg=surrogate[i].kpi/MMaq*(1.0+surrogate[i].keq*chp*surrogate[config.iHp].gamma_aq/surrogate[config.iNH4p].gamma_aq);
          }
        else if (surrogate[i].name=="HNO3")
          {
            deltaH_over_RT0=29.17;
            deltaCp0_over_R=16.83;
            surrogate[i].Kaq_inorg=surrogate[i].kpi/MMaq*(1.0+surrogate[i].keq/(chp*surrogate[config.iHp].gamma_aq*surrogate[config.iNO3m].gamma_aq));
          }
        else if (surrogate[i].name=="HCl")
          {
            deltaH_over_RT0=30.20;
            deltaCp0_over_R=19.91;
            surrogate[i].Kaq_inorg=surrogate[i].kpi/MMaq*(1.0+surrogate[i].keq/(chp*surrogate[config.iHp].gamma_aq*surrogate[config.iClm].gamma_aq));
          }
	    
      }
}
